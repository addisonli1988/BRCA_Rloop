#!/usr/bin/env python3
"""Export Q20 total nuclear fragments from queryname-sorted primary SAM pairs."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path
import re
import sys

CANONICAL_NUCLEAR = {'chr' + str(i) for i in range(1,23)} | {'chrX','chrY'}
MITOCHONDRIAL = {'chrM','MT'}
CIGAR_TOKEN = re.compile(r'(\d+)([MIDNSHP=X])')

class InputError(ValueError):
    pass

def contig_class(chrom):
    if chrom in CANONICAL_NUCLEAR:
        return 'canonical_nuclear'
    if chrom in MITOCHONDRIAL:
        return 'mitochondrial'
    return 'other_contigs'

def cigar_info(cigar, mapped):
    if cigar == '*':
        if mapped:
            raise InputError('Mapped SAM record has no CIGAR.')
        return dict(reference_length=0, softclip_total=0, softclip_left=0, softclip_right=0, aligned_bases=0)
    tokens = CIGAR_TOKEN.findall(cigar)
    if not tokens or ''.join(n + op for n, op in tokens) != cigar or any(int(n) <= 0 for n, _ in tokens):
        raise InputError('Invalid SAM CIGAR: ' + cigar)
    reference = sum(int(n) for n, op in tokens if op in 'MDN=X')
    if mapped and reference <= 0:
        raise InputError('Mapped CIGAR consumes no reference bases.')
    nonhard = [(int(n), op) for n, op in tokens if op != 'H']
    return dict(reference_length=reference, aligned_bases=sum(int(n) for n, op in tokens if op in 'M=X'),
                softclip_total=sum(int(n) for n, op in tokens if op == 'S'),
                softclip_left=nonhard[0][0] if nonhard and nonhard[0][1] == 'S' else 0,
                softclip_right=nonhard[-1][0] if nonhard and nonhard[-1][1] == 'S' else 0)

def parse_record(line, line_number):
    fields = line.rstrip('\r\n').split('\t')
    if len(fields) < 11:
        raise InputError('SAM line %d has fewer than 11 fields.' % line_number)
    if not fields[0] or fields[0] == '*' or any(c.isspace() for c in fields[0]):
        raise InputError('SAM line %d has an invalid QNAME.' % line_number)
    try:
        flag, position, mapq, mate_position, tlen = [int(fields[i]) for i in (1, 3, 4, 7, 8)]
    except ValueError as error:
        raise InputError('SAM line %d has a noninteger numeric field.' % line_number) from error
    if not 0 <= flag <= 65535 or not 0 <= mapq <= 255 or position < 0 or mate_position < 0:
        raise InputError('SAM numeric field is out of range at line %d.' % line_number)
    mapped = not flag & 4
    if mapped and (fields[2] == '*' or position < 1):
        raise InputError('Mapped record lacks reference/position at line %d.' % line_number)
    info = cigar_info(fields[5], mapped)
    return dict(qname=fields[0], flag=flag, chrom=fields[2], position=position,
                mapq=mapq, rnext=fields[6], pnext=mate_position, tlen=tlen,
                start=position - 1 if mapped else None,
                end=position - 1 + info['reference_length'] if mapped else None,
                **info)

def pair_decision(records, min_mapq, max_fragment, min_aligned_bases=0, canonical_only=False, reject_unknown_mapq=False):
    """Return (one rejection reason, None) or (None, accepted BED6 tuple)."""
    if len(records) != 2:
        return 'record_count_not_two', None
    a, b = records
    if any(r['flag'] & 2304 for r in records):
        return 'non_primary_record', None
    if any(not r['flag'] & 1 for r in records):
        return 'not_paired', None
    r1 = [r for r in records if r['flag'] & 64 and not r['flag'] & 128]
    r2 = [r for r in records if r['flag'] & 128 and not r['flag'] & 64]
    if len(r1) != 1 or len(r2) != 1:
        return 'r1_r2_identity_invalid', None
    if any(r['flag'] & 4 for r in records):
        return 'unmapped_record', None
    if any(r['flag'] & 8 for r in records):
        return 'mate_unmapped_flag', None
    if any(r['flag'] & 512 for r in records):
        return 'qc_fail', None
    if any(r['flag'] & 1024 for r in records):
        return 'duplicate', None
    if any(not r['flag'] & 2 for r in records):
        return 'not_proper_pair', None
    if any(r['mapq'] < min_mapq for r in records):
        return 'low_mate_mapq', None
    if a['chrom'] != b['chrom']:
        return 'different_contigs', None
    if bool(a['flag'] & 16) == bool(b['flag'] & 16):
        return 'same_orientation', None
    if bool(a['flag'] & 32) != bool(b['flag'] & 16) or bool(b['flag'] & 32) != bool(a['flag'] & 16):
        return 'mate_orientation_flag_inconsistent', None
    for left, right in ((a, b), (b, a)):
        mate_ref = left['chrom'] if left['rnext'] == '=' else left['rnext']
        if mate_ref != right['chrom']:
            return 'mate_reference_inconsistent', None
    if a['pnext'] != b['position'] or b['pnext'] != a['position']:
        return 'mate_position_inconsistent', None
    if a['tlen'] == 0 or a['tlen'] != -b['tlen']:
        return 'tlen_not_nonzero_opposites', None
    start, end = min(a['start'], b['start']), max(a['end'], b['end'])
    length = end - start
    if abs(a['tlen']) != length:
        return 'tlen_span_mismatch', None
    if length > max_fragment:
        return 'fragment_too_long', None
    if reject_unknown_mapq and any(r['mapq'] == 255 for r in records):
        return 'unknown_mate_mapq', None
    if any(r['aligned_bases'] < min_aligned_bases for r in records):
        return 'short_mate_aligned_bases', None
    if canonical_only and contig_class(a['chrom']) != 'canonical_nuclear':
        return 'noncanonical_contig', None
    strand = '-' if r1[0]['flag'] & 16 else '+'
    return None, (a['chrom'], start, end, a['qname'], 0, strand)


def run(stream, destination):
    destination=Path(destination)
    destination.parent.mkdir(parents=True,exist_ok=True)
    retained=0
    rejected=Counter()
    records=[]
    current=None
    overflow=False
    sort_order=None
    seen=False
    def finish(output):
        nonlocal retained
        if not records:return
        reason,fragment=('record_count_not_two',None) if overflow else pair_decision(
            records,min_mapq=20,max_fragment=2000,min_aligned_bases=50,
            canonical_only=True,reject_unknown_mapq=True)
        if reason:rejected[reason]+=1
        else:
            retained+=1
            output.write(f'{fragment[0]}\t{fragment[1]}\t{fragment[2]}\n')
    with destination.open('xb') as raw:
        with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as compressed:
            import io
            with io.TextIOWrapper(compressed,encoding='utf-8',newline='') as output:
                for number,line in enumerate(stream,1):
                    if not line.strip():raise InputError(f'Blank SAM line at {number}')
                    if line.startswith('@'):
                        if seen:raise InputError('SAM header follows alignments')
                        if line.startswith('@HD\t'):
                            fields=dict(item.split(':',1) for item in line.rstrip().split('\t')[1:] if ':' in item)
                            sort_order=fields.get('SO')
                        continue
                    if sort_order!='queryname':raise InputError('Use samtools sort -n; @HD must declare SO:queryname')
                    seen=True
                    record=parse_record(line,number)
                    if current is not None and record['qname']!=current:
                        finish(output);records=[];overflow=False
                    current=record['qname']
                    if len(records)<3:records.append(record)
                    else:overflow=True
                finish(output)
    if not seen:raise InputError('No alignments supplied')
    if not retained:raise InputError('No pairs passed the Q20 filters')
    print(json.dumps({'retained_fragments':retained,'rejected_pairs':dict(rejected)}))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True,help='New .bed.gz file')
    a=p.parse_args()
    run(sys.stdin,a.output)

if __name__=='__main__':main()
