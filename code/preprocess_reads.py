#!/usr/bin/env python3
"""Demultiplex pooled mDRIP-seq reads and process each IP or Input library."""
import argparse
from collections import Counter
from contextlib import ExitStack
import gzip
import json
from pathlib import Path
import re
import subprocess
import sys

FIXED=b'GTCGATCGCTG'
MOTIF=b'CAGCGATCGAC'
BARCODES={'S1': 'TGCA', 'S2': 'GTAC', 'S3': 'GCTA', 'S5': 'CGTA', 'S6': 'TGAC', 'S7': 'GTCA', 'S8': 'CTGA', 'S9': 'TCGA', 'S10': 'GATC', 'S11': 'AGTC', 'S12': 'CTAG', 'S13': 'TCAG', 'S14': 'CATG'}

def read_record(handle, digest=None):
    header = handle.readline()
    if not header:
        return None
    rows = (header, handle.readline(), handle.readline(), handle.readline())
    if digest is not None:
        for row in rows:
            digest.update(row)
    if not all(rows) or not all(row.endswith(b'\n') for row in rows):
        raise ValueError('Truncated FASTQ or missing final newline')
    h, s, p, q = (row.rstrip(b'\r\n') for row in rows)
    if not h.startswith(b'@') or not p.startswith(b'+') or not h[1:].split():
        raise ValueError('Invalid FASTQ header/separator')
    if not s or len(s) != len(q) or not re.fullmatch(b'[ACGTNacgtn]+', s):
        raise ValueError('Invalid FASTQ sequence or quality length')
    if not re.fullmatch(b'[!-~]+', q):
        raise ValueError('Invalid FASTQ quality byte (expected Phred+33 ASCII)')
    return rows

def read_id(record, mate):
    if mate not in (1, 2):
        raise ValueError('Mate must be 1 or 2')
    fields = record[0].rstrip(b'\r\n')[1:].split()
    token = fields[0]
    if token.endswith((b'/1', b'/2')):
        if token[-1:] != str(mate).encode():
            raise ValueError('Physical FASTQ mate suffix is reversed')
        token = token[:-2]
    if len(fields) > 1 and re.match(b'^[0-9]+:[YN]:', fields[1]):
        if fields[1].split(b':', 1)[0] != str(mate).encode():
            raise ValueError('CASAVA mate field is reversed')
    return token

def open_fastq(path):
    return gzip.open(path, 'rb') if str(path).endswith('.gz') else open(path, 'rb')

def trim_barcode(record):
    seq, qual = record[1].rstrip(b'\r\n'), record[3].rstrip(b'\r\n')
    if len(seq) <= 4:
        raise ValueError('Matched R2 has no sequence after barcode')
    # Normalize output separator; preserve physical read header and quality.
    return (record[0].rstrip(b'\r\n') + b'\n', seq[4:] + b'\n', b'+\n', qual[4:] + b'\n')

def reverse_complement(sequence):
    return sequence.translate(bytes.maketrans(b'ACGT', b'TGCA'))[::-1]

def paired_records(r1, r2):
    with open_fastq(r1) as a, open_fastq(r2) as b:
        while True:
            one, two = read_record(a), read_record(b)
            if one is None or two is None:
                if one is not None or two is not None:
                    raise ValueError('FASTQ mates differ in record count')
                return
            if read_id(one, 1) != read_id(two, 2):
                raise ValueError('FASTQ mate ID mismatch')
            yield one, two

def match_readthrough(sequence, adapter):
    """Full 15nt anywhere, or exact 11..14nt adapter prefix at physical 3prime.

    A complete internal motif with a different sample barcode is NOT trimmed.
    """
    sequence = sequence.upper()
    pos = sequence.find(adapter)
    if pos >= 0:
        return pos, 'full15'
    for length in range(14, 10, -1):
        if sequence.endswith(adapter[:length]):
            return len(sequence)-length, 'terminal_partial_' + str(length)
    return None, 'not_matched'


def demultiplex(a):
    a.output.mkdir(parents=True,exist_ok=False)
    lookup={v.encode():k for k,v in BARCODES.items()}
    counts=Counter()
    with ExitStack() as stack:
        writers={}
        for sample in [*BARCODES,'unassigned']:
            writers[sample]=[stack.enter_context(gzip.open(a.output/f'{sample}.R{mate}.fastq.gz','wb')) for mate in (1,2)]
        for one,two in paired_records(a.r1,a.r2):
            sample=lookup.get(two[1].rstrip(b'\r\n')[:4].upper(),'unassigned')
            if sample!='unassigned':two=trim_barcode(two)
            counts[sample]+=1
            for handle,record in zip(writers[sample],(one,two)):handle.writelines(record)
    print(json.dumps(dict(counts)))

def structure_trim(a, paths):
    adapter=MOTIF+reverse_complement(BARCODES[a.sample].encode())
    with ExitStack() as stack:
        writers=[stack.enter_context(gzip.open(p,'wb')) for p in paths]
        for one,two in paired_records(a.r1,a.r2):
            s1,q1,s2,q2=[r.rstrip(b'\r\n') for r in (one[1],one[3],two[1],two[3])]
            pos,_=match_readthrough(s1,adapter)
            if pos is not None:s1,q1=s1[:pos],q1[:pos]
            if min(len(s1),len(s2))<50:continue
            for handle,record,seq,qual in ((writers[0],one,s1,q1),(writers[1],two,s2,q2)):
                handle.write(record[0].rstrip(b'\r\n')+b'\n'+seq+b'\n+\n'+qual+b'\n')

def command(argv,log):
    log.write(' '.join(map(str,argv))+'\n');log.flush()
    subprocess.run(list(map(str,argv)),check=True,stdout=log,stderr=log)

def pipe(left,right,log):
    log.write(' '.join(map(str,left))+' | '+' '.join(map(str,right))+'\n');log.flush()
    with subprocess.Popen(list(map(str,left)),stdout=subprocess.PIPE,stderr=log) as process:
        try:result=subprocess.run(list(map(str,right)),stdin=process.stdout,stdout=log,stderr=log)
        finally:process.stdout.close()
        status=process.wait()
    if status or result.returncode:raise RuntimeError(f'Pipeline failed: {status}/{result.returncode}')

def library(a):
    if a.sample not in BARCODES:raise ValueError('Unknown sample barcode')
    a.output.mkdir(parents=True,exist_ok=False)
    structure=[a.output/f'structure.R{i}.fastq.gz' for i in (1,2)]
    structure_trim(a,structure)
    front=[a.output/f'front.R{i}.fastq.gz' for i in (1,2)]
    clean=[a.output/f'clean.R{i}.fastq.gz' for i in (1,2)]
    aligned=a.output/'aligned.sorted.bam'
    dedup=a.output/'deduplicated.bam'
    names=a.output/'deduplicated.queryname.bam'
    with (a.output/'commands.log').open('x',encoding='utf-8') as log:
        command([a.cutadapt,'-j',max(1,a.threads-1),'-G','^GTCGATCGCTG','--no-indels','-e','0.1',
                 '-o',front[0],'-p',front[1],*structure],log)
        command([a.cutadapt,'-j',max(1,a.threads-1),'-q','20','-m','50','--pair-filter=any',
                 '-a','AGATCGGAAGAGC','-A','AGATCGGAAGAGC','-o',clean[0],'-p',clean[1],*front],log)
        bowtie=[a.bowtie2,'--very-sensitive-local','--soft-clipped-unmapped-tlen','--fr','-I','0','-X','2000',
                '-p',max(1,a.threads-3),'-x',a.index,'--rg-id',f'{a.sample}.{a.assay}',
                '--rg',f'SM:{a.sample}','--rg',f'LB:{a.sample}.{a.assay}','--rg','PL:ILLUMINA',
                '--rg','PU:provided_pooled_files','-1',clean[0],'-2',clean[1]]
        pipe(bowtie,[a.samtools,'sort','-@','1','-m','1G','-o',aligned,'-'],log)
        picard=[a.java,'-jar',a.picard_jar] if a.picard_jar else [a.picard]
        command([*picard,'MarkDuplicates',f'I={aligned}',f'O={dedup}',f'M={a.output/"duplicate_metrics.txt"}',
                 'REMOVE_DUPLICATES=true','VALIDATION_STRINGENCY=STRICT',f'TMP_DIR={a.output}'],log)
        command([a.samtools,'sort','-n','-@',max(1,a.threads-1),'-m','1G','-o',names,dedup],log)
        pipe([a.samtools,'view','-h','-F','2304',names],
             [sys.executable,Path(__file__).with_name('filter_fragments.py'),'--output',a.output/'fragments_Q20.bed.gz'],log)

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='stage',required=True)
    d=sub.add_parser('demultiplex',help='Exact matching of the first four R2 bases')
    l=sub.add_parser('library',help='Structure trimming, alignment, deduplication and Q20 filtering')
    for q in (d,l):
        for flag in ('r1','r2','output'):q.add_argument('--'+flag,type=Path,required=True)
    d.set_defaults(func=demultiplex);l.set_defaults(func=library)
    l.add_argument('--sample',required=True,choices=sorted(BARCODES))
    l.add_argument('--assay',required=True,choices=['IP','Input'])
    l.add_argument('--index',required=True,help='Bowtie2 hg38 index prefix')
    l.add_argument('--threads',type=int,default=8)
    for flag in ('cutadapt','bowtie2','samtools','picard','java'):l.add_argument('--'+flag,default=flag)
    l.add_argument('--picard-jar',type=Path)
    a=p.parse_args();a.func(a)

if __name__=='__main__':main()
