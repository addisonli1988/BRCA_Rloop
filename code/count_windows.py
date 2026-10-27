#!/usr/bin/env python3
"""Build an independent 500-bp grid, count fragments and assemble IP/Input matrices."""
import argparse
from array import array
from contextlib import ExitStack
import csv
import itertools
from pathlib import Path
from call_regions import read_fai, IntervalIndex, bed, open_text

WIDTH=500

def need(ok,msg):
    if not ok:raise ValueError(msg)

def grid(sizes,mask,width):
    for c,length in sizes.items():
        for b in range(0,length-width+1,width):
            e=b+width
            if not mask.overlaps(c,b,e):yield f'W{width}_{c}_{b}',c,b,e,width,int(c not in ('chrX','chrY'))

def delta_arrays(sizes,width):
    return {c:array('q',[0])*((length+width-1)//width+1) for c,length in sizes.items()}

def add_fragment(deltas,c,b,e,width):
    deltas[c][b//width]+=1;deltas[c][(e-1)//width+1]-=1

def prefix_sum(deltas):
    for v in deltas.values():
        acc=0
        for i in range(len(v)):
            acc+=v[i];v[i]=acc
        need(acc==0,'Fragment delta did not balance')


def make_windows(a):
    sizes=read_fai(a.fai);mask=IntervalIndex(bed(a.blacklist,sizes))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with open_text(a.output,'xt') as output:
        output.write('window_id\tchrom\tstart\tend\tlength_bp\tautosomal\n')
        for row in grid(sizes,mask,WIDTH):output.write('\t'.join(map(str,row))+'\n')

def count(a):
    sizes=read_fai(a.fai);deltas=delta_arrays(sizes,WIDTH);n=0
    for chrom,start,end in bed(a.fragments,sizes):add_fragment(deltas,chrom,start,end,WIDTH);n+=1
    prefix_sum(deltas)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with open_text(a.windows) as windows,open_text(a.output,'xt') as output:
        output.write('window_id\tcount\n')
        for row in csv.DictReader(windows,delimiter='\t'):
            value=deltas[row['chrom']][int(row['start'])//WIDTH]
            output.write(f'{row["window_id"]}\t{value}\n')
    print(f'{n} retained nuclear fragments; fragment count is not the matrix-column sum')

def matrix(a):
    with a.metadata.open(encoding='utf-8-sig',newline='') as handle:
        samples=[r['sample'] for r in csv.DictReader(handle,delimiter='\t')]
    if not samples or len(samples)!=len(set(samples)):raise ValueError('Empty or duplicated samples')
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with ExitStack() as stack:
        readers=[csv.DictReader(stack.enter_context(open_text(a.counts_dir/f'{sample}.{a.assay}.counts.tsv.gz')),delimiter='\t') for sample in samples]
        output=stack.enter_context(open_text(a.output,'xt'))
        output.write('window_id\t'+'\t'.join(samples)+'\n')
        for rows in itertools.zip_longest(*readers):
            if any(r is None for r in rows):raise ValueError('Count files have different lengths')
            if len({r['window_id'] for r in rows})!=1:raise ValueError('Count files use different window orders')
            output.write(rows[0]['window_id']+'\t'+'\t'.join(r['count'] for r in rows)+'\n')

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='stage',required=True)
    w=sub.add_parser('windows');w.set_defaults(func=make_windows)
    for flag in ('fai','blacklist','output'):w.add_argument('--'+flag,type=Path,required=True)
    c=sub.add_parser('count');c.set_defaults(func=count)
    for flag in ('fai','fragments','windows','output'):c.add_argument('--'+flag,type=Path,required=True)
    m=sub.add_parser('matrix');m.set_defaults(func=matrix)
    for flag in ('metadata','counts-dir','output'):m.add_argument('--'+flag,type=Path,required=True)
    m.add_argument('--assay',choices=['IP','Input'],required=True)
    a=p.parse_args();a.func(a)

if __name__=='__main__':main()
