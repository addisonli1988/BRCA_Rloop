#!/usr/bin/env python3
"""Blacklist-filter Q20 fragments, call matched-Input candidates and annotate them."""
import argparse
from bisect import bisect_left
from collections import defaultdict
import gzip
from pathlib import Path
import re
import subprocess

CANONICAL={'chr'+str(i) for i in range(1,23)}|{'chrX','chrY'}
CHROMS=['chr'+str(i) for i in range(1,23)]+['chrX','chrY']

def require(ok, message):
    if not ok:
        raise ValueError(message)

def read_fai(path):
    d={}
    for line in Path(path).read_text().splitlines():
        fields=line.split('\t')
        if fields[0] in CANONICAL:
            require(fields[0] not in d and int(fields[1])>0,'Invalid FAI')
            d[fields[0]]=int(fields[1])
    require(set(d)==CANONICAL,'FAI must contain chr1..22,X,Y')
    return dict(sorted(d.items()))

class IntervalIndex:
    def __init__(self, intervals):
        grouped=defaultdict(list)
        for chrom,start,end in intervals:
            grouped[chrom].append((start,end))
        self.data={}
        for chrom,vs in grouped.items():
            merged=[]
            for start,end in sorted(vs):
                if merged and start<=merged[-1][1]:
                    merged[-1]=(merged[-1][0],max(end,merged[-1][1]))
                else:
                    merged.append((start,end))
            self.data[chrom]=([x[0] for x in merged],[x[1] for x in merged])

    def overlaps(self,chrom,start,end):
        starts,ends=self.data.get(chrom,([],[]))
        i=bisect_left(starts,end)-1
        return i>=0 and ends[i]>start

def peak_command(cfg,treatment,control,out,name):
    return [cfg['macs3'],'callpeak','-t',str(treatment),'-c',str(control),'-f','BEDPE','-g',str(cfg['effective_genome_size']),
            '-n',name,'--outdir',str(out),'-q',str(cfg['peak_qvalue']),'--keep-dup','all',
            '--min-length',str(cfg['peak_min_length']),'--max-gap',str(cfg['peak_max_gap'])]

def need(ok,message):
    if not ok: raise ValueError(message)

def merge(spans):
    out=[]
    for b,e in sorted(spans):
        need(0<=b<e,'Invalid interval')
        if out and b<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],e))
        else: out.append((b,e))
    return out

class Index:
    def __init__(self,spans):
        self.spans=merge(spans);self.starts=[v[0] for v in self.spans];self.ends=[v[1] for v in self.spans]
    def overlaps(self,b,e):
        i=bisect_left(self.starts,e)-1
        return i>=0 and self.ends[i]>b

def nearest(tss,positions,midpoint):
    if not positions:return None,[]
    i=bisect_left(positions,midpoint)
    choices=positions[max(0,i-1):min(len(positions),i+1)]
    distance=min(abs(p-midpoint) for p in choices)
    genes=sorted({g for p in choices if abs(p-midpoint)==distance for g in tss[p]})
    return distance,genes

def parse_gtf(path,sizes):
    bodies=defaultdict(list);promoters=defaultdict(list);tss=defaultdict(lambda:defaultdict(set));genes={};ntrans=0;headers=[]
    with gzip.open(path,'rt',encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):headers.append(line.rstrip());continue
            v=line.rstrip('\n').split('\t');need(len(v)==9,'Malformed GTF')
            c,_,feature,start,end,_,strand,_,attrs=v
            if c not in sizes or feature not in ('gene','transcript'):continue
            b,e=int(start)-1,int(end);need(0<=b<e<=sizes[c] and strand in ('+','-'),'GTF coordinates/build mismatch')
            a=dict(re.findall(r'(\w+) "([^"]*)"',attrs));gid=a['gene_id']
            if feature=='gene':
                need(gid not in genes,'Duplicate GTF gene ID')
                genes[gid]=dict(gene_id=gid,gene_name=a.get('gene_name',''),gene_type=a.get('gene_type',''),chrom=c,start=b,end=e,strand=strand)
                bodies[c].append((b,e))
            else:
                pos=b if strand=='+' else e-1
                promoters[c].append((max(0,pos-2000),min(sizes[c],pos+2001)))
                tss[c][pos].add(gid);ntrans+=1
    need(any('version 47' in h and 'GRCh38' in h for h in headers),'Expected GENCODE v47 GRCh38 header')
    need(genes and ntrans>0 and all(c in bodies and c in tss for c in CHROMS),'Missing chromosome annotations')
    need(all(g in genes for d in tss.values() for gs in d.values() for g in gs),'Transcript gene absent')
    return ({c:Index(promoters[c]) for c in CHROMS},{c:Index(bodies[c]) for c in CHROMS},
            tss,{c:sorted(tss[c]) for c in CHROMS},genes,dict(gtf_headers=headers,n_genes=len(genes),n_transcripts=ntrans,
            n_unique_tss=sum(len(x) for x in tss.values())))

def category(promoter,gene_body,b,e):
    p=promoter.overlaps(b,e);g=gene_body.overlaps(b,e)
    return ('promoter' if p else ('gene_body_nonpromoter' if g else 'intergenic')),int(p),int(g)


def open_text(path,mode='rt'):
    path=Path(path)
    if path.suffix=='.gz':return gzip.open(path,mode,encoding='utf-8')
    return path.open(mode,encoding='utf-8')

def bed(path,sizes):
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(('#','track','browser')):continue
            fields=line.split();chrom,start,end=fields[0],int(fields[1]),int(fields[2])
            if chrom not in sizes:continue
            require(0<=start<end<=sizes[chrom],'Invalid BED coordinates')
            yield chrom,start,end

def blacklist_filter(source,destination,sizes,mask):
    with open_text(destination,'xt') as output:
        for chrom,start,end in bed(source,sizes):
            if not mask.overlaps(chrom,start,end):output.write(f'{chrom}\t{start}\t{end}\n')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('ip','input','fai','blacklist','gtf','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--sample',required=True)
    p.add_argument('--macs3',default='macs3')
    a=p.parse_args();sizes=read_fai(a.fai);mask=IntervalIndex(bed(a.blacklist,sizes))
    a.output.mkdir(parents=True,exist_ok=False)
    ip=a.output/'IP.fragments_Q20.blacklist_filtered.bed'
    control=a.output/'Input.fragments_Q20.blacklist_filtered.bed'
    blacklist_filter(a.ip,ip,sizes,mask);blacklist_filter(a.input,control,sizes,mask)
    cfg={'macs3':a.macs3,'effective_genome_size':2913022398,'peak_qvalue':0.01,'peak_min_length':100,'peak_max_gap':50}
    subprocess.run(peak_command(cfg,ip,control,a.output,a.sample),check=True)
    peaks=sorted((c,b,e) for c,b,e in bed(a.output/f'{a.sample}_peaks.narrowPeak',sizes) if not mask.overlaps(c,b,e))
    promoters,bodies,tss,positions,genes,_=parse_gtf(a.gtf,sizes)
    destination=a.output/f'{a.sample}.candidate_regions_Q20.tsv.gz'
    with open_text(destination,'xt') as output:
        output.write('sample\tpeak_id\tchrom\tstart\tend\twidth_bp\tcategory\toverlaps_promoter\toverlaps_gene_body\tmidpoint_nearest_TSS_distance_bp\tnearest_gene_ids\tnearest_gene_names\n')
        for i,(chrom,start,end) in enumerate(peaks,1):
            label,promoter,body=category(promoters[chrom],bodies[chrom],start,end)
            distance,ids=nearest(tss[chrom],positions[chrom],(start+end-1)//2)
            row=[a.sample,f'small_Q20_{a.sample}_{i}',chrom,start,end,end-start,label,promoter,body,
                 distance,';'.join(ids),';'.join(genes[g]['gene_name'] for g in ids)]
            output.write('\t'.join('NA' if value is None else str(value) for value in row)+'\n')
    print(f'{a.sample}: {len(peaks)} candidates; {destination}')

if __name__=='__main__':main()
