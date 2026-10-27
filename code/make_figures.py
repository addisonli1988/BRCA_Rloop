#!/usr/bin/env python3
"""Draw Figures 1-3 directly from the deposited candidate, metadata and library-QC tables."""
import argparse
import csv
import gzip
import os
from pathlib import Path
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--data',type=Path,default=Path('../processed_data_for_upload'))
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();DATA=a.data;OUT=a.output;OUT.mkdir(parents=True,exist_ok=False)
os.environ.setdefault('MPLCONFIGDIR',str(OUT/'.mplconfig'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, FancyBboxPatch
from matplotlib.ticker import LogLocator, LogFormatterMathtext

def table(path):
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',encoding='utf-8-sig',newline='') as handle:return list(csv.DictReader(handle,delimiter='\t'))

metadata=table(DATA/'sample_metadata.tsv')
metadata.sort(key=lambda row:({'adjacent':0,'benign':1,'tumor':2}[row['tissue_group']],int(row['sample'][1:])))
ORDER=[r['sample'] for r in metadata]
if len(ORDER)!=12 or len(set(ORDER))!=12:raise ValueError('Expected 12 distinct analysed specimens')
GROUPS=['adjacent','benign','tumor']
GROUP_LABELS={'adjacent':'Adjacent','benign':'Benign','tumor':'Tumour'}
COLORS={'adjacent':'#007C60','benign':'#B36500','tumor':'#0072B2'}
LIGHT={'adjacent':'#E6F1ED','benign':'#F5EBDD','tumor':'#E5F0F7'}
MARKERS={'adjacent':'o','benign':'s','tumor':'^'}
PAIRS=[]
for adjacent in (r for r in metadata if r['tissue_group']=='adjacent'):
    tumour=next(r for r in metadata if r['patient_id']==adjacent['patient_id'] and r['tissue_group']=='tumor')
    PAIRS.append((adjacent['patient_id'],adjacent['sample'],tumour['sample']))
qc={(r['assay'],r['sample']):r for r in table(DATA/'library_qc_Q20.tsv')}
if len(qc)!=24:raise ValueError('Expected 24 IP/Input libraries')
nwindows={int(r['autosomal_windows_total']) for r in qc.values()}
if len(nwindows)!=1:raise ValueError('Window universes differ')
N_WINDOWS=nwindows.pop()
support={key:{'zero':r['autosomal_windows_count0'],'one_to_four':r['autosomal_windows_count1_to4'],
               'ge5':r['autosomal_windows_count_ge5']} for key,r in qc.items()}
for row in support.values():
    if sum(int(v) for v in row.values())!=N_WINDOWS:raise ValueError('Window-support counts do not sum to the autosomal universe')
regions=table(DATA/'candidate_regions_Q20.tsv.gz')
candidate={}
for m in metadata:
    sample=m['sample'];rows=[r for r in regions if r['sample']==sample]
    if not rows:raise ValueError(f'No candidates for {sample}')
    widths=[int(r['width_bp']) for r in rows]
    q25,median,q75=np.percentile(widths,[25,50,75])
    candidate[sample]={**m,'n_candidates':len(rows),'p25_width_bp':q25,'median_width_bp':median,'p75_width_bp':q75,
        'ip_fragments':int(qc['IP',sample]['n_fragments']),'input_fragments':int(qc['Input',sample]['n_fragments'])}
    for category in ('promoter','gene_body_nonpromoter','intergenic'):
        candidate[sample][category+'_percent']=100*sum(r['category']==category for r in rows)/len(rows)

def ranks(values):
    order=np.argsort(values,kind='stable');result=np.empty(len(values),dtype=float);start=0
    while start<len(values):
        end=start+1
        while end<len(values) and values[order[end]]==values[order[start]]:end+=1
        result[order[start:end]]=(start+end-1)/2+1;start=end
    return result

associations={}
for field in ('ip_fragments','input_fragments'):
    first=ranks([candidate[s]['n_candidates'] for s in ORDER]);second=ranks([candidate[s][field] for s in ORDER])
    associations['n_candidates',field]={'spearman_rho':float(np.corrcoef(first,second)[0,1])}

plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.labelsize': 8,
    'axes.titlesize': 8.5, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
    'legend.fontsize': 7.5, 'axes.linewidth': .6, 'xtick.major.width': .6,
    'ytick.major.width': .6, 'axes.spines.top': False, 'axes.spines.right': False,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'savefig.facecolor': 'white', 'hatch.linewidth': .35})

def save(fig, stem, description):
    for extension in ('png','pdf'):
        fig.savefig(OUT/f'{stem}.{extension}',dpi=600,facecolor='white')
    plt.close(fig)

def panel(ax, letter, title, x=-.075):
    ax.set_title(title, loc='left', pad=10)
    ax.text(x, 1.08, letter, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom')

def sample_axis(ax):
    ax.set_xlim(-.65, 11.65)
    ax.set_xticks(np.arange(12), ORDER)
    for cut in [3.5, 6.5]:
        ax.axvline(cut, color='#AAAAAA', lw=.5, ls=(0, (2, 2)))

def box(ax, xy, w, h, text, fc='#F4F5F6', ec='#6E7479', fontsize=8, bold=False):
    p = FancyBboxPatch(xy, w, h, boxstyle='round,pad=0.008,rounding_size=0.015',
        linewidth=.7, facecolor=fc, edgecolor=ec)
    ax.add_patch(p)
    ax.text(xy[0] + w/2, xy[1] + h/2, text, ha='center', va='center',
        fontsize=fontsize, fontweight='bold' if bold else 'normal', linespacing=1.3)
    return p

def arrow(ax, start, end, color='#5F676C'):
    ax.annotate('', xy=end, xytext=start,
        arrowprops={'arrowstyle': '-|>', 'color': color, 'lw': .8, 'shrinkA': 1.7, 'shrinkB': 1.7,
                    'mutation_scale': 8})

# Figure 1: specimen relationships and pooled mDRIP-seq workflow.
fig = plt.figure(figsize=(180/25.4, 160/25.4))
ca = fig.add_axes([.045, .495, .92, .46])
ca.set(xlim=(0, 1), ylim=(0, 1)); ca.axis('off')
ca.text(0, 1.01, 'a', fontsize=10, fontweight='bold', va='bottom')
ca.text(.035, 1.01, 'Cohort: 12 breast-tissue specimens from 8 patients', fontsize=9, fontweight='bold', va='bottom')
ca.text(.245, .92, 'Four adjacent–tumour pairs', ha='center', fontsize=8.5)
ca.text(.765, .92, 'Four independent patients', ha='center', fontsize=8.5)
for i, (patient, adjacent, tumour) in enumerate(PAIRS):
    y = .71 - i*.16
    ca.text(.015, y+.053, patient, va='center', fontsize=8)
    box(ca, (.09, y), .165, .106, f'{adjacent}  Adjacent', fc=LIGHT['adjacent'], ec=COLORS['adjacent'])
    box(ca, (.325, y), .165, .106, f'{tumour}  Tumour', fc=LIGHT['tumor'], ec=COLORS['tumor'])
    ca.plot([.264, .316], [y+.053]*2, color='#777777', lw=.8)
independent = [(r['patient_id'],r['sample'],r['tissue_group']) for r in metadata if r['patient_id'] not in {p[0] for p in PAIRS}]
for i, (patient, sample, group) in enumerate(independent):
    col, row = i%2, i//2
    x, y = .56 + col*.225, .625 - row*.245
    box(ca, (x, y), .19, .15, f'{patient}\n{sample}  {GROUP_LABELS[group]}', fc=LIGHT[group], ec=COLORS[group])
ca.text(.765, .225, '4 adjacent  |  3 benign  |  5 tumour', ha='center', fontsize=8)
ca.text(.5, .055, 'Each specimen is represented by matched IP and Input components', ha='center', fontsize=8.5)

wa = fig.add_axes([.045, .055, .92, .385])
wa.set(xlim=(0, 1), ylim=(0, 1)); wa.axis('off')
wa.text(0, 1.025, 'b', fontsize=10, fontweight='bold', va='bottom')
wa.text(.035, 1.025, 'Experimental and computational workflow', fontsize=9, fontweight='bold', va='bottom')
box(wa, (.012, .66), .267, .22, 'mDRIP-seq\nBarcode-labelled pooled processing', fontsize=8)
box(wa, (.36, .66), .27, .22, 'IP and Input libraries\nPaired-end sequencing', fontsize=8)
box(wa, (.71, .66), .275, .22, 'Barcode demultiplexing\nAdapter and quality trimming', fontsize=8)
arrow(wa, (.282, .77), (.355, .77)); arrow(wa, (.633, .77), (.705, .77))
box(wa, (.71, .21), .275, .25, 'hg38 alignment and deduplication\nQ20 fragment filtering\nBlacklist exclusion', fontsize=8)
arrow(wa, (.848, .655), (.848, .465))
box(wa, (.36, .21), .27, .25, 'IP enrichment against matched Input\nCandidate regions\nGENCODE annotation', fc='#EDF2F6', fontsize=8)
box(wa, (.012, .21), .267, .25, 'Fixed 500-bp windows\nSeparate IP and Input\nfragment-count matrices', fc='#EDF2F6', fontsize=8)
# A shared filtered-fragment node feeds two independent output representations.
wa.plot([.70, .66, .66, .145, .145], [.335, .335, .11, .11, .20], color='#5F676C', lw=.8)
arrow(wa, (.70, .335), (.635, .335))
arrow(wa, (.145, .115), (.145, .205))
wa.text(.45, .025, 'Candidate enrichment and fixed-window counts provide complementary representations',
        ha='center', va='center', fontsize=7.5)
save(fig, 'Figure_1_cohort_workflow', 'Cohort design, matched IP/Input and confirmed barcode-pooled experimental/computational processing.')

# Figure 2: library depth and autosomal 500-bp window support.
fig, axes = plt.subplots(3, 1, figsize=(180/25.4, 155/25.4))
fig.subplots_adjust(left=.115, right=.98, bottom=.09, top=.94, hspace=.72)
x = np.arange(12)
ax = axes[0]
for assay, delta, color, hatch in [('IP', -.19, '#0072B2', ''), ('Input', .19, '#C6893D', '///')]:
    vals = [int(candidate[s]['ip_fragments' if assay == 'IP' else 'input_fragments']) / 1e6 for s in ORDER]
    ax.bar(x+delta, vals, width=.35, color=color, edgecolor='#333333', linewidth=.45, hatch=hatch, label=assay, zorder=3)
ax.set_ylim(0, 20); ax.set_yticks([0, 5, 10, 15, 20]); ax.set_ylabel('Fragments (million)')
ax.grid(axis='y', color='#E0E0E0', lw=.5); ax.set_axisbelow(True)
sample_axis(ax); panel(ax, 'a', 'Filtered fragment depth')
ax.legend(frameon=False, ncol=2, loc='upper right', bbox_to_anchor=(1, 1.26), handlelength=1.5)
support_encoding = [('zero', '0', '#E4E4E4', ''), ('one_to_four', '1–4', '#8CBAD4', ''), ('ge5', '≥5', '#215B81', '///')]
for ax, assay, letter in zip(axes[1:], ['IP', 'Input'], ['b', 'c']):
    bottom = np.zeros(12)
    for col, label, color, hatch in support_encoding:
        values = np.array([100*int(support[assay, s][col])/N_WINDOWS for s in ORDER])
        ax.bar(x, values, bottom=bottom, width=.72, color=color, edgecolor='#333333', linewidth=.35, hatch=hatch, label=label)
        bottom += values
    assert np.max(np.abs(bottom-100)) < 1e-10
    ax.set_ylim(0, 100); ax.set_yticks([0, 25, 50, 75, 100]); ax.set_ylabel('Windows (%)')
    sample_axis(ax); panel(ax, letter, f'{assay}: 500-bp count support')
    if assay == 'IP':
        ax.legend(title='Raw fragments per window', title_fontsize=7.5, frameon=False, ncol=3,
            loc='upper right', bbox_to_anchor=(1, 1.34), handlelength=1.4, columnspacing=1.2)
axes[-1].set_xlabel('Specimen')
for start, end, label in [(0, 3, 'Adjacent'), (4, 6, 'Benign'), (7, 11, 'Tumour')]:
    axes[-1].text((start+end)/2, -.37, label, transform=axes[-1].get_xaxis_transform(), ha='center', fontsize=7.5)
save(fig, 'Figure_2_depth_window_support', 'Filtered-fragment library depths and raw-count support on 5,376,024 identical autosomal 500-bp windows, all 12 specimens.')

# Figure 3: individual candidates, called-interval widths and annotation, then depth associations.
fig = plt.figure(figsize=(180/25.4, 190/25.4))
aa = fig.add_axes([.085, .495, .245, .385])
ab = fig.add_axes([.405, .495, .205, .385])
ac = fig.add_axes([.675, .495, .295, .385])
for i, sample in enumerate(ORDER):
    r = candidate[sample]; group = r['tissue_group']; count = int(r['n_candidates'])
    aa.scatter(count, i, marker=MARKERS[group], s=23, color=COLORS[group], edgecolor='#222222', linewidth=.5, zorder=3)
    aa.annotate(f'{count:,}', (count, i), xytext=(4, 0), textcoords='offset points', va='center', fontsize=7.2)
    q25, med, q75 = [float(r[k]) for k in ['p25_width_bp', 'median_width_bp', 'p75_width_bp']]
    ab.plot([q25, q75], [i, i], color=COLORS[group], lw=1.3, zorder=2)
    ab.plot([q25, q25], [i-.11, i+.11], color=COLORS[group], lw=.7)
    ab.plot([q75, q75], [i-.11, i+.11], color=COLORS[group], lw=.7)
    ab.scatter(med, i, marker=MARKERS[group], s=23, color=COLORS[group], edgecolor='#222222', linewidth=.5, zorder=3)
aa.set_xscale('log', base=10); aa.set_xlim(40, 300000)
aa.set_xticks([100, 1000, 10000, 100000]); aa.xaxis.set_major_formatter(LogFormatterMathtext(base=10))
aa.xaxis.set_minor_locator(LogLocator(base=10, subs=[])); aa.set_xlabel('Candidates (log10 scale)')
ab.set_xlim(100, 355); ab.set_xticks([100, 200, 300]); ab.set_xlabel('Width (bp), median and IQR')
composition = [('promoter_percent', 'Promoter', '#45957F', ''),
    ('gene_body_nonpromoter_percent', 'Gene body (non-promoter)', '#79A9CB', '///'),
    ('intergenic_percent', 'Intergenic', '#D2D2D2', '...')]
left = np.zeros(12)
for col, label, color, hatch in composition:
    vals = np.array([float(candidate[s][col]) for s in ORDER])
    ac.barh(np.arange(12), vals, left=left, height=.68, color=color, edgecolor='#333333', linewidth=.35, hatch=hatch, label=label)
    left += vals
assert np.max(np.abs(left-100)) < 1e-8
ac.set_xlim(0, 100); ac.set_xticks([0, 50, 100]); ac.set_xlabel('Candidate records (%)')
for ax, letter, title in [(aa, 'a', 'Candidate numbers'), (ab, 'b', 'Called-interval widths'), (ac, 'c', 'Genomic composition')]:
    ax.set_ylim(11.65, -.65); ax.set_yticks(np.arange(12), ORDER if ax is aa else ['']*12)
    ax.tick_params(axis='y', length=3 if ax is aa else 0)
    for cut in [3.5, 6.5]:
        ax.axhline(cut, color='#AAAAAA', lw=.5, ls=(0, (2, 2)))
    if ax is not ac:
        ax.grid(axis='x', color='#E0E0E0', lw=.5); ax.set_axisbelow(True)
    panel(ax, letter, title, x=-.075)
composition_handles = [Patch(facecolor=color, edgecolor='#333333', linewidth=.4, hatch=hatch, label=label)
    for _, label, color, hatch in composition]
fig.legend(handles=composition_handles, frameon=False, fontsize=7.5, loc='upper center',
    bbox_to_anchor=(.54, .974), ncol=3, handlelength=1.6, handletextpad=.5, columnspacing=1.1)
ad = fig.add_axes([.11, .105, .355, .28])
ae = fig.add_axes([.62, .105, .355, .28])
label_positions = {
    'IP': {'S2': (.25, .045), 'S11': (.05, .105), 'S1': (.29, .105), 'S6': (.07, .18), 'S14': (.36, .18),
        'S13': (.06, .27), 'S9': (.29, .25), 'S7': (.29, .34), 'S10': (.63, .42), 'S12': (.39, .69), 'S3': (.63, .79), 'S5': (.84, .89)},
    'Input': {'S2': (.24, .045), 'S11': (.04, .115), 'S1': (.28, .115), 'S6': (.04, .20), 'S14': (.29, .19),
        'S13': (.29, .27), 'S9': (.05, .29), 'S7': (.29, .36), 'S10': (.38, .43), 'S12': (.44, .68), 'S3': (.66, .79), 'S5': (.84, .89)}}
for ax, field, assay, letter, max_depth in [(ad, 'ip_fragments', 'IP', 'd', 19), (ae, 'input_fragments', 'Input', 'e', 4.5)]:
    ax.set_xlim(0, max_depth); ax.set_ylim(35, 200000); ax.set_yscale('log', base=10)
    ax.set_yticks([100, 1000, 10000, 100000]); ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=[]))
    ax.set_xlabel(f'{assay} fragments (million)'); ax.set_ylabel('Candidates (log10 scale)')
    ax.grid(color='#E0E0E0', lw=.5); ax.set_axisbelow(True)
    panel(ax, letter, f'Candidate count and {assay} depth')
    rho = float(associations['n_candidates', field]['spearman_rho'])
    ax.text(.03, .96, f'Spearman ρ = {rho:.3f}', transform=ax.transAxes, va='top', fontsize=7.5,
        bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .9, 'pad': 1.5})
    for sample in ORDER:
        r = candidate[sample]; xv, yv = int(r[field])/1e6, int(r['n_candidates'])
        ax.scatter(xv, yv, s=23, marker=MARKERS[r['tissue_group']], color=COLORS[r['tissue_group']], edgecolor='#222222', lw=.5, zorder=3)
        pos = label_positions[assay][sample]
        ax.text(*pos, sample, transform=ax.transAxes, fontsize=7.5, ha='center', va='center', zorder=5,
            bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .8, 'pad': .1})
        ax.annotate('', xy=(xv,yv), xycoords='data', xytext=pos, textcoords='axes fraction',
            arrowprops={'arrowstyle': '-', 'color': '#999999', 'lw': .4, 'shrinkA': 5, 'shrinkB': 3}, zorder=1)
handles = [Line2D([], [], marker=MARKERS[g], markerfacecolor=COLORS[g], markeredgecolor='#222222',
    markeredgewidth=.5, ls='', markersize=4.7, label=GROUP_LABELS[g]) for g in GROUPS]
fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.55, .016), ncol=3,
    frameon=False, handletextpad=.5, columnspacing=1.7)
save(fig, 'Figure_3_candidate_features', 'All specimen-specific candidate counts, median/IQR called-interval widths, genomic composition and descriptive IP/Input depth associations.')


print(f'Saved Figures 1-3 in {OUT}')
