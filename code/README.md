# Breast-tissue mDRIP-seq analysis

Five scripts implement the study's Q20 total-fragment analysis and independent 500-bp windows. Use Python 3.11+; figure generation additionally requires `pip install -r requirements.txt`. External programs are Cutadapt 5.1, Bowtie2 2.5.4, SAMtools 1.20, Picard 3.4.0 and MACS3 3.0.3. Provide the hg38 Bowtie2 index/FAI, ENCODE hg38 blacklist v2 and comprehensive GENCODE v47 GRCh38 GTF (`.gtf.gz`). Run from this folder and use new output paths.

## Reads and candidate regions

```bash
python preprocess_reads.py demultiplex --r1 IP.R1.fastq.gz --r2 IP.R2.fastq.gz --output results/demux_IP
python preprocess_reads.py library --sample S1 --assay IP --r1 results/demux_IP/S1.R1.fastq.gz --r2 results/demux_IP/S1.R2.fastq.gz --index reference/hg38 --output results/S1_IP
python call_regions.py --sample S1 --ip results/S1_IP/fragments_Q20.bed.gz --input results/S1_Input/fragments_Q20.bed.gz --fai reference/hg38.fa.fai --blacklist reference/hg38-blacklist.v2.bed --gtf reference/gencode.v47.annotation.gtf.gz --output results/S1_regions
```

Repeat demultiplexing for Input and library processing for both components of each specimen. Thirteen exact four-base R2 barcodes are recognised; S8 is not included in the 12-specimen analysis. R1 read-through removal uses the exact 15-nt sample adapter or its terminal 11–14-nt prefix. No composition-based tail trimming is applied. `library` invokes `filter_fragments.py`: both primary mates must be proper, consistent, nonduplicate nuclear alignments with MAPQ≥20 (excluding 255), ≥50 aligned bases and a reference span≤2,000 bp. RNA-strand identity is not assigned.

`call_regions.py` removes whole fragments intersecting the blacklist, uses matched Input for MACS3 BEDPE calling (q=0.01, effective genome size 2,913,022,398, minimum length 100 bp, maximum gap 50 bp, `--keep-dup all`, default small-library scaling), and annotates retained candidates with promoter priority (TSS±2,000 bp). It writes the filtered IP/Input fragments needed below and one annotated candidate table per specimen. Concatenate these candidate tables with a single header to reproduce the deposited combined candidate table.

## Fixed-window counts

```bash
python count_windows.py windows --fai reference/hg38.fa.fai --blacklist reference/hg38-blacklist.v2.bed --output results/windows_500bp_hg38.tsv.gz
python count_windows.py count --fai reference/hg38.fa.fai --fragments results/S1_regions/IP.fragments_Q20.blacklist_filtered.bed --windows results/windows_500bp_hg38.tsv.gz --output results/counts/S1.IP.counts.tsv.gz
python count_windows.py matrix --metadata ../processed_data_for_upload/sample_metadata.tsv --counts-dir results/counts --assay IP --output results/IP_counts_Q20_500bp.tsv.gz
```

Repeat counts for all 12 specimens and both assays, naming count files `<sample>.<assay>.counts.tsv.gz`, then build IP and Input matrices separately. Windows are nonoverlapping 500-bp intervals on chr1–22/X/Y; short terminal and blacklist-overlapping windows are excluded. Each fragment contributes once to every window it overlaps by ≥1 bp. Matrix entries are raw integers; column sums are not fragment-library denominators.

## Figures

```bash
python make_figures.py --data ../processed_data_for_upload --output results/figures
```

This reads `sample_metadata.tsv`, `library_qc_Q20.tsv` and `candidate_regions_Q20.tsv.gz` directly. It recalculates candidate counts, width quartiles, annotation proportions and descriptive depth associations, and exports the three manuscript figures as PNG and PDF. Window-support summaries use the identical autosomal 500-bp windows; candidate counts use standard nuclear chromosomes.
