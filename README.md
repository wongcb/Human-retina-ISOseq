# Single-cell Isoform Profiling of Human Retina

This repository contains the computational workflow used in the study:

**Iso-seq profiling of novel isoform variants in human retina at single cell resolution**

The analysis integrates **PacBio Iso-Seq long-read sequencing** with **short-read single-cell RNA sequencing (scRNA-seq)** to characterize transcript isoform diversity, alternative splicing, promoter usage, disease-associated isoforms, coding potential, and splice-variant support across retinal cell types.

The workflow combines multiple tools and custom scripts, including:

- IsoSeq3
- Pigeon
- Seurat
- scisorSeqR / ScisorWiz
- ASprofile
- proActiv
- ORFanage / Pfam / HMMER
- SpliceAI
- R, Python, and shell scripts

---

## Analysis workflow

The full analysis consists of several independent modules.

```text
PacBio Iso-Seq preprocessing following the IsoSeq3 pipeline
↓
Isoform collapsing and classification with Pigeon
↓
Integration with short-read scRNA-seq in Seurat
↓
High-confidence isoform filtering
↓
Cell-type-specific promoter and splicing analysis
↓
ORF prediction and Pfam domain annotation
↓
ClinVar variants + SpliceAI prediction
↓
Splice junction, isoform, end-bias validation
↓
R-based figure and table generation
```

Because several steps depend on large intermediate files and external tools, **each module should be run separately** rather than sourcing or executing the entire workflow at once.

---

## Repository structure

```text
.
├── README.md
│
├── R/
│   └── combined_R_pipeline.R
│
├── cmd/
│   ├── split_bam.sh
│   ├── AS_profile.sh
│   ├── run_orfanage.sh
│   ├── run_pfam_hmmscan_parallel.sh
│   ├── run_end_bias.sh
│   └── scisorseqr_scisorwiz.txt
│
└── python/
    ├── split_bam.sh
    ├── AS_profile.sh
    └── ABCA4_spliceAI/
        ├── extract_ABCA4_files.py
        ├── summarize_abca4_events.py
        ├── build_abca4_spliceai_table.py
        ├── clinvar_abca4_to_vcf.py
        ├── match_spliceai_to_isoform.py
        └── filter_spliceAI_support.py

```

---

## Software requirements

### Long-read processing

- IsoSeq3
- pbmm2
- Pigeon
- samtools

### Command-line tools

- ASprofile
- ORFanage
- gffread
- HMMER: `hmmscan`, `hmmpress`
- Pfam-A HMM database
- bedtools
- awk
- perl
- SpliceAI

### R packages

- Seurat
- SingleCellExperiment
- scDblFinder
- Matrix
- tidyverse / dplyr / tidyr
- data.table
- ggplot2
- cowplot
- ggtranscript
- rtracklayer
- GenomicRanges
- proActiv
- scisorseqr
- ScisorWiz
- WGCNA
- clusterProfiler
- patchwork

### Python packages

- pandas
- pathlib
- gzip
- collections

---

## Input data

Main input raw and intermediate files include:

- PacBio HiFi reads
- PacBio Iso-Seq collapsed isoform GTF files
- Pigeon classification tables
- Isoform count matrices
- short-read scRNA-seq Seurat object
- GENCODE v44 GRCh38.p14 GFF3 annotation
- GRCh38.p14 genome FASTA
- CAGE peak annotation
- polyA motif annotation
- ClinVar ABCA4 variant table
- selected ABCA4 Sanger validation alignment files

Large raw data files can be found at The European Genome-phenome Archive (EGA): EGAD50000002227. Intermediate files are generated using different modules following the pipeline and are not included in this repository. Paths in scripts should be updated before running.

---

## Step 1 — IsoSeq3 preprocessing and Pigeon classification

Raw PacBio HiFi reads were processed using the IsoSeq3 pipeline.

Key steps include:

- primer removal
- barcode tagging
- poly(A) trimming
- concatemer removal
- barcode correction
- deduplication
- genome alignment
- isoform collapsing
- isoform classification with Pigeon

Example commands:

```bash
lima Retina1_1.ccs.bam primers.fasta Retina1_1.fl.bam --isoseq
isoseq refine Retina1_1.flt.bam primers.fasta Retina1_1.fltnc.bam
pbmm2 align --preset ISOSEQ GRCh38.genome.fa RetinaMerge.dedup.bam RetinaMerge.dedup.mapped.bam
```

Outputs include:

- filtered long-read BAM files
- collapsed isoform annotation files
- Pigeon classification tables
- filtered isoform GTF files
- isoform count matrices

---

## Step 2 — scisorSeqR and ScisorWiz processing

Command record:

```text
cmd/scisorseqr_scisorwiz.txt
```

This module records the initial scisorSeqR and ScisorWiz commands.

Main steps include:

```r
STARalign(...)
MapAndFilter(...)
GetBarcodes(...)
InfoPerLongRead(...)
IsoQuant(...)
ExonQuant(...)
DiffSplicingAnalysis(...)
triHeatmap(...)
ScisorWiz_AllInfo(...)
```

Outputs include:

- `AllInfo`
- IsoQuant results
- ExonQuant results
- differential splicing outputs
- ScisorWiz gene-level visualization outputs

---

## Step 3 — Integration with short-read scRNA-seq

Long-read isoform counts were integrated with short-read scRNA-seq data using Seurat.

Main steps include:

- loading isoform and gene count matrices
- creating Seurat objects
- matching long-read and short-read cell barcodes
- constructing isoform assays
- assigning cell type metadata
- generating UMAP and cell-type-level summary plots

Example:

```r
cur <- CreateSeuratObject(gene_list)
cur[["isoform"]] <- CreateAssayObject(counts = iso_list)
```

Outputs include:

- integrated Seurat objects
- isoform assays
- cell type annotations
- Figure 1 and supplementary overview plots

---

## Step 4 — FL2+CAGE+polyA high-confidence subset

Script:

```bash
bash cmd/split_bam.sh
```

Purpose:

This module extracts high-confidence PBids that pass:

```text
Full length reads >= 2
within CAGE peak
polyA motif present
```

The script then maps these PBids to molecule IDs and filters BAM files using `samtools view -N`.

Main outputs:

- `pbids_FL2_CAGE_polyA.txt`
- `molecule_ids_FL2_CAGE_polyA.txt`
- `*_FL2_CAGE_polyA.bam`
- `*_FL2_CAGE_polyA.bam.bai`

These BAM files are used for proActiv and promoter-related downstream analyses.

---

## Step 5 — Promoter usage analysis with proActiv

Promoter activity was quantified using proActiv.

Main steps include:

- splitting BAM files by cell type or high-confidence subset
- generating promoter annotations
- quantifying promoter activity
- identifying alternative promoter usage
- plotting major/minor promoter patterns across retinal cell types

Example:

```r
result <- proActiv(
  files = bam_files,
  promoterAnnotation = promoterAnnotation,
  condition = cell_types
)
```

Outputs include:

- promoter activity matrices
- promoter usage summary tables
- Figure 2 promoter plots
- supplementary promoter related plots

---

## Step 6 — ASprofile alternative splicing analysis

Script:

```bash
bash cmd/AS_profile.sh
```

Purpose:

This module constructs a strict FL2+CAGE+polyA GFF subset and runs ASprofile.

Main commands:

```bash
extract-as RetinaMerge.sorted.filtered_lite.FL2_CAGE_polyA_strict.gff GRCh38.p14.genome.fa > AS_output_FL2CAGEpolyAmotif.txt
perl summarize_as.pl RetinaMerge.sorted.filtered_lite.FL2_CAGE_polyA_strict.gff AS_output_FL2CAGEpolyAmotif.txt -p Retina
```

Main outputs:

- `pbids_FL2_CAGE_polyA_strict.txt`
- `RetinaMerge.sorted.filtered_lite.FL2_CAGE_polyA_strict.gff`
- `AS_output_FL2CAGEpolyAmotif.txt`
- ASprofile summary files

These outputs are used for AS event summary and Figure 3-related analyses.

---

## Step 7 — ORFanage ORF prediction

Script:

```bash
bash cmd/run_orfanage.sh
```

Purpose:

This module predicts ORFs/CDS/protein sequences for filtered Iso-Seq transcripts.

Main tools:

- ORFanage
- gffread

Main outputs:

- `RetinaMerge.filtered_lite.orfanage.gtf`
- `RetinaMerge.filtered_lite.orfanage.stats.tsv`
- `RetinaMerge.filtered_lite.transcripts.fa`
- `RetinaMerge.filtered_lite.orfanage.cds.fa`
- `RetinaMerge.filtered_lite.orfanage.protein.fa`
- `run.log`

These files are used for coding potential and translation impact analyses.

---

## Step 8 — Pfam domain annotation

Script:

```bash
bash cmd/run_pfam_hmmscan_parallel.sh
```

Purpose:

This module runs Pfam domain annotation on ORFanage protein sequences using parallel `hmmscan`.

Main outputs:

- `pfam_analysis_parallel/pfam.domtblout`
- `pfam_analysis_parallel/pfam.tblout`
- `pfam_analysis_parallel/pfam_run_summary.txt`
- chunk-level hmmscan output files

These files are used to summarize protein domain support and protein consequence categories.

---

## Step 9 — Use spliceAI to predict the effect of ABCA4 ClinVar variants on splicing events

Scripts:

```text
python/ABCA4_spliceAI/
```

Purpose:

This module evaluates whether the spliceAI predicted splicing events resulted from ClinVar ABCA4 variants support splice-altering events observed in selected Iso-Seq PBids.

Recommended run order:

```bash
# 1. Extract selected ABCA4 PBids
python python/ABCA4_spliceAI/extract_ABCA4_files.py

# 2. Build ABCA4 isoform/event tables
python python/ABCA4_spliceAI/summarize_abca4_events.py
python python/ABCA4_spliceAI/build_abca4_spliceai_table.py

# 3. Convert ClinVar ABCA4 variants to SpliceAI-compatible VCF
python python/ABCA4_spliceAI/clinvar_abca4_to_vcf.py

# 4. Run SpliceAI
spliceai -I clinvar_spliceai_input/clinvar_ABCA4.spliceai_input.vcf \
         -O clinvar_spliceai_output/clinvar_ABCA4.spliceai_output.vcf \
         -R /path/to/GRCh38.p14.genome.fa \
         -A grch38

# 5. Match SpliceAI predictions to Iso-Seq splice events
python python/ABCA4_spliceAI/match_spliceai_to_isoform.py

# 6. Filter strong support events
python python/ABCA4_spliceAI/filter_spliceAI_support.py
```

Main outputs:

- `extracted_pbids/selected_pbids.fa`
- `extracted_pbids/selected_pbids.gff`
- `ABCA4_spliceai_ready/abca4_spliceai_comparison_events.tsv`
- `clinvar_spliceai_input/clinvar_ABCA4.spliceai_input.vcf`
- `clinvar_spliceai_output/clinvar_ABCA4.spliceai_output.vcf`
- `ABCA4_spliceai_match/abca4_spliceai_variant_event_matches.tsv`
- `ABCA4_spliceai_match/abca4_spliceai_best_matches.tsv`
- `ABCA4_spliceai_match/abca4_spliceai_variant_summary.tsv`
- `ABCA4_spliceai_match/strong_support_summary/abca4_strong_support_all.tsv`
- `ABCA4_spliceai_match/strong_support_summary/abca4_very_strong_support.tsv`

These results are used by the R pipeline for ABCA4 SpliceAI support summaries and figure generation.

---

## Step 10 — Isoform validation

This module summarizes validation methods used to support isoform structures, splice junctions, and transcript completeness.

Validation evidence includes:

- FLNC read 5′/3′ end-bias analysis
- RT-PCR and Sanger validation
- RJunBase splice junction support
- bulk or short-read RNA-seq splice junction support

### 10.1 FLNC end-bias validation

Script:

```bash
bash cmd/run_end_bias.sh
```

Purpose:

This analysis evaluates whether FLNC read 5′ and 3′ ends are located near annotated TSS and TES positions, supporting full-length transcript capture.

Main tools:

- samtools
- bedtools
- awk

Main outputs:

- `tss.sorted.bed`
- `tes.sorted.bed`
- `flnc.sorted.bam`
- `reads.sorted.bed`
- `read5p_vs_tss.closest.txt`
- `read3p_vs_tes.closest.txt`
- `read_best_gene.tsv`

These files are used by the R pipeline to generate FLNC end-bias validation plots.

### 10.2 RT-PCR and Sanger validation

Purpose:

Selected isoforms and splice junctions were validated experimentally using RT-PCR followed by Sanger sequencing.

Main inputs:

- selected PBid sequences
- primer sequence table
- Sanger alignment FASTA files
- Iso-Seq GFF annotation
- reference genome annotation

Main outputs:

- merged Sanger alignment files
- primer position tables
- PBid-level validation plots
- position-level support summaries
- consensus FASTA files

These outputs are used by the R pipeline to visualize validated isoform regions and compare Sanger-supported regions with Iso-Seq transcript structures.

### 10.3 RJunBase splice junction support

Script:

```bash
python python/validation/validation_SJ_RJunBase.py
```

Purpose:

Iso-Seq splice junctions were compared with RJunBase junction annotations to evaluate whether Iso-Seq junctions are supported by external splice junction database.

Main inputs:

- `RetinaMerge_classification.filtered_lite_classification.txt`
- `RetinaMerge_classification.filtered_lite_junctions.txt`
- `detail_LS_annotation.txt`

Main outputs:

- `FSM_ISM_junction_validation_strict.csv`
- `NNC_NIC_junction_validation_strict.csv`
- `FSM_ISM_junction_validation_loose.csv`
- `NNC_NIC_junction_validation_loose.csv`
- `SJ_match_summary.csv`

### 10.4 Short-read splice junction support

Script:

```bash
python python/validation/validation_SJ_SRbulk.py
```

Purpose:

Iso-Seq splice junctions were compared with merged short-read splice junction counts to evaluate independent RNA-seq support for splice junctions.

Main inputs:

- `RetinaMerge_classification.filtered_lite_classification.txt`
- `RetinaMerge_classification.filtered_lite_junctions.txt`
- `merged_SJ_counts.tsv`

Main outputs:

- `ALL_FSM_ISM_junction_shortread_support.csv`
- `ALL_NNC_NIC_junction_shortread_support.csv`
- `FL2_FSM_ISM_junction_shortread_support.csv`
- `FL2_NNC_NIC_junction_shortread_support.csv`
- `HighConfident_FSM_ISM_junction_shortread_support.csv`
- `HighConfident_NNC_NIC_junction_shortread_support.csv`
- `SJ_shortread_support_summary.csv`

The high-confidence group is defined as isoforms with FL ≥ 2, within CAGE peak, and a valid polyA motif.

These validation outputs are used by the R pipeline to summarize external support for Iso-Seq splice junctions.

---

## Step 11 — R downstream analysis and figure generation

Main script:

```r
source("R/combined_R_pipeline.R")
```

The R script contains multiple independent modules and should be run section by section.

Main R analysis sections include:

- cell type and isoform overview
- supplementary tables
- proActiv promoter analysis
- AS event analysis
- RetNet disease isoform analysis
- translation impact analysis
- SpliceAI analysis
- isoform validation
- figure generation

Important note:

Because some helper functions and object names are reused across modules, it is recommended to run one module at a time and verify its output before running the next module.

---

## Main figure outputs

The final figure labels and output filenames are documented in `R/combined_R_pipeline.R`.

Major figure groups include:

| Figure group | Analysis |
|---|---|
| Figure 1 | cell type overview, Iso-Seq read support, and isoform summary |
| Figure 2 | promoter usage and proActiv analysis |
| Figure 3 | alternative splicing and translation impact |
| Figure 4 | alternative exon usage in retinal cells |
| Figure 5 | RetNet disease-associated isoform analysis |
| Figure 6 | ABCA4 transcript structure and SpliceAI support |
| Supplementary figures | QC, read support, CAGE/polyA, and validation analyses |

---

## Suggested run order

```text
1. IsoSeq3 and Pigeon preprocessing
2. scisorSeqR / ScisorWiz processing
3. Seurat integration and isoform assay construction
4. split_bam.sh
5. AS_profile.sh
6. run_orfanage.sh
7. run_pfam_hmmscan_parallel.sh
8. ABCA4 ClinVar-SpliceAI Python module
9. isoform validation module
10. combined_R_pipeline.R modules
```

---

## Reproducibility notes

- Run each module separately.
- Check each module's report or summary output before continuing.
- Record software versions for IsoSeq3, Pigeon, SpliceAI, ORFanage, HMMER/Pfam, ASprofile, samtools, bedtools, R, and Python.
- Keep intermediate `.tsv`, `.csv`, `.vcf`, `.gff`, `.gtf`, `.fa`, `.bam`, and report files.
- The R plotting module should be run after all external preprocessing outputs have been generated.

---

## Notes on paths

Several scripts currently use local WSL paths such as:

```text
/mnt/d/bioinformatic/isoseq/processed/from_Ray/
```

Before sharing or rerunning the workflow, update these paths or define a common project root variable.

---

## Documentation

This README serves as the main workflow documentation for the repository. It describes the overall analysis structure, required software, major input/output files, and the recommended run order for each module.

---

## Authors

Mr. Luozixian Wang  
Prof. Raymond Wong  

Centre for Eye Research Australia  
University of Melbourne
