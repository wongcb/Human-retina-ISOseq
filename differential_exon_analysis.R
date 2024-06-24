#Differential exon analysis
#do for each sample individually
# Have combined fastq.gz files for each sample in a folder called "fqFolder" in directory
#Follow ScisorSeqR pipeline until InfoPerLongRead when AllInfo files are made
#please see github at https://github.com/noush-joglekar/scisorseqr

#load in libraries
library(scisorseqr)

STARalign('fqFolder','/athena/tilgnerlab/store/hut2006/soft/src/star-mapper/2.5.2b/STAR-2.5.2b/bin/Linux_x86_64/STARlong','/athena/tilgnerlab/store/hut2006/data/seq/genomes/H.sapiens/GRCh38/wholeGenomeUnzipped/starIndex_gencode24_sequins/', 10)

MapAndFilter(annoGZ="/athena/tilgnerlab/store/hut2006/data/annotations/H.sapiens/GRCh38/gencode.v34.annotation.gtf.gz", numThreads=10, seqDir="/athena/tilgnerlab/store/hut2006/data/seq/genomes/H.sapiens/GRCh38/chromFa/", filterFullLength=TRUE, cageBed="/athena/tilgnerlab/scratch/anj2026/2020_11_30_snisorseqAnalysis/2021_04_06_ONTdataProcessing_BothRuns/updatedPhastcons/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz", polyABed="/athena/tilgnerlab/store/shardwick/genomes/human/atlas.clusters.2.0.GRCh38.96_chrNames.bed.gz", genomeVersion='hg38')

GetBarcodes(fqFolder="fqFolder/", BCClustAssignFile="../RetinaA_2_3/RetinaA23_forGetBarcodes_updated_withsample", filterReads=FALSE)

InfoPerLongRead(barcodeOutputFile="Barcodes/OutputFiltered/FilteredDeconvBC_AllFiles.csv", mapAndFilterOut="LRProcessingOutput/", minTimesIsoObserve=3, rmTmpFolder=FALSE)