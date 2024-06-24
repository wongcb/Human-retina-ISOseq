########Iso-seq profiling of novel isoform variants in human retina at single cell resolution############
########Created by Luozixian Wang, Raymond Wong##########################################################
########CERA, UNIMELB, 24/06/2024########################################################################

#need to install required packages such as isoseq3 and pbpigeon
conda activate isoseq

#step2 primer removal
lima Retina1_1.ccs.bam primers.fasta Retina1_1.fl.bam --isoseq --per-read
lima Retina1_2.ccs.bam primers.fasta Retina1_2.fl.bam --isoseq --per-read
lima Retina1_3.ccs.bam primers.fasta Retina1_3.fl.bam --isoseq --per-read
lima Retina2_1.ccs.bam primers.fasta Retina2_1.fl.bam --isoseq --per-read
lima Retina2_2.ccs.bam primers.fasta Retina2_2.fl.bam --isoseq --per-read
lima Retina2_3.ccs.bam primers.fasta Retina2_3.fl.bam --isoseq --per-read
lima Retina3_1.ccs.bam primers.fasta Retina3_1.fl.bam --isoseq --per-read
lima Retina3_2.ccs.bam primers.fasta Retina3_2.fl.bam --isoseq --per-read
lima Retina3_3.ccs.bam primers.fasta Retina3_3.fl.bam --isoseq --per-read

#step3 tag
isoseq tag Retina1_1.fl.5p--3p.bam Retina1_1.flt.bam --design T-12U-16B
isoseq tag Retina1_2.fl.5p--3p.bam Retina1_2.flt.bam --design T-12U-16B
isoseq tag Retina1_3.fl.5p--3p.bam Retina1_3.flt.bam --design T-12U-16B
isoseq tag Retina2_1.fl.5p--3p.bam Retina2_1.flt.bam --design T-12U-16B
isoseq tag Retina2_2.fl.5p--3p.bam Retina2_2.flt.bam --design T-12U-16B
isoseq tag Retina2_3.fl.5p--3p.bam Retina2_3.flt.bam --design T-12U-16B
isoseq tag Retina3_1.fl.5p--3p.bam Retina3_1.flt.bam --design T-12U-16B
isoseq tag Retina3_2.fl.5p--3p.bam Retina3_2.flt.bam --design T-12U-16B
isoseq tag Retina3_3.fl.5p--3p.bam Retina3_3.flt.bam --design T-12U-16B

#step4 refine - Remove poly(A) tails (default = 20bp polyA tail) and concatemer
isoseq refine Retina1_1.flt.bam primers.fasta Retina1_1.fltnc.bam --require-polya
isoseq refine Retina1_2.flt.bam primers.fasta Retina1_2.fltnc.bam --require-polya
isoseq refine Retina1_3.flt.bam primers.fasta Retina1_3.fltnc.bam --require-polya
isoseq refine Retina2_1.flt.bam primers.fasta Retina2_1.fltnc.bam --require-polya
isoseq refine Retina2_2.flt.bam primers.fasta Retina2_2.fltnc.bam --require-polya
isoseq refine Retina2_3.flt.bam primers.fasta Retina2_3.fltnc.bam --require-polya
isoseq refine Retina3_1.flt.bam primers.fasta Retina3_1.fltnc.bam --require-polya
isoseq refine Retina3_2.flt.bam primers.fasta Retina3_2.fltnc.bam --require-polya
isoseq refine Retina3_3.flt.bam primers.fasta Retina3_3.fltnc.bam --require-polya

#step4b merge SMRT cells
ls Retina1_1.fltnc.bam Retina1_2.fltnc.bam Retina1_3.fltnc.bam Retina2_1.fltnc.bam Retina2_2.fltnc.bam Retina2_3.fltnc.bam Retina3_1.fltnc.bam Retina3_2.fltnc.bam Retina3_3.fltnc.bam > RetinaMerge.fltnc.fofn

#step 5: Correct single cell barcodes based on an include list (percentile method) - set to 96 percentile based on kneeplot
isoseq correct --method percentile --percentile 96 --barcodes 3M-february-2018-REVERSE-COMPLEMENTED.txt.gz RetinaMerge.fltnc.fofn RetinaMerge.corrected.bam

#barcode statistics
isoseq bcstats --method percentile --percentile 96 --json RetinaMerge.bcstats.json -o RetinaMerge.bcstats.tsv RetinaMerge.corrected.bam

#step6 deduplication (use groupdedup)
samtools sort -t CB RetinaMerge.corrected.bam -o RetinaMerge.corrected.sorted.bam
isoseq groupdedup RetinaMerge.corrected.sorted.bam RetinaMerge.dedup.bam


#using pigeon to do classfication: Map reads using pbmm2 before collapsing
#step1 map reads to the reference genome 
pbmm2 align --preset ISOSEQ --sort RetinaMerge.dedup.bam GRCh38.p14.genome.fa RetinaMerge.dedup.mapped.bam

#step2 collapse into unique isoforms
isoseq collapse RetinaMerge.dedup.mapped.bam RetinaMerge.dedup.collapse.gff

#step3 sort input transcript GFF
pigeon sort RetinaMerged.dedup.mapped.collapsed.gff -o RetinaMerged.dedup.mapped.collapsed.sorted.gff

#step3b sort and index the reference files 

pigeon sort gencode.v44.chr_patch_hapl_scaff.annotation.gtf -o gencode.v44.chr_patch_hapl_scaff.annotation.sorted.gtf
pigeon index gencode.v44.chr_patch_hapl_scaff.annotation.sorted.gtf
pigeon sort refTSS_v3.3_human_coordinate.hg38.sorted.bed -o refTSS_v3.3_human_coordinate.hg38.sorted.bed
pigeon index refTSS_v3.3_human_coordinate.hg38.sorted.bed

#step4 classify isoforms
pigeon classify RetinaMerge.sorted.gff gencode.v44.chr_patch_hapl_scaff.annotation.sorted.gtf GRCh38.p14.genome.fa --fl RetinaMerge.dedup.collapse.abundance.txt  --cage-peak refTSS_v3.3_human_coordinate.hg38.sorted.bed --poly-a polyA.list.txt

#step5 filter isoforms
pigeon filter RetinaMerge_classification.txt --isoforms RetinaMerge.sorted.gff 

#step6 report gene saturation
pigeon report RetinaMerge_classification.filtered_lite_classification.txt RetinaMerge_saturation.txt

#step7 make seurat compatible input
pigeon make-seurat --dedup RetinaMerge.dedup.fasta --group RetinaMerge.dedup.collapse.group.txt -d RetinaMerge_classification.filtered_lite_classification.txt


####End of the session########################
####Author: Luozixian Wang, Raymond Wong######
####24/06/2024################################