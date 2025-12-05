########Iso-seq profiling of novel isoform variants in human retina at single cell resolution############
########created by Raymond Wong, Luozixian Wang##########################################################
########CERA, UNIMELB, 05/12/2025########################################################################

#load the packages
library(AnnotationHub)
library(BiocManager)
library(Biostrings)
library(BSgenome)
library(Seurat)
library(tidyverse)
library(harmony)
library(Matrix)
library(patchwork)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(sctransform)
library(SingleR)
library(scran)
library(SingleCellExperiment)
library(scmap)
library(singleCellNet)
library(WGCNA)
library(igraph)
library(hdWGCNA)
library(qlcMatrix)
library(magrittr)
library(biomaRt)
library(org.Hs.eg.db)
library(ggrepel)
library(plotrix)
library(gridExtra)
library(clusterProfiler)
library(topGO)
library(proActiv)
library(data.table)
library(dunn.test)
library(ggplot2)
library(ggtranscript)
library(scales)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(umap)
library(scater)
library(uwot)
library(gtable)
library(purrr)
theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 14)
organism = "org.Hs.eg.db"


#read the processed short read dataset
sr_retina <- readRDS("seuratObj_samples234_filt_logNorm_umap_celltypes_reannotRW_211121.RDS")
#prepare the long read Seurat object
#directories with the counts matrices for each replicate
sample_dirs <- 'seurat_pb_output'
#load isoform counts matrices for each sample:
iso_list <- Seurat::Read10X('seurat_output/isoforms_seurat/', gene.column = 1)
#load gene counts matrices for each sample:
gene_list <- Seurat::Read10X('seurat_output/genes_seurat/')
#create individual Seurat objects for each sample
cur <- Seurat::CreateSeuratObject(gene_list, min.cells=1)
#add a column for the original barcode
cur$barcode <- colnames(cur)
#add the iso assay to the seurat object to contain the isoform counts matrix
cur[["iso"]] <- Seurat::CreateAssayObject(counts = iso_list, min.cells=1)
cur
#merge replicates into one Seurat object
lr_retina <- cur
#check basic info about this seurat object
lr_retina

#get reverse complemented cell barcode
bc_celltype <- read.csv("bc_celltype_assignments.csv")
sequences <- bc_celltype$Barcode
reverse_complement_sequences <- reverseComplement(DNAStringSet(sequences))
new_data <- data.frame(Sequence = sequences, Reverse_Complement = as.character(reverse_complement_sequences))
bc_celltype$Barcode <- new_data$Reverse_Complement
bc_celltype
#set up a celltype metadata
lr_retina$cell_type <- NA
barcode_celltype <- bc_celltype
seurat_barcode <- lr_retina$barcode
seurat_barcode <- gsub("-1$", "", seurat_barcode)
lr_retina$barcodes_trimmed <- seurat_barcode
#barcode deconvolution
for (i in 1:nrow(barcode_celltype)) {
  barcode <- barcode_celltype$Barcode[i]
  celltype <- barcode_celltype$CellType[i]
  lr_retina$cell_type[lr_retina$barcodes_trimmed %in% barcode] <- celltype
}
print(lr_retina)
#subset with barcoded cells
Idents(lr_retina) <- lr_retina$cell_type
iso_retina_annotated <- subset(lr_retina, idents = c("Amacrine", "Astrocyte", 
                                                      "Bipolar", "Cone", 
                                                      "Fibroblast", "Horizontal", 
                                                      "Muller_Glia", "Rod", 
                                                      "RPE", "Smooth_Muscle_Cells"))
#set the read status in both datasets
sr_retina@meta.data$read_stat <- "SR"
barcodes <- iso_retina_annotated@meta.data$barcodes_trimmed
rc_barcodes <- reverseComplement(DNAStringSet(barcodes))
rc_bc <- data.frame(Reverse_Complement = as.character(rc_barcodes))
rc_bc$read_stat <- "LR"
#read_status deconvolution
for (i in 1:nrow(rc_bc)) {
  barcode <- rc_bc$Reverse_Complement[i]
  read_stat <- rc_bc$read_stat[i]
  sr_retina$read_stat[sr_retina$cell.bc.trimmed %in% barcode] <- read_stat
}
table(sr_retina$read_stat)

#load the longread data as iso assay in the seurat object
iso_counts <- iso_retina_annotated@assays$iso@counts
iso_retina <- iso_counts
colnames_split <- strsplit(colnames(iso_retina), "-", fixed = TRUE)
colnames_no_suffix <- sapply(colnames_split, `[`, 1)
suffixes <- sapply(colnames_split, function(x) if (length(x) > 1) paste0("-", x[2]) else "")
colnames_rc <- sapply(DNAStringSet(colnames_no_suffix), function(x) as.character(reverseComplement(x)))
colnames_final <- paste0(colnames_rc, suffixes)
colnames(iso_retina) <- colnames_final
cell_barcodes_iso_list <- colnames(iso_retina)
check_cell_barcodes <- sr_retina@meta.data$cell.bc
check_duplicated_barcodes <- check_cell_barcodes[duplicated(check_cell_barcodes) | duplicated(check_cell_barcodes, fromLast = TRUE)]
print(unique(check_duplicated_barcodes))
#add the iso assay to the sr_retina
merged_retina <- sr_retina
cellbc_colnames <- merged_retina@meta.data$cell.bc

#delete the replicates
full_colnames <- colnames(merged_retina)
sr_barcodes_without_prefix <- merged_retina@meta.data$cell.bc
sr_barcode_to_first_colname <- list()
for (barcode in unique(sr_barcodes_without_prefix)) {
  matched_colnames <- full_colnames[grepl(barcode, full_colnames)]
  if (length(matched_colnames) > 0) {
    sr_barcode_to_first_colname[[barcode]] <- matched_colnames[1]
  }
}
colnames_to_keep <- unname(unlist(sr_barcode_to_first_colname))
merged_retina_subset <- subset(merged_retina, cells = colnames_to_keep)
cell_barcodes_sr_seurat <- merged_retina_subset@meta.data$cell.bc

#find common barcode
common_cells <- intersect(cell_barcodes_sr_seurat, cell_barcodes_iso_list)
common_cells

#replace the rownames in iso_retina with the full name from lr_retina
cell_bc_to_full_colname <- setNames(rownames(merged_retina_subset@meta.data), merged_retina_subset@meta.data$cell.bc)
new_colnames <- sapply(colnames(iso_retina), function(bc) cell_bc_to_full_colname[bc])

#check NA
if(any(is.na(new_colnames))) {
  print("have NA，some iso_retina cells cannot find corresponding colnames in merged_retina")
} else {
  print("all colnames mapped")
  colnames(iso_retina) <- new_colnames
}

cell_barcodes_sr_seurat <- colnames(merged_retina_subset)
cell_barcodes_iso_list <- colnames(iso_retina)
common_cells <- intersect(cell_barcodes_sr_seurat, cell_barcodes_iso_list)
common_cells
merged_retina_subset <- subset(merged_retina_subset, cells = common_cells)

#add the new assay
new_assay <- CreateAssayObject(counts = iso_retina)
merged_retina_subset[["isoform"]] <- new_assay
retina_combined <- merged_retina_subset

#add a new iso_FL2 assay containing the PBids with FL>=2
iso_fl2_classification <- read.delim("RetinaMerge_classification.filtered_lite_FL2_classification.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
length(iso_fl2_classification)
head(iso_fl2_classification)
iso_fl2_ids <- iso_fl2_classification$isoform

#check the cross
iso_counts_full <- GetAssayData(retina_combined, assay = "isoform", slot = "counts")
all_iso_full <- rownames(iso_counts_full)[1:10]
all_iso_full
all_iso_base <- sub(":.*$", "", rownames(iso_counts_full))
head(cbind(full = rownames(iso_counts_full), base = all_iso_base))

keep_rows <- all_iso_base %in% iso_fl2_ids
table(keep_rows)
iso_fl2_counts <- iso_counts_full[keep_rows, ]
dim(iso_fl2_counts)

retina_combined[["iso_FL2"]] <- CreateAssayObject(counts = iso_fl2_counts)
DefaultAssay(retina_combined) <- "iso_FL2"

saveRDS(retina_combined, "iso_sr_merged_retina_FL2.RDS")
saveRDS(merged_retina_subset, "iso_sr_merged_retina.RDS")
saveRDS(merged_retina, "iso_sr_merged_raw_retina.RDS")

#annotate the cell with read status
transcript_annotation <- read.csv(
  "seurat_pb_output/RetinaMerge.annotated.info.csv",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

transcript_status <- transcript_annotation %>%
  filter(category %in% c(
    "full-splice_match",
    "incomplete-splice_match",
    "novel_not_in_catalog",
    "novel_in_catalog"
  ))

cell_status <- transcript_status %>%
  group_by(BCrev) %>%
  summarise(
    has_FSM = any(category == "full-splice_match"),
    has_ISM = any(category == "incomplete-splice_match"),
    has_NNC = any(category == "novel_not_in_catalog"),
    has_NIC = any(category == "novel_in_catalog"),
    .groups = "drop"
  ) %>%
  mutate(
    has_known = has_FSM | has_ISM,
    has_novel = has_NNC | has_NIC
  ) %>%
  mutate(
    class_all = case_when(
      has_known & !has_novel ~ "Known_only(FSM+ISM)",
      !has_known & has_novel ~ "Novel_only(NNC+NIC)",
      has_known & has_novel  ~ "Mixed",
      TRUE                   ~ "None"
    )
  )

retina_combined$iso_class_all <- NA_character_
idx_all <- match(cell_status$BCrev, retina_combined$cell.bc.trimmed)
valid <- !is.na(idx_all)
retina_combined$iso_class_all[idx_all[valid]] <- cell_status$class_all[valid]
DimPlot(retina_combined, group.by ='iso_class_all', reduction='umap', pt.size = 0.3)

#FL2 read status
fl2_isoforms <- unique(iso_fl2_classification$isoform)
transcript_status_FL2 <- transcript_status %>%
  filter(pbid %in% fl2_isoforms)

cell_status_FL2 <- transcript_status_FL2 %>%
  group_by(BCrev) %>%
  summarise(
    has_FSM_FL2 = any(category == "full-splice_match"),
    has_ISM_FL2 = any(category == "incomplete-splice_match"),
    has_NNC_FL2 = any(category == "novel_not_in_catalog"),
    has_NIC_FL2 = any(category == "novel_in_catalog"),
    .groups = "drop"
  ) %>%
  mutate(
    has_known  = has_FSM_FL2 | has_ISM_FL2,
    has_novel  = has_NNC_FL2 | has_NIC_FL2
  ) %>%
  mutate(
    class_FL2 = case_when(
      has_known & !has_novel ~ "Known_only(FSM+ISM)",
      !has_known & has_novel ~ "Novel_only(NNC+NIC)",
      has_known & has_novel  ~ "Mixed",
      TRUE                   ~ "none"
    )
  )

retina_combined$iso_class_FL2 <- NA_character_
idx_FL2 <- match(cell_status_FL2$BCrev, retina_combined$cell.bc.trimmed)
valid_FL2 <- !is.na(idx_FL2)

retina_combined$iso_class_FL2[idx_FL2[valid_FL2]] <- cell_status_FL2$class_FL2[valid_FL2]

DimPlot(retina_combined, group.by ='iso_class_FL2', reduction='umap', pt.size = 0.8)

##plot them
#figure 1E
figure_1E_umap <- DimPlot(retina_combined, group.by ='iso_class_FL2', reduction='umap', pt.size = 0.7)
ggsave("figure1E.tiff", figure_1E_umap, width = 10, height = 8, dpi = 600, bg = "WHITE")

#figure 1B
figure_1B <- DimPlot(retina_combined, group.by ='cell.type.integrated', reduction='umap', pt.size = 0.5) + 
  theme(legend.position="bottom")
ggsave("figure1B.tiff", figure_1B, width = 8, height = 8, dpi = 600, bg = "WHITE")

#figure s1A
figureS1A <- DimPlot(merged_retina, reduction = "umap", group.by = "sample.id", pt.size = 0.5)
#figure s1B
figureS1B <- DimPlot(merged_retina, group.by ='read_stat', reduction='umap', pt.size = 0.5)

combined_fS1AB <- figureS1A | figureS1B
ggsave("figureS1AB.tiff", combined_fS1AB, width = 16, height = 8, dpi = 600, bg = "WHITE")

##GO term enrichment
symbol_ensid <- read.delim('symbol_ensid.txt', header = TRUE, stringsAsFactors = FALSE)
selected_symbol_ensid <- symbol_ensid[, c("HGNC.symbol", "Gene.stable.ID")]
filtered_symbol_ensid <- selected_symbol_ensid %>%
  group_by(HGNC.symbol) %>%
  slice(1) %>%
  ungroup()

FSM_transcript_FL2 <- transcript_status_FL2 %>%
  filter(category == "full-splice_match")

ISM_transcript_FL2 <- transcript_status_FL2 %>%
  filter(category == "incomplete-splice_match")

NNC_transcript_FL2 <- transcript_status_FL2 %>%
  filter(category == "novel_not_in_catalog")

NIC_transcript_FL2 <- transcript_status_FL2 %>%
  filter(category == "novel_in_catalog")

gene_FSM_FL2 <- data.frame(gene_FSM = unique(FSM_transcript_FL2$gene),
                           stringsAsFactors = FALSE)
gene_ISM_FL2 <- data.frame(gene_ISM = unique(ISM_transcript_FL2$gene),
                           stringsAsFactors = FALSE)
gene_NNC_FL2 <- data.frame(gene_NNC = unique(NNC_transcript_FL2$gene),
                           stringsAsFactors = FALSE)
gene_NIC_FL2 <- data.frame(gene_NIC = unique(NIC_transcript_FL2$gene),
                           stringsAsFactors = FALSE)

FSM_converted_FL2 <- gene_FSM_FL2 %>%
  left_join(filtered_symbol_ensid, by = c(gene_FSM = "HGNC.symbol"))

ISM_converted_FL2 <- gene_ISM_FL2 %>%
  left_join(filtered_symbol_ensid, by = c(gene_ISM = "HGNC.symbol"))

NNC_converted_FL2 <- gene_NNC_FL2 %>%
  left_join(filtered_symbol_ensid, by = c(gene_NNC = "HGNC.symbol"))

NIC_converted_FL2 <- gene_NIC_FL2 %>%
  left_join(filtered_symbol_ensid, by = c(gene_NIC = "HGNC.symbol"))

FSM_filled_FL2 <- FSM_converted_FL2 %>%
  mutate(Gene.stable.ID = ifelse(is.na(Gene.stable.ID),
                                 gene_FSM, Gene.stable.ID))

ISM_filled_FL2 <- ISM_converted_FL2 %>%
  mutate(Gene.stable.ID = ifelse(is.na(Gene.stable.ID),
                                 gene_ISM, Gene.stable.ID))

NNC_filled_FL2 <- NNC_converted_FL2 %>%
  mutate(Gene.stable.ID = ifelse(is.na(Gene.stable.ID),
                                 gene_NNC, Gene.stable.ID))

NIC_filled_FL2 <- NIC_converted_FL2 %>%
  mutate(Gene.stable.ID = ifelse(is.na(Gene.stable.ID),
                                 gene_NIC, Gene.stable.ID))


#check duplicate
FSM_duplicated_genes_FL2 <- FSM_filled_FL2 %>%
  filter(duplicated(Gene.stable.ID) | duplicated(Gene.stable.ID, fromLast = TRUE))

ISM_duplicated_genes_FL2 <- ISM_filled_FL2 %>%
  filter(duplicated(Gene.stable.ID) | duplicated(Gene.stable.ID, fromLast = TRUE))

NNC_duplicated_genes_FL2 <- NNC_filled_FL2 %>%
  filter(duplicated(Gene.stable.ID) | duplicated(Gene.stable.ID, fromLast = TRUE))

NIC_duplicated_genes_FL2 <- NIC_filled_FL2 %>%
  filter(duplicated(Gene.stable.ID) | duplicated(Gene.stable.ID, fromLast = TRUE))

FSM_genes_FL2_for_GO <- unique(FSM_filled_FL2$Gene.stable.ID)
ISM_genes_FL2_for_GO <- unique(ISM_filled_FL2$Gene.stable.ID)
NNC_genes_FL2_for_GO <- unique(NNC_filled_FL2$Gene.stable.ID)
NIC_genes_FL2_for_GO <- unique(NIC_filled_FL2$Gene.stable.ID)

#GO enrichment
FSM_GO_FL2 <- enrichGO(
  gene          = FSM_genes_FL2_for_GO,
  ont           = "ALL",
  keyType       = "ENSEMBL",
  pvalueCutoff  = 0.05,
  OrgDb         = "org.Hs.eg.db",
  pAdjustMethod = "fdr",
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

ISM_GO_FL2 <- enrichGO(
  gene          = ISM_genes_FL2_for_GO,
  ont           = "ALL",
  keyType       = "ENSEMBL",
  pvalueCutoff  = 0.05,
  OrgDb         = "org.Hs.eg.db",
  pAdjustMethod = "fdr",
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

NNC_GO_FL2 <- enrichGO(
  gene          = NNC_genes_FL2_for_GO,
  ont           = "ALL",
  keyType       = "ENSEMBL",
  pvalueCutoff  = 0.05,
  OrgDb         = "org.Hs.eg.db",
  pAdjustMethod = "fdr",
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

NIC_GO_FL2 <- enrichGO(
  gene          = NIC_genes_FL2_for_GO,
  ont           = "ALL",
  keyType       = "ENSEMBL",
  pvalueCutoff  = 0.05,
  OrgDb         = "org.Hs.eg.db",
  pAdjustMethod = "fdr",
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

#plot GO
pFSM_dot_FL2 <- dotplot(FSM_GO_FL2, showCategory = 10) + ggtitle("FSM (FL2)")
pISM_dot_FL2 <- dotplot(ISM_GO_FL2, showCategory = 10) + ggtitle("ISM (FL2)")
pNNC_dot_FL2 <- dotplot(NNC_GO_FL2, showCategory = 10) + ggtitle("NNC (FL2)")
pNIC_dot_FL2 <- dotplot(NIC_GO_FL2, showCategory = 10) + ggtitle("NIC (FL2)")

#figure S1CDEF
print(pFSM_dot_FL2)
print(pISM_dot_FL2)
print(pNNC_dot_FL2)
print(pNIC_dot_FL2)
combined_fS1CF <- (pFSM_dot_FL2 | pISM_dot_FL2) /
  (pNNC_dot_FL2 | pNIC_dot_FL2)
ggsave("figureS1CDEF.tiff",
       combined_fS1CF,
       width = 15, height = 14,
       dpi = 600, bg = "WHITE")


#alternative promoter analysis using proActiv
#before all, run a split_bam.py to split the bam file generated by isoseq3
#following is the script --> need to run in python

#coding: utf-8

import pysam
import csv

#load the barcode-celltype
barcode_to_celltype = {}
with open("rc_bc_celltype.csv", "r") as csv_file:
  reader = csv.reader(csv_file)
for row in reader:
  barcode_to_celltype[row[0]] = row[1]

#import the original bam file
input_bam = pysam.AlignmentFile("RetinaMerge.dedup.mapped.bam", "rb")

#split bam files according to corresponding cell types
celltype_to_bamfile = {}
for celltype in set(barcode_to_celltype.values()):
  celltype_to_bamfile[celltype] = pysam.AlignmentFile(f"{celltype}.bam", "wb", header=input_bam.header)
for read in input_bam:
  if read.has_tag("CB"):
  barcode = read.get_tag("CB")
if barcode in barcode_to_celltype:
  celltype = barcode_to_celltype[barcode]

celltype_to_bamfile[celltype].write(read)
input_bam.close()
for bamfile in celltype_to_bamfile.values():
  bamfile.close()


#run proActiv
#prepare promoter annotation
promoterAnnotation.gencode.v44.subset <- preparePromoterAnnotation(file = 'gencode.v44.chr_patch_hapl_scaff.annotation.gtf' , species = 'Homo_sapiens')

#run proActiv
files <- c("RayMerged_bamfile/Amacrine_FL2.bam", 
           "RayMerged_bamfile/Astrocyte_FL2.bam", 
           "RayMerged_bamfile/Bipolar_FL2.bam", 
           "RayMerged_bamfile/Cone_FL2.bam", 
           "RayMerged_bamfile/Fibroblast_FL2.bam", 
           "RayMerged_bamfile/Horizontal_FL2.bam",
           "RayMerged_bamfile/Muller_Glia_FL2.bam",
           "RayMerged_bamfile/RPE_FL2.bam", 
           "RayMerged_bamfile/Rod_FL2.bam", 
           "RayMerged_bamfile/Smooth_Muscle_Cells_FL2.bam")
promoterAnnotation <- promoterAnnotation.gencode.v44.subset
condition <- c('Amacrine', 'Astrocyte', 'Bipolar', 'Cone', 'Fibroblast', 'Horizontal', 'Muller_Glia', 'RPE', 'Rod', 'Smooth_Muscle_Cells')
result <- proActiv(files = files,
                   promoterAnnotation = promoterAnnotation.gencode.v44.subset,
                   condition = condition,
                   genome = "BSgenome.Hsapiens.UCSC.hg38")
show(result)
#filter out promoter counts that are NA
result <- result[complete.cases(assays(result)$promoterCounts),]

#determine alternative promoter
alternativePromoters <- getAlternativePromoters(result = result, referenceCondition = "Rod")
show(alternativePromoters)

#plotting:
plots <- boxplotPromoters(result, "ENSG00000198691.14") #ensembl id for ABCA4
# Boxplot of absolute promoter activity
grid.arrange(plots[[1]], plots[[3]], nrow = 1, ncol = 2, widths = c(3, 2))
#plot ratio
rdata <- rowData(result)

#redefine a function to do dotplot
dotplotPromoters <- function(result, 
                             geneId, 
                             geneName = NULL, 
                             filterInternal = TRUE,
                             col = NULL) {
  
  rdata <- rowData(result)
  gexp <- as.matrix(assays(result)$gene)
  gexp <- gexp[, result$sampleName]
  rownames(gexp) <- rdata$geneId
  gexp <- gexp[!duplicated(gexp), ]
  gexp <- data.frame(gexp)
  genelst <- rdata$geneId
  
  geneIdx <- grep(geneId, genelst)
  if (length(geneIdx) == 0) {
    stop("Gene not found")
  }
  result <- result[geneIdx,]
  internalId <- TRUE
  if (filterInternal) {
    rdata <- rowData(result)
    internalId <- rdata$internalPromoter == FALSE
  }
  assay.abs <- assays(result)$abs[internalId,]
  assay.rel <- assays(result)$rel[internalId,]
  nonzero <- apply(assay.abs, 1, function(x) !all(x==0))
  assay.abs <- assay.abs[nonzero,,drop=FALSE]
  assay.rel <- assay.rel[nonzero,,drop=FALSE]
  gexp <- gexp[grep(geneId, rownames(gexp)),,drop=FALSE]
  if (nrow(assay.abs) == 0) {
    stop("Gene has no expressed non-internal promoters")
  }
  
  condition <- result$condition
  
  plot.abs <- generateDotplot(assay.abs, condition, 
                              main = paste0("Absolute Promoter Activity ", 
                                            geneName),
                              col = col)
  plot.rel <- generateDotplot(assay.rel, condition, 
                              main = paste0("Relative Promoter Activity ", 
                                            geneName),
                              col = col)
  plot.gexp <- generateDotplot(gexp, condition, 
                               main = paste0("Gene Expression ", geneName), 
                               promoter = FALSE,
                               col = col)
  return(list(plot.abs, plot.rel, plot.gexp))
}
generateDotplot <- function(data, 
                            condition, 
                            main, 
                            promoter = TRUE,
                            col) {
  ## Reshape to long for ggplot
  colnames(data) <- condition
  data$Feature <- rownames(data)
  data <- as.data.table(data)
  data <- melt(data, ncol(data))
  if (promoter) {
    data$Feature <- paste0('prmtr.', data$Feature)
  }
  names(data) <- c("Prmtr", "Condition", "Expression")
  data$Condition <- factor(data$Condition)
  
  Prmtr <- Expression <- Condition <- NULL
  ## Promoter and gene expression plots
  if (promoter) {
    outPlot <- ggplot(data = data, aes(x = Prmtr, y = Expression, fill = Condition)) + 
      geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.8) + 
      theme_light() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
            legend.text = element_text(size = 12)) + 
      ggtitle(main)
  } else {
    outPlot <- ggplot(data = data, aes(x = Prmtr, y = Expression, fill = Condition)) + 
      geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.8) + 
      theme_light() +
      ggtitle(main)
  }
  if (!is.null(col)) {
    outPlot <- outPlot + 
      scale_fill_manual(values = col[seq_len(length(unique(condition)))])
  }
  return(outPlot)
}

#regenerate the dotplot
ASPH_promoter_activity_dot <- dotplotPromoters(result, "ENSG00000198363", filterInternal = F)
EXOC1_promoter_activity_dot <- dotplotPromoters(result, "ENSG00000090989", filterInternal = F)
AHI1_promoter_activity_dot <- dotplotPromoters(result, "ENSG00000135541", filterInternal = F)
#figure 2C-E
grid.arrange(EXOC1_promoter_activity_dot[[1]], EXOC1_promoter_activity_dot[[2]], nrow = 1, ncol = 2, widths = c(2, 2))
grid.arrange(USH2A_promoter_activity_dot[[1]], USH2A_promoter_activity_dot[[2]], nrow = 1, ncol = 2, widths = c(2, 2))
grid.arrange(AHI1_promoter_activity_dot[[1]], AHI1_promoter_activity_dot[[2]], nrow = 1, ncol = 2, widths = c(2, 2))


#Create a long dataframe summarizing cell line and promoter class
pdata1 <- data.frame(cellLine = rep(c('Amacrine', 'Astrocyte', 'Bipolar', 'Cone', 'Fibroblast', 'Horizontal', 'Muller_Glia', 'RPE', 'Rod', 'Smooth_Muscle_Cells'), each = nrow(rdata)), promoterClass = as.factor(c(rdata$Amacrine.class, rdata$Astrocyte.class, rdata$Bipolar.class, rdata$Cone.class, rdata$Fibroblast.class, rdata$Horizontal.class, rdata$Muller_Glia.class, rdata$RPE.class, rdata$Rod.class, rdata$Smooth_Muscle_Cells.class)))
ggplot(na.omit(pdata1)) +
  geom_bar(aes(x = cellLine, fill = promoterClass)) + 
  xlab('Cell Lines') + ylab('Count') +  labs(fill = 'Promoter Category') +
  ggtitle('Categorization of Promoters')

#Because many genes have many annotated promoters, we collapse promoters from the 5th position and onward into one group for simplicity
pdata2 <- as_tibble(rdata) %>%
  mutate(promoterPosition = ifelse(promoterPosition > 5, 5, promoterPosition)) %>%
  filter(Amacrine.class %in% c('Major', 'Minor'))
ggplot(pdata2) +
  geom_bar(aes(x = promoterPosition, fill = as.factor(Amacrine.class)), position = 'fill') +
  xlab(expression(Promoter ~ Position ~ "5'" %->% "3'")) + ylab('Percentage') + 
  labs(fill = 'Promoter Category') + ggtitle('Major/Minor Promoter Proportion in Amacrine') + 
  scale_y_continuous(breaks = seq(0,1, 0.25), labels = paste0(seq(0,100,25),'%')) +
  scale_x_continuous(breaks = seq(1,5), labels = c('1','2','3','4','>=5'))

pdata2 <- as_tibble(rdata) %>%
  mutate(promoterPosition = ifelse(promoterPosition > 5, 5, promoterPosition)) %>%
  filter(Astrocyte.class %in% c('Major', 'Minor'))
ggplot(pdata2) +
  geom_bar(aes(x = promoterPosition, fill = as.factor(Astrocyte.class)), position = 'fill') +
  xlab(expression(Promoter ~ Position ~ "5'" %->% "3'")) + ylab('Percentage') + 
  labs(fill = 'Promoter Category') + ggtitle('Major/Minor Promoter Proportion in Astrocyte') + 
  scale_y_continuous(breaks = seq(0,1, 0.25), labels = paste0(seq(0,100,25),'%')) +
  scale_x_continuous(breaks = seq(1,5), labels = c('1','2','3','4','>=5'))

pdata2 <- as_tibble(rdata) %>%
  mutate(promoterPosition = ifelse(promoterPosition > 5, 5, promoterPosition)) %>%
  filter(Bipolar.class %in% c('Major', 'Minor'))
ggplot(pdata2) +
  geom_bar(aes(x = promoterPosition, fill = as.factor(Bipolar.class)), position = 'fill') +
  xlab(expression(Promoter ~ Position ~ "5'" %->% "3'")) + ylab('Percentage') + 
  labs(fill = 'Promoter Category') + ggtitle('Major/Minor Promoter Proportion in Bipolar') + 
  scale_y_continuous(breaks = seq(0,1, 0.25), labels = paste0(seq(0,100,25),'%')) +
  scale_x_continuous(breaks = seq(1,5), labels = c('1','2','3','4','>=5'))

pdata2 <- as_tibble(rdata) %>%
  mutate(promoterPosition = ifelse(promoterPosition > 5, 5, promoterPosition)) %>%
  filter(Cone.class %in% c('Major', 'Minor'))
ggplot(pdata2) +
  geom_bar(aes(x = promoterPosition, fill = as.factor(Cone.class)), position = 'fill') +
  xlab(expression(Promoter ~ Position ~ "5'" %->% "3'")) + ylab('Percentage') + 
  labs(fill = 'Promoter Category') + ggtitle('Major/Minor Promoter Proportion in Cone') + 
  scale_y_continuous(breaks = seq(0,1, 0.25), labels = paste0(seq(0,100,25),'%')) +
  scale_x_continuous(breaks = seq(1,5), labels = c('1','2','3','4','>=5'))

pdata2 <- as_tibble(rdata) %>%
  mutate(promoterPosition = ifelse(promoterPosition > 5, 5, promoterPosition)) %>%
  filter(Fibroblast.class %in% c('Major', 'Minor'))
ggplot(pdata2) +
  geom_bar(aes(x = promoterPosition, fill = as.factor(Fibroblast.class)), position = 'fill') +
  xlab(expression(Promoter ~ Position ~ "5'" %->% "3'")) + ylab('Percentage') + 
  labs(fill = 'Promoter Category') + ggtitle('Major/Minor Promoter Proportion in Fibroblast') + 
  scale_y_continuous(breaks = seq(0,1, 0.25), labels = paste0(seq(0,100,25),'%')) +
  scale_x_continuous(breaks = seq(1,5), labels = c('1','2','3','4','>=5'))

pdata2 <- as_tibble(rdata) %>%
  mutate(promoterPosition = ifelse(promoterPosition > 5, 5, promoterPosition)) %>%
  filter(Horizontal.class %in% c('Major', 'Minor'))
ggplot(pdata2) +
  geom_bar(aes(x = promoterPosition, fill = as.factor(Horizontal.class)), position = 'fill') +
  xlab(expression(Promoter ~ Position ~ "5'" %->% "3'")) + ylab('Percentage') + 
  labs(fill = 'Promoter Category') + ggtitle('Major/Minor Promoter Proportion in Horizontal') + 
  scale_y_continuous(breaks = seq(0,1, 0.25), labels = paste0(seq(0,100,25),'%')) +
  scale_x_continuous(breaks = seq(1,5), labels = c('1','2','3','4','>=5'))

pdata2 <- as_tibble(rdata) %>%
  mutate(promoterPosition = ifelse(promoterPosition > 5, 5, promoterPosition)) %>%
  filter(Muller_Glia.class %in% c('Major', 'Minor'))
ggplot(pdata2) +
  geom_bar(aes(x = promoterPosition, fill = as.factor(Muller_Glia.class)), position = 'fill') +
  xlab(expression(Promoter ~ Position ~ "5'" %->% "3'")) + ylab('Percentage') + 
  labs(fill = 'Promoter Category') + ggtitle('Major/Minor Promoter Proportion in Muller_Glia') + 
  scale_y_continuous(breaks = seq(0,1, 0.25), labels = paste0(seq(0,100,25),'%')) +
  scale_x_continuous(breaks = seq(1,5), labels = c('1','2','3','4','>=5'))

pdata2 <- as_tibble(rdata) %>%
  mutate(promoterPosition = ifelse(promoterPosition > 5, 5, promoterPosition)) %>%
  filter(RPE.class %in% c('Major', 'Minor'))
ggplot(pdata2) +
  geom_bar(aes(x = promoterPosition, fill = as.factor(RPE.class)), position = 'fill') +
  xlab(expression(Promoter ~ Position ~ "5'" %->% "3'")) + ylab('Percentage') + 
  labs(fill = 'Promoter Category') + ggtitle('Major/Minor Promoter Proportion in RPE') + 
  scale_y_continuous(breaks = seq(0,1, 0.25), labels = paste0(seq(0,100,25),'%')) +
  scale_x_continuous(breaks = seq(1,5), labels = c('1','2','3','4','>=5'))

pdata2 <- as_tibble(rdata) %>%
  mutate(promoterPosition = ifelse(promoterPosition > 5, 5, promoterPosition)) %>%
  filter(Rod.class %in% c('Major', 'Minor'))
ggplot(pdata2) +
  geom_bar(aes(x = promoterPosition, fill = as.factor(Rod.class)), position = 'fill') +
  xlab(expression(Promoter ~ Position ~ "5'" %->% "3'")) + ylab('Percentage') + 
  labs(fill = 'Promoter Category') + ggtitle('Major/Minor Promoter Proportion in Rod') + 
  scale_y_continuous(breaks = seq(0,1, 0.25), labels = paste0(seq(0,100,25),'%')) +
  scale_x_continuous(breaks = seq(1,5), labels = c('1','2','3','4','>=5'))

pdata2 <- as_tibble(rdata) %>%
  mutate(promoterPosition = ifelse(promoterPosition > 5, 5, promoterPosition)) %>%
  filter(Smooth_Muscle_Cells.class %in% c('Major', 'Minor'))
ggplot(pdata2) +
  geom_bar(aes(x = promoterPosition, fill = as.factor(Smooth_Muscle_Cells.class)), position = 'fill') +
  xlab(expression(Promoter ~ Position ~ "5'" %->% "3'")) + ylab('Percentage') + 
  labs(fill = 'Promoter Category') + ggtitle('Major/Minor Promoter Proportion in Smooth_Muscle_Cells') + 
  scale_y_continuous(breaks = seq(0,1, 0.25), labels = paste0(seq(0,100,25),'%')) +
  scale_x_continuous(breaks = seq(1,5), labels = c('1','2','3','4','>=5'))

#Get active major promoters of cell types
majorPromoter <- as_tibble(rdata) %>% group_by(geneId) %>% 
  mutate(promoterCount = n()) %>% filter(Amacrine.class == 'Major') 
pdata3 <- data.frame(proActiv = majorPromoter$Amacrine.mean,
                     geneExp = majorPromoter$Amacrine.gene.mean,
                     promoterCount = majorPromoter$promoterCount)
ggplot(pdata3, aes(x = geneExp, y = proActiv)) + 
  geom_point(aes(colour = promoterCount), alpha = 0.5) +
  ggtitle('Major Promoter Activity vs. Gene Expression in Amacrine') + 
  xlab('Average Gene Expression') + ylab('Average Major Promoter Activity') +
  labs(colour = 'Number of \n Annotated Promoters') +
  geom_abline(slope = 1, intercept = 0, colour = 'red', linetype = 'dashed')

majorPromoter <- as_tibble(rdata) %>% group_by(geneId) %>% 
  mutate(promoterCount = n()) %>% filter(Astrocyte.class == 'Major') 
pdata3 <- data.frame(proActiv = majorPromoter$Astrocyte.mean,
                     geneExp = majorPromoter$Astrocyte.gene.mean,
                     promoterCount = majorPromoter$promoterCount)
ggplot(pdata3, aes(x = geneExp, y = proActiv)) + 
  geom_point(aes(colour = promoterCount), alpha = 0.5) +
  ggtitle('Major Promoter Activity vs. Gene Expression in Astrocyte') + 
  xlab('Average Gene Expression') + ylab('Average Major Promoter Activity') +
  labs(colour = 'Number of \n Annotated Promoters') +
  geom_abline(slope = 1, intercept = 0, colour = 'red', linetype = 'dashed')

majorPromoter <- as_tibble(rdata) %>% group_by(geneId) %>% 
  mutate(promoterCount = n()) %>% filter(Bipolar.class == 'Major') 
pdata3 <- data.frame(proActiv = majorPromoter$Bipolar.mean,
                     geneExp = majorPromoter$Bipolar.gene.mean,
                     promoterCount = majorPromoter$promoterCount)
ggplot(pdata3, aes(x = geneExp, y = proActiv)) + 
  geom_point(aes(colour = promoterCount), alpha = 0.5) +
  ggtitle('Major Promoter Activity vs. Gene Expression in Bipolar') + 
  xlab('Average Gene Expression') + ylab('Average Major Promoter Activity') +
  labs(colour = 'Number of \n Annotated Promoters') +
  geom_abline(slope = 1, intercept = 0, colour = 'red', linetype = 'dashed')

majorPromoter <- as_tibble(rdata) %>% group_by(geneId) %>% 
  mutate(promoterCount = n()) %>% filter(Cone.class == 'Major') 
pdata3 <- data.frame(proActiv = majorPromoter$Cone.mean,
                     geneExp = majorPromoter$Cone.gene.mean,
                     promoterCount = majorPromoter$promoterCount)
ggplot(pdata3, aes(x = geneExp, y = proActiv)) + 
  geom_point(aes(colour = promoterCount), alpha = 0.5) +
  ggtitle('Major Promoter Activity vs. Gene Expression in Cone') + 
  xlab('Average Gene Expression') + ylab('Average Major Promoter Activity') +
  labs(colour = 'Number of \n Annotated Promoters') +
  geom_abline(slope = 1, intercept = 0, colour = 'red', linetype = 'dashed')

majorPromoter <- as_tibble(rdata) %>% group_by(geneId) %>% 
  mutate(promoterCount = n()) %>% filter(Fibroblast.class == 'Major') 
pdata3 <- data.frame(proActiv = majorPromoter$Fibroblast.mean,
                     geneExp = majorPromoter$Fibroblast.gene.mean,
                     promoterCount = majorPromoter$promoterCount)
ggplot(pdata3, aes(x = geneExp, y = proActiv)) + 
  geom_point(aes(colour = promoterCount), alpha = 0.5) +
  ggtitle('Major Promoter Activity vs. Gene Expression in Fibroblast') + 
  xlab('Average Gene Expression') + ylab('Average Major Promoter Activity') +
  labs(colour = 'Number of \n Annotated Promoters') +
  geom_abline(slope = 1, intercept = 0, colour = 'red', linetype = 'dashed')

majorPromoter <- as_tibble(rdata) %>% group_by(geneId) %>% 
  mutate(promoterCount = n()) %>% filter(Horizontal.class == 'Major') 
pdata3 <- data.frame(proActiv = majorPromoter$Horizontal.mean,
                     geneExp = majorPromoter$Horizontal.gene.mean,
                     promoterCount = majorPromoter$promoterCount)
ggplot(pdata3, aes(x = geneExp, y = proActiv)) + 
  geom_point(aes(colour = promoterCount), alpha = 0.5) +
  ggtitle('Major Promoter Activity vs. Gene Expression in Horizontal') + 
  xlab('Average Gene Expression') + ylab('Average Major Promoter Activity') +
  labs(colour = 'Number of \n Annotated Promoters') +
  geom_abline(slope = 1, intercept = 0, colour = 'red', linetype = 'dashed')

majorPromoter <- as_tibble(rdata) %>% group_by(geneId) %>% 
  mutate(promoterCount = n()) %>% filter(Muller_Glia.class == 'Major') 
pdata3 <- data.frame(proActiv = majorPromoter$Muller_Glia.mean,
                     geneExp = majorPromoter$Muller_Glia.gene.mean,
                     promoterCount = majorPromoter$promoterCount)
ggplot(pdata3, aes(x = geneExp, y = proActiv)) + 
  geom_point(aes(colour = promoterCount), alpha = 0.5) +
  ggtitle('Major Promoter Activity vs. Gene Expression in Muller_Glia') + 
  xlab('Average Gene Expression') + ylab('Average Major Promoter Activity') +
  labs(colour = 'Number of \n Annotated Promoters') +
  geom_abline(slope = 1, intercept = 0, colour = 'red', linetype = 'dashed')

majorPromoter <- as_tibble(rdata) %>% group_by(geneId) %>% 
  mutate(promoterCount = n()) %>% filter(RPE.class == 'Major') 
pdata3 <- data.frame(proActiv = majorPromoter$RPE.mean,
                     geneExp = majorPromoter$RPE.gene.mean,
                     promoterCount = majorPromoter$promoterCount)
ggplot(pdata3, aes(x = geneExp, y = proActiv)) + 
  geom_point(aes(colour = promoterCount), alpha = 0.5) +
  ggtitle('Major Promoter Activity vs. Gene Expression in RPE') + 
  xlab('Average Gene Expression') + ylab('Average Major Promoter Activity') +
  labs(colour = 'Number of \n Annotated Promoters') +
  geom_abline(slope = 1, intercept = 0, colour = 'red', linetype = 'dashed')

majorPromoter <- as_tibble(rdata) %>% group_by(geneId) %>% 
  mutate(promoterCount = n()) %>% filter(Rod.class == 'Major') 
pdata3 <- data.frame(proActiv = majorPromoter$Rod.mean,
                     geneExp = majorPromoter$Rod.gene.mean,
                     promoterCount = majorPromoter$promoterCount)
ggplot(pdata3, aes(x = geneExp, y = proActiv)) + 
  geom_point(aes(colour = promoterCount), alpha = 0.5) +
  ggtitle('Major Promoter Activity vs. Gene Expression in Rod') + 
  xlab('Average Gene Expression') + ylab('Average Major Promoter Activity') +
  labs(colour = 'Number of \n Annotated Promoters') +
  geom_abline(slope = 1, intercept = 0, colour = 'red', linetype = 'dashed')

majorPromoter <- as_tibble(rdata) %>% group_by(geneId) %>% 
  mutate(promoterCount = n()) %>% filter(Smooth_Muscle_Cells.class == 'Major') 
pdata3 <- data.frame(proActiv = majorPromoter$Smooth_Muscle_Cells.mean,
                     geneExp = majorPromoter$Smooth_Muscle_Cells.gene.mean,
                     promoterCount = majorPromoter$promoterCount)
ggplot(pdata3, aes(x = geneExp, y = proActiv)) + 
  geom_point(aes(colour = promoterCount), alpha = 0.5) +
  ggtitle('Major Promoter Activity vs. Gene Expression in Smooth_Muscle_Cells') + 
  xlab('Average Gene Expression') + ylab('Average Major Promoter Activity') +
  labs(colour = 'Number of \n Annotated Promoters') +
  geom_abline(slope = 1, intercept = 0, colour = 'red', linetype = 'dashed')


#proportion of genes using single or multiple promoters 
rowDataContent = rowData(result)
#subset to different cell type
cellTypeMatrices <- list()
cellTypes <- c("Amacrine", "Astrocyte", "Bipolar", "Cone", "Fibroblast", 
               "Horizontal", "Muller_Glia", "RPE", "Rod", "Smooth_Muscle_Cells")
for (cellType in cellTypes) {
  columnName <- paste0(cellType, ".class")
  nonNAindices <- !is.na(rowDataContent[[columnName]]) & rowDataContent[[columnName]] != "NA"
  cellTypeData <- rowDataContent[nonNAindices, ]
  cellTypeClassData <- cellTypeData[, c("geneId", paste0(cellType, ".class"))]
  cellTypeMatrices[[cellType]] <- cellTypeClassData
}
for (cellType in cellTypes) {
  varName <- paste0("df_proactiv_", cellType)
  assign(varName, cellTypeMatrices[[cellType]])
}

#calculate the promoter counts
cellTypePromoterStats <- list()

for (cellType in cellTypes) {
  currentDF <- get(paste0("df_proactiv_", cellType))
  validRows <- currentDF[[paste0(cellType, ".class")]] %in% c("Major", "Minor")
  filteredDF <- currentDF[validRows, ]
  promoterInfo <- aggregate(filteredDF[[paste0(cellType, ".class")]], by = list(filteredDF$geneId), FUN = function(x) length(unique(x)))
  colnames(promoterInfo) <- c("geneId", "PromoterCount")
  singlePromoterCount <- sum(promoterInfo$PromoterCount == 1)
  multiplePromoterCount <- sum(promoterInfo$PromoterCount > 1)
  ratio <- singlePromoterCount / multiplePromoterCount
  cellTypePromoterStats[[cellType]] <- list("SinglePromoterCount" = singlePromoterCount,
                                            "MultiplePromoterCount" = multiplePromoterCount,
                                            "Ratio" = ratio)
}
cellTypePromoterStats

num_cell_types <- length(cellTypePromoterStats)
num_rows <- 2
num_cols <- ceiling(num_cell_types / num_rows)
par(mfrow=c(num_rows, num_cols), mar=c(2, 2, 2, 2))

for (cellType in names(cellTypePromoterStats)) {
  counts <- cellTypePromoterStats[[cellType]]
  total_count <- counts$SinglePromoterCount + counts$MultiplePromoterCount
  single_percentage <- counts$SinglePromoterCount / total_count * 100
  multiple_percentage <- counts$MultiplePromoterCount / total_count * 100
  labels <- sprintf("%1.1f%%", c(single_percentage, multiple_percentage))
  pie_values <- c(counts$SinglePromoterCount, counts$MultiplePromoterCount)
  pie(pie_values, labels = labels, col = c('blue', 'red'), main = paste(cellType), radius = 1)
}

par(mfrow=c(1, 1), mar=c(5, 4, 4, 2) + 0.1)

plot.new()
legend("center", legend = c("Single Promoter", "Multiple Promoters"), fill = colors, cex = 0.8, box.lty = 1)

#save the result
saveRDS(result, "proActiv_final.RDS")


##combine the plots together
#figure 2B
cell_types <- c(
  "Amacrine", "Astrocyte", "Bipolar", "Cone", "Fibroblast",
  "Horizontal", "Muller_Glia", "RPE", "Rod", "Smooth_Muscle_Cells"
)

theme_big <- theme(
  text = element_text(size = 41), 
  axis.text = element_text(size = 37),
  axis.title = element_text(size = 41),
  plot.title = element_text(size = 41, face = "bold"),
  legend.title = element_text(size = 41),
  legend.text  = element_text(size = 37),
  strip.text = element_text(size = 41)
)

make_promoter_plot <- function(ct, show_x = FALSE, show_y = FALSE) {
  class_col <- paste0(ct, ".class")
  
  df <- as_tibble(rdata) %>%
    mutate(promoterPosition = ifelse(promoterPosition > 5, 5, promoterPosition)) %>%
    filter(.data[[class_col]] %in% c("Major", "Minor"))
  
  df$PromoterClass <- df[[class_col]]
  ct_label <- gsub("_", " ", ct)
  
  x_lab <- if (show_x) expression(Promoter ~ Position ~ "5'" %->% "3'") else NULL
  y_lab <- if (show_y) "Percentage" else NULL
  
  ggplot(df) +
    geom_bar(aes(x = promoterPosition, fill = PromoterClass), position = "fill") +
    scale_y_continuous(
      breaks = seq(0, 1, 0.25),
      labels = paste0(seq(0, 100, 25), "%"),
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_x_continuous(
      breaks = 1:5,
      labels = c("1","2","3","4",">=5")
    ) +
    labs(
      x = x_lab,  
      y = y_lab,
      fill = "Promoter Category",
      title = ct_label     
    ) +
    theme_big +
    theme(
      plot.title = element_text(size = 32, hjust = 0.5)
    )
}

n_ct <- length(cell_types)
ncol_layout <- 5
nrow_layout <- ceiling(n_ct / ncol_layout)

show_x_flags <- logical(n_ct)

for (row in 1:nrow_layout) {
  idx_in_row <- which(ceiling(seq_len(n_ct) / ncol_layout) == row)
  if (length(idx_in_row) > 0) {
    mid_idx <- idx_in_row[ceiling(length(idx_in_row) / 2)]
    show_x_flags[mid_idx] <- TRUE
  }
}

promoter_plots <- lapply(seq_along(cell_types), function(i) {
  ct <- cell_types[i]
  
  row_i <- ceiling(i / ncol_layout)
  col_i <- ((i - 1) %% ncol_layout) + 1
  
  show_x <- show_x_flags[i]
  show_y <- (col_i == 1) 
  
  make_promoter_plot(ct, show_x = show_x, show_y = show_y)
})

promoter_combined <- wrap_plots(
  promoter_plots,
  ncol = ncol_layout,
  guides = "collect"
) +
  plot_annotation(
    title = "Major and Minor Promoter Distribution Across Cell Types",
    theme = theme(
      plot.title = element_text(size = 60, face = "bold", hjust = 0.5)
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.justification = "center"
  )

promoter_combined

ggsave(
  "promoter_category_figure2B.tiff",
  plot  = promoter_combined,
  device = "tiff",
  type   = "cairo", 
  width = 40, height = 20, units = "in",
  dpi   = 600,
  compression = "lzw", bg = "WHITE"
)

#figure supplementary 6
theme_big_scatter <- theme(
  text = element_text(size = 26), 
  axis.text = element_text(size = 21),
  axis.title = element_text(size = 26),
  plot.title = element_text(size = 18, face = "bold"),
  legend.title = element_text(size = 26),
  legend.text  = element_text(size = 22),
  strip.text = element_text(size = 26), 
  legend.key.height = unit(1.0, "cm"),
  legend.key.width  = unit(0.4, "cm"),
)
make_major_scatter <- function(ct) {
  class_col      <- paste0(ct, ".class")
  proActiv_col   <- paste0(ct, ".mean")
  geneExp_col    <- paste0(ct, ".gene.mean")
  
  df <- as_tibble(rdata) %>%
    group_by(geneId) %>%
    mutate(promoterCount = n()) %>%
    filter(.data[[class_col]] == "Major")
  df$proActiv <- df[[proActiv_col]]
  df$geneExp  <- df[[geneExp_col]]
  
  ct_label <- gsub("_", " ", ct)
  
  ggplot(df, aes(x = geneExp, y = proActiv)) +
    geom_point(aes(colour = promoterCount), alpha = 0.5) +
    ggtitle(paste("Major Promoter Activity vs. Gene Expression in", ct_label)) +
    xlab("Average Gene Expression") +
    ylab("Average Major Promoter Activity") +
    labs(colour = "Number of \n Annotated Promoters") +
    geom_abline(slope = 1, intercept = 0,
                colour = "red", linetype = "dashed") + 
    theme_big_scatter
}
major_scatter_plots <- lapply(cell_types, make_major_scatter)

major_scatter_combined <- wrap_plots(
  major_scatter_plots,
  ncol = 5
) &
  theme(
    legend.position = "right",
    legend.justification = "center"
  )

major_scatter_combined
ggsave(
  "promoter_expression_supplementary6.tiff",
  plot  = major_scatter_combined,
  device = "tiff",
  type   = "cairo", 
  width = 45, height = 15, units = "in",
  dpi   = 600,
  compression = "lzw", bg = "WHITE"
)

#figure2A
library(plotrix)
library(dplyr)
library(tidyr)
library(ggplot2)

df_pie <- imap_dfr(cellTypePromoterStats, ~{
  tibble(
    cellType = .y,
    SinglePromoterCount   = .x$SinglePromoterCount,
    MultiplePromoterCount = .x$MultiplePromoterCount
  )
}) %>%
  pivot_longer(
    cols = c(SinglePromoterCount, MultiplePromoterCount),
    names_to  = "PromoterCategory",
    values_to = "Count"
  ) %>%
  mutate(
    PromoterCategory = case_when(
      PromoterCategory == "SinglePromoterCount"   ~ "Single promoter",
      PromoterCategory == "MultiplePromoterCount" ~ "Multiple promoters",
      TRUE ~ PromoterCategory
    )
  ) %>%
  group_by(cellType) %>%
  mutate(
    Percentage = Count / sum(Count) * 100,
    label = sprintf("%.1f%%", Percentage)
  ) %>%
  ungroup()

pie_cols <- c(
  "Single promoter"    = "#4F81BD",
  "Multiple promoters" = "#F39C12"
)

p_pie <- ggplot(df_pie, aes(x = 1, y = Count, fill = PromoterCategory)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ cellType, nrow = 2, scales = "free_y") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 5
  ) +
  scale_fill_manual(values = pie_cols, name = "Promoter category") +
  theme_void(base_size = 20) +
  theme(
    strip.text      = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 12),
    plot.margin     = margin(15, 15, 15, 15)
  )

p_pie
ggsave(
  "PromoterPie_AllCellTypes.tiff",
  plot  = p_pie,
  device = "tiff",
  type   = "cairo", 
  width = 10, height = 5, units = "in",
  dpi   = 1200,
  compression = "lzw", bg = "WHITE"
)

#figure2CDE
library(gridExtra)
library(grid)

FS <- list(
  base         = 28,
  axis_title   = 30,
  title        = 32,
  bigtitle     = 36,
  legend_title = 30,
  legend_point = 15
)

decorate_plot <- function(p, gene_name) {
  p +
    labs(
      title = gene_name,
      x = "Promoter ID"
    ) +
    theme(
      plot.title = element_text(size = FS$title, face = "bold", hjust = 0.5),
      
      axis.title.x = element_text(size = FS$axis_title, face = "bold",
                                  margin = margin(t = 70)),
      axis.title.y = element_text(size = FS$axis_title, face = "bold",
                                  margin = margin(r = 6)),
      
      axis.text.x  = element_text(size = FS$base, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = FS$base),
      
      legend.title = element_text(size = FS$legend_title, face = "bold"),
      legend.text  = element_text(size = FS$base)
    ) +
    guides(
      colour = guide_legend(override.aes = list(size = FS$legend_point)),
      fill   = guide_legend(override.aes = list(size = FS$legend_point))
    )
}

p1 <- decorate_plot(AHI1_promoter_activity_dot[[2]],  "AHI1")
p2 <- decorate_plot(EXOC1_promoter_activity_dot[[2]], "EXOC1")
p3 <- decorate_plot(USH2A_promoter_activity_dot[[2]], "USH2A")

row_plots <- arrangeGrob(p1, p2, p3, nrow = 1)
final_plot <- arrangeGrob(
  row_plots,
  top = textGrob(
    "Relative promoter activity",
    gp = gpar(fontsize = FS$bigtitle, fontface = "bold")
  )
)
grid.newpage()
grid.draw(final_plot)

tiff(
  filename   = "relative_promoter_activity.tiff",
  width      = 35,
  height     = 20,
  units      = "in",
  res        = 600,
  compression = "lzw"
)
grid.draw(final_plot)
dev.off()
#figure 2C-E
grid.arrange(EXOC1_promoter_activity_dot[[1]], EXOC1_promoter_activity_dot[[2]], nrow = 1, ncol = 2, widths = c(2, 2))
grid.arrange(USH2A_promoter_activity_dot[[1]], USH2A_promoter_activity_dot[[2]], nrow = 1, ncol = 2, widths = c(2, 2))
grid.arrange(AHI1_promoter_activity_dot[[1]], AHI1_promoter_activity_dot[[2]], nrow = 1, ncol = 2, widths = c(2, 2))


#figure 3A-F
#load the general Rdata
#read the gene list generated from retnet
retnet_list <- read.table("RetNet_category.txt", header = TRUE, sep = "\t")
#remove the duplicated gene
retnet_disease_unique <- lapply(retnet_list, unique)
#extract the unique gene list
retnet_genes <- as.vector(unlist(retnet_list))
retnet_gene_unique <- unique(retnet_genes)
retnet_gene_unique <- data.frame(retnet_gene_unique)
#find the corresponding pbid
retnet_pbid_FL2 <- transcript_status_FL2 %>%
  filter(gene %in% retnet_gene_unique$retnet_gene_unique) %>%
  select(gene, pbid) %>%
  distinct()

head(retnet_pbid_FL2)

mapped_retnet_FL2 <- unique(retnet_pbid_FL2$gene)
num_mapped_retnet_FL2 <- length(mapped_retnet_FL2)
print(num_mapped_retnet_FL2)

#set the barcode-celltype information first
lr_metadata <- retina_combined@meta.data %>%
  dplyr::select(celltype = cell.type.integrated) %>%
  rownames_to_column(var = "barcode")
cell_annotation_data <- as.data.frame(lr_metadata$celltype)
colnames(cell_annotation_data) <- "celltype"
rownames(cell_annotation_data) <- lr_metadata$barcode
lr_celltype <- retina_combined@meta.data$cell.type.integrated
lr_celltype <- factor(lr_celltype, levels = unique(lr_celltype))
barcode_data <- data.frame(
  barcode = colnames(retina_combined),
  cell.type.integrated = as.character(retina_combined@meta.data$cell.type.integrated)
)
barcode_data <- barcode_data[order(barcode_data$cell.type.integrated), ]
sorted_barcodes <- barcode_data$barcode


#extract the count matrix from seurat object for heatmap
iso_counts <- retina_combined@assays$isoform@counts
#keep only PBID in rownames
rownames(iso_counts) <- sub(":.*$", "", rownames(iso_counts))
selected_pbid_FL2 <- as.vector(retnet_pbid_FL2$pbid)
valid_pbid_FL2 <- selected_pbid_FL2[selected_pbid_FL2 %in% rownames(iso_counts)]
sub_iso_matrix_FL2 <- iso_counts[valid_pbid_FL2, sorted_barcodes]

rna_counts <- retina_combined@assays$RNA@counts

selected_gene_FL2 <- unique(retnet_pbid_FL2$gene)
valid_gene_FL2 <- selected_gene_FL2[selected_gene_FL2 %in% rownames(rna_counts)]
sub_rna_matrix_FL2 <- rna_counts[valid_gene_FL2, sorted_barcodes]
regular_rna_matrix_FL2 <- as.matrix(sub_rna_matrix_FL2)

#normal heatmap cannot fit in near 300 genes, thus subset them according to diseases
CORD_list <- subset(retnet_list, select = "CORD")
colnames(CORD_list) <- "gene"
unique_CORD <- distinct(CORD_list, gene, .keep_all = TRUE)
CSNB_list <- subset(retnet_list, select = "CSNB")
colnames(CSNB_list) <- "gene"
unique_CSNB <- distinct(CSNB_list, gene, .keep_all = TRUE)
LCA_list <- subset(retnet_list, select = "LCA")
colnames(LCA_list) <- "gene"
unique_LCA <- distinct(LCA_list, gene, .keep_all = TRUE)
MD_list <- subset(retnet_list, select = "MD")
colnames(MD_list) <- "gene"
unique_MD <- distinct(MD_list, gene, .keep_all = TRUE)
RP_list <- subset(retnet_list, select = "RP")
colnames(RP_list) <- "gene"
unique_RP <- distinct(RP_list, gene, .keep_all = TRUE)
US_list <- subset(retnet_list, select = "Usher_syndrome")
colnames(US_list) <- "gene"
unique_US <- distinct(US_list, gene, .keep_all = TRUE)
AMD_list <- subset(retnet_list, select = "AMD")
colnames(AMD_list) <- "gene"
unique_AMD <- distinct(AMD_list, gene, .keep_all = TRUE)


#example for RP: Retinitis pigmentosa
#extract RP-related isoforms (FL2 filtered)
RP_pbid_FL2 <- transcript_status_FL2 %>%
  filter(gene %in% unique_RP$gene) %>%      # keep RP genes only
  mutate(status = case_when(
    category %in% c("full-splice_match", "incomplete-splice_match") ~ "known",
    category %in% c("novel_in_catalog", "novel_not_in_catalog") ~ "novel",
    TRUE ~ "other"
  )) %>%
  select(gene, pbid, status) %>%
  distinct()

#keep only known/novel isoforms
RP_pbid_FL2 <- RP_pbid_FL2 %>%
  filter(status %in% c("known", "novel"))

#create transcript_status label (e.g., RHO_known / RHO_novel)
RP_pbid_FL2$transcript_status <- paste(RP_pbid_FL2$gene,
                                       RP_pbid_FL2$status,
                                       sep = "_")
RP_pbid_FL2$gene <- NULL
RP_pbid_FL2$status <- NULL

head(RP_pbid_FL2)
length(unique(RP_pbid_FL2$transcript_status))


#extract count matrix for RP FL2 isoforms
RP_selected_pbid_FL2 <- as.vector(RP_pbid_FL2$pbid)
RP_valid_pbid_FL2 <- intersect(RP_selected_pbid_FL2, rownames(iso_counts))

sub_RP_iso_matrix_FL2 <- iso_counts[RP_valid_pbid_FL2, , drop = FALSE]

#add pbid column and merge with transcript_status
sub_RP_iso_df_FL2 <- as.data.frame(sub_RP_iso_matrix_FL2)
sub_RP_iso_df_FL2$pbid <- rownames(sub_RP_iso_df_FL2)

sub_RP_iso_df_FL2 <- left_join(sub_RP_iso_df_FL2,
                               RP_pbid_FL2,
                               by = "pbid")

unique_RP_iso_count_FL2 <- length(unique(sub_RP_iso_df_FL2$transcript_status))
print(paste("Number of unique transcript_status (RP, FL2):",
            unique_RP_iso_count_FL2))

#sum counts for isoforms with same transcript_status
RP_grouped_df_FL2 <- sub_RP_iso_df_FL2 %>%
  group_by(transcript_status) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")

RP_final_matrix_FL2 <- as.matrix(RP_grouped_df_FL2[, -1])
rownames(RP_final_matrix_FL2) <- RP_grouped_df_FL2$transcript_status


#sort matrix columns by celltype order
common_barcodes_RP <- intersect(sorted_barcodes,
                                colnames(RP_final_matrix_FL2))
RP_final_sorted_matrix_FL2 <- RP_final_matrix_FL2[, common_barcodes_RP, drop = FALSE]

#create annotation vector: known vs novel
RP_row_names_FL2 <- rownames(RP_final_sorted_matrix_FL2)
RP_cluster_col_FL2 <- ifelse(grepl("novel", RP_row_names_FL2), "novel",
                             ifelse(grepl("known", RP_row_names_FL2), "known", NA))
RP_TS_col_FL2 <- data.frame(transcript_status = RP_cluster_col_FL2,
                            row.names = RP_row_names_FL2)
head(RP_TS_col_FL2)


#build Seurat object for dotplot
RP_seurat_FL2 <- CreateSeuratObject(counts = RP_final_sorted_matrix_FL2)

#add celltype information
RP_seurat_FL2[["celltype"]] <- barcode_data$cell.type.integrated[
  match(colnames(RP_seurat_FL2), barcode_data$barcode)
]

#clean rownames for Seurat
RP_new_row_names_FL2 <- gsub(pattern = "_",
                             replacement = "-",
                             x = rownames(RP_TS_col_FL2))

rownames(RP_TS_col_FL2) <- RP_new_row_names_FL2

RP_TS_col_FL2_df <- data.frame(
  gene = rownames(RP_TS_col_FL2),
  transcript_status = RP_TS_col_FL2$transcript_status,
  row.names = NULL
)

#order RP transcript_status (novel first)
RP_sorted_genes_df_FL2 <- RP_TS_col_FL2_df[
  order(RP_TS_col_FL2_df$transcript_status, decreasing = TRUE), ]


#normalize and scale
RP_seurat_FL2 <- NormalizeData(RP_seurat_FL2)
RP_seurat_FL2 <- ScaleData(RP_seurat_FL2,
                           features = rownames(RP_seurat_FL2))

Idents(RP_seurat_FL2) <- "celltype"
clusters_FL2 <- SplitObject(RP_seurat_FL2,
                            split.by = "celltype")


#identify novel transcript_status expressed in >1% cells per cluster
RP_novel_gene_FL2 <- subset(RP_sorted_genes_df_FL2,
                            transcript_status == "novel")

expressed_RP_novel_list_FL2 <- list()
for (celltype in names(clusters_FL2)) {
  cluster <- clusters_FL2[[celltype]]
  expressed_cells_per_RP_novel <- rowSums(
    cluster@assays$RNA@counts[RP_novel_gene_FL2$gene, ] > 0
  )
  expressed_RP_novel <- names(
    expressed_cells_per_RP_novel[
      expressed_cells_per_RP_novel > (ncol(cluster) * 0.01)
    ])
  expressed_RP_novel_list_FL2[[celltype]] <- expressed_RP_novel
}
all_expressed_RP_novel_list_FL2 <- unique(unlist(expressed_RP_novel_list_FL2))
all_expressed_RP_novel_list_FL2


#identify known transcript_status expressed in >1% cells
RP_known_gene_FL2 <- subset(RP_sorted_genes_df_FL2,
                            transcript_status == "known")

expressed_RP_known_list_FL2 <- list()
for (celltype in names(clusters_FL2)) {
  cluster <- clusters_FL2[[celltype]]
  expressed_cells_per_RP_known <- rowSums(
    cluster@assays$RNA@counts[RP_known_gene_FL2$gene, ] > 0
  )
  expressed_RP_known <- names(
    expressed_cells_per_RP_known[
      expressed_cells_per_RP_known > (ncol(cluster) * 0.01)
    ])
  expressed_RP_known_list_FL2[[celltype]] <- expressed_RP_known
}
all_expressed_RP_known_list_FL2 <- unique(unlist(expressed_RP_known_list_FL2))
all_expressed_RP_known_list_FL2


#dotPlot for novel isoforms
RP_novel_dotplot_FL2 <- DotPlot(RP_seurat_FL2,
                                features = all_expressed_RP_novel_list_FL2,
                                dot.scale = 10,
                                scale.by = "size",
                                col.min = 0,
                                scale.min = 1) +
  coord_flip() +
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gradientn(colours = c("#330066", "#336699",
                                    "#66CC66", "#FFCC33"))

#dotPlot for known isoforms
RP_known_dotplot_FL2 <- DotPlot(RP_seurat_FL2,
                                features = all_expressed_RP_known_list_FL2,
                                dot.scale = 10,
                                scale.by = "size",
                                col.min = 0,
                                scale.min = 1) +
  coord_flip() +
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gradientn(colours = c("#617381", "#76AB9B",
                                    "#CE5A5B", "#EA7C5B"))

grid.arrange(RP_novel_dotplot_FL2,
             RP_known_dotplot_FL2,
             nrow = 1)


#dotPlot for all novel isoforms
all.RP_novel_dotplot_FL2 <- DotPlot(RP_seurat_FL2,
                                    features = RP_novel_gene_FL2$gene,
                                    dot.scale = 10,
                                    scale.by = "size",
                                    col.min = 0,
                                    scale.min = 1) +
  coord_flip() +
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gradientn(colours = c("#330066", "#336699",
                                    "#66CC66", "#FFCC33"))

#dotPlot for all known isoforms
all.RP_known_dotplot_FL2 <- DotPlot(RP_seurat_FL2,
                                    features = RP_known_gene_FL2$gene,
                                    dot.scale = 10,
                                    scale.by = "size",
                                    col.min = 0,
                                    scale.min = 1) +
  coord_flip() +
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gradientn(colours = c("#617381", "#76AB9B",
                                    "#CE5A5B", "#EA7C5B"))

grid.arrange(all.RP_novel_dotplot_FL2,
             all.RP_known_dotplot_FL2,
             nrow = 1)

#code for single disease runs okay
#define a function to comply all steps included
library(dplyr)
library(Seurat)
library(ggplot2)
library(gridExtra)

make_disease_isoform_plots_FL2 <- function(
    disease_col,
    retnet_list,
    transcript_status_FL2,
    iso_counts,
    barcode_data,
    sorted_barcodes,
    min_frac = 0.01
) {
  message("=== Processing disease: ", disease_col, " ===")
  
  #prepare disease gene list from RetNet
  if (!disease_col %in% colnames(retnet_list)) {
    stop("Column ", disease_col, " not found in retnet_list")
  }
  
  disease_list <- retnet_list[, disease_col, drop = FALSE]
  colnames(disease_list) <- "gene"
  disease_genes <- distinct(disease_list, gene, .keep_all = TRUE)
  
  #build disease-specific isoform table (FL2 filtered)
  disease_pbid_FL2 <- transcript_status_FL2 %>%
    filter(gene %in% disease_genes$gene) %>% 
    mutate(status = case_when(
      category %in% c("full-splice_match", "incomplete-splice_match") ~ "known",
      category %in% c("novel_in_catalog", "novel_not_in_catalog") ~ "novel",
      TRUE ~ "other"
    )) %>%
    select(gene, pbid, status) %>%
    distinct() %>%
    filter(status %in% c("known", "novel"))
  
  if (nrow(disease_pbid_FL2) == 0) {
    warning("No FL2 isoforms found for disease ", disease_col)
    return(NULL)
  }
  
  disease_pbid_FL2$transcript_status <- paste(
    disease_pbid_FL2$gene,
    disease_pbid_FL2$status,
    sep = "_"
  )
  disease_pbid_FL2$gene <- NULL
  disease_pbid_FL2$status <- NULL
  
  message("  - unique transcript_status: ",
          length(unique(disease_pbid_FL2$transcript_status)))
  
  #extract isoform count matrix (PBID-based)
  iso_counts2 <- iso_counts
  #make sure rownames are PBID only
  rownames(iso_counts2) <- sub(":.*$", "", rownames(iso_counts2))
  
  selected_pbid <- as.vector(disease_pbid_FL2$pbid)
  valid_pbid <- intersect(selected_pbid, rownames(iso_counts2))
  
  if (length(valid_pbid) == 0) {
    warning("No matching PBID in iso_counts for disease ", disease_col)
    return(NULL)
  }
  
  sub_iso_matrix <- iso_counts2[valid_pbid, , drop = FALSE]
  
  sub_iso_df <- as.data.frame(sub_iso_matrix)
  sub_iso_df$pbid <- rownames(sub_iso_df)
  
  sub_iso_df <- left_join(sub_iso_df, disease_pbid_FL2, by = "pbid")
  
  unique_ts <- length(unique(sub_iso_df$transcript_status))
  message("  - transcript_status after join: ", unique_ts)
  
  #sum isoforms with same transcript_status
  grouped_df <- sub_iso_df %>%
    group_by(transcript_status) %>%
    summarise(across(where(is.numeric), sum), .groups = "drop")
  
  final_matrix <- as.matrix(grouped_df[, -1])
  rownames(final_matrix) <- grouped_df$transcript_status
  
  #order columns by cell type
  common_barcodes <- intersect(sorted_barcodes, colnames(final_matrix))
  final_sorted_matrix <- final_matrix[, common_barcodes, drop = FALSE]
  
  row_names <- rownames(final_sorted_matrix)
  cluster_col <- ifelse(grepl("novel", row_names), "novel",
                        ifelse(grepl("known", row_names), "known", NA))
  TS_col <- data.frame(transcript_status = cluster_col,
                       row.names = row_names)
  
  #build Seurat object
  seu <- CreateSeuratObject(counts = final_sorted_matrix)
  
  seu[["celltype"]] <- barcode_data$cell.type.integrated[
    match(colnames(seu), barcode_data$barcode)
  ]
  
  #clean rownames for Seurat
  new_row_names <- gsub(pattern = "_", replacement = "-", x = rownames(TS_col))
  rownames(TS_col) <- new_row_names
  
  TS_col_df <- data.frame(
    gene = rownames(TS_col),
    transcript_status = TS_col$transcript_status,
    row.names = NULL
  )
  
  sorted_genes_df <- TS_col_df[
    order(TS_col_df$transcript_status, decreasing = TRUE), ]
  
  seu <- NormalizeData(seu)
  seu <- ScaleData(seu, features = rownames(seu))
  
  Idents(seu) <- "celltype"
  clusters <- SplitObject(seu, split.by = "celltype")
  
  #find novel / known transcript_status expressed in > min_frac
  novel_gene_df <- subset(sorted_genes_df, transcript_status == "novel")
  known_gene_df <- subset(sorted_genes_df, transcript_status == "known")
  
  expr_mat <- seu@assays$RNA@counts
  
  #helper to compute "expressed in > min_frac cells per cluster"
  get_expressed_features <- function(gene_df, clusters, min_frac) {
    expressed_list <- list()
    for (celltype in names(clusters)) {
      cluster <- clusters[[celltype]]
      genes <- intersect(gene_df$gene, rownames(cluster@assays$RNA@counts))
      if (length(genes) == 0) next
      mat <- cluster@assays$RNA@counts[genes, , drop = FALSE]
      expressed_cells_per_gene <- rowSums(mat > 0)
      keep <- expressed_cells_per_gene > (ncol(cluster) * min_frac)
      expressed_list[[celltype]] <- names(expressed_cells_per_gene[keep])
    }
    unique(unlist(expressed_list))
  }
  
  all_expressed_novel <- get_expressed_features(novel_gene_df, clusters, min_frac)
  all_expressed_known <- get_expressed_features(known_gene_df, clusters, min_frac)
  
  #build dotplots
  novel_dot <- DotPlot(
    seu,
    features = all_expressed_novel,
    dot.scale = 10,
    scale.by = "size",
    col.min = 0,
    scale.min = 1
  ) +
    coord_flip() +
    ggtitle(paste0(disease_col, " novel (filtered)")) +
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_color_gradientn(colours = c("#330066", "#336699",
                                      "#66CC66", "#FFCC33"))
  
  known_dot <- DotPlot(
    seu,
    features = all_expressed_known,
    dot.scale = 10,
    scale.by = "size",
    col.min = 0,
    scale.min = 1
  ) +
    coord_flip() +
    ggtitle(paste0(disease_col, " known (filtered)")) +
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_color_gradientn(colours = c("#617381", "#76AB9B",
                                      "#CE5A5B", "#EA7C5B"))
  
  all_novel_dot <- DotPlot(
    seu,
    features = novel_gene_df$gene,
    dot.scale = 10,
    scale.by = "size",
    col.min = 0,
    scale.min = 1
  ) +
    coord_flip() +
    ggtitle(paste0(disease_col, " all novel")) +
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_color_gradientn(colours = c("#330066", "#336699",
                                      "#66CC66", "#FFCC33"))
  
  all_known_dot <- DotPlot(
    seu,
    features = known_gene_df$gene,
    dot.scale = 10,
    scale.by = "size",
    col.min = 0,
    scale.min = 1
  ) +
    coord_flip() +
    ggtitle(paste0(disease_col, " all known")) +
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_color_gradientn(colours = c("#617381", "#76AB9B",
                                      "#CE5A5B", "#EA7C5B"))
  
  #return infor
  invisible(list(
    disease = disease_col,
    seurat = seu,
    plots = list(
      novel_filtered = novel_dot,
      known_filtered = known_dot,
      novel_all = all_novel_dot,
      known_all = all_known_dot
    ),
    matrices = list(
      final_sorted_matrix = final_sorted_matrix
    ),
    data_frames = list(
      disease_pbid_FL2 = disease_pbid_FL2,
      sub_iso_df = sub_iso_df,
      grouped_df = grouped_df,
      TS_col = TS_col,
      sorted_genes_df = sorted_genes_df,
      expressed_novel_features = all_expressed_novel,
      expressed_known_features = all_expressed_known
    )
  ))
}

#run the function using each disease of interest
res_RP <- make_disease_isoform_plots_FL2(
  "RP", retnet_list,
  transcript_status_FL2,
  iso_counts,
  barcode_data,
  sorted_barcodes
)
res_CORD <- make_disease_isoform_plots_FL2(
  "CORD", retnet_list,
  transcript_status_FL2,
  iso_counts,
  barcode_data,
  sorted_barcodes
)
res_MD <- make_disease_isoform_plots_FL2(
  "MD", retnet_list,
  transcript_status_FL2,
  iso_counts,
  barcode_data,
  sorted_barcodes
)
res_CSNB <- make_disease_isoform_plots_FL2(
  "CSNB", retnet_list,
  transcript_status_FL2,
  iso_counts,
  barcode_data,
  sorted_barcodes
)
res_LCA <- make_disease_isoform_plots_FL2(
  "LCA", retnet_list,
  transcript_status_FL2,
  iso_counts,
  barcode_data,
  sorted_barcodes
)
res_US <- make_disease_isoform_plots_FL2(
  "Usher_syndrome", retnet_list,
  transcript_status_FL2,
  iso_counts,
  barcode_data,
  sorted_barcodes
)

#save the plot
RP_novel_all_plot <- res_RP$plots$novel_all
CORD_novel_all_plot <- res_CORD$plots$novel_all
MD_novel_all_plot   <- res_MD$plots$novel_all
CSNB_novel_all_plot <- res_CSNB$plots$novel_all
LCA_novel_all_plot  <- res_LCA$plots$novel_all
US_novel_all_plot   <- res_US$plots$novel_all

dot_theme_big_axis_small_legend <- theme(
  axis.title.x = element_text(size = 22, face = "bold"),
  axis.title.y = element_text(size = 22, face = "bold"),
  axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
  axis.text.y  = element_text(size = 16),
  legend.position = "right",
  legend.justification = "center",
  legend.title = element_text(size = 14),
  legend.text  = element_text(size = 14),
  plot.title = element_text(size = 26, hjust = 0.5, face = "bold")
)

RP_novel_all_plot  <- RP_novel_all_plot  + dot_theme_big_axis_small_legend
CORD_novel_all_plot <- CORD_novel_all_plot + dot_theme_big_axis_small_legend
MD_novel_all_plot   <- MD_novel_all_plot   + dot_theme_big_axis_small_legend
CSNB_novel_all_plot <- CSNB_novel_all_plot + dot_theme_big_axis_small_legend
LCA_novel_all_plot  <- LCA_novel_all_plot  + dot_theme_big_axis_small_legend
US_novel_all_plot   <- US_novel_all_plot   + dot_theme_big_axis_small_legend

p_figure3 <- (CORD_novel_all_plot / MD_novel_all_plot) |
  (CSNB_novel_all_plot / LCA_novel_all_plot / US_novel_all_plot) |
  RP_novel_all_plot

ggsave(
  filename = "Figure3new.Dotplot.IRDgenes.tiff",
  plot     = p_figure3,
  width    = 22,
  height   = 24,
  dpi      = 600,
  device   = "tiff",
  compression = "lzw",
  bg       = "WHITE"
)


tiff(
  filename = "Figure3new.Dotplot.IRDgenes.tiff",
  width = 15000, 
  height = 16000,
  res = 600,
  compression = "lzw" 
)

(
  res_CORD$plots$novel_all / res_MD$plots$novel_all
) |
  (
    res_CSNB$plots$novel_all /
      res_LCA$plots$novel_all /
      res_US$plots$novel_all
  ) |
  res_RP$plots$novel_all

dev.off()


#ggtranscript to plot the structure of transcripts
#read the annotation file
#read the file for filtering
pb_anno <- rtracklayer::import("RetinaMerge.sorted.filtered_lite_FL2.gff")
class(pb_anno)
pb_anno <- pb_anno %>% dplyr::as_tibble()
class(pb_anno)
junction_anno <- read.table("RetinaMerge_classification.filtered_lite_FL2_junctions.txt", header = TRUE)

reference_anno <- rtracklayer::import("GRCh38p14/GRCh38.p14.v44/gencode.v44.chr_patch_hapl_scaff.annotation.gff3")
class(reference_anno)
reference_anno <- reference_anno %>% dplyr::as_tibble()
class(reference_anno)

#read the file for filtering
pb_infor <- read.csv("RetinaMerge_classification.filtered_lite_FL2_classification.txt", header = TRUE, sep = "\t")
#previously only select for PBids within cage peak, then I found can filter out those NA polyA motifs, so change the code but retain the name as it is used in the following codes to reduce work
#after filter with both within cage peak and with non-NA polyA motif, the real name for this object should be pbid_cage_polyA or pbid_FL
#pbid_cage <- pb_infor %>% filter(within_cage_peak == TRUE)
#FL
pbid_fl <- pb_infor %>%
  filter(within_cage_peak == TRUE)
pbid_cage <- pb_infor %>% filter(within_cage_peak == TRUE)
#read the rds file
retina_combined <- readRDS("iso_sr_merged_retina_FL2.RDS")
#first read the RDS file
#then extract the count matrix and cell type information
#then count the unique PBids for each cell type
#finally will introduce back to the ggtranscript for plotting
seurat_list <- SplitObject(retina_combined, split.by = "cell.type.integrated")
iso_counts_list <- lapply(seurat_list, function(x) {
  GetAssayData(x, assay = "isoform", slot = "counts")
})
iso_counts_list <- lapply(iso_counts_list, function(mat) {
  rownames(mat) <- sub(":.*$", "", rownames(mat))
  return(mat)
})
expressed_genes_by_type <- lapply(iso_counts_list, function(counts) {
  genes <- rownames(counts)[rowSums(counts > 0) > 0]
  return(unique(genes))
})
transcript_cell_type_df <- do.call(rbind, lapply(names(expressed_genes_by_type), function(cell_type) {
  data.frame(
    transcript_id = expressed_genes_by_type[[cell_type]],
    cell_type = cell_type,
    stringsAsFactors = FALSE
  )
}))

#ggtranscript plot of raw exon map
#PBid for ABCA4
ABCA4_filtered_cage <- pbid_fl %>%
  filter(associated_gene == "ABCA4")

ABCA4_PBid_of_interest <- ABCA4_filtered_cage$isoform

ABCA4_annotation_from_gtf <- pb_anno %>%
  filter(
    !is.na(transcript_id),
    transcript_id %in% ABCA4_PBid_of_interest
  ) %>%
  select(
    seqnames,
    start,
    end,
    strand,
    type,
    gene_id,
    transcript_id
  )

ABCA4_annotation_from_gtf$gene_id <- "ABCA4"

#Reference GENCODE annotation for ABCA4
gene_of_interest <- "ABCA4"

ABCA4_reference_annotation_from_gtf <- reference_anno %>%
  filter(
    !is.na(gene_name),
    gene_name %in% gene_of_interest
  ) %>%
  select(
    seqnames,
    start,
    end,
    strand,
    type,
    gene_name,
    transcript_name
  )

#color for cell types
colors <- brewer.pal(12, "Paired")
color_map <- setNames(
  colors[1:11],
  c("Amacrine", "Bipolar", "Cone", "Muller_Glia", "Rod",
    "RPE", "Smooth Muscle Cells", "Horizontal",
    "Astrocyte", "Fibroblast", "zReference")
)

cell_type_label_map <- c(
  "zReference"          = "Reference",
  "Amacrine"            = "Amacrine",
  "Bipolar"             = "Bipolar",
  "Cone"                = "Cone",
  "Muller_Glia"         = "Muller glia",
  "Rod"                 = "Rod",
  "RPE"                 = "RPE",
  "Smooth Muscle Cells" = "Smooth muscle cells",
  "Horizontal"          = "Horizontal",
  "Astrocyte"           = "Astrocyte",
  "Fibroblast"          = "Fibroblast"
)

#write the exon map
ABCA4_exons_plot <- data.frame()

for (cell_type in names(expressed_genes_by_type)) {
  transcript_ids_of_interest <- expressed_genes_by_type[[cell_type]]
  
  temp_df <- ABCA4_annotation_from_gtf %>%
    filter(
      transcript_id %in% transcript_ids_of_interest,
      type == "exon"
    ) %>%
    mutate(
      cell_type          = cell_type,
      transcript_to_plot = transcript_id
    )
  
  ABCA4_exons_plot <- rbind(ABCA4_exons_plot, temp_df)
}

#Reference exons
ABCA4_reference_exons <- ABCA4_reference_annotation_from_gtf %>%
  filter(type == "exon") %>%
  mutate(
    gene_id            = gene_of_interest,
    transcript_id      = transcript_name,
    cell_type          = "zReference",
    transcript_to_plot = transcript_id
  ) %>%
  select(
    seqnames, start, end, strand, type,
    gene_id, transcript_id, cell_type, transcript_to_plot
  )

#combine PBid + reference
ABCA4_combined_plot <- rbind(ABCA4_exons_plot, ABCA4_reference_exons)

#Cell type order (reference first)
ABCA4_combined_plot <- ABCA4_combined_plot %>%
  mutate(
    cell_type = factor(
      cell_type,
      levels = c(
        "zReference",
        "Amacrine", "Bipolar", "Cone", "Muller_Glia", "Rod",
        "RPE", "Smooth Muscle Cells", "Horizontal",
        "Astrocyte", "Fibroblast"
      )
    )
  )

#Unique track id and PBID label (one row per cell_type × PBID)
ABCA4_combined_plot <- ABCA4_combined_plot %>%
  mutate(
    track_id = paste(cell_type, transcript_to_plot, sep = "__"),
    label    = transcript_to_plot
  )

#Order tracks: by cell_type then PBID, then reverse (top→bottom)
track_levels <- ABCA4_combined_plot %>%
  arrange(cell_type, transcript_to_plot) %>%
  distinct(track_id) %>%
  pull(track_id)

ABCA4_combined_plot <- ABCA4_combined_plot %>%
  mutate(
    track_id = factor(track_id, levels = rev(track_levels))
  )

track_levels <- levels(ABCA4_combined_plot$track_id)  # save for later

#ggtranscript plot
ABCA4_main_plot <- ABCA4_combined_plot %>%
  ggplot(aes(
    xstart = start,
    xend   = end,
    y      = track_id,
    fill   = cell_type
  )) +
  geom_range(
    height = 0.85,
    color  = "black"
  ) +
  geom_intron(
    data = to_intron(ABCA4_combined_plot, "track_id"),
    aes(strand = strand),
    arrow.size = 0.06,
    arrow.min.intron.length = 200,
    arrow.gap = 0.01
  ) +
  scale_fill_manual(values = color_map) +
  scale_y_discrete(
    limits = track_levels,
    labels = function(x) ABCA4_combined_plot$label[
      match(x, ABCA4_combined_plot$track_id)
    ]
  ) +
  labs(
    title = "ABCA4",
    x     = "Genomic position",
    y     = "Isoform (PBID)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title   = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 10),
    legend.position = "none"
  )

#color bar assignment
#one row per track_id with its cell_type
strip_df <- ABCA4_combined_plot %>%
  distinct(track_id, cell_type) %>%
  mutate(
    track_id = factor(track_id, levels = track_levels),
    cell_type_label = cell_type_label_map[as.character(cell_type)]
  ) %>%
  arrange(track_id)

#Numeric y positions (1...n) in the same order as main plot
strip_df <- strip_df %>%
  mutate(y_num = as.numeric(track_id))

#Identify contiguous blocks of the same cell_type
strip_df <- strip_df %>%
  arrange(y_num) %>%
  mutate(
    cell_type_chr = as.character(cell_type),
    block_id = cumsum(
      cell_type_chr != dplyr::lag(cell_type_chr, default = first(cell_type_chr))
    )
  )

#For each block, compute its vertical center for the label
label_df <- strip_df %>%
  group_by(cell_type, cell_type_label, block_id) %>%
  summarise(
    y = mean(y_num),
    .groups = "drop"
  )

#Color bar: one tile per row (track)
ABCA4_strip_plot <- ggplot(strip_df, aes(x = 1, y = track_id, fill = cell_type)) +
  geom_tile(width = 0.6, height = 1) +
  scale_fill_manual(values = color_map) +
  geom_text(
    data = label_df,
    aes(x = 1.5, y = y, label = cell_type_label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4
  ) +
  scale_y_discrete(limits = track_levels) +
  coord_cartesian(xlim = c(0.7, 2.2)) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 0)
  )

#combine main plot + color bar
ABCA4_with_bar <- plot_grid(
  ABCA4_main_plot,
  ABCA4_strip_plot,
  nrow = 1,
  rel_widths = c(10, 1),
  align = "h"
)

print(ABCA4_with_bar)

ggsave(
  filename   = "ABCA4_ggtranscript_FL2_CAGE_with_bar.tiff",
  plot       = ABCA4_with_bar,
  device     = "tiff",
  width      = 20,
  height     = 20,
  dpi        = 600,
  bg         = "WHITE"
)

#define a function to plot other genes
plot_exon_map_with_bar <- function(gene_of_interest,
                                   pbid_fl,
                                   pb_anno,
                                   reference_anno,
                                   expressed_genes_by_type) {
  
  #PBid for gene
  filtered_cage <- pbid_fl %>%
    dplyr::filter(associated_gene == gene_of_interest)
  
  PBid_of_interest <- filtered_cage$isoform
  
  annotation_from_gtf <- pb_anno %>%
    dplyr::filter(
      !is.na(transcript_id),
      transcript_id %in% PBid_of_interest
    ) %>%
    dplyr::select(
      seqnames,
      start,
      end,
      strand,
      type,
      gene_id,
      transcript_id
    )
  
  annotation_from_gtf$gene_id <- gene_of_interest
  
  #Reference GENCODE annotation for gene
  reference_annotation_from_gtf <- reference_anno %>%
    dplyr::filter(
      !is.na(gene_name),
      gene_name %in% gene_of_interest
    ) %>%
    dplyr::select(
      seqnames,
      start,
      end,
      strand,
      type,
      gene_name,
      transcript_name
    )
  
  #color for cell types
  colors <- RColorBrewer::brewer.pal(12, "Paired")
  color_map <- setNames(
    colors[1:11],
    c("Amacrine", "Bipolar", "Cone", "Muller_Glia", "Rod",
      "RPE", "Smooth Muscle Cells", "Horizontal",
      "Astrocyte", "Fibroblast", "zReference")
  )
  
  cell_type_label_map <- c(
    "zReference"          = "Reference",
    "Amacrine"            = "Amacrine",
    "Bipolar"             = "Bipolar",
    "Cone"                = "Cone",
    "Muller_Glia"         = "Muller glia",
    "Rod"                 = "Rod",
    "RPE"                 = "RPE",
    "Smooth Muscle Cells" = "Smooth muscle cells",
    "Horizontal"          = "Horizontal",
    "Astrocyte"           = "Astrocyte",
    "Fibroblast"          = "Fibroblast"
  )
  
  #write the exon map
  exons_plot <- data.frame()
  
  for (cell_type in names(expressed_genes_by_type)) {
    transcript_ids_of_interest <- expressed_genes_by_type[[cell_type]]
    
    temp_df <- annotation_from_gtf %>%
      dplyr::filter(
        transcript_id %in% transcript_ids_of_interest,
        type == "exon"
      ) %>%
      dplyr::mutate(
        cell_type          = cell_type,
        transcript_to_plot = transcript_id
      )
    
    exons_plot <- rbind(exons_plot, temp_df)
  }
  
  #Reference exons
  reference_exons <- reference_annotation_from_gtf %>%
    dplyr::filter(type == "exon") %>%
    dplyr::mutate(
      gene_id            = gene_of_interest,
      transcript_id      = transcript_name,
      cell_type          = "zReference",
      transcript_to_plot = transcript_id
    ) %>%
    dplyr::select(
      seqnames, start, end, strand, type,
      gene_id, transcript_id, cell_type, transcript_to_plot
    )
  
  #combine PBid + reference
  combined_plot <- rbind(exons_plot, reference_exons)
  
  # Cell type order (reference first)
  combined_plot <- combined_plot %>%
    dplyr::mutate(
      cell_type = factor(
        cell_type,
        levels = c(
          "zReference",
          "Amacrine", "Bipolar", "Cone", "Muller_Glia", "Rod",
          "RPE", "Smooth Muscle Cells", "Horizontal",
          "Astrocyte", "Fibroblast"
        )
      )
    )
  
  #Unique track id and PBID label (one row per cell_type × PBID)
  combined_plot <- combined_plot %>%
    dplyr::mutate(
      track_id = paste(cell_type, transcript_to_plot, sep = "__"),
      label    = transcript_to_plot
    )
  
  #Order tracks: by cell_type then PBID, then reverse (top→bottom)
  track_levels <- combined_plot %>%
    dplyr::arrange(cell_type, transcript_to_plot) %>%
    dplyr::distinct(track_id) %>%
    dplyr::pull(track_id)
  
  combined_plot <- combined_plot %>%
    dplyr::mutate(
      track_id = factor(track_id, levels = rev(track_levels))
    )
  
  track_levels <- levels(combined_plot$track_id)  # save for later
  
  #ggtranscript plot
  main_plot <- combined_plot %>%
    ggplot2::ggplot(ggplot2::aes(
      xstart = start,
      xend   = end,
      y      = track_id,
      fill   = cell_type
    )) +
    ggtranscript::geom_range(
      height = 0.85,
      color  = "black"
    ) +
    ggtranscript::geom_intron(
      data = ggtranscript::to_intron(combined_plot, "track_id"),
      ggplot2::aes(strand = strand),
      arrow.size = 0.06,
      arrow.min.intron.length = 200,
      arrow.gap = 0.01
    ) +
    ggplot2::scale_fill_manual(values = color_map) +
    ggplot2::scale_y_discrete(
      limits = track_levels,
      labels = function(x) combined_plot$label[
        match(x, combined_plot$track_id)
      ]
    ) +
    ggplot2::labs(
      title = gene_of_interest,
      x     = "Genomic position",
      y     = "Isoform (PBID)"
    ) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title.x = ggplot2::element_text(size = 16),
      axis.title.y = ggplot2::element_text(size = 16),
      axis.text.x  = ggplot2::element_text(size = 12),
      axis.text.y  = ggplot2::element_text(size = 10),
      legend.position = "none"
    )
  
  #color bar assignment
  #one row per track_id with its cell_type
  strip_df <- combined_plot %>%
    dplyr::distinct(track_id, cell_type) %>%
    dplyr::mutate(
      track_id       = factor(track_id, levels = track_levels),
      cell_type_label = cell_type_label_map[as.character(cell_type)]
    ) %>%
    dplyr::arrange(track_id)
  
  #Numeric y positions (1...n) in the same order as main plot
  strip_df <- strip_df %>%
    dplyr::mutate(y_num = as.numeric(track_id))
  
  #Identify contiguous blocks of the same cell_type
  strip_df <- strip_df %>%
    dplyr::arrange(y_num) %>%
    dplyr::mutate(
      cell_type_chr = as.character(cell_type),
      block_id = cumsum(
        cell_type_chr != dplyr::lag(cell_type_chr, default = first(cell_type_chr))
      )
    )
  
  #For each block, compute its vertical center for the label
  label_df <- strip_df %>%
    dplyr::group_by(cell_type, cell_type_label, block_id) %>%
    dplyr::summarise(
      y = mean(y_num),
      .groups = "drop"
    )
  
  #Color bar: one tile per row (track)
  strip_plot <- ggplot2::ggplot(strip_df, ggplot2::aes(x = 1, y = track_id, fill = cell_type)) +
    ggplot2::geom_tile(width = 0.6, height = 1) +
    ggplot2::scale_fill_manual(values = color_map) +
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(x = 1.5, y = y, label = cell_type_label),
      inherit.aes = FALSE,
      hjust = 0,
      size = 4
    ) +
    ggplot2::scale_y_discrete(limits = track_levels) +
    ggplot2::coord_cartesian(xlim = c(0.7, 2.2)) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 0)
    )
  
  #combine main plot + color bar
  with_bar <- cowplot::plot_grid(
    main_plot,
    strip_plot,
    nrow = 1,
    rel_widths = c(10, 1),
    align = "h"
  )
  
  return(with_bar)
}


#figure 4B
#run and save the result files
ABCA4_with_bar <- plot_exon_map_with_bar(
  gene_of_interest       = "ABCA4",
  pbid_fl                = pbid_fl,
  pb_anno                = pb_anno,
  reference_anno         = reference_anno,
  expressed_genes_by_type = expressed_genes_by_type
)

print(ABCA4_with_bar)

ggsave(
  filename   = "figure4B_ABCA4_ggtranscript_FL2_CAGE_with_bar.tiff",
  plot       = ABCA4_with_bar,
  device     = "tiff",
  width      = 16,
  height     = 20,
  dpi        = 600,
  bg         = "WHITE",
  compression = "lzw"
)

#figure S8, S9
#other genes
EXOC1_with_bar <- plot_exon_map_with_bar(
  gene_of_interest       = "EXOC1",
  pbid_fl                = pbid_fl,
  pb_anno                = pb_anno,
  reference_anno         = reference_anno,
  expressed_genes_by_type = expressed_genes_by_type
)

print(EXOC1_with_bar)

ggsave(
  filename   = "figureS9_EXOC1_ggtranscript_FL2_CAGE_with_bar.tiff",
  plot       = EXOC1_with_bar,
  device     = "tiff",
  width      = 20,
  height     = 20,
  dpi        = 600,
  bg         = "WHITE",
  compression = "lzw"
)

ASPH_with_bar <- plot_exon_map_with_bar(
  gene_of_interest       = "ASPH",
  pbid_fl                = pbid_fl,
  pb_anno                = pb_anno,
  reference_anno         = reference_anno,
  expressed_genes_by_type = expressed_genes_by_type
)

print(ASPH_with_bar)

ggsave(
  filename   = "figureS8_ASPH_ggtranscript_FL2_CAGE_with_bar.tiff",
  plot       = ASPH_with_bar,
  device     = "tiff",
  width      = 20,
  height     = 20,
  dpi        = 600,
  bg         = "WHITE",
  compression = "lzw"
)


##define a function to identify and mark the novel exons
plot_exon_map_with_bar_marknovel <- function(gene_of_interest,
                                             pbid_fl,
                                             pb_anno,
                                             reference_anno,
                                             expressed_genes_by_type,
                                             return_novel_df = FALSE) {
  
  #load PB isoforms & reference annotation
  filtered_cage <- pbid_fl %>%
    dplyr::filter(associated_gene == gene_of_interest)
  PBid_of_interest <- filtered_cage$isoform
  
  annotation_from_gtf <- pb_anno %>%
    dplyr::filter(
      !is.na(transcript_id),
      transcript_id %in% PBid_of_interest
    ) %>%
    dplyr::select(
      seqnames, start, end, strand, type,
      gene_id, transcript_id
    )
  annotation_from_gtf$gene_id <- gene_of_interest
  
  reference_annotation_from_gtf <- reference_anno %>%
    dplyr::filter(
      !is.na(gene_name),
      gene_name %in% gene_of_interest
    ) %>%
    dplyr::select(
      seqnames, start, end, strand, type,
      gene_name, transcript_name
    )
  
  #color palette and cell type labels
  colors <- RColorBrewer::brewer.pal(12, "Paired")
  color_map <- setNames(
    colors[1:11],
    c("Amacrine", "Bipolar", "Cone", "Muller_Glia", "Rod",
      "RPE", "Smooth Muscle Cells", "Horizontal",
      "Astrocyte", "Fibroblast", "zReference")
  )
  
  cell_type_label_map <- c(
    "zReference"          = "Reference",
    "Amacrine"            = "Amacrine",
    "Bipolar"             = "Bipolar",
    "Cone"                = "Cone",
    "Muller_Glia"         = "Muller glia",
    "Rod"                 = "Rod",
    "RPE"                 = "RPE",
    "Smooth Muscle Cells" = "Smooth muscle cells",
    "Horizontal"          = "Horizontal",
    "Astrocyte"           = "Astrocyte",
    "Fibroblast"          = "Fibroblast"
  )
  
  #build PB exon table by cell type
  exons_pb <- data.frame()
  for (cell_type in names(expressed_genes_by_type)) {
    transcript_ids_of_interest <- expressed_genes_by_type[[cell_type]]
    
    temp_df <- annotation_from_gtf %>%
      dplyr::filter(
        transcript_id %in% transcript_ids_of_interest,
        type == "exon"
      ) %>%
      dplyr::mutate(
        cell_type          = cell_type,
        transcript_to_plot = transcript_id
      )
    
    exons_pb <- rbind(exons_pb, temp_df)
  }
  
  #build reference exon table
  reference_exons <- reference_annotation_from_gtf %>%
    dplyr::filter(type == "exon") %>%
    dplyr::mutate(
      gene_id            = gene_of_interest,
      transcript_id      = transcript_name,
      cell_type          = "zReference",
      transcript_to_plot = transcript_id
    ) %>%
    dplyr::select(
      seqnames, start, end, strand, type,
      gene_id, transcript_id, cell_type, transcript_to_plot
    )
  
  #merge PB + reference and set track order
  combined_plot <- rbind(exons_pb, reference_exons)
  
  combined_plot <- combined_plot %>%
    dplyr::mutate(
      cell_type = factor(
        cell_type,
        levels = c(
          "zReference",
          "Amacrine", "Bipolar", "Cone", "Muller_Glia", "Rod",
          "RPE", "Smooth Muscle Cells", "Horizontal",
          "Astrocyte", "Fibroblast"
        )
      ),
      track_id = paste(cell_type, transcript_to_plot, sep = "__"),
      label    = transcript_to_plot
    )
  
  track_levels <- combined_plot %>%
    dplyr::arrange(cell_type, transcript_to_plot) %>%
    dplyr::distinct(track_id) %>%
    dplyr::pull(track_id)
  
  combined_plot <- combined_plot %>%
    dplyr::mutate(track_id = factor(track_id, levels = rev(track_levels)))
  
  track_levels <- levels(combined_plot$track_id)
  
  #identify novel exons using union of all reference exon coordinates
  reference_union <- reference_exons %>%
    dplyr::distinct(seqnames, start, end, strand) %>%
    dplyr::mutate(in_ref = TRUE)
  
  exon_pb_status <- exons_pb %>%
    dplyr::left_join(
      reference_union,
      by = c("seqnames", "start", "end", "strand")
    ) %>%
    dplyr::mutate(
      diff_type = dplyr::if_else(!is.na(in_ref), "known_exon", "novel_exon"),
      track_id  = factor(
        paste(cell_type, transcript_to_plot, sep = "__"),
        levels = track_levels
      )
    )
  
  novel_exons <- exon_pb_status %>% 
    dplyr::filter(diff_type == "novel_exon")
  
  #main exon structure plot with novel exon markers
  main_plot <- combined_plot %>%
    ggplot2::ggplot(ggplot2::aes(
      xstart = start,
      xend   = end,
      y      = track_id,
      fill   = cell_type
    )) +
    ggtranscript::geom_range(
      height = 0.85,
      color  = "black"
    ) +
    ggtranscript::geom_intron(
      data = ggtranscript::to_intron(combined_plot, "track_id"),
      ggplot2::aes(strand = strand),
      arrow.size = 0.06,
      arrow.min.intron.length = 200,
      arrow.gap = 0.01
    ) +
    ggplot2::geom_point(
      data = novel_exons,
      ggplot2::aes(
        x = (start + end) / 2,
        y = track_id
      ),
      inherit.aes = FALSE,
      shape = 17,
      size  = 2,
      color = "red",
      stroke = 1,
      position = ggplot2::position_nudge(y = 0.4)
    ) +
    ggplot2::scale_fill_manual(values = color_map) +
    ggplot2::scale_y_discrete(
      limits = track_levels,
      labels = function(x) combined_plot$label[
        match(x, combined_plot$track_id)
      ]
    ) +
    ggplot2::labs(
      title = gene_of_interest,
      x     = "Genomic position",
      y     = "Isoform (PBID)"
    ) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title.x = ggplot2::element_text(size = 16),
      axis.title.y = ggplot2::element_text(size = 16),
      axis.text.x  = ggplot2::element_text(size = 12),
      axis.text.y  = ggplot2::element_text(size = 10),
      legend.position = "none"
    )
  
  #right-side color bar for cell types
  strip_df <- combined_plot %>%
    dplyr::distinct(track_id, cell_type) %>%
    dplyr::mutate(
      track_id        = factor(track_id, levels = track_levels),
      cell_type_label = cell_type_label_map[as.character(cell_type)]
    ) %>%
    dplyr::arrange(track_id)
  
  strip_df <- strip_df %>%
    dplyr::mutate(y_num = as.numeric(track_id)) %>%
    dplyr::arrange(y_num) %>%
    dplyr::mutate(
      cell_type_chr = as.character(cell_type),
      block_id = cumsum(
        cell_type_chr != dplyr::lag(cell_type_chr, default = first(cell_type_chr))
      )
    )
  
  label_df <- strip_df %>%
    dplyr::group_by(cell_type, cell_type_label, block_id) %>%
    dplyr::summarise(y = mean(y_num), .groups = "drop")
  
  strip_plot <- ggplot2::ggplot(strip_df, ggplot2::aes(x = 1, y = track_id, fill = cell_type)) +
    ggplot2::geom_tile(width = 0.6, height = 1) +
    ggplot2::scale_fill_manual(values = color_map) +
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(x = 1.5, y = y, label = cell_type_label),
      inherit.aes = FALSE,
      hjust = 0,
      size = 4
    ) +
    ggplot2::scale_y_discrete(limits = track_levels) +
    ggplot2::coord_cartesian(xlim = c(0.7, 2.2)) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 0)
    )
  
  with_bar <- cowplot::plot_grid(
    main_plot,
    strip_plot,
    nrow = 1,
    rel_widths = c(10, 1),
    align = "h"
  )
  
  if (return_novel_df) {
    return(list(
      plot        = with_bar,
      novel_exons = novel_exons,
      exon_status = exon_pb_status
    ))
  } else {
    return(with_bar)
  }
}

#figure S7
#run and save figures
res_ABCA4 <- plot_exon_map_with_bar_marknovel(
  gene_of_interest       = "ABCA4",
  pbid_fl                = pbid_fl,
  pb_anno                = pb_anno,
  reference_anno         = reference_anno,
  expressed_genes_by_type = expressed_genes_by_type,
  return_novel_df        = TRUE
)

res_ABCA4$plot 
head(res_ABCA4$novel_exons)

ggsave(
  filename   = "figureS7_ABCA4_ggtranscript_FL2_CAGE_withbar_diffexon.tiff",
  plot       = res_ABCA4$plot,
  device     = "tiff",
  width      = 16,
  height     = 20,
  dpi        = 600,
  bg         = "WHITE",
  compression = "lzw"
)

#statistics about the diff exon
novel_exons <- res_ABCA4$novel_exons
exon_status <- res_ABCA4$exon_status
total_novel_exons <- novel_exons %>%
  distinct(seqnames, start, end, strand) %>%
  nrow()

total_novel_exons

diff_exon_summary_unique <- exon_status %>%
  dplyr::distinct(seqnames, start, end, strand, diff_type) %>%
  dplyr::count(diff_type, name = "n_unique_exons")

diff_exon_summary_unique

diff_exon_summary_occurrence <- exon_status %>%
  dplyr::count(diff_type, name = "n_exon_occurrences")

diff_exon_summary_occurrence

novel_exons_by_celltype_unique <- novel_exons %>%
  dplyr::distinct(cell_type, seqnames, start, end, strand, .keep_all = TRUE) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(
    n_unique_novel_exons = dplyr::n(),
    .groups = "drop"
  )

novel_exons_by_celltype_unique

novel_exons_by_celltype_occurrence <- novel_exons %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(
    n_novel_exon_occurrences = dplyr::n(),
    .groups = "drop"
  )

novel_exons_by_celltype_occurrence

novel_exons_by_celltype_summary <- novel_exons_by_celltype_unique %>%
  dplyr::left_join(
    novel_exons_by_celltype_occurrence,
    by = "cell_type"
  )

novel_exons_by_celltype_summary

novel_exon_recurrence <- novel_exons %>%
  dplyr::group_by(cell_type, seqnames, start, end, strand) %>%
  dplyr::summarise(
    n_transcripts = dplyr::n_distinct(transcript_id),
    transcripts   = paste(sort(unique(transcript_id)), collapse = ";"),
    .groups = "drop"
  ) %>%
  dplyr::arrange(cell_type, seqnames, start, end)

novel_exon_recurrence

write.csv(
  novel_exon_recurrence,
  file = "ABCA4_novel_exons_by_celltype_recurrence.csv",
  row.names = FALSE
)


##plot the splice junctions of PBids used for alphafold3
gene_of_interest <- "ABCA4"

selected_isoforms <- c(
  "PB.5877.497",
  "PB.5877.603",
  "PB.5877.707",
  "PB.5877.738"
)

exons_for_junction <- ABCA4_annotation_from_gtf %>%
  dplyr::filter(
    type == "exon",
    transcript_id %in% selected_isoforms
  ) %>%
  dplyr::mutate(
    transcript_to_plot = transcript_id
  )

exons_for_junction <- exons_for_junction %>%
  dplyr::mutate(
    transcript_to_plot = factor(
      transcript_to_plot,
      levels = rev(unique(selected_isoforms))
    )
  )

junction_for_junction <- junction_anno %>%
  dplyr::filter(
    !is.na(isoform),
    isoform %in% selected_isoforms
  ) %>%
  dplyr::select(
    isoform,
    genomic_start_coord,
    genomic_end_coord,
    strand,
    junction_category
  ) %>%
  dplyr::rename(
    transcript_id = isoform,
    start         = genomic_start_coord,
    end           = genomic_end_coord
  ) %>%
  dplyr::mutate(
    transcript_to_plot = factor(
      transcript_id,
      levels = levels(exons_for_junction$transcript_to_plot)
    )
  )

pbid_levels <- levels(exons_for_junction$transcript_to_plot)
n_iso <- length(pbid_levels)
palette_base <- RColorBrewer::brewer.pal(max(3, min(8, n_iso)), "Set2")
pbid_color_map <- setNames(palette_base[seq_len(n_iso)], pbid_levels)

ABCA4_junction_plot <- exons_for_junction %>%
  ggplot2::ggplot(ggplot2::aes(
    xstart = start,
    xend   = end,
    y      = transcript_to_plot
  )) +
  ggtranscript::geom_range(
    ggplot2::aes(fill = transcript_to_plot),
    color  = "black",
    height = 0.85
  ) +
  ggtranscript::geom_intron(
    data = ggtranscript::to_intron(exons_for_junction, "transcript_to_plot"),
    ggplot2::aes(strand = strand),
    arrow.min.intron.length = 200
  ) +
  ggtranscript::geom_junction(
    data          = junction_for_junction,
    junction.y.max = 0.5
  ) +
  ggtranscript::geom_junction_label_repel(
    data          = junction_for_junction,
    ggplot2::aes(label = junction_category),
    junction.y.max = 0.5,
    size = 6,
    max.overlaps = Inf
  ) +
  ggplot2::scale_fill_manual(
    values = pbid_color_map,
    name   = "Isoform (PBid)"
  ) +
  ggplot2::labs(
    title = "ABCA4 selected isoform junctions",
    x     = "Genomic position",
    y     = "Isoform (PBid)"
  ) +
  ggplot2::theme_minimal(base_size = 16) +
  ggplot2::theme(
    plot.title   = ggplot2::element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = ggplot2::element_text(size = 16),
    axis.title.y = ggplot2::element_text(size = 16),
    axis.text.x  = ggplot2::element_text(size = 12),
    axis.text.y  = ggplot2::element_text(size = 10)
  )

#figure 4D
ggsave(
  filename   = "figure4D_ABCA4_Junctions.tiff",
  plot       = ABCA4_junction_plot,
  width      = 26,
  height     = 8,
  dpi        = 600,
  bg         = "WHITE",
  compression = "lzw"
)

#ASprofile
#run the ASprofile pipeline in cmd
awk -F'\t' '
NR==FNR {
    if (FNR == 1) next  
    if ($22+0 >= 2) ids[$1] = 1
    next
}

{
    if (match($0, /transcript_id "([^"]+)"/, m)) {
        tid = m[1]
        if (tid in ids) {
            print
        }
    }
}
' RetinaMerge_classification.filtered_lite_FL2_classification.txt \
RetinaMerge.sorted.filtered_lite.gff \
> RetinaMerge.sorted.filtered_lite_FL2.gff

./ASprofile.b-1.0.4/extract-as RetinaMerge.sorted.filtered_lite_FL2.gff GRCh38.p14.genome.fa > AS_output.txt

perl ./ASprofile.b-1.0.4/summarize_as.pl RetinaMerge.sorted.filtered_lite_FL2.gff AS_output.txt -p Retina

#plot the result
library(readr)
library(ggplot2)
#read the file
nr_AS <- read.table("Retina.as.nr", header = TRUE, sep = "\t")
AS_summary <- read.table("Retina.as.summary", header = TRUE, sep = "\t")
AS_full <- read.table("AS_output.txt", header = FALSE, sep = "\t")

AS_stat <- aggregate(V9 ~ V2, data = AS_full, FUN = function(x) length(unique(x)))

AS_PBid <- ggplot(AS_stat, aes(x = V9, y = V2)) +
  geom_bar(stat = "identity", fill = "#77A381", orientation = "y") +
  labs(title = "Number of Unique Transcripts per Event Type",
       x = "Count of Unique Transcripts",
       y = "Event Type") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_minimal(base_size = 18)

AS_PBid

#figure 5A
ggsave(
  "ASevent_figure5A.tiff",
  plot  = AS_PBid,
  device = "tiff",
  type   = "cairo",
  width = 8, height = 10, units = "in",
  dpi   = 600,
  compression = "lzw",
  bg = "WHITE"
)

#supplementary figure 5A
class_file <- "RetinaMerge_classification.filtered_lite_FL2_classification.txt"
df <- read_tsv(class_file, show_col_types = FALSE)
df <- df %>%
  mutate(
    sqanti_class = dplyr::case_when(
      structural_category == "full-splice_match"       ~ "FSM",
      structural_category == "incomplete-splice_match" ~ "ISM",
      structural_category == "novel_in_catalog"        ~ "NIC",
      structural_category == "novel_not_in_catalog"    ~ "NNC",
      TRUE                                             ~ "Others"
    )
  )

head(df$within_cage_peak)

per_class <- df %>%
  group_by(sqanti_class) %>%
  summarise(
    n_total = n(),
    n_cage  = sum(within_cage_peak == TRUE, na.rm = TRUE),
    cage_pct = 100 * n_cage / n_total,
    .groups = "drop"
  )

per_class <- per_class %>%
  filter(sqanti_class %in% c("FSM", "ISM", "NIC", "NNC", "Others"))
known_novel <- df %>%
  mutate(group = dplyr::case_when(
    sqanti_class %in% c("FSM", "ISM") ~ "Known",
    sqanti_class %in% c("NIC", "NNC") ~ "Novel",
    TRUE                              ~ NA_character_
  )) %>%
  filter(!is.na(group)) %>%
  group_by(group) %>%
  summarise(
    n_total = n(),
    n_cage  = sum(within_cage_peak == TRUE, na.rm = TRUE),
    cage_pct = 100 * n_cage / n_total,
    .groups = "drop"
  )

df_plot <- bind_rows(
  per_class %>%
    transmute(label = sqanti_class, cage_pct),
  known_novel %>%
    transmute(label = group, cage_pct)
)

# Set factor order so x positions are 1..7 exactly as we want
df_plot$label <- factor(
  df_plot$label,
  levels = c("FSM", "ISM", "NIC", "NNC", "Others", "Known", "Novel")
)

p_supp5A <- ggplot(df_plot, aes(x = label, y = cage_pct)) +
  geom_col(fill = "aquamarine2", width = 0.7) +
  # vertical line between "Others" (5) and "Known" (6)
  geom_vline(xintercept = 5.5, size = 1) +
  scale_y_continuous(
    limits = c(0, 80),
    expand = expansion(mult = c(0, 0.05)),
    name = "CAGE detected (%)"
  ) +
  xlab(NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

print(p_supp5A)
ggsave("CAGE_detection_barplot.png", p_supp5A, width = 5, height = 4, dpi = 600, bg = "WHITE")

#supplementary figure 5B
subcat_stats <- df %>%
  filter(!is.na(subcategory)) %>%
  group_by(subcategory) %>%
  summarise(
    n_total = n(),
    n_cage  = sum(within_cage_peak == TRUE, na.rm = TRUE),
    cage_pct = 100 * n_cage / n_total,
    .groups = "drop"
  )

subcat_stats <- subcat_stats %>%
  arrange(desc(cage_pct))

subcat_stats$subcategory <- factor(
  subcat_stats$subcategory,
  levels = subcat_stats$subcategory
)

p_supp5B <- ggplot(subcat_stats, aes(x = subcategory, y = cage_pct)) +
  geom_col(fill = "aquamarine2", width = 0.7) +
  scale_y_continuous(
    name = "CAGE detected (%)",
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.05))
  ) +
  xlab(NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

print(p_supp5B)

ggsave("CAGE_detection_by_subcategory.png", p_supp5B,
       width = 7, height = 4.5, dpi = 600, bg = "WHITE")


#supplementary figure 5C
df <- df %>%
  mutate(
    sqanti_class = case_when(
      structural_category == "full-splice_match" ~ "FSM",
      structural_category == "incomplete-splice_match" ~ "ISM",
      structural_category == "novel_in_catalog" ~ "NIC",
      structural_category == "novel_not_in_catalog" ~ "NNC",
      TRUE ~ "Others"
    )
  )

polyA_class <- df %>%
  group_by(sqanti_class) %>%
  summarise(
    n_total = n(),
    n_polyA = sum(!is.na(polyA_motif) & polyA_motif != "", na.rm = TRUE),
    polyA_pct = 100 * n_polyA / n_total,
    .groups = "drop"
  ) %>%
  filter(sqanti_class %in% c("FSM", "ISM", "NIC", "NNC", "Others"))
polyA_known_novel <- df %>%
  mutate(group = case_when(
    sqanti_class %in% c("FSM", "ISM") ~ "Known",
    sqanti_class %in% c("NIC", "NNC") ~ "Novel",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group)) %>%
  group_by(group) %>%
  summarise(
    n_total = n(),
    n_polyA = sum(!is.na(polyA_motif) & polyA_motif != "", na.rm = TRUE),
    polyA_pct = 100 * n_polyA / n_total,
    .groups = "drop"
  )

df_plot <- bind_rows(
  polyA_class %>%
    transmute(label = sqanti_class, polyA_pct),
  polyA_known_novel %>%
    transmute(label = group, polyA_pct)
)

df_plot$label <- factor(
  df_plot$label,
  levels = c("FSM", "ISM", "NIC", "NNC", "Others", "Known", "Novel")
)
p_supp5C <- ggplot(df_plot, aes(x = label, y = polyA_pct)) +
  geom_col(fill = "salmon", width = 0.7) +
  geom_vline(xintercept = 5.5, size = 1) +
  scale_y_continuous(
    limits = c(0, 60),
    expand = expansion(mult = c(0, 0.05)),
    name = "PolyA motif detected (%)"
  ) +
  xlab(NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

print(p_supp5C)

ggsave("PolyA_motif_barplot.png", p_supp5C, width = 5, height = 4, dpi = 600, bg = "WHITE")


#validation with a bulk retina dataset
#SJ match ratio
df <- read_csv("SJ_loose_match_summary.csv")
print(df)
df_1bp <- df %>% filter(tol_bp == 1)
df_plot <- tibble(
  Category = c("FSM / ISM (known)", "NIC / NNC (novel)"),
  MatchRate = c(df_1bp$FSM_ISM_rate, df_1bp$NNC_NIC_rate)
)
#plot the SJ match ratio
p_SJ <- ggplot(df_plot, aes(x = Category, y = MatchRate, fill = Category)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = percent(MatchRate, accuracy = 0.1)),
            vjust = -0.5, size = 5) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.05)) +
  labs(
    title = "Validation of splice junctions",
    subtitle = "Iso-Seq junctions validated against RJunBase",
    x = NULL, y = "Matched junctions (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12)
  )
print(p_SJ)
#save the plot
#supplmentary 2B
ggsave("SJ_matchrate_barplot.png", p_SJ, width = 8, height = 6, dpi = 600, bg = 'WHITE')

#transcript match ratio
known_file <- "sc_known_2FL2bulk_known/RetinaMerge.FSM_ISM_FL2_CAGE_classification.txt"
novel_file <- "sc_novel_2FL2bulk_novel/RetinaMerge.NOVEL_FL2_CAGE_classification.txt"

known <- read_tsv(known_file, show_col_types = FALSE)
novel <- read_tsv(novel_file, show_col_types = FALSE)
colnames(known)
colnames(novel)

known <- known %>%
  mutate(
    mapped = !is.na(associated_transcript) & grepl("ENST", associated_transcript)
  )

novel <- novel %>%
  mutate(
    mapped = !is.na(associated_transcript) & grepl("^trans", associated_transcript)
  )

summary_known <- known %>%
  summarise(
    dataset       = "Known (FSM/ISM)",
    total_isoform = n(),
    mapped        = sum(mapped),
    unmapped      = total_isoform - mapped,
    prop_mapped   = mapped / total_isoform
  )

summary_novel <- novel %>%
  summarise(
    dataset       = "Novel (NIC/NNC)",
    total_isoform = n(),
    mapped        = sum(mapped),
    unmapped      = total_isoform - mapped,
    prop_mapped   = mapped / total_isoform
  )

summary_known
summary_novel

#plot the match ratio
df_plot_prop <- bind_rows(summary_known, summary_novel) %>%
  select(dataset, prop_mapped)

p_bulk_match <- ggplot(df_plot_prop, aes(x = dataset, y = prop_mapped, fill = dataset)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = percent(prop_mapped, accuracy = 0.1)),
            vjust = -0.5, size = 5) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.05)) +
  labs(
    title = "Ratio of Iso-Seq isoforms matched with bulk evidence",
    x = NULL,
    y = "Mapped isoforms (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    plot.title = element_text(face = "bold")
  )

print(p_bulk_match)
#save the plot
#supplementary 2A
ggsave("transcript_matchrate_barplot.png", p_bulk_match, width = 8, height = 6, dpi = 600, bg = 'WHITE')

#plot for transcript1884 and its identical transcript PB.32297.5
similar_SLC6A11_transcriptid <- c("PB.32297.5", "transcript1884.chr3.nnic")
SLC6A11_final_plot_rescaled_selected <- SLC6A11_final_plot_rescaled[SLC6A11_final_plot_rescaled$transcript_id %in% similar_SLC6A11_transcriptid, ]
#construct the modified table
SLC6A11_final_plot_rescaled_selected_exons <- SLC6A11_final_plot_rescaled_selected %>% dplyr::filter(type == "exon") 
SLC6A11_final_plot_rescaled_selected_introns <- SLC6A11_final_plot_rescaled_selected %>% dplyr::filter(type == "intron") 
#plotting
SLC6A11_rescaled_selected_nosplit <- SLC6A11_final_plot_rescaled_selected_exons %>% 
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_to_plot
  )) +
  geom_range(
    aes(fill = cell_type)
  ) +
  geom_intron(
    data = SLC6A11_final_plot_rescaled_selected_introns,
    aes(strand = strand), 
    arrow.min.intron.length = 300
  )
SLC6A11_rescaled_selected_nosplit
#supplementary 2C
ggsave("SLC6A11_nosplit_rescaled_identical.tiff", plot = SLC6A11_rescaled_selected_nosplit, width = 12, height = 6, dpi = 1200, limitsize = FALSE, bg = "white")


#pigeon report was generated using the scripts provided on SQANTI3 GitHub (https://github.com/ConesaLab/SQANTI3/tree/master/utilities/report_pigeon)

####end of the session########################
####author: Raymond Wong, Luozixian Wang######
####05/12/2025################################
