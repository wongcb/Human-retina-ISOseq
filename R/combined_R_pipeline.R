########Iso-seq profiling of novel isoform variants in human retina at single cell resolution############
########created by Luozixian Wang, Raymond Wong##########################################################
########CERA, UNIMELB, 18/06/2026########################################################################

#this combined script contain several analysis modules
#better to run every module individually to avoid crash in names of files, functions, and parameters

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
library(dplyr)
library(tidyr)
library(tibble)
library(grid)
library(AnnotationDbi)
library(scDblFinder)
library(openxlsx)
library(readr)
library(stringr)
library(forcats)
theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 14)
organism = "org.Hs.eg.db"

#the long reads and short reads were preprocessed separately
#generally, long reads were preprocessed following the isoseq3 workflow
#short reads were preprocessed by cellranger

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

#supplementary figure 1A
figureS1A <- DimPlot(merged_retina, reduction = "umap", group.by = "sample.id", pt.size = 0.5)
#supplementary figure 1B
figureS1B <- DimPlot(merged_retina, group.by ='read_stat', reduction='umap', pt.size = 0.5)

combined_fS1AB <- figureS1A | figureS1B
ggsave("figureS1AB.tiff", combined_fS1AB, width = 16, height = 8, dpi = 600, bg = "WHITE")

#GO term enrichment analysis
symbol_ensid <- read.delim('D:/bioinformatic/isoseq/processed/from_Ray/revision/redo_FL2CAGEpolyA/GO/symbol_ensid.txt', header = TRUE, stringsAsFactors = FALSE)
selected_symbol_ensid <- symbol_ensid[, c("HGNC.symbol", "Gene.stable.ID")]
filtered_symbol_ensid <- selected_symbol_ensid %>%
  group_by(HGNC.symbol) %>%
  slice(1) %>%
  ungroup()

FSM_transcript_FL2_cage_motif <- transcript_status_FL2_cage_motif %>%
  filter(category == "full-splice_match")

ISM_transcript_FL2_cage_motif <- transcript_status_FL2_cage_motif %>%
  filter(category == "incomplete-splice_match")

NNC_transcript_FL2_cage_motif <- transcript_status_FL2_cage_motif %>%
  filter(category == "novel_not_in_catalog")

NIC_transcript_FL2_cage_motif <- transcript_status_FL2_cage_motif %>%
  filter(category == "novel_in_catalog")

gene_FSM_FL2_cage_motif <- data.frame(gene_FSM = unique(FSM_transcript_FL2_cage_motif$gene),
                                      stringsAsFactors = FALSE)
gene_ISM_FL2_cage_motif <- data.frame(gene_ISM = unique(ISM_transcript_FL2_cage_motif$gene),
                                      stringsAsFactors = FALSE)
gene_NNC_FL2_cage_motif <- data.frame(gene_NNC = unique(NNC_transcript_FL2_cage_motif$gene),
                                      stringsAsFactors = FALSE)
gene_NIC_FL2_cage_motif <- data.frame(gene_NIC = unique(NIC_transcript_FL2_cage_motif$gene),
                                      stringsAsFactors = FALSE)

FSM_converted_FL2_cage_motif <- gene_FSM_FL2_cage_motif %>%
  left_join(filtered_symbol_ensid, by = c(gene_FSM = "HGNC.symbol"))

ISM_converted_FL2_cage_motif <- gene_ISM_FL2_cage_motif %>%
  left_join(filtered_symbol_ensid, by = c(gene_ISM = "HGNC.symbol"))

NNC_converted_FL2_cage_motif <- gene_NNC_FL2_cage_motif %>%
  left_join(filtered_symbol_ensid, by = c(gene_NNC = "HGNC.symbol"))

NIC_converted_FL2_cage_motif <- gene_NIC_FL2_cage_motif %>%
  left_join(filtered_symbol_ensid, by = c(gene_NIC = "HGNC.symbol"))

FSM_filled_FL2_cage_motif <- FSM_converted_FL2_cage_motif %>%
  mutate(Gene.stable.ID = ifelse(is.na(Gene.stable.ID),
                                 gene_FSM, Gene.stable.ID))

ISM_filled_FL2_cage_motif <- ISM_converted_FL2_cage_motif %>%
  mutate(Gene.stable.ID = ifelse(is.na(Gene.stable.ID),
                                 gene_ISM, Gene.stable.ID))

NNC_filled_FL2_cage_motif <- NNC_converted_FL2_cage_motif %>%
  mutate(Gene.stable.ID = ifelse(is.na(Gene.stable.ID),
                                 gene_NNC, Gene.stable.ID))

NIC_filled_FL2_cage_motif <- NIC_converted_FL2_cage_motif %>%
  mutate(Gene.stable.ID = ifelse(is.na(Gene.stable.ID),
                                 gene_NIC, Gene.stable.ID))


#check duplicate
FSM_duplicated_genes_FL2_cage_motif <- FSM_filled_FL2_cage_motif %>%
  filter(duplicated(Gene.stable.ID) | duplicated(Gene.stable.ID, fromLast = TRUE))

ISM_duplicated_genes_FL2_cage_motif <- ISM_filled_FL2_cage_motif %>%
  filter(duplicated(Gene.stable.ID) | duplicated(Gene.stable.ID, fromLast = TRUE))

NNC_duplicated_genes_FL2_cage_motif <- NNC_filled_FL2_cage_motif %>%
  filter(duplicated(Gene.stable.ID) | duplicated(Gene.stable.ID, fromLast = TRUE))

NIC_duplicated_genes_FL2_cage_motif <- NIC_filled_FL2_cage_motif %>%
  filter(duplicated(Gene.stable.ID) | duplicated(Gene.stable.ID, fromLast = TRUE))

FSM_genes_FL2_cage_motif_for_GO <- unique(FSM_filled_FL2_cage_motif$Gene.stable.ID)
ISM_genes_FL2_cage_motif_for_GO <- unique(ISM_filled_FL2_cage_motif$Gene.stable.ID)
NNC_genes_FL2_cage_motif_for_GO <- unique(NNC_filled_FL2_cage_motif$Gene.stable.ID)
NIC_genes_FL2_cage_motif_for_GO <- unique(NIC_filled_FL2_cage_motif$Gene.stable.ID)

#GO enrichment
FSM_GO_FL2_cage_motif <- enrichGO(
  gene          = FSM_genes_FL2_cage_motif_for_GO,
  ont           = "ALL",
  keyType       = "ENSEMBL",
  pvalueCutoff  = 0.05,
  OrgDb         = "org.Hs.eg.db",
  pAdjustMethod = "fdr",
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

ISM_GO_FL2_cage_motif <- enrichGO(
  gene          = ISM_genes_FL2_cage_motif_for_GO,
  ont           = "ALL",
  keyType       = "ENSEMBL",
  pvalueCutoff  = 0.05,
  OrgDb         = "org.Hs.eg.db",
  pAdjustMethod = "fdr",
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

NNC_GO_FL2_cage_motif <- enrichGO(
  gene          = NNC_genes_FL2_cage_motif_for_GO,
  ont           = "ALL",
  keyType       = "ENSEMBL",
  pvalueCutoff  = 0.05,
  OrgDb         = "org.Hs.eg.db",
  pAdjustMethod = "fdr",
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

NIC_GO_FL2_cage_motif <- enrichGO(
  gene          = NIC_genes_FL2_cage_motif_for_GO,
  ont           = "ALL",
  keyType       = "ENSEMBL",
  pvalueCutoff  = 0.05,
  OrgDb         = "org.Hs.eg.db",
  pAdjustMethod = "fdr",
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

#plot GO
pFSM_dot_FL2_cage_motif <- dotplot(FSM_GO_FL2_cage_motif, showCategory = 10) + ggtitle("FSM (FL2_cage_motif)")
pISM_dot_FL2_cage_motif <- dotplot(ISM_GO_FL2_cage_motif, showCategory = 10) + ggtitle("ISM (FL2_cage_motif)")
pNNC_dot_FL2_cage_motif <- dotplot(NNC_GO_FL2_cage_motif, showCategory = 10) + ggtitle("NNC (FL2_cage_motif)")
pNIC_dot_FL2_cage_motif <- dotplot(NIC_GO_FL2_cage_motif, showCategory = 10) + ggtitle("NIC (FL2_cage_motif)")

#supplementary figure 1C-F
print(pFSM_dot_FL2_cage_motif)
print(pISM_dot_FL2_cage_motif)
print(pNNC_dot_FL2_cage_motif)
print(pNIC_dot_FL2_cage_motif)
combined_fS1CF <- (pFSM_dot_FL2_cage_motif | pISM_dot_FL2_cage_motif) /
  (pNNC_dot_FL2_cage_motif | pNIC_dot_FL2_cage_motif)
ggsave("figureS1CDEF.tiff",
       combined_fS1CF,
       width = 15, height = 14,
       dpi = 600, bg = "WHITE", compression = "lzw")


#alternative promoter analysis using proActiv
#before all, run a split_bam.py in python to split the bam file generated by isoseq3
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
#figure2A
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

#supplementary figure 6A
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

#supplementary figure 6B
target_genes_v2 <- c("GPX4", "GUK1", "PPP2R2B")
min_gene_celltype_total_count_v2 <- 10
n_boot_v2 <- 1000
filterInternal_v2 <- FALSE

celltype_order <- c(
  "Amacrine", "Astrocyte", "Bipolar", "Cone", "Fibroblast",
  "Horizontal", "Muller_Glia", "RPE", "Rod", "Smooth_Muscle_Cells"
)

relative_fig_width_v2 <- 32
relative_fig_height_v2 <- 20

ci_fig_width_v2 <- 22
ci_fig_height_v2 <- 7

out_dir_v2 <- "Figure2CDE_BootstrapCI_v2_rescreen_celltypeTotal10"
dir.create(out_dir_v2, showWarnings = FALSE)
dir.create(file.path(out_dir_v2, "tables"), showWarnings = FALSE)

clean_ensembl <- function(x) {
  sub("\\..*$", "", as.character(x))
}

clean_celltype <- function(x) {
  gsub("_FL2_CAGE_polyA$", "", as.character(x))
}

make_promoter_color_map <- function(x) {
  promoter_ids <- unique(as.character(x))
  n_promoters <- length(promoter_ids)
  
  if (n_promoters <= 0) {
    return(character(0))
  }
  
  if (n_promoters <= 8) {
    cols <- RColorBrewer::brewer.pal(max(3, n_promoters), "Set2")[seq_len(n_promoters)]
  } else if (n_promoters <= 12) {
    cols <- RColorBrewer::brewer.pal(n_promoters, "Set3")
  } else {
    cols <- grDevices::rainbow(n_promoters)
  }
  
  setNames(cols, promoter_ids)
}

empty_plot_v2 <- function(title_text = "") {
  ggplot() +
    theme_void() +
    labs(title = title_text) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18)
    )
}

rdata_v2 <- as.data.frame(rowData(result)) %>%
  tibble::rownames_to_column("feature_id") %>%
  dplyr::mutate(
    geneId_full = as.character(geneId),
    geneId_clean = clean_ensembl(geneId),
    promoterId_clean = paste0("prmtr.", promoterId)
  )

symbol_to_ens_v2 <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = target_genes_v2,
  column = "ENSEMBL",
  keytype = "SYMBOL",
  multiVals = "list"
)

gene_map_v2 <- lapply(target_genes_v2, function(sym) {
  
  ens_candidates <- as.character(symbol_to_ens_v2[[sym]])
  ens_candidates <- ens_candidates[!is.na(ens_candidates)]
  
  matched <- rdata_v2 %>%
    dplyr::filter(geneId_clean %in% ens_candidates) %>%
    dplyr::distinct(gene_symbol = sym, geneId_full, geneId_clean)
  
  if (nrow(matched) == 0) {
    stop(paste0(
      "Cannot find ", sym, " in result. Mapped Ensembl IDs: ",
      paste(ens_candidates, collapse = ", ")
    ))
  }
  
  matched %>%
    dplyr::slice(1)
  
}) %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(factor(gene_symbol, levels = target_genes_v2))

write.csv(
  gene_map_v2,
  file.path(out_dir_v2, "tables", "v2_gene_symbol_to_ensembl_mapping.csv"),
  row.names = FALSE
)

print(gene_map_v2)

#extract gene matrices from result
get_gene_matrices_v2 <- function(result,
                                 gene_id_clean,
                                 gene_symbol,
                                 filterInternal = FALSE) {
  
  rdata <- as.data.frame(rowData(result)) %>%
    tibble::rownames_to_column("feature_id") %>%
    dplyr::mutate(
      geneId_full = as.character(geneId),
      geneId_clean = clean_ensembl(geneId),
      promoterId_clean = paste0("prmtr.", promoterId)
    )
  
  gene_idx <- which(rdata$geneId_clean == gene_id_clean)
  
  if (length(gene_idx) == 0) {
    stop(paste0("Gene not found: ", gene_symbol, " / ", gene_id_clean))
  }
  
  tmp_rdata <- rdata[gene_idx, , drop = FALSE]
  tmp_counts <- assays(result)$promoterCounts[gene_idx, , drop = FALSE]
  
  keep_idx <- rep(TRUE, nrow(tmp_rdata))
  
  if (filterInternal) {
    keep_idx <- tmp_rdata$internalPromoter == FALSE
  }
  
  tmp_rdata <- tmp_rdata[keep_idx, , drop = FALSE]
  tmp_counts <- tmp_counts[keep_idx, , drop = FALSE]
  
  promoter_ids <- paste0("prmtr.", tmp_rdata$promoterId)
  rownames(tmp_counts) <- promoter_ids
  
  list(
    rdata = tmp_rdata,
    counts = tmp_counts,
    promoter_ids = promoter_ids
  )
}

#retain celltypes directly from result
get_gene_celltype_support_from_result_v2 <- function(result,
                                                     gene_id_clean,
                                                     gene_symbol,
                                                     filterInternal = FALSE) {
  
  gene_obj <- get_gene_matrices_v2(
    result = result,
    gene_id_clean = gene_id_clean,
    gene_symbol = gene_symbol,
    filterInternal = filterInternal
  )
  
  tmp_counts <- gene_obj$counts
  
  support_df <- as.data.frame(tmp_counts) %>%
    tibble::rownames_to_column("promoterId") %>%
    tidyr::pivot_longer(
      cols = -promoterId,
      names_to = "celltype_raw",
      values_to = "promoter_count"
    ) %>%
    dplyr::mutate(
      gene_symbol = gene_symbol,
      geneId_clean = gene_id_clean,
      celltype = clean_celltype(celltype_raw),
      promoter_count = as.numeric(promoter_count)
    ) %>%
    dplyr::group_by(gene_symbol, geneId_clean, celltype) %>%
    dplyr::summarise(
      total_promoter_count = sum(promoter_count[promoter_count > 0], na.rm = TRUE),
      n_promoters_detected = sum(promoter_count > 0, na.rm = TRUE),
      retained_celltype = total_promoter_count >= min_gene_celltype_total_count_v2,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      celltype = factor(celltype, levels = celltype_order)
    ) %>%
    dplyr::arrange(
      factor(gene_symbol, levels = target_genes_v2),
      celltype
    )
  
  support_df
}

celltype_support_v2 <- lapply(seq_len(nrow(gene_map_v2)), function(i) {
  get_gene_celltype_support_from_result_v2(
    result = result,
    gene_id_clean = gene_map_v2$geneId_clean[i],
    gene_symbol = gene_map_v2$gene_symbol[i],
    filterInternal = filterInternal_v2
  )
}) %>%
  dplyr::bind_rows()

write.csv(
  celltype_support_v2,
  file.path(out_dir_v2, "tables", "v2_rescreened_gene_celltype_total_count_support.csv"),
  row.names = FALSE
)

print(celltype_support_v2)

keep_celltypes_by_gene_v2 <- lapply(target_genes_v2, function(gene_use) {
  celltype_support_v2 %>%
    dplyr::filter(
      gene_symbol == gene_use,
      retained_celltype
    ) %>%
    dplyr::pull(celltype) %>%
    as.character()
})

names(keep_celltypes_by_gene_v2) <- target_genes_v2

keep_celltype_table_v2 <- dplyr::bind_rows(
  lapply(names(keep_celltypes_by_gene_v2), function(gene_use) {
    tibble(
      gene_symbol = gene_use,
      retained_celltype = keep_celltypes_by_gene_v2[[gene_use]]
    )
  })
)

write.csv(
  keep_celltype_table_v2,
  file.path(out_dir_v2, "tables", "v2_retained_celltypes_used_for_CI_and_Figure2CDE.csv"),
  row.names = FALSE
)

print(keep_celltype_table_v2)

#re-calculate CI for retained celltypes
bootstrap_promoter_usage_from_counts <- function(count_vec, n_boot = 1000) {
  
  promoter_names <- names(count_vec)
  count_vec <- as.numeric(count_vec)
  names(count_vec) <- promoter_names
  
  # Only remove NA and zero-count promoters.
  # Do NOT use promoter_count >= 10 here.
  count_vec <- count_vec[!is.na(count_vec)]
  count_vec <- count_vec[count_vec > 0]
  
  total_count <- sum(count_vec, na.rm = TRUE)
  
  if (length(count_vec) == 0 || is.na(total_count) || total_count <= 0) {
    return(NULL)
  }
  
  if (length(count_vec) == 1) {
    return(data.frame(
      promoterId = names(count_vec)[1],
      mean_usage = 1,
      median_usage = 1,
      ci_lower = 1,
      ci_upper = 1,
      stringsAsFactors = FALSE
    ))
  }
  
  probs <- count_vec / total_count
  
  boot_mat <- sapply(seq_len(n_boot), function(i) {
    as.vector(rmultinom(1, size = total_count, prob = probs))
  })
  
  boot_mat <- as.matrix(boot_mat)
  
  if (nrow(boot_mat) != length(count_vec)) {
    boot_mat <- matrix(boot_mat, nrow = length(count_vec))
  }
  
  rownames(boot_mat) <- names(count_vec)
  colnames(boot_mat) <- paste0("boot_", seq_len(ncol(boot_mat)))
  
  boot_df <- as.data.frame(t(boot_mat))
  boot_df$total <- rowSums(boot_df)
  
  boot_long <- boot_df %>%
    tibble::rownames_to_column("boot_id") %>%
    tidyr::pivot_longer(
      cols = -c(boot_id, total),
      names_to = "promoterId",
      values_to = "count"
    ) %>%
    dplyr::mutate(
      usage = ifelse(total > 0, count / total, NA_real_)
    )
  
  boot_long %>%
    dplyr::group_by(promoterId) %>%
    dplyr::summarise(
      mean_usage = mean(usage, na.rm = TRUE),
      median_usage = median(usage, na.rm = TRUE),
      ci_lower = quantile(usage, 0.025, na.rm = TRUE),
      ci_upper = quantile(usage, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

get_promoter_ci_for_gene_rescreened_v2 <- function(result,
                                                   gene_id_clean,
                                                   gene_symbol,
                                                   filterInternal = FALSE,
                                                   n_boot = 1000) {
  
  gene_obj <- get_gene_matrices_v2(
    result = result,
    gene_id_clean = gene_id_clean,
    gene_symbol = gene_symbol,
    filterInternal = filterInternal
  )
  
  tmp_counts <- gene_obj$counts
  
  keep_celltypes <- keep_celltypes_by_gene_v2[[gene_symbol]]
  
  if (is.null(keep_celltypes) || length(keep_celltypes) == 0) {
    warning(paste0("No retained celltypes for ", gene_symbol))
    return(tibble())
  }
  
  all_celltypes_raw <- colnames(tmp_counts)
  
  out_list <- list()
  idx <- 1
  
  for (ct_raw in all_celltypes_raw) {
    
    ct_clean <- clean_celltype(ct_raw)
    
    if (!ct_clean %in% keep_celltypes) next
    
    count_vec_all <- tmp_counts[, ct_raw]
    names(count_vec_all) <- rownames(tmp_counts)
    
    count_vec_detected <- count_vec_all[
      !is.na(count_vec_all) &
        count_vec_all > 0
    ]
    
    total_count_detected <- sum(count_vec_detected, na.rm = TRUE)
    n_promoters_detected <- length(count_vec_detected)
    
    # Safety check: this should already be true because of keep_celltypes
    if (total_count_detected < min_gene_celltype_total_count_v2) next
    if (n_promoters_detected == 0) next
    
    ci_df <- bootstrap_promoter_usage_from_counts(
      count_vec_detected,
      n_boot = n_boot
    )
    
    if (is.null(ci_df) || nrow(ci_df) == 0) next
    
    ci_df <- ci_df %>%
      dplyr::mutate(
        gene_symbol = gene_symbol,
        geneId_clean = gene_id_clean,
        celltype = ct_clean,
        total_count_detected = total_count_detected,
        n_promoters_detected = n_promoters_detected
      )
    
    out_list[[idx]] <- ci_df
    idx <- idx + 1
  }
  
  dplyr::bind_rows(out_list)
}

ci_list_v2 <- lapply(seq_len(nrow(gene_map_v2)), function(i) {
  
  gene_symbol_use <- gene_map_v2$gene_symbol[i]
  gene_id_clean_use <- gene_map_v2$geneId_clean[i]
  
  message("Re-calculating CI with rescreened celltypes for ", gene_symbol_use)
  
  get_promoter_ci_for_gene_rescreened_v2(
    result = result,
    gene_id_clean = gene_id_clean_use,
    gene_symbol = gene_symbol_use,
    filterInternal = filterInternal_v2,
    n_boot = n_boot_v2
  )
})

names(ci_list_v2) <- gene_map_v2$gene_symbol

ci_all_v2 <- dplyr::bind_rows(ci_list_v2)

write.csv(
  ci_all_v2,
  file.path(out_dir_v2, "tables", "v2_bootstrap_CI_rescreened_celltypeTotal10.csv"),
  row.names = FALSE
)

for (gene_use in target_genes_v2) {
  write.csv(
    ci_list_v2[[gene_use]],
    file.path(out_dir_v2, "tables", paste0("celltype_proActiv_bootstrap_CI_", gene_use, "_v2_rescreened.csv")),
    row.names = FALSE
  )
}

#plot CI
#supplementary figure 6B
plot_ci_rescreened_v2 <- function(ci_df, gene_label) {
  
  if (is.null(ci_df) || nrow(ci_df) == 0) {
    return(empty_plot_v2(gene_label))
  }
  
  ci_df <- ci_df %>%
    dplyr::mutate(
      celltype = clean_celltype(celltype),
      celltype = factor(celltype, levels = celltype_order)
    )
  
  color_map <- make_promoter_color_map(ci_df$promoterId)
  
  ggplot(
    ci_df,
    aes(
      x = celltype,
      y = mean_usage,
      color = promoterId,
      group = promoterId
    )
  ) +
    geom_point(
      position = position_dodge(width = 0.5),
      size = 3
    ) +
    geom_errorbar(
      aes(ymin = ci_lower, ymax = ci_upper),
      position = position_dodge(width = 0.5),
      width = 0.20,
      linewidth = 0.7
    ) +
    scale_color_manual(values = color_map) +
    scale_x_discrete(
      drop = TRUE,
      expand = expansion(add = 0.5)
    ) +
    scale_y_continuous(
      labels = scales::percent_format(),
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25)
    ) +
    labs(
      title = gene_label,
      x = "Cell type",
      y = "Promoter usage",
      color = "Promoter"
    ) +
    theme_cowplot() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 13, face = "bold"),
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11, face = "bold")
    ) +
    guides(
      color = guide_legend(nrow = 2)
    )
}

p_ci_list_v2 <- lapply(target_genes_v2, function(gene_use) {
  plot_ci_rescreened_v2(
    ci_list_v2[[gene_use]],
    gene_use
  )
})

names(p_ci_list_v2) <- target_genes_v2

celltype_bootstrap_plot_v2 <- plot_grid(
  plotlist = p_ci_list_v2,
  nrow = 1,
  rel_widths = rep(1, length(p_ci_list_v2)),
  align = "h"
)

celltype_bootstrap_plot_v2

ggsave(
  filename = file.path(out_dir_v2, "celltype_proActiv_bootstrap_CI_examples_subset_celltypes_v2.tiff"),
  plot = celltype_bootstrap_plot_v2,
  device = "tiff",
  type = "cairo",
  width = 20,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "WHITE"
)

ggsave(
  filename = "celltype_proActiv_bootstrap_CI_examples_subset_celltypes_v2.tiff",
  plot = celltype_bootstrap_plot_v2,
  device = "tiff",
  type = "cairo",
  width = 20,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "WHITE"
)

#supplementary figure 6B relative dotplot
FS <- list(
  base         = 24,
  axis_title   = 27,
  title        = 30,
  bigtitle     = 36,
  legend_title = 27,
  legend_point = 14
)

decorate_plot_v2 <- function(p, gene_name) {
  p +
    labs(
      title = gene_name,
      x = "Promoter ID",
      y = "Relative promoter activity"
    ) +
    scale_y_continuous(
      labels = scales::percent_format(),
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    theme(
      plot.title = element_text(size = FS$title, face = "bold", hjust = 0.5),
      axis.title.x = element_text(
        size = FS$axis_title,
        face = "bold",
        margin = margin(t = 55)
      ),
      axis.title.y = element_text(
        size = FS$axis_title,
        face = "bold",
        margin = margin(r = 6)
      ),
      axis.text.x = element_text(
        size = FS$base,
        angle = 45,
        hjust = 1,
        vjust = 1,
        margin = margin(t = 8)
      ),
      axis.text.y = element_text(size = FS$base),
      legend.title = element_text(size = FS$legend_title, face = "bold"),
      legend.text = element_text(size = FS$base),
      legend.position = "bottom",
      plot.margin = margin(10, 15, 10, 15)
    ) +
    guides(
      colour = guide_legend(override.aes = list(size = FS$legend_point)),
      fill = guide_legend(override.aes = list(size = FS$legend_point), nrow = 2)
    )
}

subset_dotplot_to_rescreened_celltypes_v2 <- function(p, gene_symbol) {
  
  keep_celltypes <- keep_celltypes_by_gene_v2[[gene_symbol]]
  
  if (is.null(keep_celltypes) || length(keep_celltypes) == 0) {
    warning(paste0("No retained celltypes for ", gene_symbol))
    return(p)
  }
  
  p_sub <- p
  
  if (!"Condition" %in% colnames(p_sub$data)) {
    stop(
      paste0(
        "Cannot find 'Condition' column in dotplot data for ",
        gene_symbol,
        ". Available columns: ",
        paste(colnames(p_sub$data), collapse = ", ")
      )
    )
  }
  
  p_sub$data <- p_sub$data %>%
    dplyr::mutate(
      Condition_clean = clean_celltype(Condition)
    ) %>%
    dplyr::filter(
      Condition_clean %in% keep_celltypes
    ) %>%
    dplyr::mutate(
      Condition = factor(Condition_clean, levels = celltype_order)
    ) %>%
    dplyr::select(-Condition_clean)
  
  p_sub
}

dotplots_v2 <- list()

for (i in seq_len(nrow(gene_map_v2))) {
  
  gene_symbol_use <- gene_map_v2$gene_symbol[i]
  gene_id_clean_use <- gene_map_v2$geneId_clean[i]
  
  message("Generating Figure2CDE relative plot for ", gene_symbol_use)
  
  dotplots_v2[[gene_symbol_use]] <- dotplotPromoters(
    result = result,
    geneId = gene_id_clean_use,
    geneName = gene_symbol_use,
    filterInternal = filterInternal_v2
  )
}

p_rel_list_v2 <- lapply(target_genes_v2, function(gene_use) {
  dotplots_v2[[gene_use]][[2]] %>%
    subset_dotplot_to_rescreened_celltypes_v2(gene_use) %>%
    decorate_plot_v2(gene_use)
})

names(p_rel_list_v2) <- target_genes_v2

row_plots_v2 <- arrangeGrob(
  grobs = p_rel_list_v2,
  nrow = 1,
  widths = rep(1, length(p_rel_list_v2))
)

final_relative_plot_v2 <- arrangeGrob(
  row_plots_v2,
  top = textGrob(
    "Relative promoter activity",
    gp = gpar(fontsize = FS$bigtitle, fontface = "bold")
  )
)

grid.newpage()
grid.draw(final_relative_plot_v2)

tiff(
  filename = file.path(out_dir_v2, "relative_promoter_activity_figure2CDE_v2.tiff"),
  width = 36,
  height = 20,
  units = "in",
  res = 600,
  compression = "lzw"
)
grid.draw(final_relative_plot_v2)
dev.off()

tiff(
  filename = "relative_promoter_activity_figure2CDE_v2.tiff",
  width = 36,
  height = 20,
  units = "in",
  res = 600,
  compression = "lzw"
)
grid.draw(final_relative_plot_v2)
dev.off()


#figure 5A-E
#using dot plot and heatmap for two style
#but the final version use heatmap
retina_combined <- readRDS("iso_sr_merged_retina_FL2.RDS")
retnet_list <- read.table("RetNet_category.txt", header = TRUE, sep = "\t")
transcript_annotation <- read.csv(
  "RetinaMerge.annotated.info.csv",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

transcript_status <- transcript_annotation %>%
  filter(category %in% c(
    "full-splice_match",
    "incomplete-splice_match",
    "novel_not_in_catalog",
    "novel_in_catalog"
  ))
#FL2 + CAGE support read status
pb_infor_FL2_cage <- pb_infor_FL2[pb_infor_FL2$within_cage_peak == TRUE, ]
nrow(pb_infor_FL2_cage)
table(pb_infor_FL2_cage$structural_category)
iso_fl2_classification <- pb_infor_FL2_cage
length(iso_fl2_classification)
head(iso_fl2_classification)
fl2_isoforms <- unique(iso_fl2_classification$isoform)
transcript_status_FL2 <- transcript_status %>%
  filter(pbid %in% fl2_isoforms)

#FL2 + CAGE + polyA motif support read status
pb_infor_FL2_cage_motif <- pb_infor_FL2[
  pb_infor_FL2$within_cage_peak == TRUE & !is.na(pb_infor_FL2$polyA_motif),
]
nrow(pb_infor_FL2_cage_motif)
table(pb_infor_FL2_cage_motif$structural_category)
iso_fl2_cage_motif_classification <- pb_infor_FL2_cage_motif
length(iso_fl2_cage_motif_classification)
head(iso_fl2_cage_motif_classification)
fl2_cage_motif_isoforms <- unique(iso_fl2_cage_motif_classification$isoform)
transcript_status_FL2_cage_motif <- transcript_status %>%
  filter(pbid %in% fl2_cage_motif_isoforms)

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

#find the corresponding pbid for FL2 + CAGE + polyA motif
retnet_pbid_FL2_cage_motif <- transcript_status_FL2_cage_motif %>%
  filter(gene %in% retnet_gene_unique$retnet_gene_unique) %>%
  select(gene, pbid) %>%
  distinct()

head(retnet_pbid_FL2_cage_motif)

mapped_retnet_FL2_cage_motif <- unique(retnet_pbid_FL2_cage_motif$gene)
num_mapped_retnet_FL2_cage_motif <- length(mapped_retnet_FL2_cage_motif)
print(num_mapped_retnet_FL2_cage_motif)

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

selected_pbid_FL2_cage_motif <- as.vector(retnet_pbid_FL2_cage_motif$pbid)
valid_pbid_FL2_cage_motif <- selected_pbid_FL2_cage_motif[
  selected_pbid_FL2_cage_motif %in% rownames(iso_counts)
]
sub_iso_matrix_FL2_cage_motif <- iso_counts[valid_pbid_FL2_cage_motif, sorted_barcodes]

rna_counts <- retina_combined@assays$RNA@counts

selected_gene_FL2 <- unique(retnet_pbid_FL2$gene)
valid_gene_FL2 <- selected_gene_FL2[selected_gene_FL2 %in% rownames(rna_counts)]
sub_rna_matrix_FL2 <- rna_counts[valid_gene_FL2, sorted_barcodes]
regular_rna_matrix_FL2 <- as.matrix(sub_rna_matrix_FL2)

selected_gene_FL2_cage_motif <- unique(retnet_pbid_FL2_cage_motif$gene)
valid_gene_FL2_cage_motif <- selected_gene_FL2_cage_motif[
  selected_gene_FL2_cage_motif %in% rownames(rna_counts)
]
sub_rna_matrix_FL2_cage_motif <- rna_counts[valid_gene_FL2_cage_motif, sorted_barcodes]
regular_rna_matrix_FL2_cage_motif <- as.matrix(sub_rna_matrix_FL2_cage_motif)

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
  
  message("novel_gene_df nrow: ", nrow(novel_gene_df))
  message("known_gene_df nrow: ", nrow(known_gene_df))
  message("all_expressed_novel length: ", length(all_expressed_novel))
  message("all_expressed_known length: ", length(all_expressed_known))
  
  
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

#new function for FL2 + CAGE + polyA motif dotplot and heatmap
make_disease_isoform_plots_FL2_cage_motif <- function(
    disease_col,
    retnet_list,
    transcript_status_FL2_cage_motif,
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
  
  #build disease-specific isoform table (FL2 + CAGE + polyA motif filtered)
  disease_pbid_FL2_cage_motif <- transcript_status_FL2_cage_motif %>%
    filter(gene %in% disease_genes$gene) %>% 
    mutate(status = case_when(
      category %in% c("full-splice_match", "incomplete-splice_match") ~ "known",
      category %in% c("novel_in_catalog", "novel_not_in_catalog") ~ "novel",
      TRUE ~ "other"
    )) %>%
    select(gene, pbid, status) %>%
    distinct() %>%
    filter(status %in% c("known", "novel"))
  
  if (nrow(disease_pbid_FL2_cage_motif) == 0) {
    warning("No FL2 + CAGE + polyA motif isoforms found for disease ", disease_col)
    return(NULL)
  }
  
  disease_pbid_FL2_cage_motif$transcript_status <- paste(
    disease_pbid_FL2_cage_motif$gene,
    disease_pbid_FL2_cage_motif$status,
    sep = "_"
  )
  disease_pbid_FL2_cage_motif$gene <- NULL
  disease_pbid_FL2_cage_motif$status <- NULL
  
  message("  - unique transcript_status: ",
          length(unique(disease_pbid_FL2_cage_motif$transcript_status)))
  
  #extract isoform count matrix (PBID-based)
  iso_counts2 <- iso_counts
  #make sure rownames are PBID only
  rownames(iso_counts2) <- sub(":.*$", "", rownames(iso_counts2))
  
  selected_pbid <- as.vector(disease_pbid_FL2_cage_motif$pbid)
  valid_pbid <- intersect(selected_pbid, rownames(iso_counts2))
  
  if (length(valid_pbid) == 0) {
    warning("No matching PBID in iso_counts for disease ", disease_col)
    return(NULL)
  }
  
  sub_iso_matrix <- iso_counts2[valid_pbid, , drop = FALSE]
  
  sub_iso_df <- as.data.frame(sub_iso_matrix)
  sub_iso_df$pbid <- rownames(sub_iso_df)
  
  sub_iso_df <- left_join(sub_iso_df, disease_pbid_FL2_cage_motif, by = "pbid")
  
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
  
  #prepare heatmap matrix
  avg_expr <- AverageExpression(
    seu,
    assays = "RNA",
    group.by = "celltype",
    slot = "data"
  )$RNA
  
  all_novel_heatmap_mat <- avg_expr[
    intersect(novel_gene_df$gene, rownames(avg_expr)),
    ,
    drop = FALSE
  ]
  
  all_known_heatmap_mat <- avg_expr[
    intersect(known_gene_df$gene, rownames(avg_expr)),
    ,
    drop = FALSE
  ]
  
  novel_heatmap_mat <- avg_expr[
    intersect(all_expressed_novel, rownames(avg_expr)),
    ,
    drop = FALSE
  ]
  
  known_heatmap_mat <- avg_expr[
    intersect(all_expressed_known, rownames(avg_expr)),
    ,
    drop = FALSE
  ]
  
  scale_by_row <- function(mat) {
    if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
      return(mat)
    }
    mat_scaled <- t(scale(t(as.matrix(mat))))
    mat_scaled[is.na(mat_scaled)] <- 0
    return(mat_scaled)
  }
  
  all_novel_heatmap_mat_scaled <- scale_by_row(all_novel_heatmap_mat)
  all_known_heatmap_mat_scaled <- scale_by_row(all_known_heatmap_mat)
  novel_heatmap_mat_scaled <- scale_by_row(novel_heatmap_mat)
  known_heatmap_mat_scaled <- scale_by_row(known_heatmap_mat)
  
  make_heatmap_plot <- function(
    heatmap_mat,
    plot_title
  ) {
    if (is.null(heatmap_mat) || nrow(heatmap_mat) == 0 || ncol(heatmap_mat) == 0) {
      warning("Empty heatmap matrix for ", plot_title)
      return(NULL)
    }
    
    heatmap_df <- as.data.frame(heatmap_mat)
    heatmap_df$feature <- rownames(heatmap_df)
    
    heatmap_long <- heatmap_df %>%
      pivot_longer(
        cols = -feature,
        names_to = "celltype",
        values_to = "scaled_expression"
      )
    
    p <- ggplot(heatmap_long, aes(x = celltype, y = feature, fill = scaled_expression)) +
      geom_tile() +
      ggtitle(plot_title) +
      theme_cowplot() +
      theme(
        legend.position = "right",
        legend.justification = "center",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold")
      ) +
      scale_fill_gradient2(
        low = "#330066",
        mid = "white",
        high = "#FFCC33",
        midpoint = 0
      )
    
    return(p)
  }
  
  novel_heatmap <- make_heatmap_plot(
    novel_heatmap_mat_scaled,
    paste0(disease_col, " novel heatmap (filtered)")
  )
  
  known_heatmap <- make_heatmap_plot(
    known_heatmap_mat_scaled,
    paste0(disease_col, " known heatmap (filtered)")
  )
  
  all_novel_heatmap <- make_heatmap_plot(
    all_novel_heatmap_mat_scaled,
    paste0(disease_col, " all novel heatmap")
  )
  
  all_known_heatmap <- make_heatmap_plot(
    all_known_heatmap_mat_scaled,
    paste0(disease_col, " all known heatmap")
  )
  
  #return infor
  invisible(list(
    disease = disease_col,
    seurat = seu,
    plots = list(
      novel_filtered = novel_dot,
      known_filtered = known_dot,
      novel_all = all_novel_dot,
      known_all = all_known_dot,
      novel_heatmap_filtered = novel_heatmap,
      known_heatmap_filtered = known_heatmap,
      novel_heatmap_all = all_novel_heatmap,
      known_heatmap_all = all_known_heatmap
    ),
    matrices = list(
      final_sorted_matrix = final_sorted_matrix,
      avg_expr = avg_expr,
      all_novel_heatmap_mat = all_novel_heatmap_mat,
      all_known_heatmap_mat = all_known_heatmap_mat,
      novel_heatmap_mat = novel_heatmap_mat,
      known_heatmap_mat = known_heatmap_mat
    ),
    data_frames = list(
      disease_pbid_FL2_cage_motif = disease_pbid_FL2_cage_motif,
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

#run the new function using each disease of interest for FL2 + CAGE + polyA motif
res_RP_polyA <- make_disease_isoform_plots_FL2_cage_motif(
  "RP", retnet_list,
  transcript_status_FL2_cage_motif,
  iso_counts,
  barcode_data,
  sorted_barcodes
)
res_CORD_polyA <- make_disease_isoform_plots_FL2_cage_motif(
  "CORD", retnet_list,
  transcript_status_FL2_cage_motif,
  iso_counts,
  barcode_data,
  sorted_barcodes
)
res_MD_polyA <- make_disease_isoform_plots_FL2_cage_motif(
  "MD", retnet_list,
  transcript_status_FL2_cage_motif,
  iso_counts,
  barcode_data,
  sorted_barcodes
)
res_CSNB_polyA <- make_disease_isoform_plots_FL2_cage_motif(
  "CSNB", retnet_list,
  transcript_status_FL2_cage_motif,
  iso_counts,
  barcode_data,
  sorted_barcodes
)
res_LCA_polyA <- make_disease_isoform_plots_FL2_cage_motif(
  "LCA", retnet_list,
  transcript_status_FL2_cage_motif,
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

dot_theme_big_axis_small_legend <- theme(
  axis.title.x = element_text(size = 28, face = "bold"),
  axis.title.y = element_text(size = 28, face = "bold"),
  axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
  axis.text.y  = element_text(size = 18),
  legend.position = "right",
  legend.justification = "center",
  legend.title = element_text(size = 18),
  legend.text  = element_text(size = 16),
  plot.title = element_text(size = 26, hjust = 0.5, face = "bold")
)

RP_novel_all_plot  <- RP_novel_all_plot  + dot_theme_big_axis_small_legend
CORD_novel_all_plot <- CORD_novel_all_plot + dot_theme_big_axis_small_legend
MD_novel_all_plot   <- MD_novel_all_plot   + dot_theme_big_axis_small_legend
CSNB_novel_all_plot <- CSNB_novel_all_plot + dot_theme_big_axis_small_legend
LCA_novel_all_plot  <- LCA_novel_all_plot  + dot_theme_big_axis_small_legend

RP_novel_all_plot_polyA   <- res_RP_polyA$plots$novel_all
CORD_novel_all_plot_polyA <- res_CORD_polyA$plots$novel_all
MD_novel_all_plot_polyA   <- res_MD_polyA$plots$novel_all
CSNB_novel_all_plot_polyA <- res_CSNB_polyA$plots$novel_all
LCA_novel_all_plot_polyA  <- res_LCA_polyA$plots$novel_all

RP_novel_all_plot_polyA   <- RP_novel_all_plot_polyA   + dot_theme_big_axis_small_legend
CORD_novel_all_plot_polyA <- CORD_novel_all_plot_polyA + dot_theme_big_axis_small_legend
MD_novel_all_plot_polyA   <- MD_novel_all_plot_polyA   + dot_theme_big_axis_small_legend
CSNB_novel_all_plot_polyA <- CSNB_novel_all_plot_polyA + dot_theme_big_axis_small_legend
LCA_novel_all_plot_polyA  <- LCA_novel_all_plot_polyA  + dot_theme_big_axis_small_legend

#save the FL2 + CAGE + polyA motif heatmap
RP_novel_all_heatmap_polyA   <- res_RP_polyA$plots$novel_heatmap_all
CORD_novel_all_heatmap_polyA <- res_CORD_polyA$plots$novel_heatmap_all
MD_novel_all_heatmap_polyA   <- res_MD_polyA$plots$novel_heatmap_all
CSNB_novel_all_heatmap_polyA <- res_CSNB_polyA$plots$novel_heatmap_all
LCA_novel_all_heatmap_polyA  <- res_LCA_polyA$plots$novel_heatmap_all

heatmap_theme_big_axis_small_legend <- theme(
  axis.title.x = element_text(size = 28, face = "bold"),
  axis.title.y = element_text(size = 28, face = "bold"),
  axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
  axis.text.y  = element_text(size = 18),
  legend.position = "right",
  legend.justification = "center",
  legend.title = element_text(size = 18),
  legend.text  = element_text(size = 16),
  plot.title = element_text(size = 26, hjust = 0.5, face = "bold")
)

heatmap_color_clean <- scale_fill_gradient2(
  low = "#4292C6",
  mid = "white",
  high = "#D95F5F",
  midpoint = 0.6
)

if (!is.null(RP_novel_all_heatmap_polyA)) {
  RP_novel_all_heatmap_polyA <- RP_novel_all_heatmap_polyA +
    heatmap_color_clean +
    heatmap_theme_big_axis_small_legend
}
if (!is.null(CORD_novel_all_heatmap_polyA)) {
  CORD_novel_all_heatmap_polyA <- CORD_novel_all_heatmap_polyA +
    heatmap_color_clean +
    heatmap_theme_big_axis_small_legend
}
if (!is.null(MD_novel_all_heatmap_polyA)) {
  MD_novel_all_heatmap_polyA <- MD_novel_all_heatmap_polyA +
    heatmap_color_clean +
    heatmap_theme_big_axis_small_legend
}
if (!is.null(CSNB_novel_all_heatmap_polyA)) {
  CSNB_novel_all_heatmap_polyA <- CSNB_novel_all_heatmap_polyA +
    heatmap_color_clean +
    heatmap_theme_big_axis_small_legend
}
if (!is.null(LCA_novel_all_heatmap_polyA)) {
  LCA_novel_all_heatmap_polyA <- LCA_novel_all_heatmap_polyA +
    heatmap_color_clean +
    heatmap_theme_big_axis_small_legend
}

#figure 5 generation A-E
#collect the disease results
disease_results_polyA <- list(
  RP    = res_RP_polyA,
  CORD  = res_CORD_polyA,
  MD    = res_MD_polyA,
  CSNB  = res_CSNB_polyA,
  LCA   = res_LCA_polyA
)

#clean the plot title and remove the celltype label
clean_subplot <- function(p, disease_name) {
  if (is.null(p)) return(NULL)
  
  p +
    labs(
      title = disease_name,
      x = NULL,
      y = NULL
    ) +
    theme(
      plot.title = element_text(
        size = 30,
        hjust = 0.5,
        face = "bold"
      ),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(
        size = 16,
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(size = 16),
      legend.position = "right",
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16)
    )
}

#keep the original heatmap color
heatmap_color_clean <- scale_fill_gradient2(
  low = "#4292C6",
  mid = "white",
  high = "#D95F5F",
  midpoint = 0.6
)

#clean the heatmap title and remove the celltype label
clean_heatmap_subplot <- function(p, disease_name) {
  if (is.null(p)) return(NULL)
  
  p +
    labs(
      title = disease_name,
      x = NULL,
      y = NULL
    ) +
    theme(
      plot.title = element_text(
        size = 30,
        hjust = 0.5,
        face = "bold"
      ),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(
        size = 16,
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(size = 16),
      legend.position = "right",
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16)
    ) +
    heatmap_color_clean
}

#make the revised all novel heatmap
RP_all_novel_heatmap    <- clean_heatmap_subplot(res_RP_polyA$plots$novel_heatmap_all, "RP")
CORD_all_novel_heatmap  <- clean_heatmap_subplot(res_CORD_polyA$plots$novel_heatmap_all, "CORD")
MD_all_novel_heatmap    <- clean_heatmap_subplot(res_MD_polyA$plots$novel_heatmap_all, "MD")
CSNB_all_novel_heatmap  <- clean_heatmap_subplot(res_CSNB_polyA$plots$novel_heatmap_all, "CSNB")
LCA_all_novel_heatmap   <- clean_heatmap_subplot(res_LCA_polyA$plots$novel_heatmap_all, "LCA")

#combine the revised all novel heatmap
p_all_novel_heatmap_revision <- 
  (CORD_all_novel_heatmap / MD_all_novel_heatmap) |
  (CSNB_all_novel_heatmap / LCA_all_novel_heatmap) |
  RP_all_novel_heatmap

#save the revised all novel heatmap
ggsave(
  filename = "Figure5_all_novel_heatmap_FL2_cage_polyA_revised.tiff",
  plot     = p_all_novel_heatmap_revision,
  width    = 30,
  height   = 24,
  dpi      = 600,
  device   = "tiff",
  compression = "lzw",
  bg       = "WHITE"
)

ggsave(
  filename = "Figure5_all_novel_heatmap_FL2_cage_polyA_revised.pdf",
  plot     = p_all_novel_heatmap_revision,
  width    = 30,
  height   = 24,
  bg       = "WHITE"
)

#make the revised all novel dotplot
RP_novel_all_plot_polyA    <- clean_subplot(res_RP_polyA$plots$novel_all, "RP")
CORD_novel_all_plot_polyA  <- clean_subplot(res_CORD_polyA$plots$novel_all, "CORD")
MD_novel_all_plot_polyA    <- clean_subplot(res_MD_polyA$plots$novel_all, "MD")
CSNB_novel_all_plot_polyA  <- clean_subplot(res_CSNB_polyA$plots$novel_all, "CSNB")
LCA_novel_all_plot_polyA   <- clean_subplot(res_LCA_polyA$plots$novel_all, "LCA")

#combine the revised all novel dotplot
p_all_novel_dotplot_revision <- 
  (CORD_novel_all_plot_polyA / MD_novel_all_plot_polyA) |
  (CSNB_novel_all_plot_polyA / LCA_novel_all_plot_polyA) |
  RP_novel_all_plot_polyA

#save the revised all novel dotplot
ggsave(
  filename = "Figure5_all_novel_dotplot_FL2_cage_polyA_revised.tiff",
  plot     = p_all_novel_dotplot_revision,
  width    = 30,
  height   = 24,
  dpi      = 600,
  device   = "tiff",
  compression = "lzw",
  bg       = "WHITE"
)

ggsave(
  filename = "Figure5_all_novel_dotplot_FL2_cage_polyA_revised.pdf",
  plot     = p_all_novel_dotplot_revision,
  width    = 30,
  height   = 24,
  bg       = "WHITE"
)

#make the PBID export table for each disease
make_pbid_export_table <- function(
    res_obj,
    disease_name,
    barcode_data
) {
  
  if (is.null(res_obj)) return(NULL)
  
  pbid_df <- res_obj$data_frames$disease_pbid_FL2_cage_motif
  sub_iso_df <- res_obj$data_frames$sub_iso_df
  
  if (is.null(pbid_df) || nrow(pbid_df) == 0) return(NULL)
  if (is.null(sub_iso_df) || nrow(sub_iso_df) == 0) return(NULL)
  
  #prepare the basic PBID information
  pbid_export <- pbid_df %>%
    mutate(
      disease = disease_name,
      gene = sub("_(known|novel)$", "", transcript_status),
      status = sub("^.*_(known|novel)$", "\\1", transcript_status),
      gene_status = paste(gene, status, sep = "-")
    ) %>%
    select(
      disease,
      gene_status,
      gene,
      status,
      pbid,
      transcript_status
    ) %>%
    distinct()
  
  #find the barcode columns
  barcode_cols <- intersect(barcode_data$barcode, colnames(sub_iso_df))
  
  if (length(barcode_cols) == 0) {
    stop("No barcode columns were found in sub_iso_df.")
  }
  
  #prepare the long count table
  pbid_count_long <- sub_iso_df %>%
    select(pbid, transcript_status, all_of(barcode_cols)) %>%
    pivot_longer(
      cols = all_of(barcode_cols),
      names_to = "barcode",
      values_to = "count"
    )
  
  #sum the PBID counts by cell type
  pbid_celltype_count <- pbid_count_long %>%
    left_join(
      barcode_data %>%
        select(
          barcode,
          celltype = cell.type.integrated
        ),
      by = "barcode"
    ) %>%
    filter(!is.na(celltype)) %>%
    group_by(pbid, transcript_status, celltype) %>%
    summarise(
      count = sum(count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      celltype_count_col = paste0("celltype_", celltype, "_count")
    ) %>%
    select(pbid, transcript_status, celltype_count_col, count) %>%
    pivot_wider(
      names_from = celltype_count_col,
      values_from = count,
      values_fill = 0
    )
  
  #calculate the total count across all cell types
  celltype_count_cols <- grep("^celltype_.*_count$", colnames(pbid_celltype_count), value = TRUE)
  
  pbid_celltype_count <- pbid_celltype_count %>%
    mutate(
      total_count = rowSums(across(all_of(celltype_count_cols)), na.rm = TRUE)
    )
  
  #merge the basic PBID information with cell type counts
  pbid_export <- pbid_export %>%
    left_join(
      pbid_celltype_count,
      by = c("pbid", "transcript_status")
    ) %>%
    relocate(
      disease,
      gene_status,
      gene,
      status,
      pbid,
      transcript_status,
      total_count
    )
  
  return(pbid_export)
}

#extract the PBID information from all diseases
pbid_export_list <- lapply(
  names(disease_results_polyA),
  function(disease_name) {
    make_pbid_export_table(
      res_obj = disease_results_polyA[[disease_name]],
      disease_name = disease_name,
      barcode_data = barcode_data
    )
  }
)

names(pbid_export_list) <- names(disease_results_polyA)

#remove empty disease results
pbid_export_list <- pbid_export_list[!sapply(pbid_export_list, is.null)]

#make the excel workbook
wb <- createWorkbook()

for (sheet_name in names(pbid_export_list)) {
  
  df <- pbid_export_list[[sheet_name]]
  
  #clean the sheet name for excel
  clean_sheet_name <- substr(gsub("[\\[\\]\\*\\?/\\\\:]", "_", sheet_name), 1, 31)
  
  addWorksheet(wb, clean_sheet_name)
  writeData(wb, clean_sheet_name, df)
  
  freezePane(wb, clean_sheet_name, firstRow = TRUE)
  addFilter(wb, clean_sheet_name, row = 1, cols = seq_len(ncol(df)))
  setColWidths(wb, clean_sheet_name, cols = seq_len(ncol(df)), widths = "auto")
  
  header_style <- createStyle(
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    border = "Bottom"
  )
  
  addStyle(
    wb,
    sheet = clean_sheet_name,
    style = header_style,
    rows = 1,
    cols = seq_len(ncol(df)),
    gridExpand = TRUE
  )
}

#save the excel file
saveWorkbook(
  wb,
  file = "Figure5_FL2_cage_polyA_PBID_information_by_disease.xlsx",
  overwrite = TRUE
)

#save the combined PBID table
pbid_export_all_diseases <- bind_rows(pbid_export_list)

write.csv(
  pbid_export_all_diseases,
  file = "Figure5_FL2_cage_polyA_PBID_information_all_diseases.csv",
  row.names = FALSE
)


#ggtranscript to plot the structure of transcripts
#read the annotation file
#read the file for filtering
pb_anno <- rtracklayer::import("RetinaMerge.sorted.filtered_lite.FL2_CAGE_polyA.gff")
class(pb_anno)
pb_anno <- pb_anno %>% dplyr::as_tibble()
class(pb_anno)
junction_anno <- read.table("RetinaMerge_classification.filtered_lite_junctions.high_confident.txt", header = TRUE)
reference_anno <- rtracklayer::import("GRCh38p14/GRCh38.p14.v44/gencode.v44.chr_patch_hapl_scaff.annotation.gff3")
class(reference_anno)
reference_anno <- reference_anno %>% dplyr::as_tibble()
class(reference_anno)
#read the file for filtering
pb_infor <- read.csv("RetinaMerge_classification.filtered_lite_classification.high_confident.txt", header = TRUE, sep = "\t")
#previously only select for PBids within cage peak, then I found can filter out those NA polyA motifs, so change the code but retain the name as it is used in the following codes to reduce work
#after filter with both within cage peak and with non-NA polyA motif, the real name for this object should be pbid_cage_polyA or pbid_FL
#pbid_cage <- pb_infor %>% filter(within_cage_peak == TRUE)
#FL
pbid_fl <- pb_infor %>%
  filter(within_cage_peak == TRUE)
pbid_cage <- pb_infor %>% filter(within_cage_peak == TRUE)
#read the rds file
retina_combined <- readRDS("D:/bioinformatic/isoseq/processed/from_Ray/revision/redo_FL2CAGEpolyA/ggtranscript/iso_sr_merged_retina.RDS")
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
  filename   = "ABCA4_ggtranscript_FL2_CAGE_polyA_with_bar.tiff",
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


#figure 6B
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
  filename   = "figure6B_ABCA4_ggtranscript_FL2_CAGE_polyA_with_bar.tiff",
  plot       = ABCA4_with_bar,
  device     = "tiff",
  width      = 16,
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

#supplementary figure 8
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
  filename   = "figureS7_ABCA4_ggtranscript_FL2_CAGE_polyA_withbar_diffexon.tiff",
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

#figure 6D
ggsave(
  filename   = "figure6D_ABCA4_Junctions.tiff",
  plot       = ABCA4_junction_plot,
  width      = 26,
  height     = 8,
  dpi        = 600,
  bg         = "WHITE",
  compression = "lzw"
)

#ASprofile
#run the ASprofile pipeline in cmd
#plot the result
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

#figure 3A
ggsave(
  "ASevent_figure3A.tiff",
  plot  = AS_PBid,
  device = "tiff",
  type   = "cairo",
  width = 8, height = 10, units = "in",
  dpi   = 600,
  compression = "lzw",
  bg = "WHITE"
)

#supplementary figure 5A
#use the FL2 subset to retain CAGE and polyA motif information
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


#validation with a bulk retina dataset and RJunBase
#SJ match ratio using RJunBase
df <- read_csv("SJ_loose_match_summary.csv")
print(df)
df_1bp <- df %>% filter(tol_bp == 1)
df_plot <- tibble(
  Category = c("FSM / ISM (known)", "NIC / NNC (novel)"),
  MatchRate = c(df_1bp$FSM_ISM_rate, df_1bp$NNC_NIC_rate)
)
#plot the SJ match ratio
#supplementary figure 2D
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
#supplmentary figure 2D
ggsave("SJ_matchrate_barplot.png", p_SJ, width = 8, height = 6, dpi = 600, bg = 'WHITE')

#supplementary figure 2E
#short read support
df <- read_csv("short_read_support/SJ_shortread_support_summary.csv")
print(df)
#keep HighConfident strict support only
df_hc <- df %>%
  filter(set_name == "HighConfident") %>%
  mutate(
    Category = case_when(
      group == "FSM_ISM" ~ "FSM / ISM (known)",
      group == "NNC_NIC" ~ "NIC / NNC (novel)",
      TRUE ~ group
    ),
    MatchRate = support_rate_strict
  ) %>%
  select(Category, MatchRate, n_junctions, n_supported_strict)

#keep category order consistent
df_hc <- df_hc %>%
  mutate(
    Category = factor(
      Category,
      levels = c("FSM / ISM (known)", "NIC / NNC (novel)")
    )
  )

print(df_hc)

#plot the HighConfident SJ strict support ratio
p_SJ_HC <- ggplot(df_hc, aes(x = Category, y = MatchRate, fill = Category)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(
    aes(label = percent(MatchRate, accuracy = 0.1)),
    vjust = -0.5,
    size = 5
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.05)
  ) +
  labs(
    title = "Validation of high-confidence splice junctions",
    subtitle = "High-confidence Iso-Seq junctions validated by short-read support",
    x = NULL,
    y = "Supported junctions (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

print(p_SJ_HC)

#save the plot
ggsave(
  "HighConfident_SJ_strict_support_barplot.png",
  p_SJ_HC,
  width = 8,
  height = 6,
  dpi = 600,
  bg = "WHITE"
)

#supplementary figure 2B
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
#supplementary figure 2B
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
#supplementary figure 2C
ggsave("SLC6A11_nosplit_rescaled_identical.tiff", plot = SLC6A11_rescaled_selected_nosplit, width = 12, height = 6, dpi = 1200, limitsize = FALSE, bg = "white")

#save supplementary tables
#paths
pb_infor_path <- "original_from_Ray/ray pigeon output/classify.retina.merge.GRCh38.p14.with.cage.polya.FL/RetinaMerge_classification.filtered_lite_classification.txt"
pb_infor_FL2_path <- "RetinaMerge_classification.filtered_lite_FL2_classification.txt"
scRetina_path <- "iso_sr_merged_retina.RDS"

#load data
pb_infor <- read.csv(pb_infor_path, header = TRUE, sep = "\t")
pb_infor_FL2 <- read.csv(pb_infor_FL2_path, header = TRUE, sep = "\t")
scRetina <- readRDS(scRetina_path)

#count matrix
iso.mat <- GetAssayData(scRetina, assay = "isoform", slot = "counts")
rownames(iso.mat) <- sub(":.*$", "", rownames(iso.mat))

#metadata
meta <- scRetina@meta.data
meta <- meta[colnames(iso.mat), ]

#cell type counts
ct <- meta$cell.type.integrated
H_ct <- sparse.model.matrix(~ 0 + ct)
colnames(H_ct) <- sub("^ct", "", colnames(H_ct))

M_celltype <- iso.mat %*% H_ct

#sample counts
sid <- meta$sample.id
H_sid <- sparse.model.matrix(~ 0 + sid)
colnames(H_sid) <- sub("^sid", "", colnames(H_sid))

M_sample <- iso.mat %*% H_sid

#combined count table
M_all <- cbind(M_celltype, M_sample)

pb_count <- as.data.frame.matrix(M_all)
pb_count$PBid <- rownames(M_all)

#annotation table
anno_cols <- c(
  "isoform",
  "chrom", "strand", "length",
  "structural_category", "subcategory",
  "associated_gene",
  "within_cage_peak", "polyA_motif"
)

pb_anno <- pb_infor[, anno_cols]
colnames(pb_anno)[1] <- "PBid"

pb_anno_FL2 <- pb_infor_FL2[, anno_cols]
colnames(pb_anno_FL2)[1] <- "PBid"

#merge tables
matrix_all <- merge(pb_anno, pb_count, by = "PBid", all = FALSE)
matrix_FL2 <- merge(pb_anno_FL2, pb_count, by = "PBid", all = FALSE)

write.csv(matrix_all, "PBid_matrix_all.csv", row.names = FALSE)
write.csv(matrix_FL2, "PBid_matrix_FL2.csv", row.names = FALSE)

#donor support
sample_cols <- c("Retina2", "Retina3", "Retina4")

matrix_FL2_2sample <- matrix_FL2[
  rowSums(matrix_FL2[, sample_cols] > 0) >= 2,
]

matrix_all_2sample <- matrix_all[
  rowSums(matrix_all[, sample_cols] > 0) >= 2,
]

matrix_FL2_3sample <- matrix_FL2[
  rowSums(matrix_FL2[, sample_cols] > 0) >= 3,
]

nrow(matrix_FL2)
nrow(matrix_FL2_2sample)
table(matrix_FL2_2sample$structural_category)

nrow(matrix_all)
nrow(matrix_all_2sample)
table(matrix_all_2sample$structural_category)

nrow(matrix_FL2)
nrow(matrix_FL2_3sample)
table(matrix_FL2_3sample$structural_category)

#doublet rate
DefaultAssay(scRetina) <- "RNA"
sce <- as.SingleCellExperiment(scRetina, assay = "RNA")
sce <- scDblFinder(sce, samples = colData(sce)$sample.id)

scRetina$scDblFinder.class <- colData(sce)$scDblFinder.class
scRetina$scDblFinder.score <- colData(sce)$scDblFinder.score
scRetina$scDblFinder.weighted <- colData(sce)$scDblFinder.weighted

table(scRetina$scDblFinder.class)
prop.table(table(scRetina$scDblFinder.class))

#CAGE and polyA motif
matrix_FL2_cage <- matrix_FL2[
  matrix_FL2$within_cage_peak == TRUE,
]

matrix_FL2_cage_motif <- matrix_FL2[
  matrix_FL2$within_cage_peak == TRUE & !is.na(matrix_FL2$polyA_motif),
]

nrow(matrix_FL2_cage)
table(matrix_FL2_cage$structural_category)

with(matrix_FL2, table(within_cage_peak, is.na(polyA_motif)))

nrow(matrix_FL2_cage_motif)
table(matrix_FL2_cage_motif$structural_category)



#orfanage summary
#use orfanage to estimate CDS region of all isoforms detected in this study
#paths
orf_stats_file <- "RetinaMerge.filtered_lite.orfanage.stats.tsv"
class_file     <- "RetinaMerge_classification.filtered_lite_classification.txt"
out_dir        <- "orfanage_analysis"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#Helper functions

safe_num <- function(x) {
  x <- as.character(x)
  x[x %in% c("-", "", "NA", "NaN", "nan")] <- NA
  suppressWarnings(as.numeric(x))
}

save_plot_pdf_png <- function(p, file_base, width = 10, height = 7, dpi = 300) {
  ggsave(paste0(file_base, ".pdf"), p, width = width, height = height, useDingbats = FALSE)
  ggsave(paste0(file_base, ".png"), p, width = width, height = height, dpi = dpi)
}

#Read ORFanage stats (robust reader)

message("Reading ORFanage stats ...")

orf <- read.table(
  orf_stats_file,
  header = TRUE,
  sep = "",
  quote = "",
  fill = TRUE,
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = TRUE
)

message("Dimensions of ORFanage stats: ", paste(dim(orf), collapse = " x "))
message("Detected ORFanage columns:")
print(colnames(orf))

if (nrow(orf) == 0) {
  stop("ORFanage stats file was read as 0 rows.")
}

required_cols <- c("query_id", "template_id", "notes")
missing_required <- setdiff(required_cols, colnames(orf))
if (length(missing_required) > 0) {
  stop(
    "Missing required columns in ORFanage stats: ",
    paste(missing_required, collapse = ", "),
    "\nDetected columns are:\n",
    paste(colnames(orf), collapse = ", ")
  )
}

#numeric-like columns if present
num_cols <- c(
  "num_templates", "query_len", "template_len", "union_len", "pass",
  "len_match", "len_inframe", "len_outframe", "len_extra", "len_missing",
  "length_pi", "match_length_pi", "inframe_length_pi", "alignment_match",
  "start_match", "stop_match", "pi"
)
num_cols <- intersect(num_cols, colnames(orf))
for (cc in num_cols) {
  orf[[cc]] <- safe_num(orf[[cc]])
}

#force key text columns
orf$query_id    <- as.character(orf$query_id)
orf$template_id <- as.character(orf$template_id)
orf$notes       <- as.character(orf$notes)

if ("segment" %in% colnames(orf)) {
  orf$segment <- as.character(orf$segment)
} else {
  orf$segment <- NA_character_
}

#Choose the best ORFanage template row per query_id

message("Selecting best ORFanage row per query_id ...")

orf2 <- orf %>%
  mutate(
    has_template = !is.na(template_id) & template_id != "-",
    is_best_gtf  = !is.na(notes) & notes == "best_gtf",
    is_keep_cds  = !is.na(notes) & notes == "keep_cds",
    is_pass      = !is.na(pass) & pass == 1,
    
    pi                = ifelse(is.na(pi), -Inf, pi),
    match_length_pi   = ifelse(is.na(match_length_pi), -Inf, match_length_pi),
    inframe_length_pi = ifelse(is.na(inframe_length_pi), -Inf, inframe_length_pi),
    
    query_len      = ifelse(is.na(query_len), 0, query_len),
    template_len   = ifelse(is.na(template_len), 0, template_len),
    len_match      = ifelse(is.na(len_match), 0, len_match),
    len_inframe    = ifelse(is.na(len_inframe), 0, len_inframe),
    len_outframe   = ifelse(is.na(len_outframe), 0, len_outframe),
    len_extra      = ifelse(is.na(len_extra), 0, len_extra),
    len_missing    = ifelse(is.na(len_missing), 0, len_missing),
    alignment_match= ifelse(is.na(alignment_match), 0, alignment_match),
    start_match    = ifelse(is.na(start_match), 0, start_match),
    stop_match     = ifelse(is.na(stop_match), 0, stop_match),
    
    rank_best_gtf = ifelse(is_best_gtf, 1L, 0L),
    rank_keep_cds = ifelse(is_keep_cds, 1L, 0L),
    rank_pass     = ifelse(is_pass, 1L, 0L),
    rank_template = ifelse(has_template, 1L, 0L)
  )

best_orf <- orf2 %>%
  arrange(
    query_id,
    desc(rank_best_gtf),
    desc(rank_keep_cds),
    desc(rank_pass),
    desc(rank_template),
    desc(pi),
    desc(match_length_pi),
    desc(inframe_length_pi),
    desc(len_match),
    desc(len_inframe),
    desc(alignment_match),
    desc(start_match),
    desc(stop_match)
  ) %>%
  group_by(query_id) %>%
  slice(1) %>%
  ungroup()

write.csv(best_orf, file.path(out_dir, "ORFanage_best_template_per_query.csv"), row.names = FALSE)

#Read SQANTI classification
cls <- fread(class_file, sep = "\t", header = TRUE, data.table = FALSE)
cls$isoform <- as.character(cls$isoform)

wanted_cols <- c(
  "isoform", "structural_category", "associated_gene", "associated_transcript",
  "FL", "subcategory", "RTS_stage", "all_canonical", "coding", "predicted_NMD"
)
wanted_cols <- intersect(wanted_cols, colnames(cls))
cls2 <- cls[, wanted_cols, drop = FALSE]

if ("FL" %in% colnames(cls2)) {
  cls2$FL <- safe_num(cls2$FL)
}

#Merge SQANTI + ORFanage
df <- cls2 %>%
  left_join(best_orf, by = c("isoform" = "query_id"))

#Derive protein-level consequence features
df <- df %>%
  mutate(
    has_orf = !is.na(query_len) & query_len > 0,
    has_template = !is.na(template_id) & template_id != "-",
    
    query_template_ratio = ifelse(
      has_orf & !is.na(template_len) & template_len > 0,
      query_len / template_len,
      NA_real_
    ),
    
    frac_inframe = ifelse(has_orf & query_len > 0, len_inframe / query_len, NA_real_),
    frac_outframe = ifelse(has_orf & query_len > 0, len_outframe / query_len, NA_real_),
    frac_extra = ifelse(has_orf & query_len > 0, len_extra / query_len, NA_real_),
    frac_missing = ifelse(!is.na(template_len) & template_len > 0, len_missing / template_len, NA_real_)
  )

#ORF status
df <- df %>%
  mutate(
    orf_status = case_when(
      !has_orf ~ "No predicted ORF",
      has_orf & !has_template ~ "Predicted ORF, no matched template",
      has_orf & has_template ~ "Predicted ORF, matched template",
      TRUE ~ "Other"
    )
  )

#Protein size effect
df <- df %>%
  mutate(
    protein_size_effect = case_when(
      !has_orf ~ "No predicted ORF",
      has_orf & !has_template ~ "ORF without matched template",
      !is.na(query_template_ratio) & query_template_ratio < 0.50 ~ "Severely truncated (<50%)",
      !is.na(query_template_ratio) & query_template_ratio >= 0.50 & query_template_ratio < 0.90 ~ "Truncated (50-90%)",
      !is.na(query_template_ratio) & query_template_ratio >= 0.90 & query_template_ratio <= 1.10 ~ "Near reference length (90-110%)",
      !is.na(query_template_ratio) & query_template_ratio > 1.10 ~ "Extended (>110%)",
      TRUE ~ "Other"
    )
  )

#Frame effect
df <- df %>%
  mutate(
    frame_effect = case_when(
      !has_orf ~ "No predicted ORF",
      has_orf & !has_template ~ "ORF without matched template",
      frac_outframe > 0 ~ "Out-of-frame change present",
      frac_outframe == 0 & frac_inframe > 0 & (len_extra > 0 | len_missing > 0) ~ "In-frame alteration only",
      frac_outframe == 0 & len_extra == 0 & len_missing == 0 ~ "Reference-like CDS",
      TRUE ~ "Other"
    )
  )

#Simplified protein consequence for manuscript / response
df <- df %>%
  mutate(
    protein_consequence_simple = case_when(
      !has_orf ~ "No predicted ORF",
      has_orf & !has_template ~ "Novel ORF without matched template",
      frac_outframe > 0 ~ "Frameshift / out-of-frame alteration",
      !is.na(query_template_ratio) & query_template_ratio < 0.90 ~ "Predicted truncation",
      !is.na(query_template_ratio) & query_template_ratio > 1.10 ~ "Predicted extension",
      !is.na(query_template_ratio) &
        query_template_ratio >= 0.90 & query_template_ratio <= 1.10 &
        frac_outframe == 0 & (len_extra > 0 | len_missing > 0) ~ "In-frame altered CDS",
      !is.na(query_template_ratio) &
        query_template_ratio >= 0.90 & query_template_ratio <= 1.10 &
        frac_outframe == 0 & len_extra == 0 & len_missing == 0 ~ "Reference-like CDS",
      TRUE ~ "Other"
    )
  )

df <- df %>%
  mutate(
    likely_protein_altering = protein_consequence_simple %in% c(
      "Novel ORF without matched template",
      "Frameshift / out-of-frame alteration",
      "Predicted truncation",
      "Predicted extension",
      "In-frame altered CDS"
    )
  )

#Save full merged table
write.csv(df, file.path(out_dir, "SQANTI_ORFanage_merged_full.csv"), row.names = FALSE)

#Summary tables
summary_overall <- df %>%
  summarise(
    n_isoforms = n(),
    n_with_orf = sum(has_orf, na.rm = TRUE),
    frac_with_orf = mean(has_orf, na.rm = TRUE),
    n_with_template = sum(has_template, na.rm = TRUE),
    frac_with_template = mean(has_template, na.rm = TRUE),
    n_likely_protein_altering = sum(likely_protein_altering, na.rm = TRUE),
    frac_likely_protein_altering = mean(likely_protein_altering, na.rm = TRUE)
  )

write.csv(summary_overall, file.path(out_dir, "summary_overall.csv"), row.names = FALSE)

if ("structural_category" %in% colnames(df)) {
  summary_by_class <- df %>%
    group_by(structural_category) %>%
    summarise(
      n_isoforms = n(),
      n_with_orf = sum(has_orf, na.rm = TRUE),
      frac_with_orf = mean(has_orf, na.rm = TRUE),
      n_with_template = sum(has_template, na.rm = TRUE),
      frac_with_template = mean(has_template, na.rm = TRUE),
      n_likely_protein_altering = sum(likely_protein_altering, na.rm = TRUE),
      frac_likely_protein_altering = mean(likely_protein_altering, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_isoforms))
  
  write.csv(summary_by_class, file.path(out_dir, "summary_by_structural_category.csv"), row.names = FALSE)
}

summary_consequence <- df %>%
  count(protein_consequence_simple, sort = TRUE) %>%
  mutate(frac = n / sum(n))

write.csv(summary_consequence, file.path(out_dir, "summary_protein_consequence_overall.csv"), row.names = FALSE)

if ("structural_category" %in% colnames(df)) {
  summary_consequence_by_class <- df %>%
    count(structural_category, protein_consequence_simple) %>%
    group_by(structural_category) %>%
    mutate(frac_within_class = n / sum(n)) %>%
    ungroup()
  
  write.csv(summary_consequence_by_class, file.path(out_dir, "summary_protein_consequence_by_class.csv"), row.names = FALSE)
}

if ("associated_gene" %in% colnames(df)) {
  gene_candidates <- df %>%
    filter(!is.na(associated_gene), associated_gene != "") %>%
    group_by(associated_gene) %>%
    summarise(
      n_isoforms = n(),
      n_with_orf = sum(has_orf, na.rm = TRUE),
      n_likely_protein_altering = sum(likely_protein_altering, na.rm = TRUE),
      frac_likely_protein_altering = mean(likely_protein_altering, na.rm = TRUE),
      n_no_orf = sum(protein_consequence_simple == "No predicted ORF", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_likely_protein_altering), desc(frac_likely_protein_altering), desc(n_isoforms))
  
  write.csv(gene_candidates, file.path(out_dir, "gene_level_protein_altering_candidates.csv"), row.names = FALSE)
}


#FL>=2 subset analysis
df_fl2 <- df %>%
  filter(!is.na(FL) & FL >= 2)

#summary
summary_fl2 <- df_fl2 %>%
  summarise(
    n_isoforms = n(),
    n_with_orf = sum(has_orf, na.rm = TRUE),
    frac_with_orf = mean(has_orf, na.rm = TRUE),
    n_with_template = sum(has_template, na.rm = TRUE),
    frac_with_template = mean(has_template, na.rm = TRUE),
    n_likely_protein_altering = sum(likely_protein_altering, na.rm = TRUE),
    frac_likely_protein_altering = mean(likely_protein_altering, na.rm = TRUE)
  )

write.csv(summary_fl2,
          file.path(out_dir, "summary_FL2_overall.csv"),
          row.names = FALSE)

#consequence summary
summary_consequence_fl2 <- df_fl2 %>%
  count(protein_consequence_simple, sort = TRUE) %>%
  mutate(frac = n / sum(n))

write.csv(summary_consequence_fl2,
          file.path(out_dir, "summary_FL2_protein_consequence.csv"),
          row.names = FALSE)

#by structural category
summary_by_class_fl2 <- df_fl2 %>%
  group_by(structural_category) %>%
  summarise(
    n_isoforms = n(),
    n_with_orf = sum(has_orf, na.rm = TRUE),
    frac_with_orf = mean(has_orf, na.rm = TRUE),
    n_likely_protein_altering = sum(likely_protein_altering, na.rm = TRUE),
    frac_likely_protein_altering = mean(likely_protein_altering, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_by_class_fl2,
          file.path(out_dir, "summary_FL2_by_structural_category.csv"),
          row.names = FALSE)



#PBid matrix
pbid_matrix_file <- "PBid_matrix_allCT.csv"
pbm <- fread(pbid_matrix_file, data.table = FALSE)

colnames(pbm) <- make.names(colnames(pbm), unique = TRUE)

if (!("PBid" %in% colnames(pbm))) {
  stop("Column 'PBid' not found in PBid_matrix_all.csv")
}

pbm$PBid <- as.character(pbm$PBid)

meta_cols_pbm <- c(
  "PBid", "chrom", "strand", "length", "structural_category", "subcategory",
  "associated_gene", "within_cage_peak", "polyA_motif"
)

celltype_cols <- setdiff(colnames(pbm), meta_cols_pbm)

for (cc in celltype_cols) {
  pbm[[cc]] <- suppressWarnings(as.numeric(pbm[[cc]]))
}

#RetNet disease gene file
retnet_raw <- fread(
  retnet_file,
  sep = "\t",
  header = TRUE,
  fill = TRUE,
  quote = "\"",
  data.table = FALSE
)

gene_universe <- unique(c(
  as.character(pbm$associated_gene),
  if ("associated_gene" %in% colnames(df)) as.character(df$associated_gene) else character(0)
))
gene_universe <- gene_universe[!is.na(gene_universe) & gene_universe != ""]

retnet_text <- unlist(retnet_raw, use.names = FALSE)
retnet_text <- as.character(retnet_text)
retnet_text[is.na(retnet_text)] <- ""

retnet_tokens <- unlist(strsplit(gsub("[\"',;:/()]", " ", retnet_text), "\\s+"))
retnet_tokens <- retnet_tokens[nzchar(retnet_tokens)]

gene_universe_upper <- toupper(gene_universe)
retnet_tokens_upper <- unique(toupper(retnet_tokens))

disease_gene_set <- unique(gene_universe[gene_universe_upper %in% retnet_tokens_upper])
disease_gene_set <- sort(unique(disease_gene_set))

writeLines(disease_gene_set, file.path(out_dir, "RetNet_gene_set_detected_in_dataset.txt"))
message("Detected RetNet genes present in dataset: ", length(disease_gene_set))

#11.3 Merge PBid matrix into df

pbm_keep <- pbm[, c("PBid", "associated_gene", celltype_cols), drop = FALSE]

df2 <- df %>%
  left_join(pbm_keep, by = c("isoform" = "PBid"), suffix = c("", ".pbm"))

if ("associated_gene.pbm" %in% colnames(df2)) {
  df2 <- df2 %>%
    mutate(
      associated_gene_final = case_when(
        !is.na(associated_gene) & associated_gene != "" ~ associated_gene,
        !is.na(associated_gene.pbm) & associated_gene.pbm != "" ~ associated_gene.pbm,
        TRUE ~ NA_character_
      )
    )
} else {
  df2 <- df2 %>%
    mutate(associated_gene_final = associated_gene)
}

#Paths
pfam_domtblout_file <- "pfam_analysis_parallel/pfam.domtblout"
merged_orf_file     <- "orfanage_analysis/SQANTI_ORFanage_merged_full.csv"
retnet_gene_file    <- "orfanage_analysis/RetNet_gene_set_detected_in_dataset.txt"

out_dir_pfam <- "pfam_analysis_parallel/pfam_integrated_results_v2"
dir.create(out_dir_pfam, showWarnings = FALSE, recursive = TRUE)

safe_num <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "nan", "-", "Inf", "-Inf")] <- NA
  suppressWarnings(as.numeric(x))
}

save_plot_pdf_jpg <- function(p, file_base, width = 10, height = 7, dpi = 300) {
  ggsave(paste0(file_base, ".pdf"), p, width = width, height = height, useDingbats = FALSE)
  ggsave(paste0(file_base, ".jpg"), p, width = width, height = height, dpi = dpi)
}

mode_string <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  tb <- sort(table(x), decreasing = TRUE)
  names(tb)[1]
}

collapse_intervals_width <- function(start, end) {
  if (length(start) == 0 || length(end) == 0) return(0L)
  o <- order(start, end)
  start <- as.integer(start[o])
  end   <- as.integer(end[o])
  
  ms <- integer(0)
  me <- integer(0)
  
  cur_s <- start[1]
  cur_e <- end[1]
  
  if (length(start) > 1) {
    for (i in 2:length(start)) {
      if (start[i] <= cur_e + 1L) {
        cur_e <- max(cur_e, end[i])
      } else {
        ms <- c(ms, cur_s)
        me <- c(me, cur_e)
        cur_s <- start[i]
        cur_e <- end[i]
      }
    }
  }
  
  ms <- c(ms, cur_s)
  me <- c(me, cur_e)
  
  sum(me - ms + 1L)
}

#overlap fraction relative to shorter interval
interval_overlap_frac <- function(a1, a2, b1, b2) {
  ov <- max(0, min(a2, b2) - max(a1, b1) + 1)
  shorter <- min(a2 - a1 + 1, b2 - b1 + 1)
  if (shorter <= 0) return(0)
  ov / shorter
}

#keep best hit among strongly overlapping hits on same protein
collapse_overlapping_hits <- function(df_one, overlap_thr = 0.50) {
  if (nrow(df_one) <= 1) return(df_one)
  
  df_one <- df_one %>%
    arrange(i_evalue, desc(domain_score), ali_from, ali_to)
  
  keep_idx <- integer(0)
  
  for (i in seq_len(nrow(df_one))) {
    if (length(keep_idx) == 0) {
      keep_idx <- i
      next
    }
    
    ov_fracs <- vapply(
      keep_idx,
      function(j) interval_overlap_frac(
        df_one$ali_from[i], df_one$ali_to[i],
        df_one$ali_from[j], df_one$ali_to[j]
      ),
      numeric(1)
    )
    
    if (all(ov_fracs < overlap_thr, na.rm = TRUE)) {
      keep_idx <- c(keep_idx, i)
    }
  }
  
  df_one[sort(keep_idx), , drop = FALSE]
}

make_ordered_factor <- function(x, by) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unknown"
  by[is.na(by) | is.nan(by)] <- 0
  lev <- x[order(by, x)]
  factor(x, levels = unique(lev))
}

#Read ORFanage merged table
orf <- fread(merged_orf_file, data.table = FALSE)

if (!"isoform" %in% colnames(orf)) stop("Column 'isoform' not found in merged ORFanage table.")
orf$isoform <- as.character(orf$isoform)

if (!"associated_gene_final" %in% colnames(orf)) {
  if ("associated_gene" %in% colnames(orf)) {
    orf$associated_gene_final <- ifelse(
      !is.na(orf$associated_gene) & orf$associated_gene != "",
      as.character(orf$associated_gene),
      NA_character_
    )
  } else {
    orf$associated_gene_final <- NA_character_
  }
}

#Read Pfam domtblout
dom_lines <- readLines(pfam_domtblout_file, warn = FALSE)
dom_lines <- dom_lines[!grepl("^#", dom_lines)]

cat("Non-comment lines in domtblout: ", length(dom_lines), "\n")
if (length(dom_lines) == 0) {
  stop("No non-comment lines found in domtblout.")
}

parse_domtbl_line <- function(x) {
  parts <- strsplit(trimws(x), "\\s+")[[1]]
  if (length(parts) < 23) return(NULL)
  
  desc <- if (length(parts) > 22) paste(parts[23:length(parts)], collapse = " ") else ""
  
  data.frame(
    target_name      = parts[1],
    target_accession = parts[2],
    tlen             = parts[3],
    query_name       = parts[4],
    query_accession  = parts[5],
    qlen             = parts[6],
    full_evalue      = parts[7],
    full_score       = parts[8],
    full_bias        = parts[9],
    domain_num       = parts[10],
    domain_of        = parts[11],
    c_evalue         = parts[12],
    i_evalue         = parts[13],
    domain_score     = parts[14],
    domain_bias      = parts[15],
    hmm_from         = parts[16],
    hmm_to           = parts[17],
    ali_from         = parts[18],
    ali_to           = parts[19],
    env_from         = parts[20],
    env_to           = parts[21],
    acc              = parts[22],
    description      = desc,
    stringsAsFactors = FALSE
  )
}

pfam_list <- lapply(dom_lines, parse_domtbl_line)
pfam_list <- pfam_list[!vapply(pfam_list, is.null, logical(1))]

pfam_raw <- bind_rows(pfam_list)

cat("Rows parsed into pfam_raw: ", nrow(pfam_raw), "\n")
cat("Distinct isoforms in pfam_raw: ", dplyr::n_distinct(pfam_raw$query_name), "\n")

num_cols <- c(
  "tlen","qlen","full_evalue","full_score","full_bias",
  "domain_num","domain_of",
  "c_evalue","i_evalue","domain_score","domain_bias",
  "hmm_from","hmm_to","ali_from","ali_to","env_from","env_to","acc"
)
for (cc in num_cols) pfam_raw[[cc]] <- safe_num(pfam_raw[[cc]])

pfam_raw <- pfam_raw %>%
  mutate(
    isoform = trimws(as.character(query_name)),
    pfam_name = as.character(target_name),
    pfam_id = sub("\\..*$", "", as.character(target_accession)),
    hmm_cov = ifelse(!is.na(tlen) & tlen > 0, (hmm_to - hmm_from + 1) / tlen, NA_real_),
    ali_len = ali_to - ali_from + 1
  )

write.csv(pfam_raw, file.path(out_dir_pfam, "pfam_raw_hits.csv"), row.names = FALSE)

#only use i-Evalue threshold
pfam_filt <- pfam_raw %>%
  filter(
    !is.na(i_evalue),
    i_evalue <= 1e-3
  )

write.csv(pfam_filt, file.path(out_dir_pfam, "pfam_filtered_hits_relaxed.csv"), row.names = FALSE)

#Collapse overlapping redundant hits within each isoform
pfam_nr <- pfam_filt %>%
  group_by(isoform) %>%
  group_modify(~ collapse_overlapping_hits(.x, overlap_thr = 0.50)) %>%
  ungroup()

write.csv(pfam_nr, file.path(out_dir_pfam, "pfam_nonredundant_hits.csv"), row.names = FALSE)

#summarise Pfam at isoform level
pfam_iso <- pfam_nr %>%
  group_by(isoform) %>%
  summarise(
    n_pfam_domains_nr = n(),
    n_unique_pfam_nr = n_distinct(pfam_id),
    pfam_names_nr = paste(pfam_name[order(ali_from, ali_to)], collapse = ";"),
    pfam_ids_nr = paste(pfam_id[order(ali_from, ali_to)], collapse = ";"),
    pfam_architecture_nr = paste(pfam_name[order(ali_from, ali_to)], collapse = " | "),
    pfam_aa_covered_nr = collapse_intervals_width(ali_from, ali_to),
    best_i_evalue = min(i_evalue, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(pfam_iso, file.path(out_dir_pfam, "pfam_isoform_summary.csv"), row.names = FALSE)

#merge Pfam back into ORFanage table
df <- orf %>%
  left_join(pfam_iso, by = "isoform") %>%
  mutate(
    n_pfam_domains_nr = ifelse(is.na(n_pfam_domains_nr), 0L, as.integer(n_pfam_domains_nr)),
    n_unique_pfam_nr  = ifelse(is.na(n_unique_pfam_nr), 0L, as.integer(n_unique_pfam_nr)),
    pfam_aa_covered_nr = ifelse(is.na(pfam_aa_covered_nr), 0L, as.integer(pfam_aa_covered_nr)),
    has_pfam_nr = n_pfam_domains_nr > 0,
    protein_len_aa = ifelse(!is.na(query_len) & query_len > 0, floor(query_len / 3), NA_real_),
    pfam_frac_covered_nr = ifelse(
      !is.na(protein_len_aa) & protein_len_aa > 0,
      pfam_aa_covered_nr / protein_len_aa,
      NA_real_
    ),
    pfam_names_nr = ifelse(is.na(pfam_names_nr), "", pfam_names_nr),
    pfam_ids_nr = ifelse(is.na(pfam_ids_nr), "", pfam_ids_nr),
    pfam_architecture_nr = ifelse(is.na(pfam_architecture_nr), "", pfam_architecture_nr)
  )

write.csv(df, file.path(out_dir_pfam, "SQANTI_ORFanage_Pfam_merged_full.csv"), row.names = FALSE)

#build same-gene reference-like domain baseline
message("Building gene-wise reference-like domain baseline ...")

gene_baseline <- df %>%
  filter(
    !is.na(associated_gene_final), associated_gene_final != "",
    protein_consequence_simple == "Reference-like CDS"
  ) %>%
  group_by(associated_gene_final) %>%
  summarise(
    n_ref_like = n(),
    baseline_median_n_domains = median(n_pfam_domains_nr, na.rm = TRUE),
    baseline_median_frac_cov = median(pfam_frac_covered_nr, na.rm = TRUE),
    baseline_mode_architecture = mode_string(pfam_architecture_nr),
    .groups = "drop"
  )

write.csv(gene_baseline, file.path(out_dir_pfam, "pfam_gene_reference_like_baseline.csv"), row.names = FALSE)

df <- df %>%
  left_join(gene_baseline, by = "associated_gene_final") %>%
  mutate(
    pfam_architecture_impact = case_when(
      !has_orf ~ "No predicted ORF",
      has_orf & !has_pfam_nr ~ "Predicted protein, no Pfam domain",
      protein_consequence_simple == "Reference-like CDS" ~ "Reference-like architecture",
      !is.na(baseline_mode_architecture) &
        baseline_mode_architecture != "" &
        pfam_architecture_nr == baseline_mode_architecture ~ "Architecture similar to reference-like",
      !is.na(baseline_median_n_domains) &
        n_pfam_domains_nr < baseline_median_n_domains ~ "Reduced domain architecture",
      !is.na(baseline_median_n_domains) &
        n_pfam_domains_nr > baseline_median_n_domains ~ "Expanded / altered architecture",
      has_pfam_nr ~ "Altered domain architecture",
      TRUE ~ "Other"
    )
  )

write.csv(df, file.path(out_dir_pfam, "SQANTI_ORFanage_Pfam_with_architecture_impact.csv"), row.names = FALSE)

#RetNet subset
retnet_genes <- readLines(retnet_gene_file, warn = FALSE)
retnet_genes <- unique(retnet_genes[nzchar(retnet_genes)])

df_retnet <- df %>%
  filter(!is.na(associated_gene_final), associated_gene_final %in% retnet_genes)

write.csv(df_retnet, file.path(out_dir_pfam, "RetNet_Pfam_merged_full.csv"), row.names = FALSE)

#Summary tables

message("Writing summary tables ...")

#proteins with predicted ORF only
df_orf <- df %>% filter(has_orf)

pfam_overall <- df_orf %>%
  summarise(
    n_isoforms_with_orf = n(),
    n_with_pfam = sum(has_pfam_nr, na.rm = TRUE),
    frac_with_pfam = mean(has_pfam_nr, na.rm = TRUE),
    median_n_domains = median(n_pfam_domains_nr, na.rm = TRUE),
    median_frac_cov = median(pfam_frac_covered_nr, na.rm = TRUE)
  )
write.csv(pfam_overall, file.path(out_dir_pfam, "pfam_overall_summary.csv"), row.names = FALSE)

pfam_fl2 <- df_orf %>%
  filter(!is.na(FL) & FL >= 2) %>%
  summarise(
    n_isoforms_with_orf = n(),
    n_with_pfam = sum(has_pfam_nr, na.rm = TRUE),
    frac_with_pfam = mean(has_pfam_nr, na.rm = TRUE),
    median_n_domains = median(n_pfam_domains_nr, na.rm = TRUE),
    median_frac_cov = median(pfam_frac_covered_nr, na.rm = TRUE)
  )
write.csv(pfam_fl2, file.path(out_dir_pfam, "pfam_FL2_summary.csv"), row.names = FALSE)

pfam_by_consequence <- df %>%
  group_by(protein_consequence_simple) %>%
  summarise(
    n_isoforms = n(),
    n_with_pfam = sum(has_pfam_nr, na.rm = TRUE),
    frac_with_pfam = mean(has_pfam_nr, na.rm = TRUE),
    median_n_domains = median(n_pfam_domains_nr, na.rm = TRUE),
    median_frac_cov = median(pfam_frac_covered_nr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_isoforms))
write.csv(pfam_by_consequence, file.path(out_dir_pfam, "pfam_by_protein_consequence.csv"), row.names = FALSE)

pfam_architecture_summary <- df %>%
  count(pfam_architecture_impact, sort = TRUE) %>%
  mutate(frac = n / sum(n))
write.csv(pfam_architecture_summary, file.path(out_dir_pfam, "pfam_architecture_impact_summary.csv"), row.names = FALSE)

#RetNet summaries
df_retnet_orf <- df_retnet %>% filter(has_orf)

pfam_retnet_overall <- df_retnet_orf %>%
  summarise(
    n_isoforms_with_orf = n(),
    n_with_pfam = sum(has_pfam_nr, na.rm = TRUE),
    frac_with_pfam = mean(has_pfam_nr, na.rm = TRUE),
    median_n_domains = median(n_pfam_domains_nr, na.rm = TRUE),
    median_frac_cov = median(pfam_frac_covered_nr, na.rm = TRUE)
  )
write.csv(pfam_retnet_overall, file.path(out_dir_pfam, "RetNet_pfam_overall_summary.csv"), row.names = FALSE)

pfam_retnet_fl2 <- df_retnet_orf %>%
  filter(!is.na(FL) & FL >= 2) %>%
  summarise(
    n_isoforms_with_orf = n(),
    n_with_pfam = sum(has_pfam_nr, na.rm = TRUE),
    frac_with_pfam = mean(has_pfam_nr, na.rm = TRUE),
    median_n_domains = median(n_pfam_domains_nr, na.rm = TRUE),
    median_frac_cov = median(pfam_frac_covered_nr, na.rm = TRUE)
  )
write.csv(pfam_retnet_fl2, file.path(out_dir_pfam, "RetNet_pfam_FL2_summary.csv"), row.names = FALSE)

pfam_retnet_by_consequence <- df_retnet %>%
  group_by(protein_consequence_simple) %>%
  summarise(
    n_isoforms = n(),
    n_with_pfam = sum(has_pfam_nr, na.rm = TRUE),
    frac_with_pfam = mean(has_pfam_nr, na.rm = TRUE),
    median_n_domains = median(n_pfam_domains_nr, na.rm = TRUE),
    median_frac_cov = median(pfam_frac_covered_nr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_isoforms))
write.csv(pfam_retnet_by_consequence, file.path(out_dir_pfam, "RetNet_pfam_by_protein_consequence.csv"), row.names = FALSE)

pfam_retnet_architecture_summary <- df_retnet %>%
  count(pfam_architecture_impact, sort = TRUE) %>%
  mutate(frac = n / sum(n))
write.csv(pfam_retnet_architecture_summary, file.path(out_dir_pfam, "RetNet_pfam_architecture_impact_summary.csv"), row.names = FALSE)

#Top Pfam domains among altered RetNet isoforms
top_retnet_domains <- df_retnet %>%
  filter(likely_protein_altering, pfam_ids_nr != "" & !is.na(pfam_ids_nr)) %>%
  select(isoform, associated_gene_final, pfam_ids_nr) %>%
  separate_rows(pfam_ids_nr, sep = ";") %>%
  filter(pfam_ids_nr != "") %>%
  count(pfam_ids_nr, sort = TRUE, name = "n_isoforms")
write.csv(top_retnet_domains, file.path(out_dir_pfam, "RetNet_top_pfam_domains.csv"), row.names = FALSE)

#Example genes
example_genes <- c("ABCA4", "CNGB1", "PDE6A", "GNAT1", "SAG")
gene_examples <- df_retnet %>%
  filter(associated_gene_final %in% example_genes) %>%
  select(
    isoform, associated_gene_final, structural_category, FL,
    protein_consequence_simple, likely_protein_altering,
    n_pfam_domains_nr, pfam_frac_covered_nr, pfam_names_nr,
    pfam_architecture_nr, pfam_architecture_impact
  ) %>%
  arrange(associated_gene_final, desc(likely_protein_altering), desc(FL), isoform)
write.csv(gene_examples, file.path(out_dir_pfam, "RetNet_example_gene_domain_table.csv"), row.names = FALSE)

#per-gene impact summary
retnet_gene_impact_summary <- df_retnet %>%
  group_by(associated_gene_final) %>%
  summarise(
    n_isoforms = n(),
    n_with_pfam = sum(has_pfam_nr, na.rm = TRUE),
    frac_with_pfam = mean(has_pfam_nr, na.rm = TRUE),
    n_reduced_arch = sum(pfam_architecture_impact == "Reduced domain architecture", na.rm = TRUE),
    n_altered_arch = sum(pfam_architecture_impact %in% c("Altered domain architecture", "Expanded / altered architecture"), na.rm = TRUE),
    n_no_pfam = sum(pfam_architecture_impact == "Predicted protein, no Pfam domain", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_reduced_arch + n_altered_arch), desc(n_isoforms))
write.csv(retnet_gene_impact_summary, file.path(out_dir_pfam, "RetNet_gene_pfam_impact_summary.csv"), row.names = FALSE)

#summarize the orfanage and pfam mapping results
#Paths
pfam_merged_file <- "pfam_analysis_parallel/pfam_integrated_results_v2/SQANTI_ORFanage_Pfam_with_architecture_impact.csv"
pfam_hits_file   <- "pfam_analysis_parallel/pfam_integrated_results_v2/pfam_nonredundant_hits.csv"
retnet_gene_file <- "orfanage_analysis/RetNet_gene_set_detected_in_dataset.txt"

out_dir <- "translation_impact_v3"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#plot style
base_size_main   <- 15
base_size_axis   <- 13
base_size_text   <- 4.3
base_size_legend <- 12

theme_paper <- function() {
  theme_bw(base_size = base_size_main) +
    theme(
      plot.title      = element_text(face = "bold", size = 17),
      axis.title      = element_text(size = 14),
      axis.text       = element_text(size = 12, colour = "black"),
      axis.text.y     = element_text(size = 12, colour = "black"),
      legend.title    = element_text(size = 12),
      legend.text     = element_text(size = 11),
      strip.text      = element_text(size = 12, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

safe_num <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "nan", "-", "Inf", "-Inf")] <- NA
  suppressWarnings(as.numeric(x))
}

save_plot_pdf_png <- function(p, file_base, width = 10, height = 7, dpi = 300) {
  ggsave(paste0(file_base, ".pdf"), p, width = width, height = height, useDingbats = FALSE)
  ggsave(paste0(file_base, ".png"), p, width = width, height = height, dpi = dpi)
}

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Cannot find file: ", path)
}

clean_bool <- function(x) {
  if (is.logical(x)) return(x)
  x <- as.character(x)
  x <- trimws(tolower(x))
  x %in% c("true", "t", "1", "yes", "y")
}

clean_nonempty <- function(x) {
  !is.na(x) & trimws(as.character(x)) != ""
}

#Read input tables
df <- fread(pfam_merged_file, data.table = FALSE)
pfam_hits <- fread(pfam_hits_file, data.table = FALSE)
pbm <- fread(pbid_matrix_file, data.table = FALSE)

df$isoform <- as.character(df$isoform)
pfam_hits$isoform <- as.character(pfam_hits$isoform)

if (!"PBid" %in% colnames(pbm)) {
  stop("Column 'PBid' not found in PBid_matrix_allCT.csv")
}
pbm$PBid <- as.character(pbm$PBid)

#Basic checks
if (!"associated_gene_final" %in% colnames(df)) {
  if ("associated_gene" %in% colnames(df)) {
    df$associated_gene_final <- df$associated_gene
  } else {
    df$associated_gene_final <- NA_character_
  }
}

if (!"protein_consequence_simple" %in% colnames(df)) {
  stop("Column 'protein_consequence_simple' not found in ", pfam_merged_file)
}

if (!"pfam_architecture_impact" %in% colnames(df)) {
  stop("Column 'pfam_architecture_impact' not found in ", pfam_merged_file)
}

num_cols_df <- intersect(
  c("FL", "query_len", "template_len", "n_pfam_domains_nr",
    "n_unique_pfam_nr", "pfam_aa_covered_nr", "pfam_frac_covered_nr"),
  colnames(df)
)
for (cc in num_cols_df) df[[cc]] <- safe_num(df[[cc]])

num_cols_hits <- intersect(
  c("ali_from", "ali_to", "i_evalue", "domain_score"),
  colnames(pfam_hits)
)
for (cc in num_cols_hits) pfam_hits[[cc]] <- safe_num(pfam_hits[[cc]])

#PBid matrix metadata
colnames(pbm) <- make.names(colnames(pbm), unique = TRUE)
meta_cols_pbm <- c(
  "PBid", "chrom", "strand", "length", "structural_category", "subcategory",
  "associated_gene", "within_cage_peak", "polyA_motif"
)
meta_cols_pbm <- intersect(meta_cols_pbm, colnames(pbm))

celltype_cols <- setdiff(colnames(pbm), meta_cols_pbm)

numeric_like <- vapply(celltype_cols, function(cc) {
  vals <- suppressWarnings(as.numeric(pbm[[cc]]))
  mean(is.na(vals)) < 0.95
}, logical(1))

celltype_cols <- celltype_cols[numeric_like]

for (cc in celltype_cols) {
  pbm[[cc]] <- suppressWarnings(as.numeric(pbm[[cc]]))
  pbm[[cc]][is.na(pbm[[cc]])] <- 0
}

expr_mat <- as.matrix(pbm[, celltype_cols, drop = FALSE])
storage.mode(expr_mat) <- "numeric"

row_total <- rowSums(expr_mat)
row_max <- apply(expr_mat, 1, max)
dominant_idx <- if (ncol(expr_mat) > 0) max.col(expr_mat, ties.method = "first") else rep(NA_integer_, nrow(pbm))
dominant_celltype <- if (length(celltype_cols) > 0) celltype_cols[dominant_idx] else rep(NA_character_, nrow(pbm))
dominant_celltype[row_total == 0] <- NA_character_

pbm_meta <- pbm %>%
  transmute(
    isoform = PBid,
    associated_gene_pbm = if ("associated_gene" %in% colnames(pbm)) as.character(associated_gene) else NA_character_,
    within_cage_peak = if ("within_cage_peak" %in% colnames(pbm)) clean_bool(within_cage_peak) else NA,
    polyA_motif = if ("polyA_motif" %in% colnames(pbm)) as.character(polyA_motif) else NA_character_,
    dominant_celltype = dominant_celltype
  )

#Merge metadata
df2 <- df %>%
  left_join(pbm_meta, by = "isoform") %>%
  mutate(
    associated_gene_final = case_when(
      clean_nonempty(associated_gene_final) ~ associated_gene_final,
      clean_nonempty(associated_gene_pbm) ~ associated_gene_pbm,
      TRUE ~ NA_character_
    )
  )

#High-confidence subset
#FL >= 2 + within_cage_peak + polyA motif present
if (!"FL" %in% colnames(df2)) stop("Column 'FL' not found in main table.")

df_hc <- df2 %>%
  mutate(
    is_high_conf = !is.na(FL) & FL >= 2 &
      !is.na(within_cage_peak) & within_cage_peak &
      clean_nonempty(polyA_motif)
  ) %>%
  filter(is_high_conf)


#RetNet subset
retnet_genes <- readLines(retnet_gene_file, warn = FALSE)
retnet_genes <- unique(retnet_genes[clean_nonempty(retnet_genes)])

df_hc_retnet <- df_hc %>%
  filter(!is.na(associated_gene_final), associated_gene_final %in% retnet_genes)

#Consequence ordering
consequence_levels <- c(
  "Reference-like CDS",
  "In-frame altered CDS",
  "Predicted truncation",
  "Frameshift / out-of-frame alteration",
  "Predicted extension",
  "Novel ORF without matched template",
  "No predicted ORF",
  "Other"
)

df_hc$protein_consequence_simple <- factor(
  ifelse(clean_nonempty(df_hc$protein_consequence_simple), df_hc$protein_consequence_simple, "Other"),
  levels = consequence_levels
)

df_hc_retnet$protein_consequence_simple <- factor(
  ifelse(clean_nonempty(df_hc_retnet$protein_consequence_simple), df_hc_retnet$protein_consequence_simple, "Other"),
  levels = consequence_levels
)

#figure 3B
#all high-confidence isoforms
d1_df <- df_hc %>%
  count(protein_consequence_simple, name = "n") %>%
  filter(!is.na(protein_consequence_simple)) %>%
  mutate(frac = n / sum(n))

p_d1 <- ggplot(d1_df, aes(x = frac, y = fct_rev(protein_consequence_simple))) +
  geom_col(fill = "grey35") +
  geom_text(
    aes(label = paste0(comma(n), " (", percent(frac, accuracy = 0.1), ")")),
    hjust = -0.05, size = base_size_text
  ) +
  scale_x_continuous(
    labels = percent,
    limits = c(0, max(d1_df$frac) * 1.24 + 1e-6)
  ) +
  labs(
    title = "Predicted translational consequences in all high-confidence isoforms",
    x = "Fraction of isoforms",
    y = NULL
  ) +
  theme_paper()

save_plot_pdf_png(
  p_d1,
  file.path(out_dir, "panel_D1_all_highconf_protein_consequence"),
  width = 11.5, height = 6.2
)

write.csv(
  d1_df,
  file.path(out_dir, "panel_D1_all_highconf_protein_consequence_source_data.csv"),
  row.names = FALSE
)

#figure 3C
#high-confidence RetNet isoforms
d2_df <- df_hc_retnet %>%
  count(protein_consequence_simple, name = "n") %>%
  filter(!is.na(protein_consequence_simple)) %>%
  mutate(frac = n / sum(n))

p_d2 <- ggplot(d2_df, aes(x = frac, y = fct_rev(protein_consequence_simple))) +
  geom_col(fill = "grey35") +
  geom_text(
    aes(label = paste0(comma(n), " (", percent(frac, accuracy = 0.1), ")")),
    hjust = -0.05, size = base_size_text
  ) +
  scale_x_continuous(
    labels = percent,
    limits = c(0, max(d2_df$frac) * 1.24 + 1e-6)
  ) +
  labs(
    title = "Predicted translational consequences in high-confidence RetNet isoforms",
    x = "Fraction of isoforms",
    y = NULL
  ) +
  theme_paper()

save_plot_pdf_png(
  p_d2,
  file.path(out_dir, "panel_D2_highconf_RetNet_protein_consequence"),
  width = 11.5, height = 6.2
)

write.csv(
  d2_df,
  file.path(out_dir, "panel_D2_highconf_RetNet_protein_consequence_source_data.csv"),
  row.names = FALSE
)

#figure 3D
#Pfam impact by protein consequence
df_hc_e1 <- df_hc %>%
  mutate(
    pfam_group = case_when(
      pfam_architecture_impact %in% c("Reference-like architecture", "Architecture similar to reference-like") ~
        "Reference-like / similar",
      pfam_architecture_impact %in% c("Reduced domain architecture") ~
        "Reduced",
      pfam_architecture_impact %in% c("Altered domain architecture", "Expanded / altered architecture") ~
        "Altered / expanded",
      pfam_architecture_impact %in% c("Predicted protein, no Pfam domain") ~
        "No Pfam domain",
      pfam_architecture_impact %in% c("No predicted ORF") ~
        "No predicted ORF",
      TRUE ~ "Other"
    ),
    pfam_group = factor(
      pfam_group,
      levels = c(
        "Reference-like / similar",
        "Reduced",
        "Altered / expanded",
        "No Pfam domain",
        "No predicted ORF",
        "Other"
      )
    ),
    protein_consequence_simple = factor(protein_consequence_simple, levels = consequence_levels)
  )

e1_df <- df_hc_e1 %>%
  count(protein_consequence_simple, pfam_group, name = "n") %>%
  filter(!is.na(protein_consequence_simple), !is.na(pfam_group)) %>%
  group_by(protein_consequence_simple) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

p_e1 <- ggplot(e1_df, aes(x = protein_consequence_simple, y = frac, fill = pfam_group)) +
  geom_col(position = "fill", width = 0.86) +
  geom_text(
    data = e1_df %>% filter(frac >= 0.08),
    aes(label = percent(frac, accuracy = 1)),
    position = position_fill(vjust = 0.5),
    size = 3.8
  ) +
  scale_y_continuous(labels = percent) +
  labs(
    title = "Pfam architecture impact across high-confidence isoforms",
    x = NULL,
    y = "Fraction within protein consequence class",
    fill = "Pfam impact"
  ) +
  theme_paper() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 12),
    legend.position = "right"
  )

save_plot_pdf_png(
  p_e1,
  file.path(out_dir, "panel_E1_highconf_pfam_impact_by_protein_consequence"),
  width = 12.5, height = 7.2
)

write.csv(
  e1_df,
  file.path(out_dir, "panel_E1_highconf_pfam_impact_by_protein_consequence_source_data.csv"),
  row.names = FALSE
)


#descriptive statistic
df_test <- bind_rows(
  df_hc %>% mutate(subset = "All high-confidence"),
  df_hc_retnet %>% mutate(subset = "RetNet high-confidence")
) %>%
  mutate(
    likely_protein_altering = protein_consequence_simple %in% c(
      "Novel ORF without matched template",
      "Frameshift / out-of-frame alteration",
      "Predicted truncation",
      "Predicted extension",
      "In-frame altered CDS"
    )
  )

test_tbl <- table(df_test$subset, df_test$likely_protein_altering)
fisher_res <- fisher.test(test_tbl)

capture.output(
  list(contingency_table = test_tbl, fisher_test = fisher_res),
  file = file.path(out_dir, "D1_D2_likely_protein_altering_fisher_test.txt")
)


#Save full high-confidence isoform table
write.csv(
  df_hc,
  file.path(out_dir, "df_hc_high_confidence_isoforms.csv"),
  row.names = FALSE
)

saveRDS(
  df_hc,
  file.path(out_dir, "df_hc_high_confidence_isoforms.rds")
)


#PSI bootstrap analysis following scisorseqr exon logic
#bootstrap 95% CIs
#pooled summaries

#extract the position of the selected exons
events <- tribble(
  ~event_id,      ~gene_id,              ~target_exon,
  "ASPH_exon",    "ENSG00000198363.18",  "chr8_61682439_61682483_-",
  "EXOC1_exon",   "ENSG00000090989.18",  "chr4_55866872_55866892_+"
)

#input files
allinfo_file <- "AllInfo_IncompleteReads_NoSampleNumber"
donor_run_file <- "donor_run.csv"

#helper functions
split_chain <- function(x) {
  if (is.na(x) || x == "" || x == "NA") return(character(0))
  y <- unlist(strsplit(x, ";%;", fixed = TRUE))
  y <- trimws(y)
  y[y != ""]
}

parse_exon_string <- function(exon_str) {
  x <- unlist(strsplit(exon_str, "_", fixed = TRUE))
  if (length(x) != 4) {
    return(tibble(
      chr = NA_character_,
      start = NA_real_,
      end = NA_real_,
      strand = NA_character_
    ))
  }
  tibble(
    chr = x[1],
    start = as.numeric(x[2]),
    end = as.numeric(x[3]),
    strand = x[4]
  )
}

#match scisorseqr logic:
#read start = start coord of first exon in ExonChain
#read end = end coord of last exon in ExonChain
get_read_span_from_exon_chain <- function(exon_chain) {
  exons <- split_chain(exon_chain)
  if (length(exons) == 0) {
    return(tibble(
      read_chr = NA_character_,
      read_start = NA_real_,
      read_end = NA_real_,
      read_strand = NA_character_
    ))
  }
  
  first_exon <- parse_exon_string(exons[1])
  last_exon  <- parse_exon_string(exons[length(exons)])
  
  tibble(
    read_chr = first_exon$chr,
    read_start = first_exon$start,
    read_end = last_exon$end,
    read_strand = first_exon$strand
  )
}

#whether target exon is present in the read
is_included_exon <- function(exon_chain, target_exon) {
  exons <- split_chain(exon_chain)
  target_exon %in% exons
}

#whether read spans the exon coordinate interval
#exactly following ExonCounting.R (scisorseqr package) idea: s >= start & e <= end
is_spanning_exon <- function(exon_chain, target_exon) {
  span <- get_read_span_from_exon_chain(exon_chain)
  target <- parse_exon_string(target_exon)
  
  if (any(is.na(span)) || any(is.na(target))) return(FALSE)
  
  same_chr <- identical(as.character(span$read_chr), as.character(target$chr))
  same_strand <- identical(as.character(span$read_strand), as.character(target$strand))
  
  if (!same_chr || !same_strand) return(FALSE)
  
  (target$start >= span$read_start) && (target$end <= span$read_end)
}

bootstrap_psi_scisorseqr <- function(spanning_vec, included_vec, n_boot = 1000) {
  #spanning_vec and included_vec should be aligned logical vectors over reads
  keep <- !is.na(spanning_vec) & !is.na(included_vec)
  spanning_vec <- spanning_vec[keep]
  included_vec <- included_vec[keep]
  
  #only spanning reads contribute to PSI
  if (length(spanning_vec) == 0 || sum(spanning_vec) == 0) {
    return(tibble(
      mean_psi = NA_real_,
      median_psi = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_
    ))
  }
  
  dat <- tibble(
    spanning = spanning_vec,
    included = included_vec
  ) %>%
    filter(spanning)
  
  if (nrow(dat) == 0) {
    return(tibble(
      mean_psi = NA_real_,
      median_psi = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_
    ))
  }
  
  #each spanning read is either included or excluded
  inc_status <- dat$included
  
  if (all(inc_status)) {
    return(tibble(
      mean_psi = 1,
      median_psi = 1,
      ci_lower = 1,
      ci_upper = 1
    ))
  }
  
  if (all(!inc_status)) {
    return(tibble(
      mean_psi = 0,
      median_psi = 0,
      ci_lower = 0,
      ci_upper = 0
    ))
  }
  
  boot_psi <- replicate(n_boot, {
    sampled <- sample(inc_status, size = length(inc_status), replace = TRUE)
    mean(sampled)
  })
  
  tibble(
    mean_psi = mean(boot_psi, na.rm = TRUE),
    median_psi = median(boot_psi, na.rm = TRUE),
    ci_lower = quantile(boot_psi, 0.025, na.rm = TRUE),
    ci_upper = quantile(boot_psi, 0.975, na.rm = TRUE)
  )
}

summarise_one_group_scisorseqr <- function(df_group, n_boot = 1000) {
  total <- sum(df_group$spanning, na.rm = TRUE)
  inclusion <- sum(df_group$spanning & df_group$included, na.rm = TRUE)
  exclusion <- total - inclusion
  psi <- ifelse(total > 0, inclusion / total, NA_real_)
  
  boot <- bootstrap_psi_scisorseqr(
    spanning_vec = df_group$spanning,
    included_vec = df_group$included,
    n_boot = n_boot
  )
  
  bind_cols(
    tibble(
      inclusion = inclusion,
      exclusion = exclusion,
      total = total,
      PSI = psi
    ),
    boot
  )
}

run_one_event_scisorseqr <- function(allinfo, event_id, gene_id, target_exon, n_boot = 1000) {
  cat("\n====================================================\n")
  cat("Running event:", event_id, "\n")
  cat("Gene:", gene_id, "\n")
  cat("Target exon:", target_exon, "\n")
  
  df_gene <- allinfo %>%
    filter(Gene == gene_id)
  
  cat("Reads in gene:", nrow(df_gene), "\n")
  
  if (nrow(df_gene) == 0) {
    stop(paste("No reads found for gene:", gene_id))
  }
  
  event_reads <- df_gene %>%
    mutate(
      event_id = event_id,
      target_exon = target_exon
    )
  
  event_reads$included <- vapply(
    event_reads$ExonChain,
    function(x) is_included_exon(x, target_exon = target_exon),
    logical(1)
  )
  
  event_reads$spanning <- vapply(
    event_reads$ExonChain,
    function(x) is_spanning_exon(x, target_exon = target_exon),
    logical(1)
  )
  
  event_reads$call <- dplyr::case_when(
    event_reads$spanning & event_reads$included ~ "included",
    event_reads$spanning & !event_reads$included ~ "excluded",
    TRUE ~ "not_spanning"
  )
  
  write.csv(
    event_reads,
    paste0(event_id, "_read_level_scisorseqr_style_classification.csv"),
    row.names = FALSE
  )
  
  #donor x celltype
  summary_donor <- event_reads %>%
    group_by(event_id, Donor, Celltype) %>%
    group_modify(~ summarise_one_group_scisorseqr(.x, n_boot = n_boot)) %>%
    ungroup()
  
  #pooled by celltype
  summary_pooled <- event_reads %>%
    group_by(event_id, Celltype) %>%
    group_modify(~ summarise_one_group_scisorseqr(.x, n_boot = n_boot)) %>%
    ungroup()
  
  #pooled overall
  summary_overall <- event_reads %>%
    group_by(event_id) %>%
    group_modify(~ summarise_one_group_scisorseqr(.x, n_boot = n_boot)) %>%
    ungroup()
  
  write.csv(
    summary_donor,
    paste0(event_id, "_psi_bootstrap_by_donor_celltype_scisorseqr_style.csv"),
    row.names = FALSE
  )
  write.csv(
    summary_pooled,
    paste0(event_id, "_psi_bootstrap_pooled_by_celltype_scisorseqr_style.csv"),
    row.names = FALSE
  )
  write.csv(
    summary_overall,
    paste0(event_id, "_psi_bootstrap_overall_scisorseqr_style.csv"),
    row.names = FALSE
  )
  
  #donor-level plot
  p_donor <- ggplot(summary_donor, aes(x = Celltype, y = mean_psi, color = Donor)) +
    geom_point(position = position_dodge(width = 0.45), size = 3) +
    geom_errorbar(
      aes(ymin = ci_lower, ymax = ci_upper),
      position = position_dodge(width = 0.45),
      width = 0.18,
      linewidth = 0.6
    ) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    labs(
      title = paste0(event_id, " donor-level PSI (scisorseqr style)"),
      subtitle = paste0("Target exon: ", target_exon),
      x = "Cell type",
      y = "PSI"
    ) +
    theme_bw(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  ggsave(
    paste0(event_id, "_PSI_bootstrap_donor_level_scisorseqr_style.tiff"),
    p_donor,
    width = 12, height = 6, units = "in",
    dpi = 600, compression = "lzw"
  )
  
  #pooled plot
  p_pooled <- ggplot(summary_pooled, aes(x = Celltype, y = mean_psi)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.18, linewidth = 0.6) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    labs(
      title = paste0(event_id, " pooled PSI (scisorseqr style)"),
      subtitle = paste0("Target exon: ", target_exon),
      x = "Cell type",
      y = "PSI"
    ) +
    theme_bw(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  ggsave(
    paste0(event_id, "_PSI_bootstrap_pooled_scisorseqr_style.tiff"),
    p_pooled,
    width = 10, height = 6, units = "in",
    dpi = 600, compression = "lzw"
  )
  
  list(
    event_reads = event_reads,
    summary_donor = summary_donor,
    summary_pooled = summary_pooled,
    summary_overall = summary_overall
  )
}

#read inputs
donor_run <- read.csv(donor_run_file, stringsAsFactors = FALSE)
allinfo <- read.delim(
  allinfo_file,
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE
)

colnames(allinfo) <- c(
  "ReadName", "Gene", "Celltype", "Barcode", "UMI",
  "IntronChain", "TSS", "PolyA", "ExonChain",
  "Status", "NumberofIntrons"
)

allinfo <- allinfo %>%
  mutate(
    Run = sub("/.*$", "", ReadName)
  ) %>%
  left_join(donor_run, by = "Run")

cat("Total reads:", nrow(allinfo), "\n")
cat("Reads with donor assigned:", sum(!is.na(allinfo$Donor)), "\n")

#quick gene ID check
cat("\nGene IDs detected for ASPH-like query:\n")
print(unique(allinfo$Gene[grepl("ENSG00000198363", allinfo$Gene)]))

cat("\nGene IDs detected for EXOC1-like query:\n")
print(unique(allinfo$Gene[grepl("ENSG00000090989", allinfo$Gene)]))

#run all events
results_list <- pmap(
  events,
  function(event_id, gene_id, target_exon) {
    run_one_event_scisorseqr(
      allinfo = allinfo,
      event_id = event_id,
      gene_id = gene_id,
      target_exon = target_exon,
      n_boot = 1000
    )
  }
)

names(results_list) <- events$event_id
saveRDS(results_list, "PSI_bootstrap_results_scisorseqr_style.rds")

#combined outputs
all_donor <- bind_rows(lapply(results_list, function(x) x$summary_donor))
all_pooled <- bind_rows(lapply(results_list, function(x) x$summary_pooled))
all_overall <- bind_rows(lapply(results_list, function(x) x$summary_overall))

write.csv(
  all_donor,
  "ALL_events_PSI_bootstrap_by_donor_celltype_scisorseqr_style.csv",
  row.names = FALSE
)
write.csv(
  all_pooled,
  "ALL_events_PSI_bootstrap_pooled_by_celltype_scisorseqr_style.csv",
  row.names = FALSE
)
write.csv(
  all_overall,
  "ALL_events_PSI_bootstrap_overall_scisorseqr_style.csv",
  row.names = FALSE
)

#combined donor plot
p_all_donor <- ggplot(all_donor, aes(x = Celltype, y = mean_psi, color = Donor)) +
  geom_point(position = position_dodge(width = 0.45), size = 2.5) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    position = position_dodge(width = 0.45),
    width = 0.18,
    linewidth = 0.5
  ) +
  facet_wrap(~ event_id, ncol = 1) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(
    title = "Donor-level PSI bootstrap for manuscript exons",
    subtitle = "scisorseqr-style exon PSI: inclusion / spanning reads",
    x = "Cell type",
    y = "PSI"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(
  "ALL_events_PSI_bootstrap_donor_level_scisorseqr_style.tiff",
  p_all_donor,
  width = 12, height = 10, units = "in",
  dpi = 600, compression = "lzw"
)

#combined pooled plot
p_all_pooled <- ggplot(all_pooled, aes(x = Celltype, y = mean_psi)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.18, linewidth = 0.5) +
  facet_wrap(~ event_id, ncol = 1) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(
    title = "Pooled PSI bootstrap for manuscript exons",
    subtitle = "scisorseqr-style exon PSI: inclusion / spanning reads",
    x = "Cell type",
    y = "PSI"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(
  "ALL_events_PSI_bootstrap_pooled_scisorseqr_style.tiff",
  p_all_pooled,
  width = 10, height = 10, units = "in",
  dpi = 600, compression = "lzw"
)


#only show selected cell types for each event
#because this section indicate the significant differential exon usage between cell type pairs
#so retain cell types within comparison
plot_celltype_keep <- list(
  ASPH_exon  = c("amacrine", "bipolar", "horizontal", "mullerglia", "rod"),
  EXOC1_exon = c("amacrine", "mullerglia", "rod")
)

#original display order
manuscript_celltype_order <- c(
  "amacrine", "bipolar", "horizontal", "mullerglia", "rod"
)

#helper function to standardise cell type names for matching
normalise_celltype <- function(x) {
  x %>%
    as.character() %>%
    tolower() %>%
    gsub("[_ -]", "", .)
}

#add relaxed matching column
all_donor_plot <- all_donor %>%
  mutate(
    Celltype_match = normalise_celltype(Celltype)
  )

all_pooled_plot <- all_pooled %>%
  mutate(
    Celltype_match = normalise_celltype(Celltype)
  )

#convert keep list into dataframe
plot_celltype_keep_df <- bind_rows(
  lapply(names(plot_celltype_keep), function(ev) {
    tibble(
      event_id = ev,
      Celltype_match = normalise_celltype(plot_celltype_keep[[ev]])
    )
  })
)

#filter donor-level summary for manuscript cell types only
all_donor_manuscript <- all_donor_plot %>%
  inner_join(
    plot_celltype_keep_df,
    by = c("event_id", "Celltype_match")
  ) %>%
  mutate(
    Celltype_order_match = factor(
      Celltype_match,
      levels = normalise_celltype(manuscript_celltype_order)
    )
  ) %>%
  arrange(event_id, Celltype_order_match) %>%
  mutate(
    Celltype = factor(Celltype, levels = unique(Celltype))
  )

#filter pooled summary for manuscript cell types only
all_pooled_manuscript <- all_pooled_plot %>%
  inner_join(
    plot_celltype_keep_df,
    by = c("event_id", "Celltype_match")
  ) %>%
  mutate(
    Celltype_order_match = factor(
      Celltype_match,
      levels = normalise_celltype(manuscript_celltype_order)
    )
  ) %>%
  arrange(event_id, Celltype_order_match) %>%
  mutate(
    Celltype = factor(Celltype, levels = unique(Celltype))
  )

#save filtered summary tables
write.csv(
  all_donor_manuscript,
  "ALL_events_PSI_bootstrap_by_donor_celltype_MANUSCRIPT_CELLTYPES.csv",
  row.names = FALSE
)

write.csv(
  all_pooled_manuscript,
  "ALL_events_PSI_bootstrap_pooled_by_celltype_MANUSCRIPT_CELLTYPES.csv",
  row.names = FALSE
)

#combined donor-level plot, horizontal layout
p_all_donor_manuscript <- ggplot(
  all_donor_manuscript,
  aes(x = Celltype, y = mean_psi, color = Donor)
) +
  geom_point(
    position = position_dodge(width = 0.45),
    size = 2.8
  ) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    position = position_dodge(width = 0.45),
    width = 0.18,
    linewidth = 0.5
  ) +
  facet_wrap(
    ~ event_id,
    nrow = 1,
    scales = "free_x"
  ) +
  scale_y_continuous(
    labels = percent_format(),
    limits = c(0, 1)
  ) +
  labs(
    title = "Donor-level PSI bootstrap for manuscript exons",
    subtitle = "Only cell types shown in the manuscript are displayed",
    x = "Cell type",
    y = "PSI"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  "ALL_events_PSI_bootstrap_donor_level_MANUSCRIPT_CELLTYPES_horizontal.tiff",
  p_all_donor_manuscript,
  width = 12,
  height = 5.5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  "ALL_events_PSI_bootstrap_donor_level_MANUSCRIPT_CELLTYPES_horizontal.pdf",
  p_all_donor_manuscript,
  width = 12,
  height = 5.5,
  units = "in"
)

#combined pooled plot, horizontal layout
p_all_pooled_manuscript <- ggplot(
  all_pooled_manuscript,
  aes(x = Celltype, y = mean_psi)
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.18,
    linewidth = 0.6
  ) +
  facet_wrap(
    ~ event_id,
    nrow = 1,
    scales = "free_x"
  ) +
  scale_y_continuous(
    labels = percent_format(),
    limits = c(0, 1)
  ) +
  labs(
    title = "Pooled PSI bootstrap for manuscript exons",
    subtitle = "Only cell types shown in the manuscript are displayed",
    x = "Cell type",
    y = "PSI"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold")
  )

ggsave(
  "ALL_events_PSI_bootstrap_pooled_MANUSCRIPT_CELLTYPES_horizontal.tiff",
  p_all_pooled_manuscript,
  width = 10,
  height = 5.5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  "ALL_events_PSI_bootstrap_pooled_MANUSCRIPT_CELLTYPES_horizontal.pdf",
  p_all_pooled_manuscript,
  width = 10,
  height = 5.5,
  units = "in"
)

#generate individual plots
#figure 4D-F
event_ids_to_plot <- c("ASPH_exon", "EXOC1_exon")

for (ev in event_ids_to_plot) {
  pooled_one <- all_pooled_manuscript %>%
    filter(event_id == ev) %>%
    arrange(Celltype_order_match) %>%
    mutate(
      Celltype = factor(Celltype, levels = unique(Celltype))
    )
  
  p_pooled_one <- ggplot(
    pooled_one,
    aes(x = Celltype, y = mean_psi)
  ) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.18,
      linewidth = 0.6
    ) +
    scale_y_continuous(
      labels = percent_format(),
      limits = c(0, 1)
    ) +
    labs(
      title = paste0(ev, " pooled PSI"),
      subtitle = "Manuscript-selected cell types only",
      x = "Cell type",
      y = "PSI"
    ) +
    theme_bw(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  ggsave(
    paste0(ev, "_PSI_bootstrap_pooled.tiff"),
    p_pooled_one,
    width = 6.5,
    height = 5,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
  
  ggsave(
    paste0(ev, "_PSI_bootstrap_pooled.pdf"),
    p_pooled_one,
    width = 6.5,
    height = 5,
    units = "in"
  )
}


#download ABCA4 variants from clinVar and run with splicAI to estimate the effects of variants to splicing events
#paths
base_dir <- "D:/bioinformatic/isoseq/processed/from_Ray/revision/splicAI"
events_file <- file.path(base_dir, "ABCA4_spliceai_ready", "abca4_spliceai_comparison_events.tsv")
matches_file <- file.path(base_dir, "ABCA4_spliceai_match", "abca4_spliceai_variant_event_matches.tsv")
isoform_summary_file <- file.path(base_dir, "ABCA4_event_summary_v44", "abca4_isoform_event_summary.tsv")
clinvar_included_file <- file.path(base_dir, "clinvar_spliceai_input", "clinvar_ABCA4.included.tsv")
matrix_file_FL2_CAGE_polyA <- file.path(base_dir, "ABCA4_export", "ABCA4_PBid_matrix_FL2_CAGE_polyAmotif.csv")

out_dir_base <- file.path(base_dir, "ABCA4_support_overall_results_v2")
dir.create(out_dir_base, recursive = TRUE, showWarnings = FALSE)

#define functions
safe_num <- function(x) suppressWarnings(as.numeric(x))

clean_variant_name <- function(x) {
  x %>%
    str_replace("^NM_000350\\.3\\(ABCA4\\):", "") %>%
    str_replace("^ABCA4:", "")
}

first_non_na <- function(x) {
  y <- x[!is.na(x) & x != ""]
  if (length(y) == 0) return(NA_character_)
  y[[1]]
}

first_non_na_num <- function(x) {
  y <- x[!is.na(x)]
  if (length(y) == 0) return(NA_real_)
  y[[1]]
}

event_key_fun <- function(df) {
  df %>%
    mutate(
      event_key = paste(
        transcript_id,
        event_type,
        event_subtype,
        pb_exon_index,
        observed_start,
        observed_end,
        sep = "||"
      )
    )
}

get_col_if_exists <- function(df, colname) {
  if (colname %in% colnames(df)) return(as.character(df[[colname]]))
  rep(NA_character_, nrow(df))
}

get_variant_id_final <- function(df) {
  coalesce(
    get_col_if_exists(df, "variant_ID"),
    get_col_if_exists(df, "ID"),
    get_col_if_exists(df, "VariationID")
  )
}

standardize_matrix <- function(df, dataset_label) {
  if (!"PBid" %in% colnames(df)) {
    colnames(df)[1] <- "PBid"
  }
  df
}

prepare_matrix_metadata <- function(pbid_mat, dataset_label) {
  pbid_mat <- standardize_matrix(pbid_mat, dataset_label)
  
  if (dataset_label == "ALL_PBID") {
    celltype_cols <- setdiff(colnames(pbid_mat), "PBid")
    
    for (cc in intersect(celltype_cols, colnames(pbid_mat))) {
      pbid_mat[[cc]] <- safe_num(pbid_mat[[cc]])
    }
    
    pbid_meta <- pbid_mat %>%
      transmute(
        transcript_id = PBid,
        summary_count = NA_real_
      ) %>%
      distinct(transcript_id, .keep_all = TRUE)
    
    pbid_celltype_long <- pbid_mat %>%
      select(PBid, all_of(celltype_cols)) %>%
      pivot_longer(
        cols = all_of(celltype_cols),
        names_to = "cell_type_matrix",
        values_to = "cell_count"
      ) %>%
      mutate(cell_count = safe_num(cell_count)) %>%
      filter(!is.na(cell_count), cell_count > 0) %>%
      rename(transcript_id = PBid)
  } else {
    pbid_meta_cols <- c(
      "PBid", "chrom", "strand", "length", "structural_category", "subcategory",
      "associated_gene", "within_cage_peak", "polyA_motif", "summary count"
    )
    
    candidate_expr_cols <- setdiff(colnames(pbid_mat), pbid_meta_cols)
    celltype_cols <- candidate_expr_cols[!grepl("^Retina\\d+$", candidate_expr_cols)]
    
    for (cc in intersect(celltype_cols, colnames(pbid_mat))) {
      pbid_mat[[cc]] <- safe_num(pbid_mat[[cc]])
    }
    
    pbid_meta <- pbid_mat %>%
      transmute(
        transcript_id = PBid,
        chrom = if ("chrom" %in% colnames(pbid_mat)) chrom else NA_character_,
        strand = if ("strand" %in% colnames(pbid_mat)) strand else NA_character_,
        pbid_length = if ("length" %in% colnames(pbid_mat)) length else NA_real_,
        structural_category = if ("structural_category" %in% colnames(pbid_mat)) structural_category else NA_character_,
        subcategory = if ("subcategory" %in% colnames(pbid_mat)) subcategory else NA_character_,
        associated_gene = if ("associated_gene" %in% colnames(pbid_mat)) associated_gene else NA_character_,
        within_cage_peak = if ("within_cage_peak" %in% colnames(pbid_mat)) within_cage_peak else NA,
        polyA_motif = if ("polyA_motif" %in% colnames(pbid_mat)) polyA_motif else NA_character_,
        summary_count = if ("summary count" %in% colnames(pbid_mat)) `summary count` else NA_real_
      ) %>%
      distinct(transcript_id, .keep_all = TRUE)
    
    pbid_celltype_long <- pbid_mat %>%
      select(PBid, all_of(celltype_cols)) %>%
      pivot_longer(
        cols = all_of(celltype_cols),
        names_to = "cell_type_matrix",
        values_to = "cell_count"
      ) %>%
      mutate(cell_count = safe_num(cell_count)) %>%
      filter(!is.na(cell_count), cell_count > 0) %>%
      rename(transcript_id = PBid)
  }
  
  pbid_celltype_summary <- pbid_celltype_long %>%
    group_by(transcript_id) %>%
    summarise(
      expressed_cell_types = paste(sort(unique(cell_type_matrix)), collapse = ";"),
      n_expressed_cell_types = n_distinct(cell_type_matrix),
      dominant_cell_type = {
        tmp <- cur_data()
        max_n <- max(tmp$cell_count, na.rm = TRUE)
        paste(sort(unique(tmp$cell_type_matrix[tmp$cell_count == max_n])), collapse = ";")
      },
      dominant_cell_count = max(cell_count, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(
    pbid_mat = pbid_mat,
    pbids_in_dataset = unique(pbid_mat$PBid),
    celltype_cols = celltype_cols,
    pbid_meta = pbid_meta,
    pbid_celltype_long = pbid_celltype_long,
    pbid_celltype_summary = pbid_celltype_summary
  )
}

run_support_analysis <- function(
    dataset_label,
    matrix_file,
    out_dir,
    matches_raw,
    events_raw,
    isoform_meta,
    n_total_variants_included_common
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  pbid_mat <- read_csv(matrix_file, show_col_types = FALSE)
  prep <- prepare_matrix_metadata(pbid_mat, dataset_label)
  
  pbids_in_dataset <- prep$pbids_in_dataset
  pbid_meta <- prep$pbid_meta
  pbid_celltype_long <- prep$pbid_celltype_long
  pbid_celltype_summary <- prep$pbid_celltype_summary
  celltype_cols <- prep$celltype_cols
  
  matches_ds <- matches_raw %>%
    filter(transcript_id %in% pbids_in_dataset) %>%
    mutate(
      variant_short = if ("variant_Name" %in% colnames(.)) clean_variant_name(variant_Name) else NA_character_,
      variant_id_final = get_variant_id_final(.)
    ) %>%
    left_join(pbid_meta, by = "transcript_id") %>%
    left_join(isoform_meta, by = "transcript_id") %>%
    left_join(pbid_celltype_summary, by = "transcript_id") %>%
    mutate(
      cell_type_final = coalesce(dominant_cell_type, isoform_cell_type_label, cell_type)
    ) %>%
    distinct()
  
  events_ds <- events_raw %>%
    filter(transcript_id %in% pbids_in_dataset) %>%
    left_join(pbid_meta, by = "transcript_id") %>%
    left_join(isoform_meta, by = "transcript_id") %>%
    left_join(pbid_celltype_summary, by = "transcript_id") %>%
    mutate(
      cell_type_final = coalesce(dominant_cell_type, isoform_cell_type_label, cell_type)
    ) %>%
    distinct()
  
  matches_support_ds <- matches_ds %>%
    mutate(
      is_terminal = event_type %in% c("alternative_first_exon", "alternative_last_exon"),
      
      exact_AG_obsA = !is.na(diff_AG_obs_acceptor) & diff_AG_obs_acceptor == 0,
      exact_DG_obsD = !is.na(diff_DG_obs_donor) & diff_DG_obs_donor == 0,
      exact_AL_refA = !is.na(diff_AL_ref_acceptor) & diff_AL_ref_acceptor == 0,
      exact_DL_refD = !is.na(diff_DL_ref_donor) & diff_DL_ref_donor == 0,
      
      general_acceptor_component =
        (exact_AL_refA & !is.na(DS_AL) & DS_AL > 0) |
        (exact_AG_obsA & !is.na(DS_AG) & DS_AG > 0),
      
      general_donor_component =
        (exact_DL_refD & !is.na(DS_DL) & DS_DL > 0) |
        (exact_DG_obsD & !is.na(DS_DG) & DS_DG > 0),
      
      general_support = case_when(
        event_type == "alt_acceptor" ~ general_acceptor_component,
        event_type == "alt_donor" ~ general_donor_component,
        event_type == "novel_exon" ~ general_acceptor_component | general_donor_component,
        event_type %in% c("partial_intron_retention", "intron_retention", "exon_skipping") ~
          general_acceptor_component | general_donor_component,
        TRUE ~ FALSE
      ),
      
      support_acceptor_component =
        (exact_AL_refA & !is.na(DS_AL) & DS_AL >= 0.5) |
        (exact_AG_obsA & !is.na(DS_AG) & DS_AG >= 0.5),
      
      support_donor_component =
        (exact_DL_refD & !is.na(DS_DL) & DS_DL >= 0.5) |
        (exact_DG_obsD & !is.na(DS_DG) & DS_DG >= 0.5),
      
      support_05 = case_when(
        event_type == "alt_acceptor" ~ support_acceptor_component,
        event_type == "alt_donor" ~ support_donor_component,
        event_type == "novel_exon" ~ support_acceptor_component | support_donor_component,
        event_type %in% c("partial_intron_retention", "intron_retention", "exon_skipping") ~
          support_acceptor_component | support_donor_component,
        TRUE ~ FALSE
      ),
      
      support_mode = case_when(
        !support_05 ~ "none",
        support_acceptor_component & support_donor_component ~ "full_splicing_support",
        support_acceptor_component ~ "acceptor_support",
        support_donor_component ~ "donor_support",
        support_05 ~ "other_support",
        TRUE ~ "none"
      ),
      
      support_label = support_mode
    )
  
  events_main_ds <- events_ds %>%
    filter(!event_type %in% c("alternative_first_exon", "alternative_last_exon")) %>%
    distinct(event_key, .keep_all = TRUE)
  
  main_matches_ds <- matches_support_ds %>%
    filter(!is_terminal) %>%
    distinct()
  
  general_support_tbl_ds <- main_matches_ds %>%
    filter(general_support) %>%
    distinct()
  
  support_tbl_ds <- main_matches_ds %>%
    filter(support_05) %>%
    distinct()
  
  best_general_per_event_ds <- main_matches_ds %>%
    filter(general_support) %>%
    group_by(event_key) %>%
    arrange(desc(match_score), desc(maxDS), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
  
  best_support_per_event_ds <- main_matches_ds %>%
    filter(support_05) %>%
    mutate(
      support_rank = case_when(
        support_mode == "full_splicing_support" ~ 3,
        support_mode %in% c("acceptor_support", "donor_support") ~ 2,
        support_mode == "other_support" ~ 1,
        TRUE ~ 0
      )
    ) %>%
    group_by(event_key) %>%
    arrange(desc(support_rank), desc(match_score), desc(maxDS), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
  
  event_support_summary_ds <- events_main_ds %>%
    select(
      event_key, transcript_id, event_type, event_subtype,
      observed_start, observed_end,
      cell_type_final, dominant_cell_type, expressed_cell_types,
      n_expressed_cell_types, summary_count, major_event, priority_for_followup
    ) %>%
    distinct() %>%
    left_join(
      best_general_per_event_ds %>%
        select(
          event_key,
          general_support,
          general_variant_id_final = variant_id_final,
          general_variant_ID = variant_ID,
          general_VariationID = VariationID,
          general_variant_Name = variant_Name,
          general_variant_short = variant_short,
          general_match_score = match_score,
          general_maxDS = maxDS
        ),
      by = "event_key"
    ) %>%
    left_join(
      best_support_per_event_ds %>%
        select(
          event_key,
          support_05,
          support_mode,
          support_label,
          variant_id_final,
          variant_ID,
          VariationID,
          variant_Name,
          variant_short,
          match_score,
          maxDS
        ),
      by = "event_key"
    ) %>%
    mutate(
      general_support = replace_na(general_support, FALSE),
      support_mode_best = replace_na(support_mode, "none")
    )
  
  variant_overall_summary_ds <- main_matches_ds %>%
    filter(support_05) %>%
    mutate(
      variant_key = coalesce(as.character(variant_id_final), as.character(variant_ID), as.character(VariationID))
    ) %>%
    group_by(variant_key, variant_ID, VariationID, variant_Name, variant_short, VariantType) %>%
    summarise(
      n_supported_matches = n(),
      n_unique_events = n_distinct(event_key),
      n_unique_pbids = n_distinct(transcript_id),
      support_modes = paste(sort(unique(support_mode)), collapse = ","),
      matched_event_types = paste(sort(unique(event_type)), collapse = ","),
      .groups = "drop"
    ) %>%
    arrange(desc(n_unique_pbids), desc(n_unique_events))
  
  n_total_events <- n_distinct(events_main_ds$event_key)
  n_general_events <- sum(event_support_summary_ds$general_support, na.rm = TRUE)
  n_acceptor_events <- sum(event_support_summary_ds$support_mode_best == "acceptor_support", na.rm = TRUE)
  n_donor_events <- sum(event_support_summary_ds$support_mode_best == "donor_support", na.rm = TRUE)
  n_full_events <- sum(event_support_summary_ds$support_mode_best == "full_splicing_support", na.rm = TRUE)
  n_other_events <- sum(event_support_summary_ds$support_mode_best == "other_support", na.rm = TRUE)
  n_supported_events <- sum(event_support_summary_ds$support_mode_best != "none", na.rm = TRUE)
  
  pct_general_events <- if (n_total_events > 0) 100 * n_general_events / n_total_events else 0
  pct_acceptor_events <- if (n_total_events > 0) 100 * n_acceptor_events / n_total_events else 0
  pct_donor_events <- if (n_total_events > 0) 100 * n_donor_events / n_total_events else 0
  pct_full_events <- if (n_total_events > 0) 100 * n_full_events / n_total_events else 0
  pct_other_events <- if (n_total_events > 0) 100 * n_other_events / n_total_events else 0
  pct_supported_events <- if (n_total_events > 0) 100 * n_supported_events / n_total_events else 0
  
  n_variants_general <- general_support_tbl_ds %>%
    mutate(variant_key = coalesce(as.character(variant_id_final), as.character(variant_ID), as.character(VariationID))) %>%
    summarise(n = n_distinct(variant_key, na.rm = TRUE)) %>%
    pull(n)
  
  n_variants_supported <- support_tbl_ds %>%
    mutate(variant_key = coalesce(as.character(variant_id_final), as.character(variant_ID), as.character(VariationID))) %>%
    summarise(n = n_distinct(variant_key, na.rm = TRUE)) %>%
    pull(n)
  
  n_variants_acceptor <- support_tbl_ds %>%
    filter(support_mode == "acceptor_support") %>%
    mutate(variant_key = coalesce(as.character(variant_id_final), as.character(variant_ID), as.character(VariationID))) %>%
    summarise(n = n_distinct(variant_key, na.rm = TRUE)) %>%
    pull(n)
  
  n_variants_donor <- support_tbl_ds %>%
    filter(support_mode == "donor_support") %>%
    mutate(variant_key = coalesce(as.character(variant_id_final), as.character(variant_ID), as.character(VariationID))) %>%
    summarise(n = n_distinct(variant_key, na.rm = TRUE)) %>%
    pull(n)
  
  n_variants_full <- support_tbl_ds %>%
    filter(support_mode == "full_splicing_support") %>%
    mutate(variant_key = coalesce(as.character(variant_id_final), as.character(variant_ID), as.character(VariationID))) %>%
    summarise(n = n_distinct(variant_key, na.rm = TRUE)) %>%
    pull(n)
  
  figa_df <- tibble(
    category = factor(
      c(
        "ClinVar variants included",
        "Variants with DS>0 same-exon support",
        "Variants with DS>=0.5 same-exon support",
        "Variants with acceptor support",
        "Variants with donor support",
        "Variants with full splicing support"
      ),
      levels = c(
        "ClinVar variants included",
        "Variants with DS>0 same-exon support",
        "Variants with DS>=0.5 same-exon support",
        "Variants with acceptor support",
        "Variants with donor support",
        "Variants with full splicing support"
      )
    ),
    n = c(
      n_total_variants_included_common,
      n_variants_general,
      n_variants_supported,
      n_variants_acceptor,
      n_variants_donor,
      n_variants_full
    )
  )
  
  figb_df <- tibble(
    support_strength = factor(
      c("DS>0 same-exon support", "DS>=0.5 same-exon support"),
      levels = c("DS>0 same-exon support", "DS>=0.5 same-exon support")
    ),
    n = c(
      n_general_events,
      n_supported_events
    )
  )
  
  figc_df <- tibble(
    support_mode = factor(
      c("Acceptor only support", "Donor only support", "Full splicing support"),
      levels = c("Acceptor only support", "Donor only support", "Full splicing support")
    ),
    n = c(
      n_acceptor_events,
      n_donor_events,
      n_full_events
    )
  )
  
  figd_df <- variant_overall_summary_ds %>%
    filter(!is.na(variant_short), variant_short != "") %>%
    slice_head(n = 20) %>%
    mutate(variant_short = fct_reorder(variant_short, n_unique_pbids))
  
  fige_df <- event_support_summary_ds %>%
    filter(support_mode_best %in% c("acceptor_support", "donor_support", "full_splicing_support")) %>%
    mutate(
      support_mode_plot = case_when(
        support_mode_best == "acceptor_support" ~ "Acceptor only support",
        support_mode_best == "donor_support" ~ "Donor only support",
        support_mode_best == "full_splicing_support" ~ "Full splicing support"
      )
    ) %>%
    select(event_key, transcript_id, support_mode_best, support_mode_plot) %>%
    distinct() %>%
    left_join(pbid_celltype_long, by = "transcript_id", relationship = "many-to-many") %>%
    mutate(cell_type_matrix = replace_na(cell_type_matrix, "Unknown")) %>%
    distinct(event_key, transcript_id, support_mode_best, support_mode_plot, cell_type_matrix) %>%
    count(cell_type_matrix, support_mode_plot, name = "n_event_occurrence")
  
  pa <- ggplot(figa_df, aes(x = category, y = n)) +
    geom_col() +
    geom_text(aes(label = n), vjust = -0.2, size = 4) +
    labs(
      title = paste0("ClinVar variant support overview [", dataset_label, "]"),
      x = NULL,
      y = "Number of ABCA4 ClinVar variants"
    ) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  
  pb <- ggplot(figb_df, aes(x = support_strength, y = n)) +
    geom_col() +
    geom_text(aes(label = n), vjust = -0.2, size = 4) +
    labs(
      title = paste0("Support strength across non-terminal ABCA4 splicing events [", dataset_label, "]"),
      x = NULL,
      y = "Number of events"
    ) +
    theme_bw(base_size = 12)
  
  pc <- ggplot(figc_df, aes(x = support_mode, y = n)) +
    geom_col() +
    geom_text(aes(label = n), vjust = -0.2, size = 4) +
    labs(
      title = paste0("Support mode among DS>=0.5 non-terminal ABCA4 splicing events [", dataset_label, "]"),
      x = NULL,
      y = "Number of events"
    ) +
    theme_bw(base_size = 12)
  
  pd <- ggplot(figd_df, aes(x = variant_short, y = n_unique_pbids)) +
    geom_col() +
    coord_flip() +
    labs(
      title = paste0("Top clinical ABCA4 variants supporting multiple PBIDs [", dataset_label, "]"),
      x = NULL,
      y = "Number of unique PBIDs supported"
    ) +
    theme_bw(base_size = 11)
  
  pe <- ggplot(fige_df, aes(x = fct_reorder(cell_type_matrix, n_event_occurrence, .fun = sum), y = n_event_occurrence, fill = support_mode_plot)) +
    geom_col() +
    coord_flip() +
    labs(
      title = paste0("Cell-type occurrence of DS>=0.5 supported PBID events [", dataset_label, "]"),
      subtitle = "Based on PBID matrix expression > 0",
      x = NULL,
      y = "Number of PBID-event occurrences",
      fill = "Support mode"
    ) +
    theme_bw(base_size = 12)
  
  heat_df <- event_support_summary_ds %>%
    filter(support_mode_best != "none") %>%
    mutate(pbid_event = paste0(transcript_id, " | ", event_type)) %>%
    select(pbid_event, variant_short, maxDS) %>%
    filter(!is.na(variant_short), variant_short != "") %>%
    distinct()
  
  top_pbids <- heat_df %>%
    count(pbid_event, sort = TRUE) %>%
    slice_head(n = 30) %>%
    pull(pbid_event)
  
  top_vars <- heat_df %>%
    count(variant_short, sort = TRUE) %>%
    slice_head(n = 15) %>%
    pull(variant_short)
  
  heat_df2 <- heat_df %>%
    filter(pbid_event %in% top_pbids, variant_short %in% top_vars)
  
  if (nrow(heat_df2) > 0) {
    p_heat <- ggplot(heat_df2, aes(x = variant_short, y = fct_rev(pbid_event), fill = maxDS)) +
      geom_tile() +
      labs(
        title = paste0("Mode-supported PBID landscape [", dataset_label, "]"),
        x = "ClinVar variant",
        y = "PBID event",
        fill = "maxDS"
      ) +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
    
    ggsave(file.path(out_dir, "Supplementary_heatmap_mode_supported_events.pdf"), p_heat, width = 10, height = 9)
    ggsave(file.path(out_dir, "Supplementary_heatmap_mode_supported_events.png"), p_heat, width = 10, height = 9, dpi = 300)
  }
  
  supported_event_celltype_tbl <- event_support_summary_ds %>%
    filter(support_mode_best != "none") %>%
    select(event_key, transcript_id, event_type, support_mode_best, variant_short, maxDS) %>%
    distinct() %>%
    left_join(pbid_celltype_long, by = "transcript_id", relationship = "many-to-many") %>%
    arrange(transcript_id, cell_type_matrix)
  
  summary_table <- tibble(
    dataset = dataset_label,
    total_nonterminal_events = n_total_events,
    general_supported_events = n_general_events,
    pct_general_supported_events = pct_general_events,
    mode_supported_events = n_supported_events,
    pct_mode_supported_events = pct_supported_events,
    acceptor_supported_events = n_acceptor_events,
    pct_acceptor_supported_events = pct_acceptor_events,
    donor_supported_events = n_donor_events,
    pct_donor_supported_events = pct_donor_events,
    full_supported_events = n_full_events,
    pct_full_supported_events = pct_full_events,
    other_supported_events = n_other_events,
    pct_other_supported_events = pct_other_events,
    total_clinvar_variants_included = n_total_variants_included_common,
    variants_with_general_support = n_variants_general,
    variants_with_any_mode_support = n_variants_supported,
    variants_with_acceptor_support = n_variants_acceptor,
    variants_with_donor_support = n_variants_donor,
    variants_with_full_support = n_variants_full
  )
  
  diagnostic_lines <- c(
    paste0("ABCA4 support analysis diagnostics [", dataset_label, "]"),
    paste0(strrep("=", 38 + nchar(dataset_label))),
    "",
    paste("Matrix file:", matrix_file),
    paste("Number of PBIDs in matrix:", length(pbids_in_dataset)),
    paste("Rows in subsetted matches:", nrow(matches_ds)),
    paste("Rows in subsetted events:", nrow(events_ds)),
    paste("Unique non-terminal events:", n_total_events),
    paste("Matrix-derived cell-type columns:", paste(celltype_cols, collapse = ", "))
  )
  
  report_lines <- c(
    paste0("ABCA4 support summary [", dataset_label, "]"),
    paste0(strrep("=", 54 + nchar(dataset_label))),
    "",
    paste("Matrix file:", matrix_file),
    paste("PBID subset definition: PBIDs present in", basename(matrix_file)),
    "",
    paste("Total non-terminal PBID events:", n_total_events),
    paste("DS>0 same-exon supported events:", n_general_events, sprintf("(%.1f%%)", pct_general_events)),
    paste("DS>=0.5 same-exon supported events:", n_supported_events, sprintf("(%.1f%%)", pct_supported_events)),
    paste("Acceptor-supported events:", n_acceptor_events, sprintf("(%.1f%%)", pct_acceptor_events)),
    paste("Donor-supported events:", n_donor_events, sprintf("(%.1f%%)", pct_donor_events)),
    paste("Full-splicing-supported events:", n_full_events, sprintf("(%.1f%%)", pct_full_events)),
    paste("Other-supported events:", n_other_events, sprintf("(%.1f%%)", pct_other_events)),
    "",
    paste("ClinVar variants included:", n_total_variants_included_common),
    paste("Variants with DS>0 same-exon support:", n_variants_general),
    paste("Variants with DS>=0.5 same-exon support:", n_variants_supported),
    paste("Variants with acceptor support:", n_variants_acceptor),
    paste("Variants with donor support:", n_variants_donor),
    paste("Variants with full support:", n_variants_full),
    "",
    "Notes:",
    "- General support = exact relevant positional concordance + DS > 0.",
    "- Mode-defined support = exact relevant positional concordance + DS >= 0.5.",
    "- Terminal exon events are excluded from all main statistics.",
    paste0("- This analysis is explicitly restricted to PBIDs present in ", basename(matrix_file), ".")
  )
  
  write_tsv(general_support_tbl_ds, file.path(out_dir, "abca4_general_supported_matches.tsv"))
  write_tsv(support_tbl_ds, file.path(out_dir, "abca4_mode_supported_matches.tsv"))
  write_tsv(best_general_per_event_ds, file.path(out_dir, "abca4_best_general_support_per_event.tsv"))
  write_tsv(best_support_per_event_ds, file.path(out_dir, "abca4_best_mode_support_per_event.tsv"))
  write_tsv(event_support_summary_ds, file.path(out_dir, "abca4_event_support_summary.tsv"))
  write_tsv(variant_overall_summary_ds, file.path(out_dir, "abca4_variant_overall_support_summary.tsv"))
  
  write_tsv(figa_df, file.path(out_dir, "plot_data_FigureA.tsv"))
  write_tsv(figb_df, file.path(out_dir, "plot_data_FigureB.tsv"))
  write_tsv(figc_df, file.path(out_dir, "plot_data_FigureC.tsv"))
  write_tsv(figd_df, file.path(out_dir, "plot_data_FigureD.tsv"))
  write_tsv(fige_df, file.path(out_dir, "plot_data_FigureE.tsv"))
  
  write_tsv(supported_event_celltype_tbl, file.path(out_dir, "Supplementary_mode_supported_event_celltypes.tsv"))
  write_tsv(summary_table, file.path(out_dir, "summary_table_for_writing.tsv"))
  
  writeLines(diagnostic_lines, file.path(out_dir, "ABCA4_support_diagnostics.txt"))
  writeLines(report_lines, file.path(out_dir, "ABCA4_support_overall_summary.txt"))
  
  ggsave(file.path(out_dir, "FigureA_clinvar_support_overview.pdf"), p1a, width = 9, height = 5)
  ggsave(file.path(out_dir, "FigureA_clinvar_support_overview.png"), p1a, width = 9, height = 5, dpi = 300)
  
  ggsave(file.path(out_dir, "FigureB_support_strength.pdf"), p1b, width = 9, height = 5)
  ggsave(file.path(out_dir, "FigureB_support_strength.png"), p1b, width = 9, height = 5, dpi = 300)
  
  ggsave(file.path(out_dir, "FigureC_support_mode_distribution.pdf"), p1c, width = 8, height = 6)
  ggsave(file.path(out_dir, "FigureC_support_mode_distribution.png"), p1c, width = 8, height = 6, dpi = 300)
  
  ggsave(file.path(out_dir, "FigureD_top_variants_supporting_pbids.pdf"), p1d, width = 9, height = 7)
  ggsave(file.path(out_dir, "FigureD_top_variants_supporting_pbids.png"), p1d, width = 9, height = 7, dpi = 300)
  
  ggsave(file.path(out_dir, "FigureE_celltype_distribution_mode_supported_events.pdf"), p1e, width = 8, height = 6)
  ggsave(file.path(out_dir, "FigureE_celltype_distribution_mode_supported_events.png"), p1e, width = 8, height = 6, dpi = 300)
  
  list(
    dataset_label = dataset_label,
    summary_table = summary_table,
    event_support_summary = event_support_summary_ds,
    general_support_tbl = general_support_tbl_ds,
    support_tbl = support_tbl_ds
  )
}
#use FigureC from this part as the figure 6E

#read common input files once
events_raw <- read_tsv(events_file, show_col_types = FALSE)
matches_raw <- read_tsv(matches_file, show_col_types = FALSE)
isoform_summary_raw <- read_tsv(isoform_summary_file, show_col_types = FALSE)
clinvar_included <- read_tsv(clinvar_included_file, show_col_types = FALSE)

#clean common tables

num_cols_matches <- c(
  "pb_exon_index", "ref_exon_index",
  "observed_start", "observed_end",
  "observed_donor_pos", "observed_acceptor_pos",
  "reference_donor_pos", "reference_acceptor_pos",
  "variant_POS",
  "DS_AG", "DS_AL", "DS_DG", "DS_DL",
  "DP_AG", "DP_AL", "DP_DG", "DP_DL",
  "pred_AG_pos", "pred_AL_pos", "pred_DG_pos", "pred_DL_pos",
  "maxDS",
  "diff_AG_obs_acceptor", "diff_DG_obs_donor",
  "diff_AL_ref_acceptor", "diff_DL_ref_donor",
  "score_AG_obs_acceptor", "score_DG_obs_donor",
  "score_AL_ref_acceptor", "score_DL_ref_donor",
  "match_score"
)

for (cc in intersect(num_cols_matches, colnames(matches_raw))) {
  matches_raw[[cc]] <- safe_num(matches_raw[[cc]])
}

for (cc in intersect(c("pb_exon_index", "ref_exon_index", "observed_start", "observed_end"), colnames(events_raw))) {
  events_raw[[cc]] <- safe_num(events_raw[[cc]])
}

if ("cell_type" %in% colnames(matches_raw)) {
  matches_raw <- matches_raw %>% mutate(cell_type = na_if(cell_type, ""))
}
if ("cell_type" %in% colnames(events_raw)) {
  events_raw <- events_raw %>% mutate(cell_type = na_if(cell_type, ""))
}
if ("cell_type" %in% colnames(isoform_summary_raw)) {
  isoform_summary_raw <- isoform_summary_raw %>% mutate(cell_type = na_if(cell_type, ""))
}

events_raw <- event_key_fun(events_raw)
matches_raw <- event_key_fun(matches_raw)

isoform_meta <- isoform_summary_raw %>%
  group_by(transcript_id) %>%
  summarise(
    isoform_cell_type_label = first_non_na(cell_type),
    major_event = first_non_na(major_event),
    priority_for_followup = first_non_na(priority_for_followup),
    exon_count = first_non_na_num(exon_count),
    junction_count = first_non_na_num(junction_count),
    .groups = "drop"
  )

clinvar_included2 <- clinvar_included %>%
  mutate(
    variant_id_final = coalesce(
      get_col_if_exists(., "ID"),
      get_col_if_exists(., "VariationID")
    )
  )

n_total_variants_included_common <- clinvar_included2 %>%
  summarise(n = n_distinct(variant_id_final, na.rm = TRUE)) %>%
  pull(n)

#run the analysis
dataset_configs <- tibble::tribble(
  ~dataset_label,            ~matrix_file,                       ~out_dir,
  "FL2_CAGE_polyA",         matrix_file_FL2_CAGE_polyA,         file.path(out_dir_base, "FL2_CAGE_polyA")
)

analysis_results <- purrr::pmap(
  dataset_configs,
  function(dataset_label, matrix_file, out_dir) {
    run_support_analysis(
      dataset_label = dataset_label,
      matrix_file = matrix_file,
      out_dir = out_dir,
      matches_raw = matches_raw,
      events_raw = events_raw,
      isoform_meta = isoform_meta,
      n_total_variants_included_common = n_total_variants_included_common
    )
  }
)
#use FigureC from this part as the figure 6E
names(analysis_results) <- dataset_configs$dataset_label

#save the combined summary
combined_summary_table <- bind_rows(lapply(analysis_results, function(x) x$summary_table))
write_tsv(combined_summary_table, file.path(out_dir_base, "combined_dataset_summary.tsv"))
for (nm in names(analysis_results)) {
  cat("\n==============================\n")
  cat("Dataset:", nm, "\n")
  cat("==============================\n")
  print(analysis_results[[nm]]$summary_table)
  print(count(analysis_results[[nm]]$event_support_summary, general_support))
  print(count(analysis_results[[nm]]$event_support_summary, support_mode_best))
}
print(combined_summary_table)


#select full_splicing_support PBids to show on the transcript map
#settings
base_dir <- "D:/bioinformatic/isoseq/processed/from_Ray"
out_dir <- file.path(base_dir, "revision/splicAI/ABCA4_spliceAI_visual_selectedPBids")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

target_pbids <- c(
  "PB.5877.586",
  "PB.5877.634",
  "PB.5877.740",
  "PB.5877.800",
  "PB.5877.898"
)

pbid_label_map <- tibble(
  transcript_id = target_pbids,
  display_group = c(
    "Rod",
    "Rod+MG",
    "Bipolar+Rod",
    "RPE+Cone+Rod",
    "Amacrine+SMC+Rod"
  )
)

min_ds_to_plot <- 0.20
max_events_per_transcript_class <- 6

zoom_region <- NULL

fig_width  <- 30
fig_height <- 12
fig_dpi    <- 300

#define functions
safe_import_gff <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  rtracklayer::import(path) %>% as_tibble()
}

safe_fread <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  data.table::fread(path) %>% as_tibble()
}

save_plot_pdf_png <- function(plot_obj, prefix, out_dir,
                              width = 30, height = 12, dpi = 300) {
  ggsave(
    filename = file.path(out_dir, paste0(prefix, ".pdf")),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    limitsize = FALSE
  )
  ggsave(
    filename = file.path(out_dir, paste0(prefix, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    limitsize = FALSE
  )
}

make_introns_from_exons <- function(exon_df) {
  #exon_df contain: transcript_to_plot, start, end, strand
  exon_df %>%
    arrange(transcript_to_plot, start, end) %>%
    group_by(transcript_to_plot) %>%
    mutate(
      intron_start = lag(end),
      intron_end   = start
    ) %>%
    ungroup() %>%
    filter(!is.na(intron_start), !is.na(intron_end), intron_end > intron_start) %>%
    transmute(
      transcript_to_plot,
      strand,
      start = intron_start,
      end   = intron_end
    )
}

#file paths
pb_gff_path <- file.path(
  base_dir,
  "original_from_Ray/ray_pigeon_output/classify.retina.merge.GRCh38.p14.with.cage.polya.FL/RetinaMerge.sorted.filtered_lite.gff"
)

ref_gff_path <- file.path(
  base_dir,
  "GRCh38p14/GRCh38.p14.v44/gencode.v44.chr_patch_hapl_scaff.annotation.gff3"
)

splice_match_path <- file.path(
  base_dir,
  "revision/splicAI/ABCA4_support_overall_results_v2/FL2_CAGE_polyA/abca4_mode_supported_matches.tsv"
)

#load data
pb_anno <- safe_import_gff(pb_gff_path)
reference_anno <- safe_import_gff(ref_gff_path)
splice_match <- safe_fread(splice_match_path)

#build exon annotation
pb_exons <- pb_anno %>%
  filter(
    !is.na(transcript_id),
    transcript_id %in% target_pbids,
    type == "exon"
  ) %>%
  select(seqnames, start, end, strand, transcript_id) %>%
  left_join(pbid_label_map, by = "transcript_id") %>%
  mutate(
    transcript_to_plot = paste0(display_group, " | ", transcript_id),
    track_type = "PBid"
  )

if (nrow(pb_exons) == 0) {
  stop("No PB exons found for selected PBids.")
}

ref_exons <- reference_anno %>%
  filter(
    !is.na(gene_name),
    gene_name == "ABCA4",
    type == "exon",
    transcript_name == "ABCA4-201"
  ) %>%
  transmute(
    seqnames = seqnames,
    start = start,
    end = end,
    strand = strand,
    transcript_id = "ABCA4-201",
    display_group = "Reference",
    transcript_to_plot = "Reference | ABCA4-201",
    track_type = "Reference"
  )

if (nrow(ref_exons) == 0) {
  stop("Reference ABCA4-201 exons not found.")
}

plot_exons <- bind_rows(ref_exons, pb_exons)

#plot range
full_plot_start <- min(plot_exons$start, na.rm = TRUE)
full_plot_end   <- max(plot_exons$end,   na.rm = TRUE)

if (is.null(zoom_region)) {
  plot_start <- full_plot_start
  plot_end   <- full_plot_end
} else {
  plot_start <- zoom_region[1]
  plot_end   <- zoom_region[2]
}

plot_exons <- plot_exons %>%
  filter(end >= plot_start, start <= plot_end)

if (nrow(plot_exons) == 0) {
  stop("plot_exons became empty after plot range filtering.")
}

#prepare SpliceAI data
splice_sel <- splice_match %>%
  filter(transcript_id %in% target_pbids) %>%
  left_join(pbid_label_map, by = "transcript_id") %>%
  mutate(
    transcript_to_plot = paste0(display_group, " | ", transcript_id)
  )

if (nrow(splice_sel) == 0) {
  stop("No matching rows found in SpliceAI table for selected PBids.")
}

if (!"variant_ID" %in% colnames(splice_sel)) {
  splice_sel <- splice_sel %>%
    mutate(variant_ID = variant_short)
}

splice_long <- splice_sel %>%
  select(
    transcript_id,
    transcript_to_plot,
    display_group,
    variant_short,
    variant_ID,
    variant_POS,
    pred_AG_pos, pred_AL_pos, pred_DG_pos, pred_DL_pos,
    DS_AG, DS_AL, DS_DG, DS_DL
  ) %>%
  pivot_longer(
    cols = c(pred_AG_pos, pred_AL_pos, pred_DG_pos, pred_DL_pos),
    names_to = "pred_type",
    values_to = "pred_pos"
  ) %>%
  mutate(
    ds = case_when(
      pred_type == "pred_AG_pos" ~ DS_AG,
      pred_type == "pred_AL_pos" ~ DS_AL,
      pred_type == "pred_DG_pos" ~ DS_DG,
      pred_type == "pred_DL_pos" ~ DS_DL,
      TRUE ~ NA_real_
    ),
    splice_class = case_when(
      pred_type == "pred_AG_pos" ~ "AG",
      pred_type == "pred_AL_pos" ~ "AL",
      pred_type == "pred_DG_pos" ~ "DG",
      pred_type == "pred_DL_pos" ~ "DL",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(pred_pos), !is.na(ds), ds >= min_ds_to_plot) %>%
  filter(
    (variant_POS >= plot_start & variant_POS <= plot_end) |
      (pred_pos >= plot_start & pred_pos <= plot_end)
  )

if (nrow(splice_long) == 0) {
  stop("No SpliceAI events remained after filtering. Try lowering min_ds_to_plot.")
}

splice_long <- splice_long %>%
  group_by(transcript_to_plot, splice_class, variant_short, variant_POS, pred_pos) %>%
  slice_max(order_by = ds, n = 1, with_ties = FALSE) %>%
  ungroup()

splice_long <- splice_long %>%
  group_by(transcript_to_plot, splice_class) %>%
  arrange(desc(ds), .by_group = TRUE) %>%
  slice_head(n = max_events_per_transcript_class) %>%
  ungroup() %>%
  mutate(
    score_label = paste0(splice_class, "=", sprintf("%.2f", ds))
  )

#variant track
variant_track_name <- "ClinVar variants"

variant_ref <- splice_long %>%
  distinct(variant_short, variant_ID, variant_POS) %>%
  mutate(
    transcript_to_plot = variant_track_name,
    variant_label = variant_short
  )

#Y mapping
pbid_rows <- pbid_label_map %>%
  mutate(transcript_to_plot = paste0(display_group, " | ", transcript_id)) %>%
  pull(transcript_to_plot)

y_levels <- c(
  variant_track_name,
  "Reference | ABCA4-201",
  pbid_rows
)

y_map <- tibble(
  transcript_to_plot = y_levels,
  y = rev(seq_along(y_levels))
)

plot_exons <- plot_exons %>%
  left_join(y_map, by = "transcript_to_plot")

variant_ref <- variant_ref %>%
  left_join(y_map, by = "transcript_to_plot")

splice_long <- splice_long %>%
  left_join(y_map, by = "transcript_to_plot")

if (any(is.na(plot_exons$y))) stop("Some exons failed to map to y positions.")
if (any(is.na(variant_ref$y))) stop("Some variant track rows failed to map to y positions.")
if (any(is.na(splice_long$y))) stop("Some splice events failed to map to y positions.")

#event offsets
splice_long <- splice_long %>%
  mutate(
    y_offset = case_when(
      splice_class == "AG" ~  0.22,
      splice_class == "AL" ~  0.10,
      splice_class == "DG" ~ -0.10,
      splice_class == "DL" ~ -0.22,
      TRUE ~ 0
    ),
    y_event = y + y_offset
  )

#introns
plot_introns <- make_introns_from_exons(plot_exons) %>%
  left_join(y_map, by = "transcript_to_plot")

#aesthetics
splice_color_map <- c(
  "AG" = "#1f78b4",
  "AL" = "#6baed6",
  "DG" = "#d95f02",
  "DL" = "#fdae6b"
)

splice_shape_map <- c(
  "AG" = 24,
  "AL" = 25,
  "DG" = 24,
  "DL" = 25
)

#exon box thickness
exon_half_height <- 0.10

#plot
p_main <- ggplot() +
  #vertical lines for ClinVar positions
  geom_vline(
    data = variant_ref,
    aes(xintercept = variant_POS),
    color = "grey82",
    linewidth = 0.35,
    linetype = "dashed"
  ) +
  
  #exon boxes
  geom_rect(
    data = plot_exons,
    aes(
      xmin = start,
      xmax = end,
      ymin = y - exon_half_height,
      ymax = y + exon_half_height
    ),
    fill = "black",
    color = "black"
  ) +
  
  #intron lines
  geom_segment(
    data = plot_introns,
    aes(
      x = start,
      xend = end,
      y = y,
      yend = y
    ),
    color = "grey40",
    linewidth = 0.6
  ) +
  
  #ClinVar variant points
  geom_point(
    data = variant_ref,
    aes(x = variant_POS, y = y),
    color = "red3",
    size = 4
  ) +
  
  #ClinVar variant labels
  ggrepel::geom_text_repel(
    data = variant_ref,
    aes(x = variant_POS, y = y, label = variant_label),
    color = "red4",
    size = 5.0,
    direction = "y",
    nudge_y = 0.45,
    box.padding = 0.45,
    point.padding = 0.25,
    segment.size = 0.30,
    min.segment.length = 0,
    max.overlaps = 300
  ) +
  
  #spliceAI predicted sites on PBids
  geom_point(
    data = splice_long,
    aes(
      x = pred_pos,
      y = y_event,
      fill = splice_class,
      shape = splice_class
    ),
    size = 4.5,
    color = "black",
    stroke = 0.5
  ) +
  
  #spliceAI score labels
  ggrepel::geom_text_repel(
    data = splice_long,
    aes(
      x = pred_pos,
      y = y_event,
      label = score_label,
      color = splice_class
    ),
    size = 4.8,
    box.padding = 0.40,
    point.padding = 0.20,
    segment.size = 0.28,
    min.segment.length = 0,
    max.overlaps = 300
  ) +
  
  scale_fill_manual(values = splice_color_map, drop = FALSE) +
  scale_color_manual(values = splice_color_map, drop = FALSE) +
  scale_shape_manual(values = splice_shape_map, drop = FALSE) +
  
  scale_y_continuous(
    breaks = y_map$y,
    labels = y_map$transcript_to_plot,
    expand = expansion(mult = c(0.08, 0.12))
  ) +
  
  coord_cartesian(
    xlim = c(plot_start, plot_end),
    clip = "off"
  ) +
  
  labs(
    title = "ABCA4 selected Iso-Seq PBids with SpliceAI-supported splice events",
    subtitle = "Top track: ClinVar variant IDs; transcript rows: predicted acceptor/donor gain/loss scores",
    x = "Genomic position",
    y = NULL,
    fill = "SpliceAI class",
    color = "SpliceAI class",
    shape = "SpliceAI class"
  ) +
  
  theme_bw(base_size = 18) +
  theme(
    text = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    plot.title = element_text(size = 22, face = "bold"),
    plot.subtitle = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(15, 25, 15, 15)
  )


#Save outputs

#figure 6F
save_plot_pdf_png(
  plot_obj = p_main,
  prefix = "ABCA4_selectedPBids_spliceAI_overlay",
  out_dir = out_dir,
  width = fig_width,
  height = fig_height,
  dpi = fig_dpi
)

#statistics
write.csv(plot_exons,   file.path(out_dir, "plot_exons_selectedPBids.csv"), row.names = FALSE)
write.csv(plot_introns, file.path(out_dir, "plot_introns_selectedPBids.csv"), row.names = FALSE)
write.csv(splice_long,  file.path(out_dir, "spliceAI_events_selectedPBids.csv"), row.names = FALSE)
write.csv(variant_ref,  file.path(out_dir, "spliceAI_variants_reference_track.csv"), row.names = FALSE)

run_summary <- tibble(
  item = c(
    "output_directory",
    "n_target_pbids",
    "n_plot_exons",
    "n_plot_introns",
    "n_splice_events",
    "n_unique_variants",
    "plot_start",
    "plot_end",
    "min_ds_to_plot",
    "max_events_per_transcript_class"
  ),
  value = c(
    out_dir,
    length(target_pbids),
    nrow(plot_exons),
    nrow(plot_introns),
    nrow(splice_long),
    nrow(variant_ref),
    plot_start,
    plot_end,
    min_ds_to_plot,
    max_events_per_transcript_class
  )
)

write.csv(run_summary, file.path(out_dir, "run_summary.csv"), row.names = FALSE)


#perform RT-PCR to amplify specific isoforms from the total RNA extracted from the human retina
#send cDNAs for Sanger sequencing
#align the sequencing results back to the transcript
#validate the selected isoforms
#path
base_dir <- "D:/bioinformatic/isoseq/processed/from_Ray"
pb_gff_path <- file.path(
  base_dir,
  "original_from_Ray/ray_pigeon_output/classify.retina.merge.GRCh38.p14.with.cage.polya.FL/RetinaMerge.sorted.filtered_lite.gff"
)
ref_gff_path <- file.path(
  base_dir,
  "GRCh38p14/GRCh38.p14.v44/gencode.v44.chr_patch_hapl_scaff.annotation.gff3"
)
reference_gene_name <- "ABCA4"
reference_transcript_name <- "ABCA4-201"
fasta_dir <- file.path(
  base_dir,
  "revision/RT_PCR/merged_mapping_fasta"
)
primer_seq_csv_path <- file.path(
  fasta_dir,
  "primer_sequences.csv"
)
out_dir <- file.path(
  base_dir,
  "revision/RT_PCR/merged_mapping_plots"
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
target_pbids <- c(
  "PB.5877.945",
  "PB.5877.897"
)

fasta_suffix <- "_map_result.fasta"

fig_width <- 16
fig_height <- 6.6
fig_dpi <- 300

plot_pad_bp <- 80

zoom_to_mapping_region <- TRUE

use_compressed_x <- TRUE
compress_gap_threshold_bp <- 500
compressed_gap_width <- 120
compress_feature_flank_bp <- 20

reference_exon_half_height <- 0.22
pbid_exon_half_height <- 0.24
sanger_half_height <- 0.18

show_primers <- TRUE
primer_arrow_size <- 0.22
primer_label_size <- 3.2

#FASTA file
read_fasta_as_strings <- function(path) {
  if (!file.exists(path)) {
    stop("FASTA file not found: ", path)
  }
  
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[lines != ""]
  
  header_idx <- which(startsWith(lines, ">"))
  
  if (length(header_idx) == 0) {
    stop("No FASTA headers found in: ", path)
  }
  
  res <- vector("list", length(header_idx))
  
  for (i in seq_along(header_idx)) {
    start_i <- header_idx[i]
    end_i <- if (i < length(header_idx)) {
      header_idx[i + 1] - 1
    } else {
      length(lines)
    }
    
    id <- sub("^>", "", lines[start_i])
    seq_lines <- lines[(start_i + 1):end_i]
    seq <- paste(seq_lines, collapse = "")
    seq <- toupper(gsub("\\s+", "", seq))
    
    res[[i]] <- tibble(
      id = id,
      seq = seq,
      nchar = nchar(seq)
    )
  }
  
  df <- bind_rows(res)
  
  if (length(unique(df$nchar)) != 1) {
    stop(
      "Not all sequences have the same alignment length in: ",
      path,
      "\nThis script requires an aligned multi-FASTA."
    )
  }
  
  df
}

#Annotation
get_pbid_exons <- function(pb_gff_path, pbid) {
  if (!file.exists(pb_gff_path)) {
    stop("PBid GFF not found: ", pb_gff_path)
  }
  
  gff <- rtracklayer::import(pb_gff_path) %>%
    as_tibble()
  
  exons <- gff %>%
    filter(
      type == "exon",
      !is.na(transcript_id),
      as.character(transcript_id) == pbid
    ) %>%
    transmute(
      seqnames = as.character(seqnames),
      start = as.numeric(start),
      end = as.numeric(end),
      strand = as.character(strand),
      transcript_id = as.character(transcript_id)
    )
  
  if (nrow(exons) == 0) {
    stop("No exon annotation found for ", pbid, " in PBid GFF.")
  }
  
  exons
}

get_reference_exons <- function(ref_gff_path,
                                gene_name_target = "ABCA4",
                                transcript_name_target = "ABCA4-201") {
  if (!file.exists(ref_gff_path)) {
    stop("Reference GFF not found: ", ref_gff_path)
  }
  
  gff <- rtracklayer::import(ref_gff_path) %>%
    as_tibble()
  
  exons <- gff %>%
    filter(
      type == "exon",
      !is.na(gene_name),
      gene_name == gene_name_target,
      !is.na(transcript_name),
      transcript_name == transcript_name_target
    ) %>%
    transmute(
      seqnames = as.character(seqnames),
      start = as.numeric(start),
      end = as.numeric(end),
      strand = as.character(strand),
      transcript_id = as.character(transcript_name)
    )
  
  if (nrow(exons) == 0) {
    stop(
      "No reference exon annotation found for ",
      gene_name_target,
      " / ",
      transcript_name_target,
      " in reference GFF."
    )
  }
  
  exons
}

add_tx_coordinates <- function(exons) {
  strand_now <- unique(exons$strand)
  
  if (length(strand_now) != 1) {
    stop("More than one strand found for this transcript.")
  }
  
  if (strand_now == "-") {
    exons <- exons %>%
      arrange(desc(start), desc(end))
  } else {
    exons <- exons %>%
      arrange(start, end)
  }
  
  exons %>%
    mutate(
      exon_width = end - start + 1,
      tx_start = cumsum(lag(exon_width, default = 0)) + 1,
      tx_end = tx_start + exon_width - 1,
      exon_order = row_number()
    )
}

make_introns_from_exons <- function(exons_tx) {
  if (is.null(exons_tx) || nrow(exons_tx) <= 1) {
    return(tibble(start = numeric(), end = numeric()))
  }
  
  exons_tx %>%
    arrange(exon_order) %>%
    mutate(
      intron_start = lag(end),
      intron_end = start
    ) %>%
    filter(!is.na(intron_start), !is.na(intron_end)) %>%
    transmute(
      start = pmin(intron_start, intron_end),
      end = pmax(intron_start, intron_end)
    )
}

#coordinate conversion
tx_interval_to_genomic_blocks <- function(tx_start, tx_end, exons_tx) {
  strand_now <- unique(exons_tx$strand)
  
  res <- lapply(seq_len(nrow(exons_tx)), function(i) {
    ex <- exons_tx[i, ]
    
    ov_start <- max(tx_start, ex$tx_start)
    ov_end <- min(tx_end, ex$tx_end)
    
    if (ov_start > ov_end) {
      return(NULL)
    }
    
    if (strand_now == "+") {
      g_start <- ex$start + (ov_start - ex$tx_start)
      g_end <- ex$start + (ov_end - ex$tx_start)
    } else {
      g_start <- ex$end - (ov_end - ex$tx_start)
      g_end <- ex$end - (ov_start - ex$tx_start)
    }
    
    tibble(
      seqnames = ex$seqnames,
      start = min(g_start, g_end),
      end = max(g_start, g_end),
      strand = ex$strand,
      tx_start = ov_start,
      tx_end = ov_end,
      exon_order = ex$exon_order
    )
  })
  
  bind_rows(res)
}

tx_interval_to_one_genomic_range <- function(tx_start, tx_end, exons_tx) {
  blocks <- tx_interval_to_genomic_blocks(
    tx_start = tx_start,
    tx_end = tx_end,
    exons_tx = exons_tx
  )
  
  if (nrow(blocks) == 0) {
    return(tibble(
      seqnames = NA_character_,
      start = NA_real_,
      end = NA_real_,
      strand = NA_character_
    ))
  }
  
  tibble(
    seqnames = paste(unique(blocks$seqnames), collapse = ";"),
    start = min(blocks$start, na.rm = TRUE),
    end = max(blocks$end, na.rm = TRUE),
    strand = paste(unique(blocks$strand), collapse = ";")
  )
}

#compressed x-axis
#this step is to better visualize the results
#because exon is quite shorter than introns especially for ABCA4, so shorten the introns for better visualization
merge_intervals_for_plot <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(tibble(start = numeric(), end = numeric()))
  }
  
  x <- df %>%
    transmute(
      start = as.numeric(start),
      end = as.numeric(end)
    ) %>%
    filter(!is.na(start), !is.na(end)) %>%
    mutate(
      start = pmin(start, end),
      end = pmax(start, end)
    ) %>%
    arrange(start, end)
  
  if (nrow(x) == 0) {
    return(tibble(start = numeric(), end = numeric()))
  }
  
  out <- list()
  cur_start <- x$start[1]
  cur_end <- x$end[1]
  k <- 1
  
  if (nrow(x) >= 2) {
    for (i in 2:nrow(x)) {
      s <- x$start[i]
      e <- x$end[i]
      
      if (s <= cur_end + 1) {
        cur_end <- max(cur_end, e)
      } else {
        out[[k]] <- tibble(start = cur_start, end = cur_end)
        k <- k + 1
        cur_start <- s
        cur_end <- e
      }
    }
  }
  
  out[[k]] <- tibble(start = cur_start, end = cur_end)
  
  bind_rows(out)
}

build_compressed_axis <- function(feature_df,
                                  flank_bp = 20,
                                  gap_threshold_bp = 500,
                                  compressed_gap_width = 120) {
  features <- feature_df %>%
    transmute(
      start = as.numeric(start) - flank_bp,
      end = as.numeric(end) + flank_bp
    ) %>%
    filter(!is.na(start), !is.na(end)) %>%
    mutate(
      start = pmin(start, end),
      end = pmax(start, end)
    ) %>%
    arrange(start, end)
  
  features_merged <- merge_intervals_for_plot(features)
  
  if (nrow(features_merged) == 0) {
    stop("No features available for compressed axis.")
  }
  
  segments <- vector("list", nrow(features_merged))
  
  for (i in seq_len(nrow(features_merged))) {
    orig_start <- features_merged$start[i]
    orig_end <- features_merged$end[i]
    
    if (i == 1) {
      plot_start <- orig_start
    } else {
      previous_orig_end <- features_merged$end[i - 1]
      previous_plot_end <- segments[[i - 1]]$plot_end
      
      original_gap <- orig_start - previous_orig_end
      
      if (original_gap > gap_threshold_bp) {
        plot_gap <- compressed_gap_width
      } else {
        plot_gap <- original_gap
      }
      
      plot_start <- previous_plot_end + plot_gap
    }
    
    feature_width <- orig_end - orig_start
    plot_end <- plot_start + feature_width
    
    segments[[i]] <- tibble(
      orig_start = orig_start,
      orig_end = orig_end,
      plot_start = plot_start,
      plot_end = plot_end
    )
  }
  
  bind_rows(segments)
}

map_genomic_to_compressed <- function(x, axis_tbl) {
  sapply(x, function(xx) {
    if (is.na(xx)) {
      return(NA_real_)
    }
    
    inside_idx <- which(xx >= axis_tbl$orig_start & xx <= axis_tbl$orig_end)
    
    if (length(inside_idx) >= 1) {
      i <- inside_idx[1]
      return(axis_tbl$plot_start[i] + (xx - axis_tbl$orig_start[i]))
    }
    
    before_idx <- which(axis_tbl$orig_end < xx)
    after_idx <- which(axis_tbl$orig_start > xx)
    
    if (length(before_idx) == 0) {
      i <- 1
      return(axis_tbl$plot_start[i] - (axis_tbl$orig_start[i] - xx))
    }
    
    if (length(after_idx) == 0) {
      i <- nrow(axis_tbl)
      return(axis_tbl$plot_end[i] + (xx - axis_tbl$orig_end[i]))
    }
    
    left_i <- max(before_idx)
    right_i <- min(after_idx)
    
    orig_gap_start <- axis_tbl$orig_end[left_i]
    orig_gap_end <- axis_tbl$orig_start[right_i]
    plot_gap_start <- axis_tbl$plot_end[left_i]
    plot_gap_end <- axis_tbl$plot_start[right_i]
    
    frac <- (xx - orig_gap_start) / (orig_gap_end - orig_gap_start)
    
    plot_gap_start + frac * (plot_gap_end - plot_gap_start)
  })
}

add_plot_x <- function(df, axis_tbl) {
  if (is.null(df) || nrow(df) == 0) {
    return(df)
  }
  
  df %>%
    mutate(
      plot_start = map_genomic_to_compressed(start, axis_tbl),
      plot_end = map_genomic_to_compressed(end, axis_tbl),
      plot_xmin = pmin(plot_start, plot_end),
      plot_xmax = pmax(plot_start, plot_end)
    )
}

#Primer alignment
reverse_complement_dna <- function(seq) {
  seq <- toupper(gsub("\\s+", "", seq))
  comp <- chartr("ACGTNacgtn", "TGCANtgcan", seq)
  paste(rev(strsplit(comp, "")[[1]]), collapse = "")
}

read_pbid_sequence_from_alignment_fasta <- function(fasta_path, pbid) {
  aln_df <- read_fasta_as_strings(fasta_path)
  
  ref_idx <- which(stringr::str_detect(aln_df$id, stringr::fixed(pbid)))
  
  if (length(ref_idx) != 1) {
    stop(
      "Expected exactly one PBid reference sequence containing ",
      pbid,
      " in ",
      fasta_path,
      ", but found ",
      length(ref_idx)
    )
  }
  
  seq <- aln_df$seq[ref_idx]
  seq <- toupper(gsub("-", "", seq))
  seq
}

find_primer_in_pbid_sequence <- function(pbid_seq, primer_seq) {
  pbid_seq <- toupper(gsub("\\s+", "", pbid_seq))
  primer_seq <- toupper(gsub("\\s+", "", primer_seq))
  primer_rc <- reverse_complement_dna(primer_seq)
  
  hit_fwd <- gregexpr(primer_seq, pbid_seq, fixed = TRUE)[[1]]
  hit_rev <- gregexpr(primer_rc, pbid_seq, fixed = TRUE)[[1]]
  
  hit_fwd <- hit_fwd[hit_fwd > 0]
  hit_rev <- hit_rev[hit_rev > 0]
  
  hits <- tibble()
  
  if (length(hit_fwd) > 0) {
    hits <- bind_rows(
      hits,
      tibble(
        match_orientation = "primer_sequence",
        tx_start = as.numeric(hit_fwd),
        tx_end = as.numeric(hit_fwd + nchar(primer_seq) - 1)
      )
    )
  }
  
  if (length(hit_rev) > 0) {
    hits <- bind_rows(
      hits,
      tibble(
        match_orientation = "reverse_complement",
        tx_start = as.numeric(hit_rev),
        tx_end = as.numeric(hit_rev + nchar(primer_rc) - 1)
      )
    )
  }
  
  hits
}

make_primer_positions_from_sequences <- function(primer_seq_csv_path,
                                                 fasta_dir,
                                                 fasta_suffix,
                                                 pb_gff_path,
                                                 out_csv_path) {
  if (!file.exists(primer_seq_csv_path)) {
    stop("Primer sequence CSV not found: ", primer_seq_csv_path)
  }
  
  primer_df <- read.csv(primer_seq_csv_path, stringsAsFactors = FALSE) %>%
    as_tibble()
  
  required_cols <- c("pbid", "primer_name", "sequence", "direction")
  missing_cols <- setdiff(required_cols, colnames(primer_df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Primer sequence CSV is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  primer_df <- primer_df %>%
    mutate(
      pbid = as.character(pbid),
      primer_name = as.character(primer_name),
      sequence = toupper(gsub("\\s+", "", as.character(sequence))),
      direction = toupper(as.character(direction)),
      direction = case_when(
        direction %in% c("F", "FORWARD", "+") ~ "F",
        direction %in% c("R", "REVERSE", "-") ~ "R",
        TRUE ~ direction
      )
    )
  
  res <- lapply(seq_len(nrow(primer_df)), function(i) {
    row <- primer_df[i, ]
    
    pbid <- row$pbid
    primer_name <- row$primer_name
    primer_seq <- row$sequence
    
    fasta_path <- file.path(fasta_dir, paste0(pbid, fasta_suffix))
    
    pbid_seq <- read_pbid_sequence_from_alignment_fasta(
      fasta_path = fasta_path,
      pbid = pbid
    )
    
    pbid_exons <- get_pbid_exons(
      pb_gff_path = pb_gff_path,
      pbid = pbid
    )
    
    pbid_exons_tx <- add_tx_coordinates(pbid_exons)
    
    hits <- find_primer_in_pbid_sequence(
      pbid_seq = pbid_seq,
      primer_seq = primer_seq
    )
    
    if (nrow(hits) == 0) {
      warning("No primer hit found: ", primer_name, " in ", pbid)
      
      return(tibble(
        pbid = pbid,
        primer_name = primer_name,
        sequence = primer_seq,
        direction = row$direction,
        plot_direction = NA_character_,
        match_orientation = NA_character_,
        tx_start = NA_real_,
        tx_end = NA_real_,
        seqnames = NA_character_,
        start = NA_real_,
        end = NA_real_,
        strand = NA_character_,
        n_hits = 0L
      ))
    }
    
    n_hits <- nrow(hits)
    
    hits %>%
      rowwise() %>%
      mutate(
        genomic_range = list(
          tx_interval_to_one_genomic_range(
            tx_start = tx_start,
            tx_end = tx_end,
            exons_tx = pbid_exons_tx
          )
        )
      ) %>%
      ungroup() %>%
      tidyr::unnest(genomic_range) %>%
      mutate(
        pbid = pbid,
        primer_name = primer_name,
        sequence = primer_seq,
        direction = row$direction,
        n_hits = n_hits,
        plot_direction = case_when(
          strand == "-" & direction == "F" ~ "R",
          strand == "-" & direction == "R" ~ "F",
          TRUE ~ direction
        )
      ) %>%
      select(
        pbid,
        primer_name,
        sequence,
        direction,
        plot_direction,
        match_orientation,
        tx_start,
        tx_end,
        seqnames,
        start,
        end,
        strand,
        n_hits
      )
  })
  
  primer_positions <- bind_rows(res)
  
  write.csv(
    primer_positions,
    out_csv_path,
    row.names = FALSE
  )
  
  primer_positions
}

#alignment summary
summarise_alignment_against_pbid <- function(aln_df, pbid) {
  ref_idx <- which(str_detect(aln_df$id, fixed(pbid)))
  
  if (length(ref_idx) != 1) {
    stop(
      "Expected exactly one reference sequence containing ",
      pbid,
      ", but found ",
      length(ref_idx),
      ".\nPlease make sure the aligned FASTA contains one PBid reference sequence."
    )
  }
  
  ref_seq <- aln_df$seq[ref_idx]
  read_df <- aln_df[-ref_idx, , drop = FALSE]
  
  if (nrow(read_df) == 0) {
    stop("No Sanger read sequence found in alignment for ", pbid)
  }
  
  ref_chars <- strsplit(ref_seq, "")[[1]]
  aln_len <- length(ref_chars)
  
  ref_non_gap <- ref_chars != "-"
  tx_pos_vec <- cumsum(ref_non_gap)
  tx_pos_vec[!ref_non_gap] <- NA_integer_
  
  read_chars_list <- lapply(read_df$seq, function(x) {
    strsplit(x, "")[[1]]
  })
  names(read_chars_list) <- read_df$id
  
  read_ranges <- lapply(read_chars_list, function(chars) {
    non_gap <- which(chars != "-")
    
    if (length(non_gap) == 0) {
      return(c(NA_integer_, NA_integer_))
    }
    
    c(min(non_gap), max(non_gap))
  })
  
  pos_summary <- lapply(seq_len(aln_len), function(col_i) {
    ref_base <- ref_chars[col_i]
    tx_pos <- tx_pos_vec[col_i]
    
    if (is.na(tx_pos)) {
      return(NULL)
    }
    
    active_bases <- c()
    deletion_n <- 0L
    covered_read_ids <- c()
    match_read_ids <- c()
    mismatch_read_ids <- c()
    deletion_read_ids <- c()
    
    for (read_id in names(read_chars_list)) {
      chars <- read_chars_list[[read_id]]
      rg <- read_ranges[[read_id]]
      
      if (is.na(rg[1]) || col_i < rg[1] || col_i > rg[2]) {
        next
      }
      
      b <- chars[col_i]
      
      if (b == "-") {
        deletion_n <- deletion_n + 1L
        covered_read_ids <- c(covered_read_ids, read_id)
        deletion_read_ids <- c(deletion_read_ids, read_id)
      } else if (b %in% c("A", "C", "G", "T", "N")) {
        active_bases <- c(active_bases, b)
        covered_read_ids <- c(covered_read_ids, read_id)
        
        if (
          b %in% c("A", "C", "G", "T") &&
          ref_base %in% c("A", "C", "G", "T") &&
          b == ref_base
        ) {
          match_read_ids <- c(match_read_ids, read_id)
        }
        
        if (
          b %in% c("A", "C", "G", "T") &&
          ref_base %in% c("A", "C", "G", "T") &&
          b != ref_base
        ) {
          mismatch_read_ids <- c(mismatch_read_ids, read_id)
        }
      }
    }
    
    support_n <- length(active_bases)
    match_n <- length(match_read_ids)
    mismatch_n <- length(mismatch_read_ids)
    conflict_n <- mismatch_n + deletion_n
    
    consensus_base <- "N"
    
    if (support_n > 0) {
      base_tbl <- sort(table(active_bases), decreasing = TRUE)
      consensus_base <- names(base_tbl)[1]
    }
    
    status <- case_when(
      match_n > 0 ~ "match",
      match_n == 0 & conflict_n > 0 ~ "conflict",
      TRUE ~ "no_support"
    )
    
    tibble(
      tx_pos = tx_pos,
      ref_base = ref_base,
      support_n = support_n,
      match_n = match_n,
      mismatch_n = mismatch_n,
      deletion_n = deletion_n,
      conflict_n = conflict_n,
      consensus_base = consensus_base,
      status = status,
      covered_read_ids = paste(unique(covered_read_ids), collapse = ";"),
      match_read_ids = paste(unique(match_read_ids), collapse = ";"),
      mismatch_read_ids = paste(unique(mismatch_read_ids), collapse = ";"),
      deletion_read_ids = paste(unique(deletion_read_ids), collapse = ";")
    )
  })
  
  bind_rows(pos_summary) %>%
    arrange(tx_pos)
}

make_status_runs <- function(pos_summary) {
  x <- pos_summary %>%
    filter(status %in% c("match", "conflict")) %>%
    arrange(tx_pos)
  
  if (nrow(x) == 0) {
    return(tibble(
      run_id = integer(),
      status = character(),
      tx_start = numeric(),
      tx_end = numeric(),
      support_n_max = numeric(),
      conflict_n_max = numeric()
    ))
  }
  
  x %>%
    mutate(
      prev_status = lag(status),
      prev_tx_pos = lag(tx_pos),
      is_new_run = case_when(
        row_number() == 1 ~ TRUE,
        status != prev_status ~ TRUE,
        tx_pos != prev_tx_pos + 1 ~ TRUE,
        TRUE ~ FALSE
      ),
      run_id = cumsum(is_new_run)
    ) %>%
    group_by(run_id, status) %>%
    summarise(
      tx_start = min(tx_pos),
      tx_end = max(tx_pos),
      support_n_max = max(support_n, na.rm = TRUE),
      conflict_n_max = max(conflict_n, na.rm = TRUE),
      .groups = "drop"
    )
}

#consensus FASTA
make_consensus_sequence <- function(pos_summary) {
  paste(pos_summary$consensus_base, collapse = "")
}

write_consensus_fasta <- function(pos_summary, pbid, out_path) {
  consensus_seq <- make_consensus_sequence(pos_summary)
  
  starts <- seq(1, nchar(consensus_seq), by = 80)
  ends <- pmin(starts + 79, nchar(consensus_seq))
  
  wrapped <- paste(
    substring(consensus_seq, starts, ends),
    collapse = "\n"
  )
  
  writeLines(
    c(
      paste0(">", pbid, "_merged_sanger_consensus"),
      wrapped
    ),
    out_path
  )
}

#final validation interval
endpoint_min_consecutive_match <- 5L

make_validation_interval <- function(pos_summary,
                                     min_consecutive_match = endpoint_min_consecutive_match) {
  
  x <- pos_summary %>%
    arrange(tx_pos) %>%
    mutate(is_match = status == "match")
  
  match_x <- x %>%
    filter(is_match) %>%
    arrange(tx_pos)
  
  if (nrow(match_x) == 0) {
    return(tibble(
      tx_start = numeric(),
      tx_end = numeric(),
      validated_n = integer(),
      total_n = integer(),
      endpoint_min_consecutive_match = integer()
    ))
  }
  
  match_runs <- match_x %>%
    mutate(
      prev_tx_pos = lag(tx_pos),
      is_new_run = case_when(
        row_number() == 1 ~ TRUE,
        tx_pos != prev_tx_pos + 1 ~ TRUE,
        TRUE ~ FALSE
      ),
      run_id = cumsum(is_new_run)
    ) %>%
    group_by(run_id) %>%
    summarise(
      tx_start = min(tx_pos),
      tx_end = max(tx_pos),
      run_len = n(),
      .groups = "drop"
    )
  
  endpoint_runs <- match_runs %>%
    filter(run_len >= min_consecutive_match)
  
  if (nrow(endpoint_runs) == 0) {
    warning(
      "No continuous match run >= ",
      min_consecutive_match,
      " bp. Falling back to all matched positions."
    )
    
    return(tibble(
      tx_start = min(match_x$tx_pos, na.rm = TRUE),
      tx_end = max(match_x$tx_pos, na.rm = TRUE),
      validated_n = nrow(match_x),
      total_n = nrow(pos_summary),
      endpoint_min_consecutive_match = min_consecutive_match
    ))
  }
  
  tibble(
    tx_start = min(endpoint_runs$tx_start, na.rm = TRUE),
    tx_end = max(endpoint_runs$tx_end, na.rm = TRUE),
    validated_n = nrow(match_x),
    total_n = nrow(pos_summary),
    endpoint_min_consecutive_match = min_consecutive_match
  )
}

#final plot
plot_merged_mapping <- function(pbid,
                                ref_exons_tx,
                                pbid_exons_tx,
                                mapping_blocks,
                                pos_summary,
                                out_prefix,
                                primer_positions_all = NULL) {
  
  ref_exons_tx <- ref_exons_tx %>%
    arrange(exon_order) %>%
    mutate(
      ref_exon_label = paste0("exon", exon_order)
    )
  
  ref_introns_plot <- make_introns_from_exons(ref_exons_tx)
  pbid_introns_plot <- make_introns_from_exons(pbid_exons_tx)
  
  validation_interval <- make_validation_interval(pos_summary)
  
  if (nrow(validation_interval) == 0) {
    warning("No matched Sanger-supported region found for ", pbid)
    validation_blocks <- tibble()
  } else {
    validation_blocks <- tx_interval_to_genomic_blocks(
      tx_start = validation_interval$tx_start[1],
      tx_end = validation_interval$tx_end[1],
      exons_tx = pbid_exons_tx
    ) %>%
      mutate(
        status = "validated",
        validated_n = validation_interval$validated_n[1],
        total_n = validation_interval$total_n[1]
      )
  }
  
  primer_plot_raw <- tibble()
  
  if (
    isTRUE(show_primers) &&
    !is.null(primer_positions_all) &&
    is.data.frame(primer_positions_all) &&
    nrow(primer_positions_all) > 0
  ) {
    primer_plot_raw <- primer_positions_all %>%
      filter(pbid == !!pbid, n_hits > 0)
  }
  
  if (isTRUE(zoom_to_mapping_region)) {
    coord_pool <- c(
      validation_blocks$start,
      validation_blocks$end,
      pbid_exons_tx$start,
      pbid_exons_tx$end,
      primer_plot_raw$start,
      primer_plot_raw$end
    )
  } else {
    coord_pool <- c(
      validation_blocks$start,
      validation_blocks$end,
      ref_exons_tx$start,
      ref_exons_tx$end,
      pbid_exons_tx$start,
      pbid_exons_tx$end,
      primer_plot_raw$start,
      primer_plot_raw$end
    )
  }
  
  coord_pool <- coord_pool[is.finite(coord_pool)]
  
  plot_start_raw <- min(coord_pool, na.rm = TRUE) - plot_pad_bp
  plot_end_raw <- max(coord_pool, na.rm = TRUE) + plot_pad_bp
  
  ref_exons_plot <- ref_exons_tx %>%
    filter(end >= plot_start_raw, start <= plot_end_raw)
  
  ref_introns_plot <- ref_introns_plot %>%
    filter(end >= plot_start_raw, start <= plot_end_raw)
  
  pbid_exons_plot <- pbid_exons_tx %>%
    filter(end >= plot_start_raw, start <= plot_end_raw)
  
  pbid_introns_plot <- pbid_introns_plot %>%
    filter(end >= plot_start_raw, start <= plot_end_raw)
  
  validation_blocks_plot <- validation_blocks %>%
    filter(end >= plot_start_raw, start <= plot_end_raw)
  
  primer_plot <- primer_plot_raw %>%
    filter(end >= plot_start_raw, start <= plot_end_raw)
  
  if (isTRUE(use_compressed_x)) {
    axis_features <- bind_rows(
      ref_exons_plot %>% select(start, end),
      pbid_exons_plot %>% select(start, end),
      validation_blocks_plot %>% select(start, end),
      primer_plot %>% select(start, end)
    )
    
    axis_tbl <- build_compressed_axis(
      feature_df = axis_features,
      flank_bp = compress_feature_flank_bp,
      gap_threshold_bp = compress_gap_threshold_bp,
      compressed_gap_width = compressed_gap_width
    )
    
    ref_exons_plot <- add_plot_x(ref_exons_plot, axis_tbl)
    ref_introns_plot <- add_plot_x(ref_introns_plot, axis_tbl)
    pbid_exons_plot <- add_plot_x(pbid_exons_plot, axis_tbl)
    pbid_introns_plot <- add_plot_x(pbid_introns_plot, axis_tbl)
    validation_blocks_plot <- add_plot_x(validation_blocks_plot, axis_tbl)
    
    if (nrow(primer_plot) > 0) {
      primer_plot <- add_plot_x(primer_plot, axis_tbl)
    }
    
    x_min <- min(axis_tbl$plot_start, na.rm = TRUE)
    x_max <- max(axis_tbl$plot_end, na.rm = TRUE)
    
    x_breaks <- axis_tbl$plot_start
    x_labels <- format(round(axis_tbl$orig_start), scientific = FALSE)
    
    x_lab <- "Compressed genomic position"
  } else {
    ref_exons_plot <- ref_exons_plot %>%
      mutate(plot_xmin = start, plot_xmax = end)
    
    ref_introns_plot <- ref_introns_plot %>%
      mutate(plot_xmin = start, plot_xmax = end)
    
    pbid_exons_plot <- pbid_exons_plot %>%
      mutate(plot_xmin = start, plot_xmax = end)
    
    pbid_introns_plot <- pbid_introns_plot %>%
      mutate(plot_xmin = start, plot_xmax = end)
    
    validation_blocks_plot <- validation_blocks_plot %>%
      mutate(plot_xmin = start, plot_xmax = end)
    
    if (nrow(primer_plot) > 0) {
      primer_plot <- primer_plot %>%
        mutate(plot_xmin = start, plot_xmax = end)
    }
    
    x_min <- plot_start_raw
    x_max <- plot_end_raw
    
    x_breaks <- waiver()
    x_labels <- waiver()
    
    x_lab <- "Genomic position"
  }
  
  ref_introns_plot <- ref_introns_plot %>%
    mutate(
      plot_xmin = pmax(plot_xmin, x_min),
      plot_xmax = pmin(plot_xmax, x_max)
    ) %>%
    filter(plot_xmax > plot_xmin)
  
  pbid_introns_plot <- pbid_introns_plot %>%
    mutate(
      plot_xmin = pmax(plot_xmin, x_min),
      plot_xmax = pmin(plot_xmax, x_max)
    ) %>%
    filter(plot_xmax > plot_xmin)
  
  if (nrow(primer_plot) > 0) {
    primer_plot <- primer_plot %>%
      mutate(
        final_direction = ifelse(
          !is.na(plot_direction),
          plot_direction,
          direction
        ),
        arrow_x = ifelse(final_direction == "F", plot_xmin, plot_xmax),
        arrow_xend = ifelse(final_direction == "F", plot_xmax, plot_xmin),
        label_x = (plot_xmin + plot_xmax) / 2
      )
  }
  
  if (nrow(validation_blocks_plot) > 0) {
    validation_box_plot <- tibble(
      plot_xmin = min(validation_blocks_plot$plot_xmin, na.rm = TRUE),
      plot_xmax = max(validation_blocks_plot$plot_xmax, na.rm = TRUE)
    )
  } else {
    validation_box_plot <- tibble()
  }
  
  y_ref <- 2
  y_pbid <- 1
  
  #primer position
  primer_y <- y_pbid + 0.35
  primer_label_y <- y_pbid + 0.43
  primer_label_size_use <- 5
  
  validation_ymin <- y_pbid - pbid_exon_half_height - 0.10
  validation_ymax <- y_pbid + pbid_exon_half_height + 0.10
  
  p <- ggplot() +
    geom_segment(
      data = ref_introns_plot,
      aes(x = plot_xmin, xend = plot_xmax, y = y_ref, yend = y_ref),
      linewidth = 0.45,
      color = "grey55"
    ) +
    geom_rect(
      data = ref_exons_plot,
      aes(
        xmin = plot_xmin,
        xmax = plot_xmax,
        ymin = y_ref - reference_exon_half_height,
        ymax = y_ref + reference_exon_half_height
      ),
      fill = "grey25",
      color = "grey25"
    ) +
    
    #reference exon labels
    geom_text(
      data = ref_exons_plot,
      aes(
        x = (plot_xmin + plot_xmax) / 2,
        y = y_ref + reference_exon_half_height + 0.18,
        label = ref_exon_label
      ),
      size = 5,
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      color = "grey20"
    ) +
    
    geom_segment(
      data = pbid_introns_plot,
      aes(x = plot_xmin, xend = plot_xmax, y = y_pbid, yend = y_pbid),
      linewidth = 0.45,
      color = "grey55"
    ) +
    geom_rect(
      data = pbid_exons_plot,
      aes(
        xmin = plot_xmin,
        xmax = plot_xmax,
        ymin = y_pbid - pbid_exon_half_height,
        ymax = y_pbid + pbid_exon_half_height
      ),
      fill = "grey55",
      color = "grey25"
    )
  
  if (nrow(validation_box_plot) > 0) {
    p <- p +
      geom_rect(
        data = validation_box_plot,
        aes(
          xmin = plot_xmin,
          xmax = plot_xmax,
          ymin = validation_ymin,
          ymax = validation_ymax
        ),
        fill = NA,
        color = "#2ca25f",
        linewidth = 1.15
      )
  }
  
  if (nrow(primer_plot) > 0) {
    p <- p +
      geom_segment(
        data = primer_plot,
        aes(
          x = arrow_x,
          xend = arrow_xend,
          y = primer_y,
          yend = primer_y
        ),
        linewidth = 0.9,
        color = "orange3",
        arrow = arrow(
          length = unit(primer_arrow_size, "cm"),
          type = "closed"
        )
      ) +
      geom_text(
        data = primer_plot,
        aes(
          x = label_x,
          y = primer_label_y,
          label = primer_name
        ),
        size = primer_label_size_use,
        color = "orange4",
        vjust = 0,
        hjust = 0.5
      )
  }
  
  p <- p +
    scale_y_continuous(
      breaks = c(y_ref, y_pbid),
      labels = c(
        paste0("Reference | ", reference_transcript_name),
        paste0("PBid | ", pbid)
      ),
      limits = c(0.35, 2.85)
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels
    ) +
    coord_cartesian(
      xlim = c(x_min, x_max),
      clip = "off"
    ) +
    labs(
      title = paste0("Sanger-validated region: ", pbid),
      x = x_lab,
      y = NULL
    ) +
    theme_bw(base_size = 16) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 14, margin = margin(r = 12)),
      axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
      axis.title.x = element_text(size = 16, margin = margin(t = 10)),
      plot.title = element_text(face = "bold", size = 18),
      legend.position = "none",
      plot.margin = margin(45, 20, 20, 35)
    )
  
  ggsave(
    paste0(out_prefix, "_exon_labeled_connected_validation.pdf"),
    p,
    width = fig_width,
    height = fig_height,
    units = "in",
    limitsize = FALSE
  )
  
  ggsave(
    paste0(out_prefix, "_exon_labeled_connected_validation.png"),
    p,
    width = fig_width,
    height = fig_height,
    units = "in",
    dpi = fig_dpi,
    limitsize = FALSE
  )
  
  p
}

#main function for one PBid
run_one_pbid <- function(pbid) {
  message("Processing ", pbid)
  
  fasta_path <- file.path(fasta_dir, paste0(pbid, fasta_suffix))
  
  aln_df <- read_fasta_as_strings(fasta_path)
  
  pbid_exons <- get_pbid_exons(
    pb_gff_path = pb_gff_path,
    pbid = pbid
  )
  
  pbid_exons_tx <- add_tx_coordinates(pbid_exons)
  
  ref_exons <- get_reference_exons(
    ref_gff_path = ref_gff_path,
    gene_name_target = reference_gene_name,
    transcript_name_target = reference_transcript_name
  )
  
  ref_exons_tx <- add_tx_coordinates(ref_exons)
  
  pos_summary <- summarise_alignment_against_pbid(
    aln_df = aln_df,
    pbid = pbid
  )
  
  status_runs <- make_status_runs(pos_summary)
  
  if (nrow(status_runs) == 0) {
    stop("No matched or conflicted covered regions found for ", pbid)
  }
  
  mapping_blocks <- pmap_dfr(
    status_runs,
    function(run_id,
             status,
             tx_start,
             tx_end,
             support_n_max,
             conflict_n_max) {
      tx_interval_to_genomic_blocks(
        tx_start = tx_start,
        tx_end = tx_end,
        exons_tx = pbid_exons_tx
      ) %>%
        mutate(
          status = status,
          support_n_max = support_n_max,
          conflict_n_max = conflict_n_max
        )
    }
  )
  
  out_prefix <- file.path(out_dir, paste0(pbid, "_merged_mapping"))
  
  write.csv(
    pos_summary,
    paste0(out_prefix, "_position_summary.csv"),
    row.names = FALSE
  )
  
  write.csv(
    status_runs,
    paste0(out_prefix, "_status_runs.csv"),
    row.names = FALSE
  )
  
  write.csv(
    mapping_blocks,
    paste0(out_prefix, "_mapping_blocks.csv"),
    row.names = FALSE
  )
  
  write.csv(
    pbid_exons_tx,
    paste0(out_prefix, "_pbid_exons_with_tx_coordinates.csv"),
    row.names = FALSE
  )
  
  write.csv(
    ref_exons_tx,
    paste0(out_prefix, "_reference_exons_with_tx_coordinates.csv"),
    row.names = FALSE
  )
  
  write_consensus_fasta(
    pos_summary = pos_summary,
    pbid = pbid,
    out_path = paste0(out_prefix, "_consensus.fasta")
  )
  
  plot_merged_mapping(
    pbid = pbid,
    ref_exons_tx = ref_exons_tx,
    pbid_exons_tx = pbid_exons_tx,
    mapping_blocks = mapping_blocks,
    pos_summary = pos_summary,
    out_prefix = out_prefix,
    primer_positions_all = primer_positions_all
  )
}

#primer positions
primer_positions_csv_out <- file.path(
  fasta_dir,
  "primer_positions.csv"
)

primer_positions_all <- make_primer_positions_from_sequences(
  primer_seq_csv_path = primer_seq_csv_path,
  fasta_dir = fasta_dir,
  fasta_suffix = fasta_suffix,
  pb_gff_path = pb_gff_path,
  out_csv_path = primer_positions_csv_out
)

write.csv(
  primer_positions_all,
  file.path(out_dir, "primer_positions_used_for_plot.csv"),
  row.names = FALSE
)

#run selected PBids
#supplementary figure 9B-C
plots <- lapply(target_pbids, run_one_pbid)



#generate plots for Alphafold3
#plot PBid transcript with CDS red box, reference exon labels, and protein residue axis
target_pbids <- c(
  "PB.5877.897",
  "PB.5877.945"
)

orfinder_cds_files <- tibble(
  pbid = c(
    "PB.5877.897",
    "PB.5877.945"
  ),
  cds_path = c(
    file.path(fasta_dir, "PB.5877.897.ORF8.cds"),
    file.path(fasta_dir, "PB.5877.945.ORF6.cds")
  )
)

cds_box_color <- "#de2d26"
cds_box_linewidth <- 1.15

#Protein residue labels shown in AlphaFold3 PAE heatmap
protein_residue_breaks_by_pbid <- list(
  "PB.5877.945" = c(1, 33, 66, 99, 132, 165, 198),
  "PB.5877.897" = c(1, 21, 42, 63, 84, 105, 126, 147, 168)
)

#clean DNA sequence
clean_dna_sequence <- function(x) {
  x <- toupper(as.character(x))
  x <- gsub("\\s+", "", x)
  x <- gsub("-", "", x)
  x <- gsub("[^ACGTN]", "", x)
  x
}

#read ORFinder .cds FASTA
read_orfinder_cds_fasta <- function(path, pbid) {
  if (!file.exists(path)) {
    stop("ORFinder CDS file not found: ", path)
  }
  
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[lines != ""]
  
  header <- lines[startsWith(lines, ">")]
  header <- ifelse(length(header) > 0, sub("^>", "", header[1]), NA_character_)
  
  seq_lines <- lines[!startsWith(lines, ">")]
  cds_seq <- paste(seq_lines, collapse = "")
  cds_seq <- clean_dna_sequence(cds_seq)
  
  if (nchar(cds_seq) == 0) {
    stop("No CDS sequence found in: ", path)
  }
  
  tibble(
    pbid = pbid,
    cds_path = path,
    cds_header = header,
    cds_sequence = cds_seq,
    cds_length = nchar(cds_seq)
  )
}

#read all ORFinder CDS files
cds_sequence_df <- purrr::pmap_dfr(
  orfinder_cds_files,
  function(pbid, cds_path) {
    read_orfinder_cds_fasta(
      path = cds_path,
      pbid = pbid
    )
  }
)

write.csv(
  cds_sequence_df,
  file.path(out_dir, "orfinder_cds_sequences_used_for_plot.csv"),
  row.names = FALSE
)

print(cds_sequence_df)

#find CDS sequence inside PBid transcript sequence
find_cds_in_pbid_sequence <- function(pbid_seq, cds_seq) {
  pbid_seq <- clean_dna_sequence(pbid_seq)
  cds_seq <- clean_dna_sequence(cds_seq)
  cds_rc <- reverse_complement_dna(cds_seq)
  
  hit_fwd <- gregexpr(cds_seq, pbid_seq, fixed = TRUE)[[1]]
  hit_rev <- gregexpr(cds_rc, pbid_seq, fixed = TRUE)[[1]]
  
  hit_fwd <- hit_fwd[hit_fwd > 0]
  hit_rev <- hit_rev[hit_rev > 0]
  
  hits <- tibble()
  
  if (length(hit_fwd) > 0) {
    hits <- bind_rows(
      hits,
      tibble(
        match_orientation = "cds_sequence",
        tx_start = as.numeric(hit_fwd),
        tx_end = as.numeric(hit_fwd + nchar(cds_seq) - 1)
      )
    )
  }
  
  if (length(hit_rev) > 0) {
    hits <- bind_rows(
      hits,
      tibble(
        match_orientation = "reverse_complement",
        tx_start = as.numeric(hit_rev),
        tx_end = as.numeric(hit_rev + nchar(cds_rc) - 1)
      )
    )
  }
  
  hits
}

#convert CDS transcript coordinates into genomic/exon blocks
make_cds_positions_from_sequences <- function(cds_sequence_df,
                                              fasta_dir,
                                              fasta_suffix,
                                              pb_gff_path,
                                              target_pbids = NULL,
                                              out_csv_path = NULL) {
  
  required_cols <- c("pbid", "cds_sequence")
  missing_cols <- setdiff(required_cols, colnames(cds_sequence_df))
  
  if (length(missing_cols) > 0) {
    stop(
      "cds_sequence_df is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  cds_df <- cds_sequence_df %>%
    mutate(
      pbid = as.character(pbid),
      cds_sequence = clean_dna_sequence(cds_sequence)
    )
  
  if (!is.null(target_pbids)) {
    cds_df <- cds_df %>%
      filter(pbid %in% target_pbids)
  }
  
  if (nrow(cds_df) == 0) {
    stop("No CDS records left after filtering target_pbids.")
  }
  
  res <- lapply(seq_len(nrow(cds_df)), function(i) {
    row <- cds_df[i, ]
    
    pbid <- row$pbid
    cds_seq <- row$cds_sequence
    
    fasta_path <- file.path(fasta_dir, paste0(pbid, fasta_suffix))
    
    pbid_seq <- read_pbid_sequence_from_alignment_fasta(
      fasta_path = fasta_path,
      pbid = pbid
    )
    
    hits <- find_cds_in_pbid_sequence(
      pbid_seq = pbid_seq,
      cds_seq = cds_seq
    )
    
    if (nrow(hits) == 0) {
      warning("No CDS hit found in PBid transcript sequence: ", pbid)
      
      return(tibble(
        pbid = pbid,
        cds_length = nchar(cds_seq),
        match_orientation = NA_character_,
        tx_start_cds = NA_real_,
        tx_end_cds = NA_real_,
        seqnames = NA_character_,
        start = NA_real_,
        end = NA_real_,
        strand = NA_character_,
        tx_start = NA_real_,
        tx_end = NA_real_,
        exon_order = NA_integer_,
        n_hits = 0L
      ))
    }
    
    if (nrow(hits) > 1) {
      warning(
        "Multiple CDS hits found for ",
        pbid,
        ". The first hit will be used for plotting."
      )
    }
    
    hit <- hits[1, ]
    
    pbid_exons <- get_pbid_exons(
      pb_gff_path = pb_gff_path,
      pbid = pbid
    )
    
    pbid_exons_tx <- add_tx_coordinates(pbid_exons)
    
    cds_blocks <- tx_interval_to_genomic_blocks(
      tx_start = hit$tx_start,
      tx_end = hit$tx_end,
      exons_tx = pbid_exons_tx
    ) %>%
      mutate(
        pbid = pbid,
        cds_length = nchar(cds_seq),
        match_orientation = hit$match_orientation,
        tx_start_cds = hit$tx_start,
        tx_end_cds = hit$tx_end,
        n_hits = nrow(hits)
      ) %>%
      select(
        pbid,
        cds_length,
        match_orientation,
        tx_start_cds,
        tx_end_cds,
        seqnames,
        start,
        end,
        strand,
        tx_start,
        tx_end,
        exon_order,
        n_hits
      )
    
    cds_blocks
  })
  
  cds_positions <- bind_rows(res)
  
  if (!is.null(out_csv_path)) {
    write.csv(cds_positions, out_csv_path, row.names = FALSE)
  }
  
  cds_positions
}

#convert protein residue numbers into transcript CDS positions
make_protein_residue_positions <- function(pbid,
                                           cds_positions_all,
                                           pbid_exons_tx,
                                           residue_breaks_by_pbid) {
  
  cds_blocks <- cds_positions_all %>%
    filter(pbid == !!pbid, n_hits > 0)
  
  if (nrow(cds_blocks) == 0) {
    return(tibble())
  }
  
  if (!pbid %in% names(residue_breaks_by_pbid)) {
    warning("No protein residue breaks provided for ", pbid)
    return(tibble())
  }
  
  residue_breaks <- residue_breaks_by_pbid[[pbid]]
  
  cds_tx_start <- min(cds_blocks$tx_start_cds, na.rm = TRUE)
  cds_tx_end <- max(cds_blocks$tx_end_cds, na.rm = TRUE)
  cds_len <- cds_tx_end - cds_tx_start + 1
  
  max_residue_possible <- floor(cds_len / 3)
  
  residue_breaks <- residue_breaks[
    residue_breaks >= 1 &
      residue_breaks <= max_residue_possible
  ]
  
  if (length(residue_breaks) == 0) {
    warning("No valid residue breaks left for ", pbid)
    return(tibble())
  }
  
  residue_df <- tibble(
    pbid = pbid,
    residue = residue_breaks,
    #Use the first nucleotide of each codon for plotting
    tx_pos = cds_tx_start + (residue_breaks - 1) * 3 + 1
  )
  
  residue_pos <- lapply(seq_len(nrow(residue_df)), function(i) {
    row <- residue_df[i, ]
    
    g <- tx_interval_to_one_genomic_range(
      tx_start = row$tx_pos,
      tx_end = row$tx_pos,
      exons_tx = pbid_exons_tx
    )
    
    bind_cols(row, g)
  }) %>%
    bind_rows()
  
  residue_pos
}

#plot PBid with CDS red box and ABCA4-201 exon labels
plot_pbid_with_cds_box <- function(pbid,
                                   ref_exons_tx,
                                   pbid_exons_tx,
                                   cds_positions_all,
                                   out_prefix,
                                   primer_positions_all = NULL,
                                   show_primers_on_cds_plot = FALSE) {
  
  ref_exons_tx <- ref_exons_tx %>%
    arrange(exon_order) %>%
    mutate(
      ref_exon_label = paste0("exon", exon_order)
    )
  
  ref_introns_plot <- make_introns_from_exons(ref_exons_tx)
  pbid_introns_plot <- make_introns_from_exons(pbid_exons_tx)
  
  cds_blocks <- cds_positions_all %>%
    filter(pbid == !!pbid, n_hits > 0)
  
  if (nrow(cds_blocks) == 0) {
    warning("No CDS blocks available for plotting: ", pbid)
  }
  
  #Protein residue positions for matching AlphaFold3 PAE heatmap axis labels
  residue_plot_raw <- make_protein_residue_positions(
    pbid = pbid,
    cds_positions_all = cds_positions_all,
    pbid_exons_tx = pbid_exons_tx,
    residue_breaks_by_pbid = protein_residue_breaks_by_pbid
  )
  
  #zoom plot of the ABCA4-201 reference to the PBid/CDS region
  #otherwise the reference track will force the plot to show the whole transcript.
  coord_pool <- c(
    cds_blocks$start,
    cds_blocks$end,
    pbid_exons_tx$start,
    pbid_exons_tx$end,
    residue_plot_raw$start,
    residue_plot_raw$end
  )
  
  coord_pool <- coord_pool[is.finite(coord_pool)]
  
  if (length(coord_pool) == 0) {
    stop("No finite coordinates available for plotting: ", pbid)
  }
  
  plot_start_raw <- min(coord_pool, na.rm = TRUE) - plot_pad_bp
  plot_end_raw <- max(coord_pool, na.rm = TRUE) + plot_pad_bp
  
  #only keep the ABCA4-201 reference exons overlapping the PBid/CDS region
  ref_exons_plot <- ref_exons_tx %>%
    filter(end >= plot_start_raw, start <= plot_end_raw)
  
  ref_introns_plot <- ref_introns_plot %>%
    filter(end >= plot_start_raw, start <= plot_end_raw)
  
  pbid_exons_plot <- pbid_exons_tx %>%
    filter(end >= plot_start_raw, start <= plot_end_raw)
  
  pbid_introns_plot <- pbid_introns_plot %>%
    filter(end >= plot_start_raw, start <= plot_end_raw)
  
  cds_blocks_plot <- cds_blocks %>%
    filter(end >= plot_start_raw, start <= plot_end_raw)
  
  residue_plot <- residue_plot_raw %>%
    filter(end >= plot_start_raw, start <= plot_end_raw)
  
  if (isTRUE(use_compressed_x)) {
    axis_features <- bind_rows(
      ref_exons_plot %>% select(start, end),
      pbid_exons_plot %>% select(start, end),
      cds_blocks_plot %>% select(start, end),
      residue_plot %>% select(start, end)
    )
    
    axis_tbl <- build_compressed_axis(
      feature_df = axis_features,
      flank_bp = compress_feature_flank_bp,
      gap_threshold_bp = compress_gap_threshold_bp,
      compressed_gap_width = compressed_gap_width
    )
    
    ref_exons_plot <- add_plot_x(ref_exons_plot, axis_tbl)
    ref_introns_plot <- add_plot_x(ref_introns_plot, axis_tbl)
    pbid_exons_plot <- add_plot_x(pbid_exons_plot, axis_tbl)
    pbid_introns_plot <- add_plot_x(pbid_introns_plot, axis_tbl)
    cds_blocks_plot <- add_plot_x(cds_blocks_plot, axis_tbl)
    
    if (nrow(residue_plot) > 0) {
      residue_plot <- add_plot_x(residue_plot, axis_tbl) %>%
        mutate(
          plot_x = (plot_xmin + plot_xmax) / 2
        )
    }
    
    x_min <- min(axis_tbl$plot_start, na.rm = TRUE)
    x_max <- max(axis_tbl$plot_end, na.rm = TRUE)
    
    x_breaks <- axis_tbl$plot_start
    x_labels <- format(round(axis_tbl$orig_start), scientific = FALSE)
    
    x_lab <- "Compressed genomic position"
  } else {
    ref_exons_plot <- ref_exons_plot %>%
      mutate(plot_xmin = start, plot_xmax = end)
    
    ref_introns_plot <- ref_introns_plot %>%
      mutate(plot_xmin = start, plot_xmax = end)
    
    pbid_exons_plot <- pbid_exons_plot %>%
      mutate(plot_xmin = start, plot_xmax = end)
    
    pbid_introns_plot <- pbid_introns_plot %>%
      mutate(plot_xmin = start, plot_xmax = end)
    
    cds_blocks_plot <- cds_blocks_plot %>%
      mutate(plot_xmin = start, plot_xmax = end)
    
    if (nrow(residue_plot) > 0) {
      residue_plot <- residue_plot %>%
        mutate(
          plot_xmin = start,
          plot_xmax = end,
          plot_x = (plot_xmin + plot_xmax) / 2
        )
    }
    
    x_min <- plot_start_raw
    x_max <- plot_end_raw
    
    x_breaks <- waiver()
    x_labels <- waiver()
    
    x_lab <- "Genomic position"
  }
  
  ref_introns_plot <- ref_introns_plot %>%
    mutate(
      plot_xmin = pmax(plot_xmin, x_min),
      plot_xmax = pmin(plot_xmax, x_max)
    ) %>%
    filter(plot_xmax > plot_xmin)
  
  pbid_introns_plot <- pbid_introns_plot %>%
    mutate(
      plot_xmin = pmax(plot_xmin, x_min),
      plot_xmax = pmin(plot_xmax, x_max)
    ) %>%
    filter(plot_xmax > plot_xmin)
  
  #draw CDS as one connected red box across the PBid row
  if (nrow(cds_blocks_plot) > 0) {
    cds_box_plot <- tibble(
      plot_xmin = min(cds_blocks_plot$plot_xmin, na.rm = TRUE),
      plot_xmax = max(cds_blocks_plot$plot_xmax, na.rm = TRUE)
    )
  } else {
    cds_box_plot <- tibble()
  }
  
  y_ref <- 2
  y_pbid <- 1
  
  cds_ymin <- y_pbid - pbid_exon_half_height - 0.10
  cds_ymax <- y_pbid + pbid_exon_half_height + 0.10
  
  #protein residue tick and label positions
  residue_tick_ymin <- cds_ymin - 0.06
  residue_tick_ymax <- cds_ymin - 0.18
  residue_label_y <- cds_ymin - 0.28
  
  p <- ggplot() +
    geom_segment(
      data = ref_introns_plot,
      aes(x = plot_xmin, xend = plot_xmax, y = y_ref, yend = y_ref),
      linewidth = 0.45,
      color = "grey55"
    ) +
    geom_rect(
      data = ref_exons_plot,
      aes(
        xmin = plot_xmin,
        xmax = plot_xmax,
        ymin = y_ref - reference_exon_half_height,
        ymax = y_ref + reference_exon_half_height
      ),
      fill = "grey25",
      color = "grey25"
    ) +
    
    #reference exon labels
    geom_text(
      data = ref_exons_plot,
      aes(
        x = (plot_xmin + plot_xmax) / 2,
        y = y_ref + reference_exon_half_height + 0.18,
        label = ref_exon_label
      ),
      size = 5,
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      color = "grey20"
    ) +
    geom_segment(
      data = pbid_introns_plot,
      aes(x = plot_xmin, xend = plot_xmax, y = y_pbid, yend = y_pbid),
      linewidth = 0.45,
      color = "grey55"
    ) +
    geom_rect(
      data = pbid_exons_plot,
      aes(
        xmin = plot_xmin,
        xmax = plot_xmax,
        ymin = y_pbid - pbid_exon_half_height,
        ymax = y_pbid + pbid_exon_half_height
      ),
      fill = "grey55",
      color = "grey25"
    )
  
  if (nrow(cds_box_plot) > 0) {
    p <- p +
      geom_rect(
        data = cds_box_plot,
        aes(
          xmin = plot_xmin,
          xmax = plot_xmax,
          ymin = cds_ymin,
          ymax = cds_ymax
        ),
        fill = NA,
        color = cds_box_color,
        linewidth = cds_box_linewidth
      )
  }
  
  #protein residue ticks and labels matching AlphaFold3 PAE heatmap
  if (nrow(residue_plot) > 0) {
    p <- p +
      geom_segment(
        data = residue_plot,
        aes(
          x = plot_x,
          xend = plot_x,
          y = residue_tick_ymin,
          yend = residue_tick_ymax
        ),
        linewidth = 0.5,
        color = "grey20"
      ) +
      geom_text(
        data = residue_plot,
        aes(
          x = plot_x,
          y = residue_label_y,
          label = residue
        ),
        size = 4.5,
        angle = 0,
        hjust = 0.5,
        vjust = 1,
        color = "grey20"
      )
  }
  
  p <- p +
    scale_y_continuous(
      breaks = c(y_ref, y_pbid),
      labels = c(
        paste0("Reference | ", reference_transcript_name),
        paste0("PBid | ", pbid)
      ),
      limits = c(0.15, 2.85)
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels
    ) +
    coord_cartesian(
      xlim = c(x_min, x_max),
      clip = "off"
    ) +
    labs(
      title = paste0("CDS region on PBid: ", pbid),
      x = x_lab,
      y = NULL
    ) +
    theme_bw(base_size = 16) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 14, margin = margin(r = 12)),
      axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
      axis.title.x = element_text(size = 16, margin = margin(t = 10)),
      plot.title = element_text(face = "bold", size = 18),
      legend.position = "none",
      plot.margin = margin(45, 20, 35, 35)
    )
  
  ggsave(
    paste0(out_prefix, "_cds_red_box_exon_labeled_residue_axis.pdf"),
    p,
    width = fig_width,
    height = fig_height,
    units = "in",
    limitsize = FALSE
  )
  
  ggsave(
    paste0(out_prefix, "_cds_red_box_exon_labeled_residue_axis.png"),
    p,
    width = fig_width,
    height = fig_height,
    units = "in",
    dpi = fig_dpi,
    limitsize = FALSE
  )
  
  p
}

#run one PBid
run_one_pbid_cds_plot <- function(pbid) {
  message("Plotting CDS region for ", pbid)
  
  pbid_exons <- get_pbid_exons(
    pb_gff_path = pb_gff_path,
    pbid = pbid
  )
  
  pbid_exons_tx <- add_tx_coordinates(pbid_exons)
  
  ref_exons <- get_reference_exons(
    ref_gff_path = ref_gff_path,
    gene_name_target = reference_gene_name,
    transcript_name_target = reference_transcript_name
  )
  
  ref_exons_tx <- add_tx_coordinates(ref_exons)
  
  out_prefix <- file.path(out_dir, paste0(pbid, "_cds_region"))
  
  plot_pbid_with_cds_box(
    pbid = pbid,
    ref_exons_tx = ref_exons_tx,
    pbid_exons_tx = pbid_exons_tx,
    cds_positions_all = cds_positions_all,
    out_prefix = out_prefix,
    primer_positions_all = NULL,
    show_primers_on_cds_plot = FALSE
  )
}

#generate CDS positions and plots
cds_positions_csv_out <- file.path(
  out_dir,
  "cds_positions_used_for_plot.csv"
)

cds_positions_all <- make_cds_positions_from_sequences(
  cds_sequence_df = cds_sequence_df,
  fasta_dir = fasta_dir,
  fasta_suffix = fasta_suffix,
  pb_gff_path = pb_gff_path,
  target_pbids = target_pbids,
  out_csv_path = cds_positions_csv_out
)

cds_plots <- lapply(target_pbids, run_one_pbid_cds_plot)

#FLNC 5'/3' end bias analysis
#path
vis_dir <- "visualization"
dir.create(vis_dir, showWarnings = FALSE, recursive = TRUE)
#end-offset analysis: read 5'->TSS and 3'->TES distances
#read bedtools closest outputs
#format: BED(A) + BED(B) + signed distance (last column)
t5 <- fread("read5p_vs_tss.closest.txt", header = FALSE)
t3 <- fread("read3p_vs_tes.closest.txt", header = FALSE)

#extract absolute distance from the last column
t5[, dist := abs(get(paste0("V", ncol(t5))))]
t3[, dist := abs(get(paste0("V", ncol(t3))))]

#summarize end offsets within fixed windows
summ_within <- function(dt, label){
  data.table(
    metric = label,
    n = nrow(dt),
    within_50bp  = mean(dt$dist <= 50),
    within_100bp = mean(dt$dist <= 100),
    within_200bp = mean(dt$dist <= 200),
    median_bp = median(dt$dist),
    p95_bp = as.numeric(quantile(dt$dist, 0.95))
  )
}

summary_table <- rbind(
  summ_within(t5, "5prime_to_TSS"),
  summ_within(t3, "3prime_to_TES")
)

#save summary statistics
fwrite(
  summary_table,
  file.path(vis_dir, "end_offset_summary.tsv"),
  sep = "\t"
)

#distance distribution (histogram, log10 scale)
p_dist <- rbind(
  data.table(type = "5' end → TSS", dist = t5$dist),
  data.table(type = "3' end → TES", dist = t3$dist)
)

g1 <- ggplot(p_dist, aes(x = dist)) +
  geom_histogram(bins = 200) +
  scale_x_continuous(trans = "log10") +
  facet_wrap(~type, ncol = 1, scales = "free_y") +
  labs(
    x = "Absolute distance to annotated end (bp, log10 scale)",
    y = "Read count",
    title = "PacBio FLNC read end offsets relative to annotated TSS/TES"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(vis_dir, "end_offset_hist_log10.png"),
  g1, width = 7, height = 6, dpi = 300
)

#cumulative distribution function (CDF)
g2 <- ggplot(p_dist, aes(x = dist, color = type)) +
  stat_ecdf(linewidth = 1) +
  scale_x_continuous(trans = "log10") +
  labs(
    x = "Absolute distance to annotated end (bp, log10 scale)",
    y = "Empirical cumulative distribution",
    title = "CDF of read end offsets relative to transcript boundaries"
  ) +
  theme_bw(base_size = 12)

#supplementary figure 2A
ggsave(
  file.path(vis_dir, "end_offset_cdf_log10.png"),
  g2, width = 7, height = 4.5, dpi = 300
)

#gene-body relative position density curve
#input: read, chr, read_start, read_end, strand, gene_id, gene_name, gene_start, gene_end
rb <- fread(
  "read_best_gene.tsv",
  col.names = c("read", "chr", "rs", "re", "strand", "gene", "gs", "ge")
)

#compute read midpoint
rb[, mid := floor((rs + re) / 2)]

#compute gene length and filter invalid genes
rb[, gene_len := ge - gs]
rb <- rb[gene_len > 0]

#compute relative position along gene body
#0 = TSS side, 1 = TES side (strand-aware)
rb[strand == "+", rel := (mid - gs) / gene_len]
rb[strand == "-", rel := (ge - mid) / gene_len]

#keep reads mapping within gene body
rb <- rb[rel >= 0 & rel <= 1]

#bin relative positions and compute density
nbin <- 100L
rb[, bin := pmin(nbin - 1L, pmax(0L, floor(rel * nbin)))]

dens <- rb[, .N, by = bin][order(bin)]
dens[, rel_mid := (bin + 0.5) / nbin]
dens[, density := N / sum(N)]

#save density table
fwrite(
  dens,
  file.path(vis_dir, "gene_body_density.tsv"),
  sep = "\t"
)

#plot gene-body relative position density
g3 <- ggplot(dens, aes(x = rel_mid, y = density)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Relative position along gene body (0 = TSS, 1 = TES)",
    y = "Density of read midpoints",
    title = "Gene-body relative position density of PacBio FLNC reads"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(vis_dir, "gene_body_density_curve.png"),
  g3, width = 7, height = 4.5, dpi = 300
)


#pigeon report was generated using the scripts provided on SQANTI3 GitHub (https://github.com/ConesaLab/SQANTI3/tree/master/utilities/report_pigeon)

####end of the session########################
####author: Luozixian Wang, Raymond Wong######
####18/06/2026################################
