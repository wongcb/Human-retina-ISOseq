########Iso-seq profiling of novel isoform variants in human retina at single cell resolution############
########created by Luozixian Wang, Raymond Wong##########################################################
########CERA, UNIMELB, 20/06/2024########################################################################

#load the packages
library(Biostrings)
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
theme_set(theme_cowplot())
set.seed(12345)
library(Matrix)
library(SingleR)
library(scran)
library(SingleCellExperiment)
library(scmap)
library(singleCellNet)
library(WGCNA)
library(igraph)
library(hdWGCNA)
enableWGCNAThreads(nThreads = 14)
library(qlcMatrix)
library(magrittr)
library(biomaRt)
library(org.Hs.eg.db)
library(ggrepel)
library(plotrix)
library(gridExtra)
library(AnnotationHub)
library(clusterProfiler)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(topGO)
library(proActiv)
library(data.table)
library(dunn.test)
library(ggplot2)
library(ggtranscript)

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

#load the long read data as iso assay in the seurat object
iso_counts <- iso_retina_annotated@assays$iso@counts
iso_retina <- iso_counts
colnames_split <- strsplit(colnames(iso_retina), "-", fixed = TRUE)
colnames_no_suffix <- sapply(colnames_split, `[`, 1)
suffixes <- sapply(colnames_split, function(x) if (length(x) > 1) paste0("-", x[2]) else "")
colnames_rc <- sapply(DNAStringSet(colnames_no_suffix), function(x) as.character(reverseComplement(x)))
colnames_final <- paste0(colnames_rc, suffixes)
colnames(iso_retina) <- colnames_final
cell_barcodes_iso_list <- colnames(iso_retina)
check_cell_barcodes <- merged_retina_subset@meta.data$cell.bc
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

#read the gene/transcript annotation
transcript_annotation <- read.csv("seurat_pb_output/RetinaMerge.annotated.info.csv", header = TRUE, sep = "\t")
transcript_status <- subset(transcript_annotation, category %in% c("full-splice_match", "incomplete-splice_match", "novel_not_in_catalog", "novel_in_catalog"))

#subset the transcript status list
FSM_transcript <- transcript_status %>% filter(category == "full-splice_match")
ISM_transcript <- transcript_status %>% filter(category == "incomplete-splice_match")
NNC_transcript <- transcript_status %>% filter(category == "novel_not_in_catalog")
NIC_transcript <- transcript_status %>% filter(category == "novel_in_catalog")

gene_FSM <- data.frame(gene_FSM = unique(FSM_transcript$gene), stringsAsFactors = FALSE)
gene_ISM <- data.frame(gene_ISM = unique(ISM_transcript$gene), stringsAsFactors = FALSE)
gene_NNC <- data.frame(gene_NNC = unique(NNC_transcript$gene), stringsAsFactors = FALSE)
gene_NIC <- data.frame(gene_NIC = unique(NIC_transcript$gene), stringsAsFactors = FALSE)

#give the status of each cell
#subset the NNC_transcript and NIC_transcript to only BCrev
FSM_BCrev <- FSM_transcript %>%
  select(BCrev)
ISM_BCrev <- ISM_transcript %>%
  select(BCrev)
NNC_BCrev <- NNC_transcript %>%
  select(BCrev)
NIC_BCrev <- NIC_transcript %>%
  select(BCrev)

FSM_BCrev$status <- "FSM"
ISM_BCrev$status <- "ISM"
NNC_BCrev$status <- "NNC"
NIC_BCrev$status <- "NIC"

unique_FSM <- FSM_BCrev %>%
  distinct(BCrev, .keep_all = TRUE)
unique_ISM <- ISM_BCrev %>%
  distinct(BCrev, .keep_all = TRUE)
unique_NNC <- NNC_BCrev %>%
  distinct(BCrev, .keep_all = TRUE)
unique_NIC <- NIC_BCrev %>%
  distinct(BCrev, .keep_all = TRUE)

intersected_data <- unique_ISM %>%
  inner_join(unique_NNC, by = "BCrev") %>%
  inner_join(unique_NIC, by = "BCrev")

intersected_data <- intersected_data %>%
  select(BCrev) %>%
  mutate(status = "NOVEL")

df_FSM_filtered <- unique_FSM %>%
  filter(!(BCrev %in% intersected_data$BCrev))
df_ISM_filtered <- unique_ISM %>%
  filter(!(BCrev %in% intersected_data$BCrev))
df_NNC_filtered <- unique_NNC %>%
  filter(!(BCrev %in% intersected_data$BCrev))
df_NIC_filtered <- unique_NIC %>%
  filter(!(BCrev %in% intersected_data$BCrev))
cell_status <- bind_rows(intersected_data, df_FSM_filtered, df_ISM_filtered, df_NNC_filtered, df_NIC_filtered)

retina_combined$novel_status <- "Novel Transcripts"
fsm_barcodes <- cell_status$BCrev[cell_status$status == "FSM"]
retina_combined$novel_status[match(fsm_barcodes, retina_combined$cell.bc.trimmed)] <- "Known Transcripts"

intersected_barcodes <- cell_status$BCrev[cell_status$status == "NOVEL"]
retina_combined$novel_status[match(intersected_barcodes, retina_combined$cell.bc.trimmed)] <- "Novel Transcripts"

for (i in 1:nrow(cell_status)) {
  barcode <- cell_status$BCrev[i]
  status <- cell_status$status[i]
  retina_combined$novel_status[retina_combined$cell.bc.trimmed %in% barcode] <- status
}

##plot figures
#figure 1E
DimPlot(retina_combined, group.by ='novel_status', reduction='umap', pt.size = 0.3)

#Supplementary figure 1A
DimPlot(retina_combined, reduction = "umap", group.by = "batch", pt.size = 0.3)

#Supplementary figure 1B
DimPlot(retina_combined, group.by ='read_stat', reduction='umap', pt.size = 0.3)

#Figure 1B
DimPlot(retina_combined, group.by ='cell.type.integrated', reduction='umap', pt.size = 0.3)

#Figure 1C
retina_marker<- read.csv(file = './Retinamarker_reorder.csv')

DotPlot(retina_combined, features = retina_marker, dot.scale = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=16)) +
  theme(axis.text.y = element_text(size=16)) +
    scale_color_viridis(option="plasma", direction = -1)

#Figure 4A
VlnPlot(retina_combined, feature = 'ABCA4', pt.size = 0, ncol = 1)

##GO term enrichment
symbol_ensid <- read.delim('symbol_ensid.txt', header = TRUE, stringsAsFactors = FALSE)
selected_symbol_ensid <- symbol_ensid[, c("HGNC.symbol", "Gene.stable.ID")]
filtered_symbol_ensid <- selected_symbol_ensid %>%
  group_by(HGNC.symbol) %>%
  slice(1) %>%
  ungroup()
NNC_converted <- gene_NNC %>%
  left_join(filtered_symbol_ensid, by = c(gene_NNC = "HGNC.symbol"))
NIC_converted <- gene_NIC %>%
  left_join(filtered_symbol_ensid, by = c(gene_NIC = "HGNC.symbol"))
FSM_converted <- gene_FSM %>%
  left_join(filtered_symbol_ensid, by = c(gene_FSM = "HGNC.symbol"))
ISM_converted <- gene_ISM %>%
  left_join(filtered_symbol_ensid, by = c(gene_ISM = "HGNC.symbol"))

FSM_filled <- FSM_converted %>%
  mutate(Gene.stable.ID = ifelse(is.na(Gene.stable.ID), gene_FSM, Gene.stable.ID))
ISM_filled <- ISM_converted %>%
  mutate(Gene.stable.ID = ifelse(is.na(Gene.stable.ID), gene_ISM, Gene.stable.ID))
NNC_filled <- NNC_converted %>%
  mutate(Gene.stable.ID = ifelse(is.na(Gene.stable.ID), gene_NNC, Gene.stable.ID))
NIC_filled <- NIC_converted %>%
  mutate(Gene.stable.ID = ifelse(is.na(Gene.stable.ID), gene_NIC, Gene.stable.ID))

#check duplicate
NNC_duplicated_genes <- NNC_filled %>%
  filter(duplicated(Gene.stable.ID) | duplicated(Gene.stable.ID, fromLast = TRUE))
print(NNC_duplicated_genes)
NIC_duplicated_genes <- NIC_filled %>%
  filter(duplicated(Gene.stable.ID) | duplicated(Gene.stable.ID, fromLast = TRUE))
print(NIC_duplicated_genes)
FSM_duplicated_genes <- FSM_filled %>%
  filter(duplicated(Gene.stable.ID) | duplicated(Gene.stable.ID, fromLast = TRUE))
print(FSM_duplicated_genes)
ISM_duplicated_genes <- ISM_filled %>%
  filter(duplicated(Gene.stable.ID) | duplicated(Gene.stable.ID, fromLast = TRUE))
print(ISM_duplicated_genes)

#GO enrichment
NNC_GO <- enrichGO(gene = NNC_filled$Gene.stable.ID, 
                   ont ="ALL", 
                   keyType = "ENSEMBL", 
                   pvalueCutoff = 0.05,
                   OrgDb = "org.Hs.eg.db", 
                   pAdjustMethod = "fdr", 
                   qvalueCutoff = 0.2, 
                   readable = TRUE)
NIC_GO <- enrichGO(gene = NIC_filled$Gene.stable.ID, 
                   ont ="ALL", 
                   keyType = "ENSEMBL", 
                   pvalueCutoff = 0.05,
                   OrgDb = "org.Hs.eg.db", 
                   pAdjustMethod = "fdr", 
                   qvalueCutoff = 0.2, 
                   readable = TRUE)
FSM_GO <- enrichGO(gene = FSM_filled$Gene.stable.ID, 
                   ont ="ALL", 
                   keyType = "ENSEMBL", 
                   pvalueCutoff = 0.05,
                   OrgDb = "org.Hs.eg.db", 
                   pAdjustMethod = "fdr", 
                   qvalueCutoff = 0.2, 
                   readable = TRUE)
ISM_GO <- enrichGO(gene = ISM_filled$Gene.stable.ID, 
                   ont ="ALL", 
                   keyType = "ENSEMBL", 
                   pvalueCutoff = 0.05,
                   OrgDb = "org.Hs.eg.db", 
                   pAdjustMethod = "fdr", 
                   qvalueCutoff = 0.2, 
                   readable = TRUE)


#plot GO
pNNC_dot <- dotplot(NNC_GO, showCategory = 10)
pNIC_dot <- dotplot(NIC_GO, showCategory = 10)
pFSM_dot <- dotplot(FSM_GO, showCategory = 10)
pISM_dot <- dotplot(ISM_GO, showCategory = 10)

#Supplementary figure 1C-F
print(pNNC_dot)
print(pNIC_dot)
print(pFSM_dot)
print(pISM_dot)


#alternative promoter analysis using proActiv
#before all, run a split_bam.py to split the bam file generated by isoseq3

# coding: utf-8
import pysam
import csv
#load the barcode-celltype
barcode_to_celltype = {}
with open("rc_bc_celltype.csv", "r") as csv_file:
    reader = csv.reader(csv_file)
    for row in reader:
        barcode_to_celltype[row[0]] = row[1]
#import the original bam file
input_bam = pysam.AlignmentFile("Retina_merged.dedup.mapped.bam", "rb")
#split bam files according to corresponding cell types
celltype_to_bamfile = {}
for celltype in set(barcode_to_celltype.values()):
    celltype_to_bamfile[celltype] = pysam.AlignmentFile(f"{celltype}.bam", "wb", header=input_bam.header)
for read in input_bam:
    #find the tag of CB
    if read.has_tag("CB"):
        barcode = read.get_tag("CB")
        if barcode in barcode_to_celltype:
            celltype = barcode_to_celltype[barcode]
            celltype_to_bamfile[celltype].write(read)
input_bam.close()
for bamfile in celltype_to_bamfile.values():
    bamfile.close()


#run proActiv
files <- c("/RayMerged/Amacrine.bam", 
           "/RayMerged/Astrocyte.bam", 
           "/RayMerged/Bipolar.bam", 
           "/RayMerged/Cone.bam", 
           "/RayMerged/Fibroblast.bam", 
           "/RayMerged/Horizontal.bam",
           "/RayMerged/Muller_Glia.bam",
           "/RayMerged/RPE.bam", 
           "/RayMerged/Rod.bam", 
           "/RayMerged/Smooth_Muscle_Cells.bam")
promoterAnnotation <- promoterAnnotation.gencode.v44.subset
condition <- c('Amacrine', 'Astrocyte', 'Bipolar', 'Cone', 'Fibroblast', 'Horizontal', 'Muller_Glia', 'RPE', 'Rod', 'Smooth_Muscle_Cells')
result <- proActiv(files = files,
                   promoterAnnotation = promoterAnnotation.gencode.v43.subset,
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

#re-generate the dotplot
EXOC1_promoter_activity_dot <- dotplotPromoters(result, "ENSG00000090989", filterInternal = F)
USH2A_promoter_activity_dot <- dotplotPromoters(result, "ENSG00000042781", filterInternal = F)
AHI1_promoter_activity_dot <- dotplotPromoters(result, "ENSG00000135541", filterInternal = F)
#figure 2C-E
grid.arrange(EXOC1_promoter_activity_dot[[1]], EXOC1_promoter_activity_dot[[2]], nrow = 1, ncol = 2, widths = c(2, 2))
grid.arrange(USH2A_promoter_activity_dot[[1]], USH2A_promoter_activity_dot[[2]], nrow = 1, ncol = 2, widths = c(2, 2))
grid.arrange(AHI1_promoter_activity_dot[[1]], AHI1_promoter_activity_dot[[2]], nrow = 1, ncol = 2, widths = c(2, 2))

## Create a long dataframe summarizing cell line and promoter class
pdata1 <- data.frame(cellLine = rep(c('Amacrine', 'Astrocyte', 'Bipolar', 'Cone', 'Fibroblast', 'Horizontal', 'Muller_Glia', 'RPE', 'Rod', 'Smooth_Muscle_Cells'), each = nrow(rdata)), promoterClass = as.factor(c(rdata$Amacrine.class, rdata$Astrocyte.class, rdata$Bipolar.class, rdata$Cone.class, rdata$Fibroblast.class, rdata$Horizontal.class, rdata$Muller_Glia.class, rdata$RPE.class, rdata$Rod.class, rdata$Smooth_Muscle_Cells.class)))
ggplot(na.omit(pdata1)) +
  geom_bar(aes(x = cellLine, fill = promoterClass)) + 
  xlab('Cell Lines') + ylab('Count') +  labs(fill = 'Promoter Category') +
  ggtitle('Categorization of Promoters')

## Because many genes have many annotated promoters, we collapse promoters from the 5th position and onward into one group for simplicity

# plot Figure 2B
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


# Figure 2A: proportion of genes using single or multiple promoters 
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

#Supplementary figure 5
## Get active major promoters of cell types
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


#Figure 3A-F
#read the gene list generated from retnet
retnet_list <- read.table("RetNet_category.txt", header = TRUE, sep = "\t")
#remove the duplicated gene
retnet_disease_unique <- lapply(retnet_list, unique)
#extract the unique gene list
retnet_genes <- as.vector(unlist(retnet_list))
retnet_gene_unique <- unique(retnet_genes)
retnet_gene_unique <- data.frame(retnet_gene_unique)
#find the corresponding pbid
retnet_pbid <- transcript_annotation[transcript_annotation$gene %in% 
                                       retnet_gene_unique$retnet_gene_unique, c("gene", "pbid")]
head(retnet_pbid)

retnet_pbid <- retnet_pbid[!duplicated(retnet_pbid), ]
rownames(retnet_pbid) <- NULL
head(retnet_pbid)

mapped_retnet <- unique(retnet_pbid$gene)
num_mapped_retnet <- length(mapped_retnet)
print(num_mapped_retnet)

#extract the count matrix from seurat object for heatmap
iso_counts <- retina_combined@assays$iso@counts
selected_pbid <- as.vector(retnet_pbid$pbid)
valid_pbid <- selected_pbid[selected_pbid %in% row.names(iso_counts)]
sub_iso_matrix <- iso_counts[valid_pbid,]

rna_counts <- retina_combined@assays$RNA@counts
selected_gene <- as.vector(retnet_pbid$gene)
selected_gene_unique <- unique(selected_gene)
valid_gene <- selected_gene_unique[selected_gene_unique %in% row.names(rna_counts)]
sub_rna_matrix <- rna_counts[valid_gene,]
regular_rna_matrix <- as.matrix(sub_rna_matrix)

#subset genes according to diseases
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

#RP: Retinitis pigmentosa
RP_pbid <- transcript_annotation[transcript_annotation$gene %in% 
                                   unique_RP$gene, c("gene", "pbid", "transcript", "category")]
head(RP_pbid)
RP_pbid <- RP_pbid %>%
  mutate(status = case_when(
    category == "full-splice_match" ~ "known",                                     
    category %in% c("incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog") ~ "novel", 
    TRUE ~ "other"                                                    
  ))
RP_pbid$transcript  <- NULL
RP_pbid <- RP_pbid[!duplicated(RP_pbid), ]
rownames(RP_pbid) <- NULL
head(RP_pbid)
RP_pbid$transcript_status <- paste(RP_pbid$gene, RP_pbid$status, sep = "_")
RP_pbid$gene <- NULL
RP_pbid$status <- NULL
RP_selected_pbid <- as.vector(RP_pbid$pbid)
RP_valid_pbid <- RP_selected_pbid[RP_selected_pbid %in% row.names(iso_counts)]
sub_RP_iso_matrix <- iso_counts[RP_valid_pbid,]
sub_RP_iso_df <- as.data.frame(sub_RP_iso_matrix)
sub_RP_iso_df$pbid <- rownames(sub_RP_iso_df)
sub_RP_iso_df <- left_join(sub_RP_iso_df, RP_pbid, by="pbid")
unique_RP_iso_count <- length(unique(sub_RP_iso_df$transcript_status))
print(paste("Number of unique values:", unique_RP_iso_count))
RP_grouped_df <- sub_RP_iso_df %>%
  group_by(transcript_status) %>%
  summarise(across(where(is.numeric), sum))
RP_final_matrix <- as.matrix(RP_grouped_df[, -1])
rownames(RP_final_matrix) <- RP_grouped_df$transcript_status
#sort the cells according to celltype
RP_final_sorted_matrix <- RP_final_matrix[, sorted_barcodes]
ncol(RP_final_matrix)
length(sorted_barcodes)
all(sorted_barcodes %in% colnames(RP_final_matrix))
mismatched_barcodes <- sorted_barcodes[!sorted_barcodes %in% colnames(RP_final_matrix)]
print(mismatched_barcodes)
#then split the transcript status
RP_row_names <- rownames(RP_final_sorted_matrix)
RP_sorting_vector <- ifelse(grepl("novel", RP_row_names), 1, 2)
RP_final_sorted_matrix <- RP_final_sorted_matrix[order(RP_sorting_vector), ]

RP_cluster_col <- ifelse(grepl("novel", RP_row_names), "novel", 
                         ifelse(grepl("known", RP_row_names), "known", NA))
RP_TS_col <- data.frame(transcript_status = RP_cluster_col, row.names = RP_row_names)
head(RP_TS_col)

#organise into seurat objects for dotplot
RP_seurat <- CreateSeuratObject(counts = RP_final_sorted_matrix)
RP_seurat[["celltype"]] <- barcode_data$cell.type.integrated[match(colnames(RP_seurat), barcode_data$barcode)]
RP_TS_col_new <- RP_TS_col
RP_seurat[["transcript_status"]] <- RP_TS_col_new$transcript_status[rownames(RP_seurat@meta.data)]
RP_seurat <- NormalizeData(RP_seurat)
RP_seurat.all.genes <- rownames(RP_seurat)
RP_seurat <- ScaleData(RP_seurat, features = RP_seurat.all.genes)
RP_new_row_names <- gsub(pattern = "_", replacement = "-", x = rownames(RP_TS_col))
rownames(RP_TS_col_new) <- RP_new_row_names
RP_TS_col_new <- data.frame(gene = rownames(RP_TS_col_new), RP_TS_col_new, row.names = NULL)
RP_sorted_genes_df <- RP_TS_col_new[order(RP_TS_col_new$transcript_status, decreasing = TRUE), ]
Idents(RP_seurat) <- "celltype"
clusters <- SplitObject(RP_seurat, split.by = "celltype") 

RP_novel_gene <- subset(RP_sorted_genes_df, transcript_status == "novel")

all.RP_novel_dotplot <- DotPlot(RP_seurat, features = RP_novel_gene$gene, dot.scale = 10, scale.by = 'size', col.min = 0, scale.min = 1) + 
  coord_flip() + 
  theme(legend.position = "bottom", legend.justification = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 8), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))+
  ggtitle("Retinits pigmentosa" )+
  theme(plot.title = element_text(size=25, hjust = 0.5))  

#CORD
CORD_pbid <- transcript_annotation[transcript_annotation$gene %in% 
                                   unique_CORD$gene, c("gene", "pbid", "transcript", "category")]
head(CORD_pbid)
CORD_pbid <- CORD_pbid %>%
  mutate(status = case_when(
    category == "full-splice_match" ~ "known",                                     
    category %in% c("incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog") ~ "novel", 
    TRUE ~ "other"                                                    
  ))
CORD_pbid$transcript  <- NULL
CORD_pbid <- CORD_pbid[!duplicated(CORD_pbid), ]
rownames(CORD_pbid) <- NULL
head(CORD_pbid)
CORD_pbid$transcript_status <- paste(CORD_pbid$gene, CORD_pbid$status, sep = "_")
CORD_pbid$gene <- NULL
CORD_pbid$status <- NULL
CORD_selected_pbid <- as.vector(CORD_pbid$pbid)
CORD_valid_pbid <- CORD_selected_pbid[CORD_selected_pbid %in% row.names(iso_counts)]
sub_CORD_iso_matrix <- iso_counts[CORD_valid_pbid,]
sub_CORD_iso_df <- as.data.frame(sub_CORD_iso_matrix)
sub_CORD_iso_df$pbid <- rownames(sub_CORD_iso_df)
sub_CORD_iso_df <- left_join(sub_CORD_iso_df, CORD_pbid, by="pbid")
unique_CORD_iso_count <- length(unique(sub_CORD_iso_df$transcript_status))
print(paste("Number of unique values:", unique_CORD_iso_count))
CORD_grouped_df <- sub_CORD_iso_df %>%
  group_by(transcript_status) %>%
  summarise(across(where(is.numeric), sum))
CORD_final_matrix <- as.matrix(CORD_grouped_df[, -1])
rownames(CORD_final_matrix) <- CORD_grouped_df$transcript_status
#sort the cells according to celltype
CORD_final_sorted_matrix <- CORD_final_matrix[, sorted_barcodes]
ncol(CORD_final_matrix)
length(sorted_barcodes)
all(sorted_barcodes %in% colnames(CORD_final_matrix))
mismatched_barcodes <- sorted_barcodes[!sorted_barcodes %in% colnames(CORD_final_matrix)]
print(mismatched_barcodes)
#then split the transcript status
CORD_row_names <- rownames(CORD_final_sorted_matrix)
CORD_sorting_vector <- ifelse(grepl("novel", CORD_row_names), 1, 2)
CORD_final_sorted_matrix <- CORD_final_sorted_matrix[order(CORD_sorting_vector), ]

CORD_cluster_col <- ifelse(grepl("novel", CORD_row_names), "novel", 
                         ifelse(grepl("known", CORD_row_names), "known", NA))
CORD_TS_col <- data.frame(transcript_status = CORD_cluster_col, row.names = CORD_row_names)
head(CORD_TS_col)

#try to organise into seurat objects for dotplot
CORD_seurat <- CreateSeuratObject(counts = CORD_final_sorted_matrix)
CORD_seurat[["celltype"]] <- barcode_data$cell.type.integrated[match(colnames(CORD_seurat), barcode_data$barcode)]
CORD_TS_col_new <- CORD_TS_col
CORD_seurat[["transcript_status"]] <- CORD_TS_col_new$transcript_status[rownames(CORD_seurat@meta.data)]
CORD_seurat <- NormalizeData(CORD_seurat)
CORD_seurat.all.genes <- rownames(CORD_seurat)
CORD_seurat <- ScaleData(CORD_seurat, features = CORD_seurat.all.genes)
CORD_new_row_names <- gsub(pattern = "_", replacement = "-", x = rownames(CORD_TS_col))
rownames(CORD_TS_col_new) <- CORD_new_row_names
CORD_TS_col_new <- data.frame(gene = rownames(CORD_TS_col_new), CORD_TS_col_new, row.names = NULL)
CORD_sorted_genes_df <- CORD_TS_col_new[order(CORD_TS_col_new$transcript_status, decreasing = TRUE), ]
Idents(CORD_seurat) <- "celltype"
clusters <- SplitObject(CORD_seurat, split.by = "celltype")

CORD_novel_gene <- subset(CORD_sorted_genes_df, transcript_status == "novel")

all.CORD_novel_dotplot <- DotPlot(CORD_seurat, features = CORD_novel_gene$gene, dot.scale = 10, scale.by = 'size', col.min = 0, scale.min = 1) + 
  coord_flip() + 
  theme(legend.position = "bottom", legend.justification = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 8), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) +
  ggtitle("CORD")+
  theme(plot.title = element_text(size=25, hjust = 0.5))  

#MD
MD_pbid <- transcript_annotation[transcript_annotation$gene %in% 
                                     unique_MD$gene, c("gene", "pbid", "transcript", "category")]
head(MD_pbid)
MD_pbid <- MD_pbid %>%
  mutate(status = case_when(
    category == "full-splice_match" ~ "known",                                     
    category %in% c("incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog") ~ "novel", 
    TRUE ~ "other"                                                    
  ))
MD_pbid$transcript  <- NULL
MD_pbid <- MD_pbid[!duplicated(MD_pbid), ]
rownames(MD_pbid) <- NULL
head(MD_pbid)
MD_pbid$transcript_status <- paste(MD_pbid$gene, MD_pbid$status, sep = "_")
MD_pbid$gene <- NULL
MD_pbid$status <- NULL
MD_selected_pbid <- as.vector(MD_pbid$pbid)
MD_valid_pbid <- MD_selected_pbid[MD_selected_pbid %in% row.names(iso_counts)]
sub_MD_iso_matrix <- iso_counts[MD_valid_pbid,]
sub_MD_iso_df <- as.data.frame(sub_MD_iso_matrix)
sub_MD_iso_df$pbid <- rownames(sub_MD_iso_df)
sub_MD_iso_df <- left_join(sub_MD_iso_df, MD_pbid, by="pbid")
unique_MD_iso_count <- length(unique(sub_MD_iso_df$transcript_status))
print(paste("Number of unique values:", unique_MD_iso_count))
MD_grouped_df <- sub_MD_iso_df %>%
  group_by(transcript_status) %>%
  summarise(across(where(is.numeric), sum))
MD_final_matrix <- as.matrix(MD_grouped_df[, -1])
rownames(MD_final_matrix) <- MD_grouped_df$transcript_status
#sort the cells according to celltype
MD_final_sorted_matrix <- MD_final_matrix[, sorted_barcodes]
ncol(MD_final_matrix)
length(sorted_barcodes)
all(sorted_barcodes %in% colnames(MD_final_matrix))
mismatched_barcodes <- sorted_barcodes[!sorted_barcodes %in% colnames(MD_final_matrix)]
print(mismatched_barcodes)
#then split the transcript status
MD_row_names <- rownames(MD_final_sorted_matrix)
MD_sorting_vector <- ifelse(grepl("novel", MD_row_names), 1, 2)
MD_final_sorted_matrix <- MD_final_sorted_matrix[order(MD_sorting_vector), ]

MD_cluster_col <- ifelse(grepl("novel", MD_row_names), "novel", 
                           ifelse(grepl("known", MD_row_names), "known", NA))
MD_TS_col <- data.frame(transcript_status = MD_cluster_col, row.names = MD_row_names)
head(MD_TS_col)

#try to organise into seurat objects for dotplot
MD_seurat <- CreateSeuratObject(counts = MD_final_sorted_matrix)
MD_seurat[["celltype"]] <- barcode_data$cell.type.integrated[match(colnames(MD_seurat), barcode_data$barcode)]
MD_TS_col_new <- MD_TS_col
MD_seurat[["transcript_status"]] <- MD_TS_col_new$transcript_status[rownames(MD_seurat@meta.data)]
MD_seurat <- NormalizeData(MD_seurat)
MD_seurat.all.genes <- rownames(MD_seurat)
MD_seurat <- ScaleData(MD_seurat, features = MD_seurat.all.genes)
MD_new_row_names <- gsub(pattern = "_", replacement = "-", x = rownames(MD_TS_col))
rownames(MD_TS_col_new) <- MD_new_row_names
MD_TS_col_new <- data.frame(gene = rownames(MD_TS_col_new), MD_TS_col_new, row.names = NULL)
MD_sorted_genes_df <- MD_TS_col_new[order(MD_TS_col_new$transcript_status, decreasing = TRUE), ]
Idents(MD_seurat) <- "celltype"
clusters <- SplitObject(MD_seurat, split.by = "celltype")

MD_novel_gene <- subset(MD_sorted_genes_df, transcript_status == "novel")

all.MD_novel_dotplot <- DotPlot(MD_seurat, features = MD_novel_gene$gene, dot.scale = 10, scale.by = 'size', col.min = 0, scale.min = 1) + 
  coord_flip() + 
  theme(legend.position = "bottom", legend.justification = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 8), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))+
  ggtitle("Inherited macular degeneration")+
  theme(plot.title = element_text(size=25, hjust = 0.5))  


#CSNB
CSNB_pbid <- transcript_annotation[transcript_annotation$gene %in% 
                                     unique_CSNB$gene, c("gene", "pbid", "transcript", "category")]
head(CSNB_pbid)
CSNB_pbid <- CSNB_pbid %>%
  mutate(status = case_when(
    category == "full-splice_match" ~ "known",                                     
    category %in% c("incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog") ~ "novel", 
    TRUE ~ "other"                                                    
  ))
CSNB_pbid$transcript  <- NULL
CSNB_pbid <- CSNB_pbid[!duplicated(CSNB_pbid), ]
rownames(CSNB_pbid) <- NULL
head(CSNB_pbid)
CSNB_pbid$transcript_status <- paste(CSNB_pbid$gene, CSNB_pbid$status, sep = "_")
CSNB_pbid$gene <- NULL
CSNB_pbid$status <- NULL
CSNB_selected_pbid <- as.vector(CSNB_pbid$pbid)
CSNB_valid_pbid <- CSNB_selected_pbid[CSNB_selected_pbid %in% row.names(iso_counts)]
sub_CSNB_iso_matrix <- iso_counts[CSNB_valid_pbid,]
sub_CSNB_iso_df <- as.data.frame(sub_CSNB_iso_matrix)
sub_CSNB_iso_df$pbid <- rownames(sub_CSNB_iso_df)
sub_CSNB_iso_df <- left_join(sub_CSNB_iso_df, CSNB_pbid, by="pbid")
unique_CSNB_iso_count <- length(unique(sub_CSNB_iso_df$transcript_status))
print(paste("Number of unique values:", unique_CSNB_iso_count))
CSNB_grouped_df <- sub_CSNB_iso_df %>%
  group_by(transcript_status) %>%
  summarise(across(where(is.numeric), sum))
CSNB_final_matrix <- as.matrix(CSNB_grouped_df[, -1])
rownames(CSNB_final_matrix) <- CSNB_grouped_df$transcript_status
#sort the cells acCSNBing to celltype
CSNB_final_sorted_matrix <- CSNB_final_matrix[, sorted_barcodes]
ncol(CSNB_final_matrix)
length(sorted_barcodes)
all(sorted_barcodes %in% colnames(CSNB_final_matrix))
mismatched_barcodes <- sorted_barcodes[!sorted_barcodes %in% colnames(CSNB_final_matrix)]
print(mismatched_barcodes)
#then split the transcript status
CSNB_row_names <- rownames(CSNB_final_sorted_matrix)
CSNB_sorting_vector <- ifelse(grepl("novel", CSNB_row_names), 1, 2)
CSNB_final_sorted_matrix <- CSNB_final_sorted_matrix[order(CSNB_sorting_vector), ]

CSNB_cluster_col <- ifelse(grepl("novel", CSNB_row_names), "novel", 
                           ifelse(grepl("known", CSNB_row_names), "known", NA))
CSNB_TS_col <- data.frame(transcript_status = CSNB_cluster_col, row.names = CSNB_row_names)
head(CSNB_TS_col)

#try to organise into seurat objects for dotplot
CSNB_seurat <- CreateSeuratObject(counts = CSNB_final_sorted_matrix)
CSNB_seurat[["celltype"]] <- barcode_data$cell.type.integrated[match(colnames(CSNB_seurat), barcode_data$barcode)]
CSNB_TS_col_new <- CSNB_TS_col
CSNB_seurat[["transcript_status"]] <- CSNB_TS_col_new$transcript_status[rownames(CSNB_seurat@meta.data)]
CSNB_seurat <- NormalizeData(CSNB_seurat)
CSNB_seurat.all.genes <- rownames(CSNB_seurat)
CSNB_seurat <- ScaleData(CSNB_seurat, features = CSNB_seurat.all.genes)
CSNB_new_row_names <- gsub(pattern = "_", replacement = "-", x = rownames(CSNB_TS_col))
rownames(CSNB_TS_col_new) <- CSNB_new_row_names
CSNB_TS_col_new <- data.frame(gene = rownames(CSNB_TS_col_new), CSNB_TS_col_new, row.names = NULL)
CSNB_sorted_genes_df <- CSNB_TS_col_new[order(CSNB_TS_col_new$transcript_status, decreasing = TRUE), ]
Idents(CSNB_seurat) <- "celltype"
clusters <- SplitObject(CSNB_seurat, split.by = "celltype")

CSNB_novel_gene <- subset(CSNB_sorted_genes_df, transcript_status == "novel")

all.CSNB_novel_dotplot <- DotPlot(CSNB_seurat, features = CSNB_novel_gene$gene, dot.scale = 10, scale.by = 'size', col.min = 0, scale.min = 1, scale.max = 30) + 
  coord_flip() + 
  theme(legend.position = "bottom", legend.justification = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 8), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))+
  ggtitle("CSNB" )+
  theme(plot.title = element_text(size=25, hjust = 0.5))  


#LCA
LCA_pbid <- transcript_annotation[transcript_annotation$gene %in% 
                                     unique_LCA$gene, c("gene", "pbid", "transcript", "category")]
head(LCA_pbid)
LCA_pbid <- LCA_pbid %>%
  mutate(status = case_when(
    category == "full-splice_match" ~ "known",                                     
    category %in% c("incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog") ~ "novel", 
    TRUE ~ "other"                                                    
  ))
LCA_pbid$transcript  <- NULL
LCA_pbid <- LCA_pbid[!duplicated(LCA_pbid), ]
rownames(LCA_pbid) <- NULL
head(LCA_pbid)
LCA_pbid$transcript_status <- paste(LCA_pbid$gene, LCA_pbid$status, sep = "_")
LCA_pbid$gene <- NULL
LCA_pbid$status <- NULL
LCA_selected_pbid <- as.vector(LCA_pbid$pbid)
LCA_valid_pbid <- LCA_selected_pbid[LCA_selected_pbid %in% row.names(iso_counts)]
sub_LCA_iso_matrix <- iso_counts[LCA_valid_pbid,]
sub_LCA_iso_df <- as.data.frame(sub_LCA_iso_matrix)
sub_LCA_iso_df$pbid <- rownames(sub_LCA_iso_df)
sub_LCA_iso_df <- left_join(sub_LCA_iso_df, LCA_pbid, by="pbid")
unique_LCA_iso_count <- length(unique(sub_LCA_iso_df$transcript_status))
print(paste("Number of unique values:", unique_LCA_iso_count))
LCA_grouped_df <- sub_LCA_iso_df %>%
  group_by(transcript_status) %>%
  summarise(across(where(is.numeric), sum))
LCA_final_matrix <- as.matrix(LCA_grouped_df[, -1])
rownames(LCA_final_matrix) <- LCA_grouped_df$transcript_status
#sort the cells acLCAing to celltype
LCA_final_sorted_matrix <- LCA_final_matrix[, sorted_barcodes]
ncol(LCA_final_matrix)
length(sorted_barcodes)
all(sorted_barcodes %in% colnames(LCA_final_matrix))
mismatched_barcodes <- sorted_barcodes[!sorted_barcodes %in% colnames(LCA_final_matrix)]
print(mismatched_barcodes)
#then split the transcript status
LCA_row_names <- rownames(LCA_final_sorted_matrix)
LCA_sorting_vector <- ifelse(grepl("novel", LCA_row_names), 1, 2)
LCA_final_sorted_matrix <- LCA_final_sorted_matrix[order(LCA_sorting_vector), ]

LCA_cluster_col <- ifelse(grepl("novel", LCA_row_names), "novel", 
                           ifelse(grepl("known", LCA_row_names), "known", NA))
LCA_TS_col <- data.frame(transcript_status = LCA_cluster_col, row.names = LCA_row_names)
head(LCA_TS_col)

#try to organise into seurat objects for dotplot
LCA_seurat <- CreateSeuratObject(counts = LCA_final_sorted_matrix)
LCA_seurat[["celltype"]] <- barcode_data$cell.type.integrated[match(colnames(LCA_seurat), barcode_data$barcode)]
LCA_TS_col_new <- LCA_TS_col
LCA_seurat[["transcript_status"]] <- LCA_TS_col_new$transcript_status[rownames(LCA_seurat@meta.data)]
LCA_seurat <- NormalizeData(LCA_seurat)
LCA_seurat.all.genes <- rownames(LCA_seurat)
LCA_seurat <- ScaleData(LCA_seurat, features = LCA_seurat.all.genes)
LCA_new_row_names <- gsub(pattern = "_", replacement = "-", x = rownames(LCA_TS_col))
rownames(LCA_TS_col_new) <- LCA_new_row_names
LCA_TS_col_new <- data.frame(gene = rownames(LCA_TS_col_new), LCA_TS_col_new, row.names = NULL)
LCA_sorted_genes_df <- LCA_TS_col_new[order(LCA_TS_col_new$transcript_status, decreasing = TRUE), ]
Idents(LCA_seurat) <- "celltype"
clusters <- SplitObject(LCA_seurat, split.by = "celltype")

LCA_novel_gene <- subset(LCA_sorted_genes_df, transcript_status == "novel")

all.LCA_novel_dotplot <- DotPlot(LCA_seurat, features = LCA_novel_gene$gene, dot.scale = 10, scale.by = 'size', col.min = 0, scale.min = 1) + 
  coord_flip() + 
  theme(legend.position = "bottom", legend.justification = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 8), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))+
  ggtitle("LCA" )+
  theme(plot.title = element_text(size=25, hjust = 0.5))  


#US
US_pbid <- transcript_annotation[transcript_annotation$gene %in% 
                                    unique_US$gene, c("gene", "pbid", "transcript", "category")]
head(US_pbid)
US_pbid <- US_pbid %>%
  mutate(status = case_when(
    category == "full-splice_match" ~ "known",                                     
    category %in% c("incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog") ~ "novel", 
    TRUE ~ "other"                                                    
  ))
US_pbid$transcript  <- NULL
US_pbid <- US_pbid[!duplicated(US_pbid), ]
rownames(US_pbid) <- NULL
head(US_pbid)
US_pbid$transcript_status <- paste(US_pbid$gene, US_pbid$status, sep = "_")
US_pbid$gene <- NULL
US_pbid$status <- NULL
US_selected_pbid <- as.vector(US_pbid$pbid)
US_valid_pbid <- US_selected_pbid[US_selected_pbid %in% row.names(iso_counts)]
sub_US_iso_matrix <- iso_counts[US_valid_pbid,]
sub_US_iso_df <- as.data.frame(sub_US_iso_matrix)
sub_US_iso_df$pbid <- rownames(sub_US_iso_df)
sub_US_iso_df <- left_join(sub_US_iso_df, US_pbid, by="pbid")
unique_US_iso_count <- length(unique(sub_US_iso_df$transcript_status))
print(paste("Number of unique values:", unique_US_iso_count))
US_grouped_df <- sub_US_iso_df %>%
  group_by(transcript_status) %>%
  summarise(across(where(is.numeric), sum))
US_final_matrix <- as.matrix(US_grouped_df[, -1])
rownames(US_final_matrix) <- US_grouped_df$transcript_status
#sort the cells acUSing to celltype
US_final_sorted_matrix <- US_final_matrix[, sorted_barcodes]
ncol(US_final_matrix)
length(sorted_barcodes)
all(sorted_barcodes %in% colnames(US_final_matrix))
mismatched_barcodes <- sorted_barcodes[!sorted_barcodes %in% colnames(US_final_matrix)]
print(mismatched_barcodes)
#then split the transcript status
US_row_names <- rownames(US_final_sorted_matrix)
US_sorting_vector <- ifelse(grepl("novel", US_row_names), 1, 2)
US_final_sorted_matrix <- US_final_sorted_matrix[order(US_sorting_vector), ]

US_cluster_col <- ifelse(grepl("novel", US_row_names), "novel", 
                          ifelse(grepl("known", US_row_names), "known", NA))
US_TS_col <- data.frame(transcript_status = US_cluster_col, row.names = US_row_names)
head(US_TS_col)

#try to organise into seurat objects for dotplot
US_seurat <- CreateSeuratObject(counts = US_final_sorted_matrix)
US_seurat[["celltype"]] <- barcode_data$cell.type.integrated[match(colnames(US_seurat), barcode_data$barcode)]
US_TS_col_new <- US_TS_col
US_seurat[["transcript_status"]] <- US_TS_col_new$transcript_status[rownames(US_seurat@meta.data)]
US_seurat <- NormalizeData(US_seurat)
US_seurat.all.genes <- rownames(US_seurat)
US_seurat <- ScaleData(US_seurat, features = US_seurat.all.genes)
US_new_row_names <- gsub(pattern = "_", replacement = "-", x = rownames(US_TS_col))
rownames(US_TS_col_new) <- US_new_row_names
US_TS_col_new <- data.frame(gene = rownames(US_TS_col_new), US_TS_col_new, row.names = NULL)
US_sorted_genes_df <- US_TS_col_new[order(US_TS_col_new$transcript_status, decreasing = TRUE), ]
Idents(US_seurat) <- "celltype"
clusters <- SplitObject(US_seurat, split.by = "celltype")

US_novel_gene <- subset(US_sorted_genes_df, transcript_status == "novel")

all.LCA_novel_dotplot <- DotPlot(LCA_seurat, features = LCA_novel_gene$gene, dot.scale = 10, scale.by = 'size', col.min = 0, scale.min = 1) + 
  coord_flip() + 
  theme(legend.position = "bottom", legend.justification = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 8), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))+
  ggtitle("LCA" )+
  theme(plot.title = element_text(size=25, hjust = 0.5))  

all.US_novel_dotplot <- DotPlot(US_seurat, features = US_novel_gene$gene, dot.scale = 10, scale.by = 'size', col.min = 0, scale.min = 1) + 
  coord_flip() + 
  theme(legend.position = "bottom", legend.justification = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 8), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))+
  ggtitle("Usher Syndrome" )+
  theme(plot.title = element_text(size=25, hjust = 0.5))  

#patchwork Figure 3A-F

jpeg(filename = "Figure3.resize.Dotplot.IRDgenes.jpg", width = 6000, height = 7200, res = 300)
 (all.CORD_novel_dotplot / all.MD_novel_dotplot) | (all.CSNB_novel_dotplot / all.LCA_novel_dotplot / all.US_novel_dotplot ) | all.RP_novel_dotplot
dev.off()

#ggtranscript to plot the structure of transcripts
#read the annotation file
#read the file for filtering
pb_anno <- rtracklayer::import("RetinaMerge.sorted.filtered_lite.gff")
class(pb_anno)
pb_anno <- pb_anno %>% dplyr::as_tibble()
class(pb_anno)

reference_anno <- rtracklayer::import("gencode.v44.chr_patch_hapl_scaff.annotation.gff3")
class(reference_anno)
reference_anno <- reference_anno %>% dplyr::as_tibble()
class(reference_anno)

#read the file for filtering
pb_infor <- read.csv("RetinaMerge_classification.filtered_lite_classification.txt", header = TRUE, sep = "\t")
#FL
pbid_fl <- pb_infor %>%
  filter(within_cage_peak == TRUE, !is.na(polyA_motif))
pbid_cage <- pb_infor %>% filter(within_cage_peak == TRUE)

#extract the count matrix and cell type information
#count the unique PBids for each cell type
#finally will introduce back to the ggtranscript for plotting
seurat_list <- SplitObject(retina_combined, split.by = "cell.type.integrated")
iso_counts_list <- lapply(seurat_list, function(x) {
  GetAssayData(x, assay = "isoform", slot = "counts")
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

#Figure 4B
#plot ABCA4
ABCA4_filtered_cage <- pbid_fl %>% 
  filter(associated_gene == "ABCA4")
# filter your gtf for the gene of interest
ABCA4_PBid_of_interest <- ABCA4_filtered_cage$isoform
ABCA4_annotation_from_gtf <- pb_anno %>%
  dplyr::filter(
    !is.na(transcript_id),
    transcript_id %in% ABCA4_PBid_of_interest
  )
# extract the required annotation columns
ABCA4_annotation_from_gtf <- ABCA4_annotation_from_gtf %>% 
  dplyr::select(
    seqnames,
    start,
    end,
    strand,
    type,
    gene_id,
    transcript_id
  )
ABCA4_annotation_from_gtf %>% head()
#unique(ABCA4_annotation_from_gtf$transcript_id)
ABCA4_annotation_from_gtf$gene_id
ABCA4_annotation_from_gtf$gene_id <- "ABCA4"

#reference
gene_of_interest <- "ABCA4"
ABCA4_reference_annotation_from_gtf <- reference_anno %>%
  dplyr::filter(
    !is.na(gene_name),
    gene_name %in% gene_of_interest
  )
ABCA4_reference_annotation_from_gtf <- ABCA4_reference_annotation_from_gtf %>% 
  dplyr::select(
    seqnames,
    start,
    end,
    strand,
    type,
    gene_name,
    transcript_name
  )

#plot for each cell type
colors <- RColorBrewer::brewer.pal(12, "Paired")
color_map <- setNames(colors, c("Bipolar", "Rod", "Amacrine", "Muller_Glia", "Fibroblast", 
                                "Horizontal", "RPE", "Cone", "Astrocyte", "Smooth Muscle Cells", "zReference"))

temp_df <- data.frame()
ABCA4_exons_plot <- data.frame()
for (cell_type in names(expressed_genes_by_type)) {
  transcript_ids_of_interest <- expressed_genes_by_type[[cell_type]]
  temp_df <- ABCA4_annotation_from_gtf %>%
    filter(transcript_id %in% transcript_ids_of_interest, type == "exon") %>%
    mutate(cell_type = cell_type, transcript_to_plot = paste(cell_type, transcript_id, sep = "."))
  ABCA4_exons_plot <- rbind(ABCA4_exons_plot, temp_df)
}
ABCA4_reference_exons <- ABCA4_reference_annotation_from_gtf %>%
  filter(type == "exon")
ABCA4_reference_exons <- ABCA4_reference_exons %>% 
  mutate(transcript_id = transcript_name, gene_id = gene_name, cell_type = "zReference", transcript_to_plot = paste(cell_type, transcript_id, sep = "."))
ABCA4_reference_exons$gene_name <- NULL
ABCA4_reference_exons$transcript_name <- NULL
ABCA4_combined_plot <- rbind(ABCA4_exons_plot, ABCA4_reference_exons)
ABCA4_together_plot <- ABCA4_combined_plot %>%
  ggplot(aes(xstart = start, xend = end, y = transcript_to_plot)) +
  geom_range(aes(fill = cell_type), color = color_map[cell_type]) + 
  geom_intron(
    data = to_intron(ABCA4_combined_plot, "transcript_to_plot"),
    aes(strand = strand), 
    arrow.min.intron.length = 1000
  ) +
  scale_fill_manual(values = color_map) + 
  labs(title = "ABCA4") +
  theme_minimal()
print(ABCA4_together_plot)

ABCA4plots <- list()
heights <- numeric()
ABCA4_reference_rescaled <- shorten_gaps(
  ABCA4_reference_exons, 
  to_intron(ABCA4_reference_exons, "transcript_id"), 
  group_var = "transcript_id"
)
ABCA4_reference_plot <- ABCA4_reference_rescaled %>%
  dplyr::filter(type == "exon") %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
  geom_range() +
  geom_intron(
    data = ABCA4_reference_rescaled %>% dplyr::filter(type == "intron"), 
    arrow.min.intron.length = 200
  )
for (cell_type in names(expressed_genes_by_type)) {
  transcript_ids_of_interest <- expressed_genes_by_type[[cell_type]]
  ABCA4_exons <- ABCA4_annotation_from_gtf %>%
    filter(transcript_id %in% transcript_ids_of_interest, type == "exon") %>%
    mutate(cell_type = cell_type)
  if (nrow(ABCA4_exons) > 0) {
    num_transcripts <- length(unique(ABCA4_exons$transcript_id))
    heights <- c(heights, num_transcripts)
    
    ABCA4_rescaled <- shorten_gaps(
      ABCA4_exons,
      to_intron(ABCA4_exons, "transcript_id"),
      group_var = "transcript_id"
    )
    plot <- ABCA4_rescaled %>%
      filter(type == "exon") %>%
      ggplot(aes(xstart = start, xend = end, y = transcript_id)) +
      geom_range(aes(fill = cell_type), color = color_map[cell_type]) + 
      geom_intron(
        data = ABCA4_rescaled %>% filter(type == "intron"),
        arrow.min.intron.length = 200
      ) +
      scale_fill_manual(values = color_map) + 
      labs(title = cell_type) +
      theme_minimal()
    ABCA4plots[[cell_type]] <- plot
  } else {
    ABCA4plots[[cell_type]] <- NULL
  }
}
#add reference
ABCA4plots[["reference"]] <- ABCA4_reference_plot
heights <- c(heights, 7)

ABCA4_valid_plots <- Filter(Negate(is.null), ABCA4plots)
ABCA4_valid_heights <- heights[!sapply(ABCA4plots, is.null)]
ABCA4_rescaled_plot <- wrap_plots(ABCA4_valid_plots, ncol = 1) + 
  plot_layout(heights = ABCA4_valid_heights / sum(ABCA4_valid_heights))
print(ABCA4_rescaled_plot)


#Supplementary figure 7: plot ASPH
ASPH_filtered_cage <- pbid_fl %>% 
  filter(associated_gene == "ASPH")
# filter your gtf for the gene of interest
ASPH_PBid_of_interest <- ASPH_filtered_cage$isoform
ASPH_annotation_from_gtf <- pb_anno %>%
  dplyr::filter(
    !is.na(transcript_id),
    transcript_id %in% ASPH_PBid_of_interest
  )
# extract the required annotation columns
ASPH_annotation_from_gtf <- ASPH_annotation_from_gtf %>% 
  dplyr::select(
    seqnames,
    start,
    end,
    strand,
    type,
    gene_id,
    transcript_id
  )
ASPH_annotation_from_gtf %>% head()
#unique(ASPH_annotation_from_gtf$transcript_id)
ASPH_annotation_from_gtf$gene_id
ASPH_annotation_from_gtf$gene_id <- "ASPH"

#reference
gene_of_interest <- "ASPH"
ASPH_reference_annotation_from_gtf <- reference_anno %>%
  dplyr::filter(
    !is.na(gene_name),
    gene_name %in% gene_of_interest
  )
ASPH_reference_annotation_from_gtf <- ASPH_reference_annotation_from_gtf %>% 
  dplyr::select(
    seqnames,
    start,
    end,
    strand,
    type,
    gene_name,
    transcript_name
  )

#plot for each cell type
colors <- RColorBrewer::brewer.pal(12, "Paired")
color_map <- setNames(colors, c("Bipolar", "Rod", "Amacrine", "Muller_Glia", "Fibroblast", 
                                "Horizontal", "RPE", "Cone", "Astrocyte", "Smooth Muscle Cells", "zReference"))

temp_df <- data.frame()
ASPH_exons_plot <- data.frame()
for (cell_type in names(expressed_genes_by_type)) {
  transcript_ids_of_interest <- expressed_genes_by_type[[cell_type]]
  temp_df <- ASPH_annotation_from_gtf %>%
    filter(transcript_id %in% transcript_ids_of_interest, type == "exon") %>%
    mutate(cell_type = cell_type, transcript_to_plot = paste(cell_type, transcript_id, sep = "."))
  ASPH_exons_plot <- rbind(ASPH_exons_plot, temp_df)
}
ASPH_reference_exons <- ASPH_reference_annotation_from_gtf %>%
  filter(type == "exon")
ASPH_reference_exons <- ASPH_reference_exons %>% 
  mutate(transcript_id = transcript_name, gene_id = gene_name, cell_type = "zReference", transcript_to_plot = paste(cell_type, transcript_id, sep = "."))
ASPH_reference_exons$gene_name <- NULL
ASPH_reference_exons$transcript_name <- NULL
ASPH_combined_plot <- rbind(ASPH_exons_plot, ASPH_reference_exons)
ASPH_together_plot <- ASPH_combined_plot %>%
  ggplot(aes(xstart = start, xend = end, y = transcript_to_plot)) +
  geom_range(aes(fill = cell_type), color = color_map[cell_type]) + 
  geom_intron(
    data = to_intron(ASPH_combined_plot, "transcript_to_plot"),
    aes(strand = strand), 
    arrow.min.intron.length = 1000
  ) +
  scale_fill_manual(values = color_map) + 
  labs(title = "ASPH") +
  theme_minimal()
print(ASPH_together_plot)

ASPHplots <- list()
heights <- numeric()
ASPH_reference_rescaled <- shorten_gaps(
  ASPH_reference_exons, 
  to_intron(ASPH_reference_exons, "transcript_id"), 
  group_var = "transcript_id"
)
ASPH_reference_plot <- ASPH_reference_rescaled %>%
  dplyr::filter(type == "exon") %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
  geom_range() +
  geom_intron(
    data = ASPH_reference_rescaled %>% dplyr::filter(type == "intron"), 
    arrow.min.intron.length = 200
  )
for (cell_type in names(expressed_genes_by_type)) {
  transcript_ids_of_interest <- expressed_genes_by_type[[cell_type]]
  ASPH_exons <- ASPH_annotation_from_gtf %>%
    filter(transcript_id %in% transcript_ids_of_interest, type == "exon") %>%
    mutate(cell_type = cell_type)
  if (nrow(ASPH_exons) > 0) {
    num_transcripts <- length(unique(ASPH_exons$transcript_id))
    heights <- c(heights, num_transcripts)
    
    ASPH_rescaled <- shorten_gaps(
      ASPH_exons,
      to_intron(ASPH_exons, "transcript_id"),
      group_var = "transcript_id"
    )
    plot <- ASPH_rescaled %>%
      filter(type == "exon") %>%
      ggplot(aes(xstart = start, xend = end, y = transcript_id)) +
      geom_range(aes(fill = cell_type), color = color_map[cell_type]) + 
      geom_intron(
        data = ASPH_rescaled %>% filter(type == "intron"),
        arrow.min.intron.length = 200
      ) +
      scale_fill_manual(values = color_map) + 
      labs(title = cell_type) +
      theme_minimal()
    ASPHplots[[cell_type]] <- plot
  } else {
    ASPHplots[[cell_type]] <- NULL
  }
}
#add reference
ASPHplots[["reference"]] <- ASPH_reference_plot
heights <- c(heights, 30)

ASPH_valid_plots <- Filter(Negate(is.null), ASPHplots)
ASPH_valid_heights <- heights[!sapply(ASPHplots, is.null)]
ASPH_rescaled_plot <- wrap_plots(ASPH_valid_plots, ncol = 1) + 
  plot_layout(heights = ASPH_valid_heights / sum(ASPH_valid_heights))
print(ASPH_rescaled_plot)

#Supplementary figure 8: plot EXOC1
EXOC1_filtered_cage <- pbid_fl %>% 
  filter(associated_gene == "EXOC1")
# filter your gtf for the gene of interest
EXOC1_PBid_of_interest <- EXOC1_filtered_cage$isoform
EXOC1_annotation_from_gtf <- pb_anno %>%
  dplyr::filter(
    !is.na(transcript_id),
    transcript_id %in% EXOC1_PBid_of_interest
  )
# extract the required annotation columns
EXOC1_annotation_from_gtf <- EXOC1_annotation_from_gtf %>% 
  dplyr::select(
    seqnames,
    start,
    end,
    strand,
    type,
    gene_id,
    transcript_id
  )
EXOC1_annotation_from_gtf %>% head()
#unique(EXOC1_annotation_from_gtf$transcript_id)
EXOC1_annotation_from_gtf$gene_id
EXOC1_annotation_from_gtf$gene_id <- "EXOC1"

#reference
gene_of_interest <- "EXOC1"
EXOC1_reference_annotation_from_gtf <- reference_anno %>%
  dplyr::filter(
    !is.na(gene_name),
    gene_name %in% gene_of_interest
  )
EXOC1_reference_annotation_from_gtf <- EXOC1_reference_annotation_from_gtf %>% 
  dplyr::select(
    seqnames,
    start,
    end,
    strand,
    type,
    gene_name,
    transcript_name
  )

#plot for each cell type
colors <- RColorBrewer::brewer.pal(12, "Paired")
color_map <- setNames(colors, c("Bipolar", "Rod", "Amacrine", "Muller_Glia", "Fibroblast", 
                                "Horizontal", "RPE", "Cone", "Astrocyte", "Smooth Muscle Cells", "zReference"))

temp_df <- data.frame()
EXOC1_exons_plot <- data.frame()
for (cell_type in names(expressed_genes_by_type)) {
  transcript_ids_of_interest <- expressed_genes_by_type[[cell_type]]
  temp_df <- EXOC1_annotation_from_gtf %>%
    filter(transcript_id %in% transcript_ids_of_interest, type == "exon") %>%
    mutate(cell_type = cell_type, transcript_to_plot = paste(cell_type, transcript_id, sep = "."))
  EXOC1_exons_plot <- rbind(EXOC1_exons_plot, temp_df)
}
EXOC1_reference_exons <- EXOC1_reference_annotation_from_gtf %>%
  filter(type == "exon")
EXOC1_reference_exons <- EXOC1_reference_exons %>% 
  mutate(transcript_id = transcript_name, gene_id = gene_name, cell_type = "zReference", transcript_to_plot = paste(cell_type, transcript_id, sep = "."))
EXOC1_reference_exons$gene_name <- NULL
EXOC1_reference_exons$transcript_name <- NULL
EXOC1_combined_plot <- rbind(EXOC1_exons_plot, EXOC1_reference_exons)
EXOC1_together_plot <- EXOC1_combined_plot %>%
  ggplot(aes(xstart = start, xend = end, y = transcript_to_plot)) +
  geom_range(aes(fill = cell_type), color = color_map[cell_type]) + 
  geom_intron(
    data = to_intron(EXOC1_combined_plot, "transcript_to_plot"),
    aes(strand = strand), 
    arrow.min.intron.length = 1000
  ) +
  scale_fill_manual(values = color_map) + 
  labs(title = "EXOC1") +
  theme_minimal()
print(EXOC1_together_plot)

EXOC1plots <- list()
heights <- numeric()
EXOC1_reference_rescaled <- shorten_gaps(
  EXOC1_reference_exons, 
  to_intron(EXOC1_reference_exons, "transcript_id"), 
  group_var = "transcript_id"
)
EXOC1_reference_plot <- EXOC1_reference_rescaled %>%
  dplyr::filter(type == "exon") %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
  geom_range() +
  geom_intron(
    data = EXOC1_reference_rescaled %>% dplyr::filter(type == "intron"), 
    arrow.min.intron.length = 200
  )
for (cell_type in names(expressed_genes_by_type)) {
  transcript_ids_of_interest <- expressed_genes_by_type[[cell_type]]
  EXOC1_exons <- EXOC1_annotation_from_gtf %>%
    filter(transcript_id %in% transcript_ids_of_interest, type == "exon") %>%
    mutate(cell_type = cell_type)
  if (nrow(EXOC1_exons) > 0) {
    num_transcripts <- length(unique(EXOC1_exons$transcript_id))
    heights <- c(heights, num_transcripts)
    
    EXOC1_rescaled <- shorten_gaps(
      EXOC1_exons,
      to_intron(EXOC1_exons, "transcript_id"),
      group_var = "transcript_id"
    )
    plot <- EXOC1_rescaled %>%
      filter(type == "exon") %>%
      ggplot(aes(xstart = start, xend = end, y = transcript_id)) +
      geom_range(aes(fill = cell_type), color = color_map[cell_type]) + 
      geom_intron(
        data = EXOC1_rescaled %>% filter(type == "intron"),
        arrow.min.intron.length = 200
      ) +
      scale_fill_manual(values = color_map) + 
      labs(title = cell_type) +
      theme_minimal()
    EXOC1plots[[cell_type]] <- plot
  } else {
    EXOC1plots[[cell_type]] <- NULL
  }
}
#add reference
EXOC1plots[["reference"]] <- EXOC1_reference_plot
heights <- c(heights, 8)

EXOC1_valid_plots <- Filter(Negate(is.null), EXOC1plots)
EXOC1_valid_heights <- heights[!sapply(EXOC1plots, is.null)]
EXOC1_rescaled_plot <- wrap_plots(EXOC1_valid_plots, ncol = 1) + 
  plot_layout(heights = EXOC1_valid_heights / sum(EXOC1_valid_heights))
print(EXOC1_rescaled_plot)

#Figure 4D: plot splice junction for ABCA4 variants
ABCA4_exons_2 <- ABCA4_annotation_from_gtf %>% dplyr::filter(type == "exon")
ABCA4_alphafold3_exons <- ABCA4_exons_2 %>% dplyr::filter(
  transcript_id %in% c("PB.5877.603", "PB.5877.497", "PB.5877.707", "PB.5877.738"))
ABCA4_alphafold3_exons <- ABCA4_alphafold3exons %>%
  mutate(
    celltype = case_when(
      transcript_id %in% c("PB.5877.603") ~ "Rod",
      transcript_id %in% c("PB.5877.497") ~ "Cone",
      transcript_id %in% c("PB.5877.707") ~ "Smooth Muscle Cells",
      transcript_id %in% c("PB.5877.738") ~ "RPE"
    ),
    transcript_to_plot = paste(celltype, transcript_id, sep = ".")
  )

ABCA4_reference_exons_2 <- ABCA4_reference_annotation_from_gtf %>% dplyr::filter(type == "exon")
ABCA4_201_exons <- ABCA4_reference_exons_2 %>% dplyr::filter(transcript_name == "ABCA4-201")
ABCA4_201_exons <- ABCA4_201_exons %>% rename(transcript_id = transcript_name)

ABCA4_alphafold3_junction <- ABCA4_junction %>% dplyr::filter(
  isoform %in% c("PB.5877.603", "PB.5877.497", "PB.5877.707", "PB.5877.738"))
ABCA4_alphafold3_junction <- ABCA4_alphafold3_junction %>% rename(transcript_id = isoform)
transcript_info <- ABCA4_alphafold3exons %>%
  select(transcript_id, celltype, transcript_to_plot) %>%
  distinct()
ABCA4_alphafold3junction <- ABCA4_alphafold3_junction %>%
  left_join(transcript_info, by = "transcript_id")
ABCA4_alphafold3junction <- ABCA4_alphafold3_junction %>% rename(start = genomic_start_coord)
ABCA4_alphafold3junction <- ABCA4_alphafold3_junction %>% rename(end = genomic_end_coord)

ABCA4_alphafold3exons %>%
  ggplot(aes(xstart = start, xend = end, y = transcript_to_plot)) +
  geom_range(aes(fill = celltype), color = color_map[celltype]) + 
  geom_intron(data = to_intron(ABCA4_alphafold3exons, "transcript_to_plot")) + 
  geom_junction(data = ABCA4_alphafold3junction, junction.y.max = 0.5) + 
  geom_junction_label_repel(data = ABCA4_alphafold3junction,aes(label = junction_category), junction.y.max = 0.5) +
  scale_fill_manual(values = color_map) + 
  labs(title = celltype)

#Figure 1-H, supplementary figure 2-4:
#Pigeon report was generated using the scripts provided on SQANTI3 GitHub (https://github.com/ConesaLab/SQANTI3/tree/master/utilities/report_pigeon)

####end of the session########################
####Author: Luozixian Wang, Raymond Wong######
####20/06/2024################################




