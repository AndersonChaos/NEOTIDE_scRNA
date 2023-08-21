library(tidyverse)
library(ggpubr)
library(readxl)
library(SeuratDisk)
library(rstatix)
library(Seurat)

# Data Preprocess
IIT_id <- c("P343", "P438", "P519", "P529", "P531", "P547", "P549", "P586", "P589", "P590", "P591")
Epi <- readRDS("ALL_LUAD/Epi.rds")
Epi <- subset(Epi, subset = sampleID %in% IIT_id)

doublet.cell.IDs <- unique(
  c(
    WhichCells(Epi, expression = PTPRC > 0 & EPCAM > 0)
  )
)
length(doublet.cell.IDs)
Epi <- Epi[, !Epi$cellID %in% doublet.cell.IDs]
sampleIDs <- unique(Epi$sampleID)
copykat <- data.frame()
for (id in sampleIDs) {
  print(id)
  if (file.exists(paste0("copykat_ALL_LUAD/", id, "_copykat_prediction.txt")) == FALSE) next
  copy_res <- read_table(
    paste0("copykat_ALL_LUAD/", id, "_copykat_prediction.txt"),
    col_types = cols(.default = col_character())
  )
  copy_res$cell.names <- str_replace(copy_res$cell.names, "\\.", "-")
  copykat <- bind_rows(copykat, copy_res)
}
meta <- Epi@meta.data %>% dplyr::select(-starts_with("copykat"))
meta <- meta %>% left_join(CNA_copykat, by = c("cellID" = "cellID"))
rownames(meta) <- meta$cellID
head(meta)
meta <- meta %>%
  mutate(
    copykat.pred = ifelse(is.na(copykat.pred), "not.defined", copykat.pred)
  )
Epi@meta.data <- meta
Epi <- subset(Epi, subset = sampleID %in% IIT_id)
all.genes <- rownames(Epi)
Epi <- Epi %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData(feature = all.genes)

Epi <- Epi %>%
  RunPCA(features = VariableFeatures(object = Epi)) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.6) %>%
  RunUMAP(reduction = "pca", dims = 1:20)


DimPlot(Epi, reduction = "umap", group.by = "sampleID", label = TRUE, repel = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(Epi, reduction = "umap", group.by = "cell.type", pt.size = 0.5)
ggsave("copykat_pred.pdf", width = 10, height = 10)
DimPlot(Epi, pt.size = 0.5, label = TRUE, label.size = 10) + NoLegend()

saveRDS(Epi, "ALL_LUAD/Epi.rds")
Epi_ID <- c("Normal_E1_SCGB3A2", "Normal_E2_SFTPC", "Tumor_E3_CAPS", "Tumor_E4_CDK4", "Tumor_E5_CEACAM6", "Normal_E6_RGS1", "Normal_E7_KRT5", "Normal_E8_MUC5B", "Tumor_E9_TMSB10", "Tumor_E10_CXCL14", "Normal_E11_AGER", "Tumor_E12_HCN1")

Idents(Epi) <- Epi$seurat_clusters
names(Epi_ID) <- levels(Epi)
Epi <- RenameIdents(Epi, Epi_ID)
Epi$cell.type <- Idents(Epi)


#########################################
#                                       #
#            Tumor Epithelial           #
#                                       #
#########################################
Tumor_Epi <- Epi %>% subset(subset = seurat_clusters %in% c(2, 3, 4, 8, 9, 11))
Tumor_Epi <- subset(Tumor_Epi, subset = sampleID %in% IIT_id)

Idents(Tumor_Epi) <- Tumor_Epi$Treatment
Tumor_Epi_meta <- Tumor_Epi@meta.data
saved_id <- Tumor_Epi_meta %>% count(sampleID) %>% filter(n >= 10) %>% pull(sampleID)
Tumor_Epi <- subset(Tumor_Epi, subset = sampleID %in% saved_id)
Tumor_Epi <- Tumor_Epi %>% subset(sampleID %in% IIT_id)
Tumor_Epi_meta <- Tumor_Epi@meta.data
Tumor_Epi_meta <- Tumor_Epi_meta %>%
  mutate(
    NOD_ID = case_when(
      sampleID == "P343" ~ "NOD01",
      sampleID == "P438" ~ "NOD02",
      sampleID == "P589" ~ "NOD03",
      sampleID == "P590" ~ "NOD04",
      sampleID == "P591" ~ "NOD05",
      sampleID == "P529" ~ "NOD06",
      sampleID == "P531" ~ "NOD07",
      sampleID == "P519" ~ "NOD08",
      sampleID == "P547" ~ "NOD09",
      sampleID == "P549" ~ "NOD10",
      sampleID == "P586" ~ "NOD11",
      .default = "Real_World"
    ),
    response = case_when(
      NOD_ID %in% c("NOD02", "NOD05", "NOD09") ~ "Highly Resistant",
      NOD_ID %in% c("NOD06", "NOD01", "NOD08", "NOD07") ~ "Moderate response",
      NOD_ID %in% c("NOD10", "NOD11", "NOD03", "NOD04") ~ "Immune Sensitive",
    )
  )

all.genes <- rownames(Tumor_Epi)
Tumor_Epi <- Tumor_Epi %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData(feature = all.genes)

Tumor_Epi <- Tumor_Epi %>%
  RunPCA(features = VariableFeatures(object = Tumor_Epi)) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.4) %>%
  RunUMAP(reduction = "pca", dims = 1:20)

DimPlot(Tumor_Epi, pt.size = 0.5, label = TRUE, label.size = 5) + NoLegend()
DimPlot(Tumor_Epi, pt.size = 0.5, label = TRUE, label.size = 5, group.by = "copykat.pred") + NoLegend()
DimPlot(Tumor_Epi, pt.size = 0.5, label = TRUE, label.size = 5, group.by = "sampleID") + NoLegend()
Idents(Tumor_Epi) <- Tumor_Epi$response
saveRDS(Tumor_Epi, "ALL_LUAD/Tumor_Epi.rds")
Tumor_Cluster_Marker <- FindAllMarkers(Tumor_Epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
  group_by(cluster)
  #  %>% top_n(n = 100, wt = avg_log2FC)
View(Tumor_Cluster_Marker)
write_csv(Tumor_Cluster_Marker, "Tumor_Cluster_Marker.csv")
ADC_list <- c("ERBB2", "TACSTD2", "NECTIN4", "FOLR1")
Checkpoint <- c("CD274", "NT5E", "CD47", "PDCD1LG2", "LAG3", "BTLA")
TAA <- c("CEACAM5", "CEACAM6")
VlnPlot(Tumor_Epi_IIT, features = c(ADC_list, TAA, Checkpoint), pt.size = 0, group.by = "response")
VlnPlot(Tumor_Epi_IIT, features = c(ADC_list, TAA, Checkpoint), pt.size = 0)
Tumor_Exp_selected <- Tumor_Epi_IIT[["RNA"]]@data[c(TAA, ADC_list, Checkpoint), ]

Tumor_Exp_selected <- Tumor_Exp_selected %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "cellID", values_to = "expression") %>%
  mutate(
    sampleID = str_split(cellID, "_", simplify = TRUE)[, 1]
  )
Tumor_Exp_selected <- Tumor_Exp_selected %>%
  left_join(Sample)
write_csv(Tumor_Exp_selected, "cellratio/Tumor_Exp_selected.csv")



###############################
#                             #
#         CNV score           #
#                             #
###############################

sampleIDs <- unique(Epi$sampleID)
copykat_CNA_res <- data.frame()
for (id in sampleIDs) {
  print(id)
  if (file.exists(paste0("copykat_ALL_LUAD/", id, "_copykat_CNA_results.txt")) == FALSE) next
  copy_res <- read_table(
    paste0("copykat_ALL_LUAD/", id, "_copykat_CNA_results.txt")
  )
  copy_res <- pivot_longer(copy_res, cols = -c("chrom", "chrompos", "abspos"), names_to = "cellID", values_to = "copy_num")
  copy_res$cellID <- str_replace(copy_res$cellID, "\\.", "-")
  copykat_CNA_res <- bind_rows(copykat_CNA_res, copy_res)
}

CNA_score <- copykat_CNA_res %>% group_by(cellID) %>% summarise(sum = sum(abs(copy_num)))
CNA_copykat <- CNA_score %>% left_join(copykat, by = c("cellID" = "cell.names"))
CNA_copykat
ggboxplot(CNA_copykat, x = "copykat.pred", y = "sum")
head(meta)

ggboxplot(meta, x = "seurat_clusters", y = "sum")

# Epi UMAP data save

UMAP_data <- Epi@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column(var = "cellID") %>%
  as_tibble()
cell_meta <- Epi@meta.data %>% as_tibble()
UMAP_data <- left_join(UMAP_data, cell_meta, by = "cellID") %>%
  dplyr::rename(CNVscore = sum) %>%
  left_join(Sample, by = "sampleID")
write_csv(UMAP_data, "cellratio/Epi_UMAP.csv")

########################################
#                                      #
#        CellphoneDB preparing         #
#                                      #
########################################

Myeloid_part_IIT <- subset(Myeloid_part, subset = sampleID %in% IIT_id)
cellphoneDB_data <- merge(Myeloid_part_IIT, Tumor_Epi_IIT)
table(cellphoneDB_data$cell.type, cellphoneDB_data$response)
SaveH5Seurat(cellphoneDB_data, filename = "ALL_LUAD/Myeloid_Epi.h5seurat", overwrite = TRUE)
Convert("ALL_LUAD/Myeloid_Epi.h5seurat", dest = "h5ad", overwrite = TRUE)
meta_epi <- cellphoneDB_data@meta.data[, c("cellID", "cell.type")]
write_csv(meta_epi, "ALL_LUAD/Myeloid_Epi_meta.csv")

# UMAP plot
set3_color <- RColorBrewer::brewer.pal(12, "Set3")
set3_color <- c(set3_color, "#CCAC69", "#F1C3A2", "#E0C68C")
set3_color[9] <- "#D6D4BF"
set3_color[2] <- "#F7CE5C"
DimPlot(Epi, pt.size = 0.5, cols = set3_color) + theme(
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1, linetype = "solid"),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()) + labs(x = "UMAP1", y = "UMAP2", title = "Epithelial Cells")
ggsave("cellratio/Epithelial.pdf", width = 8, height = 6)