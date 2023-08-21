library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(SeuratDisk)
library(rstatix)
library(Seurat)
library(harmony)
library(readxl)
library(ggsci)
library(msigdbr)
# Read in the data

Myeloid_part <- readRDS("ALL_LUAD/Myeloid_part.rds")
# Preprocess the data and run dimension reduction
all.genes <- rownames(Myeloid_part)
Myeloid_part <- Myeloid_part %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData(feature = all.genes)

Myeloid_part <- Myeloid_part %>%
  RunPCA(features = VariableFeatures(object = Myeloid_part)) %>%
  RunHarmony(group.by.vars = "sampleID") %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.8) %>%
  RunUMAP(reduction = "harmony", dims = 1:20)

# Identify cell types
Myeloid_canonical <- list(
  "M" = c("CD68", "LYZ"),
  "cDC1" = c("CLEC9A", "FLT3", "IDO1"),
  "cDC2" = c("CD1C", "FCER1A", "HLA-DQA1"),
  "cDC3_LAMP3" = c("LAMP3", "CCR7", "FSCN1"),
  "pDC" = c("LILRA4", "GZMB", "IL3RA"),
  "Mφ_S100A8" = c("FCN1", "S100A9", "S100A8", "CD14", "VCAN"),
  "Mφ_MKI67" = c("MKI67", "STMN1"),
  "Mφ_CD16" = c("FCGR3A", "LST1", "LILRB2"),
  "Mφ_APOE" = c("APOE", "C1QC", "C1QA", "C1QB", "MARCO", "SPP1", "F13A1"),
  "Neu" = c("CSF3R", "RGS2", "FCGR3B"),
  "Mast" = c("KIT", "TPSAB1", "CPA3"),
  "PDL1" = c("CD274", "PDCD1LG2")
)
# Find markers for each cluster
M_marker <- FindAllMarkers(Myeloid_part, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC)
View(M_marker)

# Indentify cell types and assign them to the object
M_cluster <- c(
  "Mono_S100A9", "Mast cell", "cDC2", "Macro_C1QC",
  "Macro_C1QC", "Macro_EREG", "Macro_FBP1", "cDC2", "Macro_MKI67", "Mono_S100A9", "Mono_CD16", "Macro_APOE", "Macro_CXCL9", "Neutrophil", "pDC_LILRA4", "cDC1_CLEC9A", "Macro_CTSK", "mDC_LAMP3", "doublet", "Macro_C1QC"
)
Idents(Myeloid_part) <- Myeloid_part$seurat_clusters
names(M_cluster) <- levels(Myeloid_part)
Myeloid_part <- RenameIdents(Myeloid_part, M_cluster)
Myeloid_part$cell.type <- Idents(Myeloid_part)
Idents(Myeloid_part) <- factor(Idents(Myeloid_part), levels = c("Macro_C1QC", "Macro_APOE", "Macro_CXCL9", "Macro_EREG", "Macro_FBP1", "Macro_CTSK", "Macro_MKI67", "Mono_S100A9", "Mono_CD16", "mDC_LAMP3", "Neutrophil", "cDC1_CLEC9A", "cDC2", "pDC_LILRA4", "Mast cell"))

VlnPlot(Myeloid_part, features = c("CD3D", "CD19", "EPCAM", "PECAM1", "COL1A1", "NKG7", "CD68", "JCHAIN", "CSF3R", "CLEC9A"), pt.size = 0)


# Dotplot
DotPlot(Myeloid_part, features = c("C1QC", "CCL4", "C1QA", "APOE", "GPNMB", "APOC1", "CXCL9", "CXCL10", "GBP1", "EREG", "OLR1", "THBD", "FBP1", "MARCO", "FN1", "CTSK", "SPARC", "MMP9", "MKI67", "STMN1", "TUBB", "S100A9", "VCAN", "FCN1", "FCGR3A", "CDKN1C", "ZNF703", "CCL22", "LAMP3", "FSCN1", "FCGR3B", "CMTM2", "CSF3R", "CLEC9A", "IRF8", "CPNE3", "CD1C", "CD1A", "HLA-DQB2", "LILRA4", "GZMB", "IL3RA", "KIT", "TPSAB1", "CPA3"), scale.min=10, scale.max = 70, col.max = 1.5, col.min = -1.5) + RotatedAxis() + scale_colour_gradientn(colours = c("#052550", "#3984BB", "#95C6DF", "white", "#F3AA89", "red", "#a30018"))
set3_color <- RColorBrewer::brewer.pal(12, "Set3")
set3_color <- c(set3_color, "#CCAC69", "#F1C3A2", "#E0C68C")
set3_color[9] <- "#D6D4BF"
set3_color[2] <- "#F7CE5C"
DimPlot(Myeloid_part, pt.size = 0.5, cols = set3_color) + theme(
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1, linetype = "solid"),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()) + labs(x = "UMAP1", y = "UMAP2", title = "Myeloid Cells")
ggsave("cellratio/Epithelial.pdf", width = 8, height = 6)

saveRDS(Myeloid_part, "ALL_LUAD/Myeloid_part.rds")

cellratio <- prop.table(table(Myeloid_part$sampleID, Myeloid_part$cell.type), margin = 1) # calculate ratio of sub celltype
meta <- as.data.frame(cellratio)
meta <- meta %>% left_join(Sample, by = c("Var1" = "sampleID")) # merge with sample information
meta <- meta %>%
  rename(
    sampleID = Var1,
    Ratio = Freq,
    Celltype = Var2
  ) %>%
  mutate(
    Smoking_History = str_extract(Smoking_History, "(Y|N)"),
    EGFR = ifelse(str_detect(EGFR, "L858R"), "L858R",
      ifelse(str_detect(EGFR, "exon19del"), "exon19del", ifelse(str_detect(EGFR, "exon20ins"), "exon20ins", "other"))
    ),
    PRR = as.double(`Pathological Response Rate`),
    PRR = ifelse(is.na(PRR), 0, PRR),
    response = ifelse(PRR > median(PRR), "Responder", "Non-Responder"),
    Pre_treatment_PDL1 =
      ifelse(`PD-L1 TPS` == "Not tested", NA,
      ifelse(`PD-L1 TPS` == "<1%", "<1%",
      ifelse(`PD-L1 TPS` == "0%", "<1%", ">=1%")
        )
      )

  )  %>%
  filter(!str_detect(Histology, "SCLC") & !is.na(Pre_treatment_PDL1)) %>%
    mutate(MyGroup = ifelse(Pre_treatment_PDL1 == "<1%", ifelse(response == "Responder", "<1% & Responder", "<1% & Non-Responder"), ifelse(response == "Responder", ">=1% & Responder", ">=1% & Non-Responder")))

