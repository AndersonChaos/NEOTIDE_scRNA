library(tidyverse)
library(SeuratDisk)
library(rstatix)
library(Seurat)
library(harmony)
library(readxl)
# Read in the data
B_part <- readRDS("ALL_LUAD/B_part.rds")

# Preprocess the data and run dimension reduction

all.genes <- rownames(B_part)
B_part <- B_part %>%
    NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
    FindVariableFeatures(selection.method = "vst") %>%
    ScaleData(feature = all.genes)
B_part <- B_part %>%
    RunPCA(features = VariableFeatures(object = B_part)) %>%
    RunHarmony(group.by.vars = "sampleID") %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = 0.6) %>%
    RunUMAP(reduction = "harmony", dims = 1:20)


# Identify cell types
B_canonical <- list(
    B = c("CD79A", "CD79B", "CD19", "MS4A1"),
    `Naive B` = c("IGHM", "IGHD", "FCER2", "IL4R"), `Memory B` = c("CD27", "TNFRSF13B", "CD24"),
    `Germinal center B` = c("BCL6", "AICDA", "CXCR5"
    , "STMN1"),
    Plasma = c("XBP1", "MZB1", "IGHG1", "IGHA1", "JCHAIN")
)

B_marker <- FindAllMarkers(B_part, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
    group_by(cluster) %>%
    top_n(10, avg_log2FC) %>%
    ungroup()
View(B_marker)
# Vlnplot of canonical markers
VlnPlot(B_part, features = c("CD79A", "IGHM", "CD3D", "IGHD", "CD27", "IL4R", "BCL6", "PDE4D", "XBP1", "JCHAIN", "IGHA1", "IGHG1", "STMN1", "TNFSF9", "SELL", "DUSP4"), pt.size = 0)
VlnPlot(B_part, features = c("EPCAM", "COL1A1", "CD3D", "NKG7", "CD68", "CD19", "CD79A", "JCHAIN", "LYZ"), pt.size = 0)
# Rename the clusters
B_cluster <- c("B_memory_TNFSF9", "B_memory_PDE4D", "B_memory_LTB", "B_memory_EGR1", "B_memory_DUSP4", "B_naive", "Plasma", "Plasma", "B_memory_MAP4K4", "Plasma", "B_memory_TNFSF9", "B_memory_DUSP4", "GCB", "doublet", "Plasma")
Idents(B_part) <- B_part$seurat_clusters
names(B_cluster) <- levels(B_part)
B_part <- RenameIdents(B_part, B_cluster)
B_part$cell.type <- Idents(B_part)
saveRDS(B_part, "ALL_LUAD/B_part.rds")

# Cellratio plot

cellratio <- prop.table(table(B_part$sampleID, B_part$cell.type), margin = 1)  # calculate ratio of sub celltype
meta <- as.data.frame(cellratio)
write_csv(meta, "cellratio/B.csv")
meta <- meta %>%
    left_join(Sample, by = c(Var1 = "sampleID"))  # merge with sample information
meta <- meta %>% # Rename
  rename(
    sampleID = Var1,
    Ratio = Freq,
    Celltype = Var2
  )
meta <- meta %>%
  mutate(
    major_type = ifelse(str_detect(Celltype, "CD4"), "CD4", ifelse(str_detect(Celltype, "CD8"), "CD8", "NK")),
    Smoking_History = str_extract(Smoking_History, "(Y|N)"),
    EGFR = ifelse(str_detect(EGFR, "L858R"), "L858R",
      ifelse(str_detect(EGFR, "exon19del"), "exon19del", ifelse(str_detect(EGFR, "exon20ins"), "exon20ins", "other"))
    ),
    PRR = as.double(`Pathological Response Rate`),
    PRR = ifelse(is.na(PRR), 0, PRR),
    response = ifelse(PRR > median(PRR), "Responder", "Non-Responder"),
    isMPR = factor(`Pathological Response`, levels = c("non-MPR", "MPR")),
    Pre_treatment_PDL1 =
      ifelse(`PD-L1 TPS` == "Not tested", NA,
        ifelse(`PD-L1 TPS` == "<1%", "<1%",
          ifelse(`PD-L1 TPS` == "0%", "<1%", ">=1%")
        )
      ),
    PDL1_MPR = case_when(
      Pre_treatment_PDL1 == "<1%" & isMPR == "MPR" ~ "<1% & MPR",
      Pre_treatment_PDL1 == "<1%" & isMPR == "non-MPR" ~ "<1% & non-MPR",
      Pre_treatment_PDL1 == ">=1%" & isMPR == "MPR" ~ ">=1% & MPR",
      Pre_treatment_PDL1 == ">=1%" & isMPR == "non-MPR" ~ ">=1% & non-MPR"
    ),
    Smoke_MPR = case_when(
      Smoking_History == "Y" & isMPR == "MPR" ~ "ES & MPR",
      Smoking_History == "Y" & isMPR == "non-MPR" ~ "ES & non-MPR",
      Smoking_History == "N" & isMPR == "MPR" ~ "NS & MPR",
      Smoking_History == "N" & isMPR == "non-MPR" ~ "NS & non-MPR"
    )
  ) %>%
  filter(!str_detect(Histology, "SCLC"))
ggboxplot(
    meta,
    x = "Celltype", y = "Ratio",
    fill = "Smoke_MPR", palette = "jco") +
    stat_compare_means(aes(group = Smoke_MPR), method = "wilcox.test", label = "p.signif", size = 3, hide.ns = TRUE) +
theme_classic() +
theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

# UMAP plot
set3_color <- RColorBrewer::brewer.pal(9, "Set3")
set3_color <- c(set3_color, "#FED9A6")
set3_color[8] <- "#FED9A6"
set3_color[9] <- "#E78AC3"

Idents(scRNA) <- scRNA$major_cell_type
scRNA_IIT <- scRNA %>% subset(sampleID %in% IIT_id)
scRNA_RWC <- scRNA %>% subset(sampleID %in% IIT_id, invert = TRUE)
DimPlot(B_part, pt.size = 0.01, cols = set3_color) + theme(
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1, linetype = "solid"),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()) + labs(x = "UMAP1", y = "UMAP2", title = "Major_cell_type")
ggsave("cellratio/B_celltype_UMAP.pdf")

# Dotplot
DotPlot(B_part, features = c("TNFSF9", "ZNF331", "NR4A2", "PDE4D", "HOPX", "CRIP1", "LTB", "CD52", "MS4A1", "TNF", "EGR1", "CITED2", "DUSP4", "RGS1", "IFITM1", "MALAT1", "MAP4K4", "ARIH1", "TCL1A", "FCER2", "YBX3", "IGHG1", "MZB1", "XBP1",  "RGS13", "NEIL1", "LMO2"), scale.min=10, scale.max = 70, col.max = 1.5, col.min = -1.5, cols = c("lightgrey", "red")) + RotatedAxis() + scale_colour_gradientn(colours = c("#052550", "#3984BB", "#95C6DF", "white", "#F3AA89", "red", "#a30018"))
ggsave("cellratio/B_dotplot.pdf")