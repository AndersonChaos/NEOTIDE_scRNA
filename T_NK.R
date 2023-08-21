library(tidyverse)
library(ggsci)
library(ggpubr)
library(SeuratDisk)
library(rstatix)
library(Seurat)
library(harmony)
library(readxl)
library(GSVA)
library(msigdbr)
Sample <- read_excel("Sample.xlsx") %>%
  mutate(sampleID = str_extract(SampleID, "P[0-9]+")) %>%
  dplyr::select(
    sampleID, Gender, Age, Histology, Cycles, PD1,
    `Pathological Response`, EGFR, `Pathological Response Rate`,
    Smoking_History, `PD-L1 TPS`
  )

# Read in data, preprocess and run dimension reduction

T_part <- readRDS("ALL_LUAD/T_part.rds")

all.genes <- rownames(T_part)
T_part <- T_part %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData(feature = all.genes)

T_part <- T_part %>%
    RunPCA(features = VariableFeatures(object = T_part)) %>%
    RunHarmony(group.by.vars = "sampleID") %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = 1) %>%
    RunUMAP(reduction = "harmony", dims = 1:20)

T_canonical <- list(
  "Lineage marker" = c("CD3D", "CD8A", "CD4", "NKG7"),
  "Effector marker" = c("GZMK", "GZMB", "PRF1", "FGFBP2"),
  "Exhaustion marker" = c("HAVCR2", "LAG3", "PDCD1", "TIGIT", "CTLA4", "LAYN", "ENTPD1", "IFNG"),
  "Memory marker" = c("CCR7", "TCF7", "IL7R", "SELL"),
  "Proliferation marker" = c("MKI67", "STMN1"),
  "Treg marker" = c("FOXP3", "IL2RA")
)

DotPlot(T_part, features = T_canonical, cols = c("lightgrey", "red")) + RotatedAxis()
VlnPlot(T_part, features = c("CD4", "CD8A", "NKG7", "CD19", "CD3D", "CD68", "EPCAM", "PECAM1", "COL1A1", "IGHG1"), pt.size = 0)

FeaturePlot(T_part, features = c("CD4", "CD8A", "NKG7", "CD19", "CD3D", "CD68", "EPCAM", "PECAM1", "COL1A1"), pt.size = 0)

DimPlot(T_part,pt.size = 0.001, raster = FALSE) + theme(
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1, linetype = "solid"),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()) + labs(x = "UMAP1", y = "UMAP2", title = "T/NK")
ggsave("cellratio/T_part_UMAP.pdf", width = 8, height = 6)
DimPlot(T_part, label = TRUE, label.size = 4, group.by = "cell.type")



# Classify T/NK into CD4, CD8 and NK for further analysis

CD4 <- subset(T_part, subset = seurat_clusters %in% c(0, 2, 6, 8, 10, 13, 15))
CD8 <- subset(T_part, subset = seurat_clusters %in% c(1, 3, 4, 5, 7, 14, 17))
NK <- subset(T_part, subset = seurat_clusters %in% c(9, 11, 18))

saveRDS(T_part, "ALL_LUAD/T_part.rds")
saveRDS(CD4, "ALL_LUAD/CD4.rds")
saveRDS(CD8, "ALL_LUAD/CD8.rds")
saveRDS(NK, "ALL_LUAD/NK.rds")


# T_part cell ratio barplot
cellratio <- prop.table(table(T_part$sampleID, T_part$cell.type), margin = 1) # calculate ratio of sub celltype
meta <- as.data.frame(cellratio)
meta <- meta %>% left_join(Sample, by = c("Var1" = "sampleID")) # merge with sample information
meta <- meta %>% # Rename
  dplyr::rename(
    sampleID = Var1,
    Ratio = Freq,
    Celltype = Var2
  )
meta <- meta %>%
  mutate(
    Smoking_History = str_extract(Smoking_History, "(Y|N)"),
    EGFR = ifelse(str_detect(EGFR, "L858R"), "L858R",
      ifelse(str_detect(EGFR, "exon19del"), "exon19del", ifelse(str_detect(EGFR, "exon20ins"), "exon20ins", "other"))
    ),
    PRR = as.double(`Pathological Response Rate`),
    PRR = ifelse(is.na(PRR), 0, PRR),
    response = ifelse(PRR > median(PRR), "Responder", "Non-Responder"),
    isMPR = ifelse(`Pathological Response` %in% c("MPR", "pCR"), "MPR", "non-MPR"),
    isMPR = factor(isMPR, levels = c("non-MPR", "MPR")),
    data_from = ifelse(sampleID %in% IIT_id, "2104", "Real_World"),
    group_MPR = case_when(
      isMPR == "MPR" & data_from == "2104" ~ "MPR_2104",
      isMPR == "non-MPR" & data_from == "2104" ~ "non-MPR_2104",
      isMPR == "MPR" & data_from == "Real_World" ~ "MPR_Real_World",
      isMPR == "non-MPR" & data_from == "Real_World" ~ "non-MPR_Real_World"
    ),
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
    ),
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
  ) %>%
  filter(!str_detect(Histology, "SCLC"))
ggboxplot(
    meta,
    x = "Celltype", y = "Ratio",
    color = "response", palette = "jco", add = "jitter") +
    stat_compare_means(aes(group = response), method = "wilcox.test", label = "p.signif", size = 3, hide.ns = TRUE) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

ggboxplot(
    meta %>% filter(!is.na(Resistant)),
    x = "Celltype", y = "Ratio",
    color = "Highly_Resistant", palette = "jco", add = "jitter") +
    stat_compare_means(aes(group = Highly_Resistant), method = "wilcox.test", label = "p.signif", size = 3, hide.ns = TRUE) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
#########################################
#                                       #
#           CD4 part analysis           #
#                                       #
#########################################
CD4 <- readRDS("ALL_LUAD/CD4.rds")
all.genes <- rownames(CD4)
CD4 <- CD4 %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData(feature = all.genes)

CD4 <- CD4 %>%
    RunPCA(features = VariableFeatures(object = CD4)) %>%
    RunHarmony(group.by.vars = "sampleID") %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = 0.8) %>%
    RunUMAP(reduction = "harmony", dims = 1:20)


# Calculate CD4 Markers
CD4_marker <- FindAllMarkers(CD4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC) %>%
  ungroup()
View(CD4_marker)

# set CD4 idents names
CD4_id <- c("CD4T_Tm_ANXA1", "CD4T_Tm_CCL5", "CD4T_Treg_FOXP3", "CD4T_Th1-like_CXCL13", "CD4T_Tfh_CXCL13", "CD4T_Tn_CCR7", "CD4T_Treg_CCR8", "CD4T_Treg_FOXP3", "CD8T_Tem_CCL5", "CD8T_MAIT_KLRB1", "CD4T_Treg_FOXP3", "CD4T_Treg_CCR8", "CD8T_Tem_GZMK")

Idents(CD4) <- CD4$seurat_clusters
names(CD4_id) <- levels(CD4)
CD4 <- RenameIdents(CD4, CD4_id)
CD4$cell.type <- Idents(CD4)

saveRDS(CD4, "ALL_LUAD/CD4.rds")

# Dotplot and UMAP plot
Idents(CD4) <- factor(Idents(CD4), levels = c("CD4T_Treg_FOXP3", "CD4T_Treg_CCR8", "CD4T_Tm_ANXA1", "CD4T_Tm_CCL5", "CD4T_Th1-like_CXCL13", "CD4T_Tfh_CXCL13", "CD4T_Tn_CCR7", "CD8T_Tem_CCL5", "CD8T_MAIT_KLRB1", "CD8T_Tem_GZMK"))
DotPlot(CD4, features = CD4_feature, cols = c("lightgrey", "red"), scale.min = 10, scale.max = 70, col.min = -1.5, col.max = 1.5) + RotatedAxis() + scale_colour_gradientn(colours = c("#052550", "#3984BB", "#95C6DF", "white", "#F3AA89", "red", "#a30018"))
set3_color <- RColorBrewer::brewer.pal(12, "Set3")
set3_color <- c(set3_color, "#CCAC69", "#F1C3A2", "#E0C68C")
set3_color[9] <- "#b1b09f"
set3_color[2] <- "#F7CE5C"
DimPlot(CD4, pt.size = 0.01, cols = set3_color) + theme(
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1, linetype = "solid"),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()) + labs(x = "UMAP1", y = "UMAP2", title = "CD4_sub_cell_type")

#########################################
#                                       #
#           CD8 part analysis           #
#                                       #
#########################################

CD8 <- readRDS("ALL_LUAD/CD8.rds")
all.genes <- rownames(CD8)
CD8 <- CD8 %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData(feature = all.genes)

CD8 <- CD8 %>%
    RunPCA(features = VariableFeatures(object = CD8)) %>%
    RunHarmony(group.by.vars = "sampleID") %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = 0.8) %>%
    RunUMAP(reduction = "harmony", dims = 1:20)
# Plots
DimPlot(CD8, label = TRUE, pt.size = 0.5, label.size = 8)
VlnPlot(CD8, features = c("CD4", "CD8A", "NKG7", "CD68", "EPCAM", "LYZ", "CD79A", "CD3D", "PECAM1", "COL1A1", "JCHAIN"), pt.size = 0, group.by = 'seurat_clusters')
VlnPlot(CD8, features = c("CD4", "NEAT1", "KLRB1", "ITGAE"), pt.size = 0)
VlnPlot(CD8, features = c("FOXP3", "IL7R", "CTLA4", "PDCD1", "GZMK", "GZMB", "ITGAE"), pt.size = 0, )
DotPlot(CD8, features = c("HAVCR2", "IL7R", "CTLA4", "PDCD1", "SLC4A10", "TNFRSF9", "NKG7",  "SELL", "FGFBP2", "ISG15", "STMN1", "MKI67", "CCR7", "CD8A", "CD69", "LAYN", "ENTPD1", "ZNF683", "CXCR6", "XCL1", "ITGAE")) + RotatedAxis()

# Calculate CD8 Markers
CD8_Marker <- FindAllMarkers(CD8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC) %>%
  ungroup()
View(CD8_Marker)

# rename CD8 idents
CD8_id <- c("CD8T_Tem_ANXA1", "CD8T_Tm_NEAT1", "CD8T_Tex_CXCL13", "CD8T_Tem_TNFSF9", "CD8T_Tem_TNFSF9", "CD8T_Tem_GZMK", "CD8T_Trm_ZNF683", "CD8T_Tem_GZMK", "CD4T_Tm_NEAT1", "T_prolifer_STMN1", "T_prolifer_STMN1", "CD8T_ISG15")
Idents(CD8) <- CD8$seurat_clusters
names(CD8_id) <- levels(CD8)
CD8 <- RenameIdents(CD8, CD8_id)
CD8$cell.type <- Idents(CD8)
saveRDS(CD8, "ALL_LUAD/CD8.rds")



Idents(CD8) <- factor(Idents(CD8), levels = c("CD8T_Tem_ANXA1", "CD8T_Tem_GZMK", "CD8T_Tem_TNFSF9", "CD8T_Trm_ZNF683", "CD8T_Tex_CXCL13", "CD8T_Tm_NEAT1", "CD8T_ISG15", "CD4T_Tm_NEAT1", "T_prolifer_STMN1"))
DotPlot(CD8, features = CD8_feature, cols = c("lightgrey", "red"), scale.min = 10, scale.max = 70, col.min = -1.5, col.max = 1.5) + RotatedAxis() + scale_colour_gradientn(colours = c("#052550", "#3984BB", "#95C6DF", "white", "#F3AA89", "red", "#a30018"))
set3_color <- RColorBrewer::brewer.pal(12, "Set2")
set3_color <- c(set3_color, "#9176ff", "#F1C3A2", "#E0C68C")
set3_color[9] <- "#9176ff"
set3_color[2] <- "#F7CE5C"

DimPlot(CD8, pt.size = 0.01, cols = set3_color) + theme(
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1, linetype = "solid"),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()) + labs(x = "UMAP1", y = "UMAP2", title = "CD8_sub_cell_type")
#########################################
#                                       #
#            NK part analysis           #
#                                       #
#########################################


# Read in data, preprocess and run dimension reduction
NK <- readRDS("ALL_LUAD/NK.rds")
all.genes <- rownames(NK)
NK <- NK %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData(feature = all.genes)
NK <- NK %>%
    RunPCA(features = VariableFeatures(object = NK)) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = 1) %>%
    RunUMAP(reduction = "harmony", dims = 1:20)

# NK markers
NK_marker <- FindAllMarkers(NK, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC) %>%
  ungroup()
View(NK_marker)



# Rename NK idents
NK_id <- c("NK_XCL1", "CD8T_Teff_FGFBP2", "NK_FGFBP2", "CD8T_gdT_TRDV1", "NK_XCL1", "ILC3", "CD8T_Teff_FGFBP2", "NK_FGFBP2", "NK_XCL1", "ILC3")
Idents(NK) <- NK$seurat_clusters
names(NK_id) <- levels(NK)
NK <- RenameIdents(NK, NK_id)
NK$cell.type <- Idents(NK)
saveRDS(NK, "ALL_LUAD/NK.rds")

# Dot plot and UMAP
Idents(NK) <- factor(Idents(NK), levels = c("NK_XCL1", "NK_FGFBP2", "CD8T_Teff_FGFBP2", "CD8T_gdT_TRDV1", "ILC3"))
DotPlot(NK, features = NK_feature, cols = c("lightgrey", "red"), scale.min = 10, scale.max = 70, col.min = -1.5, col.max = 1.5) + RotatedAxis() + scale_colour_gradientn(colours = c("#052550", "#3984BB", "#95C6DF", "white", "#F3AA89", "red", "#a30018"))
set3_color <- RColorBrewer::brewer.pal(12, "Set2")
set3_color <- c(set3_color, "#9176ff", "#F1C3A2", "#E0C68C")
set3_color[9] <- "#9176ff"
set3_color[2] <- "#F7CE5C"

DimPlot(NK, pt.size = 0.05, cols = set3_color) + theme(
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1, linetype = "solid"),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()
) + labs(x = "UMAP1", y = "UMAP2", title = "NK_sub_cell_type") 


# combine CD4, CD8 and NK
CD4 <- readRDS("ALL_LUAD/CD4.rds")
CD8 <- readRDS("ALL_LUAD/CD8.rds")
NK <- readRDS("ALL_LUAD/NK.rds")
T_part <- merge(CD4, NK)
T_part <- merge(T_part, CD8)
saveRDS(T_part, "ALL_LUAD/T_part.rds")

#######################################
#                                     #
#             GSVA analysis           #
#                                     #
#######################################

T_part <- T_part_IIT

CP_geneset <- msigdbr(category = "C2") %>%
  filter(gs_name %in% c(
    "REACTOME_COSTIMULATION_BY_THE_CD28_FAMILY",
    "REACTOME_TCR_SIGNALING"
  ))
H_geneset <- msigdbr(category = "H")
geneset <- rbind(CP_geneset, H_geneset)
CP_geneset_name <- CP_geneset %>%
  pull(gs_name)
IFNG_geneset <- H_geneset %>% filter(gs_name == "HALLMARK_INTERFERON_GAMMA_RESPONSE")
H_geneset_name <- c(
    "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_GLYCOLYSIS",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_PI3K_AKT_MTOR_SIGNALING"
  )
Hallmark_name <- c(CP_geneset_name, H_geneset_name) %>% unique()
Exhaustion_marker <- list(c("PDCD1", "LAYN", "HAVCR2", "LAG3", "CTLA4", "TIGIT", "TOX", "VSIR", "BTLA", "ENTPD1"))
Cytokine_receptor <- list(c("CSF1", "CSF2", "ADAM19", "ADAM8", "ADAM12", "CD70", "IL12RB2", "IL17A", "IL17F", "IL10RA", "IL1R1", "IL1R2", "IL21R", "IL21", "IL26", "IL2RA", "IL2RB", "IL2RG", "IL32", "IL6R", "TGFB1", "TGFBR2", "TGFBR3"))
Cytotoxic <- list(c("GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "GNLY", "NKG7", "IFNG"))
Chemokine_and_receptor <-  list(c("CCR4", "CCR5", "CCR6", "CCR7", "CCR8", "CXCR3", "CXCR4", "CXCR5", "CXCR6", "CCL3", "CCL4", "CCL5", "CCL20", "CXCL13", "CXCL8", "XCL1"))
Treg_signature <- list(c("FOXP3", "IKZF2", "IKZF4", "IL2RA", "ENTPD1", "CCR4", "ICOS", "IL10RA", "TGFB1", "TIGIT", "CTLA4", "LAG3", "HAVCR2", "PDCD1", "CCR8"))
for (name in Hallmark_name) {
  print(name)
  HALLMARK <- geneset %>%
    filter(gs_name == name) %>%
    pull(human_gene_symbol) %>%
    list()
  T_part <- AddModuleScore(T_part, features = HALLMARK, name = name)
}
T_part <- AddModuleScore(T_part, features = Exhaustion_marker, name = "Exhaustion_marker")
T_part <- AddModuleScore(T_part, features = Cytokine_receptor, name = "Cytokine_and_receptor")
T_part <- AddModuleScore(T_part, features = Cytotoxic, name = "Cytotoxic")
T_part <- AddModuleScore(T_part, features = Chemokine_and_receptor, name = "Chemokine_and_receptor")
T_part <- AddModuleScore(T_part, features = Treg_signature, name = "Treg_signature")

heatmap_data <- T_part@meta.data[, c("sampleID", "cell.type", paste0(Hallmark_name, "1"),
"Exhaustion_marker1", "Cytokine_and_receptor1", "Cytotoxic1", "Chemokine_and_receptor1", "Treg_signature1")] %>%
left_join(Sample, by = c("sampleID" = "sampleID")) %>%
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
        )
)
write_csv(heatmap_data, "heatmap_data.csv")
heatmap_by_type <- heatmap_data %>%
  select(-PD1) %>%
  group_by(cell.type) %>%
  summarise(across(ends_with("1"), mean)) %>%
  ungroup()
write_csv(heatmap_by_type, "heatmap_by_type.csv")