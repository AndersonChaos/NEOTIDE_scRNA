library(tidyverse)
library(ggpubr)
library(SeuratDisk)
library(rstatix)
library(Seurat)
library(harmony)
library(readxl)

sample <- list.files("Raw_Matrix/Final_Use")
print(sample)
count <- 0
use_id <- c("P33", "P52", "P64", "P70", "P105", "P122", "P130", "P172", "P209", "P218", "P223", "P266", "P269", "P270", "P293","P345", "P348", "P390","P394", "P399", "P481", "P483", "P498", "P509", "P523", "P528", "P533","P551","P574", "P579", "P592", "P603", "P605", "P609")
IIT_id <- c("P343", "P438", "P519", "P529", "P531", "P547", "P549", "P586", "P589", "P590", "P591")
for (name in sample){
  print(name)
  raw <- Read10X(data.dir = paste0("Raw_Matrix/Final_Use/", name))
  raw <- CreateSeuratObject(counts = raw, project = name, min.cells = 3)
  raw[["percent.mt"]] <- PercentageFeatureSet(raw, pattern = "^MT-")
  raw <- subset(raw,
    subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA > 1000 & nCount_RNA < 40000
  )
  saveRDS(raw, paste0("LUAD", name, ".rds"))
  print(raw)
  raw$sampleID <- name
  raw$cellID <- paste0(raw$sampleID, "_", colnames(raw))
  raw <- RenameCells(raw, new.names = raw$cellID)
  rownames(raw@meta.data) <- raw@meta.data$cellID
      if (count == 0) {
      scRNA <- raw
      count <- count + 1
    } else {
      scRNA <- merge(scRNA, raw)
    }
}

doublet.cell.IDs <- unique(
  c(
    WhichCells(scRNA, expression = CD79A > 0 & CD3D > 0),
    WhichCells(scRNA, expression = CD79A > 0 & CD3E > 0),
    WhichCells(scRNA, expression = CD79A > 0 & CD14 > 0),

    WhichCells(scRNA, expression = CD19 > 0 & CD14 > 0),
    WhichCells(scRNA, expression = CD19 > 0 & CD3D > 0),
    WhichCells(scRNA, expression = CD19 > 0 & CD3E > 0),

    WhichCells(scRNA, expression = CD14 > 0 & CD3D > 0),
    WhichCells(scRNA, expression = CD14 > 0 & CD3E > 0),
    WhichCells(scRNA, expression = CD68 > 0 & CD3D > 0),
    WhichCells(scRNA, expression = CD68 > 0 & CD3E > 0)
  )
)
scRNA <- scRNA[, !scRNA$cellID %in% doublet.cell.IDs]
scRNA <- subset(scRNA, subset = nFeature_RNA > 600 & nFeature_RNA < 6000 & nCount_RNA > 1000 & nCount_RNA < 40000 & percent.mt < 10)
all.genes <- rownames(scRNA)
scRNA
scRNA <- scRNA %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData(features = all.genes)
scRNA <- scRNA %>%
  RunPCA(features = VariableFeatures(object = scRNA)) %>%
  RunHarmony(group.by.vars = "sampleID") %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.8) %>%
  RunUMAP(reduction = "harmony", dims = 1:20)

scMeta <- scRNA@meta.data
Sample <- read_excel("Sample.xlsx") %>%
  mutate(sampleID = str_extract(SampleID, "P[0-9]+")) %>%
  select(sampleID, Gender, Age, Histology, Cycles, PD1, `Pathological Response`, EGFR)

scMeta <- left_join(scMeta, Sample, by = c("sampleID"), keep = FALSE)
scMeta <- table(scRNA$major.cell.type, scRNA$sampleID) %>% as.data.frame() %>% mutate(total = sum(Freq))
write_csv(scMeta, "cellratio/scMeta.csv")
scMeta <- scMeta %>%
  mutate(data_from = ifelse(sampleID %in% IIT_id, "2104", "Real_World"))
scMeta %>% filter(data_from == "2104") %>% count(major.cell.type) %>% write_csv("cellratio/2104_QC.csv")
scMeta %>% filter(data_from == "Real_World") %>% count(major.cell.type) %>% write_csv("cellratio/RW_QC.csv")

scMeta$sampleID <- as.factor(scMeta$sampleID)
rownames(scMeta) <- scMeta$cellID
scRNA@meta.data <- scMeta

feature <- list(
  "im" = c("PTPRC"),
  "T" = c("TRAC", "CD3D", "CD3E", "CD3G"),
  "NK" = c("NKG7", "GNLY", "NCR1", "NCAM1"),
  "Myeloid" = c("CD68", "LYZ", "CD14"),
  "B" = c("CD79A", "CD19", "MS4A1", "IGHM", "IGHG1"),
  "Endo" = c("PECAM1", "FLT1"),
  "Fibro" = c("COL1A1", "COL6A1"),
  "Epi" = c("EPCAM", "KRT18"),
  "Neutrophil" = c("CSF3R", "FCGR3B", "S100A9"),
  "Mast" = c("KIT", "CPA3", "TPSAB1")
)
Idents(scRNA) <- scRNA$seurat_clusters
DotPlot(scRNA, features = feature, group.by = "seurat_clusters") + RotatedAxis()
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)

clusters_id <- c("T", "B", "T", "NK", "doublet", "Epi", "Myeloid", "Myeloid", "Myeloid", "Endo", "Fibro", "Fibro", "Fibro", "B", "T", "Myeloid", "T", "Epi", "Mast", "Fibro", "Epi", "Epi", "T", "Epi", "Epi", "Fibro", "B", "Myeloid", "Fibro", "Epi", "Endo", "Fibro", "B", "Epi", "Endo", "Mast", "Endo", "Myeloid")

Idents(scRNA) <- scRNA$seurat_clusters
names(clusters_id) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, clusters_id)
scRNA$major.cell.type <- Idents(scRNA)
saveRDS(scRNA, "ALL_LUAD/ALL_LUAD.rds")

DimPlot(scRNA, pt.size = 3e-06, label = TRUE, label.size = 10, group.by = "seurat_clusters")
DimPlot(scRNA, pt.size = 3e-06, label = TRUE, label.size = 10)


# Save the data
Idents(scRNA) <- scRNA$major.cell.type
Epi <- subset(scRNA, idents = "Epi")
T_part <- subset(scRNA, idents = c("T", "NK"))
B_part <- subset(scRNA, idents = "B")
Myeloid_part <- subset(scRNA, idents = c("Myeloid", "Mast"))
Stromal <- subset(scRNA, idents = c("Fibro", "Endo"))
doublets <- subset(scRNA, idents = "doublet")
All_immune <- subset(scRNA, idents = c("T", "NK", "B", "Myeloid", "Mast"))
scRNA <- subset(scRNA, idents = c("T", "NK", "B", "Myeloid", "Mast", "Epi", "Fibro", "Endo"))
DimPlot(scRNA)
saveRDS(Epi, "ALL_LUAD/Epi.rds")
saveRDS(T_part, "ALL_LUAD/T_part.rds")
saveRDS(B_part, "ALL_LUAD/B_part.rds")
saveRDS(Myeloid_part, "ALL_LUAD/Myeloid_part.rds")
saveRDS(Stromal, "ALL_LUAD/Stromal.rds")
saveRDS(scRNA, "ALL_LUAD/ALL_LUAD_filtered.rds")
saveRDS(All_immune, "ALL_LUAD/All_immune.rds")

IIT <- scRNA %>% subset(sampleID %in% IIT_id)
RWC <- scRNA %>% subset(sampleID %in% IIT_id, invert = TRUE)

Real_Matrix <- RWC[["RNA"]]@counts
Matrix::writeMM(Real_Matrix, "Real_Matrix.mtx")
IIT_Matrix <- IIT[["RNA"]]@counts
Matrix::writeMM(IIT_Matrix, "IIT_Matrix.mtx")
Real_Meta <- RWC@meta.data %>%
  left_join(Sample, by = c("sampleID")) %>%
  mutate(Smoking_History = str_extract(Smoking_History, "(Y|N)"))
IIT_Meta <- IIT@meta.data %>%
  left_join(Sample, by = c("sampleID")) %>%
  mutate(Smoking_History = str_extract(Smoking_History, "(Y|N)"))
write_delim(Real_Meta, "Real_Meta.txt", delim = "\t")
write_delim(IIT_Meta, "IIT_Meta.txt", delim = "\t")

# UMAP plot
Idents(scRNA) <- scRNA$major_cell_type
scRNA_IIT <- scRNA %>% subset(sampleID %in% IIT_id)
scRNA_RWC <- scRNA %>% subset(sampleID %in% IIT_id, invert = TRUE)
set3_color <- RColorBrewer::brewer.pal(9, "Set3")
set3_color <- c(set3_color, "#FED9A6")
set3_color[8] <- "#FED9A6"
set3_color[9] <- "#E78AC3"
DimPlot(scRNA, pt.size = 0.01, cols = set3_color) + theme(
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1, linetype = "solid"),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()) + labs(x = "UMAP1", y = "UMAP2", title = "Major_cell_type")

# Cellratio data

cellratio <- prop.table(table(scRNA$sampleID, scRNA$major.cell.type), margin = 1) # calculate ratio of sub celltype
meta <- as.data.frame(cellratio)
meta <- meta %>% left_join(Sample, by = c("Var1" = "sampleID")) # merge with sample information
meta <- meta %>% # Rename
  rename(
    sampleID = Var1,
    Ratio = Freq,
    Celltype = Var2
  )
meta <- meta %>%
  mutate(
    major_type = ifelse(str_detect(Celltype, "CD4"), "CD4", ifelse(str_detect(Celltype, "CD8"), "CD8", "NK")),
    my_type = ifelse(
      str_detect(Celltype, "Treg"), "Treg", 
      ifelse(str_detect(Celltype, "Tex"), "Tex", ifelse(str_detect(Celltype, "NK"), "NK", "Other"))),
    Smoking_History = str_extract(Smoking_History, "(Y|N)"),
    EGFR = ifelse(str_detect(EGFR, "L858R"), "L858R",
      ifelse(str_detect(EGFR, "exon19del"), "exon19del", ifelse(str_detect(EGFR, "exon20ins"), "exon20ins", "other"))
    ),
    PRR = as.double(`Pathological Response Rate`),
    PRR = ifelse(is.na(PRR), 0, PRR),
    response = ifelse(PRR > median(PRR), "Responder", "Non-Responder"),
    isMPR = factor(`Pathological Response`, levels = c("non-MPR", "MPR")),
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
  ) %>%
  filter(!str_detect(Histology, "SCLC"))