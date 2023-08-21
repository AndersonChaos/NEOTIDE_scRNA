library(Seurat)
library(ggpubr)
library(tidyverse)
T_part <- readRDS("ALL_LUAD/T_part.rds")
###
tcr_files <- list.files(path = "scTCR_clean/", pattern = "*.csv")
count <- 0
for (file in tcr_files) {
    tcr <- read_csv(paste0("scTCR_clean/", file))
    if (count == 0) {
        tcr_all <- tcr
        count <- count + 1
    } else {
        tcr_all <- rbind(tcr_all, tcr)
    }
}
T_meta <- T_part@meta.data
tcr_all <- tcr_all %>%
    distinct()
T_meta_TCR <- T_meta %>%
    left_join(tcr_all, by = c(cellID = "barcode"))
T_meta_TCR <- T_meta_TCR %>%
    as_tibble() %>%
    group_by(clonetype) %>%
    mutate(cloneSize = n()) %>%
    ungroup() %>%
    mutate(
        expanded = ifelse(cloneSize > 10, "expanded", "not-expanded")
    ) %>%
    as.data.frame()
rownames(T_meta_TCR) <- T_meta_TCR$cellID


TCR_ALL <- T_meta_TCR %>%
    mutate(
        group = ifelse(sampleID %in% IIT_id, "2104", "Real_World")
    ) %>%
    select("sampleID", "clonetype", "cell.type", "cloneSize", "cellID", "group")
TCR_ALL <- TCR_ALL %>%
    mutate(
        cell_type = case_when(
            cell.type == "CD4T_Th1-like CXCL13" ~ "CD4T_Th1-like_CXCL13",
            .default = cell.type
        )
    ) %>%
    filter(!is.na(clonetype))
TCR_ALL <- TCR_ALL %>% left_join(Sample, by = c(sampleID = "sampleID"), multiple = "first") %>%    mutate(
        isMPR = ifelse(`Pathological Response` %in% c("MPR", "pCR"), "MPR", "Non-MPR"),
        isMPR = factor(isMPR, levels = c("Non-MPR", "MPR")),
        sampleID_Smoke = case_when(
            str_detect(Smoking_History, "Y") ~ str_c(sampleID, "_ES"),
            str_detect(Smoking_History, "N") ~ str_c(sampleID, "_NS"),
        ),
        PRR = as.double(`Pathological Response Rate`),
        PRR = ifelse(is.na(PRR), 0.1, PRR),
        NOD_ID = case_when(
            sampleID == "P343" ~ "NOD01",
            sampleID == "P438" ~ "NOD02",
            sampleID == "P589"~ "NOD03",
            sampleID == "P590"~ "NOD04",
            sampleID == "P591"~ "NOD05",
            sampleID == "P529"~ "NOD06",
            sampleID == "P531"~ "NOD07",
            sampleID == "P519"~ "NOD08",
            sampleID == "P547"~ "NOD09",
            sampleID == "P549"~ "NOD10",
            sampleID == "P586"~ "NOD11",
            .default = "Real_World"
        ),
        NODID_Smoke = case_when(
            str_detect(Smoking_History, "Y") ~ str_c(NOD_ID, "_ES_", PRR),
            str_detect(Smoking_History, "N") ~ str_c(NOD_ID, "_NS_", PRR),
        ),
    )
Tex_clone_use <- TCR_ALL %>%
    filter(cell_type %in% c("CD8T_Tex_CXCL13") & cloneSize >= 10 & sampleID %in% IIT_id) %>%
    select(clonetype, cloneSize) %>%
    distinct() %>%
    pull(clonetype)
nTex_clone <- TCR_ALL %>%
    filter(cell_type %in% c("CD8T_Tex_CXCL13") & cloneSize >= 10) %>%
    select(sampleID, clonetype, cloneSize) %>%
    distinct() %>%
    count(sampleID) %>%
    rename(nTex = n)
nTh1_clone <- TCR_ALL %>%
    filter(cell_type == "CD4T_Th1-like_CXCL13" & cloneSize >= 3) %>%
    select(sampleID, clonetype, cloneSize) %>%
    distinct() %>%
    count(sampleID) %>%
    rename(nTh1 = n)
nTreg_clone <- TCR_ALL %>%
    filter(cell_type == "CD4T_Treg_CCR8" & cloneSize >= 3) %>%
    select(sampleID, clonetype, cloneSize) %>%
    distinct() %>%
    count(sampleID) %>%
    rename(nTreg = n)
group_info <- TCR_ALL %>%
    select(sampleID, group) %>%
    distinct()
clone_num <- nTex_clone %>%
    full_join(nTh1_clone, by = "sampleID") %>%
    full_join(nTreg_clone, by = "sampleID") %>%
    full_join(group_info, by = "sampleID") %>%
    mutate(
        EGFR_Mut = case_when(
            group == "2104" ~ "EGFR_Mut",
            group == "Real_World" ~ "EGFR_Mut",
            .default = "EGFR_WT_LUAD"
        ),
        EGFR_Mut = factor(EGFR_Mut, levels = c("EGFR_WT_LUAD", "EGFR_Mut")),
        nTex = ifelse(is.na(nTex), 0, nTex),
        nTh1 = ifelse(is.na(nTh1), 0, nTh1),
        nTreg = ifelse(is.na(nTreg), 0, nTreg)
    )
write_csv(clone_num, "cellratio/clone_num.csv")
write_csv(TCR_ALL, "cellratio/TCR_All.csv")