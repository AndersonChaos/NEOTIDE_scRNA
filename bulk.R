library(GSVA)
library(readxl)
library(tidyverse)
library(ggpubr)

bulk_data <- read_excel("bulk/2104_TPM.xlsx") %>% as.data.frame()
bulk_data <- read.csv("bulk/2104_TP.csv") %>% as.data.frame()
head(bulk_data)
rownames(bulk_data) <- bulk_data[, 1]
bulk_data <- bulk_data[, -1]
bulk_matrix <- bulk_data
head(bulk_matrix)
geneset <- list(
    "Tex" = c("HAVCR2", "PDCD1", "LAG3", "CXCL13", "TIGIT", "CTLA4", "ENTPD1", "LAYN"),
    "Treg" = c("CCR8", "FOXP3", "IL2RA")
)
enrich <- gsva(as.matrix(bulk_matrix), geneset, method = "ssgsea", kcdf = "Gaussian")
enrich_res <- enrich %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    mutate(
        Tex_cutoff = ifelse(Tex > median(Tex), "High", "Low"),
        Treg_cutoff = ifelse(Treg > median(Treg), "High", "Low"),
        Cut = ifelse(Tex_cutoff == "Low" & Treg_cutoff == "High", "Signature", "Others")
    )
View(enrich_res)
write_csv(enrich_res, "bulk/enrich_res_2104.csv")
