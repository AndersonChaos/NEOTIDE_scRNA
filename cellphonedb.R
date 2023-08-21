library(tidyverse)
library(ggpubr)

res_dir <- "ALL_LUAD/statistical_analysis_"
pvalues <- read_delim(paste0(res_dir, "pvalues_Myeloid_Epi.txt"), delim = "\t")
means <- read_delim(paste0(res_dir, "means_Myeloid_Epi.txt"), delim = "\t")
pvalues <- pvalues %>%
    pivot_longer(cols = -c("id_cp_interaction", "interacting_pair", "partner_a", "partner_b", "gene_a", "gene_b", "secreted", "receptor_a", "receptor_b", "annotation_strategy", "is_integrin"), names_to = "interacting_cells", values_to = "pvalue")
means <- means %>%
    pivot_longer(cols = -c("id_cp_interaction", "interacting_pair", "partner_a", "partner_b", "gene_a", "gene_b", "secreted", "receptor_a", "receptor_b", "annotation_strategy", "is_integrin"), names_to = "interacting_cells", values_to = "means")
integrated_data <- left_join(means, pvalues)
integrated_data <- integrated_data %>%
    mutate(log_p_value = -log10(pvalue + 1e-5)) # to avoid 0
write_csv(integrated_data, "cellratio/integrated_cellphoneDB.csv")

bubble_final <- read_csv("ALL_LUAD/bubble final.csv")

ggplot(bubble_final %>% mutate(means = ifelse(means > 1.5, 1.5, means))) +
    geom_point(aes_string(x = "interacting_cells", y = "interacting_pair", size = "log_p_value", color = "means")) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12)) +
    labs(x = NULL, y = NULL, size = "-log10(pvalue)") +
    scale_colour_gradientn(colours = c("#052550", "#3984BB", "#95C6DF", "#ececec", "#F3AA89", "red", "#a30018"))

ggsave("ALL_LUAD/cellphoneDB.pdf", width = 14, height = 8)
