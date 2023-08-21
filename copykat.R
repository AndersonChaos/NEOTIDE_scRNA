library(copykat)
library(Seurat)
library(tidyverse)
library(AnnoProbe)
Epi <- readRDS("ALL_LUAD/Epi.rds")
Epi <- subset(Epi, sampleID %in% IIT_id)
Stromal <- readRDS("ALL_LUAD/Stromal.rds")
set.seed(3047)
sampleIDs <- unique(Epi$sampleID)
print(sampleIDs)
for (sample in sampleIDs) {
    sample_ID <- str_extract(sample, "P[0-9]+")
    print(paste0("Processing ", sample_ID))
    sample <- subset(Epi, sampleID == sample_ID)
    sample_stromal <- subset(Stromal, sampleID == sample_ID)
    if(length(unique(sample$cellID)) < 10){
        next
    }
    stromal_size <- ncol(sample_stromal)
    if (stromal_size > 100) {
        stromal_cellID <- sample_stromal$cellID
        # randomly choose 100 cell
        stromal_cellID <- sample(stromal_cellID, 100)
        sample_stromal <- subset(sample_stromal, cellID %in% stromal_cellID)
    }
    sample_control <- merge(sample, sample_stromal)
    normal_cells <- colnames(sample_stromal)
    dfcount <- sample_control@assays$RNA@counts
    geneInfor <- annoGene(rownames(dfcount), "SYMBOL", "human")
    geneInfor <- geneInfor[with(geneInfor, order(chr, start)), c(1, 4:6)]
    geneInfor <- geneInfor[!duplicated(geneInfor[, 1]), ]
    head(geneInfor)
    geneInfor <- geneInfor %>%
        filter(chr != "chrM") %>%
        filter(chr != "chrX") %>%
        filter(chr != "chrY")
    dfcount <- dfcount[rownames(dfcount) %in% geneInfor[, 1], ]
    copykat.test <- copykat(rawmat = dfcount, id.type = "S", ngene.chr = 5, win.size = 25, KS.cut = 0.1, sam.name = paste0("copykat_ALL_LUAD/", sample_ID), distance = "euclidean",
        output.seg = FALSE, plot.genes = TRUE, n.cores = 20, norm.cell.names = normal_cells)
}

# calculate CNV score
sampleIDs <- IIT_id
copykat_CNA_res <- data.frame()
pred_res <- data.frame()
for (id in sampleIDs) {
  print(id)
  if (file.exists(paste0("copykat_ALL_LUAD/", id, "_copykat_CNA_results.txt")) == FALSE) next
  copy_res <- read_table(
    paste0("copykat_ALL_LUAD/", id, "_copykat_CNA_results.txt")
  )
  predicted <- read_table(
    paste0("copykat_ALL_LUAD/", id, "_copykat_prediction.txt")
  ) %>% filter(copykat.pred != "not.defined")
  copy_res <- pivot_longer(
      copy_res,
      cols = -c("chrom", "chrompos", "abspos"), names_to = "cellID", values_to = "copy_num"
  )
  copy_res$cellID <- str_replace(copy_res$cellID, "\\.", "-")
  copykat_CNA_res <- bind_rows(copykat_CNA_res, copy_res)
  pred_res <- bind_rows(pred_res, predicted)
}
CNVs <- pivot_wider(copykat_CNA_res, names_from = cellID, values_from = copy_num)

# Copykat plot
my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
chr <- as.numeric(CNVs$chrom) %% 2 + 1
rbPal1 <- colorRampPalette(c("black", "grey"))
CHR <- rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR, CHR)

rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
com.preN <- pred_res$copykat.pred
pred <- rbPal5(2)[as.numeric(factor(com.preN))]

cells <- rbind(pred, pred)
col_breaks <- c(seq(-1, -0.4, length = 50), seq(-0.4, -0.2, length = 150), seq(-0.2, 0.2, length = 600), seq(0.2, 0.4, length = 150), seq(0.4, 1, length = 50))

jpeg("copykat_pred.jpeg", width = 10, height = 10)
heatmap.3(t(CNVs[, 4:ncol(CNVs)]),
    dendrogram = "r", distfun = function(x) parallelDist::parDist(x, threads = 4, method = "euclidean"), hclustfun = function(x) hclust(x, method = "ward.D2"),
    ColSideColors = chr1, RowSideColors = cells, Colv = NA, Rowv = TRUE,
    notecol = "black", col = my_palette, breaks = col_breaks, key = TRUE,
    keysize = 1, density.info = "none", trace = "none",
    cexRow = 0.1, cexCol = 0.1, cex.main = 1, cex.lab = 0.1,
    symm = F, symkey = F, symbreaks = T, cex = 1, cex.main = 4, margins = c(10, 10)
)
legend("topright", paste("pred.", names(table(com.preN)), sep = ""), pch = 15, col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex = 0.6, bty = "n")
dev.off()
