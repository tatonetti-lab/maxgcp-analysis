library(stats)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
input_pheno_path <- args[1]
output_pca_path <- args[2]
output_loading <- args[3]


pheno <- fread(input_pheno_path, header=TRUE)
pheno_scale <- scale(pheno[,-1])
pca_result <- prcomp(pheno_scale, center = TRUE, scale. = TRUE)
pcs <- pca_result$x
pc_data <- data.frame(IID = pheno$'#IID', pcs)

write.table(pc_data, file = sub(".csv", ".tsv", output_pca_path), sep = "\t", row.names = FALSE, quote = FALSE)

loadings <- as.data.frame(pca_result$rotation)
fwrite(loadings, output_loading)

loadings <- cbind(phenotype = rownames(loadings), loadings)
rownames(loadings) <- NULL


pca_plot <- ggplot(pc_data, aes(x = PC1, y = PC2)) +
  geom_point() +
  xlab(paste0("PC1 (", round(summary(pca_result)$importance[2,1] * 100, 2), "% variance)")) +
  ylab(paste0("PC2 (", round(summary(pca_result)$importance[2,2] * 100, 2), "% variance)")) +
  ggtitle("PCA Plot of Phenotypes") +
  theme_minimal()

ggsave("./Documents/TLAB/Plink/sample_data/pca_plot.png", plot = pca_plot, width = 10, height = 7)
print(pca_plot)