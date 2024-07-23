library(qqman)
library(data.table)
library(gridExtra)

# Plot for Naive GWAS 
result_path <- "Documents/TLAB/Plink/sample_data/gwas/"
phenotypes <- paste0("pheno_", 0:9)
PCs <- paste0("PC", 1:10)

convert_chr_to_numeric <- function(chr) {
  if (chr == "X") return(23)
  return(as.numeric(chr))
}

naive_plot <- function(pheno, result_path) {
  file_path <- file.path(result_path, paste0("naive.", pheno, ".glm.linear"))
  gwas_results <- fread(file_path, header = TRUE)
  colnames(gwas_results)[which(names(gwas_results) == "#CHROM")] <- "CHR"
  colnames(gwas_results)[which(names(gwas_results) == "POS")] <- "BP"
  colnames(gwas_results)[which(names(gwas_results) == "P")] <- "P"
  colnames(gwas_results)[which(names(gwas_results) == "ID")] <- "SNP"
  
  gwas_results$CHR <- sapply(gwas_results$CHR, convert_chr_to_numeric)
  
  manhattan(gwas_results, main = paste("Manhattan Plot for", pheno), col = c("blue2", "coral2"),
            suggestiveline = -log10(1e-5), genomewideline = -log10(5e-8))
}

result_path <- "Documents/TLAB/Plink/sample_data/gwas"
phenotypes <- paste0("pheno_", 0:9)

png("Documents/TLAB/Plink/sample_data/naive_manhattan_plots.png", width = 20, height = 8, units = "in", res = 300)
par(mfrow = c(2, 5), oma = c(0, 0, 2, 0))

for (pheno in phenotypes) {
  naive_plot(pheno, result_path)
}

dev.off()

# Plot for PCA
result_path <- "Documents/TLAB/Plink/sample_data/pca_gwas/"
convert_chr_to_numeric <- function(chr) {
  if (chr == "X") return(23)
  return(as.numeric(chr))
}

pca_plot <- function(PC, result_path) {
  file_path <- file.path(result_path, paste0("pca.", PC, ".glm.linear"))
  gwas_results <- fread(file_path, header = TRUE)
  colnames(gwas_results)[which(names(gwas_results) == "#CHROM")] <- "CHR"
  colnames(gwas_results)[which(names(gwas_results) == "POS")] <- "BP"
  colnames(gwas_results)[which(names(gwas_results) == "P")] <- "P"
  colnames(gwas_results)[which(names(gwas_results) == "ID")] <- "SNP"
  
  gwas_results$CHR <- sapply(gwas_results$CHR, convert_chr_to_numeric)
  
  manhattan(gwas_results, main = paste("Manhattan Plot for", pheno), col = c("blue2", "coral2"),
            suggestiveline = -log10(1e-5), genomewideline = -log10(5e-8))
}


png("Documents/TLAB/Plink/sample_data/pca_manhattan_plots.png", width = 20, height = 8, units = "in", res = 300)
par(mfrow = c(2, 5), oma = c(0, 0, 2, 0))

for (PC in PCs) {
  pca_plot(PC, result_path)
}

dev.off()

# MaxGCP Plot
result_path <- "Documents/TLAB/Plink/sample_data/maxgcp_gwas/"
gcp_plot <- function(pheno, result_path) {
  file_path <- file.path(result_path, paste0("maxgcp.", pheno, ".glm.linear"))
  gwas_results <- fread(file_path, header = TRUE)
  colnames(gwas_results)[which(names(gwas_results) == "#CHROM")] <- "CHR"
  colnames(gwas_results)[which(names(gwas_results) == "POS")] <- "BP"
  colnames(gwas_results)[which(names(gwas_results) == "P")] <- "P"
  colnames(gwas_results)[which(names(gwas_results) == "ID")] <- "SNP"
  
  gwas_results$CHR <- sapply(gwas_results$CHR, convert_chr_to_numeric)
  
  manhattan(gwas_results, main = paste("Manhattan Plot for", pheno), col = c("blue2", "coral2"),
            suggestiveline = -log10(1e-5), genomewideline = -log10(5e-8))
}

png("Documents/TLAB/Plink/sample_data/maxgcp_manhattan_plots.png", width = 20, height = 8, units = "in", res = 300)
par(mfrow = c(2, 5), oma = c(0, 0, 2, 0))

for (pheno in phenotypes) {
  gcp_plot(pheno, result_path)
}

dev.off()

# MaxH Plot
result_path <- "Documents/TLAB/Plink/sample_data/maxh_gwas/"
maxh_plot <- function(pheno, result_path) {
  file_path <- file.path(result_path, paste0("maxh.", pheno, ".glm.linear"))
  gwas_results <- fread(file_path, header = TRUE)
  colnames(gwas_results)[which(names(gwas_results) == "#CHROM")] <- "CHR"
  colnames(gwas_results)[which(names(gwas_results) == "POS")] <- "BP"
  colnames(gwas_results)[which(names(gwas_results) == "P")] <- "P"
  colnames(gwas_results)[which(names(gwas_results) == "ID")] <- "SNP"
  
  gwas_results$CHR <- sapply(gwas_results$CHR, convert_chr_to_numeric)
  
  manhattan(gwas_results, main = paste("Manhattan Plot for", pheno), col = c("blue2", "coral2"),
            suggestiveline = -log10(1e-5), genomewideline = -log10(5e-8))
}

png("Documents/TLAB/Plink/sample_data/maxh_manhattan_plots.png", width = 20, height = 8, units = "in", res = 300)
par(mfrow = c(2, 5), oma = c(0, 0, 2, 0))

for (pheno in phenotypes) {
  maxh_plot(pheno, result_path)
}

dev.off()



