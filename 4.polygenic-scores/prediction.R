library(dplyr)
library(caret)
library(e1071)

#args <- commandArgs(trailingOnly = TRUE)
# data_path <- args[1]
# covariate_file <- args[2]
# phenotype_file <- args[3]
# prs_naive_file <- args[4]
# prs_pca_file <- args[5]
# prs_maxgcp_file <- args[6]
# prs_maxh_file <- args[7]

data_path <- "./Documents/TLAB/Plink/sample_data"
covariate_file <- file.path(data_path, "/covar.tsv")
phenotype_file <- file.path(data_path, "/phenotypes.tsv")
prs_naive_file <- file.path(data_path, "/prs/prs_cov/naiveprs_covariate.tsv")
prs_pca_file <- file.path(data_path, "/prs/prs_cov/pcaprs_covariate.tsv")
prs_maxgcp_file <- file.path(data_path, "/prs/prs_cov/maxgcprs_covariate.tsv")
prs_maxh_file <- file.path(data_path, "/prs/prs_cov/maxhprs_covariate.tsv")

pheno_df <- read.csv(phenotype_file, sep='\t')
prs_naive_df <- read.csv(prs_naive_file, sep='\t')
prs_pca_df <- read.csv(prs_pca_file, sep='\t')
prs_maxgcp_df <- read.csv(prs_maxgcp_file, sep='\t')
prs_maxh_df <- read.csv(prs_maxh_file, sep='\t')

#phenos <- c("pheno_0", "pheno_1", "pheno_2", "pheno_3", "pheno_4", "pheno_5", "pheno_6", "pheno_7", "pheno_8", "pheno_9")
phenos <- c("pheno_0", "pheno_1")

merge_train <- function(covar_data, pheno_data, target_pheno){
  target_pattern <- paste0("_", target_pheno)
  covar_data <- covar_data %>% select(X.IID, starts_with("covar_"), ends_with(target_pattern))
  pheno_data <- pheno_data %>% select(X.IID, all_of(target_pheno))
  data <- merge(covar_data, pheno_data, by="X.IID")
  
  return(data)
}


perform_train_test <- function(data, phenotype_col, train_ratio=0.7) {
  set.seed(66)

  train_index <- createDataPartition(data[[phenotype_col]], p=train_ratio, list=FALSE)
  train_data <- data[train_index,]
  test_data <- data[-train_index,]
  train_data <- train_data %>% select(-X.IID)
  test_data <- test_data %>% select(-X.IID)
  model <- train(as.formula(paste(phenotype_col, "~ .")), data=train_data, method="lm")
  predictions <- predict(model, test_data)
  

  predicted_classes <- ifelse(predictions > median(train_data[[phenotype_col]]), 1, 0)
  actual_classes <- ifelse(test_data[[phenotype_col]] > median(train_data[[phenotype_col]]), 1, 0)

  confusion <- confusionMatrix(as.factor(predicted_classes), as.factor(actual_classes))
  
  return(confusion)
}


performance_results <- data.frame()

for (pheno in phenos) {
  print(paste("merge data for", pheno))
  naive_data <- merge_train(prs_naive_df, pheno_df, pheno)
  pca_data <- merge(prs_pca_df, pheno_df[, c("X.IID", all_of(pheno))], by="X.IID")
  maxgcp_data <- merge_train(prs_maxgcp_df, pheno_df, pheno)
  maxh_data <- merge(prs_maxh_df, pheno_df[, c("X.IID", all_of(pheno))], by="X.IID")

  print(paste("train model for", pheno))
  naive_confusion <- perform_train_test(naive_data, pheno)
  #pca_confusion <- perform_train_test(pca_data, pheno)
  #maxgcp_confusion <- perform_train_test(maxgcp_data, pheno)
  #maxh_confusion <- perform_train_test(maxh_data, pheno)
  
  print(paste("Record Performance for", pheno))
  naive_performance <- data.frame(Model="Naive",Phenotype=pheno, Accuracy=naive_confusion$overall["Accuracy"], Sensitivity=naive_confusion$byClass["Sensitivity"], 
                                  Specificity=naive_confusion$byClass["Specificity"],Recall=naive_confusion$byClass["Recall"], Precision=naive_confusion$byClass["Precision"],
                                  F1=naive_confusion$byClass["F1"], row.names=NULL)
  #pca_performance <- data.frame(Model="PCA",Phenotype=pheno, Accuracy=pca_confusion$overall["Accuracy"], Sensitivity=pca_confusion$byClass["Sensitivity"],
  #                              Specificity=pca_confusion$byClass["Specificity"],Recall=pca_confusion$byClass["Recall"], Precision=pca_confusion$byClass["Precision"],
  #                              F1=pca_confusion$byClass["F1"],row.names=NULL)
  #maxgcp_performance <- data.frame(Model="MaxGCP",Phenotype=pheno, Accuracy=maxgcp_confusion$overall["Accuracy"], Sensitivity=maxgcp_confusion$byClass["Sensitivity"],
  #                                 Specificity=maxgcp_confusion$byClass["Specificity"],Recall=maxgcp_confusion$byClass["Recall"], Precision=maxgcp_confusion$byClass["Precision"],
  #                                 F1=maxgcp_confusion$byClass["F1"],row.names=NULL)
  #maxh_performance <- data.frame(Model="MaxH",Phenotype=pheno, Accuracy=maxh_confusion$overall["Accuracy"], Sensitivity=maxh_confusion$byClass["Sensitivity"],
  #                              Specificity=maxh_confusion$byClass["Specificity"],Recall=maxh_confusion$byClass["Recall"], Precision=maxh_confusion$byClass["Precision"],
  #                              F1=maxh_confusion$byClass["F1"],row.names=NULL)
  
  #performance_results <- rbind(performance_results, naive_performance, pca_performance, maxgcp_performance, maxh_performance)
  performance_results <- rbind(performance_results, naive_performance)
}




