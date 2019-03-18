library(dplyr)
library(Biobase)
devtools::load_all()
# Create a toy dataset for use in examples

# Feature data
set.seed(38)
n_features <- 20
feature_data <- data.frame(Split = "HILIC_pos",
                           Alignment = seq_len(n_features),
                           Mass = runif(n_features, 100, 500),
                           RetentionTime = runif(n_features, 0.5, 8),
                           Column = "HILIC", Mode = "pos",
                           Flag = NA_character_)

# Create Feature ID
round_mz <- as.numeric(feature_data$Mass) %>% round(digits = 4) %>%
  as.character() %>% gsub("[.]", "_", .)
round_rt <- as.numeric(feature_data$RetentionTime) %>% round(digits = 4) %>%
  as.character() %>% gsub("[.]", "_", .)
feature_data$Feature_ID <- paste0("HILIC_neg_", round_mz, "a", round_rt)
feature_data <- select(feature_data, "Feature_ID", everything())

rownames(feature_data) <- feature_data$Feature_ID

# Pheno data
n_samples <- 30

qc_idx <- seq(1, n_samples, length.out = 5) %>% round()
subject_ids <- rep(as.character(seq_len(n_samples/2)), 2)
group <- rep(sample(LETTERS[1:2], n_samples/2, replace = TRUE), 2)
group[qc_idx] <- "QC"
subject_ids[qc_idx] <- "QC"
qc <- ifelse(seq_len(n_samples) %in% qc_idx, "QC", "Sample")
pheno_data <- data.frame(Injection_order = seq_len(n_samples),
                         Sample_ID = paste0("Demo_", seq_len(n_samples)),
                         Subject_ID = subject_ids,
                         Group = group, QC = qc,
                         Time = factor(rep(c(1,2), each = n_samples/2)))

rownames(pheno_data) <- pheno_data$Sample_ID


# Assay data

# Random means for each feature
means <- runif(n_features, 3000, 33000)

# Normally distributed data around the mean
assay_data <- t(sapply(means, function(x) {rnorm(n_samples, x, 0.3*x)}))
# Add drift effect to the data
coefs <- runif(n_samples, 0.4, 0.9) * sample(c(-1, 1), n_samples, replace = TRUE)

# Randomly choose linear or logarithmic trend
for (i in seq_len(nrow(assay_data))) {
  if (rnorm(1) > 0) {
    assay_data[i,] <- assay_data[i,] + means[i] * coefs[i] * log(pheno_data$Injection_order)
  } else {
    assay_data[i,] <- assay_data[i,] + means[i] * coefs[i] * 0.08* pheno_data$Injection_order
  }
}
# Set random indexes to zero (for missing values)
n_missing <- 30
row_zeros <- sample(seq_len(nrow(assay_data)), n_missing, replace = TRUE)
col_zeros <- sample(seq_len(ncol(assay_data)), n_missing, replace = TRUE)
for (i in seq_len(n_missing)) {
  assay_data[row_zeros[i], col_zeros[i]] <- 0
}

# Set dimension names
rownames(assay_data) <- rownames(feature_data)
colnames(assay_data) <- rownames(pheno_data)



pheno_data <- AnnotatedDataFrame(data = pheno_data)
feature_data <- AnnotatedDataFrame(data = feature_data)

example_set <- MetaboSet(exprs = assay_data,
                         phenoData = pheno_data,
                         featureData = feature_data,
                         group_col = "Group",
                         time_col = "Time",
                         subject_col = "Subject_ID",
                         predicted = matrix(NA_real_, nrow = nrow(assay_data),
                                            ncol = ncol(assay_data),
                                            dimnames = dimnames(assay_data)))




usethis::use_data(example_set, overwrite = TRUE)