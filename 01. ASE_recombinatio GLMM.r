# ==============================================================================
# Model 1: Association between ASE and recombination class
# ==============================================================================
# Description:
#   Assess the association between the occurrence of allele-specific
#   expression (ASE) and recombination class across pig autosomes
#   (chromosomes 1-18) using binomial generalized linear mixed models.
#
# Model:
#   ASE.trait ~ Recombination + (1 | Sample)
#
# Comparisons:
#   1. HRR vs CS
#   2. Normal vs CS
#   3. CS vs combined HRR and Normal regions (non-CS)
#
# Input:
#   Tab-delimited file containing ASE.trait, Recombination and Sample.
#
# Output:
#   Fitted models, likelihood ratio tests, fixed-effect estimates,
#   odds ratios and 95% Wald confidence intervals.
#
# Usage:
#   Rscript scripts/01_ASE_recombination_GLMM.R
# ==============================================================================


# 1. Load required packages
# ------------------------------------------------------------------------------

library(lme4)


# 2. Define input and output paths
# ------------------------------------------------------------------------------

input_file <- file.path(
  "data",
  "chr1_18.all.sites_individual.recombination.combined"
)

output_dir <- file.path(
  "results",
  "01_ASE_recombination"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# 3. Load and prepare input data
# ------------------------------------------------------------------------------

ase_data <- read.table(
  input_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Check required input columns
required_columns <- c("ASE.trait", "Recombination", "Sample")

if (!all(required_columns %in% colnames(ase_data))) {
  stop("Required columns are missing from the input file.")
}

# Retain valid observations
ase_data <- subset(
  ase_data,
  ASE.trait %in% c("No.ASE", "ASE") &
    Recombination %in% c("CS", "Normal", "HRR") &
    !is.na(Sample)
)

# Define the binary response: No.ASE = 0, ASE = 1
ase_data$ASE.trait <- factor(
  ase_data$ASE.trait,
  levels = c("No.ASE", "ASE")
)

ase_data$Sample <- factor(ase_data$Sample)


# 4. Define a function for GLMM analysis
# ------------------------------------------------------------------------------

fit_glmm <- function(data, comparison, conf_level = 0.95) {

  data <- droplevels(data)

  message("Fitting model: ", comparison)

  # Full model: recombination class as a fixed effect
  full_model <- glmer(
    ASE.trait ~ Recombination + (1 | Sample),
    data = data,
    family = binomial(link = "logit")
  )

  # Null model: random intercept only
  null_model <- glmer(
    ASE.trait ~ 1 + (1 | Sample),
    data = data,
    family = binomial(link = "logit")
  )

  # Likelihood ratio test
  lrt <- anova(null_model, full_model, test = "Chisq")

  # Wald confidence intervals for fixed effects
  conf_int <- confint(
    full_model,
    parm = "beta_",
    level = conf_level,
    method = "Wald"
  )

  # Extract the recombination effect, excluding the intercept
  coef_table <- coef(summary(full_model))

  terms <- setdiff(rownames(coef_table), "(Intercept)")

  # Summarise fixed-effect estimates and odds ratios
  results <- data.frame(
    Comparison = comparison,
    Term = terms,
    Estimate = coef_table[terms, "Estimate"],
    SE = coef_table[terms, "Std. Error"],
    OR = exp(coef_table[terms, "Estimate"]),
    CI_lower = exp(conf_int[terms, 1]),
    CI_upper = exp(conf_int[terms, 2]),
    Confidence_level = conf_level,
    Wald_P = coef_table[terms, "Pr(>|z|)"],
    LRT_P = lrt[2, "Pr(>Chisq)"],
    N_observations = nobs(full_model),
    N_samples = length(unique(data$Sample)),
    Singular_fit = isSingular(full_model),
    row.names = NULL
  )

  # Return fitted models and statistical summaries
  list(
    full_model = full_model,
    null_model = null_model,
    model_summary = summary(full_model),
    confidence_intervals = conf_int,
    likelihood_ratio_test = lrt,
    results = results
  )
}


# 5. Prepare data for the three recombination comparisons
# ------------------------------------------------------------------------------

# 5.1 HRR vs CS
# Reference category: CS

data_HRR_CS <- subset(
  ase_data,
  Recombination %in% c("HRR", "CS")
)

data_HRR_CS$Recombination <- factor(
  data_HRR_CS$Recombination,
  levels = c("CS", "HRR")
)


# 5.2 Normal vs CS
# Reference category: CS

data_Normal_CS <- subset(
  ase_data,
  Recombination %in% c("Normal", "CS")
)

data_Normal_CS$Recombination <- factor(
  data_Normal_CS$Recombination,
  levels = c("CS", "Normal")
)


# 5.3 CS vs combined HRR and Normal regions
# Reference category: non-CS

data_CS_nonCS <- ase_data

data_CS_nonCS$Recombination <- ifelse(
  data_CS_nonCS$Recombination == "CS",
  "CS",
  "non-CS"
)

data_CS_nonCS$Recombination <- factor(
  data_CS_nonCS$Recombination,
  levels = c("non-CS", "CS")
)


# 6. Fit GLMMs for all comparisons
# ------------------------------------------------------------------------------

comparison_data <- list(
  HRR_vs_CS = data_HRR_CS,
  Normal_vs_CS = data_Normal_CS,
  CS_vs_nonCS = data_CS_nonCS
)

model_results <- lapply(
  names(comparison_data),
  function(comparison) {
    fit_glmm(
      data = comparison_data[[comparison]],
      comparison = comparison
    )
  }
)

names(model_results) <- names(comparison_data)


# 7. Summarise and export results
# ------------------------------------------------------------------------------

# Combine statistical results across comparisons
summary_table <- do.call(
  rbind,
  lapply(model_results, function(x) x$results)
)

rownames(summary_table) <- NULL

print(summary_table)

# Save fitted models and detailed statistical summaries
saveRDS(
  model_results,
  file = file.path(output_dir, "model_results.rds")
)

# Export fixed-effect estimates, odds ratios and confidence intervals
write.table(
  summary_table,
  file = file.path(output_dir, "model_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Record R and package versions for reproducibility
writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo.txt")
)

message("Model 1 analysis completed.")
