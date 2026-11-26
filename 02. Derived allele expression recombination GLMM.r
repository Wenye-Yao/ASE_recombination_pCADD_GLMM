# ==============================================================================
# Model 2: Recombination class as a function of derived allele expression
# ==============================================================================
# Description:
#   Assess the association between derived allele expression classes
#   and local recombination landscapes using binomial generalized
#   linear mixed models (GLMMs).
#
# Model:
#   Recombination ~ ASE.trait + (1 | Sample)
#
# Recombination comparisons:
#   1. HRR vs CS
#   2. Normal vs CS
#   3. HRR/Normal vs CS
#
# Derived allele expression comparisons:
#   1. Overexpression vs No ASE
#   2. Underexpression vs No ASE
#   3. No ASE vs ASE
#
# Input:
#   Individual-level data containing Recombination, ASE.trait and Sample.
#   Optional columns include pCADD, Tissue and functional annotation.
#
# Output:
#   Fitted models, likelihood ratio tests, fixed-effect estimates,
#   odds ratios and Wald confidence intervals.
#
# Dependencies:
#   R package: lme4
# ==============================================================================
# 1. Load required packages
# ------------------------------------------------------------------------------
library(lme4)
# 2. Define input and output paths
# ------------------------------------------------------------------------------
output_dir <- file.path( "results","02_Derived_expression_recombination")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# 3. Define a function for fitting individual GLMMs
# ------------------------------------------------------------------------------
fit_glmm <- function(data,
                     recombination_comparison,
                     expression_comparison,
                     category,
                     conf_level = 0.95) {
  data <- droplevels(data)
  message( "Fitting model: ", recombination_comparison, " | ", expression_comparison, " | ", category)

  # Check that both response and predictor contain two categories
  if (nlevels(data$Recombination) != 2 ||
      nlevels(data$ASE.trait) != 2) {
    stop("Both response and predictor must contain two categories.")}

  # Full model
  full_model <- glmer( Recombination ~ ASE.trait + (1 | Sample),data = data, family = binomial(link = "logit"))

  # Null model
  null_model <- glmer(Recombination ~ 1 + (1 | Sample),data = data,family = binomial(link = "logit"))

  # Likelihood ratio test
  lrt <- anova(null_model,full_model,test = "Chisq")

  # Wald confidence intervals for fixed effects
  conf_int <- confint(full_model,parm = "beta_",method = "Wald",level = conf_level)

  # Extract fixed-effect estimates, excluding the intercept
  coef_table <- coef(summary(full_model))

  terms <- setdiff(
    rownames(coef_table),
    "(Intercept)"
  )

  # Extract statistical results
  results <- data.frame(
    Category = category,
    Recombination_comparison = recombination_comparison,
    Expression_comparison = expression_comparison,
    Response_reference = levels(data$Recombination)[1],
    Response_event = levels(data$Recombination)[2],
    Predictor_reference = levels(data$ASE.trait)[1],
    Predictor_comparison = levels(data$ASE.trait)[2],
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


# 4. Define a function to prepare recombination comparisons
# ------------------------------------------------------------------------------

prepare_recombination <- function(data, comparison) {

  if (comparison == "HRR_vs_CS") {

    data <- subset(
      data,
      Recombination %in% c("CS", "HRR")
    )

    data$Recombination <- factor(
      data$Recombination,
      levels = c("CS", "HRR")
    )

  } else if (comparison == "Normal_vs_CS") {

    data <- subset(
      data,
      Recombination %in% c("CS", "Normal")
    )

    data$Recombination <- factor(
      data$Recombination,
      levels = c("CS", "Normal")
    )

  } else if (comparison == "HRR_Normal_vs_CS") {

    data <- subset(
      data,
      Recombination %in% c("CS", "HRR", "Normal")
    )

    # Combine HRR and Normal regions
    data$Recombination <- ifelse(
      data$Recombination == "CS",
      "CS",
      "HRR/Normal"
    )

    data$Recombination <- factor(
      data$Recombination,
      levels = c("CS", "HRR/Normal")
    )

  } else {

    stop("Invalid recombination comparison: ", comparison)

  }

  droplevels(data)
}


# 5. Define a function to prepare expression comparisons
# ------------------------------------------------------------------------------

prepare_expression <- function(data, comparison) {

  if (comparison == "Over_vs_NoASE") {

    data <- subset(
      data,
      ASE.trait %in% c("No.ASE", "overexpress")
    )

    data$ASE.trait <- factor(
      data$ASE.trait,
      levels = c("No.ASE", "overexpress")
    )

  } else if (comparison == "Under_vs_NoASE") {

    data <- subset(
      data,
      ASE.trait %in% c("No.ASE", "underexpress")
    )

    data$ASE.trait <- factor(
      data$ASE.trait,
      levels = c("No.ASE", "underexpress")
    )

  } else if (comparison == "NoASE_vs_ASE") {

    data <- subset(
      data,
      ASE.trait %in% c(
        "No.ASE",
        "overexpress",
        "underexpress"
      )
    )

    # Combine overexpression and underexpression into ASE
    data$ASE.trait <- ifelse(
      data$ASE.trait == "No.ASE",
      "No.ASE",
      "ASE"
    )

    # Preserve the original comparison direction:
    # ASE = 0; No.ASE = 1
    data$ASE.trait <- factor(
      data$ASE.trait,
      levels = c("ASE", "No.ASE")
    )

  } else {

    stop("Invalid expression comparison: ", comparison)

  }

  droplevels(data)
}


# 6. Define the main analysis function
# ------------------------------------------------------------------------------

run_models <- function(
    input_data,
    category,
    output_dir,
    recombination_comparisons = c(
      "HRR_vs_CS",
      "Normal_vs_CS",
      "HRR_Normal_vs_CS"
    ),
    conf_level = 0.95) {

  # Standardise the recombination column name
  if ("recombination" %in% colnames(input_data) &&
      !"Recombination" %in% colnames(input_data)) {

    colnames(input_data)[
      colnames(input_data) == "recombination"
    ] <- "Recombination"
  }

  # Check required input columns
  required_columns <- c(
    "Recombination",
    "ASE.trait",
    "Sample"
  )

  if (!all(required_columns %in% colnames(input_data))) {
    stop("Required columns are missing from the input data.")
  }

  # Check requested recombination comparisons
  valid_comparisons <- c(
    "HRR_vs_CS",
    "Normal_vs_CS",
    "HRR_Normal_vs_CS"
  )

  if (!all(recombination_comparisons %in% valid_comparisons)) {
    stop("Invalid recombination comparison specified.")
  }

  # Retain valid observations
  n_input <- nrow(input_data)

  input_data <- subset(
    input_data,
    Recombination %in% c("CS", "HRR", "Normal") &
      ASE.trait %in% c(
        "No.ASE",
        "overexpress",
        "underexpress"
      ) &
      !is.na(Sample)
  )

  message(
    "Excluded observations: ",
    n_input - nrow(input_data)
  )

  input_data$Sample <- factor(input_data$Sample)

  # Define expression comparisons
  expression_comparisons <- c(
    "Over_vs_NoASE",
    "Under_vs_NoASE",
    "NoASE_vs_ASE"
  )

  # Initialise lists for model objects and statistical results
  model_results <- list()
  result_tables <- list()

  # Fit models for each recombination comparison
  for (rec_comp in recombination_comparisons) {

    rec_data <- prepare_recombination(
      data = input_data,
      comparison = rec_comp
    )

    # Fit models for each expression comparison
    for (expr_comp in expression_comparisons) {

      model_data <- prepare_expression(
        data = rec_data,
        comparison = expr_comp
      )

      model_name <- paste(
        rec_comp,
        expr_comp,
        sep = "__"
      )

      # Check whether sufficient categories are available
      if (nrow(model_data) == 0 ||
          nlevels(model_data$Recombination) != 2 ||
          nlevels(model_data$ASE.trait) != 2 ||
          nlevels(model_data$Sample) < 2) {

        warning(
          "Skipping model due to insufficient data: ",
          model_name
        )

        next
      }

      # Fit full and null models
      result <- fit_glmm(
        data = model_data,
        recombination_comparison = rec_comp,
        expression_comparison = expr_comp,
        category = category,
        conf_level = conf_level
      )

      model_results[[model_name]] <- result

      result_tables[[model_name]] <- result$results
    }
  }

  # Check whether any models were successfully fitted
  if (length(result_tables) == 0) {
    stop("No valid models were fitted.")
  }

  # Combine results into a single table
  summary_table <- do.call(
    rbind,
    result_tables
  )

  rownames(summary_table) <- NULL

  # Create a separate output directory for each analysis category
  category_dir <- file.path(output_dir, category)

  dir.create(
    category_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # Save fitted models and detailed statistical summaries
  saveRDS(
    model_results,
    file = file.path(
      category_dir,
      "model_results.rds"
    )
  )

  # Export statistical results for downstream analysis and plotting
  write.table(
    summary_table,
    file = file.path(
      category_dir,
      "model_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # Record R and package versions
  writeLines(
    capture.output(sessionInfo()),
    con = file.path(
      category_dir,
      "sessionInfo.txt"
    )
  )

  message("Analysis completed: ", category)

  # Return model objects and statistical results
  invisible(
    list(
      models = model_results,
      summary = summary_table
    )
  )
}
