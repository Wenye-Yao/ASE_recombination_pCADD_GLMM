# ==============================================================================
# Model 3: Derived allele expression as a function of pCADD score
# ==============================================================================
# Description:
#   Investigate the relationship between predicted deleteriousness (pCADD score) and the direction of derived allele expression
#   using binomial generalized linear mixed models (GLMMs).
#
# Part 1:
#   Assess the association between pCADD score and the probabilityof derived allele overexpression.
#
# Part 2:
#   Assess the interaction between pCADD score and recombination class (HRR and CS).
#
# Response variable:
#   underexpress = 0 ; overexpress  = 1
#
# Random effect:
#   Sample identity
#
# Output:
#   Fitted models, AIC, likelihood ratio tests, fixed-effect estimates,
#   odds ratios, Wald confidence intervals, predicted probabilities
#   and figures.
#
# Usage:
#   Rscript scripts/03_Derived_expression_pCADD_GLMM.R
# ==============================================================================
# 1. Load required packages
# ------------------------------------------------------------------------------
library(lme4)
library(ggeffects)
library(ggplot2)
# 2. Define input and output paths
# ------------------------------------------------------------------------------

input_file_1 <- file.path("data","chr1_18.ASE.derved.allele.CADD.recombination.BT.data")

output_dir <- file.path("results","03_Derived_expression_pCADD")

output_part1 <- file.path(output_dir, "01_pCADD_effect")
output_part2 <- file.path(output_dir, "02_interaction_effect")

dir.create(output_part1, recursive = TRUE, showWarnings = FALSE)
dir.create(output_part2, recursive = TRUE, showWarnings = FALSE)

# Analysis settings
conf_level <- 0.95
options(contrasts = c("contr.treatment", "contr.poly"))
model_control <- glmerControl(optimizer = "bobyqa",optCtrl = list(maxfun = 2e5))

# 3. Define functions for data preparation and statistical analysis
# ------------------------------------------------------------------------------
# 3.1 Load and prepare input data
# ------------------------------------------------------------------------------
prepare_data <- function(input_file, interaction = FALSE) {
  data <- read.table(input_file,header = TRUE,sep = "\t",stringsAsFactors = FALSE)
  # Define required columns
  required_columns <- c("ASE.trait", "pCADD", "Sample")
  if (interaction) {required_columns <- c(required_columns, "Recombination")}
  if (!all(required_columns %in% colnames(data))) {stop("Required columns are missing from: ", input_file)}
  n_input <- nrow(data)
  
  # Convert pCADD scores to numeric values
  data$pCADD <- as.numeric(as.character(data$pCADD))
  
  # Retain valid ASE observations
  data <- subset(data, ASE.trait %in% c("underexpress", "overexpress") & is.finite(pCADD) & !is.na(Sample))

  # Define the binary response variable
  data$ASE.trait <- factor(data$ASE.trait, levels = c("underexpress", "overexpress"))
  data$Sample <- factor(data$Sample)
  # Retain HRR and CS observations for interaction analysis
  if (interaction) {
    data <- subset(data,Recombination %in% c("CS", "HRR"))
    # CS is the reference category
    data$Recombination <- factor(data$Recombination, levels = c("CS", "HRR"))}
  data <- droplevels(data)

  # Check whether sufficient categories are available
  if (nlevels(data$ASE.trait) != 2 ||
      nlevels(data$Sample) < 2 ||
      length(unique(data$pCADD)) < 2) {
    stop("Insufficient variation in the input data.")}

  if (interaction && nlevels(data$Recombination) != 2) {
    stop("Both CS and HRR observations are required.")}

  message("Input file: ", input_file)
  message("Original observations: ", n_input)
  message("Retained observations: ", nrow(data))
  message("Excluded observations: ", n_input - nrow(data))
  return(data)}

# 3.2 Fit a set of binomial GLMMs
# ------------------------------------------------------------------------------
fit_glmm_models <- function(formulas, data) {
  models <- lapply(names(formulas), function(model_name) {
      message("Fitting model: ", model_name)
      glmer(formula = formulas[[model_name]],data = data,family = binomial(link = "logit"),control = model_control)})
  names(models) <- names(formulas)
  return(models)}
# 3.3 Perform likelihood ratio tests between nested models
# ------------------------------------------------------------------------------
run_lrt <- function(models, comparisons) {
  lrt_results <- lapply(comparisons, function(comparison) {
      reduced_model <- comparison[1]
      full_model <- comparison[2]
      # Check that models contain the same observations
      if (nobs(models[[reduced_model]]) != nobs(models[[full_model]])) {
        stop("Models contain different numbers of observations.")}
      lrt <- anova(models[[reduced_model]],models[[full_model]],test = "Chisq")
      data.frame(
        Reduced_model = reduced_model,
        Full_model = full_model,
        Chi_square = lrt[2, "Chisq"],
        DF = lrt[2, "Chi Df"],
        LRT_P = lrt[2, "Pr(>Chisq)"],
        row.names = NULL )})
  results <- do.call(rbind, lrt_results)
  rownames(results) <- NULL
  return(results)}

# 3.4 Calculate model selection statistics
# ------------------------------------------------------------------------------
summarise_aic <- function(models) {
  aic_table <- data.frame(
    Model = names(models),
    AIC = vapply(models, AIC, numeric(1)),
    BIC = vapply(models, BIC, numeric(1)),
    LogLik = vapply(models,function(model) as.numeric(logLik(model)),numeric(1)),
    N_observations = vapply(models,function(model) as.numeric(nobs(model)),numeric(1)),
    row.names = NULL )

  # Difference from the model with the lowest AIC
  aic_table$Delta_AIC <- aic_table$AIC - min(aic_table$AIC)
  return(aic_table)}

# 3.5 Extract fixed-effect estimates and confidence intervals
# ------------------------------------------------------------------------------
summarise_fixed_effects <- function(models, conf_level = 0.95) {
  results <- lapply(names(models), function(model_name) {
      model <- models[[model_name]]
      coef_table <- coef(summary(model))
      # Calculate Wald confidence intervals
      critical_value <- qnorm((1 + conf_level) / 2)
      estimates <- coef_table[, "Estimate"]
      standard_errors <- coef_table[, "Std. Error"]
      ci_lower <- estimates - critical_value * standard_errors
      ci_upper <- estimates + critical_value * standard_errors
      data.frame(
        Model = model_name,
        Term = rownames(coef_table),
        Estimate = estimates,
        SE = standard_errors,
        OR = exp(estimates),
        CI_lower = exp(ci_lower),
        CI_upper = exp(ci_upper),
        Confidence_level = conf_level,
        Wald_P = coef_table[, "Pr(>|z|)"],
        N_observations = nobs(model),
        Singular_fit = if (inherits(model, "merMod")) {
          isSingular(model)
        } else {NA},
        row.names = NULL )})

  results <- do.call(rbind, results)
  rownames(results) <- NULL
  return(results)}

# 3.6 Export model objects and statistical summaries
# ------------------------------------------------------------------------------
export_results <- function(models, lrt_results,output_dir, conf_level = 0.95) {
  aic_table <- summarise_aic(models)
  fixed_effects <- summarise_fixed_effects(models,conf_level = conf_level)
  # Save fitted model objects
  saveRDS(models,file = file.path(output_dir, "model_results.rds"))
  
  # Save model selection statistics
  write.table(aic_table, file = file.path(output_dir, "model_selection.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  # Save likelihood ratio test results
  write.table( lrt_results, file = file.path(output_dir, "likelihood_ratio_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  # Save fixed-effect estimates
  write.table(fixed_effects,file = file.path(output_dir, "model_summary.tsv"),sep = "\t",quote = FALSE,row.names = FALSE)

  print(aic_table)
  print(lrt_results)
  print(fixed_effects)}

# 4. Part 1: Derived allele expression as a function of pCADD
# ==============================================================================
# 4.1 Load and prepare data
# ------------------------------------------------------------------------------

cadd_data <- prepare_data( input_file = input_file_1, interaction = FALSE)

# 4.2 Fit candidate models
# ------------------------------------------------------------------------------
C0 <- glm(ASE.trait ~ 1, data = cadd_data, family = binomial(link = "logit") ) # C0: Intercept-only GLM
# Define mixed-effects model formulas
cadd_formulas <- list(
  C1 = ASE.trait ~ 1 + (1 | Sample),  # C1: Intercept-only GLMM
  C2 = ASE.trait ~ pCADD + (1 | Sample) # C2: pCADD as a fixed effect
)
# Fit mixed-effects models
cadd_models <- fit_glmm_models(formulas = cadd_formulas, data = cadd_data)
# Combine GLM and GLMM objects
cadd_models <- c(list(C0 = C0),cadd_models)
# 4.3 Compare models
# ------------------------------------------------------------------------------
# Test whether pCADD improves model fit
cadd_comparisons <- list(c("C1", "C2"))
cadd_lrt <- run_lrt(models = cadd_models,comparisons = cadd_comparisons)
                    
# 4.4 Export statistical results
# ------------------------------------------------------------------------------
export_results( models = cadd_models, lrt_results = cadd_lrt, output_dir = output_part1, conf_level = conf_level)

# 4.5 Predict the probability of derived allele overexpression
# ------------------------------------------------------------------------------
# Predictions are based on the fixed effects of model C2, with sample-level random effects set to zero.
cadd_predictions <- as.data.frame(
  ggpredict(cadd_models$C2, terms = "pCADD", type = "fixed", ci_level = conf_level))

write.table( cadd_predictions, file = file.path(output_part1, "predictions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# 4.6 Plot predicted probabilities
# ------------------------------------------------------------------------------
plot_cadd <- ggplot(cadd_predictions,aes(x = x, y = predicted)) + 
                    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "blue", alpha = 0.2) +
                    geom_line(colour = "blue", linewidth = 1) +
                    labs(x = "pCADD score", y = "Predicted probability of\n derived allele overexpression") +
                    theme_bw(base_size = 14) +
                    theme(axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"))

print(plot_cadd)

ggsave(filename = file.path(output_part1,"pCADD_effect.pdf"),plot = plot_cadd,width = 6,height = 5,units = "in")

# 5. Part 2: Interaction between pCADD and recombination class
# ==============================================================================
# 5.1 Load and prepare HRR and CS data
# ------------------------------------------------------------------------------
interaction_data <- prepare_data(input_file = input_file_1,interaction = TRUE)

# 5.2 Fit candidate models
# ------------------------------------------------------------------------------
interaction_formulas <- list(
  M0 = ASE.trait ~ 1 + (1 | Sample), # M0: Intercept-only GLMM
  M1 = ASE.trait ~ pCADD + (1 | Sample), # M1: pCADD effect
  M2 = ASE.trait ~ Recombination + (1 | Sample), # M2: Recombination effect
  M3 = ASE.trait ~ pCADD + Recombination + (1 | Sample), # M3: Additive effects of pCADD and recombination
  M4 = ASE.trait ~ pCADD * Recombination + (1 | Sample) # M4: Interaction between pCADD and recombination )

interaction_models <- fit_glmm_models( formulas = interaction_formulas, data = interaction_data)

# 5.3 Compare models
# ------------------------------------------------------------------------------
interaction_comparisons <- list(
  # Test single-predictor models against the baseline
  c("M0", "M1"),
  c("M0", "M2"),
  c("M0", "M3"),
  c("M1", "M3"),
  c("M2", "M3"),
  # Compare the interaction model with reduced models
  c("M0", "M4"),
  c("M1", "M4"),
  c("M2", "M4"),
  # Specifically test the interaction effect
  c("M3", "M4")
)

interaction_lrt <- run_lrt( models = interaction_models, comparisons = interaction_comparisons)

# 5.4 Export statistical results
# ------------------------------------------------------------------------------
export_results(
  models = interaction_models,
  lrt_results = interaction_lrt,
  output_dir = output_part2,
  conf_level = conf_level
)
# 5.5 Predict probabilities across recombination classes
# ------------------------------------------------------------------------------
interaction_predictions <- as.data.frame(
  ggpredict(interaction_models$M4, terms = c("pCADD", "Recombination"), type = "fixed", ci_level = conf_level ))

write.table(interaction_predictions, file = file.path(output_part2, "predictions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# 5.6 Plot predicted probabilities
# ------------------------------------------------------------------------------
# Define colours explicitly to ensure consistent plotting
recombination_colours <- c("CS" = "black", "HRR" = "red")

plot_interaction <- ggplot( interaction_predictions, aes(x = x, y = predicted, colour = group, fill = group)) +
  geom_ribbon( aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = recombination_colours) +
  scale_fill_manual(values = recombination_colours) +
  labs(x = "pCADD score", y = "Predicted probability of\n derived allele overexpression") +
  theme_bw(base_size = 14) +
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"), legend.title = element_blank())
  
print(plot_interaction)
ggsave( filename = file.path( output_part2, "pCADD_recombination_interaction.pdf" ), plot = plot_interaction, width = 6, height = 5, units = "in" )

# 6. Record software versions
# ------------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), con = file.path(output_dir, "sessionInfo.txt"))

message("Model 3 analysis completed.")
