# Overview

This repository contains the R scripts used to investigate the association between recombination landscapes, allele-specific expression (ASE), and predicted deleteriousness in pig genomes.

The analyses integrate individual-level ASE observations, derived allele information, pCADD scores, and fine-scale recombination classes across pig autosomes (chromosomes 1–18).

Three statistical analyses are performed using generalized linear mixed models (GLMMs), with sample identity included as a random intercept to account for repeated observations within samples.
# Statistical analysis
## Model 1 : ASE as a function of recombination class
**Script** : 01_ASE_recombination_GLMM.r

This script investigates the association between the occurrence of ASE and local recombination class across pig autosomes. Binomial GLMMs are fitted with ASE status as the response variable, recombination class as a fixed effect, and sample identity as a random intercept.

Three comparisons are performed:
1. HRR versus CS.
2. Normal recombination regions versus CS.
3. CS versus combined HRR and normal recombination regions.

Odds ratios, Wald confidence intervals, and likelihood ratio tests are used to assess the association between ASE and recombination class.

This code is corresponding with the Fig. 1A

## Model 2: Recombination class as a function of derived allele expression class
**Script** : 02.Derived allele expression recombination GLMM.r

This script investigates the association between derived allele expression classes (overexpresion, underexpression and No.ASE) and local recombination landscapes (HRR, CS).
Binomial GLMMs are fitted with recombination class as the response variable, derived allele expression class as a fixed effect, and sample identity as a random intercept.

Three expression contrasts are considered:

1. Overexpression versus No ASE.
2. Underexpression versus No ASE.
3. No ASE versus combined overexpression and underexpression.

The analysis can be applied to different recombination comparisons and subsets of variants, including those defined by pCADD scores or functional annotations.

This code is corresponding with the Fig. 1B-D & Fig.S2A-F & Fig.3 & Fig.S5 C-E

## Model 3: Derived allele expression as a function of pCADD score
**script**: 03.Derived allele pCADD Recombination.r

This script investigates the relationship between predicted deleteriousness and the direction of derived allele expression.

The analysis consists of two parts:

1. Assessing the association between pCADD score and the probability of derived allele overexpression.

2. Assessing the interaction between pCADD score and recombination class (HRR and CS).

Candidate models are compared using the Akaike information criterion (AIC) and likelihood ratio tests. Model-based predictions and confidence intervals are generated to visualise the relationship between pCADD score and derived allele expression.

This code is corresponding with the Fig.2 & Fig.S3 & Fig.S4 
