# =============================================================================
# Arabidopsis Fruit Production Response to Nutrients and Damage
# STAT 642 - Spring 2026 - Team 11
#
# R-script counterpart to slides_condensed.qmd. Contains all code and produces
# all outputs (tables, plots, tests) from the slide deck, minus the prose.
# =============================================================================


# -----------------------------------------------------------------------------
# Setup: libraries, data, and global options
# -----------------------------------------------------------------------------

source('R/load_arabidopsis.R')

library(dplyr)
library(ggplot2)
library(MASS)
library(car)
library(ggplot2)
library(lme4)
library(lmerTest)
library(knitr)
library(glmmTMB)
library(DHARMa)
library(emmeans)

# Sum-to-zero contrasts so Type III ANOVA main-effect tests are
# averaged over factor levels rather than evaluated at a reference level.
options(contrasts = c("contr.sum", "contr.poly"))

# Coerce categorical predictors to factors up front so every model
# (Box-Cox lm, full LMM, NB GLMM) uses the same parameterization.
Arabidopsis$nutrient     <- factor(Arabidopsis$nutrient)
Arabidopsis$rack         <- factor(Arabidopsis$rack)
Arabidopsis$logp1_fruits <- log(Arabidopsis$total.fruits + 1)


# -----------------------------------------------------------------------------
# Slide 1: Distribution of total fruits (raw histogram)
# -----------------------------------------------------------------------------

print(
  ggplot(Arabidopsis, aes(x = total.fruits)) +
    geom_histogram(binwidth = 5, fill = "lightblue", color = "black") +
    labs(title = "Distribution of Total Fruits",
         x = "Total Fruits", y = "Count") +
    theme_minimal()
)


# -----------------------------------------------------------------------------
# Slide 3: Variables & experimental design table
# -----------------------------------------------------------------------------

vars_tbl <- data.frame(
  Variable = c("Total Fruits", "Region", "Nutrients", "AMD",
               "Status", "Rack", "Population", "Genotype"),
  Role = c("Response", "Fixed", "Fixed", "Fixed",
           "Fixed", "Fixed (additive)", "Random", "Random"),
  Levels = c("count per plant",
             "NL, SP, SW",
             "1 prill (low), 8 prills (high)",
             "clipped, unclipped",
             "Normal, Petri.Plate, Transplant",
             "1, 2",
             "9 levels, nested in Region",
             "24 seed families, nested in Population")
)
print(
  kable(vars_tbl, align = "lll",
        col.names = c("Variable", "Role", "Levels / Description"))
)


# -----------------------------------------------------------------------------
# Slide 5: Response transformation - Box-Cox on initial LM
# -----------------------------------------------------------------------------

# Shift by +1 before Box-Cox so zeros are admissible.
fit_init <- lm(
  (total.fruits + 1) ~ reg * nutrient * amd +
    reg * status + nutrient * status + amd * status + rack,
  data = Arabidopsis
)

bc <- boxcox(fit_init, lambda = seq(0, 1, by = 0.01))
optimal_lambda <- bc$x[which.max(bc$y)]
abline(v = optimal_lambda, lty = 2)


# -----------------------------------------------------------------------------
# Slide 6: Full LMM and backward elimination via LRT
# -----------------------------------------------------------------------------

# Full LMM (REML): treatment factors fully crossed, plus nuisance interactions
# with Status, Rack additive, and (popu), (gen) random intercepts.
fit_full <- lmer(
  logp1_fruits ~ reg * nutrient * amd +
    reg * status + nutrient * status + amd * status +
    rack + (1 | popu) + (1 | gen),
  data = Arabidopsis
)

# Step 1: drop Nutrients x Status
fit_no_nut_status <- lmer(
  logp1_fruits ~ reg * nutrient * amd +
    reg * status + amd * status +
    rack + (1 | popu) + (1 | gen),
  data = Arabidopsis
)
fit_full_ml <- update(fit_full, REML = FALSE)
fit_no_nut_status_ml <- update(fit_no_nut_status, REML = FALSE)
lrt1 <- anova(fit_no_nut_status_ml, fit_full_ml)

# Step 2: drop AMD x Status
fit_no_amd_status <- lmer(
  logp1_fruits ~ reg * nutrient * amd +
    reg * status +
    rack + (1 | popu) + (1 | gen),
  data = Arabidopsis
)
fit_no_amd_status_ml <- update(fit_no_amd_status, REML = FALSE)
lrt2 <- anova(fit_no_nut_status_ml, fit_no_amd_status_ml)

# Step 3: drop three-way Region x Nutrients x AMD
fit_no_three_way <- lmer(
  logp1_fruits ~ reg * nutrient + nutrient * amd + reg * amd +
    reg * status + rack + (1 | popu) + (1 | gen),
  data = Arabidopsis
)
fit_no_three_way_ml <- update(fit_no_three_way, REML = FALSE)
lrt3 <- anova(fit_no_three_way_ml, fit_no_amd_status_ml)

# Step 4: drop Region x AMD
fit_no_reg_amd <- lmer(
  logp1_fruits ~ reg * nutrient + nutrient * amd +
    reg * status + rack + (1 | popu) + (1 | gen),
  data = Arabidopsis
)
fit_no_reg_amd_ml <- update(fit_no_reg_amd, REML = FALSE)
lrt4 <- anova(fit_no_three_way_ml, fit_no_reg_amd_ml)

# Step 5: drop Nutrients x AMD -> leaves the reduced model.
fit_no_nut_amd <- lmer(
  logp1_fruits ~ reg * nutrient + amd +
    reg * status + rack + (1 | popu) + (1 | gen),
  data = Arabidopsis
)
fit_no_nut_amd_ml <- update(fit_no_nut_amd, REML = FALSE)
lrt5 <- anova(fit_no_reg_amd_ml, fit_no_nut_amd_ml)

# Reduced (final) model and its Type III F-tests (REML).
fit_reduced <- fit_no_nut_amd
anova_reduced <- anova(fit_reduced, type = 3)

# LRT summary table shown on slide 6.
lrts <- list(lrt1, lrt2, lrt3, lrt4, lrt5)
lrt_summary <- data.frame(
  Step = seq_along(lrts),
  Test = c(
    "Nutrients \u00d7 Status",
    "AMD \u00d7 Status",
    "Region \u00d7 Nutrients \u00d7 AMD",
    "Region \u00d7 AMD",
    "Nutrients \u00d7 AMD"
  ),
  Chisq   = round(sapply(lrts, function(x) x$Chisq[2]), 2),
  Df      = sapply(lrts, function(x) x$Df[2]),
  P_Value = signif(sapply(lrts, function(x) x$`Pr(>Chisq)`[2]), 3)
)
print(
  kable(
    lrt_summary,
    align = "clccc",
    col.names = c("Step", "Interaction dropped", "Chisq", "Df", "p"),
    row.names = FALSE
  )
)


# -----------------------------------------------------------------------------
# Slide 7: Diagnostics for the reduced LMM
# -----------------------------------------------------------------------------

# Shared lattice padding used in both diagnostic plots.
tight_pad <- list(
  par.xlab.text = list(cex = 0.7),
  par.ylab.text = list(cex = 0.7),
  par.main.text = list(cex = 0.8),
  layout.heights = list(
    top.padding        = 0.2,
    main.key.padding   = 0,
    key.axis.padding   = 0,
    axis.xlab.padding  = 0.2,
    xlab.key.padding   = 0,
    key.sub.padding    = 0,
    bottom.padding     = 0.2
  ),
  layout.widths = list(
    left.padding       = 0.2,
    key.ylab.padding   = 0,
    ylab.axis.padding  = 0.2,
    axis.key.padding   = 0,
    right.padding      = 0.2
  )
)

# Residuals vs. fitted (with smoother).
print(
  plot(fit_reduced,
       type = c("p", "smooth"),
       col.line = 1,
       main = list("Residuals vs Fitted", cex = 0.8),
       scales = list(cex = 0.6),
       par.settings = tight_pad)
)

# Normal Q-Q of residuals.
print(
  lattice::qqmath(fit_reduced,
    main = list("Normal Q-Q", cex = 0.8),
    scales = list(cex = 0.6),
    par.settings = tight_pad)
)

# Levene (Brown-Forsythe) test of residual variance across each factor.
lev_results <- data.frame(
  Factor = c("Region", "Nutrients", "Status", "AMD"),
  p_value = c(
    leveneTest(resid(fit_reduced) ~ Arabidopsis$reg)[1,3],
    leveneTest(resid(fit_reduced) ~ Arabidopsis$nutrient)[1,3],
    leveneTest(resid(fit_reduced) ~ Arabidopsis$status)[1,3],
    leveneTest(resid(fit_reduced) ~ Arabidopsis$amd)[1,3]
  )
)
lev_results$p_value <- signif(lev_results$p_value, 3)
print(kable(lev_results, align = "lc"))

# Shapiro-Wilk on scaled residuals (used inline on slide 7).
shapiro_p <- signif(
  shapiro.test(resid(fit_reduced, scaled = TRUE))$p.value, 2
)
cat("Shapiro-Wilk p-value (scaled residuals):", shapiro_p, "\n")


# -----------------------------------------------------------------------------
# Slide 8: Type III F-tests for the reduced model
# -----------------------------------------------------------------------------

tbl <- as.data.frame(anova_reduced)
tbl$Effect <- rownames(tbl)
tbl <- tbl[, c("Effect", "F value", "NumDF", "DenDF", "Pr(>F)")]
tbl$`F value` <- round(tbl$`F value`, 3)
tbl$DenDF    <- round(tbl$DenDF, 2)
tbl$`Pr(>F)` <- ifelse(
  tbl$`Pr(>F)` < 0.001,
  "<0.001",
  formatC(tbl$`Pr(>F)`, format = "f", digits = 3)
)
print(
  kable(
    tbl, row.names = FALSE,
    col.names = c("Effect", "F", "df1", "df2", "p-value"),
    align = "lcccc",
    caption = "Type III F-tests for log(Total Fruits + 1)"
  )
)


# -----------------------------------------------------------------------------
# Slide 9: Variance components + LRT for random effects (RQ3)
# -----------------------------------------------------------------------------

vc <- as.data.frame(VarCorr(fit_reduced))
vc <- vc[, c("grp", "vcov")]
vc$pct  <- 100 * vc$vcov / sum(vc$vcov)
vc$vcov <- round(vc$vcov, 3)
vc$pct  <- round(vc$pct, 1)
names(vc) <- c("Component", "Variance", "% of total")
print(
  kable(
    vc, row.names = FALSE,
    caption = "Variance components (reduced LMM)"
  )
)

# Population share referenced in the narrative / conclusions.
popu_pct <- vc$`% of total`[vc$Component == "popu"]
cat("Population share of total variance:", popu_pct, "%\n")

# LRT for each random effect via lmerTest::ranova.
rt <- as.data.frame(lmerTest::ranova(fit_reduced))
rt$Term <- rownames(rt)
rt <- rt[, c("Term", "LRT", "Df", "Pr(>Chisq)")]
rt$LRT <- round(rt$LRT, 2)
rt$`Pr(>Chisq)` <- ifelse(
  rt$`Pr(>Chisq)` < 0.001,
  "<0.001",
  formatC(rt$`Pr(>Chisq)`, format = "f", digits = 3)
)
print(
  kable(
    rt, row.names = FALSE,
    col.names = c("Term", "LRT", "Df", "p-value"),
    caption = "LRT for random effects"
  )
)


# -----------------------------------------------------------------------------
# Slide 10: Post-hoc EMMs (pairwise contrasts on the log scale)
# -----------------------------------------------------------------------------

# Nutrients within Region.
emm_nutrient <- emmeans(fit_reduced, pairwise ~ nutrient | reg)
print(
  kable(
    as.data.frame(emm_nutrient$contrasts)[,
      c("reg", "contrast", "estimate", "SE", "p.value")],
    digits = 3, row.names = FALSE,
    col.names = c("Reg", "Contrast", "Est.", "SE", "p")
  )
)

# AMD main effect.
emm_amd <- emmeans(fit_reduced, pairwise ~ amd)
print(
  kable(
    as.data.frame(emm_amd$contrasts)[,
      c("contrast", "estimate", "SE", "p.value")],
    digits = 3, row.names = FALSE,
    col.names = c("Contrast", "Est.", "SE", "p")
  )
)

# Status within Region.
emm_status <- emmeans(fit_reduced, pairwise ~ status | reg)
print(
  kable(
    as.data.frame(emm_status$contrasts)[,
      c("reg", "contrast", "estimate", "SE", "p.value")],
    digits = 3, row.names = FALSE,
    col.names = c("Reg", "Contrast", "Est.", "SE", "p")
  )
)


# -----------------------------------------------------------------------------
# Slide 11: Interaction plots
# -----------------------------------------------------------------------------

# Region x Nutrients: cell means with +/- 1 SE error bars.
print(
  Arabidopsis %>%
    group_by(reg, nutrient) %>%
    summarise(
      mean_fruits = mean(logp1_fruits),
      se = sd(logp1_fruits) / sqrt(n()),
      .groups = "drop"
    ) %>%
    ggplot(aes(x = factor(nutrient,
                          labels = c("Low (1 prill)", "High (8 prills)")),
               y = mean_fruits,
               color = reg,
               group = reg)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_fruits - se, ymax = mean_fruits + se),
                  width = 0.1, linewidth = 0.8) +
    labs(x = "Nutrient Level",
         y = "Mean log(Total Fruits + 1)",
         color = "Region") +
    theme_minimal() +
    theme(legend.position = "bottom",
          text = element_text(size = 10))
)

# Region x Status.
print(
  Arabidopsis %>%
    group_by(reg, status) %>%
    summarise(
      mean_fruits = mean(logp1_fruits),
      se = sd(logp1_fruits) / sqrt(n()),
      .groups = "drop"
    ) %>%
    ggplot(aes(x = status,
               y = mean_fruits,
               color = reg,
               group = reg)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_fruits - se, ymax = mean_fruits + se),
                  width = 0.15, linewidth = 0.8) +
    labs(x = "Plant Status",
         y = "Mean log(Total Fruits + 1)",
         color = "Region") +
    theme_minimal() +
    theme(legend.position = "bottom",
          text = element_text(size = 10),
          axis.text.x = element_text(angle = 15, hjust = 1))
)


# -----------------------------------------------------------------------------
# Slide 12: Robustness check - negative binomial GLMM on raw counts
# -----------------------------------------------------------------------------

fit_nb <- glmmTMB(
  total.fruits ~ reg * nutrient + reg * status + amd + rack +
    (1 | popu) + (1 | gen),
  data = Arabidopsis,
  family = nbinom2
)

nb_anova <- as.data.frame(car::Anova(fit_nb, type = 3))
nb_anova$Effect <- rownames(nb_anova)
rownames(nb_anova) <- NULL

nb_tbl <- nb_anova[, c("Effect", "Chisq", "Df", "Pr(>Chisq)")]
names(nb_tbl) <- c("Effect", "Chisq", "Df", "p-value")
nb_tbl$Chisq <- round(nb_tbl$Chisq, 2)
nb_tbl$`p-value` <- ifelse(
  nb_tbl$`p-value` < 0.0001,
  "<0.0001",
  formatC(nb_tbl$`p-value`, format = "f", digits = 4)
)
print(
  kable(nb_tbl, align = "lccc",
        caption = "Type III tests: NB GLMM")
)


# -----------------------------------------------------------------------------
# Backup B1: Type III F-tests for the full (pre-reduction) model
# -----------------------------------------------------------------------------

anova_full <- anova(fit_full, type = 3)
full_tbl <- as.data.frame(anova_full)
full_tbl$Effect <- rownames(full_tbl)
full_tbl <- full_tbl[, c("Effect", "F value", "NumDF", "DenDF", "Pr(>F)")]
full_tbl$`F value` <- round(full_tbl$`F value`, 3)
full_tbl$DenDF    <- round(full_tbl$DenDF, 2)
full_tbl$`Pr(>F)` <- ifelse(
  full_tbl$`Pr(>F)` < 0.001,
  "<0.001",
  formatC(full_tbl$`Pr(>F)`, format = "f", digits = 3)
)
print(
  kable(
    full_tbl, row.names = FALSE,
    col.names = c("Effect", "F", "df1", "df2", "p-value"),
    align = "lcccc",
    caption = "Type III F-tests for the full model"
  )
)


# -----------------------------------------------------------------------------
# Backup B2: AIC/BIC all-subsets search over candidate interactions
# -----------------------------------------------------------------------------

# Main-effects baseline (all main effects + rack + random effects).
base <- "logp1_fruits ~ reg + nutrient + amd + status + rack + (1 | popu) + (1 | gen)"

treatment_ints <- list(
  "Reg x Nut"       = "reg:nutrient",
  "Reg x AMD"       = "reg:amd",
  "Nut x AMD"       = "nutrient:amd",
  "Reg x Nut x AMD" = "reg:nutrient:amd"
)
nuisance_ints <- list(
  "Reg x Stat" = "reg:status",
  "Nut x Stat" = "nutrient:status",
  "AMD x Stat" = "amd:status"
)

# Hierarchy: if 3-way is in, all three 2-way must also be in.
t_combos <- do.call(expand.grid, replicate(4, c(FALSE, TRUE), simplify = FALSE))
names(t_combos) <- names(treatment_ints)
t_valid <- with(t_combos,
  !`Reg x Nut x AMD` | (`Reg x Nut` & `Reg x AMD` & `Nut x AMD`)
)
t_combos <- t_combos[t_valid, , drop = FALSE]

n_combos <- do.call(expand.grid, replicate(3, c(FALSE, TRUE), simplify = FALSE))
names(n_combos) <- names(nuisance_ints)

# Enumerate and fit all (treatment, nuisance) subset combinations with ML.
models <- list()
for (i in seq_len(nrow(t_combos))) {
  for (j in seq_len(nrow(n_combos))) {
    t_in <- names(treatment_ints)[as.logical(t_combos[i, ])]
    n_in <- names(nuisance_ints)[as.logical(n_combos[j, ])]
    included <- c(t_in, n_in)
    if (length(included) == 0) {
      label <- "Main effects only"
      form <- base
    } else {
      label <- paste(included, collapse = " + ")
      terms <- c(treatment_ints[t_in], nuisance_ints[n_in])
      form <- paste(base, paste(terms, collapse = " + "), sep = " + ")
    }
    models[[label]] <- lmer(as.formula(form), data = Arabidopsis, REML = FALSE)
  }
}

aic_bic <- data.frame(
  Model = names(models),
  AIC   = sapply(models, AIC),
  BIC   = sapply(models, BIC)
)
aic_bic$dAIC <- aic_bic$AIC - min(aic_bic$AIC)
aic_bic$dBIC <- aic_bic$BIC - min(aic_bic$BIC)
aic_bic <- aic_bic[order(aic_bic$AIC), ]

print(
  kable(
    head(aic_bic, 10), digits = 2, row.names = FALSE,
    caption = "Top 10 of 72 models, ranked by AIC"
  )
)


# -----------------------------------------------------------------------------
# Backup B3: Scale-location plot
# -----------------------------------------------------------------------------

scale_loc_plot <- function(m) {
  plot(m, sqrt(abs(resid(.))) ~ fitted(.),
       type = c("p", "smooth"),
       xlab = "Fitted",
       ylab = "Sqrt(|Residuals|)",
       main = "")
}
print(scale_loc_plot(fit_reduced))


# -----------------------------------------------------------------------------
# Backup B4: Shapiro-Wilk detail + full-size Q-Q plot
# -----------------------------------------------------------------------------

resids <- resid(fit_reduced, scaled = TRUE)
shapiro_test <- shapiro.test(resids)
print(shapiro_test)

print(lattice::qqmath(fit_reduced, main = "Normal Q-Q plot"))


# -----------------------------------------------------------------------------
# Backup B5: Variance diagnostics detail (Levene + variance ratios)
# -----------------------------------------------------------------------------

lev_results <- data.frame(
  Factor = c("Region", "Nutrients", "Status", "AMD"),
  p_value = c(
    leveneTest(resid(fit_reduced) ~ Arabidopsis$reg)[1,3],
    leveneTest(resid(fit_reduced) ~ Arabidopsis$nutrient)[1,3],
    leveneTest(resid(fit_reduced) ~ Arabidopsis$status)[1,3],
    leveneTest(resid(fit_reduced) ~ Arabidopsis$amd)[1,3]
  )
)
lev_results$p_value <- signif(lev_results$p_value, 3)
print(kable(lev_results, align = "lc"))

# Max/min within-group residual variance ratios for Region and Status.
var_reg    <- tapply(resid(fit_reduced), Arabidopsis$reg, var)
var_status <- tapply(resid(fit_reduced), Arabidopsis$status, var)
ratio_reg    <- max(var_reg) / min(var_reg)
ratio_status <- max(var_status) / min(var_status)

print(
  data.frame(
    Group = c("Region", "Status"),
    Variance_Ratio = round(c(ratio_reg, ratio_status), 2)
  ) %>% kable(align = "lc")
)


# -----------------------------------------------------------------------------
# Backup B6: Zoomed response distribution (counts < 10)
# -----------------------------------------------------------------------------

print(
  ggplot(Arabidopsis %>% filter(total.fruits < 10), aes(x = total.fruits)) +
    geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
    labs(title = "Distribution of Total Fruits (zoomed to < 10)",
         x = "Total Fruits", y = "Count") +
    theme_minimal()
)
