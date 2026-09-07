# Description -------------------------------------------------------------
# Cell type proportion analysis using binomial mixed models
# Cell type proportions were modelled using a binomial generalised linear mixed model (GLMM) 
# implemented in glmmTMB. For each cell type, the number of cells of that type per sample was 
# modelled against the total number of cells in that sample using cbind(n_cells, total - n_cells), 
# which specifies the number of successes (cells of the target type) and failures (all other cells) respectively. 
# This formulation correctly accounts for differences in sequencing depth between samples 
# (A proportion based on 10,000 cells carries more statistical weight than one based on 100 cells)
# Disease state (disease_harmonized), age (age_years), and sex (sex_predicted) 
# were included as fixed effects, allowing estimation of disease-associated compositional 
# shifts while adjusting for demographic covariates. Dataset (ic_id_dataset) and donor (ic_id_donor_overall) 
# were included as random intercepts to account for the multi-dataset structure of the data and the fact that some 
# donors contributed multiple samples. 
# Estimated marginal means and 95% confidence intervals were extracted using the emmeans package on the response (probability) scale. 
# P-values for each non-reference disease group were extracted from the Wald z-test in the model summary. 
# Multiple testing correction was applied across all cell types and contrasts using the Benjamini-Hochberg false discovery rate procedure.
# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)
library(glmmTMB)

cols25 <- c(
  "#1F51FF", "#d95f02", "#CBC3E3", "#e7298a", "#66a61e",
  "#e6ab02", "#800020", "#988558",
  "#00FFFF", "#b2df8a", "#fb9a99", "#FFEA00",
  "#FF00FF", "#5D3FD3", "#0FFF50", "#b15928",
  "#8dd3c7", "#bebada", "#80b1d3", "#fccde5",
  "#b3de69", "#FF5F1F", "#FFA82E", "#DFFF00", "#FAD5A5"
)

# Load --------------------------------------------------------------------
meta <- vroom::vroom(here::here("islet_cartography_scrna/data/paper_figures/files/obs.csv"))

# Number of samples per donor ---------------------------------------------
samples_per_donor <- meta%>%
  distinct(ic_id_donor_overall, ic_id_dataset, ic_id_platform_adjusted_sample)  |> 
  group_by(ic_id_donor_overall, ic_id_dataset)  |> 
  summarise(n_samples = n())

pdf(here::here("islet_cartography_scrna/data/meta_analysis/plots/samples_per_donor.pdf"),
    height = 3, 
    width = 5)
samples_per_donor |> 
  ggplot2::ggplot(ggplot2::aes(x = ic_id_donor_overall, y = n_samples, fill = ic_id_dataset)) +
  ggplot2::geom_bar(stat="identity") +
  scale_fill_manual(values = cols25) +
  my_theme() +
  ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = element_blank())
dev.off()

# Preprocess --------------------------------------------------------------
meta_counts <- meta |> 
  dplyr::select(ic_id_donor_overall, 
                ic_id_dataset,
                ic_id_platform_adjusted_sample,
                disease_harmonized,
                ethnicity_broad_harmonized,
                hba_1_c_percent,
                bmi,
                sex_predicted,
                age_years,
                diabetes_medication_harmonized) |>
  dplyr::distinct() |> 
  tidyr::replace_na(list("ethnicity_broad_harmonized" = "unknown")) |> 
  dplyr::mutate(
    hba_1_c_percent_imp = tidyr::replace_na(hba_1_c_percent, base::round(mean(hba_1_c_percent, na.rm = TRUE), 2)),
    bmi_imp             = tidyr::replace_na(bmi, base::round(mean(bmi, na.rm = TRUE), 2)),
    dplyr::across(tidyselect::where(is.character), base::as.factor),
    dplyr::across(tidyselect::where(is.numeric),  ~ as.numeric(base::scale(.x)))
  ) |> 
  dplyr::relocate(hba_1_c_percent_imp, .after = hba_1_c_percent) |> 
  dplyr::relocate(bmi_imp, .after = bmi)

# Get number of cells per sample --------------------------------
counts <- meta |>
  dplyr::group_by(ic_id_platform_adjusted_sample, ic_id_donor_overall, cell_type_lvl_3) |>
  dplyr::summarise(n_cells = dplyr::n(), .groups = 'drop') |>
  dplyr::group_by(ic_id_platform_adjusted_sample) |>          # total cells per sample (denominator)
  dplyr::mutate(total = sum(n_cells)) |>
  dplyr::ungroup() |>
  dplyr::left_join(y = meta_counts, by = c("ic_id_platform_adjusted_sample", "ic_id_donor_overall"))


# Which covarivates have an effect ----------------------------------------
base_formula <- cbind(n_cells, total - n_cells) ~ disease_harmonized +
  (1 | ic_id_dataset) + (1 | ic_id_donor_overall)

covariate_terms <- list(
  sex       = ~ disease_harmonized * sex_predicted,
  age       = ~ disease_harmonized * age_years,
  hba1c     = ~ disease_harmonized * hba_1_c_percent_imp,
  bmi       = ~ disease_harmonized * bmi_imp
)

tidy_glmmtmb <- function(m, covariate) {
  m |> summary() |>
    purrr::pluck("coefficients", "cond") |>
    as.data.frame() |>
    tibble::rownames_to_column("variable") |>
    dplyr::rename(log_odds = Estimate, pval = `Pr(>|z|)`, se = `Std. Error`, z = `z value`) |>
    dplyr::mutate(covariate = covariate)
}

glmm_covariates <- function(df) {
  m_base <- glmmTMB(base_formula, family = binomial, data = df)
  
  covariate_terms |>
    purrr::imap(\(term, name) {
      m <- update(m_base, as.formula(paste(". ~ . +", deparse1(term[[2]]))))
      tidy_glmmtmb(m, name) |>
        dplyr::mutate(
          lrt_p = anova(m_base, m) |> as.data.frame() |> dplyr::pull(`Pr(>Chisq)`) |> dplyr::last()
        )
    }) |>
    purrr::list_rbind(names_to = "terms") |>
    dplyr::bind_rows(tidy_glmmtmb(m_base, "base"), x = _)
}

results <- counts |>
  tidyr::nest(.by = "cell_type_lvl_3") |>
  dplyr::mutate(data= purrr::map(data, ~glmm_covariates(base::droplevels(.x)))) |> 
  tidyr::unnest("data") |> 
  dplyr::mutate(fdr = p.adjust(pval, method = "BH")) |> 
  dplyr::mutate(variable = stringr::str_remove(variable, "harmonized"),
                variable = stringr::str_remove(variable, "predicted")) 

main_effects <- results |> dplyr::filter(covariate == "base", !stringr::str_detect(variable, ":"))
interactions <- results |> dplyr::filter(stringr::str_detect(variable, ":"))
vik <- khroma::color("vik")


var_order <- c(
  # Main effects
  "(Intercept)",
  "disease_pre",
  "disease_t2d",
  # Sex
  "sex_m",
  "disease_pre:sex_m",
  "disease_t2d:sex_m",
  # Age
  "age_years",
  "disease_pre:age_years",
  "disease_t2d:age_years",
  # HbA1c
  "hba_1_c_percent_imp",
  "disease_pre:hba_1_c_percent_imp",
  "disease_t2d:hba_1_c_percent_imp",
  # BMI
  "bmi_imp",
  "disease_pre:bmi_imp",
  "disease_t2d:bmi_imp"
)

res_plot <- results |> 
  dplyr::mutate(variable = factor(variable, levels = var_order)) |> 
  dplyr::filter(!variable == "(Intercept)")

res_plot |>
  ggplot(aes(y = variable, x = cell_type_lvl_3, fill = log_odds)) +
  geom_tile() +
  ggplot2::geom_text(data = res_plot |> 
                       dplyr::filter(fdr <= 0.05), 
                       aes(label = paste0(round(fdr,3))), stat = "identity", size = 2, colour = "black") +
  ggplot2::scale_fill_gradientn(
    colors = vik(256),
    limits = c(-2, 2),
    values = scales::rescale(c(-2, 0, 2)),
    oob = scales::squish, name = "Log Odds ratio"
  ) + 
  theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
                 panel.spacing = unit(0.1, "lines"),
                 axis.title = ggplot2::element_blank())
















formulas <- list("base" = cbind(n_cells, total - n_cells) ~ disease_harmonized + (1 | ic_id_dataset) + (1 | ic_id_donor_overall),
                 "sex" = cbind(n_cells, total - n_cells) ~ disease_harmonized * sex_predicted + (1 | ic_id_dataset) + (1 | ic_id_donor_overall),
                 "age" = cbind(n_cells, total - n_cells) ~ disease_harmonized * age_years + (1 | ic_id_dataset) + (1 | ic_id_donor_overall),
                 "ethnicity" = cbind(n_cells, total - n_cells) ~ disease_harmonized * ethnicity_broad_harmonized + (1 | ic_id_dataset) + (1 | ic_id_donor_overall),
                 "hba1c" = cbind(n_cells, total - n_cells) ~ disease_harmonized * hba_1_c_percent_imp + (1 | ic_id_dataset) + (1 | ic_id_donor_overall),
                 "bmi" = cbind(n_cells, total - n_cells) ~ disease_harmonized * bmi_imp + (1 | ic_id_dataset) + (1 | ic_id_donor_overall))

glmm_covariates <- function(df, formula){
  res <- glmmTMB(.x,
          family = binomial, 
          data = df) |> 
    summary() |> 
    pluck("coefficients") |> 
    pluck("cond") |> 
    as.data.frame() |> 
    tibble::rownames_to_column("interactions")
  
  main_effects   <- results |> dplyr::filter(covariate == "base", !stringr::str_detect(variables, ":"))
  interactions   <- results |> dplyr::filter(stringr::str_detect(variables, ":"))}

test <- map(formulas, ~ glmmTMB(.x,
                       family = binomial, 
                       data = df) |> 
      summary() |> 
        pluck("coefficients") |> 
        pluck("cond") |> 
        as.data.frame() |> 
        tibble::rownames_to_column("variables")) |> 
  purrr::list_rbind(names_to = "covariate")
test

fit <- glmmTMB(formulas[["bmi"]],
               family = binomial, 
               data = df)

summary(fit)

fit_additive <- glmmTMB(
  cbind(n_cells, total - n_cells) ~ disease_harmonized + age_years + sex_predicted +
    (1 | ic_id_dataset) + (1 | ic_id_donor_overall),
  family = binomial, data = df
)

fit_interaction <- glmmTMB(
  cbind(n_cells, total - n_cells) ~ disease_harmonized * sex_predicted + age_years +
    (1 | ic_id_dataset) + (1 | ic_id_donor_overall),
  family = binomial, data = df
)

# Test whether interaction improves fit
anova(fit_additive, fit_interaction)

# Inspect interaction effects
summary(fit_interaction)

# Marginal means by disease and sex
emmeans::emmeans(fit_interaction, ~ disease_harmonized | sex_predicted, type = "response")
# Function ----------------------------------------------------------------
prop_fit <- function(df, formula = cbind(n_cells, total - n_cells) ~ disease_harmonized + age_years + sex_predicted + (1 | ic_id_dataset) + (1 | ic_id_donor_overall)){
  fit <- glmmTMB(
    formula,
    family = binomial,
    data   = df
  )
  
  summ  <- summary(fit)
  emm   <- summary(emmeans::emmeans(fit, ~ disease_harmonized, type = "response"))
  cond  <- summ$coefficients$cond
  pvals <- cond[grepl("^disease_harmonized", rownames(cond)), "Pr(>|z|)"]
  
  data.frame(
    state  = emm$disease_harmonized,
    pvalue = c(NA, pvals),
    prob   = emm$prob,
    lower  = emm$asymp.LCL,
    upper  = emm$asymp.UCL
  )
}

# Run ---------------------------------------------------------------------
test <- counts |> 
  tidyr::nest(.by = "cell_type_lvl_2") |>
  dplyr::mutate(data = purrr::map(data, base::droplevels)) |>
  dplyr::filter(purrr::map_vec(data, ~ nlevels(dplyr::pull(.x, disease_harmonized)) > 1)) |>
  dplyr::mutate(prop = purrr::map(data, prop_fit)) |> 
  dplyr::select(-data) |> 
  tidyr::unnest(prop) |> 
  dplyr::mutate(qval = stats::p.adjust(pvalue, method = "BH"))

test |> 
  ggplot2::ggplot(ggplot2::aes(x = state, y = prob)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::facet_wrap(~cell_type_lvl_2, scales = "free_y") +
  ggplot2::geom_text(aes(label = paste0("q=", round(qval,3))), stat = "identity", size = 2, vjust = 1.5, colour = "white") +
  my_theme()










prop_fit <- function(df){
  fit <- glmmTMB(cbind(n_cells, total - n_cells) ~ disease_harmonized + age_years + sex_predicted,
                 family = binomial,
                 data = df)
  summ <- summary(fit)
  emm <- summary(emmeans::emmeans(fit, ~ disease_harmonized, type = "response"))
  
  cond    <- summ$coefficients$cond
  pvals   <- cond[grepl("^disease_harmonized", rownames(cond)), "Pr(>|z|)"]
  
  return(data.frame(
    state  = emm$disease_harmonized,
    pvalue = c(NA, pvals),
    prob   = emm$prob,
    lower  = emm$asymp.LCL,
    upper  = emm$asymp.UCL
  ))
}

# Number of cells per donor
counts <- meta |>
  dplyr::group_by(ic_id_donor_overall, cell_type_lvl_2) |>
  dplyr::summarise(n_cells = dplyr::n(), .groups = 'drop') |>
  dplyr::left_join(y = meta_counts, relationship = "many-to-many") |>
  dplyr::group_by(ic_id_donor_overall) |> 
  dplyr::mutate(total = sum(n_cells)) |>
  dplyr::ungroup()

test <- counts |> 
  tidyr::nest(.by = c("cell_type_lvl_2")) |>
  dplyr::mutate(data = purrr::map(data, ~ base::droplevels(.x))) |> # drop unused levels
  dplyr::filter(purrr::map_vec(data, ~ dplyr::pull(.x, "disease_harmonized") |> nlevels() > 1)) |>  # filter datasets with only one disease
  dplyr::mutate(prop = purrr::map(data, prop_fit)) |> 
  dplyr::select(-data) |> 
  tidyr::unnest('prop')

# Perform analysis --------------------------------------------------------
df <- counts |> 
  dplyr::filter(ic_id_dataset == "ic_13" & cell_type_lvl_2 == "alpha")
fit <- glmmTMB(cbind(n_cells, total - n_cells) ~ disease_harmonized + age_years + sex_predicted,
               family = binomial,
               data = df)
summ <- summary(fit)
emm <- summary(emmeans::emmeans(fit, ~ disease_harmonized, type = "response"))


library(glmmTMB)
test <- counts |> 
  tidyr::nest(.by = c("ic_id_dataset", "cell_type_lvl_2")) |>
  dplyr::mutate(data = purrr::map(data, ~ base::droplevels(.x))) |> # drop unused levels
  dplyr::filter(purrr::map_vec(data, ~ dplyr::pull(.x, "disease_harmonized") |> nlevels() > 1)) |>  # filter datasets with only one disease
  dplyr::mutate(prop = purrr::map(data, prop_fit)) |> 
  dplyr::select(-data) |> 
  tidyr::unnest('prop')

# Proportions become very small, look into this
test |> 
  ggplot2::ggplot(ggplot2::aes(x = state, y = prob)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::facet_wrap(~cell_type_lvl_2, scales = "free_y")


test_celltypes <- colnames(counts)[1:15]

for (test_celltype in test_celltypes) {
  combined$celltype <- combined[,test_celltype]
  fit <- glmmTMB(cbind(celltype, total - celltype) ~ disease_harmonized_hba1c + age + gender + (1 | study), family = binomial, data = combined)
  summ <- summary(fit)
  emm <- emmeans(fit, ~ disease_harmonized_hba1c, type = "response")
  emm <- summary(emm)
  probs <- emm$prob
  LCL <- emm$asymp.LCL
  UCL <- emm$asymp.UCL
  pvals <- c(1, summ$coefficients$cond[, "Pr(>|z|)"][2:3])
  tmp <- data.frame(celltype = rep(test_celltype, 3), state = c("ND", "PRE", "T2D"), pvalue = as.vector(pvals), prob = as.vector(probs), lcl = as.vector(LCL), ucl = as.vector(UCL))
  if (test_celltype == test_celltypes[1]) {
    res <- tmp
  } else {
    res <- rbind(res, tmp)
  }
}
res$star <- ifelse(res$pvalue < 0.05, "*", "")
res$y_star <- res$ucl + (res$ucl - res$prob) * 0.1
ggplot(data = res, aes(x = state, y = prob, fill = celltype)) + geom_bar(stat = "identity") + geom_text(aes(y = y_star, label = star), vjust = 0) + geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.2) + facet_wrap(~ celltype, scales = "free_y") + theme_minimal()

# samples per donor -------------------------------------------------------
samples_per_donor <- meta%>%
  distinct(ic_id_donor_overall, ic_id_platform_adjusted_sample)  |> 
  group_by(ic_id_donor_overall)  |> 
  summarise(n_samples = n())

samples_per_donor |> 
  ggplot2::ggplot(ggplot2::aes(x = forcats::fct_reorder(ic_id_donor_overall, n_samples), y = n_samples)) +
  ggplot2::geom_bar(stat="identity") +
  my_theme() +
  ggplot2::theme(axis.text.x = ggplot2::element_blank())

# Calculate cell type proportions -----------------------------------------
# per dataset
# per donor
donor_props <- meta %>%
  dplyr::group_by(ic_id_donor_overall, manual_annotation) |> 
  dplyr::summarise(n_cells = n(), .groups = 'drop')  |> 
  dplyr::group_by(ic_id_donor_overall) |> 
  dplyr::mutate(proportion = (n_cells / sum(n_cells))*100,
                log_proportion = log(proportion))  |> 
  dplyr::ungroup() |> 
  dplyr::left_join(y = meta |> dplyr::select(ic_id_donor_overall, disease_harmonized, 
                                             ethnicity_broad_harmonized, 
                                             hba_1_c_percent, bmi, 
                                             sex_predicted, age_years) |> 
                     dplyr::distinct())

# # Check model assumptions -------------------------------------------------
# donor_beta <- donor_props |> 
#   dplyr::filter(manual_annotation == "beta")
# # Does the sex effect change depending on disease status?
# ## model ----
# model1 <- lm(proportion ~ disease_harmonized + sex_predicted + disease_harmonized:sex_predicted, data = donor_beta)
# ## Normality ----
# # residuals vs fitted
# plot(model1, which = 1)
# # qqplot
# plot(model1, which = 2)
# # residual vs fitted at a scale
# plot(model1, which = 3)
# # spread of cooks distance
# plot(model1, which = 4)
# # None of my points come close to having high leverage
# plot(model1, which = 5)
# # Tests that data comes from a normal distribtion
# # H0: Data is normally distributed
# # HA: The data is not normally distribtued
# # W should be close to 1
# # W is the fraction between how well the weighted sorted data aligns with what you'd expect from a normal distribution
# # and the total variability in the data
# # W = (∑ aᵢ x₍ᵢ₎)² / ∑(xᵢ - x̄)²
# # x₍ᵢ₎ = ordered data values (sorted from smallest to largest)
# # aᵢ = special weights (coefficients) calculated from the normal distribution
# # x̄ = mean of your data
# # ∑(xᵢ - x̄)² = sum of squared deviations from the mean
# shapiro.test(donor_beta$proportion)
# # Seems like data is not normally distributed, which is it not:
# plot(density(donor_beta$proportion))
# hist(donor_beta$proportion, breaks = 30, main = "Proportions beta cells")
# ## Heteroskedasticity ----
# # In a standard linear model, the variance of the residuals are assumed to be constant (i.e. independent) over the values of the response (fitted values)
# # These results could indicate that there is some non constant variance of the residuals (heteroscedasticity).
# # Breusch-Pagan
# car::ncvTest(model1)
# ## Multicolinearity ----
# # is there correlation between predictor values?
# # no - correlation above 4
# car::vif(model1, type = 'predictor')
# ## Outliers ----
# # Test if there are outliers
# car::outlierTest(model1)
# ## Homogeneity of variance ----
# # It seems like we have heterogeneity of variance
# car::leveneTest(residuals(model1) ~ factor(fitted(model1)))
# ## Look at infuential / high leverage points ----
# # X-axis (Hat-Values): Leverage — how extreme a point's predictor values are
# # Y-axis (Standardized Residuals): How far the residual is from zero
# # Bubble size: Cook's D — influence on regression coefficients
# # Dashed lines: Thresholds for concerning values
# # 244, 245, 246: Negative residuals, low leverage → don't worry
# # 256, 267: Higher leverage but still low Cook's D → minimal influence
# # 506: High leverage AND moderately large residual, but still acceptable (Cook's D ≈ 0.064)
# car::influencePlot(model1)
# summary(model1)
# 
# # Plot results ------------------------------------------------------------
# summary(model1) |> broom::tidy() |> 
#   dplyr::mutate(estimate)
# 
# broom::tidy(model1) %>%
#   mutate(
#     estimate = round(estimate, 3),
#     std.error = round(std.error, 3),
#     statistic = round(statistic, 3),
#     p.value = round(p.value, 4),
#     significance = case_when(
#       p.value < 0.05 ~ "*",
#       TRUE ~ "ns"
#     )
#   ) %>%
#   knitr::kable()
# 
# broom::glance(model1) %>%
#   dplyr::select(r.squared, adj.r.squared, sigma, statistic, p.value, df, df.residual, nobs)
# 
# plot(density(donor_beta$proportion), main = "Proportions of beta cells")
# 
# ggplot(donor_beta, aes(x = disease_harmonized, y = proportion, fill = disease_harmonized)) +
#   geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
#   labs(
#     title = "Beta cell proportion by disease status",
#     x = "Disease group",
#     y = "Proportion (%)",
#     fill = "Disease"
#   ) +
#   scale_fill_manual(values = disease_color) +
#   my_theme() +
#   theme(legend.position = "none")
# 
# ggplot(donor_beta, aes(x = disease_harmonized, y = proportion, fill = sex_predicted)) +
#   geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
#   labs(
#     title = "Beta cell proportion by disease and sex",
#     x = "Disease group",
#     y = "Proportion (%)",
#     fill = "Sex"
#   ) +
#   scale_fill_manual(values = c("f" = "#7F77DD", "m" = "#3266ad"),
#                     labels = c("f" = "Female", "m" = "Male")) +
#   theme_minimal() +
#   theme(legend.position = "top")
# 
# donor_beta %>%
#   group_by(disease_harmonized, sex_predicted) %>%
#   summarise(
#     mean = mean(proportion, na.rm = TRUE),
#     sd = sd(proportion, na.rm = TRUE),
#     se = sd(proportion, na.rm = TRUE) / sqrt(n()),
#     .groups = 'drop'
#   ) %>%
#   mutate(
#     sex_predicted = ifelse(sex_predicted == "f", "Female", "Male"),
#     disease_harmonized = factor(disease_harmonized, 
#                                 levels = c("nd", "pre", "t2d"),
#                                 labels = c("Control", "Pre-diabetes", "Type 2 diabetes"))
#   ) %>%
#   ggplot(aes(x = disease_harmonized, y = mean, fill = sex_predicted)) +
#   geom_col(position = "dodge", width = 0.7, alpha = 0.8) +
#   geom_errorbar(
#     aes(ymin = mean - se, ymax = mean + se),
#     position = position_dodge(0.7),
#     width = 0.2,
#     linewidth = 0.5
#   ) +
#   labs(
#     title = "Mean beta cell proportion ± SE",
#     x = "Disease group",
#     y = "Proportion (%)",
#     fill = "Sex"
#   ) +
#   scale_fill_manual(values = c("Female" = "#7F77DD", "Male" = "#3266ad")) +
#   theme_minimal() +
#   theme(legend.position = "top", panel.grid.major.x = element_blank())
# 
# donor_beta %>%
#   group_by(disease_harmonized, sex_predicted) %>%
#   summarise(
#     n = n(),
#     mean = round(mean(proportion, na.rm = TRUE), 2),
#     sd = round(sd(proportion, na.rm = TRUE), 2),
#     min = round(min(proportion, na.rm = TRUE), 2),
#     max = round(max(proportion, na.rm = TRUE), 2),
#     .groups = 'drop'
#   ) %>%
#   arrange(disease_harmonized, sex_predicted) %>%
#   knitr::kable()
# 
# 

# Beta --------------------------------------------------------------------
## Data ----
donor_beta <- donor_props |> 
  dplyr::filter(manual_annotation == "beta")
## model ----
model1 <- lm(proportion ~ disease_harmonized + sex_predicted + disease_harmonized:sex_predicted, data = donor_beta)
test <- shapiro.test(donor_beta$proportion) |> 
  broom::tidy()
ncv <- car::ncvTest(model1)
lv <- car::leveneTest(residuals(model1) ~ factor(fitted(model1))) |> 
  broom::tidy() |> 
  dplyr::select(statistic, p.value) |> 
  dplyr::mutate(method= "Levene's test")
car::outlierTest(model1)
test <- test |> 
  dplyr::bind_rows(base::data.frame("statistic" = ncv$ChiSquare, "p.value" = ncv$p, "method" = "Non-Constant Error Variance")) |> 
  dplyr::bind_rows(lv)

test %>% 
  knitr::kable()

broom::tidy(model1) %>%
  mutate(
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    statistic = round(statistic, 3),
    p.value = round(p.value, 4),
    significance = case_when(
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  knitr::kable()

broom::glance(model1) %>%
  knitr::kable()

# Our intercept is 22.3
# slope is 11.9
# std.error: how far does our coefficient estimate vary from the actual data
# statistic: t-statistic, how many std does our coefficient estimate deviate from 0
# residual standard error  (sigma) the average amount our response (proportion) will differ from the true regression line.
# r.squared: how well does our model explain the data. around 6.9% of variance found in the response variable (proportion)
# can be explain the by the predictor variable (disease + sex)

ggplot(donor_beta, aes(x = disease_harmonized, y = proportion, fill = disease_harmonized)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
  labs(
    title = "Beta cell proportion by disease status",
    x = "Disease group",
    y = "Proportion (%)",
    fill = "Disease"
  ) +
  scale_fill_manual(values = disease_color) +
  my_theme() +
  theme(legend.position = "none")

ggplot(donor_beta, aes(x = disease_harmonized, y = proportion, fill = sex_predicted)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
  labs(
    title = "Beta cell proportion by disease and sex",
    x = "Disease group",
    y = "Proportion (%)",
    fill = "Sex"
  ) +
  scale_fill_manual(values = c("f" = "#C2563A", "m" = "#3F7F93"),
                    labels = c("f" = "F", "m" = "M")) +
  my_theme() +
  theme(legend.position = "none")
par(mfrow=c(2,2))
plot(model1)


# vs hb11c ----------------------------------------------------------------
## Data ----
donor_beta <- donor_props |> 
  dplyr::filter(manual_annotation == "beta")
## model ----
model1 <- lm(proportion ~ disease_harmonized + hba_1_c_percent + disease_harmonized:hba_1_c_percent, data = donor_beta)
test <- shapiro.test(donor_beta$proportion) |> 
  broom::tidy()
ncv <- car::ncvTest(model1)
lv <- car::leveneTest(residuals(model1) ~ factor(fitted(model1))) |> 
  broom::tidy() |> 
  dplyr::select(statistic, p.value) |> 
  dplyr::mutate(method= "Levene's test")
car::outlierTest(model1)
test <- test |> 
  dplyr::bind_rows(base::data.frame("statistic" = ncv$ChiSquare, "p.value" = ncv$p, "method" = "Non-Constant Error Variance")) |> 
  dplyr::bind_rows(lv)

broom::tidy(model1) %>%
  mutate(
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    statistic = round(statistic, 3),
    p.value = round(p.value, 4),
    significance = case_when(
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  knitr::kable()

broom::glance(model1) %>%
  knitr::kable()

# Our intercept is 22.3
# slope is 11.9
# std.error: how far does our coefficient estimate vary from the actual data
# statistic: t-statistic, how many std does our coefficient estimate deviate from 0
# residual standard error  (sigma) the average amount our response (proportion) will differ from the true regression line.
# r.squared: how well does our model explain the data. around 6.9% of variance found in the response variable (proportion)
# can be explain the by the predictor variable (disease + sex)

ggplot(donor_beta, aes(x = disease_harmonized, y = proportion, fill = disease_harmonized)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
  labs(
    title = "Beta cell proportion by disease status",
    x = "Disease group",
    y = "Proportion (%)",
    fill = "Disease"
  ) +
  scale_fill_manual(values = disease_color) +
  my_theme() +
  theme(legend.position = "none")

ggplot(donor_beta, aes(y = hba_1_c_percent, x = proportion, color = disease_harmonized)) +
  geom_point() +
  scale_color_manual(values = disease_color) +
  my_theme() +
  theme(legend.position = "none")

par(mfrow=c(2,2))
plot(model1)


# Ethnicity ---------------------------------------------------------------
## Data ----
donor_beta <- donor_props |> 
  dplyr::filter(manual_annotation == "beta") |> 
  dplyr::select(proportion, disease_harmonized, ethnicity_broad_harmonized) |> 
  tidyr::drop_na()
## model ----
model1 <- lm(proportion ~ disease_harmonized + ethnicity_broad_harmonized + disease_harmonized:ethnicity_broad_harmonized, data = donor_beta)
test <- shapiro.test(donor_beta$proportion) |> 
  broom::tidy()
ncv <- car::ncvTest(model1)
lv <- car::leveneTest(residuals(model1) ~ factor(fitted(model1))) |> 
  broom::tidy() |> 
  dplyr::select(statistic, p.value) |> 
  dplyr::mutate(method= "Levene's test")
car::outlierTest(model1)
test <- test |> 
  dplyr::bind_rows(base::data.frame("statistic" = ncv$ChiSquare, "p.value" = ncv$p, "method" = "Non-Constant Error Variance")) |> 
  dplyr::bind_rows(lv)

broom::tidy(model1) %>%
  mutate(
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    statistic = round(statistic, 3),
    p.value = round(p.value, 4),
    significance = case_when(
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  knitr::kable()

broom::glance(model1) %>%
  knitr::kable()


