# Description -------------------------------------------------------------
# Proportion of cell types

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

# Load --------------------------------------------------------------------
meta <- vroom::vroom(here::here("islet_cartography_scrna/data/paper_figures/files/obs.csv"))


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

# Check model assumptions -------------------------------------------------
donor_beta <- donor_props |> 
  dplyr::filter(manual_annotation == "beta")
# Does the sex effect change depending on disease status?
## model ----
model1 <- lm(proportion ~ disease_harmonized + sex_predicted + disease_harmonized:sex_predicted, data = donor_beta)
## Normality ----
# residuals vs fitted
plot(model1, which = 1)
# qqplot
plot(model1, which = 2)
# residual vs fitted at a scale
plot(model1, which = 3)
# spread of cooks distance
plot(model1, which = 4)
# None of my points come close to having high leverage
plot(model1, which = 5)
# Tests that data comes from a normal distribtion
# H0: Data is normally distributed
# HA: The data is not normally distribtued
# W should be close to 1
# W is the fraction between how well the weighted sorted data aligns with what you'd expect from a normal distribution
# and the total variability in the data
# W = (∑ aᵢ x₍ᵢ₎)² / ∑(xᵢ - x̄)²
# x₍ᵢ₎ = ordered data values (sorted from smallest to largest)
# aᵢ = special weights (coefficients) calculated from the normal distribution
# x̄ = mean of your data
# ∑(xᵢ - x̄)² = sum of squared deviations from the mean
shapiro.test(donor_beta$proportion)
# Seems like data is not normally distributed, which is it not:
plot(density(donor_beta$proportion))
hist(donor_beta$proportion, breaks = 30)
## Heteroskedasticity ----
# In a standard linear model, the variance of the residuals are assumed to be constant (i.e. independent) over the values of the response (fitted values)
# These results could indicate that there is some non constant variance of the residuals (heteroscedasticity).
# Breusch-Pagan
ncvTest(model1)
## Multicolinearity ----
# is there correlation between predictor values?
# no - correlation above 4
vif(model1, type = 'predictor')
## Outliers ----
# Test if there are outliers
car::outlierTest(model1)
## Homogeneity of variance ----
# It seems like we have heterogeneity of variance
car::leveneTest(residuals(model1) ~ factor(fitted(model1)))
## Look at infuential / high leverage points ----
# X-axis (Hat-Values): Leverage — how extreme a point's predictor values are
# Y-axis (Standardized Residuals): How far the residual is from zero
# Bubble size: Cook's D — influence on regression coefficients
# Dashed lines: Thresholds for concerning values
# 244, 245, 246: Negative residuals, low leverage → don't worry
# 256, 267: Higher leverage but still low Cook's D → minimal influence
# 506: High leverage AND moderately large residual, but still acceptable (Cook's D ≈ 0.064)
influencePlot(model1)
summary(model1)

# Plot results ------------------------------------------------------------
summary(model1) |> broom::tidy() |> 
  dplyr::mutate(estimate)

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
  dplyr::select(r.squared, adj.r.squared, sigma, statistic, p.value, df, df.residual, nobs)

ggplot(donor_beta, aes(x = disease_harmonized, y = proportion, fill = disease_harmonized)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
  labs(
    title = "Beta cell proportion by disease status",
    x = "Disease group",
    y = "Proportion (%)",
    fill = "Disease"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(donor_beta, aes(x = disease_harmonized, y = proportion, fill = sex_predicted)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
  labs(
    title = "Beta cell proportion by disease and sex",
    x = "Disease group",
    y = "Proportion (%)",
    fill = "Sex"
  ) +
  scale_fill_manual(values = c("f" = "#7F77DD", "m" = "#3266ad"),
                    labels = c("f" = "Female", "m" = "Male")) +
  theme_minimal() +
  theme(legend.position = "top")

donor_beta %>%
  group_by(disease_harmonized, sex_predicted) %>%
  summarise(
    mean = mean(proportion, na.rm = TRUE),
    sd = sd(proportion, na.rm = TRUE),
    se = sd(proportion, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  mutate(
    sex_predicted = ifelse(sex_predicted == "f", "Female", "Male"),
    disease_harmonized = factor(disease_harmonized, 
                                levels = c("nd", "pre", "t2d"),
                                labels = c("Control", "Pre-diabetes", "Type 2 diabetes"))
  ) %>%
  ggplot(aes(x = disease_harmonized, y = mean, fill = sex_predicted)) +
  geom_col(position = "dodge", width = 0.7, alpha = 0.8) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    position = position_dodge(0.7),
    width = 0.2,
    linewidth = 0.5
  ) +
  labs(
    title = "Mean beta cell proportion ± SE",
    x = "Disease group",
    y = "Proportion (%)",
    fill = "Sex"
  ) +
  scale_fill_manual(values = c("Female" = "#7F77DD", "Male" = "#3266ad")) +
  theme_minimal() +
  theme(legend.position = "top", panel.grid.major.x = element_blank())

donor_beta %>%
  group_by(disease_harmonized, sex_predicted) %>%
  summarise(
    n = n(),
    mean = round(mean(proportion, na.rm = TRUE), 2),
    sd = round(sd(proportion, na.rm = TRUE), 2),
    min = round(min(proportion, na.rm = TRUE), 2),
    max = round(max(proportion, na.rm = TRUE), 2),
    .groups = 'drop'
  ) %>%
  arrange(disease_harmonized, sex_predicted) %>%
  knitr::kable()



# Beta --------------------------------------------------------------------
## Data ----
donor_beta <- donor_props |> 
  dplyr::filter(manual_annotation == "beta")
## model ----
model1 <- lm(proportion ~ disease_harmonized + sex_predicted + disease_harmonized:sex_predicted, data = donor_beta)
test <- shapiro.test(donor_beta$proportion) |> 
  broom::tidy()
ncv <- ncvTest(model1)
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
