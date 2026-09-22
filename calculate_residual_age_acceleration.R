#####################################################################
# Calculation of residual age acceleration (AgeAccel)               #
#####################################################################

# External reference: age acceleration measured relative to publicly available population (n = 17,726)
# Test set (n = 3,502) of hold-out samples to avoid bias with respect to age clock refitting

# Author: David Martínez-Enguita (2025)

pack_R <- c("dplyr", "tidyr", "purrr", "broom",
            "readr", "effectsize", "ggplot2")
for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

set.seed(777)

## Calculation of epigenetic age acceleration residuals

# Load cohort-specific epigenetic age predictions
pred_age_COVUM <- read.csv("~/Documents/COVUM_Project/Results/pred_age_COVUM_all_year_adjusted_sel_clocks.csv")
pred_age_COVUM_sel <- pred_age_COVUM[, c("True_Age", "Horvath", "Hannum", "Horvath_skin", "PhenoAge", "Zhang", "NCAE_Age")]

# Load epigenetic age predictions for DNAm age clock refitting data
pred_refitted <- read.csv("~/Documents/COVUM_Project/Results/pred_allmethods_refitted_models.csv")
pred_refitted <- pred_refitted[, c("age", "Horvath", "Hannum", "Horvath_skin", "PhenoAge", "Zhang", "NCAE_Age")]

# Match age of test set to cohort before estimation of coefficients:
# Divide into 5-year bins and sample the test set so its age histogram matches the specific cohort
range_cohort <- range(pred_age_COVUM$True_Age)
test_trim <- pred_refitted %>%
  filter(age >= range_cohort[1], age <= range_cohort[2])

bin_width <- 5
bins <- seq(floor(range_cohort[1]/bin_width)*bin_width,
            ceiling(range_cohort[2]/bin_width)*bin_width,
            by = bin_width)

# Count number of bins in cohort
target_n <- pred_age_COVUM %>%
  mutate(age_bin = cut(True_Age, bins, right = FALSE)) %>%
  count(age_bin, name = "n_target")

# Sample without replacement from test set to match n_target
test_match <- test_trim %>%
  mutate(age_bin = cut(age, bins, right = FALSE)) %>%
  inner_join(target_n, by = "age_bin") %>%
  group_modify(~ {
    size <- min(.y$n_target, nrow(.x))
    .x[sample(nrow(.x), size), ]
  }) %>%
  ungroup() %>%
  select(-n_target)

# Fit linear regression model in the test set
# EpiAge = beta0 + beta1 x chronological age
clock_cols <- c("Horvath", "Hannum", "Horvath_skin", "PhenoAge", "Zhang", "NCAE_Age")

coefs_test <- map_dfr(clock_cols, function(cl){
  fit <- lm(reformulate("age", response = cl), data = test_match)
  tidy(fit) %>% 
    select(term, beta = estimate) %>%
    mutate(clock = cl, .before = 1)
}) %>%
  pivot_wider(names_from = term,
              values_from = beta,
              names_glue  = "{ifelse(term == '(Intercept)', 'beta0', 'beta1')}")

# Store fitted coefficients per age clock
write.csv(coefs_test, "~/Documents/COVUM_Project/Results/linear_regression_clock_coeffs_WB_test.csv", 
          row.names = FALSE)

# Calculate age acceleration residuals in cohort
ageaccel_out <- pred_age_COVUM %>%
  pivot_longer(all_of(clock_cols),
               names_to = "clock",
               values_to = "epi_age") %>%
  left_join(coefs_test, by = "clock") %>% 
  mutate(predicted = beta0 + beta1 * True_Age,
         AgeAccel  = epi_age - predicted)

ageaccel_wide <- ageaccel_out %>% 
  pivot_wider(id_cols = c(Sample_Name, Sample_Group, Sex, Timepoint,
                          Study_ID, True_Age),
              names_from = clock,
              values_from = c(epi_age, predicted, AgeAccel), 
              names_glue = "{.value}_{clock}")

# Check for correlation of age acceleration with age
ageaccel_out %>%
  group_by(clock) %>%
  summarise(corr = cor(AgeAccel, True_Age))

# Store age acceleration residuals
write.csv(ageaccel_wide, "~/Documents/COVUM_Project/Results/AgeAccel_residuals_COVUM_all.csv", 
          row.names = FALSE)

## Determine whether patients show larger epigenetic age acceleration than controls

# Keep only the two groups of interest and AgeAccel columns
ageaccel_sub <- ageaccel_wide %>%
  filter(Sample_Group %in% c("PCC", "PCC_control")) %>%
  filter(Timepoint != 24)

clock_cols   <- grep("^AgeAccel_", names(ageaccel_sub), value = TRUE)

ageaccel_long <- ageaccel_sub %>%
  select(Sample_Name, Sample_Group, all_of(clock_cols)) %>%
  pivot_longer(all_of(clock_cols),
               names_to   = "clock",
               names_prefix = "AgeAccel_",
               values_to  = "AgeAccel")

# Run a Welch's two-sample t-test for each clock
results_tbl <- ageaccel_long %>%
  group_by(clock) %>%
  summarise(
    n_PCC = sum(Sample_Group == "PCC"),
    n_control = sum(Sample_Group == "PCC_control"),
    mean_PCC = mean(AgeAccel[Sample_Group == "PCC"]),
    sd_PCC = sd(AgeAccel[Sample_Group == "PCC"]),
    mean_control = mean(AgeAccel[Sample_Group == "PCC_control"]),
    sd_control = sd(AgeAccel[Sample_Group == "PCC_control"]),
    
    # Welch t-test (unequal variances)
    t_out = list(t.test(AgeAccel ~ Sample_Group, var.equal = FALSE)),
    p_value = t_out[[1]]$p.value,
    diff = mean_PCC - mean_control,
    ci_low = t_out[[1]]$conf.int[1],
    ci_high = t_out[[1]]$conf.int[2],
    
    # Cohen's d for effect size (Hedges' correction)
    d = cohens_d(AgeAccel, Sample_Group, adjust = TRUE)$Cohens_d
  ) %>%
  ungroup() %>%
  select(-t_out)

compare_out <- as.data.frame(results_tbl)

# Store statistics for age acceleration comparison between groups
write.csv(compare_out, "~/Documents/COVUM_Project/Results/stats_signif_PCC_vs_PCC_controls_AgeAccel.csv", 
          row.names = FALSE)

## Boxplot of age acceleration residuals per age clock and sample group

# Prepare input
plot_long <- ageaccel_wide %>%
  filter(Sample_Group %in% c("PCC", "PCC_control"),
         Timepoint != 24) %>%
  pivot_longer(
    starts_with("AgeAccel_"),
    names_to    = "clock",
    names_prefix = "AgeAccel_",
    values_to   = "AgeAccel_years")

# Add z-score column
plot_long <- plot_long %>%
  group_by(clock) %>%
  mutate(AgeAccel_z = (AgeAccel_years - mean(AgeAccel_years)) /
           sd(AgeAccel_years)) %>%
  ungroup()

# Boxplot of AgeAccel residuals
pdf(file = "~/Documents/COVUM_Project/Results/AgeAccel_residuals_boxplot_PCC_vs_control.pdf", width = 8, height = 4.5)
ggplot(plot_long,
       aes(x = clock,
           y = AgeAccel_years,
           fill = Sample_Group)) +
  geom_boxplot(outlier.shape = NA, 
               width = 0.6, 
               position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Age acceleration residual (years)",
       x = "Epigenetic clock",
       fill = "") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "top")
dev.off()

# Boxplot of z-scores
pdf(file = "~/Documents/COVUM_Project/Results/AgeAccel_zscore_boxplot_PCC_vs_control.pdf", width = 8, height = 4.5)
ggplot(plot_long,
       aes(x = clock,
           y = AgeAccel_z,
           fill = Sample_Group)) +
  geom_boxplot(outlier.shape = NA, 
               width = 0.6, 
               position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Age acceleration (z-score)",
       x = "Epigenetic clock",
       fill = "") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "top")
dev.off()


# COVUM manuscript figure
# 1) Divide by age clock
# 2) Filter for time points 3, 6 and 12
# 3) Calculate the statistical significance of the shown comparisons between PCC and PCC_control as the P-value of a Mann-Whitney U test
# 4) Plot the age acceleration values (as Z-scores).

library(dplyr)
library(ggplot2)
library(ggsignif)
library(readr)
library(rlang)


make_epiage_plot <- function(dat,
                             y_var = c("AgeAccel_years", "AgeAccel_z"),
                             y_limits = NULL,
                             save_stats_csv = NULL) {
  
  y_var <- match.arg(y_var)
  v <- sym(y_var)
  
  if (is.null(y_limits)) {
    y_limits <- if (y_var == "AgeAccel_years") c(-20, 30) else c(-3.5, 3.5)
  }
  
  clocks_order <- c("Horvath","Hannum","Horvath_skin","PhenoAge",
                    "Zhang","NCAE_Age", "mean_EpiAge")
  
  # Facet labels
  clock_labels <- c(
    Horvath      = "Horvath",
    Hannum       = "Hannum",
    Horvath_skin = "Horvath skin & blood",
    PhenoAge     = "PhenoAge",
    Zhang        = "Zhang",
    NCAE_Age     = "NCAE-Age", 
    mean_EpiAge  = "mean_EpiAge"
  )
  
  # Base data (keep rows that have either selected outcome or, for mean_EpiAge in z mode, years)
  df <- dat %>%
    filter(Timepoint %in% c(3, 6, 12),
           Sample_Group %in% c("PCC", "PCC_control"),
           clock %in% clocks_order) %>%
    mutate(
      Group = recode(Sample_Group, PCC = "PCC+", PCC_control = "PCC-"),
      Group = factor(Group, levels = c("PCC+", "PCC-")),
      Timepoint_f = factor(Timepoint, levels = c(3, 6, 12),
                           labels = c("3m", "6m", "1y")),
      # use z everywhere except force years for mean_EpiAge when plotting z
      y_plot = if (y_var == "AgeAccel_z") if_else(clock == "mean_EpiAge", AgeAccel_years, AgeAccel_z)
      else AgeAccel_years,
      clock_f = recode(clock, !!!clock_labels)
    ) %>%
    filter(is.finite(y_plot)) %>%                # drop rows without the needed outcome
    mutate(clock_f = factor(clock_f, levels = unname(clock_labels)[match(clocks_order, names(clock_labels))]))
  
  # Per–clock × timepoint Mann–Whitney on y_plot (PCC+ vs PCC-)
  stats_signif <- df %>%
    group_by(clock_f, Timepoint_f) %>%
    group_modify(~{
      x <- .x %>% filter(Group == "PCC+") %>% pull(y_plot)
      y <- .x %>% filter(Group == "PCC-") %>% pull(y_plot)
      
      W <- p <- NA_real_
      if (length(x) && length(y) && any(is.finite(x)) && any(is.finite(y))) {
        wt <- wilcox.test(x, y, exact = FALSE)
        W  <- as.numeric(wt$statistic)
        p  <- wt$p.value
      }
      
      tibble(
        n_PCC        = sum(is.finite(x)),
        n_PCC_ctrl   = sum(is.finite(y)),
        median_PCC   = median(x, na.rm = TRUE),
        median_ctrl  = median(y, na.rm = TRUE),
        mean_PCC     = mean(x, na.rm = TRUE),
        mean_ctrl    = mean(y, na.rm = TRUE),
        diff_medians = median(x, na.rm = TRUE) - median(y, na.rm = TRUE),
        W = W, p = p
      )
    }) %>%
    ungroup() %>%
    mutate(
      n1 = n_PCC, n2 = n_PCC_ctrl,
      U1 = ifelse(is.na(W), NA_real_, W - n1 * (n1 + 1) / 2),
      U2 = ifelse(is.na(U1), NA_real_, n1 * n2 - U1),
      U  = pmin(U1, U2, na.rm = FALSE),
      r_rb = 1 - (2 * U) / (n1 * n2),
      prob_superiority = U / (n1 * n2),
      p_adj = p.adjust(p, method = "BH"),
      label = case_when(
        is.na(p) ~ "NA",
        p < 1e-4 ~ "****",
        p < 1e-3 ~ "***",
        p < 1e-2 ~ "**",
        p < 0.05 ~ "*",
        TRUE     ~ "ns"
      )
    )
  
  # bracket heights from y_plot
  y_max <- df %>%
    group_by(clock_f, Timepoint_f) %>%
    summarise(y = max(y_plot, na.rm = TRUE), .groups = "drop")
  
  stats_for_plot <- stats_signif %>%
    left_join(y_max, by = c("clock_f", "Timepoint_f")) %>%
    mutate(
      xnum = as.numeric(Timepoint_f),
      xmin = xnum - 0.22, xmax = xnum + 0.22,
      y_position = pmin(y_limits[2] - 0.5, y + diff(y_limits) * 0.08),
      annotations = as.character(label)
    )
  
  if (!is.null(save_stats_csv)) readr::write_csv(stats_signif, save_stats_csv)
  
  # Axis label
  y_lab <- if (y_var == "AgeAccel_z") "EpiAge acceleration (z-score; mean panel in years)"
  else "EpiAge acceleration (years)"
  
  # Plot
  p <- ggplot(df, aes(Timepoint_f, y_plot)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_boxplot(
      aes(group = interaction(Timepoint_f, Group), fill = Timepoint_f),
      width = 0.62,
      position = position_dodge2(width = 0.7, preserve = "single"),
      outlier.shape = 21,
      outlier.size  = 1.6,
      outlier.fill  = "white",
      outlier.color = "black",
      outlier.stroke = 0.3,
      alpha = 0.9,
      coef = 1.5                  # IQR rule for defining outliers (default)
    ) +
    ggsignif::geom_signif(
      data = stats_for_plot,
      aes(xmin = xmin, xmax = xmax, y_position = y_position, annotations = annotations),
      inherit.aes = FALSE, manual = TRUE, na.rm = TRUE,
      tip_length = 0.01, textsize = 3, size = 0.25
    ) +
    facet_wrap(~ clock_f, nrow = 2, ncol = 4, drop = FALSE) +
    coord_cartesian(ylim = y_limits, expand = FALSE) +
    scale_y_continuous(breaks = pretty(y_limits, n = 7)) +
    scale_fill_manual(
      values = c("3m" = "#8C510A", "6m" = "#5AB4AC", "1y" = "#DFC27D"),
      name = NULL, labels = c("3 months", "6 months", "1 year")
    ) +
    labs(x = NULL, y = y_lab) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "right",
          strip.background = element_rect(fill = "grey95"),
          strip.text = element_text(face = "bold"))
  
  list(plot = p, stats_signif = stats_signif)
}

res <- make_epiage_plot(plot_long_years, y_var = "AgeAccel_years", y_limits = c(-20, 30))
p <- res$plot
stats_signif <- res$stats_signif

ggsave("EpiAge_residual_accel_boxplot_years.pdf", p, width = 12, height = 6.5, dpi = 300)

