# =============================================================================
# EHR COHORT ANALYSIS PIPELINE — DEMONSTRATION
# Script 2: Descriptive Statistics and Preliminary Survival Analysis
#
# Produces:
#   - Table 1: baseline characteristics by exposure group (PPI vs H2RA)
#   - Time-to-event summary by exposure group
#   - Kaplan-Meier survival curves with a log-rank test
#
# Inputs:  data/processed/analysis_data.Rds  (produced by script 01)
# Outputs: results/table1_baseline_characteristics.csv, results/km_plot.pdf
# =============================================================================


# =============================================================================
# 0. PACKAGES AND DATA
# =============================================================================

library(tidyverse)
library(survival)
library(survminer)   # ggsurvplot() for the KM plot

dir.create("results", showWarnings = FALSE)

analysis_data <- readRDS("data/processed/analysis_data.Rds")

analysis_data <- analysis_data %>%
  mutate(
    gender = factor(gender),
    ethnicity = factor(ethnicity,
                       levels = c("White", "South-Asian", "Black", "Mixed", "Other")),
    ses_person = factor(ses_person, levels = c("1", "2", "3", "4", "5"))
  )


# =============================================================================
# 1. BASELINE CHARACTERISTICS TABLE (TABLE 1)
# =============================================================================
# Continuous variables: mean (SD). Binary/categorical variables: n (%).

cat("=== OVERALL COHORT ===\n")
cat("Total patients:", nrow(analysis_data), "\n")
cat("Total deaths:  ", sum(analysis_data$died), "\n")
cat("Mortality (%): ", round(100 * mean(analysis_data$died), 1), "%\n\n")

cat("=== BY EXPOSURE GROUP ===\n")
analysis_data %>%
  group_by(ppi) %>%
  summarise(n = n(), deaths = sum(died), pct_died = round(100 * deaths / n, 1),
            .groups = "drop") %>%
  print()

build_table1 <- function(dat) {

  n_overall <- nrow(dat)
  n_ppi     <- sum(dat$ppi == "PPI")
  n_h2ra    <- sum(dat$ppi == "H2RA")

  header <- tibble(Variable = "N", Overall = as.character(n_overall),
                    PPI = as.character(n_ppi), H2RA = as.character(n_h2ra))

  stats_by_group <- function(d, var_name, type) {
    v <- d[[var_name]]
    if (type == "continuous") {
      sprintf("%.1f (%.1f)", mean(v, na.rm = TRUE), sd(v, na.rm = TRUE))
    } else if (type == "binary") {
      sprintf("%d (%.1f%%)", sum(v == 1, na.rm = TRUE), 100 * mean(v == 1, na.rm = TRUE))
    }
  }

  cont_vars   <- c("age", "bmi", "n_contacts")
  cont_labels <- c("Age at index date, years: mean (SD)",
                   "BMI, kg/m²: mean (SD)",
                   "Healthcare contacts (prior year): mean (SD)")

  cont_rows <- map2_dfr(cont_vars, cont_labels, function(v, lbl) {
    tibble(
      Variable = lbl,
      Overall  = stats_by_group(dat, v, "continuous"),
      PPI      = stats_by_group(dat[dat$ppi == "PPI", ], v, "continuous"),
      H2RA     = stats_by_group(dat[dat$ppi == "H2RA", ], v, "continuous")
    )
  })

  fup_row <- tibble(
    Variable = "Follow-up time, years: median [IQR]",
    Overall  = sprintf("%.2f [%.2f–%.2f]", median(dat$fup_days / 365.25),
                       quantile(dat$fup_days / 365.25, 0.25), quantile(dat$fup_days / 365.25, 0.75)),
    PPI      = sprintf("%.2f [%.2f–%.2f]", median(dat$fup_days[dat$ppi == "PPI"] / 365.25),
                       quantile(dat$fup_days[dat$ppi == "PPI"] / 365.25, 0.25),
                       quantile(dat$fup_days[dat$ppi == "PPI"] / 365.25, 0.75)),
    H2RA     = sprintf("%.2f [%.2f–%.2f]", median(dat$fup_days[dat$ppi == "H2RA"] / 365.25),
                       quantile(dat$fup_days[dat$ppi == "H2RA"] / 365.25, 0.25),
                       quantile(dat$fup_days[dat$ppi == "H2RA"] / 365.25, 0.75))
  )

  death_row <- tibble(
    Variable = "All-cause deaths, n (%)",
    Overall  = sprintf("%d (%.1f%%)", sum(dat$died), 100 * mean(dat$died)),
    PPI      = sprintf("%d (%.1f%%)", sum(dat$died[dat$ppi == "PPI"]), 100 * mean(dat$died[dat$ppi == "PPI"])),
    H2RA     = sprintf("%d (%.1f%%)", sum(dat$died[dat$ppi == "H2RA"]), 100 * mean(dat$died[dat$ppi == "H2RA"]))
  )

  sex_row <- tibble(
    Variable = "Female sex, n (%)",
    Overall  = sprintf("%d (%.1f%%)", sum(dat$gender == "Female"), 100 * mean(dat$gender == "Female")),
    PPI      = sprintf("%d (%.1f%%)", sum(dat$gender[dat$ppi == "PPI"] == "Female"),
                       100 * mean(dat$gender[dat$ppi == "PPI"] == "Female")),
    H2RA     = sprintf("%d (%.1f%%)", sum(dat$gender[dat$ppi == "H2RA"] == "Female"),
                       100 * mean(dat$gender[dat$ppi == "H2RA"] == "Female"))
  )

  cat_levels <- function(dat, var, levels_vec, labels_vec, section_label) {
    rows <- map2_dfr(levels_vec, labels_vec, function(lv, lb) {
      tibble(
        Variable = lb,
        Overall  = sprintf("%d (%.1f%%)", sum(dat[[var]] == lv, na.rm = TRUE),
                           100 * mean(dat[[var]] == lv, na.rm = TRUE)),
        PPI      = sprintf("%d (%.1f%%)", sum(dat[[var]][dat$ppi == "PPI"] == lv, na.rm = TRUE),
                           100 * mean(dat[[var]][dat$ppi == "PPI"] == lv, na.rm = TRUE)),
        H2RA     = sprintf("%d (%.1f%%)", sum(dat[[var]][dat$ppi == "H2RA"] == lv, na.rm = TRUE),
                           100 * mean(dat[[var]][dat$ppi == "H2RA"] == lv, na.rm = TRUE))
      )
    })
    bind_rows(tibble(Variable = section_label, Overall = "", PPI = "", H2RA = ""), rows)
  }

  eth_rows <- cat_levels(dat, "ethnicity",
                          c("White", "South-Asian", "Black", "Mixed", "Other"),
                          c("  White", "  South Asian", "  Black", "  Mixed", "  Other"),
                          "Ethnicity, n (%)")

  ses_rows <- cat_levels(dat, "ses_person", c("1", "2", "3", "4", "5"),
                          c("  1 (least deprived)", "  2", "  3", "  4", "  5 (most deprived)"),
                          "Socioeconomic status (SES), n (%)")

  cal_levels <- levels(dat$calendarperiod)
  cal_rows <- cat_levels(dat, "calendarperiod", cal_levels,
                          paste0("  ", cal_levels), "Calendar period of index date, n (%)")

  bin_vars   <- c("prior_gastric_cancer", "recent_gerd", "recent_peptic_ulcer")
  bin_labels <- c("Prior gastric cancer, n (%)",
                  "GERD in 6 months prior to index date, n (%)",
                  "Peptic ulcer in 6 months prior to index date, n (%)")

  bin_rows <- map2_dfr(bin_vars, bin_labels, function(v, lbl) {
    tibble(
      Variable = lbl,
      Overall  = sprintf("%d (%.1f%%)", sum(dat[[v]] == 1, na.rm = TRUE), 100 * mean(dat[[v]] == 1, na.rm = TRUE)),
      PPI      = sprintf("%d (%.1f%%)", sum(dat[[v]][dat$ppi == "PPI"] == 1, na.rm = TRUE),
                         100 * mean(dat[[v]][dat$ppi == "PPI"] == 1, na.rm = TRUE)),
      H2RA     = sprintf("%d (%.1f%%)", sum(dat[[v]][dat$ppi == "H2RA"] == 1, na.rm = TRUE),
                         100 * mean(dat[[v]][dat$ppi == "H2RA"] == 1, na.rm = TRUE))
    )
  })

  bmi_miss_row <- tibble(
    Variable = "BMI missing, n (%)",
    Overall  = sprintf("%d (%.1f%%)", sum(is.na(dat$bmi)), 100 * mean(is.na(dat$bmi))),
    PPI      = sprintf("%d (%.1f%%)", sum(is.na(dat$bmi[dat$ppi == "PPI"])),
                       100 * mean(is.na(dat$bmi[dat$ppi == "PPI"]))),
    H2RA     = sprintf("%d (%.1f%%)", sum(is.na(dat$bmi[dat$ppi == "H2RA"])),
                       100 * mean(is.na(dat$bmi[dat$ppi == "H2RA"])))
  )

  bind_rows(header, fup_row, death_row, cont_rows[1, ], sex_row, eth_rows, ses_rows,
            cal_rows, cont_rows[2, ], bmi_miss_row, bin_rows, cont_rows[3, ])
}

table1 <- build_table1(analysis_data)
print(table1, n = Inf)
write_csv(table1, "results/table1_baseline_characteristics.csv")
cat("\nTable 1 saved to: results/table1_baseline_characteristics.csv\n")


# =============================================================================
# 2. TIME-TO-EVENT SUMMARY BY EXPOSURE GROUP
# =============================================================================

cat("\n=== TIME-TO-EVENT SUMMARY BY EXPOSURE GROUP ===\n")

tte_summary <- analysis_data %>%
  group_by(ppi) %>%
  summarise(
    n = n(), n_events = sum(died), event_rate_pct = round(100 * n_events / n, 1),
    person_years = round(sum(fup_days) / 365.25, 0),
    crude_rate_per_1000py = round(1000 * n_events / (sum(fup_days) / 365.25), 2),
    median_fup_yrs = round(median(fup_days / 365.25), 2),
    .groups = "drop"
  )
print(tte_summary)


# =============================================================================
# 3. KAPLAN-MEIER ANALYSIS
# =============================================================================

km_fit <- survfit(Surv(fup_years, died) ~ ppi, data = analysis_data, conf.type = "log")

cat("\n=== KAPLAN-MEIER SURVIVAL ESTIMATES AT KEY TIME POINTS ===\n")
print(summary(km_fit, times = c(1, 2, 5, 10)))

logrank_test <- survdiff(Surv(fup_years, died) ~ ppi, data = analysis_data)
cat("\n=== LOG-RANK TEST ===\n")
print(logrank_test)

km_plot <- ggsurvplot(
  km_fit, data = analysis_data, conf.int = TRUE, pval = TRUE, pval.method = TRUE,
  risk.table = TRUE, risk.table.height = 0.25,
  xlab = "Time since first prescription (years)", ylab = "Survival probability",
  legend.title = "Exposure", legend.labs = c("H2RA", "PPI"),
  palette = c("#E69F00", "#0072B2"),   # colour-blind-friendly
  ggtheme = theme_bw(base_size = 12),
  xlim = c(0, max(analysis_data$fup_years)), break.time.by = 2,
  surv.median.line = "hv", tables.theme = theme_cleantable(),
  title = "Kaplan-Meier estimates of survival by exposure group"
)

pdf("results/km_plot.pdf", width = 10, height = 7)
print(km_plot)
dev.off()

png("results/km_plot.png", width = 2000, height = 1400, res = 200)
print(km_plot)
dev.off()

cat("\nKaplan-Meier plot saved to: results/km_plot.pdf and .png\n")


# =============================================================================
# 4. RE-SAVE THE ANALYSIS DATASET
# =============================================================================
# Script 03 relies on gender/ethnicity/ses_person being factors (set in
# section 0 above), so the dataset is saved back out here.

saveRDS(analysis_data, file = "data/processed/analysis_data.Rds")

# =============================================================================
# END OF SCRIPT 2
# =============================================================================
