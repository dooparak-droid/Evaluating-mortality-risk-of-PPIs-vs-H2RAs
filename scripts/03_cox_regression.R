# =============================================================================
# EHR COHORT ANALYSIS PIPELINE — DEMONSTRATION
# Script 3: Cox Regression, Model Assessment, and Figures
#
# Produces:
#   Univariable and multivariable Cox models for PPI vs H2RA
#   Functional-form checks (Martingale residuals) for continuous confounders
#   Proportional-hazards checks (Schoenfeld residuals, log-log plot)
#   Predicted survivor curves for four illustrative patient profiles
#   A cohort-flow diagram, drawn from the counts script 01 recorded
#
# Inputs: data/processed/analysis_data.Rds, data/processed/flow_counts.Rds
# Outputs: results/table2_cox_results.csv and several figures under results/
# =============================================================================


# =============================================================================
# 0. PACKAGES AND DATA
# =============================================================================

library(tidyverse)
library(survival)
library(survminer)   # ggcoxzph()
library(splines)     # ns() for natural cubic splines
library(broom)       # tidy() for clean model output
library(grid)        # for the cohort-flow diagram

analysis_data <- readRDS("data/processed/analysis_data.Rds")
flow          <- readRDS("data/processed/flow_counts.Rds")

model_vars <- c("fup_years", "died", "ppi", "age", "gender", "ethnicity",
                "ses_person", "calendarperiod", "bmi", "prior_gastric_cancer",
                "recent_gerd", "recent_peptic_ulcer", "n_contacts")

analysis_cc <- analysis_data %>% drop_na(all_of(model_vars))

cat("Original dataset:", nrow(analysis_data), "patients\n")
cat("Complete-case dataset:", nrow(analysis_cc), "patients\n")


# =============================================================================
# 1. FUNCTIONAL FORM OF CONTINUOUS CONFOUNDERS
# =============================================================================
# Method: fit a Cox model without the variable of interest, extract
# Martingale residuals, and plot against the variable with a lowess
# smoother. A linear relationship implies a linear term is appropriate;
# curvature suggests a transformation (e.g. a spline) is needed.

cox_no_age <- coxph(
  Surv(fup_years, died) ~ ppi + gender + ethnicity + ses_person +
    calendarperiod + bmi + prior_gastric_cancer + recent_gerd +
    recent_peptic_ulcer + n_contacts,
  data = analysis_cc
)
mart_no_age <- residuals(cox_no_age, type = "martingale")

cox_no_contacts <- coxph(
  Surv(fup_years, died) ~ ppi + age + gender + ethnicity + ses_person +
    calendarperiod + bmi + prior_gastric_cancer + recent_gerd +
    recent_peptic_ulcer,
  data = analysis_cc
)
mart_no_contacts <- residuals(cox_no_contacts, type = "martingale")

pdf("results/fig_martingale_functional_form.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
scatter.smooth(analysis_cc$age, mart_no_age,
                xlab = "Age at index date (years)", ylab = "Martingale residuals",
                main = "Functional form: age", col = rgb(0, 0, 0, 0.15), pch = 16, cex = 0.5,
                lpars = list(col = "#D55E00", lwd = 2))
abline(h = 0, lty = 2, col = "grey50")
scatter.smooth(analysis_cc$n_contacts, mart_no_contacts,
                xlab = "Healthcare contacts in prior year", ylab = "Martingale residuals",
                main = "Functional form: healthcare contacts", col = rgb(0, 0, 0, 0.15),
                pch = 16, cex = 0.5, lpars = list(col = "#D55E00", lwd = 2))
abline(h = 0, lty = 2, col = "grey50")
par(mfrow = c(1, 1))
dev.off()
cat("Martingale residual plot saved to: results/fig_martingale_functional_form.pdf\n")

# Inspect the plot above before finalising the model formula below.
# A natural cubic spline (ns(x, df = 3)) is used by default for both age
# and n_contacts here, as a flexible choice that will not overfit; replace
# with the raw variable for either one if its own smoother looks linear.


# =============================================================================
# 2. UNIVARIABLE COX MODEL
# =============================================================================

cox_uni <- coxph(Surv(fup_years, died) ~ ppi, data = analysis_cc)
cat("\n=== UNIVARIABLE COX MODEL ===\n")
print(summary(cox_uni))

uni_tidy <- tidy(cox_uni, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(model = "Univariable") %>%
  select(model, term, estimate, conf.low, conf.high, p.value)


# =============================================================================
# 3. MULTIVARIABLE COX MODEL
# =============================================================================
# Confounders: age, sex, ethnicity, deprivation, BMI (linear), calendar
# period, prior gastric cancer, recent GERD, recent peptic ulcer, and
# healthcare contacts. No interaction terms; no propensity-score analysis.

cox_mv <- coxph(
  Surv(fup_years, died) ~
    ppi + ns(age, df = 3) + gender + ethnicity + ses_person + bmi +
    calendarperiod + prior_gastric_cancer + recent_gerd +
    recent_peptic_ulcer + ns(n_contacts, df = 3),
  data = analysis_cc,
  ties = "efron"
)

cat("\n=== MULTIVARIABLE COX MODEL ===\n")
print(summary(cox_mv))

mv_tidy <- tidy(cox_mv, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(model = "Multivariable") %>%
  select(model, term, estimate, conf.low, conf.high, p.value)


# =============================================================================
# 4. TABLE 2 — EXPOSURE EFFECT FROM CRUDE AND MULTIVARIABLE MODELS
# =============================================================================

table2 <- bind_rows(uni_tidy, mv_tidy) %>%
  filter(grepl("ppi", term)) %>%
  mutate(
    HR = round(estimate, 2), CI_lower = round(conf.low, 2), CI_upper = round(conf.high, 2),
    p_value = ifelse(p.value < 0.001, "<0.001", round(p.value, 3)),
    HR_CI = paste0(HR, " (", CI_lower, "–", CI_upper, ")")
  ) %>%
  select(Model = model, `HR (95% CI)` = HR_CI, `p-value` = p_value)

cat("\n=== TABLE 2: EXPOSURE EFFECT — PPI vs H2RA ===\n")
print(table2)
write_csv(table2, "results/table2_cox_results.csv")


# =============================================================================
# 5. PROPORTIONAL HAZARDS ASSUMPTION
# =============================================================================
# Two complementary checks:
#   (A) Schoenfeld residuals via cox.zph() — a significant p-value
#       (conventionally < 0.05) suggests non-proportional hazards.
#   (B) A log(-log S(t)) plot for the exposure variable — parallel lines
#       support proportional hazards.
#
# If a violation is found for the exposure itself, this materially
# affects the causal estimate; a violation confined to a confounder can
# usually be handled by stratifying the baseline hazard on that variable.

ph_test <- cox.zph(cox_mv, transform = "identity")
cat("\n=== SCHOENFELD RESIDUAL TEST — PROPORTIONAL HAZARDS ===\n")
print(ph_test)

pdf("results/fig_schoenfeld_residuals.pdf", width = 12, height = 16)
par(mfrow = c(4, 3))
plot(ph_test)
par(mfrow = c(1, 1))
dev.off()

ggcoxzph_plot <- ggcoxzph(ph_test, var = "ppi", font.main = 12,
                           caption = "Scaled Schoenfeld residuals for PPI exposure")
ggsave("results/fig_schoenfeld_ppi.pdf", ggcoxzph_plot[[1]], width = 7, height = 5)

km_by_ppi <- survfit(Surv(fup_years, died) ~ ppi, data = analysis_cc)
pdf("results/fig_loglog_check.pdf", width = 7, height = 5)
plot(km_by_ppi, fun = "cloglog", col = c("#E69F00", "#0072B2"), lwd = 2,
     xlab = "log(time)", ylab = "log(-log S(t))",
     main = "Log(-log) plot: checking proportional hazards for PPI vs H2RA")
legend("topleft", legend = c("H2RA", "PPI"), col = c("#E69F00", "#0072B2"), lwd = 2)
dev.off()

# --- Decide on a final model based on the test above -------------------------
# Below are the two standard responses to a proportional-hazards violation.
# Pick whichever applies to your own cox.zph() output.

# Option A — no meaningful violation anywhere: use cox_mv unmodified.
# cox_final <- cox_mv

# Option B — violation confined to a confounder (shown here for
# calendarperiod): stratify on it, which allows the baseline hazard to
# vary by period while keeping a single hazard ratio for ppi.
cox_final <- coxph(
  Surv(fup_years, died) ~
    ppi + ns(age, df = 3) + gender + ethnicity + ses_person + bmi +
    strata(calendarperiod) + prior_gastric_cancer + recent_gerd +
    recent_peptic_ulcer + ns(n_contacts, df = 3),
  data = analysis_cc,
  ties = "efron"
)

# Option C — violation found for the exposure itself: split follow-up at
# a chosen time point and allow the hazard ratio to differ before and
# after it (uncomment and adapt if this applies):
#
# analysis_split <- survSplit(
#   Surv(fup_years, died) ~ ., cut = 5, data = analysis_cc,
#   start = "t_start", end = "t_stop", event = "died"
# ) %>%
#   mutate(ppi_early = ifelse(ppi == "PPI" & t_stop <= 5, 1, 0),
#          ppi_late  = ifelse(ppi == "PPI" & t_stop >  5, 1, 0))
# cox_final <- coxph(
#   Surv(t_start, t_stop, died) ~ ppi_early + ppi_late + ns(age, df = 3) +
#     gender + ethnicity + ses_person + bmi + calendarperiod +
#     prior_gastric_cancer + recent_gerd + recent_peptic_ulcer +
#     ns(n_contacts, df = 3),
#   data = analysis_split, ties = "efron"
# )

final_res <- tidy(cox_final, exponentiate = TRUE, conf.int = TRUE)
ppi_hr <- final_res %>% filter(term == "ppiPPI")

cat(sprintf("\nAdjusted Hazard Ratio (PPI vs H2RA): %.2f (95%% CI %.2f-%.2f), p = %.4f\n",
            ppi_hr$estimate, ppi_hr$conf.low, ppi_hr$conf.high, ppi_hr$p.value))

ph_test_final <- cox.zph(cox_final)
cat("\n--- Proportional hazards test (final model) ---\n")
print(ph_test_final$table)


# =============================================================================
# 6. ESTIMATED SURVIVOR CURVES FOR FOUR ILLUSTRATIVE PROFILES
# =============================================================================
# Four covariate profiles, crossing sex with a comorbidity/deprivation
# contrast, each shown under both PPI and H2RA:
#   (a) Man,   no comorbidities,     most deprived
#   (b) Man,   recent peptic ulcer,  least deprived
#   (c) Woman, no comorbidities,     most deprived
#   (d) Woman, recent peptic ulcer,  least deprived

cox_for_curves <- cox_mv

ref_bmi      <- median(analysis_cc$bmi, na.rm = TRUE)
ref_contacts <- median(analysis_cc$n_contacts, na.rm = TRUE)
ref_cal      <- levels(analysis_cc$calendarperiod)[
  ceiling(length(levels(analysis_cc$calendarperiod)) / 2)
]
ref_cal      <- factor(ref_cal, levels = levels(analysis_cc$calendarperiod))

base_profiles <- data.frame(
  profile_label = c("(a) Man, no comorbidities, most deprived",
                    "(b) Man, peptic ulcer, least deprived",
                    "(c) Woman, no comorbidities, most deprived",
                    "(d) Woman, peptic ulcer, least deprived"),
  age = 60,
  gender = factor(c("Male", "Male", "Female", "Female"), levels = levels(analysis_cc$gender)),
  ethnicity = factor("White", levels = levels(analysis_cc$ethnicity)),
  ses_person = factor(c("5", "1", "5", "1"), levels = levels(analysis_cc$ses_person)),
  bmi = ref_bmi,
  calendarperiod = ref_cal,
  prior_gastric_cancer = 0,
  recent_gerd = 0,
  recent_peptic_ulcer = c(0, 1, 0, 1),
  n_contacts = ref_contacts,
  ppi = factor("PPI", levels = levels(analysis_cc$ppi)),
  stringsAsFactors = FALSE
)

h2ra_profiles <- base_profiles %>% mutate(ppi = factor("H2RA", levels = levels(analysis_cc$ppi)))

# Interleave so that odd columns of surv_curves are PPI, even are H2RA
all_profiles <- bind_rows(lapply(1:4, function(i) bind_rows(base_profiles[i, ], h2ra_profiles[i, ])))

surv_curves <- survfit(cox_for_curves, newdata = all_profiles)

ppi_col  <- "#0072B2"
h2ra_col <- "#E69F00"

plot_survivor_curves <- function() {
  par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))
  for (i in 1:4) {
    idx_ppi  <- (2 * i) - 1
    idx_h2ra <- 2 * i
    plot(surv_curves$time, surv_curves$surv[, idx_ppi], type = "s", col = ppi_col, lwd = 2,
         ylim = c(0, 1), xlim = c(0, max(surv_curves$time)),
         xlab = "Time since first prescription (years)", ylab = "Estimated survival probability",
         main = base_profiles$profile_label[i], cex.main = 0.8)
    lines(surv_curves$time, surv_curves$surv[, idx_h2ra], type = "s", col = h2ra_col, lwd = 2)
    lines(surv_curves$time, surv_curves$lower[, idx_ppi], lty = 2, col = ppi_col)
    lines(surv_curves$time, surv_curves$upper[, idx_ppi], lty = 2, col = ppi_col)
    lines(surv_curves$time, surv_curves$lower[, idx_h2ra], lty = 2, col = h2ra_col)
    lines(surv_curves$time, surv_curves$upper[, idx_h2ra], lty = 2, col = h2ra_col)
    legend("bottomleft", legend = c("PPI", "H2RA"), col = c(ppi_col, h2ra_col), lwd = 2, bty = "n")
  }
}

pdf("results/fig_survivor_curves.pdf", width = 12, height = 10)
plot_survivor_curves()
dev.off()

png("results/fig_survivor_curves.png", width = 2000, height = 1700, res = 200)
plot_survivor_curves()
dev.off()

cat("Survivor curves saved to: results/fig_survivor_curves.pdf and .png\n")


# =============================================================================
# 7. COHORT-FLOW DIAGRAM
# =============================================================================
# Drawn from the counts script 01 recorded in data/processed/flow_counts.Rds,
# rather than numbers copied in by hand.

n_ppi_final  <- sum(analysis_cc$ppi == "PPI")
n_h2ra_final <- sum(analysis_cc$ppi == "H2RA")
fmt <- function(n) formatC(n, format = "d", big.mark = ",")

col_box <- "#1F3A5F"; col_excl <- "#8B1A1A"; col_alloc <- "#1F5F3A"
col_bg <- "#F5F8FC"; col_excl_bg <- "#FDF0F0"; col_alloc_bg <- "#F0FDF4"

draw_box <- function(x, y, w, h, label, sublabel = NULL, col_border = col_box,
                      col_fill = col_bg, fontsize = 9, bold = FALSE) {
  grid.roundrect(x = unit(x, "npc"), y = unit(y, "npc"), width = unit(w, "npc"),
                 height = unit(h, "npc"), r = unit(3, "mm"),
                 gp = gpar(fill = col_fill, col = col_border, lwd = 1.5))
  ytext <- if (is.null(sublabel)) y else y + h * 0.13
  grid.text(label, x = unit(x, "npc"), y = unit(ytext, "npc"),
            gp = gpar(fontsize = fontsize, fontface = if (bold) "bold" else "plain", col = col_border),
            just = "centre")
  if (!is.null(sublabel)) {
    grid.text(sublabel, x = unit(x, "npc"), y = unit(y - h * 0.13, "npc"),
              gp = gpar(fontsize = fontsize - 1, col = col_border), just = "centre")
  }
}

draw_excl_box <- function(x, y, w, h, lines) {
  grid.roundrect(x = unit(x, "npc"), y = unit(y, "npc"), width = unit(w, "npc"),
                 height = unit(h, "npc"), r = unit(2, "mm"),
                 gp = gpar(fill = col_excl_bg, col = col_excl, lwd = 1.2))
  n <- length(lines); step <- h / (n + 1)
  for (i in seq_along(lines)) {
    grid.text(lines[i], x = unit(x - w * 0.42, "npc"), y = unit(y + h / 2 - i * step, "npc"),
              just = c("left", "centre"), gp = gpar(fontsize = 7.8, col = col_excl))
  }
}

arrow_down <- function(x, y_from, y_to) {
  grid.lines(x = unit(c(x, x), "npc"), y = unit(c(y_from, y_to), "npc"),
             gp = gpar(col = "grey40", lwd = 1.3),
             arrow = arrow(length = unit(2.5, "mm"), type = "closed", angle = 20))
}

arrow_right <- function(x_from, x_to, y) {
  grid.lines(x = unit(c(x_from, x_to), "npc"), y = unit(c(y, y), "npc"),
             gp = gpar(col = col_excl, lwd = 1.2),
             arrow = arrow(length = unit(2, "mm"), type = "closed", angle = 20))
}

draw_flowchart <- function() {
grid.newpage()

cx <- 0.38; bw <- 0.46; bh <- 0.062
ex_cx <- 0.80; ex_bw <- 0.34; ex_bh <- 0.058
y1 <- 0.935; y2 <- 0.820; y3 <- 0.705; y4 <- 0.590; y5 <- 0.475; y6 <- 0.345; y7 <- 0.185

draw_box(cx, y1, bw, bh, paste0("Patients with a PPI or H2RA prescription  (n = ", fmt(flow$n_start), ")"),
         bold = TRUE)
draw_box(cx, y2, bw, bh, paste0("After excluding missing demographics  (n = ", fmt(flow$n_demog), ")"))
draw_box(cx, y3, bw, bh, paste0("After excluding < 1 year registration run-in  (n = ", fmt(flow$n_runin), ")"))
draw_box(cx, y4, bw, bh, paste0("After excluding outside study window  (n = ", fmt(flow$n_dates), ")"))
draw_box(cx, y5, bw, bh, paste0("After excluding end date before index date  (n = ", fmt(flow$n_fup), ")"))
draw_box(cx, y6, bw, bh, paste0("Final analysis cohort (complete case)  (n = ", fmt(flow$n_complete), ")"),
         col_fill = "#E8EFF8", bold = TRUE)

alloc_left_cx <- 0.19; alloc_right_cx <- 0.57; alloc_bw <- 0.30; alloc_bh <- 0.09

draw_box(alloc_left_cx, y7, alloc_bw, alloc_bh, "H2RA group", sublabel = paste0("n = ", fmt(n_h2ra_final)),
         col_border = col_alloc, col_fill = col_alloc_bg, bold = TRUE)
draw_box(alloc_right_cx, y7, alloc_bw, alloc_bh, "PPI group", sublabel = paste0("n = ", fmt(n_ppi_final)),
         col_border = col_alloc, col_fill = col_alloc_bg, bold = TRUE)

arrow_down(cx, y1 - bh / 2, y2 + bh / 2)
arrow_down(cx, y2 - bh / 2, y3 + bh / 2)
arrow_down(cx, y3 - bh / 2, y4 + bh / 2)
arrow_down(cx, y4 - bh / 2, y5 + bh / 2)
arrow_down(cx, y5 - bh / 2, y6 + bh / 2)

fork_y <- y6 - bh / 2 - 0.03
grid.lines(x = unit(c(cx, cx), "npc"), y = unit(c(y6 - bh / 2, fork_y), "npc"), gp = gpar(col = "grey40", lwd = 1.3))
grid.lines(x = unit(c(alloc_left_cx, alloc_right_cx), "npc"), y = unit(c(fork_y, fork_y), "npc"),
           gp = gpar(col = "grey40", lwd = 1.3))
grid.lines(x = unit(c(alloc_left_cx, alloc_left_cx), "npc"), y = unit(c(fork_y, y7 + alloc_bh / 2), "npc"),
           gp = gpar(col = col_alloc, lwd = 1.3), arrow = arrow(length = unit(2.5, "mm"), type = "closed", angle = 20))
grid.lines(x = unit(c(alloc_right_cx, alloc_right_cx), "npc"), y = unit(c(fork_y, y7 + alloc_bh / 2), "npc"),
           gp = gpar(col = col_alloc, lwd = 1.3), arrow = arrow(length = unit(2.5, "mm"), type = "closed", angle = 20))

excl_x_right <- cx + bw / 2
excl_mid_1 <- (y1 + y2) / 2
draw_excl_box(ex_cx, excl_mid_1, ex_bw, ex_bh,
              c("Missing year of birth, sex, or", paste0("ethnicity: n = ", fmt(flow$n_start - flow$n_demog))))
arrow_right(excl_x_right, ex_cx - ex_bw / 2, excl_mid_1)

excl_mid_2 <- (y2 + y3) / 2
draw_excl_box(ex_cx, excl_mid_2, ex_bw, ex_bh,
              c("Less than 1 year registration", paste0("prior to index date: n = ", fmt(flow$n_demog - flow$n_runin))))
arrow_right(excl_x_right, ex_cx - ex_bw / 2, excl_mid_2)

excl_mid_3 <- (y3 + y4) / 2
draw_excl_box(ex_cx, excl_mid_3, ex_bw, ex_bh,
              c("Index date outside study window", paste0("n = ", fmt(flow$n_runin - flow$n_dates))))
arrow_right(excl_x_right, ex_cx - ex_bw / 2, excl_mid_3)

excl_mid_4 <- (y4 + y5) / 2
draw_excl_box(ex_cx, excl_mid_4, ex_bw, ex_bh,
              c("End of follow-up before", paste0("index date: n = ", fmt(flow$n_dates - flow$n_fup))))
arrow_right(excl_x_right, ex_cx - ex_bw / 2, excl_mid_4)

excl_mid_5 <- (y5 + y6) / 2
draw_excl_box(ex_cx, excl_mid_5, ex_bw, ex_bh,
              c("Missing a value needed for the", paste0("model (complete-case): n = ", fmt(flow$n_fup - flow$n_complete))))
arrow_right(excl_x_right, ex_cx - ex_bw / 2, excl_mid_5)

grid.text("Study Population Flowchart", x = unit(0.5, "npc"), y = unit(0.99, "npc"),
          just = c("centre", "top"), gp = gpar(fontsize = 12, fontface = "bold", col = col_box))
}

pdf("results/fig_flowchart.pdf", width = 9, height = 10)
draw_flowchart()
dev.off()

png("results/fig_flowchart.png", width = 1800, height = 2000, res = 200)
draw_flowchart()
dev.off()

cat("Flowchart saved to: results/fig_flowchart.pdf and .png\n")


# =============================================================================
# 8. SUMMARY OF KEY RESULTS
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("SUMMARY OF KEY NUMBERS\n")
cat(strrep("=", 60), "\n")
cat("\nComplete-case N:", nrow(analysis_cc), "\n")
cat("Events (deaths):", sum(analysis_cc$died), "\n")

uni_hr <- exp(coef(cox_uni)["ppiPPI"])
uni_ci <- exp(confint(cox_uni)["ppiPPI", ])
cat(sprintf("\nUnivariable HR (PPI vs H2RA) = %.2f (95%% CI %.2f-%.2f)\n", uni_hr, uni_ci[1], uni_ci[2]))

mv_hr <- exp(coef(cox_final)["ppiPPI"])
mv_ci <- exp(confint(cox_final)["ppiPPI", ])
cat(sprintf("Final-model HR (PPI vs H2RA)  = %.2f (95%% CI %.2f-%.2f)\n", mv_hr, mv_ci[1], mv_ci[2]))

cat(sprintf("\nConcordance (final model) = %.3f (SE = %.3f)\n",
            summary(cox_final)$concordance[1], summary(cox_final)$concordance[2]))

# =============================================================================
# END OF SCRIPT 3
# =============================================================================
