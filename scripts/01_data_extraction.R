# =============================================================================
# EHR COHORT ANALYSIS PIPELINE — DEMONSTRATION
# Script 1: Data Extraction
#
# Worked example: does a prescription for a proton pump inhibitor (PPI)
# affect subsequent all-cause mortality, compared with a histamine H2-
# receptor antagonist (H2RA)? Run against synthetic data (see
# 00_generate_synthetic_data.R) — none of the numbers this produces
# describe anything real.
#
# Study design: cohort study with an active-comparator design.
# Index date:   date of first PPI or H2RA prescription (within eligibility
#               window).
# Follow-up:    index date to the first of: death, leaving the practice,
#               or STUDY_END below.
# Outcome:      all-cause mortality.
#
# This script reads the raw data files, applies eligibility criteria,
# derives confounders, and saves the analysis-ready dataset. It also
# records the size of the cohort after each exclusion step, so that
# script 03 can draw a flowchart from real numbers rather than counts
# copied in by hand.
# =============================================================================


# =============================================================================
# 0. PACKAGES, CONFIG, AND DATA IMPORT
# =============================================================================

library(tidyverse)
library(survival)   # Needed for Surv() in later scripts

# Study parameters — change these to fit a different data extract.
STUDY_START             <- as.Date("2000-01-01")
STUDY_END               <- as.Date("2019-01-01")
MIN_REGISTRATION_YEARS  <- 1

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# --- Read raw data files -----------------------------------------------------
data_files <- list.files("data/raw", full.names = TRUE, pattern = "\\.csv$")
data_names <- tolower(tools::file_path_sans_ext(basename(data_files)))
data <- lapply(data_files, read_csv, show_col_types = FALSE)
names(data) <- data_names

code_files <- list.files("data/codelists", full.names = TRUE, pattern = "\\.csv$")
code_names <- tolower(tools::file_path_sans_ext(basename(code_files)))
codes <- lapply(code_files, read_csv, show_col_types = FALSE)
names(codes) <- code_names

stopifnot(all(c("patient", "prescription", "clinic", "other",
                "ses", "death", "contactdates") %in% names(data)))


# =============================================================================
# 1. DEFINE THE ELIGIBLE COHORT AND INDEX DATE
# =============================================================================
# Eligibility criteria applied here:
#   Recorded year of birth, sex and ethnicity
#   First prescription after registration + MIN_REGISTRATION_YEARS
#   First prescription within [STUDY_START, STUDY_END]
#   Follow-up has not ended before the index date

flow <- list()  # tracks cohort size after each step, for the flowchart

# --- Step 1a: Extract all PPI and H2RA prescriptions ------------------------
ppis <- data$prescription %>%
  inner_join(codes$ppi, by = "prescode") %>%
  mutate(ppi = 1)

h2ras <- data$prescription %>%
  inner_join(codes$h2ra, by = "prescode") %>%
  mutate(ppi = 0)

all_ppi_h2ra <- bind_rows(ppis, h2ras)

# --- Step 1b: Identify each patient's FIRST PPI or H2RA prescription --------
# Where a patient has both on the same day, the PPI is retained as the
# exposure.
first_rx <- all_ppi_h2ra %>%
  arrange(id, eventdate, desc(ppi)) %>%
  filter(!duplicated(id)) %>%
  select(id, indexdate = eventdate, ppi)

flow$n_start <- nrow(first_rx)
cat("Patients with a PPI or H2RA prescription:", flow$n_start, "\n")

# --- Step 1c: Join to patient file to apply remaining eligibility criteria ---
cohort <- first_rx %>%
  left_join(
    data$patient %>% select(id, yob, gender, eth5, start, leavepractice, accept),
    by = "id"
  ) %>%
  left_join(data$death %>% filter(death == 1) %>% select(id, death_date = deathdate),
    by = "id")

cohort <- cohort %>%
  filter(!is.na(yob), !is.na(gender), !is.na(eth5), eth5 != "Unknown")

flow$n_demog <- nrow(cohort)
cat("After requiring recorded YOB, sex, and non-missing ethnicity:",
    flow$n_demog, "(excluded:", flow$n_start - flow$n_demog, ")\n")

cohort <- cohort %>%
  filter(indexdate > start + 365.25 * MIN_REGISTRATION_YEARS)

flow$n_runin <- nrow(cohort)
cat("After requiring >=", MIN_REGISTRATION_YEARS, "year(s) registration run-in:",
    flow$n_runin, "(excluded:", flow$n_demog - flow$n_runin, ")\n")

cohort <- cohort %>%
  filter(indexdate >= STUDY_START, indexdate <= STUDY_END)

flow$n_dates <- nrow(cohort)
cat("After applying study date window (", format(STUDY_START), "to", format(STUDY_END),
    "):", flow$n_dates, "(excluded:", flow$n_runin - flow$n_dates, ")\n")

cohort <- cohort %>%
  mutate(
    potential_end = pmin(STUDY_END, leavepractice, death_date, na.rm = TRUE)) %>%
  filter(potential_end >= indexdate)

flow$n_fup <- nrow(cohort)
cat("Excluded due to follow-up ending before index date:",
    flow$n_dates - flow$n_fup, "\n")

cohort <- cohort %>%
  select(id, indexdate, ppi)

cat("Eligible cohort size:", nrow(cohort), "\n")


# =============================================================================
# 2. DEFINE OUTCOME AND FOLLOW-UP TIMES
# =============================================================================
# Follow-up ends at the earliest of: date of death, date of leaving the
# practice, or STUDY_END.

follow_up <- cohort %>%
  left_join(data$patient %>% select(id, leavepractice), by = "id") %>%
  left_join(data$death %>% filter(death == 1) %>% select(id, death_date = deathdate),
             by = "id") %>%
  mutate(
    enddate  = pmin(STUDY_END, leavepractice, death_date, na.rm = TRUE),
    died     = as.integer(!is.na(death_date) & death_date == enddate),
    fup_days = as.numeric(enddate - indexdate),
    fup_days = ifelse(fup_days == 0, 0.5, fup_days),   # allow zero-day follow-up in Cox model
    fup_years = fup_days / 365.25
  ) %>%
  select(id, enddate, died, fup_days, fup_years)


# =============================================================================
# 3. DERIVE CONFOUNDERS
# =============================================================================

# --- 3a: Demographics: age, sex, ethnicity, deprivation, calendar year -------
# Date of birth is recorded as year only; 15 June in that year is used to
# approximate age.

n_periods <- 4
period_breaks <- seq(STUDY_START, STUDY_END, length.out = n_periods + 1)

demographics <- cohort %>%
  left_join(data$patient %>% select(id, yob, gender, eth5), by = "id") %>%
  left_join(data$ses %>% select(id, ses_person), by = "id") %>%
  mutate(
    dob_approx = as.Date(paste0(yob, "-06-15")),
    age = as.numeric(indexdate - dob_approx) / 365.25,
    calendarperiod = cut(as.numeric(indexdate), breaks = as.numeric(period_breaks),
                          labels = paste0(format(period_breaks[-length(period_breaks)], "%Y"),
                                           "-", format(period_breaks[-1] - 1, "%Y")),
                          include.lowest = TRUE),
    ethnicity = factor(eth5, levels = c("White", "South-Asian", "Black", "Mixed", "Other")),
    ses_person = case_when(
      ses_person == "Least Deprived (1)" ~ "1",
      ses_person == "Most Deprived (5)"  ~ "5",
      TRUE ~ as.character(ses_person)
    ),
    ses_person = factor(ses_person, levels = c("1", "2", "3", "4", "5"))
  ) %>%
  select(id, age, gender, ethnicity, ses_person, calendarperiod)


# --- 3b: BMI (most recent measurement up to 5 years before index date) ------
# BMI is preferentially calculated from weight and height; a direct BMI
# entry is used as a fallback. Implausible values are excluded.

yob_lookup <- data$patient %>% select(id, yob)

bmi_direct <- data$clinic %>%
  filter(entity == 17) %>%
  rename(bmi_date = eventdate) %>%
  inner_join(cohort, by = "id") %>%
  inner_join(data$other, by = c("id", "entity", "otherid")) %>%
  rename(bmi = data3) %>%
  filter(
    !is.na(bmi),
    between(bmi, 5, 200),
    between(as.numeric(indexdate - bmi_date), 0, 5 * 365.25)
  ) %>%
  group_by(id, bmi_date) %>%
  summarise(bmi = mean(bmi), .groups = "drop") %>%
  arrange(id, bmi_date) %>%
  group_by(id) %>%
  filter(bmi_date == max(bmi_date)) %>%
  mutate(preference = 2) %>%
  select(id, bmi, bmi_date, preference)

height_data <- data$clinic %>%
  filter(entity == 18) %>%
  rename(height_date = eventdate) %>%
  inner_join(cohort, by = "id") %>%
  inner_join(data$other, by = c("id", "entity", "otherid")) %>%
  rename(height_m = data1) %>%
  inner_join(yob_lookup, by = "id") %>%
  mutate(year_of_height = as.numeric(format(height_date, "%Y"))) %>%
  filter(
    year_of_height - yob >= 17,
    !is.na(height_m),
    between(height_m, 1.20, 2.15)
  ) %>%
  group_by(id, height_date) %>%
  summarise(height_m = mean(height_m), .groups = "drop") %>%
  arrange(id, height_date) %>%
  group_by(id) %>%
  filter(height_date == max(height_date)) %>%
  select(id, height_m, height_date)

weight_data <- data$clinic %>%
  filter(entity == 17) %>%
  rename(weight_date = eventdate) %>%
  inner_join(cohort, by = "id") %>%
  inner_join(data$other, by = c("id", "entity", "otherid")) %>%
  rename(weight_kg = data1) %>%
  filter(
    !is.na(weight_kg),
    weight_kg >= 20,
    between(as.numeric(indexdate - weight_date), 0, 5 * 365.25)
  ) %>%
  group_by(id, weight_date) %>%
  summarise(weight_kg = mean(weight_kg), .groups = "drop") %>%
  arrange(id, weight_date) %>%
  group_by(id) %>%
  filter(weight_date == max(weight_date)) %>%
  select(id, weight_kg, weight_date)

bmi_calculated <- inner_join(height_data, weight_data, by = "id") %>%
  mutate(bmi = weight_kg / height_m^2, preference = 1) %>%
  rename(bmi_date = weight_date) %>%
  select(id, bmi, bmi_date, preference)

bmi_final <- bind_rows(bmi_calculated, bmi_direct) %>%
  arrange(id, preference) %>%
  filter(!duplicated(id)) %>%
  select(id, bmi)

cat("Patients with BMI measurement:", nrow(bmi_final), "of", nrow(cohort),
    "(missing:", nrow(cohort) - nrow(bmi_final), ")\n")


# --- 3c: Prior gastric cancer -------------------------------------------------

prior_gastric_cancer <- data$clinic %>%
  inner_join(codes$gastric_cancer, by = "clincode") %>%
  arrange(id, eventdate) %>%
  filter(!duplicated(id)) %>%
  inner_join(cohort, by = "id") %>%
  filter(eventdate < indexdate) %>%
  mutate(prior_gastric_cancer = 1) %>%
  select(id, prior_gastric_cancer)


# --- 3d: GERD in the 6 months before the index date ---------------------------

recent_gerd <- data$clinic %>%
  inner_join(cohort, by = "id") %>%
  inner_join(codes$gerd, by = "clincode") %>%
  filter(eventdate >= indexdate - 180, eventdate < indexdate) %>%
  arrange(id, eventdate) %>%
  filter(!duplicated(id)) %>%
  mutate(recent_gerd = 1) %>%
  select(id, recent_gerd)


# --- 3e: Peptic ulcer in the 6 months before the index date -------------------

recent_peptic_ulcer <- data$clinic %>%
  inner_join(cohort, by = "id") %>%
  inner_join(codes$peptic_ulcer, by = "clincode") %>%
  filter(eventdate >= indexdate - 180, eventdate < indexdate) %>%
  arrange(id, eventdate) %>%
  filter(!duplicated(id)) %>%
  mutate(recent_peptic_ulcer = 1) %>%
  select(id, recent_peptic_ulcer)


# --- 3f: Number of healthcare contacts in the year before the index date -----

relevant_contact_types <- c(
  "Surgery consultation", "Followup/routine visit", "Clinic",
  "Telephone call from a patient", "Acute visit", "Home Visit",
  "Emergency Consultation"
)

n_contacts <- data$contactdates %>%
  inner_join(cohort, by = "id") %>%
  filter(eventdate >= indexdate - 365.25, eventdate < indexdate,
         contacttype %in% relevant_contact_types) %>%
  group_by(id) %>%
  summarise(n_contacts = n(), .groups = "drop")

n_contacts <- cohort %>%
  select(id) %>%
  left_join(n_contacts, by = "id") %>%
  mutate(n_contacts = replace_na(n_contacts, 0))


# =============================================================================
# 4. ASSEMBLE THE ANALYSIS DATASET
# =============================================================================

analysis_data <- cohort %>%
  inner_join(follow_up, by = "id") %>%
  left_join(demographics, by = "id") %>%
  left_join(bmi_final, by = "id") %>%
  left_join(prior_gastric_cancer, by = "id") %>%
  left_join(recent_gerd, by = "id") %>%
  left_join(recent_peptic_ulcer, by = "id") %>%
  left_join(n_contacts, by = "id") %>%
  mutate(
    across(c(prior_gastric_cancer, recent_gerd, recent_peptic_ulcer),
           ~ replace_na(.x, 0)),
    ppi = factor(ppi, levels = c(0, 1), labels = c("H2RA", "PPI")),
    died = as.integer(died)
  )

flow$n_complete <- sum(complete.cases(
  analysis_data %>%
    select(fup_years, died, ppi, age, gender, ethnicity, ses_person,
           calendarperiod, bmi, prior_gastric_cancer, recent_gerd,
           recent_peptic_ulcer, n_contacts)
))

cat("\n--- Analysis dataset dimensions ---\n")
cat("Rows:", nrow(analysis_data), "\n")
cat("Columns:", ncol(analysis_data), "\n")

cat("\n--- Deaths by exposure group ---\n")
print(
  analysis_data %>%
    group_by(ppi) %>%
    summarise(n = n(), deaths = sum(died), pct_died = round(100 * deaths / n, 1))
)


# =============================================================================
# 5. SAVE ANALYSIS DATASET AND COHORT-FLOW COUNTS
# =============================================================================

saveRDS(analysis_data, file = "data/processed/analysis_data.Rds")
saveRDS(flow, file = "data/processed/flow_counts.Rds")
cat("\nAnalysis dataset saved to: data/processed/analysis_data.Rds\n")
cat("Cohort-flow counts saved to: data/processed/flow_counts.Rds\n")

# =============================================================================
# END OF SCRIPT 1
# =============================================================================
