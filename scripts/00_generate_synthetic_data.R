# =============================================================================
# EHR COHORT ANALYSIS PIPELINE — DEMONSTRATION
# Script 0: Generate synthetic data
#
# Produces a fake patient population with prescription, clinical, and
# demographic records shaped like a typical UK primary care database export.
# Every value here is randomly generated; none of it represents real
# patients. This script exists purely so scripts 01-03 can be run end to
# end without any real dataset attached to the repository.
#
# Outputs (all under data/):
#   raw/patient.csv, raw/prescription.csv, raw/clinic.csv, raw/other.csv,
#   raw/death.csv, raw/ses.csv, raw/contactdates.csv
#   codelists/ppi.csv, codelists/h2ra.csv, codelists/gastric_cancer.csv,
#   codelists/gerd.csv, codelists/peptic_ulcer.csv
# =============================================================================

library(tidyverse)

set.seed(42)

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/codelists", recursive = TRUE, showWarnings = FALSE)

n_patients <- 5000

# Period covered by the synthetic database. Change these to widen or
# narrow the simulated data if needed.
db_start <- as.Date("1995-01-01")
db_end   <- as.Date("2020-01-01")

# --- Codelists -----------------------------------------------------------
# Toy numeric codes standing in for real clinical terminology codes.
write_csv(
  tibble(prescode = 90001:90003,
         term = c("Omeprazole", "Lansoprazole", "Esomeprazole")),
  "data/codelists/ppi.csv"
)
write_csv(
  tibble(prescode = 91001:91002,
         term = c("Ranitidine", "Cimetidine")),
  "data/codelists/h2ra.csv"
)
write_csv(
  tibble(clincode = 92001:92002,
         term = c("Gastric cancer", "Gastric cancer, recurrence")),
  "data/codelists/gastric_cancer.csv"
)
write_csv(
  tibble(clincode = 93001:93003,
         term = c("Gastro-oesophageal reflux disease", "Reflux oesophagitis", "Heartburn")),
  "data/codelists/gerd.csv"
)
write_csv(
  tibble(clincode = 94001:94002,
         term = c("Peptic ulcer", "Duodenal ulcer")),
  "data/codelists/peptic_ulcer.csv"
)


# --- Patient file ----------------------------------------------------------

patient <- tibble(
  id     = 1:n_patients,
  yob    = sample(1925:1995, n_patients, replace = TRUE),
  gender = sample(c("Male", "Female"), n_patients, replace = TRUE),
  eth5   = sample(c("White", "South-Asian", "Black", "Mixed", "Other"),
                  n_patients, replace = TRUE,
                  prob = c(0.75, 0.10, 0.08, 0.04, 0.03)),
  start  = as.Date(runif(n_patients, as.numeric(db_start), as.numeric(db_end) - 365),
                    origin = "1970-01-01"),
  accept = 1
)

leaves <- rbinom(n_patients, 1, 0.3) == 1
patient$leavepractice <- as.Date(NA)
patient$leavepractice[leaves] <- pmin(
  db_end,
  patient$start[leaves] + round(runif(sum(leaves), 365, 365 * 15))
)

write_csv(patient, "data/raw/patient.csv")


# --- Death file --------------------------------------------------------------

died <- rbinom(n_patients, 1, 0.25) == 1
death <- tibble(
  id        = 1:n_patients,
  death     = as.integer(died),
  deathdate = as.Date(NA)
)
death$deathdate[died] <- pmin(
  db_end,
  patient$start[died] + round(runif(sum(died), 365, 365 * 18))
)
write_csv(death, "data/raw/death.csv")


# --- Socioeconomic status file -----------------------------------------------

ses <- tibble(
  id         = 1:n_patients,
  ses_person = sample(c("Least Deprived (1)", "2", "3", "4", "Most Deprived (5)"),
                       n_patients, replace = TRUE)
)
write_csv(ses, "data/raw/ses.csv")


# --- Prescription file ---------------------------------------------------------
# Every patient gets exactly one PPI or H2RA prescription (roughly 4:1 in
# favour of PPI, reflecting real-world prescribing patterns), so the whole
# synthetic population forms the exposed cohort.

is_ppi <- rbinom(n_patients, 1, 0.8) == 1
presc_date <- pmin(
  db_end,
  patient$start + round(runif(n_patients, 30, 365 * 10))
)

prescription <- tibble(
  id        = 1:n_patients,
  eventdate = as.Date(presc_date, origin = "1970-01-01"),
  prescode  = ifelse(is_ppi,
                      sample(90001:90003, n_patients, replace = TRUE),
                      sample(91001:91002, n_patients, replace = TRUE))
)
write_csv(prescription, "data/raw/prescription.csv")


# --- Clinic and Other files: BMI/height/weight, and condition codes --------
# Entity 17 = weight/BMI record, entity 18 = height record, matching the
# convention the analysis scripts expect.

next_otherid <- 1
clinic_rows  <- vector("list", n_patients * 2)
other_rows   <- vector("list", n_patients * 2)
row_i <- 0

add_measurement <- function(id, entity, event_date, data1 = NA_real_) {
  row_i <<- row_i + 1
  oid <- next_otherid
  next_otherid <<- next_otherid + 1
  clinic_rows[[row_i]] <<- tibble(id = id, entity = entity, otherid = oid,
                                   eventdate = event_date, clincode = NA_integer_)
  other_rows[[row_i]]  <<- tibble(id = id, entity = entity, otherid = oid,
                                   data1 = data1, data3 = NA_real_)
}

for (i in 1:n_patients) {
  meas_date <- max(presc_date[i] - round(runif(1, 30, 365 * 3)), patient$start[i] + 30)
  add_measurement(i, 18, as.Date(meas_date, origin = "1970-01-01"),
                   data1 = round(rnorm(1, 1.68, 0.09), 2))
  add_measurement(i, 17, as.Date(meas_date, origin = "1970-01-01"),
                   data1 = round(rnorm(1, 78, 15), 1))
}

clinic_measurements <- bind_rows(clinic_rows[1:row_i])
other_measurements  <- bind_rows(other_rows[1:row_i])

# Condition codes: gastric cancer (rare, any time before index), GERD and
# peptic ulcer (more common, usually shortly before index).
condition_rows <- vector("list", n_patients)
row_j <- 0
add_condition <- function(id, clincode, event_date) {
  row_j <<- row_j + 1
  condition_rows[[row_j]] <<- tibble(id = id, entity = 1, otherid = NA_integer_,
                                      eventdate = event_date, clincode = clincode)
}

for (i in 1:n_patients) {
  if (runif(1) < 0.02) {
    add_condition(i, sample(92001:92002, 1),
                  as.Date(presc_date[i] - round(runif(1, 100, 3000)), origin = "1970-01-01"))
  }
  if (runif(1) < 0.15) {
    add_condition(i, sample(93001:93003, 1),
                  as.Date(presc_date[i] - round(runif(1, 10, 170)), origin = "1970-01-01"))
  }
  if (runif(1) < 0.08) {
    add_condition(i, sample(94001:94002, 1),
                  as.Date(presc_date[i] - round(runif(1, 10, 170)), origin = "1970-01-01"))
  }
}

conditions <- bind_rows(condition_rows[seq_len(row_j)])

clinic <- bind_rows(clinic_measurements, conditions)
write_csv(clinic, "data/raw/clinic.csv")
write_csv(other_measurements, "data/raw/other.csv")


# --- Contact dates file ----------------------------------------------------

contact_types <- c("Surgery consultation", "Followup/routine visit", "Clinic",
                    "Telephone call from a patient", "Acute visit",
                    "Home Visit", "Emergency Consultation", "Letter received")

contactdates <- map_dfr(1:n_patients, function(i) {
  n_contacts <- rpois(1, 4)
  if (n_contacts == 0) return(NULL)
  tibble(
    id          = i,
    eventdate   = as.Date(presc_date[i] - round(runif(n_contacts, 1, 365)),
                           origin = "1970-01-01"),
    contacttype = sample(contact_types, n_contacts, replace = TRUE)
  )
})
write_csv(contactdates, "data/raw/contactdates.csv")

cat("Synthetic data written to data/raw/ and data/codelists/\n")
cat("Patients:", n_patients, "\n")

# =============================================================================
# END OF SCRIPT 0
# =============================================================================
