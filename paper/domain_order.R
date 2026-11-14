# paper/domain_order.R ---------------------------------------------------------
# Display order of the reporting domains (config/domains.csv) for the patient
# card and the within-outcome correlation table. Sourced from the project root.
# The order is the oxygen-delivery sequence (uptake, transport, carrying
# capacity, tissue adequacy, acid-base consequence, renal handling), then the
# remaining systems, fixed before any card was drawn; it carries no ranking of
# the measurements. Within a domain, channels keep their order in domains.csv.
# ------------------------------------------------------------------------------
DOMAIN_ORDER <- c("neuro_sedation", "respiratory_support", "hemodynamic_support", "hematologic_support",
                  "perfusion_metabolic", "acid_base", "renal_urinary", "electrolyte_osmotic",
                  "hepatobiliary", "inflammatory_response", "glycemic_regulation")
DOM <- read.csv("config/domains.csv", stringsAsFactors = FALSE)
stopifnot(setequal(DOMAIN_ORDER, unique(DOM$domain)), !anyDuplicated(DOM$signal))
DOM <- DOM[order(match(DOM$domain, DOMAIN_ORDER), seq_len(nrow(DOM))), ]
DOM_LAB <- c(neuro_sedation = "Neurologic", respiratory_support = "Respiratory",
             hemodynamic_support = "Hemodynamic", hematologic_support = "Hematologic",
             perfusion_metabolic = "Perfusion", acid_base = "Acid-base", renal_urinary = "Renal",
             electrolyte_osmotic = "Electrolyte", hepatobiliary = "Hepatobiliary",
             inflammatory_response = "Inflammatory", glycemic_regulation = "Glycemic")
stopifnot(setequal(names(DOM_LAB), DOMAIN_ORDER))
CH_LAB <- c(mbp = "Mean blood pressure", heart_rate = "Heart rate", spo2 = "Oxygen saturation",
            resp_rate = "Respiratory rate", glucose = "Glucose", gcs_motor = "Glasgow motor",
            gcs_eyes = "Glasgow eyes", gcs_verbal = "Glasgow verbal", urine_output_rate = "Urine output rate",
            creatinine = "Creatinine", platelet = "Platelets", hemoglobin = "Hemoglobin",
            temperature = "Temperature", sodium = "Sodium", bicarbonate = "Bicarbonate", bun = "Urea nitrogen",
            wbc = "White cell count", lactate = "Lactate", bilirubin_total = "Bilirubin")
stopifnot(setequal(names(CH_LAB), DOM$signal))
