library(regmedint)
library(dplyr)
library(broom)
library(mediation)
library(survival)
library(doParallel)
library(openxlsx)
library(boot)
data_clean <- readRDS("/mnt/DATA/home/juan/system_disease/data_clean.rds")

treat <- "M_status"        
outcome <- "C34"           
time <- "time"  

mediator <- c(
  "WBC_count","CRP", "Neutrophil_count", 
  "Lymphocyte_count", "Monocyte_count", "Basophil_count", "Eosinophil_count", 
  "Platelet_count", "Albumin", "RBC_count", 
  "RBC_distribution_width", "Cholesterol", "Glucose", "HbA1c", 
  "HDL", "IGF1", "LDL", "Triglycerides", "RAR", "SIRI", "SII", "NLR", 
  "PLR", "LMR", "TyG", "NHR", "NPR", "PAR", "CALLY", "MHR", "CAR", 
  "LCR", "AISI", "IBI", "PIV"
)

## Run the above mediators one by one


covariates <- c("Age", "Sex", "Townsend_deprivation_index",
                "BMI", "Waist_to_hip_ratio", "Race", "Education", "Smoking_status",
                "Alcohol_frequency", "Family_cancer", "Oily_fish_intake", 
                "Processed_meat_intake")

formula_med <- as.formula(
  paste(mediator, "~", treat, "+", paste(covariates, collapse = " + "))
)

mod_med <- lm(formula_med, data = data_clean)


names(data_clean)[1]<-'id'
data_clean$max <- data_clean$time
tm_data <- tmerge(data1 = data_clean, data2 = data_clean, id = id, tstart =0, tstop = max)

tm_data <- tmerge(tm_data, data_clean, id = id, diag_change = event(date_disease,M_status))

tm_data$M_status[tm_data$diag_change == 1]<-0

y <- with(tm_data,Surv(tstart, tstop, C34))

formula_total <- as.formula(paste0("y ~", treat, "+",paste0(covariates,collapse = "+") ))
mod_total <- coxph(formula_total, tm_data)
summary(mod_total)

formula_out <- as.formula(
  paste0("y ~ ", mediator, " + ", treat, " + ", paste(covariates, collapse = " + "))
)

mod_out <- coxph(formula_out, data = tm_data)


coef_a <- coef(mod_med)["M_status1"]
coef_b <- coef(mod_out)["WBC_count"]
coef_c_prime <- coef(mod_out)["M_status1"]
coef_c <- coef(mod_total)["M_status1"]

p_a <- summary(mod_med)$coefficients["M_status1", "Pr(>|t|)"]
p_b <- summary(mod_out)$coefficients["WBC_count", "Pr(>|z|)"]
p_c_prime <- summary(mod_out)$coefficients["M_status1", "Pr(>|z|)"]
p_c <- summary(mod_total)$coefficients["M_status1", "Pr(>|z|)"]



mediation_boot_fn <- function(data, indices) {
  tryCatch({

    d <- data[indices, ]
    d$original_id <- d$id
    d$id <- 1:nrow(d)
    
    med_fit <- lm(formula_med, data = d)
    if (any(is.na(coef(med_fit)))) return(rep(NA, 4))
    if (!paste0(treat, "1") %in% names(coef(med_fit))) return(rep(NA, 4))
    
    d$max <- d[[time]]
    
    tm_d <- tmerge(
      data1 = d, 
      data2 = d, 
      id = id, 
      tstart = 0, 
      tstop = max
    )
    
    tm_d <- tmerge(
      tm_d, 
      d, 
      id = id, 
      diag_change = event(date_disease, M_status)
    )
    
    tm_d$M_status[tm_d$diag_change == 1] <- 0

    y <- with(tm_d, Surv(tstart, tstop, get(outcome)))
    
    formula_total_boot <- as.formula(
      paste0("y ~ ", treat, " + ", paste(covariates, collapse = " + "))
    )
    
    cox_total <- coxph(formula_total_boot, data = tm_d)
    if (any(is.na(coef(cox_total)))) return(rep(NA, 4))
    if (!paste0(treat, "1") %in% names(coef(cox_total))) return(rep(NA, 4))
    

    formula_out_boot <- as.formula(
      paste0("y ~ ", mediator, " + ", treat, " + ", paste(covariates, collapse = " + "))
    )
    
    cox_out <- coxph(formula_out_boot, data = tm_d)
    if (any(is.na(coef(cox_out)))) return(rep(NA, 4))
    if (!mediator %in% names(coef(cox_out))) return(rep(NA, 4))
    if (!paste0(treat, "1") %in% names(coef(cox_out))) return(rep(NA, 4))
    

    a <- coef(med_fit)[paste0(treat, "1")]
    b <- coef(cox_out)[mediator]
    c_prime <- coef(cox_out)[paste0(treat, "1")]
    c <- coef(cox_total)[paste0(treat, "1")]
    
    if (!is.finite(a) || !is.finite(b) || !is.finite(c_prime) || !is.finite(c)) {
      return(rep(NA, 4))
    }
    
    acme <- a * b
    ade <- c_prime
    total <- c
    prop <- if (abs(total) > 1e-10) acme / total else NA
    
    return(c(acme, ade, total, prop))
    
  }, error = function(e) {
    return(rep(NA, 4))
  })
}

set.seed(111)
boot_results <- boot(data = data_clean, statistic = mediation_boot_fn, R = 1000)

acme_ci <- boot.ci(boot_results, type = "perc", index = 1)
ade_ci <- boot.ci(boot_results, type = "perc", index = 2)
total_ci <- boot.ci(boot_results, type = "perc", index = 3)
prop_ci <- boot.ci(boot_results, type = "perc", index = 4)

calc_p <- function(boot_t, est) {
  boot_t <- boot_t[!is.na(boot_t)]
  if (length(boot_t) == 0) return(NA)
  p <- ifelse(est >= 0, mean(boot_t <= 0), mean(boot_t >= 0))
  return(2 * min(p, 1 - p))
}

acme_p <- calc_p(boot_results$t[, 1], boot_results$t0[1])
ade_p <- calc_p(boot_results$t[, 2], boot_results$t0[2])
total_p <- calc_p(boot_results$t[, 3], boot_results$t0[3])
prop_p <- calc_p(boot_results$t[, 4], boot_results$t0[4])


mediation_results <- data.frame(
  Effect_Type = c(
    "ACME (control)",
    "ACME (treated)",
    "ADE (control)",
    "ADE (treated)",
    "Total Effect",
    "Prop. Mediated (control)",
    "Prop. Mediated (treated)",
    "ACME (average)",
    "ADE (average)",
    "Prop. Mediated (average)"
  ),
  Estimate = c(
    boot_results$t0[1],  # ACME (control)
    boot_results$t0[1],  # ACME (treated) 
    boot_results$t0[2],  # ADE (control)
    boot_results$t0[2],  # ADE (treated)
    boot_results$t0[3],  # Total Effect
    boot_results$t0[4],  # Prop (control)
    boot_results$t0[4],  # Prop (treated)
    boot_results$t0[1],  # ACME (average)
    boot_results$t0[2],  # ADE (average)
    boot_results$t0[4]   # Prop (average)
  ),
  CI_Lower = c(
    acme_ci$percent[4],
    acme_ci$percent[4],
    ade_ci$percent[4],
    ade_ci$percent[4],
    total_ci$percent[4],
    prop_ci$percent[4],
    prop_ci$percent[4],
    acme_ci$percent[4],
    ade_ci$percent[4],
    prop_ci$percent[4]
  ),
  CI_Upper = c(
    acme_ci$percent[5],
    acme_ci$percent[5],
    ade_ci$percent[5],
    ade_ci$percent[5],
    total_ci$percent[5],
    prop_ci$percent[5],
    prop_ci$percent[5],
    acme_ci$percent[5],
    ade_ci$percent[5],
    prop_ci$percent[5]
  ),
  P_value = c(
    acme_p,
    acme_p,
    ade_p,
    ade_p,
    total_p,
    prop_p,
    prop_p,
    acme_p,
    ade_p,
    prop_p
  )
)


mediation_results_HR <- mediation_results
mediation_results_HR$Estimate_HR <- exp(mediation_results$Estimate)
mediation_results_HR$CI_Lower_HR <- exp(mediation_results$CI_Lower)
mediation_results_HR$CI_Upper_HR <- exp(mediation_results$CI_Upper)

path_coefficients <- data.frame(
  Path = c(
    "a: M_status → WBC_count (β)",
    "b: WBC_count → Lung Cancer (log HR)",
    "c': M_status → Lung Cancer | WBC_count (log HR)",
    "c: M_status → Lung Cancer (log HR)"
  ),
  Coefficient = c(coef_a, coef_b, coef_c_prime, coef_c),
  HR = c(NA, exp(coef_b), exp(coef_c_prime), exp(coef_c)),
  P_value = c(p_a, p_b, p_c_prime, p_c),
  Significance = c(
    ifelse(p_a < 0.001, "***", ifelse(p_a < 0.01, "**", ifelse(p_a < 0.05, "*", ""))),
    ifelse(p_b < 0.001, "***", ifelse(p_b < 0.01, "**", ifelse(p_b < 0.05, "*", ""))),
    ifelse(p_c_prime < 0.001, "***", ifelse(p_c_prime < 0.01, "**", ifelse(p_c_prime < 0.05, "*", ""))),
    ifelse(p_c < 0.001, "***", ifelse(p_c < 0.01, "**", ifelse(p_c < 0.05, "*", "")))
  )
)

analysis_info <- data.frame(
  Parameter = c(
    "Sample Size",
    "Number of Events",
    "Bootstrap Simulations",
    "Confidence Method",
    "Treatment Variable",
    "Mediator Variable",
    "Outcome Variable",
    "Time Variable",
    "Covariates",
    "Model Type (Mediator)",
    "Model Type (Outcome)"
  ),
  Value = c(
    format(nrow(data_clean), big.mark = ","),
    format(sum(data_clean$C34), big.mark = ","),
    "1000",
    "Percentile Bootstrap",
    treat,
    mediator,
    outcome,
    time,
    paste(covariates, collapse = ", "),
    "Linear Regression",
    "Cox Proportional Hazards"
  )
)


save_path <- "/mnt/DATA/home/juan/system_disease/mediation.xlsx"

wb <- createWorkbook()

addWorksheet(wb, "Mediation_Results")
writeData(wb, "Mediation_Results", mediation_results, startRow = 1)

addWorksheet(wb, "Mediation_Results_HR")
writeData(wb, "Mediation_Results_HR", mediation_results_HR, startRow = 1)

addWorksheet(wb, "Path_Coefficients")
writeData(wb, "Path_Coefficients", path_coefficients, startRow = 1)

addWorksheet(wb, "Analysis_Info")
writeData(wb, "Analysis_Info", analysis_info, startRow = 1)

saveWorkbook(wb, save_path, overwrite = TRUE)