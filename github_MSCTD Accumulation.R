library(tidyverse)
library(survival)
library(lubridate)

cancer1 <- readRDS("/mnt/DATA/home/juan/system_disease/cancer1.rds")

mapping <- list(
  M00 = c("M00", "M01", "M02", "M03"),
  M05 = c("M05", "M06", "M07", "M08", "M09", "M10", "M11", "M12", "M13", "M14"),
  M15 = c("M15", "M16", "M17", "M18", "M19"),
  M20 = c("M20", "M21", "M22", "M23", "M24", "M25"),
  M30 = c("M30", "M31", "M32", "M33", "M34", "M35", "M36"),
  M40 = c("M40", "M41", "M42", "M43"),
  M45 = c("M45", "M46", "M47", "M48", "M49"),
  M50 = c("M50", "M51", "M53", "M54"),
  M60 = c("M60", "M61", "M62", "M63"),
  M65 = c("M65", "M66", "M67", "M68"),
  M70 = c("M70", "M71", "M72", "M73", "M75", "M76", "M77", "M79"),
  M80 = c("M80", "M81", "M82", "M83", "M84", "M85"),
  M86 = c("M86", "M87", "M88", "M89", "M90"),
  M91 = c("M91", "M92", "M93", "M94"),
  M95 = c("M95", "M96", "M99")
)

for (new_col in names(mapping)) {
  group_cols <- intersect(mapping[[new_col]], colnames(cancer1)) 
  if(length(group_cols) > 0){
    cancer1[[paste0(new_col, "_combine")]] <- apply(cancer1[ , group_cols], 1, function(x) { x[!is.na(x)][1] })
  }
}

msctd_cols <- c("M00_combine", "M05_combine", "M15_combine", "M20_combine", 
                "M30_combine", "M40_combine", "M45_combine", "M50_combine",
                "M60_combine", "M65_combine", "M70_combine", "M80_combine",
                "M86_combine", "M91_combine", "M95_combine")

covariates <- c("Age", "Sex", "Townsend_deprivation_index",
                "BMI", "Waist.to.hip.ratio", "Race", "Education", 
                "Smoking_status", "Alcohol_frequency", "Family_cancer", 
                "Oily_fish_intake")

columns_to_check <- covariates


cancertype <- 'C34'


cancer1[[cancertype]][is.na(cancer1[[cancertype]])] <- 0

data1 <- cancer1[complete.cases(cancer1[, columns_to_check]), ]


data1 <- data1[!data1$follow_up_time < 1, ]


data1<-data1[!(!is.na(data1$date_cancer) &
                 !is.na(data1$M100 ) &
                 data1$date_cancer - data1$M100  < 365 ),] 


msctd_dates_list <- lapply(1:nrow(data1), function(i) {
  dates <- as.Date(unlist(data1[i, msctd_cols]), origin = "1970-01-01")
  dates <- dates[!is.na(dates)]

    dates <- dates[dates >= data1$date_min[i] & dates <= data1$date_max[i]]

    dates <- sort(dates)
  return(dates)
})


data1$msctd_dates <- msctd_dates_list
data1$total_msctd <- sapply(msctd_dates_list, length)

print(table(data1$total_msctd))


data1$msctd_at_cancer <- NA

for (i in 1:nrow(data1)) {
  if (!is.na(data1$date_cancer[i]) && data1[[cancertype]][i] == 1) {
    
    msctd_dates <- data1$msctd_dates[[i]]
    cancer_date <- data1$date_cancer[i]
    
    if (length(msctd_dates) > 0) {
      msctd_before_cancer <- sum(msctd_dates < cancer_date)
      data1$msctd_at_cancer[i] <- msctd_before_cancer
    } else {
      data1$msctd_at_cancer[i] <- 0
    }
  }
}

data1$follow_up_time <- as.numeric(
  ifelse(data1[[cancertype]] == 1,
         as.numeric(data1$date_cancer - data1$date_min) / 365.25,
         as.numeric(data1$date_max - data1$date_min) / 365.25)
)

data1$max <- data1$follow_up_time
data1$cancertype <- data1[[cancertype]]

data1 <- data1[data1$max != 0, ]


names(data1)[1] <- 'id'

tm_data <- tmerge(data1 = data1, data2 = data1, id = id, 
                  tstart = 0, tstop = max)


all_msctd_events <- data1 %>%
  select(id, date_min, max, msctd_dates) %>%
  filter(sapply(msctd_dates, length) > 0) %>% 
  rowwise() %>%
  mutate(
    msctd_times = list(as.numeric(msctd_dates - date_min) / 365.25)
  ) %>%
  ungroup() %>%
  tidyr::unnest_longer(msctd_times, indices_to = "msctd_num") %>%
  filter(msctd_times > 0 & msctd_times < max) %>%
  arrange(id, msctd_times)


if (nrow(all_msctd_events) > 0) {
  tm_data <- tmerge(tm_data, all_msctd_events, id = id,
                    msctd_event = event(msctd_times))

} else {
  tm_data$msctd_event <- 0
}

tm_data <- tm_data %>%
  arrange(id, tstart) %>%
  group_by(id) %>%
  mutate(msctd_count = cumsum(msctd_event)) %>%
  ungroup()

tm_data$msctd_category <- factor(
  case_when(
    tm_data$msctd_count == 0 ~ "0",
    tm_data$msctd_count == 1 ~ "1",
    tm_data$msctd_count == 2 ~ "2",
    tm_data$msctd_count >= 3 ~ "≥3"
  ),
  levels = c("0", "1", "2", "≥3")
)

print(table(tm_data$msctd_category))


summary_stats <- tm_data %>%
  group_by(msctd_category) %>%
  summarise(
    n_intervals = n(),
    n_patients = n_distinct(id),
    n_events = sum(cancertype),
    person_years = sum(tstop - tstart),
    incidence_rate = (n_events / person_years) * 1000,
    .groups = 'drop'
  )

print(summary_stats)

cancer_patients_group <- tm_data %>%
  filter(cancertype == 1) %>%
  group_by(id) %>%
  filter(row_number() == n()) %>%  
  ungroup() %>%
  count(msctd_category, name = "n_cancer_patients")

print(cancer_patients_group)

# Cox
y <- with(tm_data, Surv(tstart, tstop, cancertype))

model1 <- coxph(y ~ msctd_category, data = tm_data)

model2 <- coxph(y ~ msctd_category + Age + Sex, data = tm_data)

FML3 <- as.formula(paste0("y ~ msctd_category + ", 
                          paste0(covariates, collapse = " + ")))
model3 <- coxph(FML3, data = tm_data)

print(summary(model3))


tm_data <- tm_data %>%
  mutate(
    msctd_numeric = case_when(
      msctd_category == "0" ~ 0,
      msctd_category == "1" ~ 1,
      msctd_category == "2" ~ 2,
      msctd_category == "≥3" ~ 3
    )
  )


trend_model1 <- coxph(y ~ msctd_numeric, data = tm_data)
trend_summary1 <- summary(trend_model1)

hr_trend1 <- exp(coef(trend_model1))
ci_trend1 <- exp(confint(trend_model1))
p_trend1 <- trend_summary1$coefficients["msctd_numeric", "Pr(>|z|)"]


trend_model2 <- coxph(y ~ msctd_numeric + Age + Sex, data = tm_data)
trend_summary2 <- summary(trend_model2)

hr_trend2 <- exp(coef(trend_model2)["msctd_numeric"])
ci_trend2 <- exp(confint(trend_model2)["msctd_numeric", ])
p_trend2 <- trend_summary2$coefficients["msctd_numeric", "Pr(>|z|)"]


FML_trend3 <- as.formula(paste0("y ~ msctd_numeric + ", 
                                paste0(covariates, collapse = " + ")))
trend_model3 <- coxph(FML_trend3, data = tm_data)
trend_summary3 <- summary(trend_model3)

hr_trend3 <- exp(coef(trend_model3)["msctd_numeric"])
ci_trend3 <- exp(confint(trend_model3)["msctd_numeric", ])
p_trend3 <- trend_summary3$coefficients["msctd_numeric", "Pr(>|z|)"]


trend_results <- data.frame(
  Model = c("Model 1 (Unadjusted)", 
            "Model 2 (Demographics)", 
            "Model 3 (Fully adjusted)"),
  HR_per_MSCTD = c(hr_trend1, hr_trend2, hr_trend3),
  Lower_CI = c(ci_trend1[1], ci_trend2[1], ci_trend3[1]),
  Upper_CI = c(ci_trend1[2], ci_trend2[2], ci_trend3[2]),
  P_trend = c(p_trend1, p_trend2, p_trend3)
)

print(trend_results, row.names = FALSE)


