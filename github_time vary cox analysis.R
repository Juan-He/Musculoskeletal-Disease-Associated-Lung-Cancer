library(tidyverse)
library(survival)
library(stringr)
library(openxlsx)

cancer1 <- readRDS("/mnt/DATA/home/juan/system_disease/cancer1.rds")

save_path <- "/mnt/DATA/home/juan/system_disease"

save_path1 <- paste(save_path ,name, 'one year.xlsx', sep = "")

#main analysis
{
  # Tmerge ####  
  performCoxAnalysis11 <- function(cancer1, diagnosis_code,cancertype,columns_check,mul_variants) {

    
    first_char <- substr(diagnosis_code, 1, 1)
    specific_code <- paste0(first_char, "100")
    
    cancer1$DIAG_TYPE <- ifelse(!is.na(cancer1[[diagnosis_code]]), 1, 0)
    
    data1 <- cancer1
    

    data1[[cancertype]][is.na(data1[[cancertype]])]<-0
  
    data1$DIAG_TYPE <- ifelse(
      !is.na(data1[[diagnosis_code]]) & !is.na(data1$date_max) & data1[[diagnosis_code]] >= data1$date_max,
      0,
      data1$DIAG_TYPE
    )
    
    
    data1 <- data1[complete.cases(data1[, columns_to_check]), ]
    
    data1<-data1[!data1$follow_up_time<1,]

    
    data1<-data1[!(!is.na(data1$date_cancer) &
                     !is.na(data1[[specific_code]] ) &
                     data1$date_cancer - data1[[specific_code]]  < 365 ),]  

    data1$follow_up_time<-as.numeric(ifelse(data1[[cancertype]]==1,as.numeric(data1$date_cancer-data1$date_min)/365.25,as.numeric(data1$date_max-data1$date_min)/365.25))
    
    data1$max<-data1$follow_up_time
    
    data1$date_disease <- as.Date(data1[[diagnosis_code]],format = "%Y-%m-%d" )

    # 癌症变量
    data1$cancertype <- data1[[cancertype]]
    
   data1$date_disease <-as.numeric( data1$date_disease - as.numeric(data1$date_min))/365.25
    

    data1 <- data1[data1$max !=0,]
    
    names(data1)[1]<-'id'
    tm_data <- tmerge(data1 = data1, data2 = data1, id = id, tstart =0, tstop = max)
    
    tm_data <- tmerge(tm_data, data1, id = id, diag_change = event(date_disease,DIAG_TYPE))
    
    tm_data$DIAG_TYPE[tm_data$diag_change == 1]<-0
    
    # colnames(tm_data)
    
    
    y <- with(tm_data,Surv(tstart, tstop, cancertype))
    
    FML <- as.formula(paste0("y ~", paste0(mul_variants1,collapse = "+") ))
    # #summary(data1$max)
    cox <- coxph( FML, tm_data)
    sum <- summary(cox)
    
        p_value <- sum$coefficients[, 5]
    HR <- sum$coef[,2]  # exp(beta)  
    HR.confint.lower <- sum$conf.int[,"lower .95"]
    HR.confint.upper <- sum$conf.int[,"upper .95"]
    HR1 <- paste0(round(HR,2), " (", round(HR.confint.lower,2), "-", round(HR.confint.upper,2), ")")
    N  <- length(data1$id)
    
    People_with_disease <-paste(sum(data1[[cancertype]] == 1 & data1$DIAG_TYPE == 1),'/',round(sum( tm_data$tstop[tm_data$DIAG_TYPE == 1]-tm_data$tstart[tm_data$DIAG_TYPE == 1]),0))
    People_without_disease <-paste(sum(data1[[cancertype]] == 1 & data1$DIAG_TYPE == 0),'/',round(sum(  tm_data$tstop[tm_data$DIAG_TYPE == 0]-tm_data$tstart[tm_data$DIAG_TYPE == 0]),0))
    
    People_with_disease1 <-paste(sum(data1[[cancertype]] == 1 & data1$DIAG_TYPE == 1),'/',sum(data1$DIAG_TYPE == 1))
    People_without_disease1 <-paste(sum(data1[[cancertype]] == 1 & data1$DIAG_TYPE == 0),'/',sum(data1$DIAG_TYPE == 0))
    
    num2 <- ifelse(sum(data1[[cancertype]] == 1 & data1$DIAG_TYPE == 1)<10 | sum(data1$DIAG_TYPE == 1) <100 ,0,1)

        Uni_cox_model <- data.frame('Characteristics' = diagnosis_code,
                                'HR' = HR ,
                                'HR1' = HR1, 
                                'HR.confint.lower'=HR.confint.lower,
                                'HR.confint.upper'=HR.confint.upper,
                                'p_value' = p_value,
                                'People_with_disease' = People_with_disease,
                                'People_without_disease' = People_without_disease,
                                'Exposome' = People_with_disease1 ,
                                'Non-exposome' =  People_without_disease1 ,
                                'N' = N, 
                                'Num2' = num2)
    
    Uni_cox_model_first_row <- Uni_cox_model[1, ]
    
    return(Uni_cox_model_first_row)
  }
  
  results_list11 <- list()

    for (i in seq_along(cancertype)) {
    variable <- cancertype[i]
    
    result_df <- data.frame('Characteristics' = character(),
                            'HR' = numeric() ,
                            'HR1' = numeric(),
                            'HR.confint.lower'=numeric(),
                            'HR.confint.upper'=numeric(),
                            'p_value' = numeric(),
                            'People_with_disease' = numeric(),
                            'People_without_disease' = numeric(),
                            'Exposome' = numeric() ,
                            'Non-exposome' = numeric() ,
                            'N'= numeric(),
                            'Num2' = numeric())
    
    for (col_name in col) {
      temp_result <- tryCatch({
        result <- performCoxAnalysis11( cancer1, col_name, variable) 
      }, error = function(e) {
        data.frame('Characteristics' = col_name,
                   'HR' = NA,
                   'HR1' = NA,
                   'HR.confint.lower' = NA,
                   'HR.confint.upper' = NA,
                   'p_value' = NA,
                   'People_with_disease' = NA,
                   'People_without_disease' = NA,
                   'Exposome' = NA,
                   'Non-exposome' = NA,
                   'N' = NA)
      })
      result_df <- rbind(result_df, temp_result)  
    }
    
    filtered_df <- result_df %>%
      filter(!grepl("100", Characteristics) & !is.na(p_value)) %>%
      mutate(Bonferroni = p.adjust(p_value, method = "bonferroni")) %>%
      mutate(FDR = p.adjust(p_value, method = "fdr")) %>%
      select(Characteristics,FDR ,Bonferroni ) 
    
    result_df <- result_df %>%
      left_join(filtered_df, by = "Characteristics")
    
    
    results_list11[[i]] <- result_df
  }
  
  
  performCoxAnalysis22 <- function(cancer1, diagnosis_code,cancertype,columns_check,mul_variants) {

    first_char <- substr(diagnosis_code, 1, 1)
    specific_code <- paste0(first_char, "100")
    #colnames(cancer1)
    cancer1$DIAG_TYPE <- ifelse(!is.na(cancer1[[diagnosis_code]]), 1, 0)
    
    data1 <- cancer1
    
    data1[[cancertype]][is.na(data1[[cancertype]])]<-0
    
    data1$DIAG_TYPE <- ifelse(
      !is.na(data1[[diagnosis_code]]) & !is.na(data1$date_max) & data1[[diagnosis_code]] >= data1$date_max,
      0,
      data1$DIAG_TYPE
    )
    
    
    data1 <- data1[complete.cases(data1[, columns_to_check]), ]
    
    data1<-data1[!data1$follow_up_time<1,]

    
    data1<-data1[!(!is.na(data1$date_cancer) &
                     !is.na(data1[[specific_code]] ) &
                     data1$date_cancer - data1[[specific_code]]  < 365 ),]  


    data1$follow_up_time<-as.numeric(ifelse(data1[[cancertype]]==1,as.numeric(data1$date_cancer-data1$date_min)/365.25,as.numeric(data1$date_max-data1$date_min)/365.25))
    
    data1$max<-data1$follow_up_time
    
    data1$date_disease <- as.Date(data1[[diagnosis_code]],format = "%Y-%m-%d" )

    data1$cancertype <- data1[[cancertype]]

    data1$date_disease <-as.numeric( data1$date_disease - as.numeric(data1$date_min))/365.25

    data1 <- data1[data1$max !=0,]
    
    names(data1)[1]<-'id'
    tm_data <- tmerge(data1 = data1, data2 = data1, id = id, tstart =0, tstop = max)
    
    tm_data <- tmerge(tm_data, data1, id = id, diag_change = event(date_disease,DIAG_TYPE))
    
    tm_data$DIAG_TYPE[tm_data$diag_change == 1]<-0
    

    y <- with(tm_data,Surv(tstart, tstop, cancertype))
    
    FML <- as.formula(paste0("y ~", paste0(mul_variants2,collapse = "+") ))
    cox <- coxph( FML, tm_data)
    sum <- summary(cox)
    
    
    p_value <- sum$coefficients[, 5]
    HR <- sum$coef[,2]  # exp(beta)  
    HR.confint.lower <- sum$conf.int[,"lower .95"]
    HR.confint.upper <- sum$conf.int[,"upper .95"]
    HR1 <- paste0(round(HR,2), " (", round(HR.confint.lower,2), "-", round(HR.confint.upper,2), ")")
    N  <- length(data1$id)
    
    People_with_disease <-paste(sum(data1[[cancertype]] == 1 & data1$DIAG_TYPE == 1),'/',round(sum( tm_data$tstop[tm_data$DIAG_TYPE == 1]-tm_data$tstart[tm_data$DIAG_TYPE == 1]),0))
    People_without_disease <-paste(sum(data1[[cancertype]] == 1 & data1$DIAG_TYPE == 0),'/',round(sum(  tm_data$tstop[tm_data$DIAG_TYPE == 0]-tm_data$tstart[tm_data$DIAG_TYPE == 0]),0))
    
    People_with_disease1 <-paste(sum(data1[[cancertype]] == 1 & data1$DIAG_TYPE == 1),'/',sum(data1$DIAG_TYPE == 1))
    People_without_disease1 <-paste(sum(data1[[cancertype]] == 1 & data1$DIAG_TYPE == 0),'/',sum(data1$DIAG_TYPE == 0))
    
    num2 <- ifelse(sum(data1[[cancertype]] == 1 & data1$DIAG_TYPE == 1)<10 | sum(data1$DIAG_TYPE == 1) <100 ,0,1)
    Uni_cox_model <- data.frame('Characteristics' = diagnosis_code,
                                'HR' = HR ,
                                'HR1' = HR1, 
                                'HR.confint.lower'=HR.confint.lower,
                                'HR.confint.upper'=HR.confint.upper,
                                'p_value' = p_value,
                                'People_with_disease' = People_with_disease,
                                'People_without_disease' = People_without_disease,
                                'Exposome' = People_with_disease1 ,
                                'Non-exposome' =  People_without_disease1 ,
                                'N' = N, 
                                'Num2' = num2)
    Uni_cox_model_first_row <- Uni_cox_model[1, ]
    
    return(Uni_cox_model_first_row)
  }
  
  results_list22 <- list()
  for (i in seq_along(cancertype)) {
    variable <- cancertype[i]
    
    result_df <- data.frame('Characteristics' = character(),
                            'HR' = numeric() ,
                            'HR1' = numeric(),
                            'HR.confint.lower'=numeric(),
                            'HR.confint.upper'=numeric(),
                            'p_value' = numeric(),
                            'People_with_disease' = numeric(),
                            'People_without_disease' = numeric(),
                            'Exposome' = numeric() ,
                            'Non-exposome' = numeric() ,
                            'N'= numeric(),
                            'Num2' = numeric())
    
    for (col_name in col) {
      temp_result <- tryCatch({
        result <- performCoxAnalysis22( cancer1, col_name, variable) 
      }, error = function(e) {
        data.frame('Characteristics' = col_name,
                   'HR' = NA,
                   'HR1' = NA,
                   'HR.confint.lower' = NA,
                   'HR.confint.upper' = NA,
                   'p_value' = NA,
                   'People_with_disease' = NA,
                   'People_without_disease' = NA,
                   'Exposome' = NA,
                   'Non-exposome' = NA,
                   'N' = NA)
      })
      result_df <- rbind(result_df, temp_result)  
    }
    
    filtered_df <- result_df %>%
      filter(!grepl("100", Characteristics) & !is.na(p_value)) %>%
      mutate(Bonferroni = p.adjust(p_value, method = "bonferroni")) %>%
      mutate(FDR = p.adjust(p_value, method = "fdr")) %>%
      select(Characteristics,FDR ,Bonferroni ) 
    
    result_df <- result_df %>%
      left_join(filtered_df, by = "Characteristics")
    
    
    results_list22[[i]] <- result_df
  }
  
  
  
  performCoxAnalysis33 <- function(cancer1, diagnosis_code,cancertype,columns_check,mul_variants) {

    first_char <- substr(diagnosis_code, 1, 1)
    specific_code <- paste0(first_char, "100")
    #colnames(cancer1)
    cancer1$DIAG_TYPE <- ifelse(!is.na(cancer1[[diagnosis_code]]), 1, 0)
    
    data1 <- cancer1
    
    data1[[cancertype]][is.na(data1[[cancertype]])]<-0
    
    data1$DIAG_TYPE <- ifelse(
      !is.na(data1[[diagnosis_code]]) & !is.na(data1$date_max) & data1[[diagnosis_code]] >= data1$date_max,
      0,
      data1$DIAG_TYPE
    )
    
    
    data1 <- data1[complete.cases(data1[, columns_to_check]), ]
    
    data1<-data1[!data1$follow_up_time<1,]

    
    data1<-data1[!(!is.na(data1$date_cancer) &
                     !is.na(data1[[specific_code]] ) &
                     data1$date_cancer - data1[[specific_code]]  < 365 ),]  

    data1$follow_up_time<-as.numeric(ifelse(data1[[cancertype]]==1,as.numeric(data1$date_cancer-data1$date_min)/365.25,as.numeric(data1$date_max-data1$date_min)/365.25))
    
    data1$max<-data1$follow_up_time
    
    data1$date_disease <- as.Date(data1[[diagnosis_code]],format = "%Y-%m-%d" )

    data1$cancertype <- data1[[cancertype]]

    data1$date_disease <-as.numeric( data1$date_disease - as.numeric(data1$date_min))/365.25
    

    data1 <- data1[data1$max !=0,]
    
    names(data1)[1]<-'id'
    tm_data <- tmerge(data1 = data1, data2 = data1, id = id, tstart =0, tstop = max)
    
    tm_data <- tmerge(tm_data, data1, id = id, diag_change = event(date_disease,DIAG_TYPE))
    
    tm_data$DIAG_TYPE[tm_data$diag_change == 1]<-0
  
    
    y <- with(tm_data,Surv(tstart, tstop, cancertype))
    
    FML <- as.formula(paste0("y ~", paste0(mul_variants3,collapse = "+") ))
    cox <- coxph( FML, tm_data)
    sum <- summary(cox)
    
    
    p_value <- sum$coefficients[, 5]
    HR <- sum$coef[,2]  # exp(beta)  
    HR.confint.lower <- sum$conf.int[,"lower .95"]
    HR.confint.upper <- sum$conf.int[,"upper .95"]
    HR1 <- paste0(round(HR,2), " (", round(HR.confint.lower,2), "-", round(HR.confint.upper,2), ")")
    N  <- length(data1$id)
    
    People_with_disease <-paste(sum(data1[[cancertype]] == 1 & data1$DIAG_TYPE == 1),'/',round(sum( tm_data$tstop[tm_data$DIAG_TYPE == 1]-tm_data$tstart[tm_data$DIAG_TYPE == 1]),0))
    People_without_disease <-paste(sum(data1[[cancertype]] == 1 & data1$DIAG_TYPE == 0),'/',round(sum(  tm_data$tstop[tm_data$DIAG_TYPE == 0]-tm_data$tstart[tm_data$DIAG_TYPE == 0]),0))
    
    People_with_disease1 <-paste(sum(data1[[cancertype]] == 1 & data1$DIAG_TYPE == 1),'/',sum(data1$DIAG_TYPE == 1))
    People_without_disease1 <-paste(sum(data1[[cancertype]] == 1 & data1$DIAG_TYPE == 0),'/',sum(data1$DIAG_TYPE == 0))
    
    num2 <- ifelse(sum(data1[[cancertype]] == 1 & data1$DIAG_TYPE == 1)<10 | sum(data1$DIAG_TYPE == 1) <100 ,0,1)
    
    Uni_cox_model <- data.frame('Characteristics' = diagnosis_code,
                                'HR' = HR ,
                                'HR1' = HR1, 
                                'HR.confint.lower'=HR.confint.lower,
                                'HR.confint.upper'=HR.confint.upper,
                                'p_value' = p_value,
                                'People_with_disease' = People_with_disease,
                                'People_without_disease' = People_without_disease,
                                'Exposome' = People_with_disease1 ,
                                'Non-exposome' =  People_without_disease1 ,
                                'N' = N, 
                                'Num2' = num2)
    Uni_cox_model_first_row <- Uni_cox_model[1, ]
    
    return(Uni_cox_model_first_row)
  }
  
  results_list33 <- list()
  for (i in seq_along(cancertype)) {
    variable <- cancertype[i]
    
    result_df <- data.frame('Characteristics' = character(),
                            'HR' = numeric() ,
                            'HR1' = numeric(),
                            'HR.confint.lower'=numeric(),
                            'HR.confint.upper'=numeric(),
                            'p_value' = numeric(),
                            'People_with_disease' = numeric(),
                            'People_without_disease' = numeric(),
                            'Exposome' = numeric() ,
                            'Non-exposome' = numeric() ,
                            'N'= numeric(),
                            'Num2' = numeric())
    
    for (col_name in col) {
      temp_result <- tryCatch({
        result <- performCoxAnalysis33( cancer1, col_name, variable) 
      }, error = function(e) {
        data.frame('Characteristics' = col_name,
                   'HR' = NA,
                   'HR1' = NA,
                   'HR.confint.lower' = NA,
                   'HR.confint.upper' = NA,
                   'p_value' = NA,
                   'People_with_disease' = NA,
                   'People_without_disease' = NA,
                   'Exposome' = NA,
                   'Non-exposome' = NA,
                   'N' = NA)
      })
      result_df <- rbind(result_df, temp_result)  
    }
    
    filtered_df <- result_df %>%
      filter(!grepl("100", Characteristics) & !is.na(p_value)) %>%
      mutate(Bonferroni = p.adjust(p_value, method = "bonferroni")) %>%
      mutate(FDR = p.adjust(p_value, method = "fdr")) %>%
      select(Characteristics,FDR ,Bonferroni ) 
    
    result_df <- result_df %>%
      left_join(filtered_df, by = "Characteristics")
    
    
    results_list33[[i]] <- result_df
  }
  
  

  wb11 <- createWorkbook()
  
  sheet_names <- c('C34',"Non_small", "Small_cell", "Adeno", "Squamous") 
  
  add_results_to_workbook <- function(result_list, model_name, sheet_names, wb) 
    for (i in seq_along(result_list)) {
      sheet_name <- paste(model_name, sheet_names[i], sep = "_") 
      addWorksheet(wb, sheet_name) 
      writeData(wb, sheet = sheet_name, result_list[[i]],rowNames = T) 
    }
  }
  

  add_results_to_workbook(results_list11, "Model1", sheet_names, wb11)
  add_results_to_workbook(results_list22, "Model2", sheet_names, wb11)
  add_results_to_workbook(results_list33, "Model3", sheet_names, wb11)
  
  
  saveWorkbook(wb11, save_path1, overwrite = TRUE)
  
  
}
