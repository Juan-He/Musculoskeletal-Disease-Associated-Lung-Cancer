library(forestploter)
library(tidyverse)
library(survival)
library(readxl)
library(readr)
library(openxlsx)
library(data.table)
library(jstable)
library(ggplot2)

save_path <- "/mnt/DATA/home/juan/system_disease"


tm <- forest_theme(base_size = 10,  
                   
                   # Confidence interval point shape, line type/color/width
                   ci_pch = 15,   
                   ci_col = "#762a83",    
                   ci_fill = "blue",     
                   ci_alpha = 0.8,        
                   ci_lty = 1,            
                   ci_lwd = 1.5,          
                   ci_Theight = 0.2, 
                   # Reference line width/type/color   
                   refline_lwd = 1,       
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   # Vertical line width/type/color  
                   vertline_lwd = 1,              
                   vertline_lty = "dashed",
                   vertline_col = "grey20",
                   # Change summary color for filling and borders  
                   summary_fill = "yellow",       
                   summary_col = "#4575b4",
                   # Footnote font size/face/color  
                   footnote_cex = 0.6,
                   footnote_fontface = "italic",
                   footnote_col = "red",
                   text = element_text(family = "myFont")
)

cancertype
performCoxAnalysis4 <- function(cancer1, diagnosis_code,cancertype) {

  
  
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
  

  data1<-data1[!(!is.na(data1$date_cancer) &
                   !is.na(data1[[specific_code]] ) &
                   data1$date_cancer - data1[[specific_code]]  < 365 ),]  

  data1$follow_up_time<-as.numeric(ifelse(data1[[cancertype]]==1,as.numeric(data1$date_cancer-data1$date_min)/365.25,as.numeric(data1$date_max-data1$date_min)/365.25))
  
  data1$max<-data1$follow_up_time
  
  data1$date_disease <- as.Date(data1[[diagnosis_code]],format = "%Y-%m-%d" )
  data1<-data1[!data1$follow_up_time<1,]
  
  data1$cancertype <- data1[[cancertype]]

  data1$date_disease <-as.numeric( data1$date_disease - as.numeric(data1$date_min))/365.25
  
  data1 <- data1[data1$max !=0,]
  
  names(data1)[1]<-'id'
  tm_data <- tmerge(data1 = data1, data2 = data1, id = id, tstart =0, tstop = max)
  
  tm_data <- tmerge(tm_data, data1, id = id, diag_change = event(date_disease,DIAG_TYPE))
  
  tm_data$DIAG_TYPE[tm_data$diag_change == 1]<-0
  
  
  tm_data$Smoking <- recode_factor(tm_data$Smoking_status, 
                                   'Never' = 'Never',
                                   'Previous' = 'Smoker',
                                   'Current' = 'Smoker',
                                   'Prefer not to answer' = 'Prefer not to answer'
  )
  
  
  tm_data$Race1 <- recode_factor(tm_data$Race, 
                                 'White' = 'White',
                                 'Mixed' = 'Non-White',
                                 'Asian or Asian British' = 'Non-White',
                                 'Black or Black British' = 'Non-White',
                                 'Other ethnic group' = 'Non-White',
                                 'Unknown' = 'Unknown')
  
  tm_data$Alcohol <- recode_factor(tm_data$Alcohol_frequency, 
                                   'Never' = 'Never',
                                   'Special occasion' = 'Drinker',
                                   '1 to 3 times a month' = 'Drinker',
                                   '1 to 2 times a week' = 'Drinker',
                                   '3 to 4 times a week' = 'Drinker',
                                   'Almost everday' = 'Drinker',
                                   'Unknown' = 'Unknown'
  )
  
  
  tm_data$Oily_fish <- recode_factor(tm_data$Oily_fish_intake,
                                     'Never' = 'Never',
                                     '<1 time a week' = 'Oily_fish_intake',
                                     '1 time a week' = 'Oily_fish_intake',
                                     '2 to 4 times a week' = 'Oily_fish_intake',
                                     '5 to 6 times a week' = 'Oily_fish_intake',
                                     'Almost everday' = 'Oily_fish_intake',
                                     'Unknown' = 'Unknown'
  )
  tm_data$Processed_meat <- recode_factor(tm_data$Processed_meat_intake,
                                          'Never' = 'Never',
                                          '<1 time a week' = 'Processed_meat_intake',
                                          '1 time a week' = 'Processed_meat_intake',
                                          '2 to 4 times a week' = 'Processed_meat_intake',
                                          '5 to 6 times a week' = 'Processed_meat_intake',
                                          'Almost everday' = 'Processed_meat_intake',
                                          'Unknown' = 'Unknown'
  )
  
  res<-TableSubgroupMultiCox(formula = Surv(tstart, tstop, cancertype) ~ DIAG_TYPE ,
                             var_subgroups = c( 'Age_group',"Sex","Race1", 'TDI_group', 'BMI_group' ,
                                                'WHR_group', "Education" ,'Smoking',"Family_cancer",
                                                "Alcohol","Oily_fish","Processed_meat"
                             ),
                             var_cov = c( "Age", "Sex", "Townsend_deprivation_index",
                                          "BMI", "Waist.to.hip.ratio", "Race", "Education", "Smoking_status",
                                          "Alcohol_frequency", "Family_cancer", "Oily_fish_intake", 
                                          "Processed_meat_intake" ),
                             data = tm_data)
  
  plot_df <- res
  plot_df[,c(2:8)][is.na(plot_df[,c(2:8)])] <- " "
  plot_df$` ` <- paste(rep(" ", nrow(plot_df)), collapse = " ")
  plot_df[,4:6] <- apply(plot_df[,4:6],2,as.numeric)
  plot_df$HR <- paste0(plot_df$`Point Estimate`, " (",plot_df$Lower, "-", plot_df$Upper, ")")
  plot_df[, 10][plot_df[, 10] == "NA (NA-NA)"] <- ""
  row1<-which(plot_df$`P value` < 0.05) 
  row2<-which(plot_df$`P for interaction` < 0.05)
  plot_df$Variable <- gsub("_", " ", plot_df$Variable)
  
  
  names(plot_df)[3]<-"Percent (%)"
  
  
  p <- forest(
    data = plot_df[,c(1,3,10,9,7,8)],
    lower = plot_df$Lower,
    upper = plot_df$Upper,
    est = plot_df$`Point Estimate`,
    ci_column = 4,
    ref_line = 1, 
    arrow_lab = c("protective factor","risk factor"),
    xlim = c(0, 3),
    ticks_at = c(0.5, 1, 1.5, 2, 2.5),
    theme = tm
  )
  
  p <- add_border(p, part = "header")
  p <- edit_plot(p, col = 5 ,
                 row = row1,
                 gp = gpar(fontface = "bold"))
  p <- edit_plot(p, col = 6 ,
                 row = row2,
                 gp = gpar(fontface = "bold"))
  
  
  
  return(p)
}

for (col_name in col) {

  p <- performCoxAnalysis4(cancer1, col_name,cancertype)
  file_name <- paste0(save_path,col_name, ".pdf")
  ggsave(file_name, p, width = 12, height = 12, dpi = 1000)
}

