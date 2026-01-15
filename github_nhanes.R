library("tidyverse")
library(survey)
library(mice)
library(gtsummary)
library(flextable)

load('data.rda')
colnames(cancer)
cancer<-select(cancer,c("SEQN","MCQ230A","MCQ230B","MCQ230C","MCQ230D" ))
cancer <- cancer %>%
  mutate(C34 = case_when(
    MCQ230A == 23  ~ 1,
    T ~ 0
  ))
cancer <- cancer %>%
  mutate(C34pure = case_when(
    MCQ230A == 23 & is.na(MCQ230B) & is.na(MCQ230C) & is.na(MCQ230D) ~ 1,
    TRUE ~ 0  # 
  )) %>% select(c("SEQN","C34","C34pure"))
table(cancer$C34)
table(cancer$C34pure)

colnames(demo)

demo <- demo %>%
  mutate(weight = case_when(
    SDDSRVYR %in% c(8, 9) ~ 2 / 7.2 * WTMEC2YR,
    SDDSRVYR == 66 ~ 3.2 / 7.2 * WTMECPRP,
    TRUE ~ NA_real_  
  ))

sum(is.na(demo$SDMVPSU))
sum(is.na(demo$SDMVSTRA))

demo<-select(demo,c("SEQN","SDMVPSU","SDMVSTRA" ,"weight","RIDAGEYR","RIAGENDR","RIDRETH3","DMDEDUC2","INDFMPIR"))%>%filter(!is.na(weight))
demo$RIDRETH3 <- recode_factor(demo$RIDRETH3, 
                               `1` = 'Mexican American',
                               `2` = 'Other Hispanic',
                               `3` = 'Non-Hispanic White',
                               `4` = 'Non-Hispanic Black',
                               `6` = 'Other/multiracial',
                               `7` = 'Other/multiracial'
)
demo$DMDEDUC2 <- recode_factor(demo$DMDEDUC2, 
                               `1` = 'Less Than 9th Grade',
                               `2` = '9-11th Grade',
                               `3`= 'High School Grad/GED',
                               `4`= 'Some College or AA degree',
                               `5`= 'College Graduate or above'
)
demo$RIAGENDR <- ifelse(demo$RIAGENDR == 1, 'male', 'female')


colnames(bmx)
bmx <- bmx %>%
  select(c("SEQN","BMXBMI"))

colnames(rxq)
row_count <- rxq %>% filter(RXDUSE == 1 | RXDUSE == 2) %>% nrow()
rxq <- rxq %>%
  dplyr::filter(RXDUSE == 1&grepl("^M", RXDRSC1)|grepl("^M", RXDRSC2)|grepl("^M", RXDRSC3)|RXDUSE == 2 )%>%
  select(c("SEQN","RXDUSE"))

table(rxq$RXDUSE)

colnames(smq)

smq <- smq%>% mutate(smoking = case_when( SMQ020 == 1~'YES' ,SMQ040 == 1|SMQ040 == 2|SMQ040 == 3~'YES',
                                          rowSums(is.na(select(., -SEQN))) == (ncol(.) - 1)~"Unknown",
                                          T~ 'NO')) %>% select(c("SEQN","smoking"))
table(smq$smoking)
smq$smoking[is.na(smq$smoking)] <- 'Unknown'


# week-month: *4; year-month:/12
ori.alq.unit <- alq$ALQ120U
table(ori.alq.unit)
trans.unit.month <- ifelse(ori.alq.unit == 1, 4, 
                           ifelse(ori.alq.unit == 3, 1/12, 
                                  ifelse(ori.alq.unit == 7 | ori.alq.unit == 9, NA, 1)))
alq$trans.unit.month <- trans.unit.month
# View(alq)

ori.alq.quantity <- alq$ALQ120Q
trans.quantity.month <- ifelse(ori.alq.quantity >= 0,  
                               ori.alq.quantity * trans.unit.month, NA)
alq$trans.quantity.month <- trans.quantity.month

alq101 <- alq$ALQ101
trans.quantity.month.factor <- ifelse(trans.quantity.month >= 4, '≥ 4 drinks/day',
                                      ifelse(trans.quantity.month >= 1, '1-3 drinks/day',
                                             ifelse(trans.quantity.month < 1, 'Nondrinkers', 'wait')))

alq$trans.quantity.month.factor <- trans.quantity.month.factor

index.1 <- which((trans.quantity.month.factor == 'wait' | is.na(trans.quantity.month.factor)) & alq101 == 1)
trans.quantity.month.factor[index.1] <- '1-3 drinks/day'

index.nondrinker <- which(alq101 == 2)
trans.quantity.month.factor[index.nondrinker] <- 'Nondrinkers'

table(trans.quantity.month.factor)
alq$alq.group <- trans.quantity.month.factor

# alq$alq.group[rowSums(is.na(alq[, -1])) == (ncol(alq) - 1)] <- 'Unknown'
alq$alq.group[is.na(alq$alq.group)] <- 'Unknown'
table(alq$alq.group)
alq <- select(alq,c("SEQN","alq.group" ))

diet_1 <- diet_1 %>% mutate(processed_meat = case_when(
  DR1CCMTX == 7 & DR1CCMTX == 13 ~ "Yes",
  is.na(DR1CCMTX) ~ "Unknown",
  TRUE ~ "No"
))%>% select(c("SEQN","processed_meat"))
diet_total <- diet_total %>% mutate(fish = case_when(
  DRD370V == 1 ~ "No",
  DRD370V == 2 ~ "Yes",
  is.na(DRD370V)& is.na(DRD360) ~ "Unknown",
  TRUE ~ "Yes"
))%>% select(c("SEQN","fish"))

paper.data <- left_join(cancer, demo, by = "SEQN") %>%
  left_join(., rxq, by = "SEQN") %>%
  left_join(., smq, by = "SEQN") %>%
  left_join(., alq, by = "SEQN") %>%
  left_join(., diet_1, by = "SEQN") %>%
  left_join(., diet_total, by = "SEQN") %>%
  left_join(., bmx, by = "SEQN")

paper.data <- paper.data %>%
  dplyr::filter(RXDUSE == 1 | RXDUSE == 2)
paper.data$RXDUSE <- ifelse(paper.data$RXDUSE ==1,"Yes",'No')
table(paper.data$RXDUSE)
colnames(paper.data)

paper.data <- paper.data %>%rename(BMI = BMXBMI, alcohol = alq.group, age = RIDAGEYR,sex = RIAGENDR,
                                   race = RIDRETH3, education = DMDEDUC2, PIR = INDFMPIR, lung_cancer = C34)

table(paper.data$lung_cancer)

set.seed(123)

paper.data  <- paper.data %>%
  group_by(SEQN) %>%
  slice_sample(n = 1) %>%
  ungroup()  

colnames(paper.data)
missing_values <- sapply(paper.data, function(x) sum(is.na(x)))

print(missing_values)
data_imputed <- mice(paper.data, m = 5, method ="rf", ,seed = 123)

data_imputed <- complete(data_imputed, 3)

missing_values <- sapply(data_imputed, function(x) sum(is.na(x)))
print(missing_values)


data_imputed <- data_imputed %>%
  replace_na(list(
    smoking = "Unknown",  
    alcohol = "Unknown",  
    processed_meat = "Unknown",  
    fish = "Unknown"  
  ))

save(data_imputed, file = "imputed_data.rda")

colnames(paper.data)
paper.data <- paper.data %>%drop_na()


load('imputed_data.rda')
paper.data <- paper.data %>%rename(diseases = RXDUSE)
imputed_data <- data_imputed %>%rename(diseases = RXDUSE)

design <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~weight, data = paper.data,nest = TRUE)
design_imputed <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~weight, data =imputed_data ,nest = TRUE)
colnames(paper.data)

#### table1 ####

b <- tbl_svysummary(design_imputed , missing = 'no', by=diseases,
                    include = c( lung_cancer, age,sex,BMI,race,education,PIR,smoking,alcohol,processed_meat,fish) ,
                    statistic = list(all_continuous()  ~ "{mean} ({sd})", 
                                     all_categorical() ~ "{n_unweighted} ({p}%)")
                    ,digits = list(all_continuous()~2, all_categorical()~ 2),) %>%  
  add_overall( all_stat_cols() ~ "**{level}**,N = {n_unweighted}({style_percent(p)}%)" ) %>%  
  modify_header(all_stat_cols() ~ "**{level}**, N = {n_unweighted} ({style_percent(p)}%)")%>%  
  modify_spanning_header(
    stat_0 ~ NA,
    update = all_stat_cols() ~ "**MSCTDs**",
  )

b <- as_flex_table(b)	

save_as_docx(b, path = "table1.docx")

logist_model <- svyglm(lung_cancer ~ diseases , design = design_imputed , family = quasibinomial )
b1 <- tbl_regression(logist_model, exponentiate = TRUE,digits =2 )

logist_mode2 <- svyglm(lung_cancer ~ diseases+age+sex , design = design_imputed , family = quasibinomial )
b2<- tbl_regression(logist_mode2, exponentiate = TRUE, include = c(diseases),digits =2 )

logist_mode3 <- svyglm(lung_cancer ~ diseases+age+race+sex+BMI+education+PIR+smoking+alcohol+processed_meat+fish, design = design_imputed, family = quasibinomial )
b3 <- tbl_regression(logist_mode3, exponentiate = TRUE, include = c(diseases),digits =2 )

bb <- tbl_merge(
  tbls = list(b1,b2,b3),
  tab_spanner = c("**Model1**", "**Model2**", "**Model3**"))
bb <- as_flex_table(bb)	

save_as_docx(bb, path = "logist.docx")

