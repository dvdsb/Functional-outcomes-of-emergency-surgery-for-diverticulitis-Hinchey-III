
gc()
options("install.lock"=FALSE)
options(scipen = 6, digits = 4)
`%notin%` <- Negate(`%in%`)
library(Hmisc)
library(tidyverse)
library(haven)
library(dplyr)
library(MASS)
library(effects)
library(DescTools)
library(ggplot2)
library(glm2)
library(AER)
library(gtsummary)
library(knitr)
library(kableExtra)
library(rms)
library(splines)
library(ordinal)
library(lmtest)
library(janitor)
library(ggpubr)
library(scimple)
library(grid)
library(tiff)
library(naniar)

# Sample size                           
# Assume the resection group provide a distribution across the five categories as
p. <- c(0.15, 0.2, 0.3, 0.2, 0.15)
odds.ratio. <- c(0.75, 0.675, 0.5)
alpha. <- 0.05
# Lavage
n1. <- 95   # 11     95
# Resection
n2. <- 119   # 80 119

# Proportional odds model
Sample.size <- popower(p=p., odds.ratio=odds.ratio., n1=n1., n2=n2.) 

#-----------------------------------------------------------------------------------
# Read data
PRO <- read_rds("T://SSORG/LapLAV/Data/Frågeformulär/laplavfrage.rds")
Clin <- read_rds("S://RCVR/kod/LapLavQuestionnaire/LapLavA4_210908.rds")
AND_0 <- read_rds("T://SSORG/LapLAV/Data/Frågeformulär/as210914.rds")
P <- PRO %>% 
  dplyr::select(studienummer, Q1, Q3, Q79) %>%
  mutate(SET=1)
C <- Clin %>%
  dplyr::select(studienummer, Surgery, treatment, age_cent, 
                kon, Sjukhustyp, hjart_karl_sjd, diabetes, 
                iva_under_VT, VT_Index, VT_90_days, 
                Comorbidity, kol, CVD, immunosuppr, sepsis_op, age) %>%
  mutate(kon=factor(kon), 
         Sjukhustyp=factor(Sjukhustyp), 
         hjart_karl_sjd= factor(hjart_karl_sjd), 
         diabetes=factor(diabetes), 
         iva_under_VT=factor(iva_under_VT), 
         kol= factor(kol, levels=c(0, 1), labels=c("No", "Yes")), 
         immunosuppr= factor(immunosuppr, levels=c(0, 1), labels=c("No", "Yes")),
         CVD= factor(CVD, levels=c(0, 1), labels=c("No", "Yes")),
         sepsis_op= factor(sepsis_op, levels=c(0, 1), labels=c("No", "Yes")))        
PC <- full_join(P, C) %>% 
  filter(SET==1 & 
           Surgery %in% c("Hartmann´s procedure",  "Laparoscopic lavage")) %>%
  mutate( Surgery= factor(case_when(Surgery=="Hartmann´s procedure"~0, 
                                    Surgery=="Laparoscopic lavage"~1), 
                          levels=c(0, 1), labels=c("Resection", "Lavage")),
          Q79 = factor(Q79, levels = c(1, 2), labels = c("Stoma: No", "Stoma: Yes")), 
          ANA_SET=1)
# Calculating and evaluating the IPW weights                                     
#Use the same weights as was used in the article on short term complications and morbidity (Samuelsson, 2021).
#See the paper on the evaluation of it´s ability to create banalced groups.

anadata.w <- tibble(read_rds("S://RCVR/kod/LapLavQuestionnaire/LapLav_PS_w_ATE.rds")) %>%
  dplyr::select(studienummer, Surgery, w) %>%
  filter(Surgery %in% c("Hartmann´s procedure",  "Laparoscopic lavage")) %>%
  mutate( S_= factor(case_when(Surgery=="Hartmann´s procedure"~0, 
                               Surgery=="Laparoscopic lavage"~1), 
                     levels=c(0, 1), labels=c("Resection", "Lavage"))) %>%
  dplyr::select(-c(Surgery))

# Merge the propensity score with the data set
PC <- full_join(PC, anadata.w) %>% 
  filter(ANA_SET==1)
#.........................................................................
desc.0 <- PC %>%
  group_by(Surgery, Q79) %>%
  summarise(median = median(age), 
            min = min(age), 
            max = max(age),
            n=n()) %>%
  kbl() %>%
  kable_classic(full_width = F, 
                html_font = "Cambria", position = "left")

desc.1 <- PC %>% 
  dplyr::select(studienummer, Surgery, Q79,
                kon, Sjukhustyp, hjart_karl_sjd, diabetes, 
                iva_under_VT, Comorbidity, kol, CVD, immunosuppr, sepsis_op) %>%
  pivot_longer(!c(studienummer, Surgery, Q79), names_to = "variable", values_to = "Response") %>%
  mutate( Response= na_if(Response , c(999))) %>%
  group_by(variable, Surgery, Q79, Response) %>%
  summarise(n = n() ) %>%
  mutate(Percent=signif(100*(n/sum(n)), digits = 2) ,
         Numbers_perc = paste0(n, "(", signif(Percent, 2), "%)"), 
         Surg.stoma=paste(str_sub(Surgery, 1, 1), str_sub(Q79, 7))) %>%
  ungroup() %>%
  dplyr::select(variable, Surg.stoma, Numbers_perc, Response) %>%
  pivot_wider(names_from=Surg.stoma, values_from=Numbers_perc) %>%
  relocate(variable, Response, "R  No", "L  No", "R  Yes", "L  Yes" ) %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Cambria", position = "left")

# Demography and patient charactertistics                                                          

P_dem <- PRO %>% 
  dplyr::select(studienummer, Q1, Q3, Q79) %>%
  mutate(SET=1)

AND_dem <- AND_0 %>% 
  dplyr::select(studienummer, smoke_status , auditc, BMI, 
                Q16,  Q28, Q30, Q31, Q32, Q33, Q34, Q35, 
                Q116 , Q118 , civil, syssel )


AND_dem2 <- AND_dem %>%
  mutate(civil = factor(case_when(civil %in% c("Divorced", "Single household without partner", 
                                               "Unmarried", "Widow/widower")~0,
                                  civil %in% c("Domestic partner", "Single household with partner", 
                                               "Married")~1), 
                        levels=c(0, 1), labels=c("Not in a relationshop", "In a relationship")), 
         syssel = factor(case_when(syssel %in% c("Early retired/disability retired", "Long-term sick leave (> 3 months", 
                                                 "Labor market measure", "Retired - Age", "Retired - Contract", "Unemployed" )~0, 
                                   syssel %in% c("Employed")~1), 
                         levels=c(0, 1), labels=c("Retired/unemplyed", "Employed")),
         Q118 = factor(case_when(Q118 %in% c("No formal education", "Elementary school", "High school, vocational school", 
                                             "Advanced higher vocational school")~0,
                                 Q118 %in% c("University/collage", "Postgraduate education")~1), 
                       levels=c(0, 1), labels=c("No university education", "University education")), 
         
         smoke_status2 = factor(case_when(smoke_status =="non_smoker"~1, smoke_status =="former_smoker"~2,
                                          smoke_status %in% c("someday_smoker", "current_smoker", "pipe_smoker")~3),
                                levels=c(1, 2, 3), labels=c("Non-smoker", "Past smoker", "Current smoker"))) %>%
  replace_with_na(replace=list(Q32=c("Don´t know", "No Applicable"), 
                               Q33=c("Don´t know", "No Applicable"), 
                               Q34=c("Don´t know", "No Applicable")))

C_dem <- Clin %>%
  dplyr::select(studienummer, Surgery, treatment, age_cent, 
                age, CCI, kon, Sjukhustyp, hjart_karl_sjd, diabetes, 
                iva_under_VT, VT_Index, VT_90_days, 
                Comorbidity, kol, CVD, immunosuppr, sepsis_op, age) %>%
  mutate(kon=factor(kon), 
         Sjukhustyp=factor(Sjukhustyp), 
         hjart_karl_sjd= factor(hjart_karl_sjd), 
         diabetes=factor(diabetes), 
         iva_under_VT=factor(iva_under_VT), 
         kol= factor(kol, levels=c(0, 1), labels=c("No", "Yes")), 
         immunosuppr= factor(immunosuppr, levels=c(0, 1), labels=c("No", "Yes")),
         CVD= factor(CVD, levels=c(0, 1), labels=c("No", "Yes")),
         sepsis_op= factor(sepsis_op, levels=c(0, 1), labels=c("No", "Yes")))        
PC_dem <- full_join(P_dem, C_dem) %>%
  full_join(AND_dem2) %>% 
  filter( Surgery %in% c("Hartmann´s procedure",  "Laparoscopic lavage")) %>%
  mutate(  SET=replace_na(SET, 99), 
           Surgery_Dem=case_when(SET==1 & Surgery=="Hartmann´s procedure"~0, 
                                 SET==1 & Surgery=="Laparoscopic lavage"~1,
                                 SET==99 & Surgery=="Hartmann´s procedure"~2, 
                                 SET==99 & Surgery=="Laparoscopic lavage"~3), 
           Surgery_Dem= factor(Surgery_Dem, levels=c(0, 1, 2, 3), 
                               labels=c("Resection", "Lavage", "N.R Resection", "N.R Lavage")),
           Q79 = factor(Q79, levels = c(1, 2), labels = c("Stoma: No", "Stoma: Yes")), 
           ANA_SET=1)

theme_gtsummary_journal(journal = "jama")
theme_gtsummary_compact()

PC_dem2 <- PC_dem %>% 
  dplyr::select( Surgery_Dem, age, kon, 
                 hjart_karl_sjd, diabetes,  Comorbidity, kol, CVD, immunosuppr, sepsis_op, 
                 smoke_status2, auditc, BMI, Q16,  Q32, Q33, Q34, Q116, Q118, Q79,
                 civil, syssel) 

PC_demC <- PC_dem2 %>% 
  dplyr::select(Surgery_Dem, age, kon, 
                diabetes,  Comorbidity, kol, CVD, immunosuppr, sepsis_op)
TBL1C <- PC_demC %>% 
  tbl_summary( by = Surgery_Dem ,
               label =  list( kon ~ "Sex" ,  
                              age ~ "Age" ,
                              Comorbidity ~"Comorbidity", 
                              kol ~"COPD", 
                              CVD ~"Cardiovascular disease",  
                              diabetes~"Diabetes", 
                              immunosuppr~"Immunosuppressed", 
                              sepsis_op~"Sepsis at surgery"  ) ,
               
               type = list(c(kon, Comorbidity, kol, CVD,  
                             diabetes, immunosuppr, sepsis_op,) ~ "categorical" ,
                           c( age) ~ "continuous") ,
               statistic = list(all_continuous() ~ "{median} ({p25}, {p75})" , 
                                all_categorical() ~ "{n} / {N} ({p}%)") ) %>%
  add_overall()   %>%
  italicize_levels() %>%
  bold_labels()

PC_demD <- PC_dem2 %>% 
  filter(Surgery_Dem %in% c("Resection", "Lavage")) %>%
  mutate(Surgery_Dem = fct_drop(Surgery_Dem, only=c("N.R Resection", "N.R Lavage"))) %>%
  dplyr::select(Surgery_Dem,  smoke_status2, auditc, BMI, Q16,  Q32, Q33, Q34, Q116, Q118, Q79,
                civil, syssel) %>%
  mutate(Q32=fct_drop(Q32), Q33=fct_drop(Q33), Q34=fct_drop(Q34))

TBL1D <- PC_demD %>% 
  tbl_summary( by = Surgery_Dem ,
               label =  list( civil ~ "Marital status", 
                              syssel ~ "Employment status",
                              Q118~"Education",
                              Q116~"Residence",
                              smoke_status2~"Smoking status", 
                              auditc~"AUDIT-C score",
                              BMI ~"Body Mass Index(BMI)", 
                              Q79~"Stoma",
                              Q16~"Physical activity", 
                              Q32 ~"Cortisone treatment", 
                              Q33~"Rheumatoid arthritis", 
                              Q34~"Crohn's disease or ulcerative colitis"  ) ,
               type = list(c(civil, syssel, Q118, Q116, smoke_status2, Q79, 
                             Q16, Q32, Q33, Q34) ~ "categorical" ,
                           c( BMI, auditc ) ~ "continuous") ,
               statistic = list(all_continuous() ~ "{median} ({p25}, {p75})" , 
                                all_categorical() ~ "{n} / {N} ({p}%)") ) %>%
  add_overall()   %>%
  italicize_levels() %>%
  bold_labels()

PC_demCC <- PC_dem2 %>% 
  dplyr::select(age, kon, Q79,
                hjart_karl_sjd, diabetes,  Comorbidity, kol, CVD, immunosuppr, sepsis_op)
TBL1E <- PC_demCC %>% 
  tbl_summary( by = Q79 ,
               label =  list( kon ~ "Sex" ,  
                              age ~ "Age" ,
                              Comorbidity ~"Comorbidity", 
                              kol ~"COPD", 
                              CVD ~"Cardiovascular disease", 
                              hjart_karl_sjd~"Cardiovas dis", 
                              diabetes~"Diabetes", 
                              immunosuppr~"Immunosuppressed", 
                              sepsis_op~"Sepsis at surgery"  ) ,
               
               type = list(c(kon, Comorbidity, kol, CVD, hjart_karl_sjd, 
                             diabetes, immunosuppr, sepsis_op,) ~ "categorical" ,
                           c( age) ~ "continuous") ,
               statistic = list(all_continuous() ~ "{median} ({p25}, {p75})" , 
                                all_categorical() ~ "{n} / {N} ({p}%)") ) %>%
  add_overall()   %>%
  italicize_levels() %>%
  bold_labels()

#---------------------------------------------------------------------------------------------
# Questionnaire data 
as210406_komp1 <- tibble(
  read_sas("T:/SSORG/LapLAV/Data/Frågeformulär/as210406_komp1/as210406_komp1.sas7bdat",  
           NULL)) %>%
  dplyr::select(studienummer, Q4, Q5X0, Q5X1, Q5X2, Q5X3A, Q6) %>%
  mutate(fQ4 = factor(Q4, levels=c(1, 2, 3), labels=c("             Yes", 
                                                      "              No", 
                                                      "      Don´t know")),
         Q5 = case_when(Q5X0==2~1,Q5X1==2~2, Q5X2==2~3, Q5X3A==2~4), 
         fQ5 = factor(Q5, levels=c(1, 2, 3, 4), 
                      labels=c("       Not worried", 
                               "           Relapse", 
                               "Additional surgery", 
                               "            Other")), 
         fQ6 = factor(Q6, levels=c(1, 2, 3), labels=c("             Yes", 
                                                      "              No", 
                                                      "      Don´t know"))) %>%
  full_join(PC) %>% 
  dplyr::select(studienummer, Q4, Q5, Q6, 
                fQ4, fQ5, fQ6, 
                Surgery, 
                age_cent, kol, CVD, immunosuppr, sepsis_op, w, ANA_SET) %>%
  filter(ANA_SET==1) %>% dplyr::select(-ANA_SET)

# Saltin Grimby and smoking
SG_CAT <- PRO %>% 
  dplyr::select(studienummer, 
                Q16, Q17X1, Q17X2, Q17X3, Q17X4, Q17X5) %>%
  mutate(Smoking=case_when(Q17X1==2~1, Q17X2==2~2, 
                           Q17X3==2 | Q17X4==2 | Q17X5==2~3),
         fSmoking=factor(Smoking,levels=c(1, 2, 3), labels=c("      Not smoker", 
                                                             "     Past smoker", 
                                                             "  Current smoker")),
         fQ16=factor(Q16, levels=c(1, 2, 3, 4), 
                     labels=c("Mostly sedentary", 
                              "Light to moderate", 
                              "Regular heavier", 
                              "Hard training"))  ) %>%
  dplyr::select(studienummer, Q16, Smoking, fQ16, fSmoking) %>%
  full_join(PC) %>% 
  filter(ANA_SET==1) %>%
  dplyr::select(studienummer, Surgery, age_cent, kol, CVD, immunosuppr, sepsis_op, w,
                Q16, Smoking, fQ16, fSmoking) 

# Quality of life and well-being 
QoL_CONT <- PRO %>% 
  dplyr::select(studienummer, 
                Q7, Q8, Q9, Q10, Q11, Q12, Q13) %>%
  full_join(PC) %>% 
  dplyr::select(studienummer, Surgery, age_cent, kol, CVD, immunosuppr, sepsis_op, w,
                Q7, Q8, Q9, Q10, Q11, Q12, Q13, ANA_SET) %>%
  filter(ANA_SET==1) %>%
  pivot_longer(!c(studienummer, Surgery, 
                  age_cent, kol, CVD, immunosuppr, sepsis_op, w), 
               names_to = "variable", 
               values_to = "Response.n") %>% 
  mutate(variable=factor(variable, 
                         levels = c("Q7", "Q8", "Q9", "Q10", "Q11", "Q12", "Q13"),
                         labels = c("Quality of life", 
                                    "Life felt meaningful", 
                                    "Bodily stength", 
                                    "Feel depressed", 
                                    "Worry or anxiety", 
                                    "Mental well-being", 
                                    "Bodily health")  )) %>%
  replace_with_na_at( .vars = c("Response.n"),
                      condition = ~.x %in% c(777, 888, 999))

# EQ5-D-5L
EQ5D <- read_rds("S://RCVR/kod/LapLavQuestionnaire/LapLav.EQ.5D.5L210920.rds")
EQ5D <- full_join(PC, EQ5D) %>% 
  dplyr::select(studienummer, Surgery, age_cent, kol, CVD, immunosuppr, sepsis_op, w,
                MO, SC, UA, PD, AD, VAS, SCORE. , ANA_SET) %>%
  filter(ANA_SET==1)
LB <- c( "No", "Slight", "Moderate", "Severe", "Unable" )
LB2 <- c( "No", "Slight", "Moderate", "Severe", "Extreme" )
LB3 <- c( "No", "Slightly", "Moderately", "Severely", "Extremely" )
EQ5D_CAT <- EQ5D %>% dplyr::select(studienummer, Surgery, 
                                   age_cent, kol, CVD, immunosuppr, sepsis_op, w,
                                   MO, SC, UA, PD, AD) %>%
  mutate(fMO=factor(MO, levels=c(1, 2, 3, 4, 5), labels=LB),
         fSC=factor(SC, levels=c(1, 2, 3, 4, 5), labels=LB),
         fUA=factor(UA, levels=c(1, 2, 3, 4, 5), labels=LB),
         fPD=factor(PD, levels=c(1, 2, 3, 4, 5),  labels=LB2), 
         fAD=factor(AD, levels=c(1, 2, 3, 4, 5),  labels=LB3)) 
EQ5D_CONT <- EQ5D %>% 
  dplyr::select(studienummer, Surgery, 
                age_cent, kol, CVD, immunosuppr, sepsis_op, w,
                VAS, SCORE.) %>%
  pivot_longer(!c(studienummer, Surgery, 
                  age_cent, kol, CVD, immunosuppr, sepsis_op, w), 
               names_to = "variable", 
               values_to = "Response.n") %>% 
  mutate(variable=factor(variable, 
                         levels = c("SCORE.", "VAS"),
                         labels = c("EQ-5D indice", "Visual Analogue Scale")  ))  %>%
  replace_with_na_at( .vars = c("Response.n"),
                      condition = ~.x %in% c(777, 888, 999))

# LARS 
P <- PRO %>% 
  dplyr::select(studienummer, 
                Q80, Q81, Q82, 
                Q83, Q84, Q90)                

AND <- AND_0 %>%
  dplyr::select(studienummer, larssum) %>%
  mutate( LARS_cat = as.numeric(cut(larssum, breaks=c(0, 20, 29, 100), right=T)), 
          fLARS_cat = factor(LARS_cat, 
                             labels=c("          No LARS", 
                                      "       Minor LARS", 
                                      "       Major LARS")))
AND2 <- full_join(PC, AND) %>%
  full_join(P) %>%
  filter(Q79!="Stoma: Yes" & ANA_SET==1)
#---------------------------------------------------------------------------
LARS_CAT <- AND2 %>%  
  mutate(fQ80=factor(Q80, levels=c(1, 2, 3), 
                     labels=c("              No, never", 
                              "Less than once per week", 
                              "At least once per week")), 
         fQ81=factor(Q81, levels=c(1, 2, 3), 
                     labels=c("              No, never", 
                              "Less than once per week", 
                              "At least once per week")),
         fQ82=factor(Q82, levels=c(1, 2, 3, 4), 
                     labels=c( "       At least 7/24h", 
                               "        4-7 times/24h", 
                               "        1-3 times/24h", 
                               "      Less than 1/24h")),
         fQ83=factor(Q83, levels=c(1, 2, 3), 
                     labels=c("            No, never", 
                              "     Less than weekly", 
                              "At least once per week")), 
         fQ84=factor(Q84, levels=c(1, 2, 3, 4), 
                     labels=c("            No, never", 
                              "Less than once per week", 
                              "At least once per week", 
                              "                 Daily")), 
         fQ90=factor(Q90, levels=c(1, 2, 3, 4, 5), 
                     labels=c( "      Not applicable" , 
                               "                  No" , 
                               "            Slightly" , 
                               "          Moderately" , 
                               "           Very much"))) %>%
  dplyr::select(studienummer, Surgery, 
                age_cent, kol, CVD, immunosuppr, sepsis_op, w,
                Q80, Q81, Q82, Q83, Q84, Q90,
                fQ80, fQ81, fQ82, fQ83, fQ84, fQ90, fLARS_cat, LARS_cat ) %>%
  mutate(fQ84=forcats::fct_collapse(fQ84, 
                                    "At least once per week" = c("At least once per week", 
                                                                 "                 Daily")))
LARS_CONT <- AND2  %>% 
  dplyr::select(studienummer, larssum, Surgery, 
                age_cent, kol, CVD, immunosuppr, sepsis_op, w) %>%
  pivot_longer(!c(studienummer, Surgery, 
                  age_cent, kol, CVD, immunosuppr, sepsis_op, w), 
               names_to = "variable", 
               values_to = "Response.n") %>% 
  mutate(variable=factor(variable, 
                         levels = c("larssum"),
                         labels = c("LARS Score")  ))   %>%
  replace_with_na_at( .vars = c("Response.n"),
                      condition = ~.x %in% c(777, 888, 999))

# Merge all the data sets
EL_CAT <- full_join(EQ5D_CAT , LARS_CAT) %>% 
  full_join(as210406_komp1) %>% 
  full_join(SG_CAT)
EL_CONT <-  full_join(EQ5D_CONT , LARS_CONT) %>%
  full_join(QoL_CONT) 
# Estimate proportional odds model
# EQ5-D: High levels indicate more poor QRQoL
# LARS: High levels indicate more poor bowel health

# EQ5-D: High levels indicate more poor QRQoL
# LARS: High levels indicate more poor bowel health
EL_CAT2 <- EL_CAT %>%
  dplyr::select(studienummer, Surgery, 
                age_cent, kol, CVD, immunosuppr, sepsis_op, w,
                MO, SC, UA, PD, AD, fMO, fSC, fUA, fPD, fAD,
                Q80, Q81, Q82, Q83, Q84, Q90, LARS_cat, fQ80, fQ81, fQ82, fQ83, fQ84, fQ90,        
                fLARS_cat, 
                Q4, Q5, Q6, fQ4, fQ5, fQ6, 
                Q16, Smoking, fQ16, fSmoking) 

f <- EL_CAT2 %>%
  pivot_longer(cols=c(fMO, fSC, fUA, fPD, fAD,
                      fQ80, fQ81, fQ82, fQ83, fQ84, fQ90, fLARS_cat, 
                      fQ4, fQ5, fQ6, 
                      fQ16, fSmoking), 
               names_to = "Var.f", 
               values_to = "Response.f") %>% 
  mutate(Var.n= substring(Var.f,2)) 
n <- EL_CAT2 %>%
  pivot_longer(cols=c(MO, SC, UA, PD, AD,
                      Q80, Q81, Q82, Q83, Q84, Q90, LARS_cat,
                      Q4, Q5, Q6, 
                      Q16, Smoking), 
               names_to = "Var.n", 
               values_to = "Response.n") 
EL_CAT2 <- full_join(n, f) %>%
  dplyr::select(studienummer, Surgery, 
                age_cent, kol, CVD, immunosuppr, sepsis_op, w,
                Var.n, Response.n, Response.f) %>%
  replace_with_na_at( .vars = c("Response.n"),
                      condition = ~.x %in% c(777, 888, 999)) %>%
  mutate(variable=factor(Var.n,  levels = c("MO", "SC", "UA", "PD", "AD", 
                                            "Q80", "Q81", "Q82", "Q83", "Q84", "Q90", 
                                            "Q4", "Q5", "Q6", 
                                            "Q16", "Smoking",  "LARS_cat"), 
                         labels=c("Mobility", "Self-care", "Usual activities", 
                                  "Pain/discomfort", "Anxiety/depression", 
                                  "Control of flatus (wind)", 
                                  "Leakage of liquid stool", 
                                  "How often open bowels",                    
                                  "Open bowels again within 1 h",     
                                  "Rush to the toilet", 
                                  "Distress due to bowel dysfunction",
                                  "Worried after surgery", 
                                  "What worried the most", 
                                  "Satisfied with the treatment", 
                                  "Physical activity", 
                                  "Smoking", 
                                  "LARS categories"))) 

#______________________________________________________________________________
# Proportional odds regression using the propensity score weights

EL_CAT.CONT <- full_join(EL_CAT2, EL_CONT)  %>%                       
  dplyr::select(studienummer, Surgery, kol,
                age_cent, CVD,  immunosuppr, sepsis_op, w, 
                Response.n , variable) %>% 
  drop_na(c(Response.n, variable)) %>%
  group_by(variable) %>% 
  nest() %>% 
  mutate(lm_obj = map(data, ~ clm(factor(Response.n) ~  
                                    Surgery  ,
                                  data = . , weights = w ))) %>%
  mutate(lm_tidy = map(lm_obj, broom::tidy)) %>%
  ungroup() %>%
  transmute(variable, lm_tidy) %>%
  unnest(cols = c(lm_tidy)) %>%
  filter(term=="SurgeryLavage") %>%
  mutate(L=estimate-1.96*std.error, 
         U=estimate+1.96*std.error, 
         output=paste0(round(exp(estimate), digits = 2),
                       "(95%CI:",round(exp(L), digits = 2),";",round(exp(U), digits = 2),")" )) %>%
  dplyr::select(variable, output) %>%
  mutate(output_=case_when(variable %in% 
                             c("Worried after surgery", 
                               "What worried the most", 
                               "Satisfied with the treatment", "Smoking")~"", 
                           variable %notin% c("Worried after surgery", 
                                              "What worried the most", 
                                              "Satisfied with the treatment", "Smoking")~output)) 

# Odds ratio: Lavage vs Resection
#_____________________________________________________________________________

DESC_ <- EL_CAT2 %>%
  drop_na(Response.f) %>%
  group_by(variable, Surgery,  Response.f) %>%
  summarise(n = n() ) %>%
  group_by(variable, Surgery) %>%
  mutate(Percent=100*(n/sum(n)),
         Numbers_perc = paste0(n, "(", signif(Percent, 2), "%)"),
         Low = 100*scimp_wald(inpmat=n, alpha=0.05)[[4]], 
         Up = 100*scimp_wald(inpmat=n, alpha=0.05)[[5]]) %>%
  dplyr::select(variable, Surgery, Response.f, n, Percent, Numbers_perc, Low, Up)  %>%
  full_join(EL_CAT.CONT)

#___________________________________________________________________________________
# Bar plots
MDL = function(z) { A <- DESC_ %>% 
  filter(variable == z) %>%
  ungroup %>%
  drop_na() %>%
  mutate(Response.f2=factor(Response.f, ordered = TRUE), 
         Response.f2=fct_drop(Response.f2))

EST <- unique(A$output_)
PLOTT_1 <- ggplot(A, aes(x= Response.f2, y=Percent, fill=Surgery)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Low, ymax=Up), width=.2,
                position=position_dodge(.9)) + 
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) 
PLOTT_1 <- PLOTT_1+labs(title="", x="", y = " ", subtitle = paste0(z,"\n" , EST))+
  theme_classic() +
  scale_fill_manual(values=c('#999999','#E69F00')) +
  theme(legend.position="none") + 
  theme(plot.subtitle = element_text(size = 8),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 7.0, face="bold"),
        axis.title.x = element_text(color = "white", size = 1),
        axis.text.y = element_text(color = "black", size = 5), 
        axis.title.y = element_text(color = "white", size = 1)) 
PLOTT_1
}

EQ_1 <- MDL(z = "Mobility")
EQ_2 <- MDL(z = "Self-care")       
EQ_3 <- MDL(z = "Usual activities")
EQ_4 <- MDL(z = "Pain/discomfort")
EQ_5 <- MDL(z = "Anxiety/depression")  

LARS_1 <- MDL(z = "Control of flatus (wind)")
LARS_2 <- MDL(z = "Leakage of liquid stool")       
LARS_3 <- MDL(z = "How often open bowels")
LARS_4 <- MDL(z = "Open bowels again within 1 h")
LARS_5 <- MDL(z = "Rush to the toilet")  
LARS_6 <- MDL(z = "LARS categories")  
LARS_7 <- MDL(z = "Distress due to bowel dysfunction")

WOR_1 <- MDL(z = "Worried after surgery")  
WOR_2 <- MDL(z = "What worried the most")
WOR_3 <- MDL(z = "Satisfied with the treatment")

Phys_1 <- MDL(z = "Physical activity")  
Phys_2 <- MDL(z = "Smoking")  
#____________________________________________________________________________________  
# Violin plots for VAS, EQ5-D Score, LARS score and the Likert scales 

EL_CONT2 <- EL_CONT %>% full_join(EL_CAT.CONT) %>%
  na.omit()
MDLcont = function(z) { A <- EL_CONT2 %>% 
  filter(variable == z) %>%
  drop_na() 
EST <- unique(A$output_)
PLOTT_1 <- ggplot(A, aes(x= Surgery, y=Response.n, fill=Surgery)) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1, fill="white") 
PLOTT_1 <- PLOTT_1+labs(title="", x="", y = " ", subtitle = paste0(z, "\n", EST))+
  theme_classic() +
  scale_fill_manual(values=c('#999999','#E69F00')) +
  theme(legend.position="none") + 
  theme(plot.subtitle = element_text(size = 8),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 7.0, color="white"),
        axis.title.x = element_text(color = "white", size = 1),
        axis.text.y = element_text(color = "black", size = 5), 
        axis.title.y = element_text(color = "white", size = 1))
PLOTT_1
}

EQcont_1 <- MDLcont(z = "EQ-5D indice")
EQcont_2 <- MDLcont(z = "Visual Analogue Scale") 
LARScont_1 <- MDLcont(z = "LARS Score") 

QoL_1 <- MDLcont(z = "Quality of life") 
QoL_2 <- MDLcont(z = "Life felt meaningful") 
QoL_3 <- MDLcont(z = "Bodily stength") 
QoL_4 <- MDLcont(z = "Feel depressed")   
QoL_5 <- MDLcont(z = "Worry or anxiety") 
QoL_6 <- MDLcont(z = "Mental well-being")   
QoL_7 <- MDLcont(z = "Bodily health")  

#-----------------------------------------------------------------------------------------
# Combine plots
FIG_EQ5d <-  ggarrange( EQcont_1, EQ_1, EQ_2, EQ_3, EQ_4, EQ_5, EQcont_2 , 
                        ncol=4, nrow=2, common.legend = TRUE, legend="top", heights=c(1, 1))

FIG_LARS <-  ggarrange( LARS_7, LARScont_1, LARS_6, LARS_1, LARS_2, LARS_3, LARS_4 , LARS_5, 
                        ncol=4, nrow=2, common.legend = TRUE, legend="top", heights=c(1, 1))  

FIG_WOR <-  ggarrange( WOR_1, WOR_2, WOR_3, 
                       ncol=3, nrow=1, common.legend = TRUE, legend="top", heights=c(1, 1)) 

FIG_QoL <-  ggarrange( QoL_1, QoL_2, QoL_3, QoL_4, QoL_5, QoL_6, QoL_7, 
                       ncol=4, nrow=2, common.legend = TRUE, legend="top", heights=c(1, 1)) 

FIG_Phys <-  ggarrange( Phys_1, Phys_2, 
                        ncol=2, nrow=1, common.legend = TRUE, legend="top", heights=c(1, 1)) 


# Primary objective: Comparison between Laparoscopic lavage and Resection with regard to bother due to bowel dysfunction                            
# PRIMARY ANALYSIS: Statistical analysis of the effect of Lavage vs Resection on distress due to bowel dysfunction                                                       
# Primary endpoint: Bother/distress due to bowel dysfunction                               

PRIM <- PC %>%
  filter(Q79 == "Stoma: No") 
Ppp <- PRO %>% 
  dplyr::select(studienummer, Q90)

PRIM <- full_join(PRIM, Ppp) %>% 
  replace_with_na_at( .vars = c("Q90"),
                      condition = ~.x %in% c(777, 888, 999)) %>%
  mutate(fQ90=factor(Q90, levels=c(1, 2, 3, 4, 5), 
                     labels=c( "       N/A" , 
                               "        No" , 
                               "  Slightly" , 
                               "Moderately" , 
                               " Very much"))) %>%
  filter(Q79 == "Stoma: No") 

DESC_1 <- PRIM %>%
  drop_na(fQ90) %>%
  tabyl(fQ90, Surgery) %>% 
  adorn_percentages("col") %>%
  adorn_pct_formatting(digits = 2) %>%
  adorn_ns() %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Cambria", position = "left")


#The odds ratio for a tendency of responding with higher 
# categories(much distress) rather than lower categories (none or low distress)               

M0 <- clm(factor(Q90) ~  Surgery , data = PRIM  )
Unadj1 <- paste( signif(exp(M0$beta[[1]]), digits = 3),
                 "(95%CI:",signif(exp(confint(M0)[[1,1]]), digits = 3),
                 ";",signif(exp(confint(M0)[[1,2]]), digits = 3),
                 "),p=", signif(   anova(M0, type="III")[[1,3]] , digits = 2) )

M1 <- clm(factor(Q90) ~  Surgery + rcs(age_cent ,5)  +
            kol +
            CVD +
            immunosuppr +
            sepsis_op, data = PRIM , weights = w )
Adj1 <- paste( signif(exp(M1$beta[[1]]), digits = 3),
               "(95%CI:",signif(exp(confint(M1)[[1,1]]), digits = 3),
               ";",signif(exp(confint(M1)[[1,2]]), digits = 3),
               "),p=", signif(   anova(M1, type="III")[[1,3]] , digits = 2) )

# The tendency for bother due to bowel dysfunction is higher among patients 
# operated by Laparoscopic lavage compared with resection.

# Stoma function                    

ST <- read_delim("LpLvQ211012.txt",  "\t", escape_double = FALSE, trim_ws = TRUE)
ST2 <- full_join(PC, ST) %>%
  full_join(P) %>%
  filter(Q79=="Stoma: Yes" & ANA_SET==1)

lab1 <- c("No" , "At least 1/month" , "At least 1/week" , 
          "At least 3/week" , "At least daily" )
lab2 <-c("Not applicable" , "No" , "Slightly" , "Moderately" , "Very much")
lab3<-c("Not true", "Partly true", "Largely true", "Fully true")

ST2 <- ST2 %>%
  mutate(Q93=factor(Q93, levels=c(1, 2, 3), labels=c("Ileostom", "Colostomy", "Don´t know")), 
         Q94X1=factor(Q94X1, levels=c(1, 2, 3), 
                      labels=c("Had it before perforated diverticulitis", 
                               "During diagnosis and treatment of perforated diverticulitis", 
                               "After acute treatment of perforated diverticulitis")), 
         Q95=factor(Q95, levels=c(1, 2, 3), labels=c("No", "Yes", "Don´t know/Can´t remember")), 
         Q96=factor(Q96, levels=c(1, 2, 3, 4, 5), labels=lab1), 
         Q97=factor(Q97, levels=c(0, 1, 2, 3, 4), labels=lab2 ),
         Q98=factor(Q98, levels=c(1, 2, 3, 4, 5), labels=lab1),
         Q99=factor(Q99, levels=c(0, 1, 2, 3, 4), labels=lab2 ),
         Q100=factor(Q100, levels=c(1, 2, 3, 4, 5), labels=lab1),
         Q101=factor(Q101, levels=c(0, 1, 2, 3, 4), labels=lab2 ),
         Q102=factor(Q102, levels=c(0, 1, 2, 3, 4), 
                     labels=c("N/A. Not worried" , "Atleast 1/month" , "Atleast 1/week" , 
                              "Atleast 3/week" , "Atleast daily" ) ),
         Q103=factor(Q103, levels=c(1, 2, 3, 4, 9), 
                     labels=c("No" , "Slightly" ,  "Moderately" , "Very much", "Don´t know")),
         Q104=factor(Q104, levels=c(1, 2, 3, 4), 
                     labels=c("No" , "Slightly" ,  "Moderately" , "Very much")),
         Q105=factor(Q105, levels=c(1, 2, 3), labels=c("No", "Yes", "Don´t know")), 
         Q106=factor(Q106, levels=c(1, 2, 3), labels=c("No", "Yes", "Don´t know")),
         Q107=factor(Q107, levels=c(0, 1, 2, 3, 4), labels=lab2 ),
         Q108A=factor(Q108A, levels=c(1, 2, 3, 4), labels=lab3 ), 
         Q108B=factor(Q108B, levels=c(1, 2, 3, 4), labels=lab3 ),
         Q108C=factor(Q108C, levels=c(1, 2, 3, 4), labels=lab3 ), 
         Q108D=factor(Q108D, levels=c(1, 2, 3, 4), labels=lab3 ),  
         Q108E=factor(Q108E, levels=c(1, 2, 3, 4), labels=lab3 ) ) %>%
  dplyr::select(studienummer, Surgery, 
                Q93, Q94X1, Q95, Q96, Q97, Q98, Q99, 
                Q100, Q101, Q102, Q103, Q104, Q105, Q106, Q107, 
                Q108A, Q108B, Q108C, Q108D, Q108E)

ST3 <- ST2 %>%
  pivot_longer(cols=-c(studienummer, Surgery) ,
               names_to = "Var.f", 
               values_to = "Response.f") %>% 
  drop_na(Response.f) %>% 
  mutate(variable=factor(Var.f,  
                         levels = c("Q93", "Q94X1", "Q95", "Q96", "Q97", "Q98", "Q99", 
                                    "Q100", "Q101", "Q102", "Q103", "Q104", "Q105", "Q106", "Q107", 
                                    "Q108A", "Q108B", "Q108C", "Q108D", "Q108E"), 
                         labels=c("Type of stoma",
                                  "When did you get a stoma?",
                                  "Feel involved in the decision to have a stoma",
                                  "Had loud gas departures from the stoma, last month",
                                  "Bother due to loud gas departures",
                                  "Had foul-smelling gas escapes from the stoma, last month", 
                                  "Bother due to foul-smelling gas escapes",
                                  "Had a leak of stool from the stoma, last month?",
                                  "Bother due to leak of stool",
                                  "Worried about the stoma leaking, last month", 
                                  "Skin around the stoma been irritated, last month?",
                                  "Problems taking care of your stoma, last month?",
                                  "Further surgery due to problems with the stoma",
                                  "You or a doctor noticed a bulge (hernia) at the stoma",
                                  "Bother due to overall problems with the stoma", 
                                  "I can live a full life with my stoma",
                                  "I feel comfortable with my stoma",
                                  "Worried something embarrassing can happen \n during sexual activity due to my stoma",
                                  "I feel dirty and unclean due to my stoma",
                                  "I have the leisure activities and social life I desire" ))) %>%
  group_by(variable, Surgery,  Response.f) %>%
  summarise(n = n() ) %>%
  group_by(variable, Surgery) %>%
  mutate(Percent=100*(n/sum(n)),
         Numbers_perc = paste0(n, "(", signif(Percent, 2), "%)")) %>%
  dplyr::select(variable, Surgery, Response.f, Numbers_perc)  

ST4 <- ST3 %>%
  pivot_wider(names_from="Surgery", values_from="Numbers_perc") %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Cambria", position = "left")

# Presentation of secondary results in the supplement                                 

spl <- PRO %>% 
  dplyr::select(studienummer, 
                Q46 , Q47 , Q48 , Q49 , Q50 , Q51 , Q52 , Q53 , Q54,   # Urin
                Q55 , Q56 , Q57 , Q58 , Q59 , Q60 , Q61 , Q62,         # Sex all 
                Q64 , Q65 , Q66 , Q67 , Q68 , Q69 , Q70 , Q71,         # Sex men
                Q72, Q73 , Q74 , Q75 , Q76 , Q77 , Q78                # Sex women 
  ) %>%
  full_join(PC) %>% 
  filter(ANA_SET==1) 
#__________________________Urine_______________________________________________________________

lab1 <- c("          Never" , 
          "Less than half " , 
          "  Half the time" , 
          "More than half " , 
          "         Always" )
lab2 <- c("  Not applicable", lab1)
lab3 <- c("Not applicable", 
          "           No", 
          "          Yes")
lab4 <- c("   Not applicable" , 
          "              No" , 
          "At least 1/month" , 
          " At least 1/week" , 
          "At least 2-3/week" , 
          "           Daily" )
lab5 <- c("            Yes" , 
          "             No" , 
          "     Don't know")
lab6 <- c("  Not applicable" , 
          "Less than 1/day" , 
          "     Once a day" , 
          "       2-3 /day" , 
          "       4-5 /day" , 
          "  6 or more/day")
lab7 <-c("Not applicable" , 
         "           No" , 
         "     Slightly" , 
         "   Moderately" , 
         "    Very much")
H <- spl %>%
  dplyr::select(studienummer , Q46 , Q47 , Q48 , Q49 , Q50 , Q51 , Q52 , Q53 , Q54) %>%
  replace_with_na_all(condition = ~.x == 999) %>%
  replace_with_na_all(condition = ~.x == 777) %>%
  mutate(Q46 = factor(Q46, levels=c(1, 2, 3, 4, 5), labels = lab1)) %>%
  mutate(Q47 = factor(Q47, levels=c(1, 2, 3, 4, 5), labels = lab1)) %>%
  mutate(Q48 = factor(Q48, levels=c(0, 1, 2, 3, 4, 5), labels = lab2)) %>%
  mutate(Q50 = factor(Q50, levels=c(0, 1, 2, 3, 4, 5), labels = lab2)) %>%
  mutate(Q49 = factor(Q49, levels=c(0, 1, 2),  labels = lab3)) %>%
  mutate(Q51 = factor(Q51, levels=c(0, 1, 2, 3, 4, 5), labels = lab4)) %>%
  mutate(Q52 = factor(Q52, levels = c(1, 2, 3), labels = lab5)) %>%
  mutate(Q53 = factor(Q53, levels=c(0, 1, 2, 3, 4, 5), , labels = lab6)) %>%
  mutate(Q54 = factor(Q54, levels = c(0, 1, 2, 3, 4), labels = lab7))

#__________________________Sex all____________________________________
lab11<-c("Yes" , 
         "No")
lab12<- c("            No" , 
          "  Yes, Slightly" , 
          "Yes, moderatley" , 
          " Yes, very much")
lab13<- c("            Never" , 
          "  About 1-2/month" , 
          "  About 3-4/month" , 
          "    About 1/week" , 
          "More than 2/week")
lab14 <- c( "           Not applicable" , 
            "          No, not at all" , 
            "Yes, true to some extent" , 
            "       Yes, largely true" , 
            "    Yes, completely true" )
lab15 <- c("    N/A, could not\nbefore surgery" , 
           "   N/A, not tried\nafter surgery" , 
           "                  No" , 
           "Yes, more problems\nafter surgery" , 
           "     Yes, not able\nafter surgery")
I <- spl %>%
  dplyr::select(studienummer , Q55 , Q56 , Q57 , Q58 , Q59 , Q60 , Q61 , Q62) %>%
  replace_with_na_all(condition = ~.x == 999) %>%
  replace_with_na_all(condition = ~.x == 777) %>%
  mutate(Q55 = factor(Q55, levels=c(1, 2), labels = lab11)) %>%
  mutate(Q56 = factor(Q56,  levels=c(1, 2, 3, 4), labels = lab12)) %>%
  mutate(Q57 = factor(Q57,  levels=c(1, 2, 3, 4, 5), labels = lab13)) %>%
  mutate(Q58 = factor(Q58,  levels=c(1, 2, 3, 4), labels = lab12)) %>%
  mutate(Q59 = factor(Q59, levels=c(0, 1, 2, 3, 4), labels = lab14)) %>%
  mutate(Q60 = factor(Q60, levels=c(0, 9, 1, 2, 3), labels = lab15)) %>%
  mutate(Q61 = factor(Q61, levels=c(1, 2, 3, 4),  labels = lab12)) %>%
  mutate(Q62 = factor(Q62, levels=c(0, 1, 2, 3, 4), labels = lab7))

#__________________________Sex men________888=women_________________________________
lab21 <- c("N/A, could not \n before surgery" , 
           "               No" , 
           "              Yes" , 
           "       Don't know")
lab22 <- c(" Not applicable" , 
           " Not applicable" , 
           " Not applicable" , 
           "            No" , 
           " Yes, slightly" , 
           "Yes, moderately" , 
           "Yes, very much")
lab23<-c( "Not applicable" ,
          "Not applicable" , 
          "          No" , 
          "    Slightly" , 
          "  Moderately" , 
          "   Very much")
lab24 <- c("      Not applicable",
           "              Never" , 
           "Less than half\nthe time" ,  
           "More than half\nthe time" , 
           "          Always" )
J <- spl %>%
  dplyr::select(studienummer , Q64 , Q65 , Q66 , Q67 , Q68 , Q69 , Q70 , Q71) %>%
  replace_with_na_all(condition = ~.x == 999) %>%
  replace_with_na_all(condition = ~.x == 888) %>% 
  mutate(Q64 = factor(Q64, levels=c(0 ,1, 2, 3) , 
                      labels = lab21)) %>%
  mutate(Q65 = factor(Q65, levels=c(0, 8 , 9 , 1, 2, 3 , 4) , 
                      labels = lab22)) %>%
  mutate(Q66 = factor(Q66, levels=c(0, 8 , 9 , 1, 2, 3 , 4) , 
                      labels = lab22)) %>%
  mutate(Q67 = factor(Q67, levels=c( 1, 2 ,888) , labels = c("No" , "Yes" ,888))) %>%
  mutate(Q68 = factor(Q68, levels=c(0, 1 , 2 , 3 , 4) , 
                      labels = lab7)) %>%
  mutate(Q69 = factor(Q69, levels=c(0, 1 , 2 , 3 , 4) , 
                      labels = lab24)) %>%
  mutate(Q70 = factor(Q70, levels=c(0, 1 , 2) , 
                      labels = lab3)) %>%
  mutate(Q71 = factor(Q71, levels=c(0 , 9 , 1, 2, 3 , 4) , 
                      labels = lab23)) 

#__________________________Sex women________888=man_____________________________________ 

K <- spl %>%
  dplyr::select(studienummer ,  Q72, Q73 , Q74 , Q75 , Q76 , Q77 , Q78) %>%
  replace_with_na_all(condition = ~.x == 999) %>%
  replace_with_na_all(condition = ~.x == 888) %>%
  mutate(Q72 = factor(Q72, levels=c(1 , 2 , 3) , 
                      labels = lab5)) %>%
  mutate(Q73 = factor(Q73, levels=c( 0 , 1 , 2 , 3 , 4) , 
                      labels = lab24)) %>%
  mutate(Q74 = factor(Q74, levels=c(0, 1 , 2 , 3 , 4) , 
                      labels = lab7)) %>%
  mutate(Q75 = factor(Q75, levels=c(0, 1 , 2 , 3 , 4) , 
                      labels = lab7)) %>%
  mutate(Q76 = factor(Q76, levels=c(0, 1 , 2 , 3 , 4) , 
                      labels = lab7)) %>%
  mutate(Q77 = factor(Q77, levels=c(0 , 9 , 1 , 2 , 3 , 4) , 
                      labels = lab23)) %>%
  mutate(Q78 = factor(Q78, levels=c(0 , 9 , 1 , 2 , 3 , 4) , 
                      labels = lab23))
#___________________________________________________________________________________
HIJK <- H %>% 
  full_join(I) %>%
  full_join(J) %>%
  full_join(K) %>%
  pivot_longer(cols=-studienummer ,
               names_to = "Var.f", 
               values_to = "Response.f") %>% 
  drop_na(Response.f) %>% 
  mutate(variable=factor(Var.f,  
levels = c( "Q46", "Q47", "Q48", 
"Q49", "Q50", "Q51" ,
"Q52", "Q53", "Q54" ,
"Q55", "Q56", "Q57" ,
"Q58", "Q59", "Q60",
"Q61", "Q62", 
"Q64" ,"Q65", "Q66", 
"Q67" , "Q68", "Q69", 
"Q70" , "Q71", 
"Q72", "Q73" , "Q74", 
"Q75", "Q76","Q77", 
"Q78"), 
labels=c("Feeling that the bladder has\nnot been emptied, last month", 
"Need to hurry to the toilet to\nurinate, last month", 
"Urinary leakage due to not\nreaching a toilet in time", 
"Urinary leakage due to physical\nstress, last month",
"Urinary leakage, last month", 
"Have urinary dysfunction caused\nyou to refrain from activities,\nlast month", 
"Did you have problems with urinary\nleakage before your surgery for\nperforated diverticulitis", 
"How often do you change pad,\ndiaper or sanitary aid\nduring a typical day", 
"How would you feel if urinary dysfunction\nas it were during the last month were to\nremain for the rest of your life?", 

"Is intercourse part of your sex-life", 
"Has sex been important to you,\nlast 6 months",
"How often have you had intercourse\nor any other sexual activity, last 6 months",
"Have you had thoughts or longing\nfor sex, last 6 months", 
"Have you refrained from sexual activities\nout of fear of failure, last 6 months", 
"Have your ability to have an orgasm\nchanged after your surgery due to perforated\ndiverticulitis", 
"Are you satisfied with your current sex-life,\nlast 6 months", 
"How would you feel if sexual impairments\nwere to remain the same for the rest\nof your life?", 

"Have your ability to get an erection\ndeteriorated or ceased entirely after\nthe surgery for perforated diverticulitis", 
"If your ability to get an erection has\ndeteriorated or ceased entirely,did\nthis affect your self-esteem", 
"If your ability to get an erection has\ndeteriorated or ceased entirely, how would\nyou feel if the impairment were to remain the same\nfor the rest of your life?",
"Have you used technical aids or medication\nto improve or prolong your erection,\nlast 6 months", 
"If you used technical aids or medication to\nimprove or prolong your erection\nlast 6 months, did it help",
"Have you had premature ejaculation at\nsexual activity, last 6 months", 
"Have you been unable to ejaculate at\nsexual activity, last 6 months", 
"How would you feel if sexual impairments\nwere to remain the same for the\nrest of your life?", 

"Have you reached menopause?", 
"At sexual arousal, have labia, clitoris\nor the vulva felt swollen and full of blood,\nlast 6 months", 
"Have your vagina felt lubricated at\nsexual arousal, last 6 months", 
"Have you had superficial pain around vulva\nduring intercourse or comparable activity,\nlast 6 months",
"Have you had deep pain in your pelvic\nregion during intercourse or\ncomparable activity,last 6 months", 
"If you experienced pain during intercourse\nor comparable activity last 6 months, how\nwould you feel if the impairment were to remain the\nsame for the rest of your life?", 
"How would you feel if sexual impairments were\nto remain the same for the rest of your life?"))) 


DESC_ <- HIJK %>% 
  full_join(PC) %>% 
  dplyr::select(studienummer, variable, Surgery,  Response.f) %>% 
  group_by(variable, Surgery,  Response.f) %>%
  summarise(n = n() ) %>%
  group_by(variable, Surgery) %>%
  mutate(Percent=100*(n/sum(n)),
         Numbers_perc = paste0(n, "(", signif(Percent, 2), "%)"),
         Low = 100*scimp_wald(inpmat=n, alpha=0.05)[[4]], 
         Up = 100*scimp_wald(inpmat=n, alpha=0.05)[[5]]) %>%
  ungroup() %>%
  dplyr::select(variable, Surgery, Response.f, n, Percent, Numbers_perc, Low, Up)  

MDL = function(z) { A <- DESC_ %>% filter(variable == z) %>%
  drop_na() %>%
  ungroup()
PLOTT_1 <- ggplot(A, aes(x= Response.f, y=Percent, fill=Surgery)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Low, ymax=Up), width=.2,
                position=position_dodge(.9)) + 
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) 
PLOTT_1 <- PLOTT_1+labs(title="", x="", y = " ", subtitle = paste0(z))+
  theme_classic() +
  scale_fill_manual(values=c('#999999','#E69F00')) +
  theme(legend.position="none") + 
  theme(plot.subtitle = element_text(size = 6),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 7.0, face="bold"),
        axis.title.x = element_text(color = "white", size = 1),
        axis.text.y = element_text(color = "black", size = 5), 
        axis.title.y = element_text(color = "white", size = 1)) 
PLOTT_1
}
#------------------------------------------------------------------------------------
U_1 <- MDL(z = "Feeling that the bladder has\nnot been emptied, last month")
U_2 <- MDL(z = "Need to hurry to the toilet to\nurinate, last month") 
U_3 <- MDL(z = "Urinary leakage due to not\nreaching a toilet in time") 
U_4 <- MDL(z = "Urinary leakage due to physical\nstress, last month")
U_5 <- MDL(z = "Urinary leakage, last month") 
U_6 <- MDL(z = "Have urinary dysfunction caused\nyou to refrain from activities,\nlast month") 
U_7 <- MDL(z = "Did you have problems with urinary\nleakage before your surgery for\nperforated diverticulitis") 
U_8 <- MDL(z = "How often do you change pad,\ndiaper or sanitary aid\nduring a typical day") 
U_9 <- MDL(z = "How would you feel if urinary dysfunction\nas it were during the last month were to\nremain for the rest of your life?") 
#------------------------------------------------------------------------------------

S_1 <- MDL(z ="Is intercourse part of your sex-life") 
S_2 <- MDL(z = "Has sex been important to you,\nlast 6 months")
S_3 <- MDL(z = "How often have you had intercourse\nor any other sexual activity, last 6 months")
S_4 <- MDL(z = "Have you had thoughts or longing\nfor sex, last 6 months") 
S_5 <- MDL(z = "Have you refrained from sexual activities\nout of fear of failure, last 6 months") 
S_6 <- MDL(z = "Have your ability to have an orgasm\nchanged after your surgery due to perforated\ndiverticulitis") 
S_7 <- MDL(z = "Are you satisfied with your current sex-life,\nlast 6 months") 
S_8 <- MDL(z = "How would you feel if sexual impairments\nwere to remain the same for the rest\nof your life?")

#------------------------------------------------------------------------------------
SM_1 <- MDL(z ="Have your ability to get an erection\ndeteriorated or ceased entirely after\nthe surgery for perforated diverticulitis") 
SM_2 <- MDL(z ="If your ability to get an erection has\ndeteriorated or ceased entirely,did\nthis affect your self-esteem") 
SM_3 <- MDL(z ="If your ability to get an erection has\ndeteriorated or ceased entirely, how would\nyou feel if the impairment were to remain the same\nfor the rest of your life?")
SM_4 <- MDL(z ="Have you used technical aids or medication\nto improve or prolong your erection,\nlast 6 months" )
SM_5 <- MDL(z ="If you used technical aids or medication to\nimprove or prolong your erection\nlast 6 months, did it help")
SM_6 <- MDL(z ="Have you had premature ejaculation at\nsexual activity, last 6 months" )
SM_7 <- MDL(z ="Have you been unable to ejaculate at\nsexual activity, last 6 months") 
SM_8 <- MDL(z ="How would you feel if sexual impairments\nwere to remain the same for the\nrest of your life?") 
#------------------------------------------------------------------------------------
SF_1 <- MDL(z ="Have you reached menopause?") 
SF_2 <- MDL(z ="At sexual arousal, have labia, clitoris\nor the vulva felt swollen and full of blood,\nlast 6 months")
SF_3 <- MDL(z ="Have your vagina felt lubricated at\nsexual arousal, last 6 months") 
SF_4 <- MDL(z ="Have you had superficial pain around vulva\nduring intercourse or comparable activity,\nlast 6 months")
SF_5 <- MDL(z ="Have you had deep pain in your pelvic\nregion during intercourse or\ncomparable activity,last 6 months") 
SF_6 <- MDL(z ="If you experienced pain during intercourse\nor comparable activity last 6 months, how\nwould you feel if the impairment were to remain the\nsame for the rest of your life?") 
SF_7 <- MDL(z ="How would you feel if sexual impairments were\nto remain the same for the rest of your life?")
#------------------------------------------------------------------------------------ 
# Combine plots
FIG_U <-  ggarrange(U_1, U_2, U_3, U_4, U_5, U_6,
                    U_7, U_8, U_9, 
                    ncol=5, nrow=2, common.legend = TRUE, legend="top", heights=c(1, 1))
FIG_S <-  ggarrange(S_1, S_2, S_3, S_4, S_5, S_6,
                    S_7, S_8, 
                    ncol=4, nrow=2, common.legend = TRUE, legend="top", heights=c(1, 1))

FIG_SM <-  ggarrange(SM_1, SM_2, SM_3, SM_4, SM_5, SM_6,
                     SM_7, SM_8, 
                     ncol=4, nrow=2, common.legend = TRUE, legend="top", heights=c(1, 1))
FIG_SF <-  ggarrange(SF_1, SF_2, SF_3, SF_4, SF_5, SF_6,
                     SF_7, 
                     ncol=4, nrow=2, common.legend = TRUE, legend="top", heights=c(1, 1))
#_____________________________________________________________________________________________________________
#_____________________________________________________________________________________________________________

# Additional analysis 220333: Combine bother question for those with and without stoma 
#----------------------------------------------------------------------------
# Combine data

# Swedish version of the question:
#90. Om du resten av ditt liv skulle leva med dina sammantagna mag- och tarmbesvär, 
#som det varit den senaste månaden, hur skulle du uppleva det?
#0.		Inte aktuellt, jag har inte haft några tarmbesvär den senaste månaden
#1.		Det skulle inte beröra mig alls
#2.		Det skulle beröra mig lite
#3.		Det skulle beröra mig måttligt
#4.		Det skulle beröra mig mycket
lab <-c( "Not applicable" , "No" , "A little" , "Moderately" , "Very much")
PRIM2 <- PRIM %>% 
  mutate(foutcome=factor(Q90, levels=c(0, 1, 2, 3, 4), labels=lab ), 
         outcome=Q90, STOMA="N") %>% 
  dplyr::select(studienummer, outcome,  foutcome, STOMA, Surgery, 
                age_cent, kol, CVD, immunosuppr, sepsis_op,  w)    # Q90, fQ90
#------------------------------------------------------------------------------
ST2 <- full_join(PC, ST) %>%
  full_join(P) %>%
  filter(Q79=="Stoma: Yes" & ANA_SET==1)
ST2 <- ST2 %>%
  mutate(foutcome=factor(Q107, levels=c(0, 1, 2, 3, 4), labels=lab ), 
         outcome=Q107, STOMA="Y") %>% 
  dplyr::select(studienummer, outcome, foutcome, STOMA, Surgery, 
                age_cent, kol, CVD, immunosuppr, sepsis_op,  w)    # Q107

PS_ <- PRIM2 %>%
  full_join(ST2) %>%
  mutate(STOMA=factor(STOMA))
#----------------------------------------------------------------------------
# Describe
DESC_1 <- PS_ %>%
  drop_na(foutcome) %>%
  tabyl(foutcome, Surgery) %>% 
  adorn_percentages("col") %>%
  adorn_pct_formatting(digits = 2) %>%
  adorn_ns() %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Cambria", position = "left")
#---------------------------------------------------------------------------
PS_sub <- PS_ # %>%
#  filter(outcome != 0)
# Estimation
# Odds ratio: Lavage vs Resection
M0 <- clm(factor(outcome) ~  Surgery , data = PS_sub   )
Unadj1 <- paste( signif(exp(M0$beta[[1]]), digits = 3),
                 "(95%CI:",signif(exp(confint(M0)[[1,1]]), digits = 3),
                 ";",signif(exp(confint(M0)[[1,2]]), digits = 3),
                 "),p=", signif(   anova(M0, type="III")[[1,3]] , digits = 2) )


M01 <- clm(factor(outcome)  ~  Surgery , data = PS_sub , weights = w )
Adj01 <- paste( signif(exp(M01$beta[[1]]), digits = 3),
                "(95%CI:",signif(exp(confint(M01)[[1,1]]), digits = 3),
                ";",signif(exp(confint(M01)[[1,2]]), digits = 3),
                "),p=", signif(   anova(M01, type="III")[[1,3]] , digits = 2) )
#--------------------------------------------------------------------------
# Plot
DESC_ <- PS_ %>%
  drop_na(foutcome) %>%
  group_by(Surgery,  foutcome) %>%
  summarise(n = n() ) %>%
  group_by(Surgery) %>%
  mutate(Percent=100*(n/sum(n)),
         Numbers_perc = paste0(n, "(", signif(Percent, 2), "%)"),
         Low = 100*scimp_wald(inpmat=n, alpha=0.05)[[4]], 
         Up = 100*scimp_wald(inpmat=n, alpha=0.05)[[5]]) %>%
  dplyr::select(Surgery, foutcome, n, Percent, Numbers_perc, Low, Up)  

PLOTT_1 <- ggplot(DESC_, aes(x= foutcome, y=Percent, fill=Surgery)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Low, ymax=Up), width=.2,
                position=position_dodge(.9)) + 
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) 
PLOTT_1 <- PLOTT_1+labs(title="", x="", y = " ", 
                        subtitle = paste0("Distress due to bowel or stoma dysfunction","\n" , Adj01 )) +
  theme_classic() +
  scale_fill_manual(values=c('#999999','#E69F00')) +
  theme(legend.position="none") + 
  theme(plot.subtitle = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 7.0, face="bold"),
        axis.title.x = element_text(color = "white", size = 1),
        axis.text.y = element_text(color = "black", size = 5), 
        axis.title.y = element_text(color = "white", size = 1)) 
PLOTT_1

FIG_ <-  ggarrange( PLOTT_1,  ncol=1, nrow=1, common.legend = TRUE, legend="top", heights=c(1, 1))


