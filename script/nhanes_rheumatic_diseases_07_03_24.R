## Analysis script for NHANES paper ##

# Author: Gabriel P. Esteves

# In order for this script to be reproducible:
# Execute the .Rproject file, and then open this analysis script.

# Important: this script begins by loading a .Rdata object, which builds the dataframe
# utilized in the analysis. This dataframe was simply built by downloading the 
# data from NHANES cycles 2007-2018, reading these in R using the "haven" package,
# and joining all dataframes together using the "full_join" function from tidyverse.
# I decided to not include this code as it is not fundamental to the analysis and is 
# costly to run, and also would require the presence of all the NHANES datasets in the
# repository. If you have any difficulty downloading and loading NHANES datasets in R,
# check their online tutorials and sample code found in: https://wwwn.cdc.gov/nchs/nhanes/tutorials/

# load packages
library(tidyverse); library(here); library(survey); library(broom); 
library(writexl); library(dagitty); library(patchwork); library(glue); 
library(gtsummary); library(marginaleffects)

# no scientific notation
options(scipen=999)

#custom ggplot theme
theme_avp <- function() {
  theme_bw(base_size=10) +
    theme(axis.title.x = element_text(size = 9),
          axis.title.y = element_text(size = 9),
          panel.grid.minor = element_blank())
}

# load raw dataset
load(here("data/nhanes_raw_dataframe.Rdata"))

# creating and transforming variables before setting up survey design
# nutrients

rh_dxa_all <- nhanes_raw_dataframe |> 
  mutate(kcal = case_when(
    DRDINT.x == 2 ~ ((DR1TKCAL+DR2TKCAL)/2), # these lines calculate the mean between two dietary recalls, when two are available
    DRDINT.x == 1 ~ DR1TKCAL # or just pick the first one that is available, if only one is.
  )) |> 
  mutate(PTN = case_when(
    DRDINT.x == 2 ~ ((DR1TPROT+DR2TPROT)/2), # same for other nutrients.
    DRDINT.x == 1 ~ DR1TPROT
  )) |> 
  mutate(CHO = case_when(
    DRDINT.x == 2 ~ ((DR1TCARB+DR2TCARB)/2),
    DRDINT.x == 1 ~ DR1TCARB
  )) |> 
  mutate(LIP = case_when(
    DRDINT.x == 2 ~ ((DR1TTFAT+DR2TTFAT)/2),
    DRDINT.x == 1 ~ DR1TTFAT
  ))

# other relevant variables

rh_dxa_all <- rh_dxa_all |> 
  mutate(PTNgkg = PTN/BMXWT) |> 
  mutate(femur_popmean_male = 0.93) |> # mean and sd values extracted from NHANES to calculate t-scores
  mutate(femur_popsd_male = 0.144) |> 
  mutate(femur_popmean_female = 0.94) |> 
  mutate(femur_popsd_female = 0.122) |> 
  mutate(femur_t_score = case_when(
    RIAGENDR == "1" ~ ((DXXOFBMD-femur_popmean_male)/femur_popsd_male),
    RIAGENDR == "2" ~ ((DXXOFBMD-femur_popmean_female)/femur_popsd_female)
  )) |> 
  mutate(spine_popmean_male = 1.061) |> 
  mutate(spine_popsd_male = 0.115) |> 
  mutate(spine_popmean_female = 1.068) |> 
  mutate(spine_popsd_female = 0.118) |> 
  mutate(spine_t_score = case_when(
    RIAGENDR == "1" ~ ((DXXOSBMD-spine_popmean_male)/spine_popsd_male),
    RIAGENDR == "2" ~ ((DXXOSBMD-spine_popmean_female)/spine_popsd_female)
  )) |> 
  mutate(spine_osteoporosis = as.factor(case_when(
    spine_t_score >= -2.5 ~ "No",
    spine_t_score < -2.5 ~ "Yes"
  ))) |> 
  mutate(femur_osteoporosis = as.factor(case_when(
    femur_t_score >= -2.5 ~ "No",
    femur_t_score < -2.5 ~ "Yes"
  ))) |> 
  mutate(disease_cat = case_when(
    MCQ195 == "1" ~ "Osteoarthritis", # classifying rheumatic disease according to questionnaire responses
    MCQ195 == "2" ~ "Rheumatoid arthritis",
    MCQ195 == "3" ~ "Psoriatic arthritis",
    MCQ160N == "1" ~ "Gout",
    MCQ190 == "1" ~ "Rheumatoid arthritis",
    MCQ190 == "2" ~ "Osteoarthritis",
    MCQ191 == "2" ~ "Osteoarthritis",
    MCQ191 == "1" ~ "Rheumatoid arthritis",
    MCQ191 == "3" ~ "Psoriatic arthritis",
  )) 

# creating adequate sample weights
# herein we use the 24h recall day 1 weight, and divide by the number of cycles (6)

rh_dxa_all <- rh_dxa_all |> 
  mutate(weights_all = case_when(
    DRDINT.x == 2 ~ (WTDR2D.x * (1/7)),
    DRDINT.x == 1 ~ (WTDRD1.x * (1/7))), 
    weights_all_exam = WTMEC2YR * (1/7))

# calculating physical activity
# obtain METs for each domain and intensity, then sum it up, then classify

rh_dxa_all <- rh_dxa_all |> 
  # cleaning variables first
  mutate(PAQ605 = case_when(PAQ605 %in% c(7, 9) ~ NA_real_, # these are codes for non-responses within the questionnaires, important to remove them
                            TRUE ~ PAQ605),
         PAQ610 = case_when(PAQ610 %in% c(77, 99) ~ NA_real_,
                            TRUE ~ PAQ610),
         PAD615 = case_when(PAD615 %in% c(7777, 9999) ~ NA_real_,
                            TRUE ~ PAD615),
         
         PAQ620 = case_when(PAQ620 %in% c(7, 9) ~ NA_real_,
                            TRUE ~ PAQ620),
         PAQ625 = case_when(PAQ625 %in% c(77, 99) ~ NA_real_,
                            TRUE ~ PAQ625),
         PAD630 = case_when(PAD630 %in% c(7777, 9999) ~ NA_real_,
                            TRUE ~ PAD630),
         
         PAQ635 = case_when(PAQ635 %in% c(7, 9) ~ NA_real_,
                            TRUE ~ PAQ635),
         PAQ640 = case_when(PAQ640 %in% c(77, 99) ~ NA_real_,
                            TRUE ~ PAQ640),
         PAD645 = case_when(PAD645 %in% c(7777, 9999) ~ NA_real_,
                            TRUE ~ PAD645),
         
         PAQ650 = case_when(PAQ650 %in% c(7, 9) ~ NA_real_,
                            TRUE ~ PAQ650),
         PAQ655 = case_when(PAQ655 %in% c(77, 99) ~ NA_real_,
                            TRUE ~ PAQ655),
         PAD660 = case_when(PAD660 %in% c(7777, 9999) ~ NA_real_,
                            TRUE ~ PAD660),
         
         PAQ665 = case_when(PAQ665 %in% c(7, 9) ~ NA_real_,
                            TRUE ~ PAQ665),
         PAQ670 = case_when(PAQ670 %in% c(77, 99) ~ NA_real_,
                            TRUE ~ PAQ670),
         PAD675 = case_when(PAD675 %in% c(7777, 9999) ~ NA_real_,
                            TRUE ~ PAD675),
         # then calculate METs per domain and overall
         work_vig_MET = if_else(PAQ605 == 2, 0, PAQ610 * PAD615), # calculating MET-minutes based on GPAQ questionnaire scoring
         work_mod_MET = if_else(PAQ620 == 2, 0, PAQ625 * PAD630),
         transport_MET = if_else(PAQ635 == 2, 0, PAQ640 * PAD645),
         rec_vig_MET = if_else(PAQ650 == 2, 0, PAQ655 * PAD660),
         rec_mod_MET = if_else(PAQ665 == 2, 0, PAQ670 * PAD675),
         work_MET = (work_vig_MET * 8) + (work_mod_MET * 4),
         rec_MET = (rec_vig_MET * 8) + (rec_mod_MET * 4),
         MET = work_MET + rec_MET + transport_MET, # summing it up
         # classification based on METs
         MET_class = case_when(MET < 600 ~ "Very low PA",
                               MET >= 600 & MET < 1200 ~ "Moderate PA",
                               MET >= 1200 & MET < 1800 ~ "High PA",
                               MET >= 1800 ~ "Very high PA")
  )

# medication use, needs to be processed in separate due to different dataframe structure
# calculating gc use and days of gc use

# load previously joined medication dataframe

load("data/medication_df.Rdata")

# create a string with the names of glucocorticoids 
gluc_names <- c('^(DEFLAZACORT|PREDNISOLONE|DEXAMETHASONE|HYDROCORTISONE|METHYLPREDNISOLONE|PREDNISONE|CORTISONE)$')

med_df <- med_df |> 
  mutate(uses_gc = if_else(stringr::str_detect(RXDDRUG, gluc_names), 1, 0), # detect these names
         days_gc = case_when(uses_gc == 0 ~ 0,
                             RXDDAYS %in% c(77777, 99999) ~ 77777,
                             TRUE ~ RXDDAYS * uses_gc)) |> 
  select(SEQN, uses_gc, days_gc, RXDDRUG) |> 
  filter(stringr::str_detect(RXDDRUG, gluc_names)) |> 
  group_by(SEQN) |> # sometimes some participants would report using more than one GC type,
  summarise(days_gc = mean(days_gc)) |> # so we group by participant and take the mean of both medications to have one single number
  ungroup()

# merge into main and fix

rh_dxa_all <- rh_dxa_all |> full_join(med_df) |> 
  mutate(days_gc = if_else(is.na(days_gc), 0, days_gc), # if days_gc is NA, it means the person does not take GC, i.e. a vlaue of 0
         days_gc = if_else(days_gc == 77777, NA_real_, days_gc)) # non response code

# calculating PAL from METs, EEE (IOM) and class of reporting
# this calculating of reporting status follows the methodology recommended by Jessri et al: 10.1017/S0007114515004237
# and described in detail in McCrory et al: 10.1079/PHN2002387

# sd value for cutoff points 

cv_t_mtee = 0.082 # values from McCrory et al based on DLW data
cv_w_ptee = 0.177

main_sample <- rh_dxa_all |> 
  filter(!is.na(disease_cat) & DR1DRSTZ == "1" & RIDAGEYR >= 18 & DRDINT.x == 2)

all_measurements <- c(main_sample$DR1TKCAL, main_sample$DR2TKCAL)
overall_mean <- mean(all_measurements)
overall_sd <- sd(all_measurements)
cv_w_ei <- (overall_sd / overall_mean) * 100 # calculating coefficient of variance of energy in our sample

sd_cutoff = sqrt(
  cv_w_ei^2/2 + cv_t_mtee^2 + cv_w_ptee^2 # creating a cutoff based on our distribution
)

lower_cutoff = 100 - sd_cutoff
upper_cutoff = 100 + sd_cutoff # creating cutoffs

rh_dxa_all <- rh_dxa_all |> 
  mutate(BEE = if_else(RIAGENDR == 1, # calculating basal energy based on IOM formula
                       293 + (456.4 * (BMXHT/100)) + (10.12 * BMXWT) - (3.8 * RIDAGEYR),
                       247 + (401.5 * (BMXHT/100)) + (8.6 * BMXWT) - (2.67 * RIDAGEYR)),
         delta_PAL_total = ( # converting MET into PAL according to Gerrior et al: PMID 16978504
           ((8 - 1) * ((1.15/0.9) * work_vig_MET/7)/1440)/(BEE/(0.0175 * 1440 * BMXWT)) +
             ((4 - 1) * ((1.15/0.9) * work_mod_MET/7)/1440)/(BEE/(0.0175 * 1440 * BMXWT)) +
             ((4 - 1) * ((1.15/0.9) * transport_MET/7)/1440)/(BEE/(0.0175 * 1440 * BMXWT)) +
             ((8 - 1) * ((1.15/0.9) * rec_vig_MET/7)/1440)/(BEE/(0.0175 * 1440 * BMXWT)) +
             ((4 - 1) * ((1.15/0.9) * rec_mod_MET/7)/1440)/(BEE/(0.0175 * 1440 * BMXWT))
         ),
         PAL = 1.1 + delta_PAL_total,
         PA = if_else(RIAGENDR == 1,
                      case_when(PAL < 1.4 ~ 1,
                                PAL >= 1.4 & PAL < 1.6 ~ 1.12,
                                PAL >= 1.6 & PAL < 1.9 ~ 1.27,
                                PAL >= 1.9 ~ 1.54),
                      
                      case_when(PAL < 1.4 ~ 1,
                                PAL >= 1.4 & PAL < 1.6 ~ 1.14,
                                PAL >= 1.6 & PAL < 1.9 ~ 1.27,
                                PAL >= 1.9 ~ 1.45)
         ),
         EER = if_else(RIAGENDR == 1, # calculating estimated energy requirements
                       662 - (9.53 * RIDAGEYR) + (PA * ((15.91 * BMXWT) + (539.6 * (BMXHT/100)))),
                       354 - (6.91 * RIDAGEYR) + (PA * ((9.36 * BMXWT) + (726 * (BMXHT/100))))
         ),
         EI_EER_ratio = (kcal/EER)*100, # ratio between what was consumed and requirements
         EI_EER_class = case_when(EI_EER_ratio < lower_cutoff ~ "Under-reporter", # classification based on the cutoff points
                                  EI_EER_ratio >= lower_cutoff & EI_EER_ratio < upper_cutoff ~ "Plausible reporter",
                                  EI_EER_ratio >= upper_cutoff ~ "Over-reporter")
  )

# setting up survey design

rh_dxa_all_analyse <- subset(rh_dxa_all, !is.na(SDMVPSU) & !is.na(SDMVSTRA)) |> 
  dplyr::select(SEQN, PTN, CHO, LIP, kcal, PTNgkg, DXDTOLE, DXDTOBMD, DXXOFBMD, 
                DXXOSBMD, RIAGENDR, BMXHT, BMXWT, RIDAGEYR, disease_cat, 
                SDMVPSU, SDMVSTRA, weights_all, RIDRETH1, BMXBMI, femur_t_score, 
                spine_t_score, DR1DRSTZ, DXAEXSTS, DXAFMRST, DXASPNST, MET, MET_class,
                days_gc, EI_EER_class, cycle)

# this dataset can then be saved and loaded when you actually want to analyze it.
# save(rh_dxa_all_analyse, file=here("data/NHANES_RH.Rdata"))

# set survey
# set individuals without sample weights to 0, to exclude later
rh_dxa_all_analyse <- rh_dxa_all_analyse |> 
  mutate(weights_all = if_else(is.na(weights_all), 0, weights_all))

NHANES_all <- svydesign(data=rh_dxa_all_analyse, 
                        id=~SDMVPSU, 
                        strata=~SDMVSTRA, 
                        weights=~weights_all,
                        nest=TRUE)

# subsetting
# quick codebook for subset: MCQ195 = 1, 2, 3 (osteoarth, RA, psoriatic A patients); 
#                            MCQ190, MCQ191 = 1 (RA patients, different documentation);
#                            DXAEXSTS; DXASPNST; DXAFMRST  = 1 (adequate DXA exam);  

NHANES_rh <- subset(NHANES_all, weights_all != 0) # excluding individuals without sample weight

NHANES_rh <- subset(NHANES_all, !is.na(disease_cat)) # selecting RDs 

NHANES_rh <- subset(NHANES_rh, DR1DRSTZ == "1") # selecting adequate r24 data

NHANES_rh <- subset(NHANES_rh, RIDAGEYR >= 18 & RIDAGEYR < 98) # selecting adults 

NHANES_subset_lm <- subset(NHANES_rh, DXAEXSTS == "1" & !is.na(DXDTOLE))

NHANES_subset_wb <- subset(NHANES_rh, DXAEXSTS == "1" & !is.na(DXDTOBMD))

NHANES_subset_femur <- subset(NHANES_rh, DXAFMRST == "1" & !is.na(DXXOFBMD))

NHANES_subset_spine <- subset(NHANES_rh, DXASPNST == "1" & !is.na(DXXOSBMD))

# separate dfs for each outcome, for table 1

lm_dat <- as_tibble(NHANES_subset_lm$variables)

wb_dat <- as_tibble(NHANES_subset_wb$variables)

fem_dat <- as_tibble(NHANES_subset_femur$variables)

spn_dat <- as_tibble(NHANES_subset_spine$variables)

sample_dat <- reduce(list(lm_dat, wb_dat, fem_dat, spn_dat), bind_rows) |> distinct()

nrow(sample_dat) # sample size

# creating subsets per outcome
NHANES_subset_lm <- subset(NHANES_rh, DXAEXSTS == "1" & !is.na(DXDTOLE) # make sure subsets
                           & !is.na(PTN) & !is.na(CHO) & !is.na(LIP) # have available data on all covariates
                           & !is.na(RIAGENDR) & !is.na(BMXWT) & !is.na(RIDAGEYR)
                           & !is.na(MET) & !is.na(days_gc) & !is.na(EI_EER_class))

nrow(NHANES_subset_lm$variables) # subset sample size

NHANES_subset_wb <- subset(NHANES_rh, DXAEXSTS == "1" & !is.na(DXDTOBMD) # and so on..
                           & !is.na(PTN) & !is.na(CHO) & !is.na(LIP)
                           & !is.na(RIAGENDR) & !is.na(BMXWT) & !is.na(RIDAGEYR)
                           & !is.na(MET) & !is.na(days_gc) & !is.na(EI_EER_class))

nrow(NHANES_subset_wb$variables)

NHANES_subset_femur <- subset(NHANES_rh, DXAFMRST == "1" & !is.na(DXXOFBMD)
                              & !is.na(PTN) & !is.na(CHO) & !is.na(LIP)
                              & !is.na(RIAGENDR) & !is.na(BMXWT) & !is.na(RIDAGEYR)
                              & !is.na(MET) & !is.na(days_gc) & !is.na(EI_EER_class))

nrow(NHANES_subset_femur$variables)

NHANES_subset_spine <- subset(NHANES_rh, DXASPNST == "1" & !is.na(DXXOSBMD)
                              & !is.na(PTN) & !is.na(CHO) & !is.na(LIP)
                              & !is.na(RIAGENDR) & !is.na(BMXWT) & !is.na(RIDAGEYR)
                              & !is.na(MET) & !is.na(days_gc) & !is.na(EI_EER_class))

nrow(NHANES_subset_spine$variables)

## Modelling ###################################################################

# unadjusted models

# Lean mass

u_model_leanmass_protein <- svyglm(DXDTOLE ~ PTN, # protein in grams
                                 design=NHANES_subset_lm) # model fit
summary(u_model_leanmass_protein) # summary
u_model1<-broom::tidy(u_model_leanmass_protein, conf.int=T) # broom::tidy extracts model coefficients and CIs. Doing this for all models.

u_model_leanmass_proteingkg <- svyglm(DXDTOLE ~ PTNgkg, # relative protein g/kg 
                                    design=NHANES_subset_lm) 
summary(u_model_leanmass_proteingkg) 
u_model2<-broom::tidy(u_model_leanmass_proteingkg,conf.int=T)

# BMD

# whole-body

u_model_wbBMD_protein <- svyglm(DXDTOBMD ~ PTN, 
                              design=NHANES_subset_wb) 
summary(u_model_wbBMD_protein) 
u_model3<-broom::tidy(u_model_wbBMD_protein,conf.int=T)

u_model_wbBMD_proteingkg <- svyglm(DXDTOBMD~PTNgkg, 
                                 design=NHANES_subset_wb) 
summary(u_model_wbBMD_proteingkg) 
u_model4<-broom::tidy(u_model_wbBMD_proteingkg,conf.int=T)

# femur

u_model_femurBMD_protein <- svyglm(DXXOFBMD ~ PTN, design=NHANES_subset_femur) 
summary(u_model_femurBMD_protein) 
u_model5<-broom::tidy(u_model_femurBMD_protein,conf.int=T)

u_model_femurBMD_proteingkg <- svyglm(DXXOFBMD ~ PTNgkg, design=NHANES_subset_femur) 
summary(u_model_femurBMD_proteingkg) 
u_model6<-broom::tidy(u_model_femurBMD_proteingkg,conf.int=T)

# spine

u_model_spineBMD_protein <- svyglm(DXXOSBMD ~ PTN, design=NHANES_subset_spine) 
summary(u_model_spineBMD_protein)
u_model7<-broom::tidy(u_model_spineBMD_protein,conf.int=T)

u_model_spineBMD_proteingkg <- svyglm(DXXOSBMD ~ PTNgkg, 
                                    design=NHANES_subset_spine) 
summary(u_model_spineBMD_proteingkg)
u_model8<-broom::tidy(u_model_spineBMD_proteingkg,conf.int=T)

# creating a dataframe with all results and tidying 
# for later use in plots and summaries

table2_un <- reduce(list(u_model1,u_model2,u_model3,u_model4,
            u_model5,u_model6,u_model7,u_model8),
                            full_join) |> 
  filter(term=="PTN" | term=="PTNgkg") |> 
  mutate(type=c("lean mass","lean mass","wb BMD","wb BMD","femur BMD","femur BMD",
                "spine BMD","spine BMD")) |> 
  mutate(across(c(estimate, conf.low, conf.high),
                .fns = ~ if_else(term == "PTN", .x * 50, .x))) |> 
  mutate(value=paste(estimate," (",conf.low,"; ",conf.high,")",sep="")) |> 
  mutate(estimate=if_else(type == "lean mass", round(estimate, 0), round(estimate, 4))) |> 
  mutate(conf.low=if_else(type == "lean mass", round(conf.low), round(conf.low, 4))) |> 
  mutate(conf.high=if_else(type == "lean mass", round(conf.high, 0), round(conf.high, 4))) |> 
  mutate(value2=paste(conf.low,"; ",conf.high,sep="")) |> 
  mutate(p.value = round(p.value, 4)) |> 
  select(term, estimate, value2, p.value, type) 

# write_xlsx(table2_un, here("Tables/nhanes_model_summary_unadjusted.xlsx"))

# adjusted models

# Lean mass

model_leanmass_protein <- svyglm(DXDTOLE ~ PTN + CHO + LIP + RIAGENDR + BMXWT # same process
                                 + RIDAGEYR + MET + days_gc + EI_EER_class, # but now includes covariates determined by DAG
                              design=NHANES_subset_lm) 

summary(model_leanmass_protein) 
model1<-broom::tidy(model_leanmass_protein, conf.int=T)

model_leanmass_proteingkg <- svyglm(DXDTOLE ~ PTNgkg + CHO + LIP + RIAGENDR + 
                                      BMXWT + RIDAGEYR + MET + days_gc + EI_EER_class, 
                                 design=NHANES_subset_lm) 
summary(model_leanmass_proteingkg) 
model2<-broom::tidy(model_leanmass_proteingkg,conf.int=T)

# BMD

# whole-body

model_wbBMD_protein <- svyglm(DXDTOBMD ~ PTN + CHO + LIP + RIAGENDR + BMXWT + RIDAGEYR
                              + MET + days_gc + EI_EER_class, 
                              design=NHANES_subset_wb) 
summary(model_wbBMD_protein) # total body BMD and total protein
model3<-broom::tidy(model_wbBMD_protein,conf.int=T)

model_wbBMD_proteingkg <- svyglm(DXDTOBMD~PTNgkg  + CHO + LIP + RIAGENDR + BMXWT + 
                                   RIDAGEYR + MET + days_gc + EI_EER_class, 
                                 design=NHANES_subset_wb) 
summary(model_wbBMD_proteingkg) # total body BMD and protein g/kg 
model4<-broom::tidy(model_wbBMD_proteingkg,conf.int=T)

# femur

model_femurBMD_protein <- svyglm(DXXOFBMD ~ PTN + CHO + LIP + 
                                   RIAGENDR + BMXWT + RIDAGEYR
                                 + MET + days_gc + EI_EER_class, design=NHANES_subset_femur) 
summary(model_femurBMD_protein) # femur BMD and total protein
model5<-broom::tidy(model_femurBMD_protein,conf.int=T)

model_femurBMD_proteingkg <- svyglm(DXXOFBMD ~ PTNgkg + CHO + LIP +
                                      RIAGENDR + BMXWT + RIDAGEYR
                                    + MET + days_gc + EI_EER_class, design=NHANES_subset_femur) 
summary(model_femurBMD_proteingkg) # femur BMD and protein g/kg
model6<-broom::tidy(model_femurBMD_proteingkg,conf.int=T)

# spine

model_spineBMD_protein <- svyglm(DXXOSBMD ~ PTN + CHO + LIP + 
                                   RIAGENDR + BMXWT + RIDAGEYR
                                 + MET + days_gc + EI_EER_class, design=NHANES_subset_spine) 
summary(model_spineBMD_protein) # spine BMD and total protein
model7<-broom::tidy(model_spineBMD_protein,conf.int=T)

model_spineBMD_proteingkg <- svyglm(DXXOSBMD ~ PTNgkg + CHO + LIP + 
                                      RIAGENDR + BMXWT + RIDAGEYR
                                    + MET + days_gc + EI_EER_class, 
                                    design=NHANES_subset_spine) 
summary(model_spineBMD_proteingkg) # spine BMD and protein g/kg
model8<-broom::tidy(model_spineBMD_proteingkg,conf.int=T)

# building a dataframe with results

model_summary_all <- reduce(list(model1,model2,model3,model4,
                                 model5,model6,model7,model8),
                            full_join) |> 
  filter(term=="PTN" | term=="PTNgkg") |> 
  mutate(type=c("lean mass","lean mass","wb BMD","wb BMD","femur BMD","femur BMD",
                 "spine BMD","spine BMD")) |> 
  mutate(across(c(estimate, conf.low, conf.high),
                  .fns = ~ if_else(term == "PTN", .x * 50, .x))) |> 
  mutate(value=paste(estimate," (",conf.low,"; ",conf.high,")",sep="")) |> 
  mutate(estimate=round(estimate, 4)) |> 
  mutate(conf.low=round(conf.low, 4)) |> 
  mutate(conf.high=round(conf.high, 4)) |> 
  mutate(value=paste(conf.low,"; ",conf.high,sep="")) |> 
  select(term, estimate, value, p.value, type, conf.low, conf.high)

# code to save the table in Excel format. Commented out to avoid saving multiple files.
# write_xlsx(model_summary_all, here("Tables/nhanes_model_summary.xlsx"))

# standardized unadjusted coefficients

# Lean mass

model_leanmass_protein_std <- svyglm(scale(DXDTOLE) ~ scale(PTN), # using scale() on predictor and outcome variables 
                                     design=NHANES_subset_lm) # to produce standardized coefficients
summary(model_leanmass_protein_std) 
model1<-broom::tidy(model_leanmass_protein_std, conf.int=T)

model_leanmass_proteingkg_std <- svyglm(scale(DXDTOLE) ~ scale(PTNgkg), 
                                        design=NHANES_subset_lm) 
summary(model_leanmass_proteingkg_std) 
model2<-broom::tidy(model_leanmass_proteingkg_std, conf.int=T)

# BMD

# whole-body

model_wbBMD_protein_std <- svyglm(scale(DXDTOBMD) ~ scale(PTN),
                                  design=NHANES_subset_wb) 
summary(model_wbBMD_protein_std) 
model3<-broom::tidy(model_wbBMD_protein_std, conf.int=T)

model_wbBMD_proteingkg_std <- svyglm(scale(DXDTOBMD) ~ scale(PTNgkg),
                                     design=NHANES_subset_wb) 
summary(model_wbBMD_proteingkg_std) # total body BMD and protein g/kg 
model4<-broom::tidy(model_wbBMD_proteingkg_std ,conf.int=T)

# femur

model_femurBMD_protein_std <- svyglm(scale(DXXOFBMD) ~ scale(PTN), 
                                     design=NHANES_subset_femur) 
summary(model_femurBMD_protein_std) # femur BMD and total protein
model5<-broom::tidy(model_femurBMD_protein_std, conf.int=T)

model_femurBMD_proteingkg_std <- svyglm(scale(DXXOFBMD) ~ scale(PTNgkg), 
                                        design=NHANES_subset_femur) 
summary(model_femurBMD_proteingkg_std) # femur BMD and protein g/kg
model6<-broom::tidy(model_femurBMD_proteingkg_std, conf.int=T)

# spine

model_spineBMD_protein_std <- svyglm(scale(DXXOSBMD) ~ scale(PTN), 
                                     design=NHANES_subset_spine) 
summary(model_spineBMD_protein_std) # spine BMD and total protein
model7<-broom::tidy(model_spineBMD_protein_std, conf.int=T)

model_spineBMD_proteingkg_std <- svyglm(scale(DXXOSBMD) ~ scale(PTNgkg), 
                                        design=NHANES_subset_spine) 
summary(model_spineBMD_proteingkg_std) # spine BMD and protein g/kg
model8<-broom::tidy(model_spineBMD_proteingkg_std, conf.int=T)

# table summary

model_summary_all_std_unadj <- reduce(list(model1,model2,model3,model4,
                                     model5,model6,model7,model8),
                                full_join) |> 
  filter(term=="scale(PTN)" | term=="scale(PTNgkg)") |> 
  mutate(type=c("lean mass","lean mass","wb BMD","wb BMD","femur BMD","femur BMD",
                "spine BMD","spine BMD")) |> 
  mutate(estimate=round(estimate, 2)) |> 
  mutate(conf.low=round(conf.low, 2)) |> 
  mutate(conf.high=round(conf.high, 2)) |> 
  mutate(value=paste(estimate," (",conf.low,"; ",conf.high,")",sep="")) |>  # creating different formatting for later use in tables and results reporting
  mutate(value2=paste(conf.low,"; ",conf.high,sep="")) |> 
  mutate(p.value = round(p.value, 4))

# write_xlsx(model_summary_all_std_unadj, here("Tables/nhanes_model_summary_std_1.xlsx"))

# standardized adjusted coefficients

# Lean mass

model_leanmass_protein_std <- svyglm(scale(DXDTOLE) ~ scale(PTN) + CHO + LIP + 
                                       disease_cat + RIAGENDR + BMXWT + RIDAGEYR 
                                       + MET + days_gc + EI_EER_class, 
                                 design=NHANES_subset_lm) 
summary(model_leanmass_protein_std) 
model1<-broom::tidy(model_leanmass_protein_std, conf.int=T)

model_leanmass_proteingkg_std <- svyglm(scale(DXDTOLE) ~ scale(PTNgkg) + CHO + LIP + disease_cat + RIAGENDR + 
                                      BMXWT + RIDAGEYR + MET + days_gc + EI_EER_class, 
                                    design=NHANES_subset_lm) 
summary(model_leanmass_proteingkg_std) 
model2<-broom::tidy(model_leanmass_proteingkg_std, conf.int=T)

# BMD

model_wbBMD_protein_std <- svyglm(scale(DXDTOBMD) ~ scale(PTN) + CHO + LIP + disease_cat + 
                                    RIAGENDR + BMXWT + RIDAGEYR+ MET + days_gc + EI_EER_class, 
                              design=NHANES_subset_wb) 
summary(model_wbBMD_protein_std) # total body BMD and total protein
model3<-broom::tidy(model_wbBMD_protein_std, conf.int=T)

model_wbBMD_proteingkg_std <- svyglm(scale(DXDTOBMD) ~ scale(PTNgkg)  + CHO + LIP + disease_cat+ RIAGENDR + BMXWT + RIDAGEYR
                                     + MET + days_gc + EI_EER_class, 
                                 design=NHANES_subset_wb) 
summary(model_wbBMD_proteingkg_std) # total body BMD and protein g/kg 
model4<-broom::tidy(model_wbBMD_proteingkg_std ,conf.int=T)

# femur

model_femurBMD_protein_std <- svyglm(scale(DXXOFBMD) ~ scale(PTN) + CHO + LIP + disease_cat+ 
                                   RIAGENDR + BMXWT + RIDAGEYR
                                   + MET + days_gc + EI_EER_class, design=NHANES_subset_femur) 
summary(model_femurBMD_protein_std) # femur BMD and total protein
model5<-broom::tidy(model_femurBMD_protein_std, conf.int=T)

model_femurBMD_proteingkg_std <- svyglm(scale(DXXOFBMD) ~ scale(PTNgkg) + CHO + LIP + disease_cat+
                                      RIAGENDR + BMXWT + RIDAGEYR+ MET + days_gc + EI_EER_class, design=NHANES_subset_femur) 
summary(model_femurBMD_proteingkg_std) # femur BMD and protein g/kg
model6<-broom::tidy(model_femurBMD_proteingkg_std, conf.int=T)

# spine

model_spineBMD_protein_std <- svyglm(scale(DXXOSBMD) ~ scale(PTN) + CHO + LIP+ disease_cat + 
                                   RIAGENDR + BMXWT + RIDAGEYR+ MET + days_gc + EI_EER_class, design=NHANES_subset_spine) 
summary(model_spineBMD_protein_std) # spine BMD and total protein
model7<-broom::tidy(model_spineBMD_protein_std, conf.int=T)
  
model_spineBMD_proteingkg_std <- svyglm(scale(DXXOSBMD) ~ scale(PTNgkg) + CHO + LIP+ disease_cat + 
                                      RIAGENDR + BMXWT + RIDAGEYR+ MET + days_gc + EI_EER_class, 
                                    design=NHANES_subset_spine) 
summary(model_spineBMD_proteingkg_std) # spine BMD and protein g/kg
model8<-broom::tidy(model_spineBMD_proteingkg_std, conf.int=T)

# table summary

model_summary_all_std <- reduce(list(model1,model2,model3,model4,
                                 model5,model6,model7,model8),
                            full_join) |> 
  filter(term=="scale(PTN)" | term=="scale(PTNgkg)") |> 
  mutate(type=c("lean mass","lean mass","wb BMD","wb BMD","femur BMD","femur BMD",
                "spine BMD","spine BMD")) |> 
  mutate(estimate=round(estimate, 2)) |> 
  mutate(conf.low=round(conf.low, 2)) |> 
  mutate(conf.high=round(conf.high, 2)) |> 
  mutate(value2=paste(conf.low,"; ",conf.high,sep="")) |> 
  mutate(value=paste(estimate," (",conf.low,"; ",conf.high,")",sep="")) |> 
  mutate(p.value = round(p.value, 4)) 

# write_xlsx(model_summary_all_std, here("Tables/nhanes_model_summary_std_2.xlsx"))

#### Tables ####################################################################

#creating table 1 with gtsummary

table_1_summary <- sample_dat |>
 select(disease_cat, RIDAGEYR, RIAGENDR, RIDRETH1, BMXHT, BMXWT, BMXBMI, days_gc, MET, MET_class, kcal, CHO, PTN, PTNgkg, LIP,
        DXDTOLE, DXDTOBMD, DXXOFBMD, femur_t_score, DXXOSBMD, spine_t_score) |>
 tbl_summary(digits = list(all_categorical() ~ 0,
                           all_continuous() ~ 3), by = disease_cat, missing = "no", statistic = all_continuous() ~ "{mean} ({sd})") |>
 add_overall()

#save
#table_1_summary |> as_gt() |> gt::gtsave(filename = "Tables/table1_gt_07_03_24.rtf")

#### Plots #####################################################################

# tidying result dataframes to be used later

model_summary_all_std_unadj <- model_summary_all_std_unadj |> 
  mutate(term = as.factor(c("Lean mass ~ protein (g/d)",
                            "Lean mass ~ protein (g/kg)",
                            "Whole-body BMD ~ protein (g/d)",
                            "Whole-body BMD ~ protein (g/kg)",
                            "Femur BMD ~ protein (g/d)",
                            "Femur BMD ~ protein (g/kg)",
                            "Spine BMD ~ protein (g/d)",
                            "Spine BMD ~ protein (g/kg)")),
         adj = "Unadjusted coefficients") 

model_summary_all_std_unadj2 <- model_summary_all_std_unadj |> 
  mutate(gkg = rep(c("g", "g/kg"), 4)) |> 
  filter(gkg == "g")
  
model_summary_all_std <- model_summary_all_std |> 
  mutate(term = as.factor(c("Lean mass ~ protein (g/d)",
                  "Lean mass ~ protein (g/kg)",
                  "Whole-body BMD ~ protein (g/d)",
                  "Whole-body BMD ~ protein (g/kg)",
                  "Femur BMD ~ protein (g/d)",
                  "Femur BMD ~ protein (g/kg)",
                  "Spine BMD ~ protein (g/d)",
                  "Spine BMD ~ protein (g/kg)")),
         adj = "Adjusted coefficients")

model_std_coefs <- full_join(model_summary_all_std_unadj2, 
                             model_summary_all_std) |> 
  mutate(adj = fct_relevel(as.factor(adj), "Unadjusted coefficients",
                                 "Adjusted coefficients")) |> 
  mutate(term = fct_relevel(term, "Lean mass ~ protein (g/d)",
                            "Lean mass ~ protein (g/kg)",
                            "Whole-body BMD ~ protein (g/d)",
                            "Whole-body BMD ~ protein (g/kg)",
                            "Femur BMD ~ protein (g/d)",
                            "Femur BMD ~ protein (g/kg)",
                            "Spine BMD ~ protein (g/d)",
                            "Spine BMD ~ protein (g/kg)")) |> 
  mutate(term = fct_rev(term))

## prediction plots

# creating result reporting texts for plots

txt1 <- model_std_coefs |>
  filter(term=="Lean mass ~ protein (g/d)" & adj=="Adjusted coefficients") |> 
  select(estimate, p.value, conf.low, conf.high)  

txt2 <- model_std_coefs |>
  filter(term=="Lean mass ~ protein (g/kg)" & adj=="Adjusted coefficients") |> 
  select(estimate, p.value, conf.low, conf.high)  

txt3 <- model_std_coefs |>
  filter(term=="Whole-body BMD ~ protein (g/d)" & adj=="Adjusted coefficients") |> 
  select(estimate, p.value, conf.low, conf.high)  

txt4 <- model_std_coefs |>
  filter(term=="Whole-body BMD ~ protein (g/kg)" & adj=="Adjusted coefficients") |> 
  select(estimate, p.value, conf.low, conf.high)  

txt5 <- model_std_coefs |>
  filter(term=="Femur BMD ~ protein (g/d)" & adj=="Adjusted coefficients") |> 
  select(estimate, p.value, conf.low, conf.high)  

txt6 <- model_std_coefs |>
  filter(term=="Femur BMD ~ protein (g/kg)" & adj=="Adjusted coefficients") |> 
  select(estimate, p.value, conf.low, conf.high)  

txt7 <- model_std_coefs |>
  filter(term=="Spine BMD ~ protein (g/d)" & adj=="Adjusted coefficients") |> 
  select(estimate, p.value, conf.low, conf.high)  

txt8 <- model_std_coefs |>
  filter(term=="Spine BMD ~ protein (g/kg)" & adj=="Adjusted coefficients") |> 
  select(estimate, p.value, conf.low, conf.high)  

summa <- model_summary_all |> 
  mutate(estimate = if_else(type == "lean mass", round(estimate, 0), round(estimate, 4))) |> 
  mutate(conf.low = if_else(type == "lean mass", round(conf.low, 0), round(conf.low, 4))) |> 
  mutate(conf.high = if_else(type == "lean mass", round(conf.high, 0), round(conf.high, 4))) |> 
  mutate(p.value = if_else(p.value >= 0.05, round(p.value, 2), round(p.value, 5)))

txt1_an <- glue("Coefficient: {summa$estimate[1]} g lean mass per 50 g protein (95%CI {summa$conf.low[1]} to {summa$conf.high[1]})\nβ: {txt1$estimate} (95%CI {txt1$conf.low} to {txt1$conf.high}), p-value={txt1$p.value}")
txt2_an <- glue("Coefficient: {summa$estimate[2]} g lean mass per 1 g/kgBM protein (95%CI {summa$conf.low[2]} to {summa$conf.high[2]})\nβ: {txt2$estimate} (95%CI {txt2$conf.low} to {txt2$conf.high}), p-value={txt2$p.value}")
txt3_an <- glue("Coefficient: {summa$estimate[3]} gm/cm² BMD per 50 g protein (95%CI {summa$conf.low[3]} to {summa$conf.high[3]})\nβ: {txt3$estimate} (95%CI {txt3$conf.low} to {txt3$conf.high}), p-value={txt3$p.value}")
txt4_an <- glue("Coefficient: {summa$estimate[4]} gm/cm² BMD per 1 g/kgBM protein (95%CI {summa$conf.low[4]} to {summa$conf.high[4]})\nβ: {txt4$estimate} (95%CI {txt4$conf.low} to {txt4$conf.high}), p-value={txt4$p.value}")
txt5_an <- glue("Coefficient: {summa$estimate[5]} gm/cm² BMD per 50 g protein (95%CI {summa$conf.low[5]} to {summa$conf.high[5]})\nβ: {txt5$estimate} (95%CI {txt5$conf.low} to {txt5$conf.high}), p-value={txt5$p.value}")
txt6_an <- glue("Coefficient: {summa$estimate[6]} gm/cm² BMD per 1 g/kgBM protein (95%CI {summa$conf.low[6]} to {summa$conf.high[6]})\nβ: {txt6$estimate} (95%CI {txt6$conf.low} to {txt6$conf.high}), p-value={txt6$p.value}")
txt7_an <- glue("Coefficient: {summa$estimate[7]} gm/cm² BMD per 50 g protein (95%CI {summa$conf.low[7]} to {summa$conf.high[7]})\nβ: {txt7$estimate} (95%CI {txt7$conf.low} to {txt7$conf.high}), p-value={txt7$p.value}")
txt8_an <- glue("Coefficient: {summa$estimate[8]} gm/cm² BMD per 1 g/kgBM protein (95%CI {summa$conf.low[8]} to {summa$conf.high[8]})\nβ: {txt8$estimate} (95%CI {txt8$conf.low} to {txt8$conf.high}), p-value={txt8$p.value}")

## prediction plots with marginaleffects

pred1 <- predictions(model_leanmass_protein, # predicts the linear relationship between protein and outcome using the model
                   by=c("PTN"), newdata="mean") # averages using the mean of covariates

pred2 <- predictions(model_leanmass_proteingkg,
                   by=c("PTNgkg"), newdata="mean")

pred3 <- predictions(model_wbBMD_protein,
                   by=c("PTN"), newdata="mean")

pred4 <- predictions(model_wbBMD_proteingkg,
                   by=c("PTNgkg"), newdata="mean")

pred5 <- predictions(model_femurBMD_protein,
                   by=c("PTN"), newdata="mean")

pred6 <- predictions(model_femurBMD_proteingkg,
                   by=c("PTNgkg"), newdata="mean")

pred7 <- predictions(model_spineBMD_protein,
                   by=c("PTN"), newdata="mean")

pred8 <- predictions(model_spineBMD_proteingkg,
                   by=c("PTNgkg"), newdata="mean")

# creating the plots

m_p1 <- ggplot(pred1, aes(x=PTN, y=estimate)) +
  geom_point(data=NHANES_subset_lm$variables, aes(x=PTN, y=DXDTOLE, size=weights_all), alpha=.1) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  #geom_rug(alpha=.1, sides="b") +
  labs(x="Protein intake (g)", y="Lean mass (g)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt1_an, hjust=1, vjust=-0.2,
           size=2.5) +
  scale_x_continuous(limits=c(25, 250)) +
  scale_y_continuous(limits=c(20000, 90000)) +
  guides(size='none') +
  theme_avp()

m_p2 <- ggplot(pred2, aes(x=PTNgkg, y=estimate)) +
  geom_point(data=NHANES_subset_lm$variables, aes(x=PTNgkg, y=DXDTOLE, size=weights_all), alpha=.1) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  #geom_rug(alpha=.1, sides="b") +
  labs(x="Protein intake (g/kgBM)", y="Lean mass (g)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt2_an, hjust=1, vjust=-0.2,
           size=2.5) +
  scale_x_continuous(limits=c(0.25, 3), breaks=seq(0.5, 3, .5)) +
  scale_y_continuous(limits=c(20000, 90000)) +
  guides(size='none') +
  theme_avp()

m_p3 <- ggplot(pred3, aes(x=PTN, y=estimate)) +
  geom_point(data=NHANES_subset_wb$variables, aes(x=PTN, y=DXDTOBMD, size=weights_all), alpha=.1) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  #geom_rug(alpha=.1, sides="b") +
  labs(x="Protein intake (g)", y="Whole-body BMD (gm/cm²)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt3_an, hjust=1, vjust=-0.2,
           size=2.5) +
  scale_x_continuous(limits=c(25, 250)) +
  scale_y_continuous(limits=c(0.8, 1.5)) +
  guides(size='none') +
  theme_avp()

m_p4 <- ggplot(pred4, aes(x=PTNgkg, y=estimate)) +
  geom_point(data=NHANES_subset_wb$variables, aes(x=PTNgkg, y=DXDTOBMD, size=weights_all), alpha=.1) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  #geom_rug(alpha=.1, sides="b") +
  labs(x="Protein intake (g/kgBM)", y="Whole-body BMD (gm/cm²)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt4_an, hjust=1, vjust=-0.2,
           size=2.5) +
  scale_x_continuous(limits=c(0.25, 3), breaks=seq(0.5, 3, .5)) +
  scale_y_continuous(limits=c(0.8, 1.5)) +
  guides(size='none') +
  theme_avp()

m_p5 <- ggplot(pred5, aes(x=PTN, y=estimate)) +
  geom_point(data=NHANES_subset_femur$variables, aes(x=PTN, y=DXXOFBMD, size=weights_all), 
             alpha=.075) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  #geom_rug(alpha=.1, sides="b") +
  labs(x="Protein intake (g)", y="Femur BMD (gm/cm²)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt5_an, hjust=1, vjust=-0.2,
           size=2.5) +
  scale_x_continuous(limits=c(25, 250)) +
  scale_y_continuous(limits=c(0.4, 1.5), breaks=seq(0.4, 1.4, .2)) +
  guides(size='none') +
  theme_avp()

m_p6 <- ggplot(pred6, aes(x=PTNgkg, y=estimate)) +
  geom_point(data=NHANES_subset_femur$variables, aes(x=PTNgkg, y=DXXOFBMD, size=weights_all), alpha=.075) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  #geom_rug(alpha=.1, sides="b") +
  labs(x="Protein intake (g/kgBM)", y="Femur BMD (gm/cm²)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt6_an, hjust=1, vjust=-0.2,
           size=2.5) +
  scale_x_continuous(limits=c(0.25, 3), breaks=seq(0.5, 3, .5)) +
  scale_y_continuous(limits=c(0.4, 1.5), breaks=seq(0.4, 1.4, .2)) +
  guides(size='none') +
  theme_avp()

m_p7 <- ggplot(pred7, aes(x=PTN, y=estimate)) +
  geom_point(data=NHANES_subset_spine$variables, aes(x=PTN, y=DXXOSBMD, size=weights_all), alpha=.1) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  #geom_rug(alpha=.1, sides="b") +
  labs(x="Protein intake (g)", y="Spine BMD (gm/cm²)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt7_an, hjust=1, vjust=-0.2,
           size=2.5) +
  scale_x_continuous(limits=c(25, 250)) +
  scale_y_continuous(limits=c(.4, 1.6), breaks=seq(0.4, 1.6, .2)) +
  guides(size='none') +
  theme_avp()

m_p8 <- ggplot(pred8, aes(x=PTNgkg, y=estimate)) +
  geom_point(data=NHANES_subset_spine$variables, aes(x=PTNgkg, y=DXXOSBMD, size=weights_all), alpha=.1) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  #geom_rug(alpha=.1, sides="b") +
  labs(x="Protein intake (g/kgBM)", y="Spine BMD (gm/cm²)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt8_an, hjust=1, vjust=-0.2,
           size=2.5) +
  scale_x_continuous(limits=c(0.25, 3), breaks=seq(0.5, 3, .5)) +
  scale_y_continuous(limits=c(.4, 1.6), breaks=seq(0.4, 1.6, .2)) +
  guides(size='none') +
  theme_avp()

# assembling a pannel with patchwork
prp_all <- (m_p1 | m_p2) / (m_p3 | m_p4) / (m_p5 | m_p6) / (m_p7 | m_p8) & 
  plot_annotation(tag_levels = 'A')

prp_all

# ggsave(here("figures/marginal_plots_ALL_24_01_24_A.png"), units="in", width=9, height=14, dpi=600)

## DAG #########################################################################

# this is just a visual representation of the dag within R.
# the code was generated by creating the dag manually in dagitty.net
# the article features also the plot generated by the website.

dag <- dagitty::dagitty('dag { 
bb="0,0,1,1"
"Body weight" [pos="0.400,0.546"]
"GC use" [pos="0.325,0.546"]
"Lean mass/BMD" [outcome,pos="0.604,0.410"]
"Physical activity" [pos="0.444,0.301"]
Age [pos="0.355,0.301"]
Disease [pos="0.250,0.546"]
Energy [pos="0.473,0.546"]
Protein [exposure,pos="0.473,0.410"]
Sex [pos="0.250,0.301"]
"Body weight" -> "Lean mass/BMD" [pos="0.546,0.485"]
"Body weight" -> Energy
"Body weight" -> Protein
"GC use" -> "Body weight"
"GC use" -> "Lean mass/BMD" [pos="0.496,0.486"]
"GC use" -> Protein
"Physical activity" -> "Lean mass/BMD" [pos="0.533,0.342"]
"Physical activity" -> Protein
Age -> "Lean mass/BMD"
Age -> Protein
Disease -> "GC use"
Disease -> "Lean mass/BMD" [pos="0.472,0.452"]
Disease -> Protein [pos="0.355,0.450"]
Energy -> "Lean mass/BMD" [pos="0.558,0.504"]
Protein -> "Lean mass/BMD"
Protein -> Energy
Sex -> "Lean mass/BMD" [pos="0.517,0.193"]
Sex -> Protein [pos="0.381,0.384"]
}
')

plot(dag, ylim=c(-0.575, -0.230))

## Unadjusted scatterplots for supplementary material

## text df

model_summary_all_unadj <- table2_un

model_summary_all_unadj <- model_summary_all_unadj |> 
  mutate(coef = "raw")

unadj_models <- model_summary_all_std_unadj |> 
  mutate(term = rep(c("PTN", "PTNgkg"), 4),
         coef = "beta") |> 
  select(term, estimate, value2, p.value, type, coef) |> 
  full_join(model_summary_all_unadj) |> 
  mutate(p.value = case_when(p.value <= 0.0001 ~ "< 0.0001",
                             p.value > 0.0001 & p.value < 0.001 ~ "< 0.001",
                             p.value > 0.001 ~ format(round(p.value, digits=3), nsmall=3)))

unadj_models <- unadj_models |> 
  mutate(term = as.factor(rep(c("Lean mass ~ protein (g/d)",
                            "Lean mass ~ protein (g/kg)",
                            "Whole-body BMD ~ protein (g/d)",
                            "Whole-body BMD ~ protein (g/kg)",
                            "Femur BMD ~ protein (g/d)",
                            "Femur BMD ~ protein (g/kg)",
                            "Spine BMD ~ protein (g/d)",
                            "Spine BMD ~ protein (g/kg)"), 2))) |> 
  mutate(value2 = stringr::str_replace(value2, ";", " to"))

txt1 <- unadj_models |>
  filter(term=="Lean mass ~ protein (g/d)") 

txt2 <- unadj_models |>
  filter(term=="Whole-body BMD ~ protein (g/d)")

txt3 <- unadj_models |>
  filter(term=="Femur BMD ~ protein (g/d)")  

txt4 <- unadj_models |>
  filter(term=="Spine BMD ~ protein (g/d)")

txt1_an <- glue("Coefficient: {txt1$estimate[2]} g lean mass per 50 g protein (95%CI {txt1$value2[2]})\nβ: {txt1$estimate[1]} (95%CI {txt1$value2[1]}), p-value = {txt1$p.value[1]}")

txt2_an <- glue("Coefficient: {txt2$estimate[2]} gm/cm² BMD per 50 g protein (95%CI {txt2$value2[2]})\nβ: {txt2$estimate[1]} (95%CI {txt2$value2[1]}), p-value = {txt2$p.value[1]}")

txt3_an <- glue("Coefficient: {txt3$estimate[2]} gm/cm² BMD per 50 g protein (95%CI {txt3$value2[2]})\nβ: {txt3$estimate[1]} (95%CI {txt3$value2[1]}), p-value = {txt3$p.value[1]}")

txt4_an <- glue("Coefficient: {txt4$estimate[2]} gm/cm² BMD per 50 g protein (95%CI {txt4$value2[2]})\nβ: {txt4$estimate[1]} (95%CI {txt4$value2[1]}), p-value = {txt4$p.value[1]}")

# unadjusted models

model_leanmass_protein_un <- svyglm(DXDTOLE ~ PTN, 
                                     design=NHANES_subset_lm) 

model_wbBMD_protein_un <- svyglm(DXDTOBMD ~ PTN,
                                  design=NHANES_subset_wb) 

model_femurBMD_protein_un <- svyglm(DXXOFBMD ~ PTN, 
                                     design=NHANES_subset_femur) 

model_spineBMD_protein_un <- svyglm(DXXOSBMD ~ PTN, 
                                     design=NHANES_subset_spine) 

# unadjusted predictions

un_pred1 <- predictions(model_leanmass_protein_un, by=c("PTN"))
un_pred2 <- predictions(model_wbBMD_protein_un, by=c("PTN"))
un_pred3 <- predictions(model_femurBMD_protein_un, by=c("PTN"))
un_pred4 <- predictions(model_spineBMD_protein_un, by=c("PTN"))

# unadjusted plots

un_p1 <- ggplot(un_pred1, aes(x=PTN, y=estimate)) +
  geom_point(data=NHANES_subset_lm$variables, aes(x=PTN, y=DXDTOLE, size=weights_all), alpha=.1) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  labs(x="Protein intake (g)", y="Lean mass (g)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt1_an, hjust=1, vjust=-0.2,
           size=4) +
  scale_x_continuous(limits=c(25, 250)) +
  scale_y_continuous(limits=c(20000, 90000)) +
  guides(size='none') +
  theme_avp()

un_p2 <- ggplot(un_pred2, aes(x=PTN, y=estimate)) +
  geom_point(data=NHANES_subset_wb$variables, aes(x=PTN, y=DXDTOBMD, size=weights_all), alpha=.1) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  labs(x="Protein intake (g)", y="Whole-body BMD (gm/cm²)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt2_an, hjust=1, vjust=-0.2,
           size=4) +
  scale_x_continuous(limits=c(25, 250)) +
  scale_y_continuous(limits=c(0.8, 1.5)) +
  guides(size='none') +
  theme_avp()

un_p3 <- ggplot(un_pred3, aes(x=PTN, y=estimate)) +
  geom_point(data=NHANES_subset_femur$variables, aes(x=PTN, y=DXXOFBMD, size=weights_all), alpha=.1) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  labs(x="Protein intake (g)", y="Femur BMD (gm/cm²)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt3_an, hjust=1, vjust=-0.2,
           size=4) +
  scale_x_continuous(limits=c(25, 250)) +
  scale_y_continuous(limits=c(0.4, 1.5), breaks=seq(0.4, 1.4, .2)) +
  guides(size='none') +
  theme_avp()

un_p4 <- ggplot(un_pred4, aes(x=PTN, y=estimate)) +
  geom_point(data=NHANES_subset_spine$variables, aes(x=PTN, y=DXXOSBMD, size=weights_all), alpha=.1) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  labs(x="Protein intake (g)", y="Spine BMD (gm/cm²)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt4_an, hjust=1, vjust=-0.2,
           size=4) +
  scale_x_continuous(limits=c(25, 250)) +
  scale_y_continuous(limits=c(.4, 1.6), breaks=seq(0.4, 1.6, .2)) +
  guides(size='none') +
  theme_avp()

(un_p1 + un_p2 + un_p3 + un_p4) &
  plot_annotation(tag_levels = "A")
 
# ggsave(here("figures/scatter_unadjusted.png"), units="in", width=14, height=10, dpi=600, type="cairo")
