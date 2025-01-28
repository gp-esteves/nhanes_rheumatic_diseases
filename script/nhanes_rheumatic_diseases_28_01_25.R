## Analysis script for NHANES paper ##

# Author: Gabriel P. Esteves

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
library(gtsummary); library(marginaleffects); library(ggeffects);
library(splines)

# no scientific notation

options(scipen=999)
set.seed(1996)

#custom ggplot theme

theme_avp <- function() {
  theme_bw(base_size=10) +
    theme(axis.title.x = element_text(size = 9),
          axis.title.y = element_text(size = 9),
          panel.grid.minor = element_blank())
}

theme_sharp2 <- function() {
  theme_classic(base_size = 12) +
    theme(panel.grid.major.x = element_line(),
          plot.tag = element_text(face="bold"),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank())
}

# setting up a custom function for model predicitons

get_pred_df <- function(design, values, predictor, remove_outliers=TRUE) {
  if (predictor == "PTN") {
    
    if (remove_outliers == TRUE) {
      values <- subset(values, values < 250) 
    }
    
    newdata <- data.frame(   
      PTN = values,
      CHO = mean(design$variables$CHO),
      LIP = mean(design$variables$LIP),
      RIAGENDR = "1", 
      BMXWT = mean(design$variables$BMXWT),
      RIDAGEYR = mean(design$variables$RIDAGEYR),
      MET = mean(design$variables$MET),
      days_gc = mean(design$variables$days_gc),
      EI_EER_class = "Plausible reporter",
      disease_cat = "Osteoarthritis")
  } else {
    if (remove_outliers == TRUE) {
      values <- subset(values, values < 2.5) 
    }
    
    newdata <- data.frame(   
      PTNgkg = values,
      CHO = mean(design$variables$CHO),
      LIP = mean(design$variables$LIP),
      RIAGENDR = "1", 
      BMXWT = mean(design$variables$BMXWT),
      RIDAGEYR = mean(design$variables$RIDAGEYR),
      MET = mean(design$variables$MET),
      days_gc = mean(design$variables$days_gc),
      EI_EER_class = "Plausible reporter",
      disease_cat = "Osteoarthritis")
  }                       
  
  return(newdata)
}

# functions for customizing p values

remove_left_zero <- function(p_values) {
  stringr::str_replace_all(p_values, "^0", "")
}


p_format_g <- function(p.value) {
  p_value_tib <- tibble(p.value = p.value)
  
  new_p <- p_value_tib |> mutate(
    p.value = case_when(is.na(p.value) ~ "-", 
                        p.value >= 0.01 ~ remove_left_zero(format(round(p.value, 2), nsmall = 2)),
                        p.value < 0.01 & p.value > 0.001 ~ remove_left_zero(format(round(p.value, 3), nsmall = 3)),
                        p.value <= 0.001 ~ "<.001")
  )
  
  return(new_p$p.value)
}

# load dataset

load("data/nhanes_raw_dataframe.Rdata")

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
  mutate(femur_popmean_male = 1.04, # Reference values from NHANES III, Looker AC et al 1997.
         femur_popsd_male = 0.144,
         femurneck_popmean_male = 0.93,
         femurneck_popsd_male = 0.137,
         femur_popmean_female = 0.94,  
         femur_popsd_female = 0.122, 
         femurneck_popmean_female = 0.86,
         femurneck_popsd_female = 0.12,
         femur_t_score = case_when(
    RIAGENDR == "1" ~ ((DXXOFBMD-femur_popmean_male)/femur_popsd_male),
    RIAGENDR == "2" ~ ((DXXOFBMD-femur_popmean_female)/femur_popsd_female)
  )) |> 
  mutate(femurneck_t_score = case_when(
    RIAGENDR == "1" ~ ((DXXNKBMD-femurneck_popmean_male)/femurneck_popsd_male),
    RIAGENDR == "2" ~ ((DXXNKBMD-femurneck_popmean_female)/femurneck_popsd_female))) |> 
  mutate(femur_osteoporosis = as.factor(case_when(
    femur_t_score >= -2.5 ~ "No",
    femur_t_score < -2.5 ~ "Yes"
  )),
  femurneck_osteoporosis = as.factor(case_when(
    femurneck_t_score >= -2.5 ~ "No",
    femurneck_t_score < -2.5 ~ "Yes"
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
# herein we use the 24h recall day 1 weight, and divide by the number of cycles

rh_dxa_all <- rh_dxa_all |> 
  mutate(weights_all = case_when(
    DRDINT.x == 2 ~ (WTDR2D.x * (1/6)),
    DRDINT.x == 1 ~ (WTDRD1.x * (1/6))), 
    weights_all_exam = WTMEC2YR * (1/6))

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
  dplyr::select(SEQN, PTN, CHO, LIP, kcal, PTNgkg, 
                DXDTOLE, DXDTOBMD, DXXOFBMD, 
                DXXOSBMD, RIAGENDR, BMXHT, BMXWT, RIDAGEYR, disease_cat, 
                SDMVPSU, SDMVSTRA, weights_all, RIDRETH1, BMXBMI, femur_t_score, 
                femurneck_t_score, femur_osteoporosis, femurneck_osteoporosis, 
                DR1DRSTZ, DXAEXSTS, DXAFMRST, DXASPNST, DXXNKBMD, MET, MET_class,
                days_gc, EI_EER_class, cycle)

# this dataset can then be saved and loaded when you actually want to analyze it.
# save(rh_dxa_all_analyse, file="data/NHANES_RH.Rdata")

# set survey
# set individuals without sample weights to 0, to exclude later

rh_dxa_all_analyse <- rh_dxa_all_analyse |> 
  mutate(weights_all = if_else(is.na(weights_all), 0, weights_all)) |> 
  mutate(RIAGENDR = as.factor(RIAGENDR),
         EI_EER_class = as.factor(EI_EER_class),
         cycle = as.factor(cycle)) # set categorical variables to factor

NHANES_all <- svydesign(data=rh_dxa_all_analyse, 
                        id=~SDMVPSU, 
                        strata=~SDMVSTRA, 
                        weights=~weights_all,
                        nest=TRUE)

nrow(NHANES_all)

# subsetting
# quick codebook for subset: MCQ195 = 1, 2, 3 (osteoarth, RA, psoriatic A patients); 
#                            MCQ190, MCQ191 = 1 (RA patients, different documentation);
#                            DXAEXSTS; DXASPNST; DXAFMRST  = 1 (adequate DXA exam);  

NHANES_rh <- subset(NHANES_all, weights_all != 0) # excluding individuals without sample weight

NHANES_rh <- subset(NHANES_all, !is.na(disease_cat)) # selecting RDs 

NHANES_rh <- subset(NHANES_rh, DR1DRSTZ == "1") # selecting adequate r24 data

NHANES_rh <- subset(NHANES_rh, RIDAGEYR >= 18 & RIDAGEYR <= 85) # selecting adults 

nrow(NHANES_rh)

# adults older than 85 are still coded as 85 in NHANES.

NHANES_subset_lm <- subset(NHANES_rh, DXAEXSTS == "1" & !is.na(DXDTOLE))

NHANES_subset_wb <- subset(NHANES_rh, DXAEXSTS == "1" & !is.na(DXDTOBMD))

NHANES_subset_femur <- subset(NHANES_rh, DXAFMRST == "1" & !is.na(DXXOFBMD))

NHANES_subset_spine <- subset(NHANES_rh, DXASPNST == "1" & !is.na(DXXOSBMD))

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

# separate dfs for each outcome, for table 1

lm_dat <- as_tibble(NHANES_subset_lm$variables)

wb_dat <- as_tibble(NHANES_subset_wb$variables)

fem_dat <- as_tibble(NHANES_subset_femur$variables)

spn_dat <- as_tibble(NHANES_subset_spine$variables)

sample_dat <- reduce(list(lm_dat, wb_dat, fem_dat, spn_dat), bind_rows) |> distinct()

nrow(sample_dat) # sample size

## Modelling ###################################################################

# unadjusted models

# Lean mass

u_model_leanmass_protein <- svyglm(DXDTOLE ~ PTN, # protein in grams
                                 design=NHANES_subset_lm) # model fit
summary(u_model_leanmass_protein) # summary
u_model1<-broom::tidy(u_model_leanmass_protein, conf.int=T) # broom::tidy extracts model coefficients and CIs. Doing this for all models.

# BMD

# whole-body

u_model_wbBMD_protein <- svyglm(DXDTOBMD ~ PTN, 
                              design=NHANES_subset_wb) 
summary(u_model_wbBMD_protein) 
u_model2<-broom::tidy(u_model_wbBMD_protein,conf.int=T)

# femoral neck

u_model_femurneckBMD_protein <- svyglm(DXXNKBMD ~ PTN, design=NHANES_subset_femur) 
summary(u_model_femurneckBMD_protein) 
u_model3<-broom::tidy(u_model_femurneckBMD_protein,conf.int=T)

# spine

u_model_spineBMD_protein <- svyglm(DXXOSBMD ~ PTN, design=NHANES_subset_spine) 
summary(u_model_spineBMD_protein)
u_model4<-broom::tidy(u_model_spineBMD_protein,conf.int=T)

# creating a dataframe with all results and tidying 
# for later use in plots and summaries

table2_un <- reduce(list(u_model1,u_model2,u_model3,u_model4),
                            full_join) |> 
  filter(term=="PTN") |> 
  mutate(type=c("lean mass", "wb BMD", "femur BMD", "spine BMD")) |> 
  mutate(across(c(estimate, conf.low, conf.high),
                .fns = ~ if_else(term == "PTN", .x * 50, .x))) |> 
  mutate(across(c(estimate, conf.low, conf.high), .f = ~ if_else(type=="lean mass", .x / 1000, .x))) |> 
  mutate(value=paste(estimate," (",conf.low,"; ",conf.high,")",sep="")) |> 
  mutate(estimate=if_else(type == "lean mass", round(estimate, 2), round(estimate, 4))) |> 
  mutate(conf.low=if_else(type == "lean mass", round(conf.low, 2), round(conf.low, 4))) |> 
  mutate(conf.high=if_else(type == "lean mass", round(conf.high, 2), round(conf.high, 4))) |> 
  mutate(value2=paste(conf.low,"; ",conf.high,sep="")) |> 
  mutate(p.value = round(p.value, 4)) |> 
  select(term, estimate, value2, p.value, type) 

# write_xlsx(table2_un, here("Tables/nhanes_model_summary_unadjusted.xlsx"))

# adjusted models

# Lean mass

model_leanmass_protein <- svyglm(DXDTOLE ~ PTN + CHO + LIP + RIAGENDR + BMXWT # same process
                                 + RIDAGEYR + MET + days_gc + EI_EER_class + disease_cat, # but now includes covariates determined by DAG
                              design=NHANES_subset_lm) 

summary(model_leanmass_protein) 
model1<-broom::tidy(model_leanmass_protein, conf.int=T)

# visualizing partial residuals to assess model fit and non-linearity
plot(predict_response(model_leanmass_protein, "PTN [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE) +
  labs(y="Lean mass", title=NULL)

model_leanmass_proteingkg <- svyglm(DXDTOLE ~ PTNgkg + CHO + LIP + RIAGENDR + 
                                      BMXWT + RIDAGEYR + MET + days_gc + EI_EER_class + disease_cat, 
                                 design=NHANES_subset_lm) 
summary(model_leanmass_proteingkg) 
model2<-broom::tidy(model_leanmass_proteingkg,conf.int=T)

plot(predict_response(model_leanmass_proteingkg, "PTNgkg [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE) +
  labs(y="Lean mass", title=NULL)

# the blue line shows a loess fit in the residuals. the drop towards the end of the curve suggests
# a non-linear association between protein(gkg) and lean mass.
# we fit a spline to better test this.

# attempting to fit a polynomial term to assess non-linearity

model_leanmass_proteingkg_poly <- svyglm(DXDTOLE ~ poly(PTNgkg, 2) + CHO + LIP + RIAGENDR + 
                                         BMXWT + RIDAGEYR + MET + days_gc + EI_EER_class + disease_cat, 
                                       design=NHANES_subset_lm) 

summary(model_leanmass_proteingkg_poly) 
car::Anova(model_leanmass_proteingkg_poly, type="II", test="F") # we now use an F test via a type II anova table to check for significance of the polynomial term as a whole

anova(model_leanmass_proteingkg, model_leanmass_proteingkg_poly) # Likelihood ratio test provided by survey package to compare models 

plot(predict_response(model_leanmass_proteingkg_poly, "PTNgkg [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE) +
  labs(y="Lean mass", title=NULL)

plot_predictions(model_leanmass_proteingkg_poly,
                 newdata=get_pred_df(NHANES_subset_lm, NHANES_subset_lm$variables$PTNgkg,
                                     predictor="PTNgkg"),
                 by="PTNgkg") ## a version of the plot without the raw data can be seen here.

# results suggest that the polynomial model provides a better fit for predicting
# lean mass by PTNgkg.

# BMD

# whole-body

model_wbBMD_protein <- svyglm(DXDTOBMD ~ PTN + CHO + LIP + RIAGENDR + BMXWT + RIDAGEYR
                              + MET + days_gc + EI_EER_class + disease_cat, 
                              design=NHANES_subset_wb) 
summary(model_wbBMD_protein) # total body BMD and total protein
model3<-broom::tidy(model_wbBMD_protein,conf.int=T)

plot(predict_response(model_wbBMD_protein, "PTN [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE) +
  labs(y="Whole-body BMD", title=NULL)

model_wbBMD_proteingkg <- svyglm(DXDTOBMD~PTNgkg  + CHO + LIP + RIAGENDR + BMXWT + 
                                   RIDAGEYR + MET + days_gc + EI_EER_class + disease_cat, 
                                 design=NHANES_subset_wb) 
summary(model_wbBMD_proteingkg) # total body BMD and protein g/kg 
model4<-broom::tidy(model_wbBMD_proteingkg,conf.int=T)

plot(predict_response(model_wbBMD_proteingkg, "PTNgkg [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE) +
  labs(y="Whole-body BMD", title=NULL)

# both plots show some non-linearity towards the higher values of PTN and PTNgkg,
# so we attempt polynomial terms.

# polynomial WB BMD models

model_wbBMD_protein_poly <- svyglm(DXDTOBMD ~ poly(PTN, 2) + CHO + LIP + RIAGENDR + BMXWT + RIDAGEYR
                                   + MET + days_gc + EI_EER_class + disease_cat, 
                                   design=NHANES_subset_wb) 

summary(model_wbBMD_protein_poly) 
car::Anova(model_wbBMD_protein_poly, type="II", test="F")

AIC(model_wbBMD_protein, model_wbBMD_protein_poly)
anova(model_wbBMD_protein, model_wbBMD_protein_poly) 

plot(predict_response(model_wbBMD_protein_poly, terms="PTN"), 
     show_residuals=TRUE, show_residuals_line = TRUE, limit_range=TRUE) +
  labs(y="Whole-body BMD", title=NULL)

model_wbBMD_proteingkg_poly <- svyglm(DXDTOBMD ~ poly(PTNgkg, 2) + CHO + LIP + RIAGENDR + BMXWT + RIDAGEYR
                                   + MET + days_gc + EI_EER_class + disease_cat, 
                                   design=NHANES_subset_wb) 

summary(model_wbBMD_proteingkg_poly) 
car::Anova(model_wbBMD_proteingkg_poly, type="II", test="F")

AIC(model_wbBMD_proteingkg, model_wbBMD_proteingkg_poly)
anova(model_wbBMD_proteingkg, model_wbBMD_proteingkg_poly) 

plot(predict_response(model_wbBMD_proteingkg_poly, "PTNgkg [all]", back_transform = TRUE), 
     show_residuals=TRUE, show_residuals_line = TRUE) +
  labs(y="Whole-body BMD", title=NULL)

# although the association with protein remains non-significant,
# the partial residual plot and the AIC suggest better fit of the polynomial model
# in both cases.

# femur
# we modelled femur first, then femoral neck (more clinically used and relevant)
# total femur is not reported in the paper but results are extremelly similar.

model_femurBMD_protein <- svyglm(DXXOFBMD ~ PTN + CHO + LIP + 
                                   RIAGENDR + BMXWT + RIDAGEYR
                                 + MET + days_gc + EI_EER_class + disease_cat, 
                                 design=NHANES_subset_femur) 
summary(model_femurBMD_protein) # femur BMD and total protein
model5<-broom::tidy(model_femurBMD_protein,conf.int=T)

plot(predict_response(model_femurBMD_protein, "PTN [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE) +
  labs(y="Femur BMD", title=NULL)

model_femurBMD_proteingkg <- svyglm(DXXOFBMD ~ PTNgkg + CHO + LIP +
                                      RIAGENDR + BMXWT + RIDAGEYR
                                    + MET + days_gc + EI_EER_class + disease_cat, 
                                    design=NHANES_subset_femur) 
summary(model_femurBMD_proteingkg) # femur BMD and protein g/kg
model6<-broom::tidy(model_femurBMD_proteingkg,conf.int=T)

plot(predict_response(model_femurBMD_proteingkg, "PTNgkg [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE) +
  labs(y="Femur BMD", title=NULL)

## both plots show some non-linearity. the ptngkg plot even suggest a small downfall towards the end
## of the x axis, suggesting a more complex pattern. we fit third degree polynomial and natural splines.
## attempting best fits.

model_femurBMD_protein_poly <- svyglm(DXXOFBMD ~ poly(PTN, 2) + CHO + LIP + 
                                   RIAGENDR + BMXWT + RIDAGEYR
                                 + MET + days_gc + EI_EER_class + disease_cat, 
                                 design=NHANES_subset_femur) 
summary(model_femurBMD_protein_poly)
car::Anova(model_femurBMD_protein_poly, type="II")

plot(predict_response(model_femurBMD_protein_poly, "PTN [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE) +
  labs(y="Femur BMD", title=NULL)

## femur ptngkg

model_femurBMD_proteingkg_poly_2 <- svyglm(DXXOFBMD ~ poly(PTNgkg, 2) + CHO + LIP +
                                      RIAGENDR + BMXWT + RIDAGEYR
                                    + MET + days_gc + EI_EER_class + disease_cat, 
                                    design=NHANES_subset_femur) 

model_femurBMD_proteingkg_poly_3 <- svyglm(DXXOFBMD ~ poly(PTNgkg, 3) + CHO + LIP +
                                             RIAGENDR + BMXWT + RIDAGEYR
                                           + MET + days_gc + EI_EER_class+ disease_cat, 
                                           design=NHANES_subset_femur)

model_femurBMD_proteingkg_ns <- svyglm(DXXOFBMD ~ ns(PTNgkg, 3) + CHO + LIP +
                                             RIAGENDR + BMXWT + RIDAGEYR
                                           + MET + days_gc + EI_EER_class + disease_cat, 
                                       design=NHANES_subset_femur)

car::Anova(model_femurBMD_proteingkg_poly_2, type="II", test="F")
car::Anova(model_femurBMD_proteingkg_poly_3, type="II", test="F")
car::Anova(model_femurBMD_proteingkg_ns, type="II", test="F")

AIC(model_femurBMD_proteingkg,
    model_femurBMD_proteingkg_poly_2,
    model_femurBMD_proteingkg_poly_3,
    model_femurBMD_proteingkg_ns)

plot(predict_response(model_femurBMD_proteingkg_poly_2, "PTNgkg [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE) +
  labs(y="Femur BMD", title=NULL)

plot(predict_response(model_femurBMD_proteingkg_poly_3, "PTNgkg [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE) +
  labs(y="Femur BMD", title=NULL)

plot(predict_response(model_femurBMD_proteingkg_ns, "PTNgkg [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE) +
  labs(y="Femur BMD", title=NULL)

plot_predictions(model_femurBMD_proteingkg_ns,
                 newdata=get_pred_df(NHANES_subset_femur, NHANES_subset_femur$variables$PTNgkg,
                                     predictor="PTNgkg"),
                 by="PTNgkg")

## all models have statistically significant associations with ptnkg,
## with the poly_3 and ns models having better AIC values.
## however, given the partial residual plot, the poly_3 model seems to be overestimating
## the curve towards the end. the natural spline model appears to have a better fit in that sense.
## we select the NS model as our final model for FEMUR BMD and PTNgkg.

## looking at femural neck (more clinically relevant), follows same strategy as in total femur

model_fneckBMD_protein <- svyglm(DXXNKBMD ~ PTN + CHO + LIP +
                                      RIAGENDR + BMXWT + RIDAGEYR
                                    + MET + days_gc + EI_EER_class + disease_cat, 
                                 design=NHANES_subset_femur)

model_fneckBMD_protein_poly_2 <- svyglm(DXXNKBMD ~ poly(PTN, 2) + CHO + LIP +
                                    RIAGENDR + BMXWT + RIDAGEYR
                                  + MET + days_gc + EI_EER_class + disease_cat, 
                                  design=NHANES_subset_femur)

model_fneckBMD_protein_poly_3 <- svyglm(DXXNKBMD ~ poly(PTN, 3) + CHO + LIP +
                                           RIAGENDR + BMXWT + RIDAGEYR
                                         + MET + days_gc + EI_EER_class + disease_cat, 
                                        design=NHANES_subset_femur)

model_fneckBMD_protein_ns <- svyglm(DXXNKBMD ~ ns(PTN, 3) + CHO + LIP +
                                           RIAGENDR + BMXWT + RIDAGEYR
                                         + MET + days_gc + EI_EER_class + disease_cat, 
                                    design=NHANES_subset_femur)

AIC(model_fneckBMD_protein,
    model_fneckBMD_protein_poly_2,
    model_fneckBMD_protein_poly_3,
    model_fneckBMD_protein_ns)

car::Anova(model_fneckBMD_protein, type="II", test="F")
car::Anova(model_fneckBMD_protein_poly_2, type="II", test="F")
car::Anova(model_fneckBMD_protein_poly_3, type="II", test="F")
car::Anova(model_fneckBMD_protein_ns, type="II", test="F")

plot(predict_response(model_fneckBMD_protein, "PTN [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE) +
  labs(y="Femoral Neck BMD", title=NULL)

plot(predict_response(model_fneckBMD_protein_poly_2, "PTN [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE) +
  labs(y="Femoral Neck BMD", title=NULL)

plot(predict_response(model_fneckBMD_protein_poly_3, "PTN [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE) +
  labs(y="Femoral Neck BMD", title=NULL)

plot(predict_response(model_fneckBMD_protein_ns, "PTN [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE) +
  labs(y="Femoral Neck BMD", title=NULL)

# ptngkg

model_fneckBMD_proteingkg <- svyglm(DXXNKBMD ~ PTNgkg + CHO + LIP +
                                         RIAGENDR + BMXWT + RIDAGEYR
                                       + MET + days_gc + EI_EER_class + disease_cat, 
                                    design=NHANES_subset_femur)

model_fneckBMD_proteingkg_poly_2 <- svyglm(DXXNKBMD ~ poly(PTNgkg, 2) + CHO + LIP +
                                      RIAGENDR + BMXWT + RIDAGEYR
                                    + MET + days_gc + EI_EER_class + disease_cat, 
                                    design=NHANES_subset_femur)

model_fneckBMD_proteingkg_poly_3 <- svyglm(DXXNKBMD ~ poly(PTNgkg, 3) + CHO + LIP +
                                      RIAGENDR + BMXWT + RIDAGEYR
                                    + MET + days_gc + EI_EER_class + disease_cat, 
                                    design=NHANES_subset_femur)


model_fneckBMD_proteingkg_ns <- svyglm(DXXNKBMD ~ ns(PTNgkg, 3) + CHO + LIP +
                                      RIAGENDR + BMXWT + RIDAGEYR
                                    + MET + days_gc + EI_EER_class + disease_cat, 
                                    design=NHANES_subset_femur)

AIC(model_fneckBMD_proteingkg,
    model_fneckBMD_proteingkg_poly_2,
    model_fneckBMD_proteingkg_poly_3,
    model_fneckBMD_proteingkg_ns)

car::Anova(model_fneckBMD_proteingkg, type="II", test="F")
car::Anova(model_fneckBMD_proteingkg_poly_2, type="II", test="F")
car::Anova(model_fneckBMD_proteingkg_poly_3, type="II", test="F")
car::Anova(model_fneckBMD_proteingkg_ns, type="II", test="F")

plot(predict_response(model_fneckBMD_proteingkg, "PTNgkg [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE) +
  labs(y="Femoral Neck BMD", title=NULL)

plot(predict_response(model_fneckBMD_proteingkg_poly_2, "PTNgkg [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE) +
  labs(y="Femoral Neck", title=NULL)

plot(predict_response(model_fneckBMD_proteingkg_poly_3, "PTNgkg [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE) +
  labs(y="Femoral Neck BMD", title=NULL)

plot(predict_response(model_fneckBMD_proteingkg_ns, "PTNgkg [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE) +
  labs(y="Femoral Neck", title=NULL)

# spine

model_spineBMD_protein <- svyglm(DXXOSBMD ~ PTN + CHO + LIP + 
                                   RIAGENDR + BMXWT + RIDAGEYR
                                 + MET + days_gc + EI_EER_class + disease_cat, 
                                 design=NHANES_subset_spine) 
summary(model_spineBMD_protein) # spine BMD and total protein
model7<-broom::tidy(model_spineBMD_protein,conf.int=T)

plot(predict_response(model_spineBMD_protein, "PTN [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE) +
  labs(y="Spine BMD", title=NULL)

model_spineBMD_proteingkg <- svyglm(DXXOSBMD ~ PTNgkg + CHO + LIP + 
                                      RIAGENDR + BMXWT + RIDAGEYR
                                    + MET + days_gc + EI_EER_class + disease_cat, 
                                    design=NHANES_subset_spine) 
summary(model_spineBMD_proteingkg) # spine BMD and protein g/kg
model8<-broom::tidy(model_spineBMD_proteingkg,conf.int=T)

plot(predict_response(model_spineBMD_proteingkg, "PTNgkg [all]"), 
     show_residuals=TRUE, show_residuals_line = TRUE) +
  labs(y="Spine BMD", title=NULL)

## for the spine models, the linear trend appears to be accurate, and so
## we don't fit any more complex associations.

# building a dataframe with results from linear models

model_summary_all <- reduce(list(model1,model2,model3,model4,
                                 model5,model6,model7,model8),
                            full_join) |> 
  filter(term=="PTN" | term=="PTNgkg") |> 
  mutate(type=c("lean mass","lean mass","wb BMD","wb BMD","femur BMD","femur BMD",
                 "spine BMD","spine BMD")) |> 
  mutate(across(c(estimate, conf.low, conf.high),
                  .fns = ~ if_else(term == "PTN", .x * 50, .x))) |> 
  mutate(across(c(estimate, conf.low, conf.high),
                .fns = ~ if_else(type == "lean mass", round(.x/1000, 2), .x))) |>
  mutate(value=paste(estimate," (",conf.low,"; ",conf.high,")",sep="")) |> 
  mutate(estimate=round(estimate, 4)) |> 
  mutate(conf.low=round(conf.low, 4)) |> 
  mutate(conf.high=round(conf.high, 4)) |> 
  mutate(value=paste(conf.low,"; ",conf.high,sep="")) |> 
  select(term, estimate, value, p.value, type, conf.low, conf.high)

# code to save the table in Excel format. Commented out to avoid saving multiple files.
# write_xlsx(model_summary_all, here("Tables/nhanes_model_summary.xlsx"))

#### Tables ####################################################################

# creating table 1 with gtsummary

table_1_summary <- sample_dat |>
 select(disease_cat, RIDAGEYR, RIAGENDR, RIDRETH1, BMXHT, BMXWT, BMXBMI, days_gc, 
        MET, MET_class, kcal, CHO, PTN, PTNgkg, LIP,
        DXDTOLE, DXDTOBMD, DXXNKBMD, femurneck_t_score, femurneck_osteoporosis, DXXOSBMD) |>
 tbl_summary(by = disease_cat, missing = "no", statistic = all_continuous() ~ "{mean} ({sd})",
             digits = list(RIDAGEYR ~ 0, 
                           BMXWT ~ 1,
                           BMXHT ~ 1,
                           BMXBMI ~ 1,
                           days_gc ~ 0,
                           MET ~ 0,
                           kcal ~ 0,
                           PTN ~ 0,
                           PTNgkg ~ 2,
                           LIP ~ 0,
                           DXDTOLE ~ 1,
                           DXDTOBMD ~ 3,
                           DXXNKBMD ~ 3,
                           femurneck_t_score ~ 1,
                           DXXOSBMD ~ 3)) |>
  add_overall()

# save
# table_1_summary |> as_hux_table() |> writexl::write_xlsx("tables/table1.xlsx")

#### Plots #####################################################################

## prediction plots

# creating texts for plots

summa <- model_summary_all |> 
  mutate(estimate = if_else(type == "lean mass", round(estimate, 2), round(estimate, 4))) |> 
  mutate(conf.low = if_else(type == "lean mass", round(conf.low, 2), round(conf.low, 4))) |> 
  mutate(conf.high = if_else(type == "lean mass", round(conf.high, 2), round(conf.high, 4))) |> 
  mutate(p.value = if_else(p.value >= 0.05, round(p.value, 2), round(p.value, 5)))

m2_anova <- car::Anova(model_leanmass_proteingkg_poly, type="II")
m3_anova <- car::Anova(model_wbBMD_protein_poly, type="II")
m4_anova <- car::Anova(model_wbBMD_proteingkg_poly, type="II")
m5_anova <- car::Anova(model_femurBMD_protein_poly, type="II")
m5_anova_neck <- car::Anova(model_fneckBMD_protein_poly_2, type="II")
m6_anova <- car::Anova(model_femurBMD_proteingkg_ns, type="II")
m6_anova_neck <- car::Anova(model_fneckBMD_proteingkg_ns, type="II")

txt1_an <- glue("Coefficient: {summa$estimate[1]} kg lean mass per 50 g protein (95%CI {summa$conf.low[1]} to {summa$conf.high[1]})\np-value {p_format_g(summa$p.value[1])}")
txt2_an <- glue("Polynomial term p-value {p_format_g(m2_anova$`Pr(>Chisq)`[1])}")
txt3_an <- glue("Polynomial term p-value = {p_format_g(m3_anova$`Pr(>Chisq)`[1])}")
txt4_an <- glue("Polynomial term p-value = {p_format_g(m4_anova$`Pr(>Chisq)`[1])}")
txt5_an <- glue("Polynomial term p-value {p_format_g(m5_anova$`Pr(>Chisq)`[1])}")
txt5_an_neck <- glue("Polynomial term p-value {p_format_g(m5_anova_neck$`Pr(>Chisq)`[1])}")
txt6_an <- glue("Spline term p-value = {p_format_g(m6_anova$`Pr(>Chisq)`[1])}")
txt6_an_neck <- glue("Spline term p-value = {p_format_g(m6_anova_neck$`Pr(>Chisq)`[1])}")
txt7_an <- glue("Coefficient: {summa$estimate[7]} gm/cm² BMD per 50 g protein (95%CI {summa$conf.low[7]} to {summa$conf.high[7]})\np-value = {p_format_g(summa$p.value[7])}")
txt8_an <- glue("Coefficient: {summa$estimate[8]} gm/cm² BMD per 1 g/kgBM protein (95%CI {summa$conf.low[8]} to {summa$conf.high[8]})\np-value = {p_format_g(summa$p.value[8])}")

## prediction plots with marginaleffects

# get a newdata df, using original observations of protein intake, and holding
# covariates at their mean

pred_data_1 <- get_pred_df(NHANES_subset_lm, NHANES_subset_lm$variables$PTN, "PTN")
pred_data_2 <- get_pred_df(NHANES_subset_lm, NHANES_subset_lm$variables$PTNgkg, "PTNgkg")
pred_data_3 <- get_pred_df(NHANES_subset_wb, NHANES_subset_wb$variables$PTN, "PTN")
pred_data_4 <- get_pred_df(NHANES_subset_wb, NHANES_subset_wb$variables$PTNgkg, "PTNgkg")
pred_data_5 <- get_pred_df(NHANES_subset_femur, NHANES_subset_femur$variables$PTN, "PTN")
pred_data_6 <- get_pred_df(NHANES_subset_femur, NHANES_subset_femur$variables$PTNgkg, "PTNgkg")
pred_data_7 <- get_pred_df(NHANES_subset_spine, NHANES_subset_spine$variables$PTN, "PTN")
pred_data_8 <- get_pred_df(NHANES_subset_spine, NHANES_subset_spine$variables$PTNgkg, "PTNgkg")

# function to retrieve weights correctly

get_weights <- function(design, predictor) {
  df_weights <- tibble(
    ptn = design$variables$PTN,
    ptngkg = design$variables$PTNgkg,
    weights = design$variables$weights_all
    ) 
  
  if (predictor == "PTN") {
   result <- df_weights |> 
      filter(ptn < 250) |> 
      pull(weights)
  } else {
    result <- df_weights |> 
      filter(ptngkg < 2.5) |> 
      pull(weights)
  }
  
  return(result)
}

# make predictions

# model list

pred1 <- predictions(model_leanmass_protein, # predicts the linear relationship between protein and outcome using the model
                   by=c("PTN"), newdata=pred_data_1, # averages using the mean of covariates
                   wts=get_weights(NHANES_subset_lm, "PTN")) # use weights for accurate predictions

pred2 <- predictions(model_leanmass_proteingkg_poly, # same thing but here we use the polynomial iteration of the model
                   by=c("PTNgkg"), newdata=pred_data_2,
                   wts=get_weights(NHANES_subset_lm, "PTNgkg"))

pred3 <- predictions(model_wbBMD_protein_poly,
                   by=c("PTN"), newdata=pred_data_3,
                   wts=get_weights(NHANES_subset_wb, "PTN"))

pred4 <- predictions(model_wbBMD_proteingkg_poly,
                   by=c("PTNgkg"), newdata=pred_data_4,
                   wts=get_weights(NHANES_subset_wb, "PTNgkg"))

pred5 <- predictions(model_femurBMD_protein_poly,
                   by=c("PTN"), newdata=pred_data_5,
                   wts=get_weights(NHANES_subset_femur, "PTN"))

pred5_neck <- predictions(model_fneckBMD_protein_poly_2,
                          by=c("PTN"), newdata=pred_data_5,
                          wts=get_weights(NHANES_subset_femur, "PTN"))

pred6 <- predictions(model_femurBMD_proteingkg_ns,
                   by=c("PTNgkg"), newdata=pred_data_6,
                   wts=get_weights(NHANES_subset_femur, "PTNgkg"))

pred6_neck <- predictions(model_fneckBMD_proteingkg_ns,
                          by=c("PTNgkg"), newdata=pred_data_6,
                          wts=get_weights(NHANES_subset_femur, "PTNgkg"))

pred7 <- predictions(model_spineBMD_protein,
                   by=c("PTN"), newdata=pred_data_7,
                   wts=get_weights(NHANES_subset_spine, "PTN"))

pred8 <- predictions(model_spineBMD_proteingkg,
                   by=c("PTNgkg"), newdata=pred_data_8,
                   wts=get_weights(NHANES_subset_spine, "PTNgkg"))


m_p1 <- ggplot(pred1, aes(x=PTN, y=estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  labs(x="Protein intake (g/day)", y="Lean mass (kg)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt1_an, hjust=1, vjust=-0.2,
           size=2.5) +
  scale_x_continuous(limits=c(25, 250)) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-3),
                     limits=c(56000, 68000), breaks=seq(56000, 68000, 2000)) +
  guides(size='none') +
  theme_avp()

m_p2 <- ggplot(pred2, aes(x=PTNgkg, y=estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  labs(x="Protein intake (g/kgBM/day)", y="Lean mass (kg)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt2_an, hjust=1, vjust=-0.6,
           size=2.5) +
  scale_x_continuous(limits=c(0.25, 2.5), breaks=seq(0.5, 2.5, .5)) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-3),
                     limits=c(56000, 68000), breaks=seq(56000, 68000, 2000)) +
  guides(size='none') +
  theme_avp()

m_p3 <- ggplot(pred3, aes(x=PTN, y=estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  labs(x="Protein intake (g/day)", y="Whole-body BMD (gm/cm²)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt3_an, hjust=1, vjust=-0.6,
           size=2.5) +
  scale_x_continuous(limits=c(25, 250)) +
  scale_y_continuous(limits=c(1.05, 1.275), breaks=seq(1.05, 1.275, .05)) +
  guides(size='none') +
  theme_avp()

m_p4 <- ggplot(pred4, aes(x=PTNgkg, y=estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  labs(x="Protein intake (g/kgBM/day)", y="Whole-body BMD (gm/cm²)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt4_an, hjust=1,vjust=-0.6,
           size=2.5) +
  scale_x_continuous(limits=c(0.25, 2.5), breaks=seq(0.5, 2.5, .5)) +
  scale_y_continuous(limits=c(1.05, 1.275), breaks=seq(1.05, 1.275, .05)) +
  guides(size='none') +
  theme_avp()

m_p5 <- ggplot(pred5, aes(x=PTN, y=estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  labs(x="Protein intake (g/day)", y="Femur BMD (gm/cm²)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt5_an, hjust=1,vjust=-0.6,
           size=2.5) +
  scale_x_continuous(limits=c(25, 250)) +
  #scale_y_continuous(limits=c(0.4, 1.5), breaks=seq(0.4, 1.4, .2)) +
  guides(size='none') +
  theme_avp()

m_p5_neck <- ggplot(pred5_neck, aes(x=PTN, y=estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  #geom_rug(alpha=.1, sides="b") +
  labs(x="Protein intake (g/day)", y="Femoral neck BMD (gm/cm²)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt5_an_neck, hjust=1, vjust=-0.2,
           size=2.5) +
  scale_x_continuous(limits=c(25, 250)) +
  scale_y_continuous(limits=c(0.59, 1.0), breaks=seq(0.6, 1.0, 0.1)) +
  guides(size='none') +
  theme_avp()

m_p6 <- ggplot(pred6, aes(x=PTNgkg, y=estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  labs(x="Protein intake (g/kgBM/day)", y="Femur BMD (gm/cm²)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt6_an, hjust=1, vjust=-0.6,
           size=2.5) +
  scale_x_continuous(limits=c(0.25, 2.5), breaks=seq(0.5, 2.5, .5)) +
  guides(size='none') +
  theme_avp()

m_p6_neck <- ggplot(pred6_neck, aes(x=PTNgkg, y=estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  #geom_rug(alpha=.1, sides="b") +
  labs(x="Protein intake (g/kgBM/day)", y="Femoral neck BMD (gm/cm²)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt6_an_neck, hjust=1, vjust=-0.2,
           size=2.5) +
  scale_x_continuous(limits=c(0.25, 2.5), breaks=seq(0.5, 2.5, .5)) +
  scale_y_continuous(limits=c(0.59, 1.0), breaks=seq(0.6, 1.0, 0.1)) +
  guides(size='none') +
  theme_avp()

m_p7 <- ggplot(pred7, aes(x=PTN, y=estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  labs(x="Protein intake (g/day)", y="Spine BMD (gm/cm²)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt7_an, hjust=1, vjust=-0.2,
           size=2.5) +
  scale_x_continuous(limits=c(25, 250)) +
  scale_y_continuous(limits=c(0.975, 1.15), breaks=seq(0.95, 1.15, 0.05)) +
  guides(size='none') +
  theme_avp()

m_p8 <- ggplot(pred8, aes(x=PTNgkg, y=estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  labs(x="Protein intake (g/kgBM/day)", y="Spine BMD (gm/cm²)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt8_an, hjust=1, vjust=-0.2,
           size=2.5) +
  scale_x_continuous(limits=c(0.25, 2.5), breaks=seq(0.5, 2.5, .5)) +
  scale_y_continuous(limits=c(0.975, 1.15), breaks=seq(0.95, 1.15, 0.05)) +
  guides(size='none') +
  theme_avp()

# assembling a pannel with patchwork
prp_all <- (m_p1 | m_p2) / (m_p3 | m_p4) / (m_p5_neck | m_p6_neck) / (m_p7 | m_p8) & 
  plot_annotation(tag_levels = 'A')

prp_all

# ggsave(here("figures/marginal_plots_ALL.png"), units="in", width=9, height=14, dpi=600)

#### analyzing t-scores ###############################

model_tscore_protein_poly <- svyglm(femurneck_t_score ~ poly(PTN, 2) + CHO + LIP +
                                      RIAGENDR + BMXWT + RIDAGEYR
                                    + MET + days_gc + EI_EER_class + disease_cat, 
                                    design=NHANES_subset_femur)

car::Anova(model_tscore_protein_poly, type="II", test="F")

model_tscore_proteingkg_ns <- svyglm(femurneck_t_score ~ ns(PTNgkg, 3) + CHO + LIP +
                                       RIAGENDR + BMXWT + RIDAGEYR
                                     + MET + days_gc + EI_EER_class + disease_cat, 
                                     design=NHANES_subset_femur)

car::Anova(model_tscore_proteingkg_ns, type="II", test="F")

plot(predict_response(model_tscore_protein_poly, "PTN [all]"),
     show_residuals = TRUE, show_residuals_line = TRUE)

plot(predict_response(model_tscore_proteingkg_ns, "PTNgkg [all]"),
     show_residuals = TRUE, show_residuals_line = TRUE)

t_model_1 <- car::Anova(model_tscore_protein_poly, type="II")
t_model_2 <- car::Anova(model_tscore_proteingkg_ns, type="II")

txt_neck_t1 <- glue("Polynomial term p-value {p_format_g(t_model_1$`Pr(>Chisq)`[1])}")
txt_neck_t2 <- glue("Spline term p-value = {p_format_g(t_model_2$`Pr(>Chisq)`[1])}")

pred_tscore_1 <- predictions(model_tscore_protein_poly,
                             by=c("PTN"), newdata=pred_data_5,
                             wts=get_weights(NHANES_subset_femur, "PTN"))

pred_tscore_2 <- predictions(model_tscore_proteingkg_ns,
                             by=c("PTNgkg"), newdata=pred_data_6,
                             wts=get_weights(NHANES_subset_femur, "PTNgkg"))

(
  tscore_p1 <- ggplot(pred_tscore_1, aes(x=PTN, y=estimate)) +
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
    geom_line(color="blue2") +
    labs(x="Protein intake (g/day)", y="Femoral neck T-score") +
    annotate(geom="text", x=Inf, y=-Inf, label=txt_neck_t1, hjust=1, vjust=-0.2,
             size=2.5) +
    scale_x_continuous(limits=c(25, 250)) +
    scale_y_continuous(limits=c(-2.604, 0.5), breaks=seq(-2.5, 1, .5)) +
    guides(size='none') +
    theme_avp()
)

(
  tscore_p2 <- ggplot(pred_tscore_2, aes(x=PTNgkg, y=estimate)) +
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
    geom_line(color="blue2") +
    labs(x="Protein intake (g/kgBM/day)", y="Femoral neck T-score") +
    annotate(geom="text", x=Inf, y=-Inf, label=txt_neck_t2, hjust=1, vjust=-0.2,
             size=2.5) +
    scale_x_continuous(limits=c(0.25, 2.5), breaks=seq(0.5, 2.5, .5)) +
    scale_y_continuous(limits=c(-2.604, 0.5), breaks=seq(-2.5, 1, .5)) +
    guides(size='none') +
    theme_avp()
)

tscore_p1 + tscore_p2 + plot_annotation(tag_levels = 'A')

#ggsave(here("figures/femoralneck_tscore_marginals.png"), units="in",
#       width=9, height=5, dpi=600)

## predictions at specific points 

## filter dataframes for specific points

target_values_ptn <- c(30, 60, 90, 120, 150, 187.5)
target_values_ptngkg <- c(0.4, .8, 1.2, 1.6, 2, 2.5)

filter_closest <- function(df, column, outcome="bone") {
  if (column == "PTN") {
    target_values <- target_values_ptn
  } else {
    target_values <- target_values_ptngkg
  }
  
  # Function to find the closest match for each target value
  find_closest <- function(target) {
    df %>%
      dplyr::filter(abs(.data[[column]] - target) == min(abs(.data[[column]] - target))) %>%
      slice(1) # Pick the first one in case of ties
  }
  
  # Apply find_closest to each target value and combine results
  df <- do.call(rbind, lapply(target_values, find_closest)) 
  
  if (column == "PTN") {
    df <- df |> 
      select(PTN, estimate, conf.low, conf.high) |> as_tibble() 
  } else {
    df <- df |> 
      select(PTNgkg, estimate, conf.low, conf.high) |> as_tibble() 
  }
  
  if (outcome == "lean mass") {
    df <- df |> mutate(label = paste0(round((estimate)/1000, 1), " kg (95% CI ",
                        round((conf.low)/1000, 1), ", ",
                        round((conf.high)/1000, 1), ")"))
  } 
  
  if (outcome == "t-score") {
    df <- df |> mutate(label = paste0(round((estimate), 2), " (95% CI ",
                                      round((conf.low), 2), ", ",
                                      round((conf.high), 2), ")"))
    } else {
    
    df <- df |> mutate(label = paste0(round((estimate), 4), " gm/cm² (95% CI ",
                                      round((conf.low), 4), ", ",
                                      round((conf.high), 4), ")"))
    }
  
  if (column == "PTN") {
    df <- df |> mutate(predictor_level = target_values_ptn) |> 
      select(predictor_level, label, PTN, estimate, conf.low, conf.high) 
  } else {
    df <- df |> mutate(predictor_level = target_values_ptngkg)  |> 
      select(predictor_level, label, PTNgkg, estimate, conf.low, conf.high) 
    
  }

  
  return(df)
}

pred_table1 <- filter_closest(pred1, "PTN", "lean mass") |> mutate(row = "lean")
pred_table2 <- filter_closest(pred2, "PTNgkg", "lean mass") |> mutate(row = "lean")
pred_table3 <- filter_closest(pred3, "PTN") |> mutate(row = "wb")
pred_table4 <- filter_closest(pred4, "PTNgkg") |> mutate(row = "wb")
pred_table5 <- filter_closest(pred5_neck, "PTN") |> mutate(row = "fem")
pred_table6 <- filter_closest(pred6_neck, "PTNgkg") |> mutate(row = "fem")
pred_table5_t <- filter_closest(pred_tscore_1, "PTN", "t-score") |> mutate(row = "fem_t")
pred_table6_t <- filter_closest(pred_tscore_2, "PTNgkg", "t-score") |> mutate(row = "fem_t")
pred_table7 <- filter_closest(pred7, "PTN") |> mutate(row = "sp")
pred_table8 <- filter_closest(pred8, "PTNgkg") |> mutate(row = "sp")

pred_t1 <- reduce(list(pred_table1, pred_table3, 
                       pred_table5, pred_table5_t, pred_table7),
                  full_join) |> 
  select(predictor_level, label, row) |> 
  pivot_wider(names_from=predictor_level, values_from=label)

pred_t2 <- reduce(list(pred_table2, pred_table4, 
                       pred_table6, pred_table6_t, pred_table8),
                  full_join) |> 
  select(predictor_level, label, row) |> 
  pivot_wider(names_from=predictor_level, values_from=label)

pred_t1 |> writexl::write_xlsx("tables/predicted_values_models_t1.xlsx")
pred_t2 |> writexl::write_xlsx("tables/predicted_values_models_t2.xlsx")

## DAG #########################################################################

# this is just a visual representation of the dag within R.
# the code was generated by creating the dag manually in dagitty.net
# the article also features the plot generated by the website.

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

unadj_models <- model_summary_all_unadj |> 
  mutate(p.value = case_when(p.value <= 0.0001 ~ "< 0.0001",
                             p.value > 0.0001 & p.value < 0.001 ~ "< 0.001",
                             p.value > 0.001 ~ format(round(p.value, digits=3), nsmall=3)))

unadj_models <- unadj_models |> 
  mutate(term = as.factor(c("Lean mass ~ protein (g/d)",
                            "Whole-body BMD ~ protein (g/d)",
                            "Femoral neck BMD ~ protein (g/d)",
                            "Spine BMD ~ protein (g/d)"))) |> 
  mutate(value2 = stringr::str_replace(value2, ";", " to"))

txt1 <- unadj_models |>
  filter(term=="Lean mass ~ protein (g/d)")

txt2 <- unadj_models |>
  filter(term=="Whole-body BMD ~ protein (g/d)")

txt3 <- unadj_models |>
  filter(term=="Femoral neck BMD ~ protein (g/d)")  

txt4 <- unadj_models |>
  filter(term=="Spine BMD ~ protein (g/d)")

txt1_an <- glue("Coefficient: {txt1$estimate[1]} g lean mass per 50 g protein (95%CI {txt1$value2[1]}), p-value {txt1$p.value[1]}")

txt2_an <- glue("Coefficient: {txt2$estimate[1]} gm/cm² BMD per 50 g protein (95%CI {txt2$value2[1]}), p-value {txt2$p.value[1]}")

txt3_an <- glue("Coefficient: {txt3$estimate[1]} gm/cm² BMD per 50 g protein (95%CI {txt3$value2[1]}), p-value {txt3$p.value[1]}")

txt4_an <- glue("Coefficient: {txt4$estimate[1]} gm/cm² BMD per 50 g protein (95%CI {txt4$value2[1]}), p-value {txt4$p.value[1]}")

# unadjusted models

model_leanmass_protein_un <- svyglm(DXDTOLE ~ PTN, 
                                     design=NHANES_subset_lm) 

model_wbBMD_protein_un <- svyglm(DXDTOBMD ~ PTN,
                                  design=NHANES_subset_wb) 

model_femurneckBMD_protein_un <- svyglm(DXXNKBMD ~ PTN, 
                                     design=NHANES_subset_femur) 

model_spineBMD_protein_un <- svyglm(DXXOSBMD ~ PTN, 
                                     design=NHANES_subset_spine) 

# unadjusted predictions

un_pred1 <- predictions(model_leanmass_protein_un, by=c("PTN"),
                        wts=NHANES_subset_lm$variables$weights_all)
un_pred2 <- predictions(model_wbBMD_protein_un, by=c("PTN"),
                        wts=NHANES_subset_wb$variables$weights_all)
un_pred3 <- predictions(model_femurneckBMD_protein_un, by=c("PTN"),
                        wts=NHANES_subset_femur$variables$weights_all)
un_pred4 <- predictions(model_spineBMD_protein_un, by=c("PTN"),
                        wts=NHANES_subset_spine$variables$weights_all)

# unadjusted plots

un_p1 <- ggplot(un_pred1, aes(x=PTN, y=estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  labs(x="Protein intake (g)", y="Lean mass (g)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt1_an, hjust=1, vjust=-0.2,
           size=3) +
  guides(size='none') +
  theme_avp()

un_p2 <- ggplot(un_pred2, aes(x=PTN, y=estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  labs(x="Protein intake (g)", y="Whole-body BMD (gm/cm²)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt2_an, hjust=1, vjust=-0.2,
           size=3) +
  guides(size='none') +
  theme_avp()

un_p3 <- ggplot(un_pred3, aes(x=PTN, y=estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  labs(x="Protein intake (g)", y="Femoral neck BMD (gm/cm²)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt3_an, hjust=1, vjust=-0.2,
           size=3) +
  guides(size='none') +
  theme_avp()

un_p4 <- ggplot(un_pred4, aes(x=PTN, y=estimate)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.2, fill="blue1") +
  geom_line(color="blue2") +
  labs(x="Protein intake (g)", y="Spine BMD (gm/cm²)") +
  annotate(geom="text", x=Inf, y=-Inf, label=txt4_an, hjust=1, vjust=-0.2,
           size=3) +
  guides(size='none') +
  theme_avp()

(un_p1 + un_p2 + un_p3 + un_p4) &
  plot_annotation(tag_levels = "A")
 
#ggsave(here("figures/scatter_unadjusted.png"), units="in", width=14, height=10, dpi=600, type="cairo")

################################################################################