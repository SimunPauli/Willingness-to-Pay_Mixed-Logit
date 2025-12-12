source("O:/Public/4233-110918-FLOW-persondata/WTP/Willingness-to-Pay_Mixed-Logit/WTP_random_coef_function.R")



#load("O:/Public/4233-110918-FLOW-persondata/Temp/model_RE_full_inc_prichar.RData")
estimated_csv <- read.csv2("O:/Public/4233-110918-FLOW-persondata/WTP/Filer til WTP beregning/modeller 27_11_2025/Flow_RE_full_v2/Flow_RE_full_v2_estimates",
                           sep=";")
estimated_vec <- setNames(as.numeric(estimated_csv$model.estimate), estimated_csv$X)[-which(estimated_csv[,2] == 0)]
robvarcov_mat <- read.csv2("O:/Public/4233-110918-FLOW-persondata/WTP/Filer til WTP beregning/modeller 27_11_2025/Flow_RE_full_v2/Flow_RE_full_v2_robcovvar",
                           sep=";", row.names = 1) |> as.matrix()

library(openxlsx)
database <- read.csv("O:/Public/4233-110918-FLOW-persondata/WTP/Filer til WTP beregning/Flow_SP1_data_v1.csv")
### Use only SP data with response time higher than 5.4 minutes
database = subset(database,database$lowtime==0)

### Create new variables 
database$HomeChargAvail = as.numeric(database$HomeCharge < 3 )
database$Income_new = pmax(database$income_val,0)
summary(database$Income_new)
database$Income_std = (database$Income_new - 3430)/10000
summary(database$Income_std)
summary(database$RanBat_PHEV1)
summary(database$RanBat_PHEV2)

database$MissInc = (database$Income_new == 0)
database$inc_norm=(database$Income_new/3430)+database$MissInc
summary(database$inc_norm)
database$av_icv1 = 1
database$av_icv2 = 1
database$av_bev1 = 1
database$av_bev2 = 1
database$av_phev1 = 1
database$av_phev2 = 1

database$HomeChargAvail = (database$HomeCharge < 3 )

database$age1 = 1*(database$age < 31)
database$age2 = 1*(database$age > 30) * (database$age < 48)
database$age3 = 1*(database$age > 47) * (database$age < 66)
database$age4 = 1*(database$age > 65)
database$edu1 = 1*(database$edu < 4)
database$edu2 = 1*(database$edu == 4)
database$edu3 = 1*(database$edu > 4)

database$Income_std_cz_de <-  database$Income_std*(database$Germany + database$Czech)

wtp <- apollo_WTP(
  alpha_name = c("b_accel"),
  beta_name = c("b_pprice"), 
  estimated_coef = estimated_vec, 
  estimated_Sigma = robvarcov_mat, 
  database = database,
  R = 1000,
  K = NULL
)
round(wtp$wtp_mean*1000)
round(wtp$wtp_95_confidence*1000)


#wtp <- apollo_WTP(
#  alpha_name = c("mu_min", "mu_min_fem", "mu_min_age1", "mu_min_it", "mu_min_prichar"), #all parameters to be included in the numerator of WTP
#  random_alpha_normal = list("mu_min" = c("sigma_min")), #list of parameters with random coefficients (relevant for this WTP) and their respective random effects: list("beta1" = c("beta1_sigam","beta1_rho"),...)
#  interaction_alpha = list(
#    "mu_min_fem" = "Kvinde",
#    "mu_min_it" = "Italy",
#    "mu_min_age1" = "age1",
#    "mu_min_prichar" = "HomeChargAvail"
#  ), #list of parameters with interactions (relevant for this WTP) and their respective attribute: list("beta_inter1" = "inter1")
#  
#  beta_name = c("b_pprice","l_pprice", "b_pprice_miss", "b_pp_cz", "l_pprice_cz_de", "b_pp_ei"), #all parameters to be included in the denominator of WTP
#  random_beta_uniform = list(
#    "b_pprice" = c("sigma_pp")
#  ),
#  interaction_beta = list(
#    "l_pprice" = "Income_std",
#    "b_pprice_miss" = "MissInc",
#    "b_pp_cz" = "Czech",
#    "l_pprice_cz_de" = "Income_std_cz_de",
#    "b_pp_ei" = "Ireland"
#  ),
#  
#  draw_transformation = list( #if draws need be transformed. E.g. log-uniform: list("exp" = "sigma_uniform)
#    "-" = c("b_pprice","l_pprice", "b_pprice_miss", "b_pp_cz", "l_pprice_cz_de", "b_pp_ei"),
#    "exp" = c("b_pprice","l_pprice", "b_pprice_miss", "b_pp_cz", "l_pprice_cz_de", "b_pp_ei")
#  ),
#  
#  estimated_coef = estimated_vec, #all coefficients in the model. Should not included coefficient set to 0
#  estimated_Sigma = robvarcov_mat, #(robust) covariance matrix
#  database = database,
#  R = 500,
#  K = 10
#)




