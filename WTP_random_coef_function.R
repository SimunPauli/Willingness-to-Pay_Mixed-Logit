


#load("O:/Public/4233-110918-FLOW-persondata/Temp/model_RE_full_inc_prichar.RData")
estimated_csv <- read.csv2("O:/Public/4233-110918-FLOW-persondata/WTP/Filer til WTP beregning/modeller 27_11_2025/Flow_RE_full_se_v3/Flow_RE_full_se_v3_estimates",
                          sep=";")
estimated_vec <- setNames(as.numeric(estimate_csv$model.estimate), estimate_csv$X)
robvarcov_mat <- read.csv2("O:/Public/4233-110918-FLOW-persondata/WTP/Filer til WTP beregning/modeller 27_11_2025/Flow_RE_full_se_v3/Flow_RE_full_se_v3_robcovvar",
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
library(MASS)
library(doParallel)
library(foreach)

#############################################
WTP <- function(alpha_name, #all parameters to be included in the numerator of WTP
                beta_name, #all parameters to be included in the denominator of WTP
                random_coef_alpha = NULL, #list of parameters with random coefficients (relevant for this WTP) and their respective random effects: list("beta1" = c("beta1_sigam","beta1_rho"),...)
                interaction_alpha = NULL, #list of parameters with interactions (relevant for this WTP) and their respective attribute: list("beta_inter1" = "inter1")
                estimated_coef = model$estimate[-which(names(model$estimate) %in% apollo_fixed)], #all coefficients in the model. Should not included coefficient set to 0
                estimated_Sigma = model$robvarcov, #(robust) covariance matrix
                database = database,
                R = NULL,
                K = NULL
) {
  if(length(unlist(interaction_alpha)) != sum(colnames(database) %in% as.character(unlist(interaction_alpha)))) {
    stop(paste0("List of attributes/variables from interaction_alpha doesn't match column names in database"))
  }
  if(any(!(random_coef_alpha %in% alpha_name)) & any(!(names(random_coef_alpha) %in% alpha_name))) {
    stop(paste0("Element in random_coef_alpha not in alpha_name"))
  }
  
  draws <- MASS::mvrnorm(n = R,
                         mu = estimated_coef, 
                         Sigma = estimated_Sigma)
  alpha_beta_drawn <- draws[,match(c(alpha_name, beta_name),colnames(draws))]
  
  if (!is.null(random_coef_alpha)) {#if there are any random coefficients
    alpha_beta_drawn_boot <- matrix(NA, nrow = K*R, ncol = length(random_coef_alpha)) #Matrix to store matrix during paralell runs
    j <- 1
    for(k in 1:R) {
      for (i in 1:length(random_coef_alpha)) {# if there are multiple random coefficients, loop one at a time
        mean_rand_coef <- alpha_beta_drawn[k, which(colnames(alpha_beta_drawn) == names(random_coef_alpha)[i])] #
        sigma_rand_coef <- alpha_beta_drawn[k, which(colnames(alpha_beta_drawn) == random_coef_alpha[[i]][1])]
        if (length(random_coef_alpha[[i]]) != 1) { # checking for correlation coefficent (rho)
          sigma_rand_coef <- sigma_rand_coef + alpha_beta_drawn[k,which(colnames(alpha_beta_drawn) == random_coef_alpha[[i]][2])]
        }
        alpha_beta_drawn_boot[j:(j+K-1), i] <- rnorm(n = K, mean_rand_coef, abs(sigma_rand_coef)) #rnrom(n_draws, mean, sd)
      }
      j <- j + K
    }
    
    colnames(alpha_beta_drawn_boot) <- paste0(names(random_coef_alpha), "_random_coef_draws")
  }
  alpha_beta_drawn_rep <- alpha_beta_drawn[rep(seq_len(R), each = K), , drop = FALSE]
  alpha_beta_drawn_rep <- cbind(alpha_beta_drawn_rep, alpha_beta_drawn_boot)
  
  # Convert to matrix for faster access
  database_indiv_inter <- database[ #database on individual level and 
    !duplicated(database$Respondent_Serial), 
    match(as.character(unlist(interaction_alpha)), colnames(database))
  ] |>
    as.matrix()
  
  # Precompute column indices
  beta_idx <- match(beta_name, colnames(alpha_beta_drawn_rep))
  
  interaction_cols <- names(interaction_alpha)
  interaction_idx <- match(interaction_cols, colnames(alpha_beta_drawn_rep))
  random_coef_idx <- grep("_random_coef_draws$", colnames(alpha_beta_drawn_rep))
  non_interaction_idx <- setdiff(
    seq_len(ncol(alpha_beta_drawn_rep)),
    c( #idx in this vector not to be included in non_interaction_idx
      random_coef_idx,
      match(
        c(names(random_coef_alpha), unlist(random_coef_alpha)), 
        colnames(alpha_beta_drawn_rep)
        ), #These are replaced by _random_coef_draws$
      interaction_idx,
      beta_idx
    )
  )
  
  
  wtp_matrix <- matrix(NA, nrow = nrow(database_indiv_inter), ncol = R*K)
  for (i in 1:nrow(database_indiv_inter)) {
    x_i <- as.numeric(database_indiv_inter[i, ])
    wtp_i <- numeric(nrow(alpha_beta_drawn_rep))
    
    for (j in 1:(R*K)) {
      alpha <- 
        sum(alpha_beta_drawn_rep[j, random_coef_idx]) +
        sum(alpha_beta_drawn_rep[j, interaction_idx] * x_i) + #if interaction_idx is integer(0) this equal 0: sum(numeric(0)) -> 0
        sum(alpha_beta_drawn_rep[j, non_interaction_idx]) #if non_interaction_idx is integer(0) this equal 0: sum(numeric(0)) -> 0
      
      beta <- alpha_beta_drawn_rep[j, beta_idx]
      
      wtp_i[j] <- alpha / beta
    }
    
    wtp_matrix[i, ] <- wtp_i
  }
  
  wtp_dist_pos <- as.numeric(rowMeans(wtp_matrix)) #means of wtp across simulation per indivudal
  wtp_dist <- -wtp_dist_pos
  list(
    wtp_dist = wtp_dist,
    wtp_mean = mean(wtp_dist),
    wtp_90_confidence = quantile(wtp_dist, c(0.05,0.95)) #90% confidence?
  )
}
 


WTP(alpha_name = c("mu_min", "sigma_min", "mu_min_fem", "mu_min_age1", "mu_min_it", "mu_min_prichar"), #all parameters to be included in the numerator of WTP
    beta_name = "b_pprice", #all parameters to be included in the denominator of WTP
    random_coef_alpha = list("mu_min" = c("sigma_min")), #list of parameters with random coefficients (relevant for this WTP) and their respective random effects: list("beta1" = c("beta1_sigam","beta1_rho"),...)
    interaction_alpha = list(
      "mu_min_fem" = "Kvinde",
      "mu_min_age1" = "age",
      "mu_min_it" = "Italy",
      "mu_min_prichar" = "HomeCharge"
      ), #list of parameters with interactions (relevant for this WTP) and their respective attribute: list("beta_inter1" = "inter1")
    estimated_coef = estimate_csv[-which(estimate_csv[,2] == 0),2], #all coefficients in the model. Should not included coefficient set to 0
    estimated_Sigma = robvarcov_mat, #(robust) covariance matrix
    database = database,
    R = 100,
    K = 10
)















