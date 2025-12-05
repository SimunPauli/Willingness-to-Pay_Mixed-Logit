


#load("O:/Public/4233-110918-FLOW-persondata/Temp/model_RE_full_inc_prichar.RData")
estimated_csv <- read.csv2("O:/Public/4233-110918-FLOW-persondata/WTP/Filer til WTP beregning/modeller 27_11_2025/Flow_RE_full_se_v3/Flow_RE_full_se_v3_estimates",
                          sep=";")
estimated_vec <- setNames(as.numeric(estimated_csv$model.estimate), estimated_csv$X)[-which(estimated_csv[,2] == 0)]
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
apollo_WTP <- function(alpha_name, #all parameters to be included in the numerator of WTP
                beta_name, #all parameters to be included in the denominator of WTP
                random_alpha_normal = NULL, #list of parameters with random coefficients (relevant for this WTP) and their respective random effects: list("beta1" = c("beta1_sigam","beta1_rho"),...)
                interaction_alpha = NULL, #list of parameters with interactions (relevant for this WTP) and their respective attribute: list("beta_inter1" = "inter1")
                estimated_coef = model$estimate[-which(names(model$estimate) %in% apollo_fixed)], #all coefficients in the model. Should not included coefficient set to 0
                estimated_Sigma = model$robvarcov, #(robust) covariance matrix
                database = database,
                R = NULL,
                K = NULL
) {
  
  if (!is.numeric(estimated_coef) || is.null(names(estimated_coef))) stop("estimated_coef must be a named numeric vector.")
  if (!is.matrix(estimated_Sigma) || any(colnames(estimated_Sigma) != rownames(estimated_Sigma))) stop("estimated_Sigma must be a square matrix with matching row/col names.")
  if (!all(names(estimated_coef) == colnames(estimated_Sigma))) {
    # try to reorder Sigma to match coef vector
    if (all(names(estimated_coef) %in% colnames(estimated_Sigma))) {
      estimated_Sigma <- estimated_Sigma[names(estimated_coef), names(estimated_coef)]
    } else stop("Names of estimated_coef must match column/row names of estimated_Sigma.")
  }
  if (!all(unlist(interaction_alpha) %in% colnames(database))) {
    stop("Some interaction variable names from interaction_alpha are not columns of database.")
  }
  if (!all(unlist(interaction_alpha) %in% interaction_alpha)) {} # noop - kept simple
  
  # ensure alpha_name and beta_name are present in estimated_coef
  if (!all(alpha_name %in% names(estimated_coef))) stop("Some alpha_name entries not in estimated_coef.")
  if (!all(beta_name %in% names(estimated_coef))) stop("Some beta_name entries not in estimated_coef.")
  
  # random_alpha_normal validation
  if (!is.null(random_alpha_normal)) {
    if (!is.list(random_alpha_normal) || is.null(names(random_alpha_normal))) stop("random_alpha_normal must be a named list.")
    if (!all(names(random_alpha_normal) %in% alpha_name)) stop("Names of random_alpha_normal must be a subset of alpha_name.")
    # check sigma names exist
    sigma_names <- unlist(random_alpha_normal)
    if (!all(sigma_names %in% names(estimated_coef))) stop("Some sigma parameter names referenced in random_alpha_normal are not in estimated_coef.")
  }
  
  #### 2. Draw R multivariate normal draws for fixed (mean) coefficients ####
  draws <- MASS::mvrnorm(n = R,
                         mu = estimated_coef, 
                         Sigma = estimated_Sigma)
  alpha_beta_drawn <- draws[,match(c(alpha_name, beta_name),colnames(draws))]
  
  if (!is.null(random_alpha_normal)) {#if there are any random coefficients
    alpha_beta_drawn_boot <- matrix(NA, nrow = K*R, ncol = length(random_alpha_normal)) #Matrix to store matrix during paralell runs
    j <- 1
    for(k in 1:R) {
      for (i in 1:length(random_alpha_normal)) {# if there are multiple random coefficients, loop one at a time
        mean_rand_coef <- alpha_beta_drawn[k, which(colnames(alpha_beta_drawn) == names(random_alpha_normal)[i])] #
        sigma_rand_coef <- alpha_beta_drawn[k, which(colnames(alpha_beta_drawn) == random_alpha_normal[[i]][1])]
        #if (length(random_alpha_normal[[i]]) != 1) { # checking for correlation coefficent (rho)
        #  sigma_rand_coef <- sigma_rand_coef + alpha_beta_drawn[k,which(colnames(alpha_beta_drawn) == random_alpha_normal[[i]][2])]
        #}
        alpha_beta_drawn_boot[j:(j+K-1), i] <- rnorm(n = K, mean_rand_coef, abs(sigma_rand_coef)) #rnrom(n_draws, mean, sd)
      }
      j <- j + K
    }
    
    colnames(alpha_beta_drawn_boot) <- paste0(names(random_alpha_normal), "_random_coef_draws")
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
        c(names(random_alpha_normal), unlist(random_alpha_normal)), 
        colnames(alpha_beta_drawn_rep)
        ), #These are replaced by _random_coef_draws$
      interaction_idx,
      beta_idx
    )
  )
  
  
  wtp_matrix <- matrix(NA, nrow = nrow(database_indiv_inter), ncol = R*K)
  alpha_matrix <- matrix(NA, nrow = nrow(database_indiv_inter), ncol = R*K)
  beta_matrix <- matrix(NA, nrow = nrow(database_indiv_inter), ncol = R*K)
  for (i in 1:nrow(database_indiv_inter)) {
    x_i <- as.numeric(database_indiv_inter[i, ])
    wtp_i <- numeric(nrow(alpha_beta_drawn_rep))
    alpha_i <- numeric(nrow(alpha_beta_drawn_rep))
    beta_i <- numeric(nrow(alpha_beta_drawn_rep))
    
    for (j in 1:(R*K)) {
      alpha <- 
        sum(alpha_beta_drawn_rep[j, random_coef_idx]) +
        sum(alpha_beta_drawn_rep[j, interaction_idx] * x_i) + #if interaction_idx is integer(0) this equal 0: sum(numeric(0)) -> 0
        sum(alpha_beta_drawn_rep[j, non_interaction_idx]) #if non_interaction_idx is integer(0) this equal 0: sum(numeric(0)) -> 0
      
      beta <- alpha_beta_drawn_rep[j, beta_idx]
      
      wtp_i[j] <- - alpha / beta # negative
      alpha_i[j] <- alpha
      beta_i[j] <- beta
    }
    wtp_matrix[i, ] <- wtp_i #each row one individual
    alpha_matrix[i, ] <- alpha
    beta_matrix[i, ] <- beta
  }
  
  wtp_dist <- as.numeric(colMeans(wtp_matrix)) #means of wtp across simulation per indivudal
  alpha_dist <- as.numeric(colMeans(alpha_matrix))
  beta_dist <- as.numeric(colMeans(beta_matrix))
  list(
    wtp_dist = wtp_dist,
    alpha_dist = alpha_dist,
    beta_dist = beta_dist,
    wtp_mean = mean(wtp_dist),
    wtp_90_confidence = quantile(wtp_dist, c(0.05,0.95)) #90% confidence?
  )
}
 


wtp <- apollo_WTP(alpha_name = c("mu_min", "sigma_min", "mu_min_fem", "mu_min_age1", "mu_min_it", "mu_min_prichar"), #all parameters to be included in the numerator of WTP
    beta_name = "b_pprice", #all parameters to be included in the denominator of WTP
    random_alpha_normal = list("mu_min" = c("sigma_min")), #list of parameters with random coefficients (relevant for this WTP) and their respective random effects: list("beta1" = c("beta1_sigam","beta1_rho"),...)
    interaction_alpha = list(
      "mu_min_fem" = "Kvinde",
      "mu_min_age1" = "age1",
      "mu_min_it" = "Italy",
      "mu_min_prichar" = "HomeChargAvail"
      ), #list of parameters with interactions (relevant for this WTP) and their respective attribute: list("beta_inter1" = "inter1")
    estimated_coef = estimated_vec, #all coefficients in the model. Should not included coefficient set to 0
    estimated_Sigma = robvarcov_mat, #(robust) covariance matrix
    database = database,
    R = 500,
    K = 10
)

wtp$wtp_mean
hist(wtp$wtp_dist, breaks = 100)

kvinde_share <- mean(database$Kvinde == 1)
kvinde_est <- estimated_csv[which(estimated_csv[,1]=="mu_min_fem"),2]
it_share <- mean(database$Italy == 1)
it_est <- estimated_csv[which(estimated_csv[,1]=="mu_min_it"),2]
age1_share <- mean(database$age1 == 1)
age1_est <- estimated_csv[which(estimated_csv[,1]=="mu_min_age1"),2]
prichar_share <- mean(database$HomeChargAvail == 1)
prichar_est <- estimated_csv[which(estimated_csv[,1]=="mu_min_prichar"),2]
-(-1.301293269 + kvinde_share*kvinde_est + it_share*it_est + age1_share*age1_est + prichar_share*prichar_est )/(-0.501882423)


