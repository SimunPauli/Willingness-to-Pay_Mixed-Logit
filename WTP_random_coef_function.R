


#load("O:/Public/4233-110918-FLOW-persondata/Temp/model_RE_full_inc_prichar.RData")
estimate_csv <- read.csv2("O:/Public/4233-110918-FLOW-persondata/WTP/Filer til WTP beregning/BasisModel 12-11/Flow_RE_full_estimates",
                          sep=";")
estimate_vec <- setNames(as.numeric(estimate_csv$model.estimate), estimate_csv$X)
robvarcov_mat <- read.csv2("O:/Public/4233-110918-FLOW-persondata/WTP/Filer til WTP beregning/BasisModel 12-11/Flow_RE_full_robcovvar",
                           sep=";", row.names = 1) |> as.matrix()

library(openxlsx)
database <- read.csv("O:/Public/4233-110918-FLOW-persondata/WTP/Filer til WTP beregning/Flow_SP1_data_v1.csv")

library(MASS)
library(doParallel)
library(foreach)
R <- 200
K <- 2

alpha_name <- c("mu_sma", "sigma_sma", "asc_sma_not_dk", "asc_prichar_sma") #all coef on the numerator of WTP (including random, correlation and interactions coefficients)
beta_name <- "b_pprice"

database$not_dk <- database$Spain + database$Czech + database$Italy + database$Germany + database$Ireland
interaction_alpha <- list("asc_sma_not_dk" = "not_dk",
                          "asc_prichar_sma" = "HomeCharge")
random_coef_alpha <- list("mu_sma" = c("sigma_sma")) # "mean" = c("sandard deviation coef", "correlation coef"), leave out correlation if not relevant

estimated_coef <- estimate_csv[-which(estimate_csv[,2] == 0),2]
estimated_Sigma <- robvarcov_mat
#############################################
WTP <- function(alpha_name, #all parameters to be included in the numerator of WTP
                beta_name, #all parameters to be included in the denominator of WTP
                random_coef_alpha = NULL, #list of parameters with random coefficients (relevant for this WTP) and their respective random effects: list("beta1" = c("beta1_sigam","beta1_rho"),...)
                interaction_alpha = NULL, #list of parameters with interactions (relevant for this WTP) and their respective attribute: list("beta_inter1" = "inter1")
                estimated_coef = model$estimate[-which(names(model$estimate) %in% apollo_fixed)], #all coefficients in the model. Should not included coefficient set to 0
                estimated_Sigma = model$robvarcov, #(robust) covariance matrix
                database = database,
                R,
                K,
) {
  
  draws <- MASS::mvrnorm(n = R,
                         mu = estimated_coef, 
                         Sigma = estimated_Sigma)
  alpha_beta_drawn <- draws[,which(colnames(draws) %in% c(alpha_name,beta_name))]
  
  if (!is.null(random_coef_alpha)) {#if there are any random coefficients
    # Set up parallel backend
    cl <- makeCluster(detectCores() - 1) #number of CPU cores used 
    registerDoParallel(cl)
    
    alpha_beta_drawn_boot <- foreach(k = 1:R, .combine = rbind) %dopar% {
      temp_matrix <- matrix(NA, nrow = K, ncol = length(random_coef_alpha)) #Matrix to store matrix during paralell runs
      for (i in 1:length(random_coef_alpha)) {# if there are multiple random coefficients, loop one at a time
        mean_rand_coef <- alpha_beta_drawn[k, which(colnames(alpha_beta_drawn) == names(random_coef_alpha)[i])] #
        sigma_rand_coef <- alpha_beta_drawn[k, which(colnames(alpha_beta_drawn) == random_coef_alpha[[i]][1])]
        if (length(random_coef_alpha[[i]]) != 1) { # checking for correlation coefficent (rho)
          sigma_rand_coef <- sigma_rand_coef + alpha_beta_drawn[k,which(colnames(alpha_beta_drawn) == random_coef_alpha[[i]][2])]
        }
        temp_matrix[, i] <- rnorm(n = K, mean_rand_coef, abs(sigma_rand_coef)) #rnrom(n_draws, mean, sd)
      }
      temp_matrix
    }
    
    stopCluster(cl) #stop parallel backend
    alpha_beta_drawn_boot <- as.matrix(alpha_beta_drawn_boot)
    names(alpha_beta_drawn_boot) <- paste0(names(random_coef_alpha), "_random_coef_draws")
  }
  alpha_beta_drawn_rep <- alpha_beta_drawn[rep(seq_len(R), each = K), , drop = FALSE]
  alpha_beta_drawn_rep <- cbind(alpha_beta_drawn_rep, alpha_beta_drawn_boot)
  
  # Convert to matrix for faster access
  database_indiv_inter <- database[ #database on individual level and 
    !duplicated(database$Respondent_Serial), 
    which(colnames(database) %in% as.character(unlist(interaction_alpha)))
  ] |>
    as.matrix()
  
  # Precompute column indices
  interaction_cols <- names(interaction_alpha)
  interaction_idx <- match(interaction_cols, colnames(alpha_beta_drawn_rep))
  random_coef_idx <- grep("_random_coef_draws$", colnames(alpha_beta_drawn_rep))
  non_interaction_idx <- setdiff(
    seq_len(ncol(alpha_beta_drawn_rep)),
    c(
      random_coef_idx,
      match(c(names(random_coef_alpha), unlist(random_coef_alpha)), colnames(alpha_beta_drawn_rep)),
      interaction_idx,
      match(beta_name, colnames(alpha_beta_drawn_rep))
    )
  )
  beta_idx <- match(beta_name, colnames(alpha_beta_drawn_rep))
  
  
  cl <- makeCluster(detectCores() - 4)
  registerDoParallel(cl)
  #Calculate alpha for each individual (considering interactions terms)
  result <- foreach(i = 1:nrow(database_indiv_inter), .combine = rbind) %dopar% {
    x_i <- as.numeric(database_indiv_inter[i, ])
    wtp_i <- matrix(NA, nrow = nrow(alpha_beta_drawn_rep), ncol = 1)
    
    for (j in seq_len(nrow(alpha_beta_drawn_rep))) {
      alpha <- sum(alpha_beta_drawn_rep[j, interaction_idx] * x_i) +
        sum(alpha_beta_drawn_rep[j, random_coef_idx]) +
        sum(alpha_beta_drawn_rep[j, non_interaction_idx])
      beta <- alpha_beta_drawn_rep[j, beta_idx]
      wtp_i[j] <- alpha / beta
    }
    
    wtp_i
  }
  stopCluster(cl)
  wtp_dist_pos <- as.numeric(rowSums(result))
  wtp_dist <- wtp_dist_pos
  list(
    wtp_dist = wtp_dist,
    wtp_mean = mean(wtp_dist),
    wtp_90_confidence = quantile(wtp_dist, c(0.05,0.95)) #90% confidence?
  )
}
 



