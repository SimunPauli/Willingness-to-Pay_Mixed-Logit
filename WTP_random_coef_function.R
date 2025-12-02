library(MASS)
library(doParallel)
library(foreach)

load("O:/Public/4233-110918-FLOW-persondata/Temp/model_RE_full_inc_prichar.RData")
estimate_csv <- read.csv2("O:/Public/4233-110918-FLOW-persondata/WTP/Filer til WTP beregning/BasisModel 12-11/Flow_RE_full_estimates",
                          sep=";")
estimate_vec <- setNames(as.numeric(estimate_csv$model.estimate), estimate_csv$X)
robvarcov_mat <- read.csv2("O:/Public/4233-110918-FLOW-persondata/WTP/Filer til WTP beregning/BasisModel 12-11/Flow_RE_full_robcovvar",
                           sep=";", row.names = 1) |> as.matrix()



R = 200
K = 50

alpha_name <- c("mu_sma", "sigma_sma", "asc_sma_not_dk", "asc_prichar_sma") #all coef on the numerator of WTP (including random, correlation and interactions coefficients)
beta_name <- "b_pprice"

database$not_dk <- database$Spain + database$Czech + database$Italy + database$Germany + database$Ireland
interaction_alpha <- list("asc_sma_not_dk" = "not_dk",
                          "asc_prichar_sma" = "HomeChargAvail")
random_coef_alpha <- list("mu_sma" = c("sigma_sma")) # "mean" = c("sandard deviation coef", "correlation coef"), leave out correlation if not relvant



function(alpha_name, beta_name,
         random_coef_alpha, interaction_alpha,
         estimated_coef = model$estimate[-which(names(model$estimate) %in% apollo_fixed)],
         estimated_Sigma = model$robvarcov
) {
  
  draws <- MASS::mvrnorm(n = R,
                         mu = estimated_coef, 
                         Sigma = estimated_Sigma)
  alpha_beta_drawn <- draws[,which(colnames(draws) %in% c(alpha_name,beta_name))]
  
  cl <- makeCluster(detectCores() - 4)
  if (!is.null(random_coef_alpha)) {#if there are any random coefficients
    alpha_beta_drawn_boot <- matrix(NA, nrow = K * R, ncol = length(random_coef_alpha))
    colnames(alpha_beta_drawn_boot) <- names(random_coef_alpha)
    
    # Set up parallel backend
    registerDoParallel(cl)
    
    results <- foreach(k = 1:R, .combine = rbind) %dopar% {
      temp_matrix <- matrix(NA, nrow = K, ncol = length(random_coef_alpha))
      for (i in 1:length(random_coef_alpha)) {# if there are multiple random coefficients, loop one at a time
        mean_rand_coef <- alpha_beta_drawn[k, which(colnames(alpha_beta_drawn) == names(random_coef_alpha)[i])]
        sigma_rand_coef <- alpha_beta_drawn[k, which(colnames(alpha_beta_drawn) == random_coef_alpha[[i]][1])]
        if (length(random_coef_alpha[[i]]) != 1) { # checking for correlation coefficent (rho)
          sigma_rand_coef <- sigma_rand_coef + alpha_beta_drawn[,which(colnames(alpha_beta_drawn) == random_coef_alpha[[i]][2])]
        }
        temp_matrix[, i] <- rnorm(n = K, mean_rand_coef, sigma_rand_coef^2)
      }
      temp_matrix
    }
    
    stopCluster(cl)
    alpha_beta_drawn_boot <- as.data.frame(results)
    names(alpha_beta_drawn_boot) <- paste0(names(random_coef_alpha), "_random_coef_draws")
    
  }
  alpha_beta_drawn_rep <- alpha_beta_drawn[rep(seq_len(R), each = K), , drop = FALSE]
  alpha_beta_drawn_rep <- cbind(alpha_beta_drawn_rep, alpha_beta_drawn_boot)
  
  
  
  unregister_dopar <- function() {
    env <- foreach:::.foreachGlobals
    rm(list=ls(name=env), pos=env)
  }
  unregister_dopar()
  
  
  
  # Convert to matrix for faster access
  alpha_beta_drawn_rep <- as.matrix(alpha_beta_drawn_rep)
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
  result <- foreach(i = 1:nrow(database_indiv_inter), .combine = cbind) %dopar% {
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
  
  wtp_dist <- as.numeric(rowSums(result))
}
quantile(wtp_dist, c(0.1,0.9))



