library(MASS)
library(doParallel)
library(foreach)
 
#############################################
apollo_WTP <- function(alpha_name, #all parameters to be included in the numerator of WTP
                random_alpha_normal = NULL, #list of parameters with random coefficients (relevant for this WTP) and their respective random effects: list("beta1" = c("beta1_sigam","beta1_rho"),...)
                random_alpha_uniform = NULL,
                interaction_alpha = NULL, #list of parameters with interactions (relevant for this WTP) and their respective attribute: list("beta_inter1" = "inter1")
                
                beta_name, #all parameters to be included in the denominator of WTP
                random_beta_normal = NULL,
                random_beta_uniform  = NULL,
                interaction_beta = NULL,
                
                trans_param = NULL, #transform individual paramater after drawing. list("function1" = c("param_to_transform1"), "function2" = c("param_to_Transform2))
                trans_alpha = NULL, #transform sum of alpha. character vector of function c("function","function2"). function1() will be applied before function2()
                trans_beta = NULL,  #transform sum of beta character vector of function c("function","function2"). function1() will be applied before function2()
                
                estimated_coef = model$estimate[-which(names(model$estimate) %in% apollo_fixed)], #all coefficients in the model. Should not included coefficient set to 0
                estimated_Sigma = model$robvarcov, #(robust) covariance matrix
                database = database,
                R,
                K = NULL
) {
  
  if (!is.numeric(estimated_coef) || is.null(names(estimated_coef))) stop("estimated_coef must be a named numeric vector.")
  if (!is.matrix(estimated_Sigma) || any(colnames(estimated_Sigma) != rownames(estimated_Sigma))) stop("estimated_Sigma must be a square matrix with matching row/col names.")
  if (!all(names(estimated_coef) %in% colnames(estimated_Sigma))) stop("Names of estimated_coef must match column/row names of estimated_Sigma.")
  if (!all(unlist(interaction_alpha) %in% colnames(database))) stop("Some interaction variable names from interaction_alpha are not columns of database.")
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
  if(!is.numeric(K) & ((!is.null(random_alpha_normal)) | (!is.null(random_beta_normal))) ) stop("K needs to be numeric if WTP included random coefficient.")
  if(length(names(interaction_alpha)) != length(unique(names(interaction_alpha)))) stop("Parameters in interaction_alpha should be unique.")
  if(length(names(interaction_beta)) != length(unique(names(interaction_beta)))) stop("Parameters in interaction_beta should be unique.")
  
  #### 2. Draw R multivariate normal draws for fixed (mean) coefficients ####
  n_draws <- R
  draws <- MASS::mvrnorm(n = n_draws,
                         mu = estimated_coef, 
                         Sigma = estimated_Sigma)
  
  if(!is.null(K)) draws_exp <- draws[rep(seq_len(R), each = K), , drop = FALSE] # (R*K) x p. draws_exp only relevant if K exist 
  
  #extra draws for NORMAL random coefficient for ALPHA
  if (!is.null(random_alpha_normal)) {
    n_draws <- R * K
    # for each named random alpha, replace the mean value in draws_exp by rnorm using sd = sigma parameter
    for (coef_name in names(random_alpha_normal)) {
      sigma_name <- random_alpha_normal[[coef_name]][1] # first element is sigma name
      if (!sigma_name %in% colnames(draws_exp)) stop(paste0("sigma '", sigma_name, "' not found in estimated coefficients/draws."))
      # per-row sd should be the repeated sigma draws from draws_exp
      sigma_vec <- draws_exp[, sigma_name]
      mean_vec <- draws_exp[, coef_name]
      if (length(random_alpha_normal[[coef_name]]) == 2) { # checking for correlation coefficent (rho)
        rho_name <- random_alpha_normal[[coef_name]][2]
        sigma_vec <- sqrt(sigma_vec^2 + (draws_exp[, rho_name])^2)
      }
      # sample individual-level random coefficients
      new_vals <- rnorm(n = n_draws, mean = mean_vec, sd = sigma_vec)
      draws_exp[, coef_name] <- new_vals
    }
    draws <- draws_exp
  }
  
  #extra draws for NORMAL random coefficient for BETA
  if (!is.null(random_beta_normal)) {
    n_draws <- R * K
    for (coef_name in names(random_beta_normal)) {
      sigma_name <- random_beta_normal[[coef_name]][1] # first element is sigma name
      if (!sigma_name %in% colnames(draws_exp)) stop(paste0("sigma '", sigma_name, "' not found in estimated coefficients/draws."))
      sigma_vec <- draws_exp[, sigma_name]
      mean_vec <- draws_exp[, coef_name]
      if (length(random_beta_normal[[coef_name]]) == 2) { # checking for correlation coefficent (rho)
        rho_name <- random_beta_normal[[coef_name]][2]
        sigma_vec <- sqrt(sigma_vec^2 + (draws_exp[, rho_name])^2)
      }
      new_vals <- rnorm(n = n_draws, mean = mean_vec, sd = sigma_vec)
      draws_exp[, coef_name] <- new_vals
    }
    draws <- draws_exp
  }
  
  #extra draws for UNIFORM random coefficient for ALPHA
  if (!is.null(random_alpha_uniform)) {
    n_draws <- R * K
    # for each named random alpha, replace the mean value in draws_exp by rnorm using sd = sigma parameter
    for (coef_name in names(random_alpha_uniform)) {
      sigma_name <- random_alpha_uniform[[coef_name]][1] # first element is sigma name
      if (!sigma_name %in% colnames(draws_exp)) stop(paste0("sigma '", sigma_name, "' not found in estimated coefficients/draws."))
      # per-row sd should be the repeated sigma draws from draws_exp
      sigma_vec <- draws_exp[, sigma_name]
      mean_vec <- draws_exp[, coef_name]
      # sample individual-level random coefficients
      new_vals <- runif(n = n_draws, min = mean_vec-abs(sigma_vec), max = mean_vec+abs(sigma_vec))
      draws_exp[, coef_name] <- new_vals
    }
    draws <- draws_exp
  }
  
  #extra draws for UNIFORM random coefficient for BETA
  if (!is.null(random_beta_uniform)) {
    n_draws <- R * K
    # for each named random alpha, replace the mean value in draws_exp by rnorm using sd = sigma parameter
    for (coef_name in names(random_beta_uniform)) {
      sigma_name <- random_beta_uniform[[coef_name]][1] # first element is sigma name
      if (!sigma_name %in% colnames(draws_exp)) stop(paste0("sigma '", sigma_name, "' not found in estimated coefficients/draws."))
      # per-row sd should be the repeated sigma draws from draws_exp
      sigma_vec <- draws_exp[, sigma_name]
      mean_vec <- draws_exp[, coef_name]
      # sample individual-level random coefficients
      new_vals <- runif(n = n_draws, min = mean_vec-abs(sigma_vec), max = mean_vec+abs(sigma_vec))
      draws_exp[, coef_name] <- new_vals
    }
    draws <- draws_exp
  }
  
  #Transform draws
  if(!is.null(trans_param)) {
    for(i in 1:length(trans_param)) {
      trans_fun <- get(names(trans_param)[i])
      for(j in 1:length(trans_param[[i]])) {
        coef_name <- trans_param[[i]][j]
        draws[, coef_name] <- trans_fun(draws[, coef_name])
      }
    }
  }


  alpha_mat <- NULL
  beta_mat <- NULL
  if(!is.null(interaction_alpha)) {
    # Convert to matrix for faster access
    database_indiv <- database[!duplicated(database$Respondent_Serial),] |> as.matrix()
    n_indiv <- nrow(database_indiv)
    
    # alpha matrix
    X_alpha <- matrix(
      1,
      nrow = n_indiv,
      ncol = length(alpha_name),
      dimnames = list(NULL, alpha_name)
    )
    
    for(param_name in names(interaction_alpha)){
      attri_name <- unlist(interaction_alpha[param_name], use.names = F)
      X_alpha[,param_name] <- database_indiv[, attri_name]
    }
    draws_alpha <- draws[, alpha_name]
    alpha_mat <- draws_alpha%*%t(X_alpha) # column: one individual across simulation;  row: one simulation across individuals;
    
    #Transform alpha_mat
    if(!is.null(trans_alpha)) {
      for(i in 1:length(trans_alpha)) {
        trans_alpha_fun <- get(trans_alpha[i])
        alpha_mat <- trans_alpha_fun(alpha_mat)
      }
    } 
  } else {
    draws_alpha <- draws[, alpha_name]
    if(is.matrix(draws_alpha)) { 
      alpha_vec <- rowSums(draws_alpha)
    } else {
      alpha_vec <- draws_alpha
    }
    #Transform alpha_vec
    if(!is.null(trans_alpha)) {
      for(i in 1:length(trans_alpha)) {
        trans_alpha_fun <- get(trans_alpha[i])
        alpha_vec <- trans_alpha_fun(alpha_vec)
      }
    }
  }
  
  
  if(!is.null(interaction_beta)) {
    database_indiv <- database[!duplicated(database$Respondent_Serial),] |> as.matrix()
    n_indiv <- nrow(database_indiv)
    
    # beta matrix
    X_beta <- matrix(
      1,
      nrow = n_indiv,
      ncol = length(beta_name),
      dimnames = list(NULL, beta_name)
    )
    for(param_name in names(interaction_beta)){
      attri_name <- unlist(interaction_beta[param_name], use.names = F)
      X_beta[,param_name] <- database_indiv[, attri_name]
    }
    draws_beta <- draws[, beta_name]
    beta_mat <- draws_beta%*%t(X_beta) 
    
    #Transform beta_mat
    if(!is.null(trans_beta)) {
      for(i in 1:length(trans_beta)) {
        trans_beta_fun <- get(trans_beta[i])
        beta_mat <- trans_beta_fun(beta_mat)
      }
    }
    
  } else {
    draws_beta <- draws[, beta_name]
    if(is.matrix(draws_beta)) {
      beta_vec <- rowSums(draws_beta) 
    } else {
      beta_vec <- draws_beta
    }
    #Transform beta_vec
    if(!is.null(trans_beta)) {
      for(i in 1:length(trans_beta)) {
        trans_beta_fun <- get(trans_beta[i])
        beta_vec <- trans_beta_fun(beta_vec)
      }
    }
  }
  
  
  #calculate WTP matrix and take means of rows -> average from one simulation across individuals 
  if(is.null(alpha_mat) & is.null(beta_mat))  {
    wtp_dist <- alpha_vec / beta_vec
  } else {
    if(is.null(alpha_mat)) alpha_mat <- alpha_vec%*%t(matrix(1, n_indiv))
    if(is.null(beta_mat)) beta_mat <- beta_vec%*%t(matrix(1, n_indiv))
    wtp_mat <- alpha_mat / beta_mat
    wtp_dist <- as.numeric(rowMeans(wtp_mat))
    
  }
  
  
  wtp_dist_hidden <- structure(wtp_dist, class = "hidden_vector")
  list(
    wtp_dist = function() wtp_dist_hidden,
    wtp_mean = mean(wtp_dist),
    wtp_median = median(wtp_dist),
    wtp_95_confidence = quantile(wtp_dist, c(0.025,0.975)) #95% simulation interval
  )
}






















