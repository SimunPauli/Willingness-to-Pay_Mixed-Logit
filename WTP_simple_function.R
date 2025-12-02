library(MASS)

simulate_wtp_simple <- function(mu, Sigma, beta_name, alpha_name, R = 1e5, seed = 2025) {
  stopifnot(length(mu) == nrow(Sigma), nrow(Sigma) == ncol(Sigma))
  set.seed(seed)
  
  # Draw joint random coefficients ~ N(mu, Sigma)
  draws <- MASS::mvrnorm(n = R, mu = mu, Sigma = Sigma)
  
  if (!alpha_name %in% colnames(draws)) {
    stop("alpha_name not found in column names of 'draws'. Name your mu and Sigma dims.")
  }
  if (!beta_name %in% colnames(draws)) {
    stop("beta_name not found in column names of 'draws'.")
  }
  
  alpha <- draws[, alpha_name]
  beta  <- draws[, beta_name]
  
  # WTP as ratio; heavy tails possible if alpha ~ 0
  wtp <- - beta / alpha
  
  out <- list(
    wtp_vec           = wtp,
    n_draws           = length(wtp),
    share_alpha_pos   = mean(alpha > 0),          # check sign violations
    share_alpha_near0 = mean(abs(alpha) < 1e-4),  # near-zero denom indicator
    mean              = mean(wtp),
    sd                = sd(wtp),
    median            = median(wtp),
    p025              = quantile(wtp, 0.025),
    p05               = quantile(wtp, 0.05),
    p25               = quantile(wtp, 0.25),
    p75               = quantile(wtp, 0.75),
    p95               = quantile(wtp, 0.95),
    p975              = quantile(wtp, 0.975)
  )
  attr(out, "wtp_draws") <- wtp
  out
}

res <- simulate_wtp_simple(model$estimate[-which(names(model$estimate) == c("asc_med", "asc_icv"))], #Drop reference group
                           model$robvarcov, #refernce group not included in the covariance matrix
                           beta_name = "asc_sma", 
                           alpha_name = "b_pprice", 
                           R = 5000) 
res$median; res$p05; res$p95
