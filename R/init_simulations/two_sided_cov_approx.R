# PURPOSE: Initialize a simpler covariance approximation with fewer terms, and
#          compare with other approximations.

library(tidyverse)

# ------------------------------------------------------------------------------

# Define a function that returns multivariate normal data and the quantity to
# estimate for the covariance regression:
return_mvrnorm_samples <- function(n_tests, rho) {
  cov_matrix <- matrix(c(1, rho, rho, 1), nrow = 2, byrow = TRUE)
  
  MASS::mvrnorm(n_tests, c(0, 0), 
                Sigma = matrix(c(1, rho, rho, 1), nrow = 2, byrow = TRUE)) %>%
    as.data.frame() %>%
    rename(z1 = V1, z2 = V2) %>%
    mutate(z_rho = rho,
           pval1 = 2 * pnorm(-abs(z1)),
           pval2 = 2 * pnorm(-abs(z2)),
           chi_pval1 = -2 * log(pval1),
           chi_pval2 = -2 * log(pval2)) %>%
    return()
}

# Generate data then run a regression and fit polynomial without intercept through
# covariance approximation:
sim_pval_cov_data <- map_dfr(seq(-1, 1, by = 0.01),
                             function(sim_rho) {
                               return_mvrnorm_samples(1000000, sim_rho) %>%
                                 group_by(z_rho) %>%
                                 summarize(obs_cov = cov(chi_pval1, chi_pval2))
                             })

summary(lm(obs_cov ~ -1 + I(z_rho^2) + I(z_rho^4) + I(z_rho^6), 
           data = sim_pval_cov_data))

# Call:
#   lm(formula = obs_cov ~ -1 + I(z_rho^2) + I(z_rho^4) + I(z_rho^6), 
#      data = sim_pval_cov_data)
# 
#   Residuals:
#          Min         1Q     Median         3Q        Max 
#   -0.0209753 -0.0037117  0.0005219  0.0044475  0.0188766 
# 
#   Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#   I(z_rho^2) 3.902364   0.008348 467.442   <2e-16 ***
#   I(z_rho^4) 0.051520   0.027116   1.900   0.0589 .  
#   I(z_rho^6) 0.032832   0.020748   1.582   0.1151    
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
#   Residual standard error: 0.006785 on 198 degrees of freedom
#   Multiple R-squared:      1,	Adjusted R-squared:      1 
#   F-statistic: 4.651e+06 on 3 and 198 DF,  p-value: < 2.2e-16

# ------------------------------------------------------------------------------

# Now make functions for comparing the covariance approximation of this estimate
# with that of Yang (2016) and the recent GFisher pre-print:

compute_yang_cov_approx <- function(rho) {
  3.9081*(rho)^2 + 0.0313*(rho)^4 + 0.1022*(rho)^6 - 0.1378*(rho)^8 + 0.0941*(rho)^(10)
}

compute_gfisher_cov_approx <- function(rho) {
  3.9068*(rho)^2 + 0.0506*(rho)^4 + 0.0173*(rho)^6 + 0.0082*(rho)^8 + 0.0046*(rho)^(10)
}

compute_simpler_cov_approx <- function(rho) {
  3.902364*(rho)^2 + 0.051520*(rho)^4 + 0.032832*(rho)^6
}

# Compute the covariance approximation differences between each:
gfisher_yang_differences <- map_dbl(seq(-1, 1, by = 0.01),
                                  function(sim_rho) {
                                    compute_gfisher_cov_approx(sim_rho) - compute_yang_cov_approx(sim_rho)
                                  })
summary(gfisher_yang_differences)
max(abs(gfisher_yang_differences))
# [1] 0.0104

yang_simpler_differences <- map_dbl(seq(-1, 1, by = 0.01),
                                    function(sim_rho) {
                                      compute_yang_cov_approx(sim_rho) - compute_simpler_cov_approx(sim_rho)
                                    })
summary(yang_simpler_differences)
max(abs(yang_simpler_differences))
# [1] 0.011184

gfisher_simpler_differences <- map_dbl(seq(-1, 1, by = 0.01),
                                    function(sim_rho) {
                                      compute_gfisher_cov_approx(sim_rho) - compute_simpler_cov_approx(sim_rho)
                                    })
summary(gfisher_simpler_differences)
max(abs(gfisher_simpler_differences))
# [1] 0.0009217044



