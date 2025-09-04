
### This R file has all functions needed in the assignment of Group 2

#Load library
library(tidyverse)

## DATA CLEANING FUNCTIONS
#F1: Check that log-moneyness is in the interval [-0.2, 0.2]
check_log_moneyness <- function(data){
 data['log_moneyness'] <- log(data['strike']/data['underlying_price'])
 data_filtered <- data[data$log_moneyness >= -0.2 & data$log_moneyness <= 0.2, ]
 return(data_filtered)
}

#F2: Check Merton Constraints for Put options
check_merton_constraints <- function(data){
  r <- data$interest_rate
  d <- exp(-r*(data$ttm/252))
  discounted_strike <- data$strike*d
  data$payoff_t <- pmax(0, data$underlying_price - discounted_strike)
  data_no_arb <- data[data$option_price >= data$payoff_t &
                 data$option_price <= data$underlying_price, ]
  return(data_no_arb)
}

## IMPLIED VOLATILITY
#F3: Compute the B&S price for call options taken ttm as number of days
bs_call_opt<-function(S, K, ttm, sigma, rf){
  #Time to maturity expressed in Years
  ttm_y <- ttm/252
  d1<- (log(S/K) + (rf + 0.5*sigma^2)*(ttm_y))/(sigma * sqrt(ttm_y))
  d2 <- d1 - sigma * sqrt(ttm_y)
  call_price <- S*pnorm(d1) - K*exp(-rf*ttm_y)*pnorm(d2)
  return(call_price)
}

#F4: Compute Vega (i.e. derivative of the option w.r.t. volatility) 
# Exact analytical formula
compute_vega<-function(sigma, S, K, ttm, rf){
  ttm_y <- ttm / 252
  d1<- ( log(S/K) + (rf + 0.5 * sigma^2) * ttm_y)/(sigma * sqrt(ttm_y))
  vega <- S*dnorm(d1)*sqrt(ttm_y)
  return(vega)
}

#F5: Compute the error between the B&S price and the market price of an option
root_sigma <- function(par, S, K, ttm, rf, mkt_price){
  bs_price <- bs_call_opt(S = S, K = K, ttm = ttm, rf = rf, sigma = par)
  error <- bs_price - mkt_price
  return(error)
}

#F6: Newton Rapson algorithm to find the uni root of a function
newton_rapson_iv <- function(xold,S,K,ttm,rf,mkt_price,
                             max_iter=30, iter=0, cond = 10^-6){
  fold <- root_sigma(par = xold,
                     S = S,
                     K = K,
                     ttm = ttm,
                     rf = rf,
                     mkt_price = mkt_price)
  dfold <- compute_vega(sigma = xold,
                        S = S,
                        K = K,
                        ttm = ttm,
                        rf = rf)
  while(iter < max_iter | abs(fold) >= cond ){
    xnew <- xold - fold/ dfold
    xold <- xnew
    fold <- root_sigma(par = xold,
                       S = S,
                       K = K,
                       ttm = ttm,
                       rf = rf,
                       mkt_price = mkt_price)
    dfold <- compute_vega(sigma = xold,
                          S = S,
                          K = K,
                          ttm = ttm,
                          rf = rf)
    iter <- iter+1
    cat("\n",xnew, iter)
  }
  return(xnew)
}

#F7: Compute vega with the finite difference method
compute_approx_vega <- function(eps, sigma, S, K, ttm, rf){
  veg <- bs_call_opt(S, K, ttm, sigma, rf)
  veg_eps <- bs_call_opt(S, K, ttm, sigma + eps, rf)
  approx_vega <- (veg_eps - veg)/eps
  return(approx_vega)
}

## CALIBRATION
#F8: Compute the mean squared error between the B&S prices computed with a 
#certain volatility and the market prices of options
calibrationMSE_optimFUN <- function(S, K, r, ttm, par, mkt_price){
  error <- (bs_call_opt(S = S, K = K, ttm = ttm, r = r, sigma = par) - mkt_price)^2
  return(mean(error))
}

# ASIAN  OPTIONS PRICING VIA MONTE CARLO
#F9: Geometric Brownian Motion simulation for underlying asset price
simGBM_forloop <- function(S0, t0 = 0, r=0, ttm = 1, sigma = 1, mu=0, nsim=100, N=1000){
  S0_initial = rep(S0, nsim)
  deltat = (ttm - t0)/N
  path_gbm <- matrix(NA, nrow = nsim, ncol = N+1)
  path_gbm[, 1] <- S0_initial
  
  for (t in c(2:(N+1))) {
    path_gbm[,t] <- path_gbm[, t-1] * exp((r - 0.5*sigma^2)*(deltat) + 
                                            sigma * sqrt(deltat)*rnorm(nsim))
  }
  time_grid <- seq(t0, ttm, by =deltat)
  
  return(list(grid = time_grid, paths = round(path_gbm, 8)))
}

# pricing Asian put option with floating strike price
pricing_AsianPutOption <- function(S0, t0 = 0, ttm = 1, sigma = 1, r=0, 
                                   nsim = 100, N=100){
  
  sim_us <- simGBM_forloop(S0 = S0, t0 = t0, ttm = ttm, r = r,
                           sigma = sigma, nsim = nsim, N = N)$paths
  
  Kgeo_mean <- apply(sim_us, 1, function(x) exp(mean(log(x))))
  ST <- sim_us[, ncol(sim_us)]
  
  price_AsianPut <- exp(-r*(ttm - t0)) * pmax(Kgeo_mean - ST, 0)
  
  
  res <- mean(price_AsianPut)
  std_dev <- sqrt(var(price_AsianPut))
  
  LB <- res - qnorm(0.975) * sqrt(var(price_AsianPut))/sqrt(nsim)
  UB <- res + qnorm(0.975) * sqrt(var(price_AsianPut))/sqrt(nsim)
  
  return(list(exotic_price = res, std_dev = std_dev, LB = LB, UB = UB, n_simulation = nsim))
  
}
