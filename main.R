### This R file has results for the assignment of Group 2

# Load libraries needed
library(readxl)
library(ggplot2)

#Import functions
setwd('C:\\Users\\Admin\\Downloads\\Num')
source("functions.R")

#Import data
data <- read_xlsx('C:\\Users\\Admin\\Downloads\\Num\\group_2.xlsx')

#POINT 1: DATA CLEANING
#Filter only options with log-moneyness within the interval [-0.2, 0.2]
data_checked <- check_log_moneyness(data)

#Check Merton Constraints
no_arb_data <- check_merton_constraints(data_checked)

#Percentage of options in cleaned dataset
print(paste0('The % of options in the clean dataset vs. dirty is ', round(nrow(no_arb_data)/nrow(data)*100), '%' ))

#Verify Monotonicity and Convexity
plot(data$strike, data$option_price, 
     ylab = "Call price", xlab = "Strike price")

#Point 2: IMPLIED VOLATILITY
#2.1. Volatility smile using the Newton-Raphson algorithm

#Compute IVs with Newton Raphson method
xold <- 0.3 # starting guess/initial value
newton_rapson_vec <- c()
for (i in 1:nrow(no_arb_data)) {
  iv <- newton_rapson_iv(xold,
                         S = no_arb_data$underlying_price[i],
                         K = no_arb_data$strike[i],
                         ttm = no_arb_data$ttm[i],
                         rf = no_arb_data$interest_rate[i],
                         mkt_price = no_arb_data$option_price[i])
  newton_rapson_vec <- c(newton_rapson_vec, iv)
}

no_arb_data['iv_n_r'] <- newton_rapson_vec #Add column to the dataframe

#Compute IVs with uniroot method
uni_root_vec <- c()
for (i in 1:nrow(no_arb_data)) {
  iv <- uniroot(root_sigma, c(0.01, 0.8),
          S = no_arb_data$underlying_price[i],
          K = no_arb_data$strike[i],
          ttm = no_arb_data$ttm[i],
          rf = no_arb_data$interest_rate[i],
          mkt_price = no_arb_data$option_price[i])$root
  uni_root_vec <- c(uni_root_vec, iv)
}
#Add column to the dataframe
no_arb_data['iv_uniroot'] <- uni_root_vec

#Plot the two IVS
ggplot(no_arb_data, aes(x = strike)) + theme_light() +
  geom_point(aes(y = iv_n_r, color = "IVs_n_r", shape = "IVs_n_r"), alpha = 0.6) +
  geom_point(aes(y = iv_uniroot, color = "IVS_uni_root", shape = "IVS_uni_root"), alpha = 0.6) +
  xlab("K") + ylab("Implied Volatility") +
  scale_y_continuous("", limits = c(0.29, 0.42), breaks = seq(0.25, 0.45, 0.01)) +
  scale_color_manual(name = "Legend", values = c("IVs_n_r" = "blue", "IVS_uni_root" = "red")) +
  scale_shape_manual(name = "Legend", values = c("IVs_n_r" = 16, "IVS_uni_root" = 17)) +
  theme(legend.position = "top")

#2.2 Finite difference method
#Setting parameters
S <- seq(530, 1030, 5)
K <- 780 #Pick K = 780 since is the nearest ATM underlying = 780.56
r <- 0.01580662
sigma <- 0.3061677 #Volatility Newton Raphson
ttm <- 26
#Compute the exact vega
exact_vega_vect <- compute_vega(sigma= sigma,
             S =S,
             K = K,
             ttm = ttm,
             rf = r)
#setting small error
eps <- 0.0001
#Compute approximate Vega
approx_vega_vect <- compute_approx_vega(eps = eps,
                                sigma= sigma,
                                S =S,
                                K = K,
                                ttm = ttm,
                                rf = r)
#Creating the dataframe for ggplot
d_veg <- data.frame(cbind(S, exact_vega_vect, approx_vega_vect))

#Plot the two Vega
ggplot(d_veg, aes(x = S)) + theme_light() +
  geom_point(aes(y = exact_vega_vect, color = "Exact Vega", shape = "Exact Vega"), 
             alpha = 0.5, size = 2) +
  geom_line(aes(y = exact_vega_vect, color = "Exact Vega"), linewidth = 0.8) +
  
  geom_point(aes(y = approx_vega_vect, color = "Approx Vega", shape = "Approx Vega"), 
             alpha = 0.5, size = 2) +
  geom_line(aes(y = approx_vega_vect, color = "Approx Vega"), linewidth = 0.8, linetype = "dashed") +
  
  xlab("Underlying Price (S)") + ylab("Vega") +
  scale_color_manual(name = "Legend", 
                     values = c("Exact Vega" = "blue", "Approx Vega" = "red")) +
  scale_shape_manual(name = "Legend", 
                     values = c("Exact Vega" = 16, "Approx Vega" = 17)) +
  theme(legend.position = "top")

#Point 3: CALIBRATION OF VOLATILITY and PRICE OF BULL CALL SPREAD STRATEGY

calibrated_vola <- optim(par = 0.10, fn = calibrationMSE_optimFUN, lower = 0.01,
      method ="L-BFGS-B",
      S = S, K = no_arb_data$strike ,
      ttm = ttm,
      r = r, mkt_price =  no_arb_data$option_price)

calibrated_vola$par #calibrated implied volatility
calibrated_vola$convergence

#Setting parameters
S <- 780.56
K1 <- S*0.975
K2 <- S*1.025
ttm <- 60
r <- 0.01580662
sigma <- calibrated_vola$par

#Pricing the bull call spread
long_c <- bs_call_opt(S, K1, ttm, sigma, r)
short_c <- bs_call_opt(S, K2, ttm, sigma, r)
#Bull Call spread price
bull_call_spread <- long_c - short_c
bull_call_spread


###Point 4: PRICING FLOATING STRIKE ASIAN PUT OPTION
## Exact simulation of  GBM for underlying asset price 
sim <- simGBM_forloop(S0 = 780.56, t0 = 0, ttm = 65/252, r = 0.01580662,
                      sigma = calibrated_vola$par, mu=0, nsim=100, N=13)
matplot(x=sim$grid, y = t(sim$paths), type = "l", ylab = "Paths", xlab = "Grid")


## Monte Carlo pricing of floating strike Asian put option
set.seed(12345)
pricing_AsianPutOption(S0 = 780.56, sigma = calibrated_vola$par, r = 0.01580662, 
                       nsim = 10000, N = 13, ttm = 65/252)

## Compare with results in case of less number of simulations
pricing_AsianPutOption(S0 = 780.56, sigma = calibrated_vola$par, r = 0.01580662, 
                       nsim = 1000, N = 13, ttm = 65/252)

pricing_AsianPutOption(S0 = 780.56, sigma = calibrated_vola$par, r = 0.01580662, 
                       nsim = 100, N = 13, ttm = 65/252)

