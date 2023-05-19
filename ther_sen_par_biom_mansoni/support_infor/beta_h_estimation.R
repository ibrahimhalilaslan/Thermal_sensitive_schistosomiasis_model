#This script perform the selected model for mu_c
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)


#J.R. Lawson and R.A.Wilson 1980, experiment with S. mansoni 
temp_lawson_80 <- c(10,          15,           20,      25,        30,        35,       40)
rate_lawson_80 <- c(0.6844903, 0.347094,   0.9361688, 0.9384906, 1.767248, 3.79802, 7.210983)



#R.E. Purnell, 1966, experiment with S. mansoni
temp_purnell_66 <- c(12,             15,        18,            21,     24,           27,       30)
rate_purnell_66 <- c(0.9958894, 2.150529,  1.729903,       1.026659, 1.085069,  1.230383,  1.718108)



# load in data
temp <- c(temp_lawson_80, temp_purnell_66)
rate <- c(rate_lawson_80, rate_purnell_66)


# keep just a single curve
.x  <- data.frame(temp, rate)

# fit model

fit <- nls_multstart(rate~spain_1982(temp = temp, a,b,c,r0),
                     data = .x,
                     iter = c(4,4,4,4),
                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') - 1,
                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') + 1,
                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                     supp_errors = 'Y',
                     convergence_count = FALSE)
fit


# calculate additional traits
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data

new_data <- data.frame(temp = seq(min(d$temp), max(d$temp), 0.5))
preds <- augment(fit, newdata = new_data)





# R program to calculate root of an equation
fun <-function(beta) {worm+beta*mice*cercari/(beta*mice+mu)*exp(-(beta*mice+mu)*hour/24)-beta*mice*cercari/(beta*mice+mu)}


temp <- 25
worm <-  30.9 *76 * 12/100

hour <- .5
mice <- 12
mu <- predict(fit, newdata = data.frame(temp))
#time_pre_exp <- 24
#cercari <- 200*exp(-mu * time_pre_exp/24)
cercari <- 76 * 12
beta_s <- uniroot(fun, lower = 0, upper = 5)
beta_s
