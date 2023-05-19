# This script is to estimate the transmission rates in snails


# This script perform the best curve for mu_m
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)


#load data  
temp <- c(7.5,             (7.5+20)/2,         20,             (20+27.5)/2,      27.5,        (27.5+36.5)/2,     36.5 )
rate <- c(2.275331, (2.275331+2.210028)/2,   2.210028,  (2.210028+2.24239)/2,   2.24239, (2.24239+4.354702)/2, 4.354702)


# keep just a single curve
d <- data.frame(temp, rate)

# fit model

fit <- nls_multstart(rate~thomas_2017(temp = temp, a,b,c,d,e),
                     data = .x,
                     iter = c(3,3,3,3,3),
                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') - 10,
                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') + 10,
                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'thomas_2017'),
                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'thomas_2017'),
                     supp_errors = 'Y')
fit


# calculate additional traits
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data

new_data <- data.frame(temp = seq(min(d$temp), max(d$temp), 0.5))
preds <- augment(fit, newdata = new_data)



t<-1

# R program to calculate root of an equation
fun <-function(beta) {log(1-i) - m_0*beta*exp(-(n*beta+mu)*t)/(n*beta+mu)+m_0*beta/(n*beta+mu)}




temp = 5
i = 1/100 

n = 
m_0 = 10*n
mu <- predict(fit, newdata = data.frame(temp))
beta_s <- uniroot(fun, lower = 0, upper = 1)
beta_s
