# This script is to estimate the transmission rates in snails

# This script perform the best curve for mu_m
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)


#R.M. Anderson et al. 1982, experiment with S. mansoni
temp_anderson_82 <- c(5,            10,        15,       20,        25,        30,        35,      40)
rate_anderson_82 <- c(4.957653,  2.141901,  1.50839,  1.974984,  2.513615,  4.322767,  4.490178,  5.114)


#S. K. Prah and C. James, 1977 S. mansoni 
temp_prah_77_man <- c(7.5,            20,      27.5,       36.5)
rate_prah_77_man <- c(2.56078,    0.6536434, 1.499468,   3.940166)


#R.E. Purnell, 1966, experiment with S. mansoni and two time range
temp_purnell_66 <- c(12,                       14,                           16,                 18.5,                        21.5,                     24.8,                  28.6,             32.7)
rate_purnell_66 <- c((3.261705+1.637893)/2,   (2.396054+1.785148)/2,     (2.904859+1.957561)/2,     (1.520372+3.054279)/2,        (2.337589+2.812790)/2,       (2.752958+3.891444)/2,  (3.842463+8.757026)/2, (4.877587+9.721674)/2)

#load data 
temp <- c(temp_anderson_82, temp_prah_77_man, temp_purnell_66)
rate <- c(rate_anderson_82, rate_prah_77_man, rate_purnell_66)


# keep just a single curve
d <- data.frame(temp, rate)
.x <- d
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





# R program to calculate root of an equation
fun <-function(beta) {log(1-i) - m_0*beta*exp(-(n*beta+mu)*(h/24))/(n*beta+mu)+m_0*beta/(n*beta+mu)}


temp <- 15
i <- 83/100


h <- 1
n <- 30
m_0 <- 5*n
mu <- predict(fit, newdata = data.frame(temp))
beta_s <- uniroot(fun, lower = 0, upper = 1)
beta_s
