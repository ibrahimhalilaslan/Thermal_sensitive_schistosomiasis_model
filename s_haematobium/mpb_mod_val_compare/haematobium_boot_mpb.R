############################# FIND THE CONFIDENT INTERVAL AND MPB ###############
# Clear the global environment 

rm(list = ls(all = TRUE)) 

# load packages
library(boot)
library(rTPC)
library(nls.multstart)
library(broom)
library(car)
library(magrittr)
library(dplyr)
library(tidyr)
library(purrr)
library(patchwork)
library(minpack.lm)


####################### SET CONSTANT PARAMETERS ###############

#The factor of density-dependent fecundity rate
# We adjust this value to get the abundance of sanil in a site same as C. Wood at al. 2019 paper
nu <- 5.662715e-07 # 


#The maximal rate of snail invasion by miracidia
# We adjust this value to get the the prevalence of disease in snail to be 5 percent
lambda <- 0.0009592123

#Constant human population having schistosomiasis 
# We adjust this value to get the number of human in a site similar as C. Wood at al. 2019 paper also Susanne H. Sokolow, 2015 has similar value
h <- 1000 # (human)



#The reduction factor of the fecundity rate due to infection
# Get from Susanne H. Sokolow, 2015
r <- 0.5


#The number of eggs produced by a mature parasite
nu_e <- 350/3 #(350 eggs per parasite in a day) 

#The maturation rate of parasite
sigma_p <- 0.02 #(~ 50 days) maturation time of parasite 

#The death rate of parasites due to human death
mu_h <- 1/(70*365)# (70 years) life span of human

#The death rate of parasite 
mu_p <- 1/(8 * 365) # (8 years) life span of parasite


# set reduction rate for the transmission rate in humans
contact_reduction_h <- 1.542194e-08

# set the reduction rate for transmission rate in snails
contact_reduction_s <-0.32

####### START BOOTSTRAPPING 
############################# TRANSMSISION RATE IN HUMAN #################


min_temp <- 12
max_temp <- 37

# Set the number of row to bootstrap to remove becuase 
# we do not have sufficient data for haematobium and the bootstrap give error last a few bootstrap
n <- 25

#(Purnell 1966),, S. mansoni
purnell_66_temp <- c(       12,             15,           18,             21,        24,          27,           30,       33)
purnell_66_rate_1 <- c(  0.9246078,      0.810981,    1.320736,       1.092269,   0.6998873,   1.231836,    0.7901292,  0.6292689)
#normalize the data set with the maximum value 
purnell_66_rate_1 <- purnell_66_rate_1/max(purnell_66_rate_1)


purnell_66_temp <- c(       12,             15,           18,             21,        24,          27,           30,       33)
purnell_66_rate_2 <- c(  1.0535,         1.002849,     0.9284408,       1.481275,   0.7225294,   0.6554545,    1.251547,  0.7475422)
#normalize the data set with the maximum value 
purnell_66_rate_2 <- purnell_66_rate_2/max(purnell_66_rate_2)


purnell_66_temp <- c(       12,             15,           18,             21,        24,          27,           30,       33)
purnell_66_rate_3 <- c(  2.182979,      1.523512,    1.777458,       1.840471,   0.8305635,   0.7054314,      1.339291,  0)
#normalize the data set with the maximum value 
purnell_66_rate_3 <- purnell_66_rate_3/max(purnell_66_rate_3)

#(Stirewalt 1954), S. mansoni 
stirewalt_54_temp <- c(    24,         27)
stirewalt_54_rate <- c(0.09590499,  0.1552625)
#normalize the data set with the maximum value 
stirewalt_54_rate <- stirewalt_54_rate/max(stirewalt_54_rate)


#(DeWitt 1965) with S. mansoni
dewitt_55_temp_1 <- c(     15,      20,          25,      30,      35)
dewitt_55_rate_1 <- c(0.4618475, 0.8464159,  1.00836, 1.237522, 0.6192221)
#normalize the data set with the maximum value 
dewitt_55_rate_1 <- dewitt_55_rate_1/max(dewitt_55_rate_1)

dewitt_55_temp_2 <- c(     15,      20,          25,      30,         35,        40)
dewitt_55_rate_2 <- c(0.385854, 0.6910937,    0.635236,  0.9783329, 1.0652, 0.2764123)
#normalize the data set with the maximum value 
dewitt_55_rate_2 <- dewitt_55_rate_2/max(dewitt_55_rate_2)

dewitt_55_temp_3 <- c(   10,        15,      20,         25,      30,       35,        40,   45)
dewitt_55_rate_3 <- c(0.7439515,  1.76077,  2.484149,  3.50177,  3.9655, 3.421058,  0.7532138, 0)
#normalize the data set with the maximum value 
dewitt_55_rate_3 <- dewitt_55_rate_3/max(dewitt_55_rate_3)

dewitt_55_temp_4 <- c(0, 5,        10,        15,      20,        25)
dewitt_55_rate_4 <- c(0, 0,    0.1190698, 0.6479505, 1.39834, 1.50036)
#normalize the data set with the maximum value 
dewitt_55_rate_4 <- dewitt_55_rate_4/max(dewitt_55_rate_4)

#combine the data set in one row 
temp <- c(purnell_66_temp, purnell_66_temp, purnell_66_temp,stirewalt_54_temp, dewitt_55_temp_1, dewitt_55_temp_2, dewitt_55_temp_3, dewitt_55_temp_4)
rate <- c(purnell_66_rate_1, purnell_66_rate_2, purnell_66_rate_3, stirewalt_54_rate, dewitt_55_rate_1, dewitt_55_rate_2, dewitt_55_rate_3, dewitt_55_rate_4)

# keep just a single curve
d <- data.frame(temp, rate)

# fit Sharpe-Schoolfield model
d_fit <- nest(d, data = c(temp, rate)) %>%
  mutate(briere2_1999 = map(data, ~nls_multstart(rate~briere2_1999(temp = temp, tmin, tmax, a,b),
                                                 data = d,
                                                 iter = c(4,4,4,4),
                                                 start_lower = get_start_vals(d$temp, d$rate, model_name = 'briere2_1999') - 10,
                                                 start_upper = get_start_vals(d$temp, d$rate, model_name = 'briere2_1999') + 10,
                                                 lower = get_lower_lims(d$temp, d$rate, model_name = 'briere2_1999'),
                                                 upper = get_upper_lims(d$temp, d$rate, model_name = 'briere2_1999'),
                                                 supp_errors = 'Y',
                                                 convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), by = 0.1))),
         # predict over that data,
         preds =  map2(briere2_1999, new_data, ~augment(.x, newdata = .y)))


# refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(rate~briere2_1999(temp = temp, tmin, tmax, a,b),
                               data = d,
                               start = coef(d_fit$briere2_1999[[1]]),
                               lower = get_lower_lims(d$temp, d$rate, model_name = 'briere2_1999'),
                               upper = get_upper_lims(d$temp, d$rate, model_name = 'briere2_1999'),
                               weights = rep(1, times = nrow(d)))

# bootstrap using case resampling
boot1 <- Boot(fit_nlsLM, method = 'case')


# create predictions of each bootstrapped model
boot1_preds <- boot1$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min_temp, max_temp, by = 0.1))) %>%
  ungroup() %>%
  mutate(pred = briere2_1999(temp = temp, tmin, tmax, a,b))

# remove the bootstrap does not exist 
boot1_preds_beta_h <- head(boot1_preds, -n*length(seq(min_temp, max_temp, by = 0.1)))

################################ TRANSMSISION RATE IN SNAIL 

#(Prah and James 1977), Bulinus (Physopsis) globosus with S. haematobium.
prah_77_temp <- c(       5,        12,            19,     26)
prah_77_rate <- c(0.001279126, 0.04565688,  0.1084776, 0.1849867)
#normalize the data set with the maximum value 
prah_77_rate <- prah_77_rate/max(prah_77_rate)

#(Shiff 1974), Bulinus globosus infected with S. haematobium.
shiff_74_temp <- c(      25,          20,           18,         16,          13,            11) 
shiff_74_rate <- c(0.03730699,  0.003512327,  0.01022082, 0.001782906,  0.002553961,  0.001281865)
#normalize the data set with the maximum value 
shiff_74_rate <- shiff_74_rate/max(shiff_74_rate)

#Chu et al. 1966), Bulinus trancatus  infected with S. haematobium. 
chu_66_temp <- c(    10,            12,           14,         15,        20,          25,         30,        35,         38)
chu_66_rate <- c(0.002009499,  0.008145561,  0.02029471, 0.02016721, 0.04754069, 0.03577785,  0.06475144, 0.04395093, 0.05660857)
#normalize the data set with the maximum value 
chu_66_rate <- chu_66_rate/max(chu_66_rate)           

# combine the data set in one row 
temp <- c(prah_77_temp, shiff_74_temp, chu_66_temp)
rate <- c(prah_77_rate, shiff_74_rate, chu_66_rate)


# keep just a single curve
d <- data.frame(temp, rate)

# fit Sharpe-Schoolfield model to raw data
d_fit <- nest(d, data = c(temp, rate)) %>%
  mutate(flinn_1991 = map(data, ~nls_multstart(rate~flinn_1991(temp = temp, a, b, c),
                                               data = .x,
                                               iter = c(5,5,5),
                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') - 10,
                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') + 10,
                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                                               supp_errors = 'Y',
                                               convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), by = 0.1))),
         # predict over that data,
         preds =  map2(flinn_1991, new_data, ~augment(.x, newdata = .y)))

# refit model using nlsLM
fit_nlsLM2 <- nlsLM(rate~flinn_1991(temp = temp, a, b, c),
                    data = d,
                    start = coef(d_fit$flinn_1991[[1]]),
                    lower = get_lower_lims(d$temp, d$rate, model_name = 'flinn_1991'),
                    upper = get_upper_lims(d$temp, d$rate, model_name = 'flinn_1991'),
                    control = nls.lm.control(maxiter=500),
                    weights = rep(1, times = nrow(d)))


# bootstrap using residual resampling
boot3 <- Boot(fit_nlsLM2, method = 'residual')

# predict over new data
boot3_preds <- boot3$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min_temp, max_temp, by = 0.1))) %>%
  ungroup() %>%
  mutate(pred = flinn_1991(temp = temp, a, b, c))

# remove the bootstrap does not exist 
boot3_preds_beta_s <- head(boot3_preds, -n*length(seq(min_temp, max_temp, by = 0.1)))

############################## CERCARIAL RELEASE RATE #################

#Nguyen et al. 2021, S. mansoni with Biomphalaria glabrata
temp_nguyen_21 <- c(5,  9,     13,     17,     21,     25,     29,     33,     37)
rate_nguyen_21 <- c(0, 15,     260,   1350,   1970,   2900,   1970,   1770,   690)

#(Upatham 1973), B. glabrata with S. mansoni  
temp_upatham_73 <- c(10, 13, 16,   19, 22,    25,   28,  31,   34,   37, 40)
rate_upatham_73 <- c( 0,  0, 206, 114, 442, 1200, 1400, 1233, 1241, 243, 121)

#combine the data set in one row 
temp <- c(temp_nguyen_21, temp_upatham_73)
rate <- c(rate_nguyen_21, rate_upatham_73)

# keep just a single curve
d <- data.frame(temp, rate)

# fit Sharpe-Schoolfield model
d_fit <- nest(d, data = c(temp, rate)) %>%
  mutate(gaussian_1987 = map(data, ~nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                                  data = .x,
                                                  iter = c(4,4,4),
                                                  start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                                                  start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                                                  lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                                  upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                                  supp_errors = 'Y',
                                                  convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), by = 0.1))),
         # predict over that data,
         preds =  map2(gaussian_1987, new_data, ~augment(.x, newdata = .y)))


# refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(rate~gaussian_1987(temp = temp, rmax, topt, a),
                               data = d,
                               start = coef(d_fit$gaussian_1987[[1]]),
                               lower = get_lower_lims(d$temp, d$rate, model_name = 'gaussian_1987'),
                               upper = get_upper_lims(d$temp, d$rate, model_name = 'gaussian_1987'),
                               weights = rep(1, times = nrow(d)))

# bootstrap using case resampling
boot1 <- Boot(fit_nlsLM, method = 'case')

# create predictions of each bootstrapped model
boot1_preds <- boot1$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min_temp, max_temp, by = 0.1))) %>%
  ungroup() %>%
  mutate(pred = gaussian_1987(temp = temp, rmax, topt, a))

# remove the bootstrap does not exist 
boot1_preds_nu_c <- head(boot1_preds, -n*length(seq(min_temp, max_temp, by = 0.1)))

########################## FECUNDITY RATE OF SNAIL ######################

#Chester Kalinda. 2017, experiment with Bulinus globosus. 
temp_ches_17 <- c(15.5,    21.2,   25.8,    31, 35.5)
rate_ches_17 <- c(  0,    0.027,  0.037,  0.019, 0)


#El- Hassan 1974, Bulinus truncatus
temp_has_74_t <- c(15,        20,       25,       30,    35)
rate_has_74_t <- c(0.0103,   0.025,   0.0348,  0.0199,    0)


#(Shiff 1967), Bulinus globosus
temp_shiff_67 <- c(  18,      22,      25,      27)
rate_shiff_67 <- c(0.0137,  0.0447,  0.0497,  0.0391)


#(Kubirizajournal et al. 2010), Bulinus nyassanus
temp_kubir_10 <- c(  22,     25,    28,    31)
rate_kubir_10 <- c(0.0297, 0.0427, 0.0463, 0.0451)

# We assume there is no egg production at 32, 33, 34 
temp <-c(temp_ches_17, temp_has_74_t,  temp_shiff_67, temp_kubir_10,  32, 33, 34)
rate <-c(rate_ches_17, rate_has_74_t,  rate_shiff_67, rate_kubir_10,   0, 0, 0)

# keep just a single curve
d <- data.frame(temp, rate)


# fit Sharpe-Schoolfield model
d_fit <- nest(d, data = c(temp, rate)) %>%
  mutate(johnsonlewin_1946 = map(data, ~nls_multstart(rate~johnsonlewin_1946(temp = temp, r0, e, eh, topt),
                                                      data = .x,
                                                      iter = c(4,4,4,4),
                                                      start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'johnsonlewin_1946') - 1,
                                                      start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'johnsonlewin_1946') + 1,
                                                      lower = get_lower_lims(.x$temp, .x$rate, model_name = 'johnsonlewin_1946'),
                                                      upper = get_upper_lims(.x$temp, .x$rate, model_name = 'johnsonlewin_1946'),
                                                      supp_errors = 'Y',  
                                                      convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), by = 0.1))),
         # predict over that data,
         preds =  map2(johnsonlewin_1946, new_data, ~augment(.x, newdata = .y)))

# refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(rate~johnsonlewin_1946(temp = temp, r0, e, eh, topt),
                               data = d,
                               start = coef(d_fit$johnsonlewin_1946[[1]]),
                               lower = get_lower_lims(d$temp, d$rate, model_name = 'johnsonlewin_1946'),
                               upper = get_upper_lims(d$temp, d$rate, model_name = 'johnsonlewin_1946'),
                               weights = rep(1, times = nrow(d)))

# bootstrap using case resampling
boot1 <- Boot(fit_nlsLM, method = 'case')

# create predictions of each bootstrapped model
boot1_preds <- boot1$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min_temp, max_temp, by = 0.1))) %>%
  ungroup() %>%
  mutate(pred = johnsonlewin_1946(temp = temp, r0, e, eh, topt))

# remove the bootstrap does not exist 
boot1_preds_nu_s <- head(boot1_preds, -n*length(seq(min_temp, max_temp, by = 0.1)))

################################# LATENT PERIOD OF SNAIL ##################
#W. Pflüger 1984, Bulinus truncatus with S. haematobium 
temp_pfluger_84 <- c(18,       19,     20,    21,    22,     23,    25,  28,   30,    31,  32)
rate_pfluger_84 <- c(1/117,   1/110,  1/79, 1/60,   1/57,  1/56,  1/44, 1/32, 1/30, 1/31, 1/30)


#R. M. Gordon, 1934, Physopsis (bolinus) globosa  with haematobium 
temp_gordon_34_glo <- c(22,     26.3,      26.4,     26.8,       31.9,    32.1,     33,     35.2)
rate_gordon_34_glo <- c(1/67,  1/38.5,    1/39.3,   1/34.8,     1/27.5,   1/23,   1/25.5,   1/26.3)


#load data, Latent period of snail
temp <- c(temp_pfluger_84, temp_gordon_34_glo)
rate <- c(rate_pfluger_84, rate_gordon_34_glo)

# keep just a single curve
d <- data.frame(temp, rate)

# fit model
d_fit <- nest(d, data = c(temp, rate)) %>%
  mutate(gaussian_1987 = map(data, ~nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                                  data = .x,
                                                  iter = c(4,4,4),
                                                  start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                                                  start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                                                  lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                                  upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                                  supp_errors = 'Y',
                                                  convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), by = 0.1))),
         # predict over that data,
         preds =  map2(gaussian_1987, new_data, ~augment(.x, newdata = .y)))


# refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(rate~gaussian_1987(temp = temp, rmax, topt, a),
                               data = d,
                               start = coef(d_fit$gaussian_1987[[1]]),
                               lower = get_lower_lims(d$temp, d$rate, model_name = 'gaussian_1987'),
                               upper = get_upper_lims(d$temp, d$rate, model_name = 'gaussian_1987'),
                               weights = rep(1, times = nrow(d)))

# bootstrap using case resampling
boot1 <- Boot(fit_nlsLM, method = 'case')


# create predictions of each bootstrapped model
boot1_preds <- boot1$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min_temp, max_temp, by = 0.1))) %>%
  ungroup() %>%
  mutate(pred = gaussian_1987(temp = temp, rmax, topt, a))

# remove the bootstrap does not exist 
boot1_preds_sigma_s <- head(boot1_preds, -n*length(seq(min_temp, max_temp, by = 0.1)))

############################# MIRACIDIA HATCH RATE ############## 

#Miracidia hatching rate  (Nguyen et al. 2021, Schistosoma mansoni)
temp <- c(5,  9,  13, 17,  21, 22, 25,  29,  33, 37)
rate <- c(70, 60, 90, 100, 76, 114, 137, 105, 150, 135)/300

# keep just a single curve
d <- data.frame(temp, rate)

# fit Sharpe-Schoolfield model
d_fit <- nest(d, data = c(temp, rate)) %>%
  mutate(flinn_1991 = map(data, ~nls_multstart(rate~flinn_1991(temp = temp, a, b, c),
                                               data = .x,
                                               iter = c(5,5,5),
                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') - 10,
                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') + 10,
                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                                               supp_errors = 'Y')),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), by = 0.1))),
         # predict over that data,
         preds =  map2(flinn_1991, new_data, ~augment(.x, newdata = .y)))

# refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(rate~flinn_1991(temp = temp, a, b, c), 
                               data = d,
                               start = coef(d_fit$flinn_1991[[1]]),
                               lower = get_lower_lims(d$temp, d$rate, model_name = 'flinn_1991'),
                               upper = get_upper_lims(d$temp, d$rate, model_name = 'flinn_1991'),
                               weights = rep(1, times = nrow(d)))

# bootstrap using case resampling
boot1 <- Boot(fit_nlsLM, method = 'case')

# create predictions of each bootstrapped model
boot1_preds <- boot1$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min_temp, max_temp, by = 0.1))) %>%
  ungroup() %>%
  mutate(pred = flinn_1991(temp = temp, a, b, c))

# remove the bootstrap does not exist 
boot1_preds_delta_e <- head(boot1_preds, -n*length(seq(min_temp, max_temp, by = 0.1)))

##################### MORTALITY RATE OF CERCARIA ##################

#J.R. Lawson and R.A.Wilson 1980, experiment with S. mansoni 
temp_lawson_80 <- c(10,          15,           20,      25,        30,        35,       40)
rate_lawson_80 <- c(0.6844903, 0.347094,   0.9361688, 0.9384906, 1.767248, 3.79802, 7.210983)

#R.E. Purnell, 1966, experiment with S. mansoni
temp_purnell_66 <- c(12,             15,        18,            21,     24,           27,       30)
rate_purnell_66 <- c(0.9958894, 2.150529,  1.729903,       1.026659, 1.085069,  1.230383,  1.718108)

#Combine the data set in one row 
temp <- c(temp_lawson_80, temp_purnell_66)
rate <- c(rate_lawson_80, rate_purnell_66)

# keep just a single curve
d <- data.frame(temp, rate)

# fit Sharpe-Schoolfield model to raw data
d_fit <- nest(d, data = c(temp, rate)) %>%
  mutate(spain_1982 = map(data, ~nls_multstart(rate~spain_1982(temp = temp, a,b,c,r0),
                                               data = .x,
                                               iter = c(4,4,4,4),
                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') - 1,
                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') + 1,
                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                                               supp_errors = 'Y',
                                               convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), by = 0.1))),
         # predict over that data,
         preds =  map2(spain_1982, new_data, ~augment(.x, newdata = .y)))

# refit model using nlsLM
fit_nlsLM2 <- nlsLM(rate~spain_1982(temp = temp, a,b,c,r0),
                    data = d,
                    start = coef(d_fit$spain_1982[[1]]),
                    lower = get_lower_lims(d$temp, d$rate, model_name = 'spain_1982'),
                    upper = get_upper_lims(d$temp, d$rate, model_name = 'spain_1982'),
                    control = nls.lm.control(maxiter=500),
                    weights = rep(1, times = nrow(d)))


# bootstrap using residual resampling
boot3 <- Boot(fit_nlsLM2, method = 'residual')

# predict over new data
boot3_preds <- boot3$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min_temp, max_temp, by = 0.1))) %>%
  ungroup() %>%
  mutate(pred = spain_1982(temp = temp, a,b,c,r0))

# remove the bootstrap does not exist 
boot3_preds_mu_c <- head(boot3_preds, -n*length(seq(min_temp, max_temp, by = 0.1)))

######################### MORTALITY RATE OF INFECTED SNAILS #####################

# Chester Kalinda 2017, Bulinus globosu with S. haematobium
temp_chester_17 <- c(15.5,            21.2,        25.8,           31,         36)
rate_chester_17 <- c(0.00244925,  0.006901087,   0.01151302,   0.02073716,  0.04642702)

# W. Pflüger 1984, Bulinus truncatus with S. haematobium
temp_pfluger_84 <- c(18,              19,         20,          21,                 22,      23,             25,            28,        30,          31,    32)
rate_pfluger_84 <- c(0.0131776,  0.01141752,     0.01963277, 0.02130685,     0.02204905, 0.02636316,    0.02185846,   0.01590461, 0.01570362,  0.02473268, 0.0262486)

#(Chu et al. 1966) Bulinus trancatus with S. haematobium  
temp_chu_66 <- c(   9.967,          13,          16.067,         19,         22.2,         25.067,          28,        31.133,        34.2,         37.067,      39.967)
rate_chu_66 <- c(0.0058551278, 0.0085350331, 0.0061192732, 0.0195844770, 0.0085350331, 0.0032587267, 0.0004931228, 0.0061192732, 0.0058551278, 0.0968812480, 0.0915518574)

#Combine the data one curve 
temp <- c(temp_chester_17, temp_pfluger_84, temp_chu_66)    
rate <- c(rate_chester_17, rate_pfluger_84, rate_chu_66)

# keep just a single curve
d <- data.frame(temp, rate)

# fit Sharpe-Schoolfield model to raw data
d_fit <- nest(d, data = c(temp, rate)) %>%
  mutate(spain_1982 = map(data, ~nls_multstart(rate~spain_1982(temp = temp, a,b,c,r0),
                                               data = .x,
                                               iter = c(4,4,4,4),
                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') - 1,
                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') + 1,
                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                                               supp_errors = 'Y',
                                               convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), by = 0.1))),
         # predict over that data,
         preds =  map2(spain_1982, new_data, ~augment(.x, newdata = .y)))

# refit model using nlsLM
fit_nlsLM2 <- nlsLM(rate~spain_1982(temp = temp, a,b,c,r0),
                    data = d,
                    start = coef(d_fit$spain_1982[[1]]),
                    lower = get_lower_lims(d$temp, d$rate, model_name = 'spain_1982'),
                    upper = get_upper_lims(d$temp, d$rate, model_name = 'spain_1982'),
                    control = nls.lm.control(maxiter=500),
                    weights = rep(1, times = nrow(d)))

# bootstrap using residual resampling
boot3 <- Boot(fit_nlsLM2, method = 'residual')

# predict over new data
boot3_preds_1 <- boot3$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d$temp), max(d$temp), by = 0.1))) %>%
  ungroup() %>%
  mutate(pred = spain_1982(temp = temp, a,b,c,r0))
################### The second bootstrap for higher temperature ########

# chester_kalinda_2017_Bulinus_globosus 
temp_ches_17_2 <-      35.5
rate_ches_17_2 <-  0.04642702


#p_h_joubert_1986_Bul_globosus
temp_joubert_86_g <- c(     34,         36,       38,       40)
rate_joubert_86_g <- c(0.02924029, 0.1490773, 0.8022027, 3.540546)


#p_h_joubert_1986_Bul_africanus
temp_joubert_86_a <- c(     34,         36,       38,       40)
rate_joubert_86_a <- c(0.1044946, 0.5664857, 1.618022,  4.344277)

# El_hassan_1974
temp_has_74_t_2 <-     35
rate_has_74_t_2 <-  0.09622

#Combine the data set one row 
temp <- c(temp_ches_17_2,temp_joubert_86_g, temp_joubert_86_a, temp_has_74_t_2)
rate <- c(rate_ches_17_2, rate_joubert_86_g, rate_joubert_86_a, rate_has_74_t_2)

# keep just a single curve
d_2 <- data.frame(temp, rate)


# fit Sharpe-Schoolfield model to raw data
d_fit <- nest(d_2, data = c(temp, rate)) %>%
  mutate(sharpeschoollow_1981 = map(data, ~nls_multstart(rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 15),
                                                         data = .x,
                                                         iter = c(4,4,4,4),
                                                         start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981') - 10,
                                                         start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981') + 10,
                                                         lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981'),
                                                         upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981'),
                                                         supp_errors = 'Y',
                                                         convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp),  by = 0.1))),
         # predict over that data,
         preds =  map2(sharpeschoollow_1981, new_data, ~augment(.x, newdata = .y)))


# refit model using nlsLM
fit_nlsLM2 <- nlsLM(rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 15),
                    data = d_2,
                    start = coef(d_fit$sharpeschoollow_1981[[1]]),
                    lower = get_lower_lims(d_2$temp, d_2$rate, model_name = 'sharpeschoollow_1981'),
                    upper = get_upper_lims(d_2$temp, d_2$rate, model_name = 'sharpeschoollow_1981'),
                    control = nls.lm.control(maxiter=500),
                    weights = rep(1, times = nrow(d_2)))


# bootstrap using residual resampling
boot3 <- Boot(fit_nlsLM2, method = 'residual')


# predict over new data
boot3_preds_2 <- boot3$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d_2$temp), max(d_2$temp), by = 0.1))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 15))

# create an empty row to combine the result of two bootstrap 
boot_preds_mu_i <- data.frame(matrix(0, ncol = 2, nrow = length(boot1_preds_beta_h$temp)))
colnames(boot_preds_mu_i) = c("temp","pred") 

# Set the number of division for the number of bootstrap 
temp_row_length <- length(seq(min_temp, max_temp, by = 0.1))

# Set the number of first bootstrap end 
temp_row_length_1 <- length(seq(min(d$temp), max(d$temp), by = 0.1))

# Set the number of second bootstrap start
temp_row_length_2 <- length(seq(min(d_2$temp), max(d_2$temp), by = 0.1))

# start bootstrap to combine both bootsrapt 
for (k in 1:(length(boot1_preds_beta_h$temp)/temp_row_length)){
  
  temp_boot_1 <- boot3_preds_1[(1+(k-1)*temp_row_length_1):(k*temp_row_length_1),]
  temp_boot_2 <- boot3_preds_2[(1+(k-1)*temp_row_length_2):(k*temp_row_length_2),]
  
  boot_preds_mu_i[(1+(k-1)*temp_row_length):(k*temp_row_length), 1:2] <- rbind(temp_boot_1[which(temp_boot_1$temp <= max_temp & temp_boot_1$temp > min_temp-0.1), 6:7], temp_boot_2[which(temp_boot_2$temp > 2*max_temp & temp_boot_2$temp <= max_temp), 6:7])
  
}

################################ MORTALITIY OF MIRACIDIA ################

#S. K. Prah and C. James, 1977 S. haematobium
temp<- c(7.5,              20,     27.5,       36.5)
rate <- c(2.275331,      2.210028,  2.24239,  4.354702)

# keep just a single curve
d <- data.frame(temp, rate)

# fit Sharpe-Schoolfield model to raw data
d_fit <- nest(d, data = c(temp, rate)) %>%
  mutate(spain_1982 = map(data, ~nls_multstart(rate~spain_1982(temp = temp, a,b,c,r0),
                                               data = .x,
                                               iter = c(4,4,4,4),
                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') - 1,
                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') + 1,
                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                                               supp_errors = 'Y',
                                               convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), by = 0.1))),
         # predict over that data,
         preds =  map2(spain_1982, new_data, ~augment(.x, newdata = .y)))

# refit model using nlsLM
fit_nlsLM2 <- nlsLM(rate~spain_1982(temp = temp, a,b,c,r0),
                    data = d,
                    start = coef(d_fit$spain_1982[[1]]),
                    lower = get_lower_lims(d$temp, d$rate, model_name = 'spain_1982'),
                    upper = get_upper_lims(d$temp, d$rate, model_name = 'spain_1982'),
                    control = nls.lm.control(maxiter=500),
                    weights = rep(1, times = nrow(d)))


# bootstrap using residual resampling
boot3 <- Boot(fit_nlsLM2, method = 'residual')

# predict over new data
boot3_preds <- boot3$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min_temp, max_temp, by = 0.1))) %>%
  ungroup() %>%
  mutate(pred = spain_1982(temp = temp, a,b,c,r0))

# remove the bootstrap does not exist 
boot3_preds_mu_m <- head(boot3_preds, -n*length(seq(min_temp, max_temp, by = 0.1)))

############################ MORTALITY OF SNAIL ##############################

# chester_kalinda_2017_Bulinus_globosus 
temp_ches_17 <- c(15.5,           21.2,        25.8,      31)
rate_ches_17 <- c(0.00244925, 0.006901087, 0.01151302, 0.02073716)


#El- Hassan 1974, Bulinus truncatus
temp_has_74_t <- c(   10,        15,     20,       25,      30)
rate_has_74_t <- c(0.00442,   0.00218, 0.00144, 0.00518, 0.03077)

#(El-Emam and Madsen 1982), Bulinus truncatus
temp_eleman_82 <- c(    10,         18,          26,        28 )
rate_eleman_82 <- c(0.00256859, 0.004560943, 0.003184598, 0.004560943)

#(Shiff 1964), Bulinus globosus
temp_shiff_64 <- c(   18,           22,            25,        27)
rate_shiff_64 <- c(0.005314696, 0.006322683, 0.01045628, 0.01057945)

#(Kubirizajournal et al. 2010),  Bulinus nyassanus
temp_kubirizajournal_10 <- c(22,              25,          28,          31)
rate_kubirizajournal_10 <- c(0.006849573,  0.00513718, 0.02054872, 0.01027436)


# Combine the data in one row 
temp <- c(temp_ches_17, temp_has_74_t,  temp_eleman_82, temp_shiff_64,  temp_kubirizajournal_10)
rate <- c(rate_ches_17, rate_has_74_t,  rate_eleman_82, rate_shiff_64,  rate_kubirizajournal_10)

# keep just a single curve
d <- data.frame(temp, rate)

# fit model
d_fit <- nest(d, data = c(temp, rate)) %>%
  mutate(quadratic_2008 = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                                   data = .x,
                                                   iter = c(4,4,4),
                                                   start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 0.5,
                                                   start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 0.5,
                                                   lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                                   upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                                   supp_errors = 'Y',
                                                   convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), by = 0.1))),
         # predict over that data,
         preds =  map2(quadratic_2008, new_data, ~augment(.x, newdata = .y)))


# refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(rate~quadratic_2008(temp = temp, a, b, c),
                               data = d,
                               start = coef(d_fit$quadratic_2008[[1]]),
                               lower = get_lower_lims(d$temp, d$rate, model_name = 'quadratic_2008'),
                               upper = get_upper_lims(d$temp, d$rate, model_name = 'quadratic_2008'),
                               weights = rep(1, times = nrow(d)))

# bootstrap using case resampling
boot1 <- Boot(fit_nlsLM, method = 'case')

# create predictions of each bootstrapped model
boot1_preds_1 <- boot1$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d$temp), max_temp, by = 0.1))) %>%
  ungroup() %>%
  mutate(pred = quadratic_2008(temp = temp, a, b, c))


# chester_kalinda_2017_Bulinus_globosus 
temp_ches_17_2 <-      35.5
rate_ches_17_2 <-  0.04642702


#p_h_joubert_1986_Bul_globosus
temp_joubert_86_g <- c(     34,         36,       38,       40)
rate_joubert_86_g <- c(0.02924029, 0.1490773, 0.8022027, 3.540546)


#p_h_joubert_1986_Bul_africanus
temp_joubert_86_a <- c(     34,         36,       38,       40)
rate_joubert_86_a <- c(0.1044946, 0.5664857, 1.618022,  4.344277)


temp_has_74_t_2 <-     35
rate_has_74_t_2 <-  0.09622


# Combine the data set in one row 
temp <- c(temp_ches_17_2,temp_joubert_86_g,  temp_has_74_t_2)
rate <- c(rate_ches_17_2, rate_joubert_86_g,  rate_has_74_t_2)

# keep just a single curve
d_2 <- data.frame(temp, rate)


# fit Sharpe-Schoolfield model to raw data
d_fit <- nest(d_2, data = c(temp, rate)) %>%
  mutate(sharpeschoollow_1981 = map(data, ~nls_multstart(rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 15),
                                                         data = .x,
                                                         iter = c(4,4,4,4),
                                                         start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981') - 10,
                                                         start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981') + 10,
                                                         lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981'),
                                                         upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981'),
                                                         supp_errors = 'Y',
                                                         convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp),  by = 0.1))),
         # predict over that data,
         preds =  map2(sharpeschoollow_1981, new_data, ~augment(.x, newdata = .y)))


# refit model using nlsLM
fit_nlsLM2 <- nlsLM(rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 15),
                    data = d_2,
                    start = coef(d_fit$sharpeschoollow_1981[[1]]),
                    lower = get_lower_lims(d_2$temp, d_2$rate, model_name = 'sharpeschoollow_1981'),
                    upper = get_upper_lims(d_2$temp, d_2$rate, model_name = 'sharpeschoollow_1981'),
                    control = nls.lm.control(maxiter=500),
                    weights = rep(1, times = nrow(d_2)))


# bootstrap using residual resampling
boot3 <- Boot(fit_nlsLM2, method = 'residual')

# predict over new data
boot3_preds_2 <- boot3$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d_2$temp), max(d_2$temp), by = 0.1))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 15))

# create an empty row to combine the results of both bootstrap 
boot_preds_mu <- data.frame(matrix(0, ncol = 2, nrow = length(boot1_preds_beta_h$temp)))
colnames(boot_preds_mu) = c("temp","pred") 

# Set the division factor for the bootstrap 
temp_row_length <- length(seq(min_temp, max_temp, by = 0.1))

# set the first bootstrap end 
temp_row_length_1 <- length(seq(min(d$temp), max_temp, by = 0.1))

# set the second bootstrap start 
temp_row_length_2 <- length(seq(min(d_2$temp), max(d_2$temp), by = 0.1))

# do for loop to combine the bootstrap 
for (k in 1:(length(boot1_preds_beta_h$temp)/temp_row_length)){
  
  temp_boot_1 <- boot1_preds_1[(1+(k-1)*temp_row_length_1):(k*temp_row_length_1),]
  temp_boot_2 <- boot3_preds_2[(1+(k-1)*temp_row_length_2):(k*temp_row_length_2),]
  
  boot_preds_mu[(1+(k-1)*temp_row_length):(k*temp_row_length), 1:2] <- rbind(temp_boot_1[which(temp_boot_1$temp <= 37 & temp_boot_1$temp > min_temp-0.1), 5:6], temp_boot_2[which(temp_boot_2$temp > 2*max_temp & temp_boot_2$temp <= max_temp), 6:7])
  
}


################################ RUN FOR MEAN PARASITE BURDEN ################
# The number of year to run the model. 
year <- 30

# set a time span run
run_time <- seq(from = 0, to = 366*year-1, by = 1)

# Set a sample space for temperature 
sample_parameters <- seq(from = 0, to = 366*year, by = 1)


# Min temperature, (since we want to have a nice shape we pick up a small number)  
min_temp <- 12
# Max temperature, 
max_temp <- 37
# set temperature
temperature <- seq(min_temp, max_temp, by = 0.1)

# set the number of bootstrap 
num_boot <- 1

# upload necessary package for parallel code 
library(foreach)
library(doParallel)

#Setup backend to use many processors
totalCores = detectCores()

#Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1]-1) 
registerDoParallel(cluster)

# start for loop for each bootstrap 
boot_out_puts <- foreach (j = 1:num_boot, .combine = 'cbind') %dopar% {
  
  #Call the library 
  library(deSolve)
  # Set a matrix to record outputs 
  out_puts <- matrix(0, nrow = 6, ncol = length(temperature))
  
  # run the model for each temperature
  for (i in 1:length(temperature)){
    
    
    # This is mu_i function 
    preds_mu_i<- boot_preds_mu_i$pred[(1+(j-1)*length(temperature)):(length(temperature)*j)]
    
   
    # This is mu function 
    preds_mu <- boot_preds_mu$pred[(1+(j-1)*length(temperature)):(length(temperature)*j)]
    
    # Generate linearly interpolate point with temperature dependent parameter function 
    nu_s_afun <- approxfun(x = sample_parameters, y = rep(boot1_preds_nu_s$pred[(1+(j-1)*length(temperature)):(length(temperature)*j)][i], length(sample_parameters)))  
    mu_m_afun <- approxfun(x = sample_parameters, y =   rep(boot3_preds_mu_m$pred[(1+(j-1)*length(temperature)):(length(temperature)*j)][i], length(sample_parameters)))  
    mu_afun <- approxfun(x = sample_parameters, y =   rep(preds_mu[i], length(sample_parameters)))  
    sigma_s_afun <- approxfun(x = sample_parameters, y =   rep(boot1_preds_sigma_s$pred[(1+(j-1)*length(temperature)):(length(temperature)*j)][i], length(sample_parameters)))  
    mu_i_afun <- approxfun(x = sample_parameters, y =   rep(preds_mu_i[i], length(sample_parameters)))  
    nu_c_afun <- approxfun(x = sample_parameters, y =   rep(boot1_preds_nu_c$pred[(1+(j-1)*length(temperature)):(length(temperature)*j)][i], length(sample_parameters)))  
    mu_c_afun <- approxfun(x = sample_parameters, y =   rep(boot3_preds_mu_c$pred[(1+(j-1)*length(temperature)):(length(temperature)*j)][i], length(sample_parameters)))  
    delta_e_afun <- approxfun(x = sample_parameters, y =   rep(boot1_preds_delta_e$pred[(1+(j-1)*length(temperature)):(length(temperature)*j)][i], length(sample_parameters)))  
    beta_s_afun <- approxfun(x = sample_parameters, y =   rep(boot3_preds_beta_s$pred[(1+(j-1)*length(temperature)):(length(temperature)*j)][i] * contact_reduction_s, length(sample_parameters)))  
    beta_h_afun <- approxfun(x = sample_parameters, y =   rep(boot1_preds_beta_h$pred[(1+(j-1)*length(temperature)):(length(temperature)*j)][i] * contact_reduction_h, length(sample_parameters)))  
    
    
    # the diff equation solver 
    thermal_sensitive_model <- function(t, y, parms){
      with(as.list(c(y, parms)),
           {
             nu_s <- nu_s_afun(t)
             mu_m <- mu_m_afun(t)
             mu <- mu_afun(t) 
             sigma_s <- sigma_s_afun(t)
             mu_i <- mu_i_afun(t)    
             nu_c <- nu_c_afun(t)             
             mu_c <- mu_c_afun(t)
             delta_e <- delta_e_afun(t)
             beta_s <- beta_s_afun(t)
             beta_h <- beta_h_afun(t)
             
             dS <- (nu_s - (S + E + I) * nu) * (S + r * E) -  lambda * (1-exp(-(beta_s * delta_e * nu_e * h * P_m)/(mu_m * (S + E + I)))) * S - mu * S
             dE <- lambda * (1-exp(-(beta_s * delta_e * nu_e *h * P_m)/(mu_m * (S + E + I)))) * S - (sigma_s + mu_i) * E
             dI <- sigma_s * E - mu_i * I
             dP <-  beta_h * (nu_c/mu_c) * I - sigma_p * P
             dP_m <- sigma_p * P - (mu_h + mu_p) * P_m
             
             return(list(c(S = dS, E = dE, I = dI, P = dP, P_m = dP_m)))
           })
    }
    
    # Specified the parameter value. 
    parms0 <- c(nu, lambda)
    
    y0 <- c(S = (60434 + 33232)/2, E = 1260, I = 2340, P = 2, P_m = 105)
    
    ## solve the system 
    model_outputs <- ode(y = y0, times = run_time, func = thermal_sensitive_model, parms = parms0)
    
    # record the outcome 
    out_puts[1:5,i] <- model_outputs[length(model_outputs[, 1]), 2:6]
    out_puts[6,i] <-  model_outputs[length(model_outputs[, 1]),4]/sum(model_outputs[length(model_outputs[, 1]),2:4])
  }
     out_puts[5, ]
}
#Stop cluster
stopCluster(cluster)

# save as csv
write.csv(data.frame(boot_out_puts), "haematobium_mpb_boot_out_puts.csv", row.names=FALSE)

