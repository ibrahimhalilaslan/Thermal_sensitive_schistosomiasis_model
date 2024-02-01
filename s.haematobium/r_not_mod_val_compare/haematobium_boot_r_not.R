# load packages
library(boot)
library(car)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)

source("par_set_haematobium.R")

min_temp <- 12
max_temp <- 37

# Set the number of row to bootstrap to remove 
n <- 25

#(Purnell 1966),, S. mansoni
purnell_66_temp <- c(       12,             15,           18,             21,        24,          27,           30,       33)
purnell_66_rate_1 <- c(  0.9246078,      0.810981,    1.320736,       1.092269,   0.6998873,   1.231836,    0.7901292,  0.6292689)
purnell_66_rate_1 <- purnell_66_rate_1/max(purnell_66_rate_1)


purnell_66_temp <- c(       12,             15,           18,             21,        24,          27,           30,       33)
purnell_66_rate_2 <- c(  1.0535,         1.002849,     0.9284408,       1.481275,   0.7225294,   0.6554545,    1.251547,  0.7475422)
purnell_66_rate_2 <- purnell_66_rate_2/max(purnell_66_rate_2)


purnell_66_temp <- c(       12,             15,           18,             21,        24,          27,           30,       33)
purnell_66_rate_3 <- c(  2.182979,      1.523512,    1.777458,       1.840471,   0.8305635,   0.7054314,      1.339291,  0)
purnell_66_rate_3 <- purnell_66_rate_3/max(purnell_66_rate_3)



#purnell_66_rate <- purnell_66_rate/max(purnell_66_rate)


#(Stirewalt 1954), S. mansoni 
stirewalt_54_temp <- c(    24,         27)
stirewalt_54_rate <- c(0.09590499,  0.1552625)
stirewalt_54_rate <- stirewalt_54_rate/max(stirewalt_54_rate)


#(DeWitt 1965) with S. mansoni
dewitt_55_temp_1 <- c(     15,      20,          25,      30,      35)
dewitt_55_rate_1 <- c(0.4618475, 0.8464159,  1.00836, 1.237522, 0.6192221)
dewitt_55_rate_1 <- dewitt_55_rate_1/max(dewitt_55_rate_1)

dewitt_55_temp_2 <- c(     15,      20,          25,      30,         35,        40)
dewitt_55_rate_2 <- c(0.385854, 0.6910937,    0.635236,  0.9783329, 1.0652, 0.2764123)
dewitt_55_rate_2 <- dewitt_55_rate_2/max(dewitt_55_rate_2)

dewitt_55_temp_3 <- c(   10,        15,      20,         25,      30,       35,        40,   45)
dewitt_55_rate_3 <- c(0.7439515,  1.76077,  2.484149,  3.50177,  3.9655, 3.421058,  0.7532138, 0)
dewitt_55_rate_3 <- dewitt_55_rate_3/max(dewitt_55_rate_3)

dewitt_55_temp_4 <- c(0, 5,        10,        15,      20,        25)
dewitt_55_rate_4 <- c(0, 0,    0.1190698, 0.6479505, 1.39834, 1.50036)
dewitt_55_rate_4 <- dewitt_55_rate_4/max(dewitt_55_rate_4)

# load data
#The transmission rate of schisto in human
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

boot1_preds_beta_h <- head(boot1_preds, -n*length(seq(min_temp, max_temp, by = 0.1)))

################################ TRANSMSISION RATE IN SNAIL 


#(Prah and James 1977), Bulinus (Physopsis) globosus with S. haematobium.
prah_77_temp <- c(       5,        12,            19,     26)
prah_77_rate <- c(0.001279126, 0.04565688,  0.1084776, 0.1849867)
prah_77_rate <- prah_77_rate/max(prah_77_rate)

#(Shiff 1974), Bulinus globosus infected with S. haematobium.
shiff_74_temp <- c(      25,          20,           18,         16,          13,            11) 
shiff_74_rate <- c(0.03730699,  0.003512327,  0.01022082, 0.001782906,  0.002553961,  0.001281865)
shiff_74_rate <- shiff_74_rate/max(shiff_74_rate)

#Chu et al. 1966), Bulinus trancatus  infected with S. haematobium. 
chu_66_temp <- c(    10,            12,           14,         15,        20,          25,         30,        35,         38)
chu_66_rate <- c(0.002009499,  0.008145561,  0.02029471, 0.02016721, 0.04754069, 0.03577785,  0.06475144, 0.04395093, 0.05660857)
chu_66_rate <- chu_66_rate/max(chu_66_rate)           


#The disease transmission rate in snails (Anderson et al. 1982, experiment with Biomphalaria glabra and schistosoma mansoni)
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

boot3_preds_beta_s <- head(boot3_preds, -n*length(seq(min_temp, max_temp, by = 0.1)))


############################## CERCARIAL RELEASE RATE #################

#Nguyen et al. 2021, S. mansoni with Biomphalaria glabrata
temp_nguyen_21 <- c(5,  9,     13,     17,     21,     25,     29,     33,     37)
rate_nguyen_21 <- c(0, 15,     260,   1350,   1970,   2900,   1970,   1770,   690)

#(Upatham 1973), B. glabrata with S. mansoni  
temp_upatham_73 <- c(10, 13, 16,   19, 22,    25,   28,  31,   34,   37, 40)
rate_upatham_73 <- c( 0,  0, 206, 114, 442, 1200, 1400, 1233, 1241, 243, 121)


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

boot1_preds_nu_c <- head(boot1_preds, -n*length(seq(min_temp, max_temp, by = 0.1)))


########################## FECUNDITY RATE OF SNAIL ######################

#Chester Kalinda. 2017, experiment with Bulinus globosus. 
#(Note Age of maturity and hatching rate getting from El-Hassan 1974 for buliunus)
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

boot1_preds_sigma_s <- head(boot1_preds, -n*length(seq(min_temp, max_temp, by = 0.1)))


############################# MIRACIDIA HATCH RATE ############## 

# load in data
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

boot1_preds_delta_e <- head(boot1_preds, -n*length(seq(min_temp, max_temp, by = 0.1)))


##################### MORTALITY RATE OF CERCARIA ##################


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

# Data 
temp <- c(temp_ches_17_2,temp_joubert_86_g, temp_joubert_86_a, temp_has_74_t_2)
rate <- c(rate_ches_17_2, rate_joubert_86_g, rate_joubert_86_a, rate_has_74_t_2)

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


boot_preds_mu_i <- data.frame(matrix(0, ncol = 2, nrow = length(boot1_preds_beta_h$temp)))
colnames(boot_preds_mu_i) = c("temp","pred") 


temp_row_length <- length(seq(min_temp, max_temp, by = 0.1))

temp_row_length_1 <- length(seq(min(d$temp), max(d$temp), by = 0.1))

temp_row_length_2 <- length(seq(min(d_2$temp), max(d_2$temp), by = 0.1))


for (k in 1:(length(boot1_preds_beta_h$temp)/temp_row_length)){
  
  temp_boot_1 <- boot3_preds_1[(1+(k-1)*temp_row_length_1):(k*temp_row_length_1),]
  temp_boot_2 <- boot3_preds_2[(1+(k-1)*temp_row_length_2):(k*temp_row_length_2),]
  
  boot_preds_mu_i[(1+(k-1)*temp_row_length):(k*temp_row_length), 1:2] <- rbind(temp_boot_1[which(temp_boot_1$temp <= max_temp & temp_boot_1$temp > min_temp-0.1), 6:7], temp_boot_2[which(temp_boot_2$temp > 2*max_temp & temp_boot_2$temp <= max_temp), 6:7])
  
}



################################ MORTALITIY OF MIRACIDIA ################
#load data 
#S. K. Prah and C. James, 1977 S. haematobium

temp_prah_77_haem <- c(7.5,              20,     27.5,       36.5)
rate_prah_77_haem <- c(2.275331,      2.210028,  2.24239,  4.354702)


temp <- temp_prah_77_haem
rate <- rate_prah_77_haem


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


boot3_preds_mu_m <- head(boot3_preds, -n*length(seq(min_temp, max_temp, by = 0.1)))

############################ MORTALITY OF SNAIL ##############################


# This analysis are up to 31 degree temperature. 

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


# Data 
temp <- c(temp_ches_17_2,temp_joubert_86_g,  temp_has_74_t_2)
rate <- c(rate_ches_17_2, rate_joubert_86_g,  rate_has_74_t_2)

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



boot_preds_mu <- data.frame(matrix(0, ncol = 2, nrow = length(boot1_preds_beta_h$temp)))
colnames(boot_preds_mu) = c("temp","pred") 

temp_row_length <- length(seq(min_temp, max_temp, by = 0.1))

temp_row_length_1 <- length(seq(min(d$temp), max_temp, by = 0.1))

temp_row_length_2 <- length(seq(min(d_2$temp), max(d_2$temp), by = 0.1))


for (k in 1:(length(boot1_preds_beta_h$temp)/temp_row_length)){
  
  temp_boot_1 <- boot1_preds_1[(1+(k-1)*temp_row_length_1):(k*temp_row_length_1),]
  temp_boot_2 <- boot3_preds_2[(1+(k-1)*temp_row_length_2):(k*temp_row_length_2),]
  
  boot_preds_mu[(1+(k-1)*temp_row_length):(k*temp_row_length), 1:2] <- rbind(temp_boot_1[which(temp_boot_1$temp <= 37 & temp_boot_1$temp > min_temp-0.1), 5:6], temp_boot_2[which(temp_boot_2$temp > 2*max_temp & temp_boot_2$temp <= max_temp), 6:7])
  
}



########################## CALCULATE R_NOT ##################

# This is analytic representation of R not value without control. 
r_not_boot <- (lambda * (boot3_preds_beta_s$pred * contact_reduction_s) * h * boot1_preds_delta_e$pred * nu_e * (boot1_preds_beta_h$pred * contact_reduction_h) 
               * boot1_preds_nu_c$pred * boot1_preds_sigma_s$pred/
                 (boot3_preds_mu_m$pred * (mu_h + mu_p) * boot3_preds_mu_c$pred * boot_preds_mu_i$pred 
                  * (boot1_preds_sigma_s$pred + boot_preds_mu_i$pred)))^(1/2)      


temp <- boot1_preds_beta_h$temp
boot_out_puts <- data.frame(temp, r_not_boot)

# calculate bootstrapped confidence intervals
boot_conf_preds <- boot_out_puts %>% group_by(temp) %>%
  summarise(conf_lower = quantile(r_not_boot, 0.025, na.rm = TRUE),
            conf_upper = quantile(r_not_boot, 0.975, na.rm = TRUE)) %>%
  ungroup()

optim_temps <- c()


temp_row_length <- length(seq(min_temp, max_temp, by = 0.1))


for (i in 1:(length(boot1_preds_beta_h$temp)/temp_row_length)){
  index <-  which.max(r_not_boot[(1+(i-1)*temp_row_length):(temp_row_length*i)])
  optim_temps[i] <- boot_out_puts$temp[index+(i-1)*temp_row_length] 
}


lower_bound <- quantile(optim_temps, 0.025, na.rm = TRUE)
upper_bound <- quantile(optim_temps, 0.975, na.rm = TRUE)

lower_bound
# 23.6
upper_bound
# 27.9


min_temp <- 12
max_temp <- 37

# Generate a sequence of temperature    
temperature <- data.frame(temp = seq(min_temp, max_temp, by = 0.1))

# combine two piece wise function of mu_i
preds_1 <- fn_mu_i_1(temperature)$.fitted[which(temperature$temp <= 37)]
preds_2 <- fn_mu_i_2(temperature)$.fitted[which(temperature$temp > 37)]

# This is mu_i function 
preds_mu_i <- c(preds_1, preds_2)  

# This is analytic representation of R not value without control. 
r_not_best <- (abs(lambda * fn_beta_s(temperature)$.fitted * h * fn_delta_e(temperature)$.fitted * nu_e 
                   * fn_beta_h(temperature)$.fitted * fn_nu_c(temperature)$.fitted * fn_sigma_s(temperature)$.fitted/
                     (fn_mu_m(temperature)$.fitted * (mu_h + mu_p) * fn_mu_c(temperature)$.fitted * preds_mu_i * 
                        (fn_sigma_s(temperature)$.fitted + preds_mu_i))))^(1/2)      


r_not_data <- data.frame(temp, r_not_best)


# We actually do this to test Nguyen results with our model, we do simulation with increasing some of
# our parameter values as disease control strategies. 
prev_percent <- read.csv(file = 'gntd_vars_all.csv')
head(prev_percent)

sch_heamatobium_ind <- which(prev_percent$parasite_s == "S. haematobium")


sch_heamatobium <- prev_percent[sch_heamatobium_ind, ]
sch_heamatobium <- sch_heamatobium[-which(sch_heamatobium$percent_pos == 0),]
sch_heamatobium <- sch_heamatobium[-c(which(is.na(sch_heamatobium$bio01))),]
sch_heamatobium_sort <- sch_heamatobium[order(sch_heamatobium$bio01), ] 


sum(is.na(sch_heamatobium_sort$bio01))
length(sch_heamatobium_sort$bio01)

sch_heamatobium_temp_unq <- sch_heamatobium %>% group_by(bio01) %>% 
  summarize(med = median(percent_pos), na.rm = TRUE)


# plot bootstrapped CIs
pdf(file = "heamatobium_data_median.pdf", width = 5, height = 5)
ggplot() +
  geom_point(aes(bio01, med), sch_heamatobium_temp_unq, col = 'blue') +
  labs(title = "S.heamatobium", 
       x = 'Temperature (ºC)',
       y = expression("R"[0]))+
  labs(title = "S.haematobium", x = 'Temperature (ºC)',y = "Percent of positive case")+
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 17),
        axis.text.x = element_text(color="black", size=17), axis.text.y = element_text(color="black", size=17),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=14, face="bold.italic", hjust=0.5))
dev.off() 


bins_maker <- function(percentile, lenth_of_bin){
  
  out_put <- matrix(0, nrow = 2, ncol = round(length(sch_heamatobium_sort$bio01)/lenth_of_bin))
  
  for (i in 1:round(length(sch_heamatobium_sort$bio01)/lenth_of_bin)){
    index_of_bin <- (1+(i-1)*lenth_of_bin):(i*lenth_of_bin)
    out_put[1, i] <- sch_heamatobium_sort$bio01[index_of_bin[lenth_of_bin/2]]
    out_put[2, i] <- quantile(sch_heamatobium_sort$percent_pos[index_of_bin], na.rm = TRUE, probs = percentile)
  }
  return(out_put)
}

bins_data <-bins_maker(percentile = .95, lenth_of_bin = 400)
d <- data.frame(bins_data[1,], bins_data[2,])
colnames(d) <- c("temp", "prev")

# plot bootstrapped CIs
pdf(file = "haematobium_boost_r_not.pdf", width = 5, height = 5)
ggplot() +
  geom_line(aes(temp, r_not_best), r_not_data, col = 'blue', size = 1) +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot_conf_preds, fill = 'blue', alpha = 0.1) +
  geom_point(aes(temp, prev*max(r_not_best)/max(prev)), d, size = 1, alpha = 0.5, colour = "blue") + 
  labs(title = "S.haematobium", x = 'Temperature (ºC)',y = expression("R"[0]))+
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 17),
        axis.text.x = element_text(color="black", size=17), axis.text.y = element_text(color="black", size=17),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=14, face="bold.italic", hjust=0.5))
dev.off() 
         

############# Nguyen's result plot ######################
temp <- c(15.127, 15.763, 16.525, 16.949, 17.458, 17.881, 18.305, 18.814, 19.322, 19.831, 20.508, 21.186, 21.949, 22.797, 23.644, 24.153, 24.746,  25.085, 25.593,  26.102,  26.61,  27.119,  27.542,  28.136,  28.898, 29.576, 30.339,  31.441, 32.458, 33.814, 34.619)
rate <- c(0.11, 0.172, 0.282, 0.373, 0.467, 0.533, 0.621,  0.712, 0.806, 0.884, 0.953, 0.984, 0.994,  0.95, 0.84, 0.765, 0.665, 0.586,  0.502, 0.439, 0.376, 0.301, 0.254, 0.197, 0.144, 0.11, 0.072, 0.038, 0.025, 0.009, 0.009)

# keep just a single curve
.x <- data.frame(temp, rate)

nguyen <- nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                        data = .x,
                        iter = c(4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                        supp_errors = 'Y')

fn_nguyen <- function(x) augment(nguyen, newdata = x)


df <- data.frame(x = r_not_data$temp, y2 = r_not_data$r_not_best, y1 = fn_nguyen(temperature)$.fitted)

#expression("R"[0]/max(("R"[0])))

# plot bootstrapped CIs
pdf(file = "haematobium_compare_boost_r_not.pdf", width = 5, height = 5)
ggplot(data = df, aes(x = x)) +
  geom_line(aes(y = y1/max(y1), colour = "Prev. Est."),size=1) +
  geom_line(aes(y = y2/max(r_not_best), colour = "New Est."), size=1) +
  scale_colour_manual("", 
                      breaks = c("Prev. Est.", "New Est."),
                      values = c("grey", "blue")) +
  geom_point(aes(temp, prev/max(prev)), d, size = 1, alpha = 0.5, colour = "blue") + 
  geom_ribbon(aes(temp, ymin = conf_lower/max(r_not_best), ymax = conf_upper/max(r_not_best)), boot_conf_preds, fill = 'blue', alpha = 0.1) +  
  labs(title = "S.haematobium", x = '',y = "")+  
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 17),
        axis.text.x = element_text(color="black", size=17), axis.text.y = element_text(color="black", size=17),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=14, face="bold.italic", hjust=0.5)) +
  #Temperature (ºC)
  # Add mark segment
  annotate("text", x = temperature$temp[which.max(df$y1)], y = 0, label = '^', colour = "grey", size = 5) + 
  annotate("text", x = temperature$temp[which.max(df$y1)], y = 0.08,  label = round(df$x[which.max(df$y1)],1), colour = "grey", size = 5) +
  annotate("text", x = temperature$temp[which.max(df$y2)], y = 0, label = '^', colour = "blue", size = 5) +
  annotate("text", x = temperature$temp[which.max(df$y2)], y = 0.08, label = round(df$x[which.max(df$y2)],1), colour = "blue", size = 5) +
  # Add horizontal line segment
   geom_segment(aes(x = lower_bound, y = -0.02, xend = upper_bound, yend = -0.02), colour = "blue", size = 1)
dev.off()




#Calculate the mean sum square to compare the both ,model output quantitatively

new_temp_data_mss <- data.frame(temp = d$temp)

# combine two piece wise function of mu_i
preds_1 <- fn_mu_i_1(new_temp_data_mss)$.fitted[which(new_temp_data_mss$temp <= 37)]
preds_2 <- fn_mu_i_2(new_temp_data_mss)$.fitted[which(new_temp_data_mss$temp > 37)]

# This is mu_i function 
preds_mu_i <- c(preds_1, preds_2)

r_not_best_mss <- (abs(lambda * fn_beta_s(new_temp_data_mss)$.fitted * h * fn_delta_e(new_temp_data_mss)$.fitted * nu_e 
                       * fn_beta_h(new_temp_data_mss)$.fitted * fn_nu_c(new_temp_data_mss)$.fitted * fn_sigma_s(new_temp_data_mss)$.fitted/
                         (fn_mu_m(new_temp_data_mss)$.fitted * (mu_h + mu_p) * fn_mu_c(new_temp_data_mss)$.fitted * preds_mu_i * 
                            (fn_sigma_s(new_temp_data_mss)$.fitted + preds_mu_i))))^(1/2)   


nguyen_mss <- (fn_nguyen(new_temp_data_mss)$.fitted/max(fn_nguyen(new_temp_data_mss)$.fitted))*max(r_not_best_mss)

data_nss <- (d$prev/max(d$prev))*max(r_not_best_mss)

new_result <- sum((r_not_best_mss - data_nss)^2)/length(data_nss)
prev_result <- sum((nguyen_mss - data_nss)^2)/length(data_nss)






