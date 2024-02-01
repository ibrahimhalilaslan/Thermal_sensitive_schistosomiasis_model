############################# FIND THE CONFIDENT INTERVAL AND MPB ###############
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



############### SET CONSTANT PARAMETERS  #######################
min_temp <- 12
max_temp <- 37

#The factor of density-dependent fecundity rate
# We adjust this value to get the abundance of sanil in a site same as C. Wood at al. 2019 paper
nu <-7.595228e-07 # (~ 450000 snails)

#The maximal rate of snail invasion by miracidia
# We adjust this value to get the the prevalence of disease in snail to be 5 percent
lambda <- 0.001724518

#Constant human population having schistosomiasis 
# We adjust this value to get the number of human in a site similar as C. Wood at al. 2019 paper also Susanne H. Sokolow, 2015 has similar value
h <- 1000 #(human)

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

# set the reduction rate for transmission rate in human 
contact_reduction_h <- 1.557485e-08


# set the reduction rate for transmission rate in snails 
contact_reduction_s <- 0.8

####### START BOOTSTRAPPING 
# Set the number of row to bootstrap to remove becuase 
# we do not have sufficient data for haematobium and the bootstrap give error last a few bootstrap
n <- 25

############################# TRANSMSISION RATE IN HUMAN #################

#(Purnell 1966),, S. mansoni
purnell_66_temp <- c(       12,             15,           18,             21,        24,          27,           30,       33)
purnell_66_rate_1 <- c(  0.9246078,      0.810981,    1.320736,       1.092269,   0.6998873,   1.231836,    0.7901292,  0.6292689)
#normalize the data with maximum value 
purnell_66_rate_1 <- purnell_66_rate_1/max(purnell_66_rate_1)


purnell_66_temp <- c(       12,             15,           18,             21,        24,          27,           30,       33)
purnell_66_rate_2 <- c(  1.0535,         1.002849,     0.9284408,       1.481275,   0.7225294,   0.6554545,    1.251547,  0.7475422)
#normalize the data with maximum value 
purnell_66_rate_2 <- purnell_66_rate_2/max(purnell_66_rate_2)


purnell_66_temp <- c(       12,             15,           18,             21,        24,          27,           30,       33)
purnell_66_rate_3 <- c(  2.182979,      1.523512,    1.777458,       1.840471,   0.8305635,   0.7054314,      1.339291,  0)
#normalize the data with maximum value 
purnell_66_rate_3 <- purnell_66_rate_3/max(purnell_66_rate_3)

#(Stirewalt 1954), S. mansoni 
stirewalt_54_temp <- c(    24,         27)
stirewalt_54_rate <- c(0.09590499,  0.1552625)
#normalize the data with maximum value 
stirewalt_54_rate <- stirewalt_54_rate/max(stirewalt_54_rate)


#(DeWitt 1965) with S. mansoni
dewitt_55_temp_1 <- c(     15,      20,          25,      30,      35)
dewitt_55_rate_1 <- c(0.4618475, 0.8464159,  1.00836, 1.237522, 0.6192221)
#normalize the data with maximum value 
dewitt_55_rate_1 <- dewitt_55_rate_1/max(dewitt_55_rate_1)

dewitt_55_temp_2 <- c(     15,      20,          25,      30,         35,        40)
dewitt_55_rate_2 <- c(0.385854, 0.6910937,    0.635236,  0.9783329, 1.0652, 0.2764123)
#normalize the data with maximum value 
dewitt_55_rate_2 <- dewitt_55_rate_2/max(dewitt_55_rate_2)

dewitt_55_temp_3 <- c(   10,        15,      20,         25,      30,       35,        40,   45)
dewitt_55_rate_3 <- c(0.7439515,  1.76077,  2.484149,  3.50177,  3.9655, 3.421058,  0.7532138, 0)
#normalize the data with maximum value 
dewitt_55_rate_3 <- dewitt_55_rate_3/max(dewitt_55_rate_3)

dewitt_55_temp_4 <- c(0, 5,        10,        15,      20,        25)
dewitt_55_rate_4 <- c(0, 0,    0.1190698, 0.6479505, 1.39834, 1.50036)
#normalize the data with maximum value 
dewitt_55_rate_4 <- dewitt_55_rate_4/max(dewitt_55_rate_4) 

# combine the data one row 
temp <- c(purnell_66_temp, purnell_66_temp, purnell_66_temp,stirewalt_54_temp, dewitt_55_temp_1, dewitt_55_temp_2, dewitt_55_temp_3, dewitt_55_temp_4)
rate <- c(purnell_66_rate_1, purnell_66_rate_2, purnell_66_rate_3, stirewalt_54_rate, dewitt_55_rate_1, dewitt_55_rate_2, dewitt_55_rate_3, dewitt_55_rate_4)

# keep just a single curve
d <- data.frame(temp, rate)

# fit Sharpe-Schoolfield model
d_fit <- nest(d, data = c(temp, rate)) %>%
  mutate(briere2_1999 = map(data, ~nls_multstart(rate~briere2_1999(temp = temp, tmin, tmax, a,b),
                                                 data = .x,
                                                 iter = c(4,4,4,4),
                                                 start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') - 10,
                                                 start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') + 10,
                                                 lower = get_lower_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                                                 upper = get_upper_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
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

# (Stirewalt 1954) Australorbis glabratus with S. mansoni.
stirewalt_54_temp <- c(    24,         27,        32)
stirewalt_54_rate <- c(0.07515584, 0.4808139, 0.3922359)
#normalize the data with maximum value 
stirewalt_54_rate <- stirewalt_54_rate/max(stirewalt_54_rate)

#(Foster 1964), Biomphalaria pfeifferi with S. mansoni.
foster_64_temp <- c(    22.85,      24.01,       26.26,      28.07,       30.04,    31.75)
foster_64_rate <- c(0.04079299,   0.04333207,  0.2203394, 0.08210066,  0.05756868, 0.07958718)
#normalize the data with maximum value 
foster_64_rate <- foster_64_rate/max(foster_64_rate)

#(Prah and James 1977),  Biomphalaria pfeifferi with S. mansoni. 
prah_77_temp <- c(    5,             15,         20,         22.5)
prah_77_rate <- c(0.0006344901, 0.05059652, 0.1629059,  0.2249733)
#normalize the data with maximum value 
prah_77_rate <- prah_77_rate/max(prah_77_rate)

#(Anderson et al. 1982), Biomphalaria glabrata with S. mansoni 
anderson_82_temp <- c(15,               20,          25,           30,         35)
anderson_82_rate <- c(0.3658224,    0.3865833,    0.9074138,     0.487644, 0.4107836)
#normalize the data with maximum value 
anderson_82_rate <-  anderson_82_rate/max(anderson_82_rate)

#(DeWitt 1955), B. glabrata with S. mansoni. 
dewitt_55_temp <-   c(10,      25,       35,    40)
dewitt_55_rate_1 <- c(0,   0.0078523, 0.0283053, 0)
#normalize the data with maximum value 
dewitt_55_rate_1 <- dewitt_55_rate_1/max(dewitt_55_rate_1)
dewitt_55_rate_2 <- c(0,   0.032993,  0.0241807, 0)
#normalize the data with maximum value 
dewitt_55_rate_2 <- dewitt_55_rate_2/max(dewitt_55_rate_2)
dewitt_55_rate_3 <- c(0,   0.01388562,  0.03692433, 0.01246124)
#normalize the data with maximum value 
dewitt_55_rate_3 <- dewitt_55_rate_3/max(dewitt_55_rate_3)
dewitt_55_rate_4 <- c(0,   0.07405795,  0.09039597, 0.00401955)
#normalize the data with maximum value 
dewitt_55_rate_4 <- dewitt_55_rate_4/max(dewitt_55_rate_4)

#Coelho and Bezerra 2006), Biomphalari glabrata.
coelho_temp_06 <- c(15, 20, 30)
coelho_rate_06 <- c(0.0001342012, 0.007285127, 0.04153248)
#normalize the data with maximum value 
coelho_rate_06 <- coelho_rate_06/max(coelho_rate_06)

#(Upatham 1973), B. glabrata with S. mansoni. 
upatham_temp_73 <- c(10,      13,     16,           19,           22,           25,        28,          31,          34,          37,          40)
upatham_rate_73 <- c(0,        0,  0.002895729, 0.004360468,  0.009807099, 0.01299096, 0.01373957,  0.01528639,  0.02722317, 0.002743025,  0.00199986)
#normalize the data with maximum value 
upatham_rate_73 <- upatham_rate_73/max(upatham_rate_73)

#combine the data set in one row 
temp <- c(stirewalt_54_temp, foster_64_temp, prah_77_temp,  anderson_82_temp, dewitt_55_temp,   dewitt_55_temp,   dewitt_55_temp,   dewitt_55_temp,  coelho_temp_06, upatham_temp_73)
rate <- c(stirewalt_54_rate, foster_64_rate, prah_77_rate,  anderson_82_rate, dewitt_55_rate_1, dewitt_55_rate_2, dewitt_55_rate_3, dewitt_55_rate_4,  coelho_rate_06, upatham_rate_73)

d <- data.frame(temp, rate)

# fit Sharpe-Schoolfield model
d_fit <- nest(d, data = c(temp, rate)) %>%
  mutate(spain_1982 = map(data, ~nls_multstart(rate~spain_1982(temp = temp, a,b,c,r0),
                                               data = d,
                                               iter = c(4,4,4,4),
                                               start_lower = get_start_vals(d$temp, d$rate, model_name = 'spain_1982') - 1,
                                               start_upper = get_start_vals(d$temp, d$rate, model_name = 'spain_1982') + 1,
                                               lower = get_lower_lims(d$temp, d$rate, model_name = 'spain_1982'),
                                               upper = get_upper_lims(d$temp, d$rate, model_name = 'spain_1982'),
                                               supp_errors = 'Y',
                                               convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), by = 0.1))),
         # predict over that data,
         preds =  map2(spain_1982, new_data, ~augment(.x, newdata = .y)))


# refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(rate~spain_1982(temp = temp, a,b,c,r0),
                               data = d,
                               start = coef(d_fit$spain_1982[[1]]),
                               lower = get_lower_lims(d$temp, d$rate, model_name = 'spain_1982'),
                               upper = get_upper_lims(d$temp, d$rate, model_name = 'spain_1982'),
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
  mutate(pred = spain_1982(temp = temp, a,b,c,r0))

# remove the bootstrap does not exist 
boot1_preds_beta_s <- head(boot1_preds, -n*length(seq(min_temp, max_temp, by = 0.1)))

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

# McCreesh et al. 2014, Biomphalaria sudanica.
temp_mccreesh_14 <- c(13.4, 15.7,  16.7,    18.9,    20.9,    22.8,   26.7,  28.3, 29.5, 32.0)
rate_mccreesh_14 <- c(  0,    0,  0.0349,  0.0369,  0.0395,  0.0386, 0.0499,  0,     0,   0)

#C. C. Appleton, 1977, experiment with Biomphalaria pfeifferi
temp_applet_77 <- c(25,       27,      29)
rate_applet_77 <- c(0.07058,  0.06824,    0)


#El- Hassan 1974, Biomphalaria Alexandrina
temp_has_74_a <- c(15,        20,      25,   30,    35)
rate_has_74_a <- c(0.00683, 0.0229, 0.0309, 0.00824, 0)


#R. F. Sturrock, 1966, experiment with Biomphalaria pfeifferi
temp_sturr_66 <- c(19,        25,        30,    35)
rate_sturr_66 <- c(0.02720, 0.07378,    0.06568, 0)

# Michelson 1961, Biomphalaria glabrata.
temp_michel_61 <- c(5, 15,   20,     25,     30,  35)
rate_michel_61 <- c(0, 0, 0.0132, 0.0349, 0.00969, 0)


# Shiff and Garnett 1963, Biomphalaria pfeifferi
temp_shiff_63 <- c(   18,     22,   25,     27)
rate_shiff_63 <- c(0.0138, 0.0507, 0.0349, 0.036)

#combine the data set in one row 
temp <- c(temp_mccreesh_14, temp_applet_77, temp_has_74_a, temp_sturr_66, temp_michel_61, temp_shiff_63)
rate <- c(rate_mccreesh_14, rate_applet_77, rate_has_74_a, rate_sturr_66, rate_michel_61, rate_shiff_63)

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
boot1_preds_nu_s <- head(boot1_preds, -n*length(seq(min_temp, max_temp, by = 0.1)))

################################# LATENT PERIOD OF SNAIL ##################


#W. Pflüger 1981, Biomphalaria glabrata with S. mansoni, diurnally changing water temperature 
temp_pfluger_81 <- c(24.3,     21.9,   20.2,  19.8,  17.5,   17.9,   19.2,  21.1, 18.9,    25,   29.9,  33)
rate_pfluger_81 <- c(1/26,    1/35.5,  1/44, 1/50.5, 1/68,  1/57.5,  1/52,  1/38, 1/57.5, 1/17, 1/16.5, 1/16.5)

#W. Pflüger 1980, Biomphalaria glabrata with S. mansoni, constant temperature 
temp_pfluger_80 <- c(  17,      19,     22,   25,      28,     30,      31,       32,     33,     34,    35)
rate_pfluger_80 <- c(1/92.5,  1/55.5,  1/34, 1/24.8, 1/19.2,  1/16.3,  1/15.9,  1/14.9, 1/15.3, 1/15.9, 1/16)


#Foster 1964, Biomphalaria pfeifleri with s. mansoni 
temp_foster_64 <- c(18,     21,    22.85,   24.01,   26.26,  28.07,  30.04, 31.75)
rate_foster_64 <- c(1/57,  1/37,    1/32,   1/30,     1/23,   1/19,   1/18,   1/16)


# R. M. Gordon, 1934, Planorbis pfeifferi with S. mansoni
temp_gordon_34_pfei <- c(21.7,    22,       26.3,      26.6,       26.9,     27.3,    27.7,    31.8,    32.1,    32.7,   32.8,   35)
rate_gordon_34_pfei <- c(1/36,   1/33,     1/23.5,     1/23,       1/22.2,   1/22,    1/19,   1/15.4, 1/15.8,  1/17.2,  1/14.7, 1/16.4)


#combine the data set in one row 
temp <- c(temp_pfluger_81,  temp_pfluger_80, temp_foster_64, temp_gordon_34_pfei)
rate <- c(rate_pfluger_81,  rate_pfluger_80, rate_foster_64, rate_gordon_34_pfei)


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
boot1_preds_sigma_s <- head(boot1_preds, -n*length(seq(min_temp, max_temp, by = 0.1)))


############################# MIRACIDIA HATCH RATE ############## 

#(Nguyen et al. 2021, Schistosoma mansoni)
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


#combine the data set in one row 
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


#W. Pflüger 1981, Biomphalaria glabrata with S. mansoni
temp_pfluger_81 <- c(24.3,            21.9,          20.2,       19.7)
rate_pfluger_81 <- c(0.02478296,   0.02157024,   0.02852877,  0.02473841)

#W. Pflüger 1981, Biomphalaria glabrata with S. mansoni with sinual function
temp_pfluger_81_sin <- c(17.5,           17.9,         19.2,            21.1,        18.9,      25,           29.9,     33.1)
rate_pfluger_81_sin <- c(0.00524522,  0.01841793,  0.00982357,   0.008281862,  0.02503713,  0.01368767, 0.02098088, 0.03903461)


#W. Pflüger 1980, Biomphalaria glabrata with S. mansoni, constant temperature (This result from prepatent snails)
temp_pfluger_80 <- c(16,               17,            18,              19,       22,         25,            28,      30,            31,           32,         33,        34,       35)
rate_pfluger_80 <- c(0.02354884,   0.01456296,  0.01980421,      0.01161004, 0.02098088,  0.01863046,  0.01136229, 0.01292767,  0.01194658,  0.009346447, 0.01923994, 0.07677861, 0.04147427)


#R. Foster, 1964 experiment with Biomphalaria pfeifferi invasion by S. mansoni
temp_foster_64 <- c(22.85,           24.01,        26.26,      28.07)
rate_foster_64 <- c(0.006674626,  0.01155245,    0.0260108,   0.03678792)

#(Upatham 1973), B. glabrata with S. mansoni  
temp_upatham_73 <- c(      10,          13,          16,        19,          22,          25,           28,          31,          34,          37,      40)
rate_upatham_73 <- c(0.006316809, 0.009430725, 0.006047558, 0.02049171, 0.009058285, 0.003967988, 0.0004177118, 0.006226962, 0.006226962, 0.09977736,  0.0916344)

#combine the data set in one row 
temp <- c(temp_pfluger_81, temp_pfluger_81_sin, temp_pfluger_80,  temp_foster_64, temp_upatham_73)
rate <- c(rate_pfluger_81, rate_pfluger_81_sin, rate_pfluger_80,  rate_foster_64, rate_upatham_73)

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
boot3_preds_1 <- boot3$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min_temp, max_temp, by = 0.1))) %>%
  ungroup() %>%
  mutate(pred = flinn_1991(temp = temp, a, b, c))

# we need to do second bootstrap for higher temperature 
# set the minumum and maximum temperature 
min_temp_2 <- 34
max_temp_2 <- 37

#P. H. Joubert, 1986, Biomphalaria pfeifferi 
temp_joubert_86_p <- c(     34,         36,       38,       40)
rate_joubert_86_p <- c(0.1457006, 0.6661235, 1.649401,  5.821733)

#El- Hassan 1974, Biomphalaria Alexandrina
temp_has_74_a_h <- c(35,    37)
rate_has_74_a_h <- c(0.03891, 0.3)

#combine the data set in one row 
temp <- c(temp_joubert_86_p,  temp_has_74_a_h) 
rate <- c(rate_joubert_86_p,  rate_has_74_a_h)

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
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), , by = 0.1))),
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
  do(data.frame(temp = seq(min_temp_2, max_temp_2, by = 0.1))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 15))


# create an empty row to use in combination of bootstrap
boot_preds_mu_i <- data.frame(matrix(0, ncol = 2, nrow = length(boot1_preds_beta_h$temp)))
colnames(boot_preds_mu_i) = c("temp","pred") 

# set the starting point of second bootstrap 
temp_row_length_2 <- length(seq(min_temp_2, max_temp_2, by = 0.1))

# set the end point of second bootstrap 
temp_row_length_1 <- length(seq(min_temp, max_temp, by = 0.1))

# use for loop to combine the two bootstraps 
for (k in 1:(length(boot1_preds_beta_h$temp)/temp_row_length_1)){
  
  temp_boot_1 <- boot3_preds_1[(1+(k-1)*temp_row_length_1):(k*temp_row_length_1),]
  temp_boot_2 <- boot3_preds_2[(1+(k-1)*temp_row_length_2):(k*temp_row_length_2),]
  
  boot_preds_mu_i[(1+(k-1)*temp_row_length_1):(k*temp_row_length_1), 1:2] <- rbind(temp_boot_1[which(temp_boot_1$temp <= 37), 5:6], temp_boot_2[which(temp_boot_2$temp > 37), 6:7])
  
}


################################ MORTALITIY OF MIRACIDIA ################

#R.M. Anderson et al. 1982, experiment with S. mansoni
temp_anderson_82 <- c(5,            10,        15,       20,        25,        30,        35,      40)
rate_anderson_82 <- c(4.957653,  2.141901,  1.50839,  1.974984,  2.513615,  4.322767,  4.490178,  5.114)


#S. K. Prah and C. James, 1977 S. mansoni 
temp_prah_77_man <- c(7.5,            20,      27.5,       36.5)
rate_prah_77_man <- c(2.56078,    0.6536434, 1.499468,   3.940166)


#R.E. Purnell, 1966, experiment with S. mansoni and two time range
temp_purnell_66 <- c(12,                       14,                           16,                 18.5,                        21.5,                     24.8,                  28.6,             32.7)
rate_purnell_66 <- c((3.261705+1.637893)/2,   (2.396054+1.785148)/2,     (2.904859+1.957561)/2,     (1.520372+3.054279)/2,        (2.337589+2.812790)/2,       (2.752958+3.891444)/2,  (3.842463+8.757026)/2, (4.877587+9.721674)/2)

#Combine the data set in one row 
temp <- c(temp_anderson_82, temp_prah_77_man, temp_purnell_66)
rate <- c(rate_anderson_82, rate_prah_77_man, rate_purnell_66)

# keep just a single curve
d <- data.frame(temp, rate)


# fit Sharpe-Schoolfield model
d_fit <- nest(d, data = c(temp, rate)) %>%
  mutate(thomas_2017 = map(data, ~nls_multstart(rate~thomas_2017(temp = temp, a,b,c,d,e),
                                                data = .x,
                                                iter = c(3,3,3,3,3),
                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') - 10,
                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') + 10,
                                                lower = get_lower_lims(.x$temp, .x$rate, model_name = 'thomas_2017'),
                                                upper = get_upper_lims(.x$temp, .x$rate, model_name = 'thomas_2017'),
                                                supp_errors = 'Y',
                                                convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), by = 0.1))),
         # predict over that data,
         preds =  map2(thomas_2017, new_data, ~augment(.x, newdata = .y)))


# refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(rate~thomas_2017(temp = temp, a,b,c,d,e),
                               data = d,
                               start = coef(d_fit$thomas_2017[[1]]),
                               lower = get_lower_lims(d$temp, d$rate, model_name = 'thomas_2017'),
                               upper = get_upper_lims(d$temp, d$rate, model_name = 'thomas_2017'),
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
  mutate(pred = thomas_2017(temp = temp, a,b,c,d,e))


# remove the bootstrap does not exist 
boot1_preds_mu_m <- head(boot1_preds, -n*length(seq(min_temp, max_temp, by = 0.1)))

############################ MORTALITY OF SNAIL ##############################


#(McCreesh et al. 2014), Biomphalaria sudanica
temp_mccreesh_14 <- c(    13.4,        15.7,       16.7,        18.9,        20.9,        22.8,        26.7,        28.3,        29.5,      32.0)
rate_mccreesh_14 <- c(0.009404046, 0.01682591, 0.005224019, 0.007002545, 0.001578888, 0.008951098, 0.007300832,  0.02229403, 0.03898983, 0.05906629)


#C. C. Appleton, 1977, experiment with Biomphalaria pfeifferi
temp_applet_77 <- c(  25,             27,       29)
rate_applet_77 <- c(0.002832593, 0.01014083, 0.02164846)


#El- Hassan 1974, Biomphalaria Alexandrina
temp_has_74_a <- c(   10,      12.5,    15,      20,      25,      30)
rate_has_74_a <- c(0.01161, 0.01119, 0.01077, 0.00596, 0.00674, 0.01505)


#R. F. Sturrock, 1966, experiment with Biomphalaria pfeifferi
temp_sturr_66 <- c(  19,            25,        30)
rate_sturr_66 <- c(0.01665,      0.01348,   0.03842)


#Foster 1964, Biomphalaria pfeifleri
temp_foster_64 <- c(  22.85,        24.01,      26.26,    28.07)
rate_foster_64 <- c(0.006674626, 0.01155245, 0.0260108, 0.03678792)


#(El-Emam and Madsen 1982),  Biomphalaria alexandrina
temp_eleman_82 <- c(   18,         26,          28)
rate_eleman_82 <- c(0.006188814, 0.006188814, 0.009373412)


#(Shiff and Garnett 1963), Biomphalaria pfeifferi. 
temp_shiff_63 <- c(      18,                 22,            25,                27)
rate_shiff_63 <- c( 0.001650303,     0.002659904,      0.002373008,     0.004170385)

# combine the data set in one row 
temp <- c(temp_mccreesh_14, temp_applet_77,  temp_has_74_a,  temp_sturr_66, temp_foster_64, temp_eleman_82, temp_shiff_63) 
rate <- c(rate_mccreesh_14, rate_applet_77,  rate_has_74_a,  rate_sturr_66, rate_foster_64, rate_eleman_82, rate_shiff_63)

# keep just a single curve
d <- data.frame(temp, rate)

# fit Sharpe-Schoolfield model
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
fit_nlsLM <- minpack.lm::nlsLM(rate~spain_1982(temp = temp, a,b,c,r0),
                               data = d,
                               start = coef(d_fit$spain_1982[[1]]),
                               lower = get_lower_lims(d$temp, d$rate, model_name = 'spain_1982'),
                               upper = get_upper_lims(d$temp, d$rate, model_name = 'spain_1982'),
                               weights = rep(1, times = nrow(d)))

# bootstrap using case resampling
boot1 <- Boot(fit_nlsLM, method = 'case')

# create predictions of each bootstrapped model
boot1_preds_mu <- boot1$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min_temp, max_temp, by = 0.1))) %>%
  ungroup() %>%
  mutate(pred = spain_1982(temp = temp, a,b,c,r0))

# we did another bootstrap for higher temperature 
min_temp_2 <- 34
max_temp_2 <- 37

#P. H. Joubert, 1986, Biomphalaria pfeifferi 
temp_joubert_86_p <- c(     34,         36,       38,       40)
rate_joubert_86_p <- c(0.1457006, 0.6661235, 1.649401,  5.821733)


#El- Hassan 1974, Biomphalaria Alexandrina
temp_has_74_a_h <- c(35,    37)
rate_has_74_a_h <- c(0.03891, 0.3)


#Combine the data set in one curve 
temp <- c(temp_joubert_86_p,  temp_has_74_a_h) 
rate <- c(rate_joubert_86_p,  rate_has_74_a_h)

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
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), , by = 0.1))),
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
  do(data.frame(temp = seq(min_temp_2, max_temp_2, by = 0.1))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 15))

# create an empty row to save the two bootstraps 
boot_preds_mu <- data.frame(matrix(0, ncol = 2, nrow = length(boot1_preds_beta_h$temp)))
colnames(boot_preds_mu) = c("temp","pred") 

# set the starting point of second bootstrap 
temp_row_length_1 <- length(seq(min_temp, max_temp, by = 0.1))

# set the starting end of first bootstrap 
temp_row_length_2 <- length(seq(min_temp_2, max_temp_2, by = 0.1))

for (k in 1:(length(boot1_preds_beta_h$temp)/temp_row_length_1)){
  
  temp_boot_1 <- boot1_preds_mu[(1+(k-1)*temp_row_length_1):(k*temp_row_length_1),]
  temp_boot_2 <- boot3_preds_2[(1+(k-1)*temp_row_length_2):(k*temp_row_length_2),]
  
  boot_preds_mu[(1+(k-1)*temp_row_length_1):(k*temp_row_length_1), 1:2] <- rbind(temp_boot_1[which(temp_boot_1$temp <= 37), 6:7], temp_boot_2[which(temp_boot_2$temp > 37), 6:7])
  
}


######################## RUN FOR MEAN PARASITE BURDEN ################
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
num_boot <- 2

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
  
  # Set a matrix to record outputs 
  out_puts <- matrix(0, nrow = 6, ncol = length(temperature))
  
  #Call the library 
  library(deSolve)
  
  # run the model for each temperature
  for (i in 1:length(temperature)){
    
    # This is mu_i function 
    preds_mu_i<- boot_preds_mu_i$pred[(1+(j-1)*length(temperature)):(length(temperature)*j)]
    
    # This is mu function 
    preds_mu <- boot_preds_mu$pred[(1+(j-1)*length(temperature)):(length(temperature)*j)]
    
    # Generate linearly interpolate point with temperature dependent parameter function 
    nu_s_afun <- approxfun(x = sample_parameters, y = rep(boot1_preds_nu_s$pred[(1+(j-1)*length(temperature)):(length(temperature)*j)][i], length(sample_parameters)))  
    mu_m_afun <- approxfun(x = sample_parameters, y =   rep(boot1_preds_mu_m$pred[(1+(j-1)*length(temperature)):(length(temperature)*j)][i], length(sample_parameters)))  
    mu_afun <- approxfun(x = sample_parameters, y =   rep(preds_mu[i], length(sample_parameters)))  
    sigma_s_afun <- approxfun(x = sample_parameters, y =   rep(boot1_preds_sigma_s$pred[(1+(j-1)*length(temperature)):(length(temperature)*j)][i], length(sample_parameters)))  
    mu_i_afun <- approxfun(x = sample_parameters, y =   rep(preds_mu_i[i], length(sample_parameters)))  
    nu_c_afun <- approxfun(x = sample_parameters, y =   rep(boot1_preds_nu_c$pred[(1+(j-1)*length(temperature)):(length(temperature)*j)][i], length(sample_parameters)))  
    mu_c_afun <- approxfun(x = sample_parameters, y =   rep(boot3_preds_mu_c$pred[(1+(j-1)*length(temperature)):(length(temperature)*j)][i], length(sample_parameters)))  
    delta_e_afun <- approxfun(x = sample_parameters, y =   rep(boot1_preds_delta_e$pred[(1+(j-1)*length(temperature)):(length(temperature)*j)][i], length(sample_parameters)))  
    beta_s_afun <- approxfun(x = sample_parameters, y =   rep(boot1_preds_beta_s$pred[(1+(j-1)*length(temperature)):(length(temperature)*j)][i] * contact_reduction_s, length(sample_parameters)))  
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
    
    y0 <- c(S = (60434 + 33232)/2, E = 1285, I = 2340, P = 2, P_m = 105)
    
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
write.csv(data.frame(boot_out_puts), "mansoni_mpb_boot_out_puts.csv", row.names=FALSE)



