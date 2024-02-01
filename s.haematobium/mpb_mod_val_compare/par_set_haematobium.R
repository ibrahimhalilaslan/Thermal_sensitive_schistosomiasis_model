# Model parameters are set in this script 

# Clear the global environment 

rm(list = ls(all = TRUE)) 

library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)



#The factor of density-dependent fecundity rate
# We adjust this value to get the abundance of sanil in a site same as C. Wood at al. 2019 paper
nu <- 5.662715e-07 # (~ 1000000 snails)


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


#Rest of parameter in here are temperature dependent 
#and they have been fitted the best curve in the other scripts 

################ FECUNDITY RATE OF SNAILS #################

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
# Combine the data set in one row 
temp <-c(temp_ches_17, temp_has_74_t,  temp_shiff_67, temp_kubir_10,  32, 33, 34)
rate <-c(rate_ches_17, rate_has_74_t,  rate_shiff_67, rate_kubir_10,   0, 0, 0)


# keep just a single curve
.x <- data.frame(temp, rate)

fit_nu_s <- nls_multstart(rate~johnsonlewin_1946(temp = temp, r0, e, eh, topt),
                          data = .x,
                          iter = c(4,4,4,4),
                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'johnsonlewin_1946') - 1,
                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'johnsonlewin_1946') + 1,
                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'johnsonlewin_1946'),
                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'johnsonlewin_1946'),
                          supp_errors = 'Y',
                          convergence_count = FALSE)

fn_nu_s <- function(x) augment(fit_nu_s, newdata = x)

################ MORTALITY RATE OF SNAILS #################

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


#Combine the data set in one row 
temp <- c(temp_ches_17, temp_has_74_t,  temp_eleman_82, temp_shiff_64,  temp_kubirizajournal_10)
rate <- c(rate_ches_17, rate_has_74_t,  rate_eleman_82, rate_shiff_64,  rate_kubirizajournal_10)

# keep just a single curve
.x <- data.frame(temp, rate)

# fit model

fit_mu_1 <- nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                          data = .x,
                          iter = c(4,4,4),
                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 0.5,
                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 0.5,
                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                          supp_errors = 'Y',
                          convergence_count = FALSE)

fn_mu_1 <- function(x) augment(fit_mu_1, newdata = x)

################ MORTALITY RATE OF SNAILS FOR HIGHER TEMPERATURE #################

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

# Combine the data set in one row 
temp <- c(temp_ches_17_2,temp_joubert_86_g, temp_joubert_86_a, temp_has_74_t_2)
rate <- c(rate_ches_17_2, rate_joubert_86_g, rate_joubert_86_a, rate_has_74_t_2)

.x <- data.frame(temp, rate)

fit_mu_2 <- nls_multstart(rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 15),
                       data = .x,
                       iter = c(4,4,4,4),
                       start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981') - 10,
                       start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981') + 10,
                       lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981'),
                       upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981'),
                       supp_errors = 'Y',
                       convergence_count = FALSE)

fn_mu_2 <- function(x) augment(fit_mu_2, newdata = x)

################ MORTALITY RATE OF INFECTED SNAILS #################

# Chester Kalinda 2017, Bulinus globosu with S. haematobium
temp_chester_17 <- c(15.5,            21.2,        25.8,           31,         36)
rate_chester_17 <- c(0.00244925,  0.006901087,   0.01151302,   0.02073716,  0.04642702)

# W. Pflüger 1984, Bulinus truncatus with S. haematobium
temp_pfluger_84 <- c(18,              19,         20,          21,                 22,      23,             25,            28,        30,          31,    32)
rate_pfluger_84 <- c(0.0131776,  0.01141752,     0.01963277, 0.02130685,     0.02204905, 0.02636316,    0.02185846,   0.01590461, 0.01570362,  0.02473268, 0.0262486)

#(Chu et al. 1966) Bulinus trancatus with S. haematobium  
temp_chu_66 <- c(   9.967,          13,          16.067,         19,         22.2,         25.067,          28,        31.133,        34.2,         37.067,      39.967)
rate_chu_66 <- c(0.0058551278, 0.0085350331, 0.0061192732, 0.0195844770, 0.0085350331, 0.0032587267, 0.0004931228, 0.0061192732, 0.0058551278, 0.0968812480, 0.0915518574)

# Combine the data set in one row 
temp <- c(temp_chester_17, temp_pfluger_84, temp_chu_66)    
rate <- c(rate_chester_17, rate_pfluger_84, rate_chu_66)


# keep just a single curve
.x <- data.frame(temp, rate)
# fit model
fit_mu_i_1 <- nls_multstart(rate~spain_1982(temp = temp, a,b,c,r0),
                            data = .x,
                            iter = c(4,4,4,4),
                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') - 1,
                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') + 1,
                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                            supp_errors = 'Y')

fn_mu_i_1 <- function(x) augment(fit_mu_i_1, newdata = x)
# we use same mortality rate for higher temperature 
fn_mu_i_2 <- fn_mu_2

################ MORTALTIY RATE OF MIRACIDA #################


#S. K. Prah and C. James, 1977 S. haematobium
temp_prah_77_haem <- c(7.5,            20,      27.5,       36.5)
rate_prah_77_haem <- c(2.275331,    2.210028,  2.24239,  4.354702)


#Combine the data set
temp <- temp_prah_77_haem
rate <- rate_prah_77_haem


# keep just a single curve
.x <- data.frame(temp, rate)

# fit model
fit_mu_m <- nls_multstart(rate~spain_1982(temp = temp, a,b,c,r0),
                          data = .x,
                          iter = c(4,4,4,4),
                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') - 1,
                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') + 1,
                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                          supp_errors = 'Y',
                          convergence_count = FALSE)

fn_mu_m <- function(x) augment(fit_mu_m, newdata = x)

################ MORTALITY RATE OF CERCARIA #################

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
.x <- data.frame(temp, rate)

# fit model
fit_mu_c <- nls_multstart(rate~spain_1982(temp = temp, a,b,c,r0),
                          data = .x,
                          iter = c(4,4,4,4),
                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') - 1,
                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') + 1,
                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                          supp_errors = 'Y',
                          convergence_count = FALSE)

fn_mu_c <- function(x) augment(fit_mu_c, newdata = x)

################ HATCHING RATE OF MIRACIDIA #################

# Nguyen et al. 2021, experiment with Schistosoma mansoni)
temp <- c(5,  9,  13, 17,  21, 22, 25,  29,  33, 37)
rate <- c(70, 60, 90, 100, 76, 114, 137, 105, 150, 135)/300

# keep just a single curve
.x <- data.frame(temp, rate)

fit_delta_e <- nls_multstart(rate~flinn_1991(temp = temp, a, b, c),
                             data = .x,
                             iter = c(5,5,5),
                             start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') - 10,
                             start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') + 10,
                             lower = get_lower_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                             upper = get_upper_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                             supp_errors = 'Y')
fn_delta_e <- function(x) augment(fit_delta_e, newdata = x)

# PREPATENT PERIOD IN SNAILS 

#W. Pflüger 1984, Bulinus truncatus with S. haematobium 
temp_pfluger_84 <- c(18,       19,     20,    21,    22,     23,    25,  28,   30,    31,  32)
rate_pfluger_84 <- c(1/117,   1/110,  1/79, 1/60,   1/57,  1/56,  1/44, 1/32, 1/30, 1/31, 1/30)


#R. M. Gordon, 1934, Physopsis (bolinus) globosa  with haematobium 
temp_gordon_34_glo <- c(22,     26.3,      26.4,     26.8,       31.9,    32.1,     33,     35.2)
rate_gordon_34_glo <- c(1/67,  1/38.5,    1/39.3,   1/34.8,     1/27.5,   1/23,   1/25.5,   1/26.3)


#Combine the data set in one curve 
temp <- c(temp_pfluger_84, temp_gordon_34_glo)
rate <- c(rate_pfluger_84, rate_gordon_34_glo)

# keep just a single curve
.x <- data.frame(temp, rate)

# fit model
fit_sigma_s <- nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                             data = .x,
                             iter = c(4,4,4),
                             start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                             start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                             lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                             upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                             supp_errors = 'Y',
                             convergence_count = FALSE)

fn_sigma_s <- function(x) augment(fit_sigma_s, newdata = x)

################ CERCARI RELEASE RATE #################


#Nguyen et al. 2021, S. mansoni with Biomphalaria glabrata
temp_nguyen_21 <- c(5,  9,     13,     17,     21,     25,     29,     33,     37)
rate_nguyen_21 <- c(0, 15,     260,   1350,   1970,   2900,   1970,   1770,   690)

#(Upatham 1973), B. glabrata with S. mansoni  
temp_upatham_73 <- c(10, 13, 16,   19, 22,    25,   28,  31,   34,   37, 40)
rate_upatham_73 <- c( 0,  0, 206, 114, 442, 1200, 1400, 1233, 1241, 243, 121)


# Combine the data set in one curve 
temp <- c(temp_nguyen_21, temp_upatham_73)
rate <- c(rate_nguyen_21, rate_upatham_73)


# keep just a single curve
.x <- data.frame(temp, rate)

# fit model
fit_nu_c <- nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                          data = .x,
                          iter = c(4,4,4),
                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                          supp_errors = 'Y')

fn_nu_c <- function(x) augment(fit_nu_c, newdata = x)

################ DISEASES TRANSMISSION RATE IN SNAILS #################

#(Prah and James 1977), Bulinus (Physopsis) globosus with S. haematobium.
prah_77_temp <- c(       5,        12,            19,     26)
prah_77_rate <- c(0.001279126, 0.04565688,  0.1084776, 0.1849867)

# normalize the data set with the maximum value 
prah_77_rate <- prah_77_rate/max(prah_77_rate)

#(Shiff 1974), Bulinus globosus infected with S. haematobium.
shiff_74_temp <- c(      25,          20,           18,         16,          13,            11) 
shiff_74_rate <- c(0.03730699,  0.003512327,  0.01022082, 0.001782906,  0.002553961,  0.001281865)
# normalize the data set with the maximum value 
shiff_74_rate <- shiff_74_rate/max(shiff_74_rate)

#(Chu et al. 1966), Bulinus trancatus  infected with S. haematobium. 
chu_66_temp <- c(    10,            12,           14,         15,        20,          25,         30,        35,         38)
chu_66_rate <- c(0.002009499,  0.008145561,  0.02029471, 0.02016721, 0.04754069, 0.03577785,  0.06475144, 0.04395093, 0.05660857)
# normalize the data set with the maximum value 
chu_66_rate <- chu_66_rate/max(chu_66_rate)           

# combine the data set in one row
temp <- c(prah_77_temp, shiff_74_temp, chu_66_temp)
rate <- c(prah_77_rate, shiff_74_rate, chu_66_rate)

# keep just a single curve
.x <- data.frame(temp, rate)

# fit model
fit_beta_s <- nls_multstart(rate~flinn_1991(temp = temp, a, b, c),
                            data = .x,
                            iter = c(5,5,5),
                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') - 10,
                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') + 10,
                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                            supp_errors = 'Y')

# set the reduction rate for transmission rate in snails
contact_reduction_s <-0.32
fn_beta_s <- function(x) augment(fit_beta_s, newdata = x) * contact_reduction_s

############### DISEASES TRASMISSISON RATE IN SNAILS ############# 

#(Purnell 1966),, S. mansoni
purnell_66_temp <-      c(12,             15,           18,             21,        24,          27,           30,       33)
purnell_66_rate_1 <- c(0.9246078,      0.810981,    1.320736,       1.092269,   0.6998873,   1.231836,    0.7901292,  0.6292689)
# normalize the data set with the maximum value 
purnell_66_rate_1 <- purnell_66_rate_1/max(purnell_66_rate_1)


purnell_66_temp <- c(       12,             15,           18,             21,        24,          27,           30,       33)
purnell_66_rate_2 <- c(  1.0535,         1.002849,     0.9284408,       1.481275,   0.7225294,   0.6554545,    1.251547,  0.7475422)
# normalize the data set with the maximum value 
purnell_66_rate_2 <- purnell_66_rate_2/max(purnell_66_rate_2)


purnell_66_temp <- c(       12,             15,           18,             21,        24,          27,           30,       33)
purnell_66_rate_3 <- c(  2.182979,      1.523512,    1.777458,       1.840471,   0.8305635,   0.7054314,      1.339291,  0)
# normalize the data set with the maximum value 
purnell_66_rate_3 <- purnell_66_rate_3/max(purnell_66_rate_3)



#(Stirewalt 1954), S. mansoni 
stirewalt_54_temp <- c(    24,         27)
stirewalt_54_rate <- c(0.09590499,  0.1552625)
# normalize the data set with the maximum value 
stirewalt_54_rate <- stirewalt_54_rate/max(stirewalt_54_rate)


#(DeWitt 1965) with S. mansoni
dewitt_55_temp_1 <- c(     15,      20,          25,      30,      35)
dewitt_55_rate_1 <- c(0.4618475, 0.8464159,  1.00836, 1.237522, 0.6192221)
# normalize the data set with the maximum value 
dewitt_55_rate_1 <- dewitt_55_rate_1/max(dewitt_55_rate_1)

dewitt_55_temp_2 <- c(     15,      20,          25,      30,         35,        40)
dewitt_55_rate_2 <- c(0.385854, 0.6910937,    0.635236,  0.9783329, 1.0652, 0.2764123)
# normalize the data set with the maximum value 
dewitt_55_rate_2 <- dewitt_55_rate_2/max(dewitt_55_rate_2)

dewitt_55_temp_3 <- c(   10,        15,      20,         25,      30,       35,        40,   45)
dewitt_55_rate_3 <- c(0.7439515,  1.76077,  2.484149,  3.50177,  3.9655, 3.421058,  0.7532138, 0)
# normalize the data set with the maximum value 
dewitt_55_rate_3 <- dewitt_55_rate_3/max(dewitt_55_rate_3)

dewitt_55_temp_4 <- c(0, 5,        10,        15,      20,        25)
dewitt_55_rate_4 <- c(0, 0,    0.1190698, 0.6479505, 1.39834, 1.50036)
# normalize the data set with the maximum value 
dewitt_55_rate_4 <- dewitt_55_rate_4/max(dewitt_55_rate_4) 

# Combine the data set one row 
temp <- c(purnell_66_temp, purnell_66_temp, purnell_66_temp,stirewalt_54_temp, dewitt_55_temp_1, dewitt_55_temp_2, dewitt_55_temp_3, dewitt_55_temp_4)
rate <- c(purnell_66_rate_1, purnell_66_rate_2, purnell_66_rate_3, stirewalt_54_rate, dewitt_55_rate_1, dewitt_55_rate_2, dewitt_55_rate_3, dewitt_55_rate_4)


# keep just a single curve
.x <- data.frame(temp, rate)

# fit model

fit_beta_h <- nls_multstart(rate~briere2_1999(temp = temp, tmin, tmax, a,b),
                            data = .x,
                            iter = c(4,4,4,4),
                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') - 10,
                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') + 10,
                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                            supp_errors = 'Y')

# set reduction rate for the transmission rate in humans
contact_reduction_h <- 1.542194e-08
fn_beta_h <- function(x) augment(fit_beta_h, newdata = x) * contact_reduction_h

