# This script perform multiple different thermal sensitive models to find the best fit 
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)

# write function to label ggplot2 panels
label_facets_num <- function(string){
  len <- length(string)
  string = paste('(', 1:len, ') ', string, sep = '')
  return(string)
}


# (Stirewalt 1954) Australorbis glabratus with S. mansoni.
stirewalt_54_temp <- c(    24,         27,        32)
stirewalt_54_rate <- c(0.07515584, 0.4808139, 0.3922359)

stirewalt_54_rate <- stirewalt_54_rate/max(stirewalt_54_rate)

#(Foster 1964), Biomphalaria pfeifferi with S. mansoni.
foster_64_temp <- c(    22.85,      24.01,       26.26,      28.07,       30.04,    31.75)
foster_64_rate <- c(0.04079299,   0.04333207,  0.2203394, 0.08210066,  0.05756868, 0.07958718)
foster_64_rate <- foster_64_rate/max(foster_64_rate)

#(Prah and James 1977),  Biomphalaria pfeifferi with S. mansoni. 
prah_77_temp <- c(    5,             15,         20,         22.5)
prah_77_rate <- c(0.0006344901, 0.05059652, 0.1629059,  0.2249733)
prah_77_rate <- prah_77_rate/max(prah_77_rate)

#(Anderson et al. 1982), Biomphalaria glabra with S. mansoni 
anderson_82_temp <- c(15,               20,          25,           30,         35)
anderson_82_rate <- c(0.3658224,    0.3865833,    0.9074138,     0.487644, 0.4107836)
anderson_82_rate <-  anderson_82_rate/max(anderson_82_rate)

#(DeWitt 1955), B. glabrata with S. mansoni. 
dewitt_55_temp <-   c(10,      25,       35,    40)
dewitt_55_rate_1 <- c(0,   0.0078523, 0.0283053, 0)
dewitt_55_rate_1 <- dewitt_55_rate_1/max(dewitt_55_rate_1)

dewitt_55_rate_2 <- c(0,   0.032993,  0.0241807, 0)
dewitt_55_rate_2 <- dewitt_55_rate_2/max(dewitt_55_rate_2)

dewitt_55_rate_3 <- c(0,   0.01388562,  0.03692433, 0.01246124)
dewitt_55_rate_3 <- dewitt_55_rate_3/max(dewitt_55_rate_3)

dewitt_55_rate_4 <- c(0,   0.07405795,  0.09039597, 0.00401955)
dewitt_55_rate_4 <- dewitt_55_rate_4/max(dewitt_55_rate_4)

#(Coelho and Bezerra 2006), Biomphalari glabrata.

coelho_temp_06 <- c(15, 20, 30)
coelho_rate_06 <- c(0.0001342012, 0.007285127, 0.04153248)
coelho_rate_06 <- coelho_rate_06/max(coelho_rate_06)

#(Upatham 1973), B. glabrata with S. mansoni. 
upatham_temp_73 <- c(10,      13,     16,           19,           22,           25,        28,          31,          34,          37,          40)
upatham_rate_73 <- c(0,        0,  0.002895729, 0.004360468,  0.009807099, 0.01299096, 0.01373957,  0.01528639,  0.02722317, 0.002743025,  0.00199986)
upatham_rate_73 <- upatham_rate_73/max(upatham_rate_73)

#The disease transmission rate in snails (Anderson et al. 1982, experiment with Biomphalaria glabra and schistosoma mansoni)
temp <- c(stirewalt_54_temp, foster_64_temp, prah_77_temp,  anderson_82_temp, dewitt_55_temp,   dewitt_55_temp,   dewitt_55_temp,   dewitt_55_temp,  coelho_temp_06, upatham_temp_73)
rate <- c(stirewalt_54_rate, foster_64_rate, prah_77_rate,  anderson_82_rate, dewitt_55_rate_1, dewitt_55_rate_2, dewitt_55_rate_3, dewitt_55_rate_4,  coelho_rate_06, upatham_rate_73)

 
# keep just a single curve
d <- data.frame(temp, rate)
# show the data
ggplot(d, aes(temp, rate)) +
  geom_point() +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Transmission rate in snail',
       title = 'Transmission rate across temperatures') 


# fit every model formulation in rTPC
d_fits <- nest(d, data = c(temp, rate)) %>%
  mutate(boatman = map(data, ~nls_multstart(rate~boatman_2017(temp = temp, rmax, tmin, tmax, a,b),
                                            data = .x,
                                            iter = c(4,4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') - 10,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') + 10,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         briere2 = map(data, ~nls_multstart(rate~briere2_1999(temp = temp, tmin, tmax, a,b),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') - 10,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') + 10,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         
         flinn = map(data, ~nls_multstart(rate~flinn_1991(temp = temp, a, b, c),
                                          data = .x,
                                          iter = c(5,5,5),
                                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') - 10,
                                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') + 10,
                                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
         gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                             data = .x,
                                             iter = c(4,4,4),
                                             start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                                             start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                                             lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                             upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
         hinshelwood = map(data, ~nls_multstart(rate~hinshelwood_1947(temp = temp, a, e, b, eh),
                                                data = .x,
                                                iter = c(5,5,5,5),
                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'hinshelwood_1947') - 1,
                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'hinshelwood_1947') + 1,
                                                lower = get_lower_lims(.x$temp, .x$rate, model_name = 'hinshelwood_1947'),
                                                upper = get_upper_lims(.x$temp, .x$rate, model_name = 'hinshelwood_1947'),
                                                supp_errors = 'Y',
                                                convergence_count = FALSE)),
         joehnk = map(data, ~nls_multstart(rate~joehnk_2008(temp = temp, rmax, topt, a, b, c),
                                           data = .x,
                                           iter = c(4,4,4,4, 4),
                                           start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'joehnk_2008') - 10,
                                           start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'joehnk_2008') + 10,
                                           lower = get_lower_lims(.x$temp, .x$rate, model_name = 'joehnk_2008'),
                                           upper = get_upper_lims(.x$temp, .x$rate, model_name = 'joehnk_2008'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)),
         
         
         kamykowski = map(data, ~nls_multstart(rate~kamykowski_1985(temp = temp, tmin, tmax, a, b, c),
                                               data = .x,
                                               iter = c(4,4,4,4,4),
                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'kamykowski_1985') - 10,
                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'kamykowski_1985') + 10,
                                               supp_errors = 'Y',
                                               convergence_count = FALSE)),
         lactin2 = map(data, ~nls_multstart(rate~lactin2_1995(temp = temp, a, b, tmax, delta_t),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'lactin2_1995') - 10,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'lactin2_1995') + 10,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'lactin2_1995'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'lactin2_1995'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         lrf = map(data, ~nls_multstart(rate~lrf_1991(temp = temp, rmax, topt, tmin, tmax),
                                        data = d,
                                        iter = c(3,3,3,3),
                                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'lrf_1991') - 10,
                                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'lrf_1991') + 10,
                                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'lrf_1991'),
                                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'lrf_1991'),
                                        supp_errors = 'Y',
                                        convergence_count = FALSE)),
         modifiedgaussian = map(data, ~nls_multstart(rate~modifiedgaussian_2006(temp = temp, rmax, topt, a, b),
                                                     data = .x,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') - 10,
                                                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') + 10,
                                                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006'),
                                                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         oneill = map(data, ~nls_multstart(rate~oneill_1972(temp = temp, rmax, ctmax, topt, q10),
                                           data = .x,
                                           iter = c(4,4,4,4),
                                           start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') - 10,
                                           start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') + 10,
                                           
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)),
         pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref, e, eh, topt, tref = 15),
                                          data = .x,
                                          iter = c(4,4,4,4),
                                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'pawar_2018') - 10,
                                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'pawar_2018') + 10,
                                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'pawar_2018'),
                                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'pawar_2018'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
         quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                              data = .x,
                                              iter = c(4,4,4),
                                              start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 0.5,
                                              start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 0.5,
                                              lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                              upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         ratkowsky = map(data, ~nls_multstart(rate~ratkowsky_1983(temp = temp, tmin, tmax, a, b),
                                              data = .x,
                                              iter = c(4,4,4,4),
                                              start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'ratkowsky_1983') - 10,
                                              start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'ratkowsky_1983') + 10,
                                              lower = get_lower_lims(.x$temp, .x$rate, model_name = 'ratkowsky_1983'),
                                              upper = get_upper_lims(.x$temp, .x$rate, model_name = 'ratkowsky_1983'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         
         sharpeschoolfull = map(data, ~nls_multstart(rate~sharpeschoolfull_1981(temp = temp, r_tref,e,el,tl,eh,th, tref = 15),
                                                     data = .x,
                                                     iter = c(4,4,4,4,4,4),
                                                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981') - 10,
                                                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981') + 10,
                                                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981'),
                                                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         sharpeschoolhigh = map(data, ~nls_multstart(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                                     data = .x,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') - 10,
                                                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') + 10,
                                                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                                                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         sharpeschoollow = map(data, ~nls_multstart(rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 15),
                                                    data = .x,
                                                    iter = c(4,4,4,4),
                                                    start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981') - 10,
                                                    start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981') + 10,
                                                    lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981'),
                                                    upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981'),
                                                    supp_errors = 'Y',
                                                    convergence_count = FALSE)),
         spain = map(data, ~nls_multstart(rate~spain_1982(temp = temp, a,b,c,r0),
                                          data = .x,
                                          iter = c(4,4,4,4),
                                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') - 1,
                                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') + 1,
                                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
         thomas1 = map(data, ~nls_multstart(rate~thomas_2012(temp = temp, a,b,c,topt),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2012') - 1,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2012') + 2,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'thomas_2012'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'thomas_2012'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         thomas2 = map(data, ~nls_multstart(rate~thomas_2017(temp = temp, a,b,c,d,e),
                                            data = .x,
                                            iter = c(3,3,3,3,3),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') - 10,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') + 10,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'thomas_2017'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'thomas_2017'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         
         weibull = map(data, ~nls_multstart(rate~weibull_1995(temp = temp, a, topt, b, c),
                                            data = d,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'weibull_1995') - 10,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'weibull_1995') + 10,
                                            supp_errors = 'Y', 
                                            convergence_count = FALSE)))


glimpse(select(d_fits, 1:7))


d_fits$boatman[[1]]

# stack models
d_stack <- select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', boatman:weibull)

# get parameters using tidy
params <- d_stack %>%
  mutate(., est = map(fit, tidy)) %>%
  select(-fit) %>%
  unnest(est) 

# get predictions using augment
newdata <- tibble(temp = seq(min(d$temp), max(d$temp), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)

# plot
pdf(file = "multip_trans_rate_in_snail.pdf")
ggplot(d_preds, aes(temp, rate)) +
  geom_point(aes(temp, rate), d) +
  geom_line(aes(temp, .fitted), col = 'blue') +
  facet_wrap(~model_name, labeller = labeller(model_name = label_facets_num), scales = 'free', ncol = 5) +
  theme_bw(base_size = 8) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(temp = 'Temperature (ºC)',
       rate = 'Death rate of cercariae',
       title = 'Fits of every model available in rTPC') +
  geom_hline(aes(yintercept = 0), linetype = 2)
dev.off()

d_ic <- d_stack %>%
  mutate(., info = map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  select(-fit) %>%
  unnest(info) %>%
  select(model_name, sigma, AIC, AICc, BIC, df.residual)

print(n = length(d_ic$model_name), d_ic[order(d_ic$AIC),])
