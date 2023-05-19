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

#(Prah and James 1977), Bulinus (Physopsis) globosus with S. haematobium.
prah_77_temp <- c(       5,        12,            19,     26)
prah_77_rate <- c(0.001279126, 0.04565688,  0.1084776, 0.1849867)
prah_77_rate <- prah_77_rate/max(prah_77_rate)

#(Shiff 1974), Bulinus globosus infected with S. haematobium.
shiff_74_temp <- c(      25,          20,           18,         16,          13,            11) 
shiff_74_rate <- c(0.03730699,  0.003512327,  0.01022082, 0.001782906,  0.002553961,  0.001281865)
shiff_74_rate <- shiff_74_rate/max(shiff_74_rate)

#(Chu et al. 1966), Bulinus trancatus  infected with S. haematobium. 
chu_66_temp <- c(    10,            12,           14,         15,        20,          25,         30,        35,         38)
chu_66_rate <- c(0.002009499,  0.008145561,  0.02029471, 0.02016721, 0.04754069, 0.03577785,  0.06475144, 0.04395093, 0.05660857)
chu_66_rate <- chu_66_rate/max(chu_66_rate)           
                 


#The disease transmission rate in snails (Anderson et al. 1982, experiment with Biomphalaria glabra and schistosoma mansoni)
temp <- c(prah_77_temp, shiff_74_temp, chu_66_temp)
rate <- c(prah_77_rate, shiff_74_rate, chu_66_rate)

 
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
         delong = map(data, ~nls_multstart(rate~delong_2017(temp = temp, c, eb, ef, tm, ehc),
                                           data = .x,
                                           iter = c(4,4,4,4,4),
                                           start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'delong_2017') - 10,
                                           start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'delong_2017') + 10,
                                           lower = get_lower_lims(.x$temp, .x$rate, model_name = 'delong_2017'),
                                           upper = get_upper_lims(.x$temp, .x$rate, model_name = 'delong_2017'),
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
                                            data = .x,
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
