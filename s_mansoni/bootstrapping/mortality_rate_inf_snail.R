# load packages
library(boot)
library(car)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)

min_temp <- 10
max_temp <- 40


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

#Data
temp <- c(temp_pfluger_81, temp_pfluger_81_sin, temp_pfluger_80,  temp_foster_64, temp_upatham_73)
rate <- c(rate_pfluger_81, rate_pfluger_81_sin, rate_pfluger_80,  rate_foster_64, rate_upatham_73)


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

# unnest predictions of original model fit
d_preds <- select(d_fit, preds) %>%
  unnest(preds)

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

# calculate bootstrapped confidence intervals
boot3_conf_preds <- group_by(boot3_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs
pdf(file = "mortality_rate_inf_snail.pdf", width = 5, height = 5)
ggplot() +
  geom_line(aes(temp, .fitted), d_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot3_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d, size = 2, alpha = 0.5) + 
  theme(panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 17),
        axis.text.x = element_text(color="black", size=17),axis.text.y = element_text(color="black", size=17)) +   labs(x = 'Temperature (ºC)',
       y = 'Mortality rate of infected snails per day') 
dev.off() 



min_temp_2 <- 30
max_temp_2 <- 40

#P. H. Joubert, 1986, Biomphalaria pfeifferi 
temp_joubert_86_p <- c(     34,         36,       38,       40)
rate_joubert_86_p <- c(0.1457006, 0.6661235, 1.649401,  5.821733)


#El- Hassan 1974, Biomphalaria Alexandrina
temp_has_74_a_h <- c(35,    37)
rate_has_74_a_h <- c(0.03891, 0.3)


# Data 
temp <- c(temp_joubert_86_p,  temp_has_74_a_h) 
rate <- c(rate_joubert_86_p,  rate_has_74_a_h)

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


# unnest predictions of original model fit
d_preds_2 <- select(d_fit, preds) %>%
  unnest(preds)

# bootstrap using residual resampling
boot3 <- Boot(fit_nlsLM2, method = 'residual')

# predict over new data
boot3_preds <- boot3$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min_temp_2, max_temp_2, by = 0.1))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 15))


# calculate bootstrapped confidence intervals
boot3_conf_preds_2 <- group_by(boot3_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


boot_conf_preds_combine <- rbind(boot3_conf_preds[which(boot3_conf_preds$temp <= 34), , ], boot3_conf_preds_2[which(boot3_conf_preds_2$temp > 34), , ])
d_preds_combine <- rbind(d_preds[which(d_preds$temp <= 34), ], d_preds_2[which(d_preds_2$temp > 34), ])


# plot bootstrapped CIs
pdf(file = "mortality_rate_inf_snail_combine.pdf", width = 5, height = 5)
ggplot() +
  geom_line(aes(temp, .fitted), d_preds_combine, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot_conf_preds_combine, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d, size = 2, alpha = 0.5) + 
  geom_point(aes(temp, rate), d_2, size = 2, alpha = 0.5) + 
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 17),
        axis.text.x = element_text(color="black", size=17), axis.text.y = element_text(color="black", size=17),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=14, face="bold.italic", hjust=0.5))+
  labs(x = 'Temperature (ºC)',
       y = 'Mortality rate of snails per day') 
dev.off() 







