# load packages
library(boot)
library(car)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)

min_temp <- 12
max_temp <- 37


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
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))),
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
  do(data.frame(temp = seq(min(d$temp), max(d$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = spain_1982(temp = temp, a,b,c,r0))

# calculate bootstrapped confidence intervals
boot3_conf_preds <- group_by(boot3_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()

#######
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


# unnest predictions of original model fit
d_preds_2 <- select(d_fit, preds) %>%
  unnest(preds)

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


# calculate bootstrapped confidence intervals
boot3_conf_preds_2 <- group_by(boot3_preds_2, temp) %>%
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

