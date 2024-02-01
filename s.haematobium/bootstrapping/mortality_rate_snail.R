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

d <- data.frame(temp, rate)

# fit Sharpe-Schoolfield model
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

# unnest predictions
d_preds <- select(d_fit, preds) %>%
  unnest(preds)

# refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(rate~quadratic_2008(temp = temp, a, b, c),
                               data = d,
                               start = coef(d_fit$quadratic_2008[[1]]),
                               lower = get_lower_lims(d$temp, d$rate, model_name = 'quadratic_2008'),
                               upper = get_upper_lims(d$temp, d$rate, model_name = 'quadratic_2008'),
                               weights = rep(1, times = nrow(d)))

# bootstrap using case resampling
boot1 <- Boot(fit_nlsLM, method = 'case')

# look at the data
head(boot1$t)


# create predictions of each bootstrapped model
boot1_preds <- boot1$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min_temp, max_temp,  by = 0.1))) %>%
  ungroup() %>%
  mutate(pred = quadratic_2008(temp = temp, a, b, c))

#calculate the 95 percent credible interval for the optimal temperature 
opt_temp <- boot1_preds %>% group_by(iter) %>% summarise(opt_temp = temp[which.min(pred)])
quantile(opt_temp$opt_temp, c(0.025, 0.975),na.rm = TRUE)

# calculate bootstrapped confidence intervals
boot1_conf_preds <- group_by(boot1_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()

# plot bootstrapped CIs
pdf(file = "mortality_rate_snail.pdf", width = 5, height = 5)
ggplot() +
  geom_line(aes(temp, .fitted), d_preds[which(d_preds$temp <= max(d$temp)),], col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds[which(boot1_conf_preds$temp <= max(d$temp)), ,], fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d, size = 2, alpha = 0.5) + 
  theme(panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 17),
        axis.text.x = element_text(color="black", size=17),axis.text.y = element_text(color="black", size=17)) + 
  labs(x = 'Temperature (ºC)',
       y = 'Mortality rate of snails per day') 
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

boot_conf_preds_combine <- rbind(boot1_conf_preds[which(boot1_conf_preds$temp <= 34), , ], boot3_conf_preds_2[which(boot3_conf_preds_2$temp > 34), , ])
d_preds_combine <- rbind(d_preds[which(d_preds$temp <= 34), ], d_preds_2[which(d_preds_2$temp > 34), ])


# plot bootstrapped CIs
pdf(file = "mortality_rate_snail_combine.pdf", width = 5, height = 5)
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

