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



# Data 

temp <- c(temp_mccreesh_14, temp_applet_77,  temp_has_74_a,  temp_sturr_66, temp_foster_64, temp_eleman_82, temp_shiff_63) 
rate <- c(rate_mccreesh_14, rate_applet_77,  rate_has_74_a,  rate_sturr_66, rate_foster_64, rate_eleman_82, rate_shiff_63)



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

# unnest predictions
d_preds <- select(d_fit, preds) %>%
  unnest(preds)

# refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(rate~spain_1982(temp = temp, a,b,c,r0),
                               data = d,
                               start = coef(d_fit$spain_1982[[1]]),
                               lower = get_lower_lims(d$temp, d$rate, model_name = 'spain_1982'),
                               upper = get_upper_lims(d$temp, d$rate, model_name = 'spain_1982'),
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
  mutate(pred = spain_1982(temp = temp, a,b,c,r0))

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

boot_conf_preds_combine <- rbind(boot1_conf_preds[which(boot1_conf_preds$temp <= 36), , ], boot3_conf_preds_2[which(boot3_conf_preds_2$temp > 36), , ])
d_preds_combine <- rbind(d_preds[which(d_preds$temp <= 36), ], d_preds_2[which(d_preds_2$temp > 36), ])


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




###################
