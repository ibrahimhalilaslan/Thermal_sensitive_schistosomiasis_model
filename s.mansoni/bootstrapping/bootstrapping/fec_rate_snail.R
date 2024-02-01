# load packages
library(boot)
library(car)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)


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



#R. F. Sturrock, 1972, experiment with Biomphalaria glabrata
temp_sturr_72 <- c(20,        25,        30,    35)
rate_sturr_72 <- c(0.05856,	0.08750,	0.06577, 0)


# Michelson 1961, Biomphalaria glabrata.
temp_michel_61 <- c(5, 15,   20,     25,     30,  35)
rate_michel_61 <- c(0, 0, 0.0132, 0.0349, 0.00969, 0)


# Shiff and Garnett 1963, Biomphalaria pfeifferi
temp_shiff_63 <- c(   18,     22,   25,     27)
rate_shiff_63 <- c(0.0138, 0.0507, 0.0349, 0.036)


temp <- c(temp_mccreesh_14, temp_applet_77, temp_has_74_a, temp_sturr_66, temp_sturr_72, temp_michel_61, temp_shiff_63)
rate <- c(rate_mccreesh_14, rate_applet_77, rate_has_74_a, rate_sturr_66, rate_sturr_72, rate_michel_61, rate_shiff_63)

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
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))),
         # predict over that data,
         preds =  map2(gaussian_1987, new_data, ~augment(.x, newdata = .y)))

# unnest predictions
d_preds <- select(d_fit, preds) %>%
  unnest(preds)

# refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(rate~gaussian_1987(temp = temp, rmax, topt, a),
                               data = d,
                               start = coef(d_fit$gaussian_1987[[1]]),
                               lower = get_lower_lims(d$temp, d$rate, model_name = 'gaussian_1987'),
                               upper = get_upper_lims(d$temp, d$rate, model_name = 'gaussian_1987'),
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
  do(data.frame(temp = seq(min(d$temp), max(d$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = gaussian_1987(temp = temp, rmax, topt, a))

#calculate the 95 percent credible interval for the optimal temperature 
opt_temp <- boot1_preds %>% group_by(iter) %>% summarise(opt_temp = temp[which.max(pred)])
quantile(opt_temp$opt_temp, c(0.025, 0.975),na.rm = TRUE)

# calculate bootstrapped confidence intervals
boot1_conf_preds <- group_by(boot1_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()

# plot bootstrapped CIs
pdf(file = "fecundity_sanil.pdf", width = 5, height = 5)
ggplot() +
  geom_line(aes(temp, .fitted), d_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d, size = 2, alpha = 0.5) + 
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 17),
        axis.text.x = element_text(color="black", size=17), axis.text.y = element_text(color="black", size=17),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=14, face="bold.italic", hjust=0.5))+
  labs(x = 'Temperature (ºC)',
       y = 'Fecundity rate of snails per day') 
dev.off() 

###################
