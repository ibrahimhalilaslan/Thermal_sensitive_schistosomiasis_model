# load packages
library(boot)
library(car)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)


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
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))),
         # predict over that data,
         preds =  map2(johnsonlewin_1946, new_data, ~augment(.x, newdata = .y)))

# unnest predictions
d_preds <- select(d_fit, preds) %>%
  unnest(preds)

# refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(rate~johnsonlewin_1946(temp = temp, r0, e, eh, topt),
                               data = d,
                               start = coef(d_fit$johnsonlewin_1946[[1]]),
                               lower = get_lower_lims(d$temp, d$rate, model_name = 'johnsonlewin_1946'),
                               upper = get_upper_lims(d$temp, d$rate, model_name = 'johnsonlewin_1946'),
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
  mutate(pred = johnsonlewin_1946(temp = temp, r0, e, eh, topt))

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
       y = 'Fecundity rate') 
dev.off() 

