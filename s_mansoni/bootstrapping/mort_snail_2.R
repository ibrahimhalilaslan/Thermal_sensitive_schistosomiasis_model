# load packages
library(boot)
library(car)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)

min_temp <- 30
max_temp <- 40

#P. H. Joubert, 1986, Biomphalaria pfeifferi 
temp_joubert_86_p <- c(     34,         36,       38,       40)
rate_joubert_86_p <- c(0.1457006, 0.6661235, 1.649401,  5.821733)


#El- Hassan 1974, Biomphalaria Alexandrina
temp_has_74_a_h <- c(35,    37)
rate_has_74_a_h <- c(0.03891, 0.3)


# Data 
temp <- c(temp_joubert_86_p,  temp_has_74_a_h) 
rate <- c(rate_joubert_86_p,  rate_has_74_a_h)

d <- data.frame(temp, rate)

# fit Sharpe-Schoolfield model to raw data
d_fit <- nest(d, data = c(temp, rate)) %>%
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
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))),
         # predict over that data,
         preds =  map2(sharpeschoollow_1981, new_data, ~augment(.x, newdata = .y)))

# refit model using nlsLM
fit_nlsLM2 <- nlsLM(rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 15),
                    data = d,
                    start = coef(d_fit$sharpeschoollow_1981[[1]]),
                    lower = get_lower_lims(d$temp, d$rate, model_name = 'sharpeschoollow_1981'),
                    upper = get_upper_lims(d$temp, d$rate, model_name = 'sharpeschoollow_1981'),
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
  do(data.frame(temp = seq(min(d$temp), max(d$temp),  length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 15))


# calculate bootstrapped confidence intervals
boot3_conf_preds <- group_by(boot3_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()

# plot bootstrapped CIs
pdf(file = "mort_snail_rate_2.pdf", width = 5, height = 5)
ggplot() +
  geom_line(aes(temp, .fitted), d_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot3_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d, size = 2, alpha = 0.5) + 
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 17),
        axis.text.x = element_text(color="black", size=17), axis.text.y = element_text(color="black", size=17),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=14, face="bold.italic", hjust=0.5))+
   labs(x = 'Temperature (ºC)',
       y = 'Mortality of snails') 
dev.off() 



###################
