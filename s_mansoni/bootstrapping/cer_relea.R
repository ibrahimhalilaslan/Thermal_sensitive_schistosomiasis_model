# load packages
library(boot)
library(car)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)


#Nguyen et al. 2021, S. mansoni with Biomphalaria glabrata
temp_nguyen_21 <- c(5,  9,     13,     17,     21,     25,     29,     33,     37)
rate_nguyen_21 <- c(0, 15,     260,   1350,   1970,   2900,   1970,   1770,   690)

#(Upatham 1973), B. glabrata with S. mansoni  
temp_upatham_73 <- c(10, 13, 16,   19, 22,    25,   28,  31,   34,   37, 40)
rate_upatham_73 <- c( 0,  0, 206, 114, 442, 1200, 1400, 1233, 1241, 243, 121)



temp <- c(temp_nguyen_21, temp_upatham_73)
rate <- c(rate_nguyen_21, rate_upatham_73)



# keep just a single curve
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

# calculate bootstrapped confidence intervals
boot1_conf_preds <- group_by(boot1_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()

########
# plot bootstrapped CIs
pdf(file = "cercariae_release_rate.pdf", width = 5, height = 5)
ggplot() +
  geom_line(aes(temp, .fitted), d_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d, size = 2, alpha = 0.5) + 
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 17),
        axis.text.x = element_text(color="black", size=17), axis.text.y = element_text(color="black", size=17),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=14, face="bold.italic", hjust=0.5))+ 
  labs(x = 'Temperature (ºC)',
       y = 'Cercarial/sail/day') 
dev.off() 


############

