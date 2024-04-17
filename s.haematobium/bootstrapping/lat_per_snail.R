# load packages
library(boot)
library(car)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)


#W. Pflüger 1984, Bulinus truncatus with S. haematobium 
temp_pfluger_84 <- c(18,       19,     20,    21,    22,     23,    25,  28,   30,    31,  32)
rate_pfluger_84 <- c(1/117,   1/110,  1/79, 1/60,   1/57,  1/56,  1/44, 1/32, 1/30, 1/31, 1/30)




#R. M. Gordon, 1934, Physopsis (bolinus) globosa  with haematobium 
temp_gordon_34_glo <- c(22,     26.3,      26.4,     26.8,       31.9,    32.1,     33,     35.2)
rate_gordon_34_glo <- c(1/67,  1/38.5,    1/39.3,   1/34.8,     1/27.5,   1/23,   1/25.5,   1/26.3)


#load data, Latent period of snail
temp <- c(temp_pfluger_84, temp_gordon_34_glo)
rate <- c(rate_pfluger_84, rate_gordon_34_glo)

# keep just a single curve
d <- data.frame(temp, rate)

# fit model
d_fit <- nest(d, data = c(temp, rate)) %>%
  mutate(lrf_1991 = map(data, ~nls_multstart(rate~lrf_1991(temp = temp, rmax, topt, tmin, tmax),
                                                  data = .x,
                                                  iter = c(3,3,3,3),
                                                  start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'lrf_1991') - 10,
                                                  start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'lrf_1991') + 10,
                                                  lower = get_lower_lims(.x$temp, .x$rate, model_name = 'lrf_1991'),
                                                  upper = get_upper_lims(.x$temp, .x$rate, model_name = 'lrf_1991'),
                                                  supp_errors = 'Y',
                                                  convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))),
         # predict over that data,
         preds =  map2(lrf_1991, new_data, ~augment(.x, newdata = .y)))

# unnest predictions
d_preds <- select(d_fit, preds) %>%
  unnest(preds)

# refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(rate~lrf_1991(temp = temp, rmax, topt, tmin, tmax),
                               data = d,
                               start = coef(d_fit$lrf_1991[[1]]),
                               lower = get_lower_lims(d$temp, d$rate, model_name = 'lrf_1991'),
                               upper = get_upper_lims(d$temp, d$rate, model_name = 'lrf_1991'),
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
  mutate(pred = lrf_1991(temp = temp, rmax, topt, tmin, tmax))

# calculate bootstrapped confidence intervals
boot1_conf_preds <- group_by(boot1_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()

# plot bootstrapped CIs
pdf(file = "lat_per_snail.pdf", width = 5, height = 5)
ggplot() +
  geom_line(aes(temp, 1/.fitted), d_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = 1/conf_lower, ymax = 1/conf_upper), boot1_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, 1/rate), d, size = 2, alpha = 0.5) + 
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 17),
        axis.text.x = element_text(color="black", size=17), axis.text.y = element_text(color="black", size=17),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=14, face="bold.italic", hjust=0.5))+
  labs(x = 'Temperature (ºC)',
       y = 'Incubation period in snails (days)') 
dev.off() 

