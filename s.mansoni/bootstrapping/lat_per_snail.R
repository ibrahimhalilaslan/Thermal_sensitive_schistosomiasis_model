# load packages
library(boot)
library(car)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)


#W. Pflüger 1981, Biomphalaria glabrata with S. mansoni, diurnally changing water temperature 
temp_pfluger_81 <- c(24.3,     21.9,   20.2,  19.8,  17.5,   17.9,   19.2,  21.1, 18.9,    25,   29.9,  33)
rate_pfluger_81 <- c(1/26,    1/35.5,  1/44, 1/50.5, 1/68,  1/57.5,  1/52,  1/38, 1/57.5, 1/17, 1/16.5, 1/16.5)

#W. Pflüger 1980, Biomphalaria glabrata with S. mansoni, constant temperature 
temp_pfluger_80 <- c(  17,      19,     22,   25,      28,     30,      31,       32,     33,     34,    35)
rate_pfluger_80 <- c(1/92.5,  1/55.5,  1/34, 1/24.8, 1/19.2,  1/16.3,  1/15.9,  1/14.9, 1/15.3, 1/15.9, 1/16)


#Foster 1964, Biomphalaria pfeifleri with s. mansoni 
temp_foster_64 <- c(18,     21,    22.85,   24.01,   26.26,  28.07,  30.04, 31.75)
rate_foster_64 <- c(1/57,  1/37,    1/32,   1/30,     1/23,   1/19,   1/18,   1/16)


# R. M. Gordon, 1934, Planorbis pfeifferi with S. mansoni
temp_gordon_34_pfei <- c(21.7,    22,       26.3,      26.6,       26.9,     27.3,    27.7,    31.8,    32.1,    32.7,   32.8,   35)
rate_gordon_34_pfei <- c(1/36,   1/33,     1/23.5,     1/23,       1/22.2,   1/22,    1/19,   1/15.4, 1/15.8,  1/17.2,  1/14.7, 1/16.4)


#load data, Latent period of snail
temp <- c(temp_pfluger_81,  temp_pfluger_80, temp_foster_64, temp_gordon_34_pfei)
rate <- c(rate_pfluger_81,  rate_pfluger_80, rate_foster_64, rate_gordon_34_pfei)


d <- data.frame(temp, rate)


# fit Sharpe-Schoolfield model
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
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), by = 0.1))),
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
  geom_ribbon(aes(temp, ymax = 1/conf_lower, ymin = 1/conf_upper), boot1_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, 1/rate), d, size = 2, alpha = 0.5) + 
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 17),
        axis.text.x = element_text(color="black", size=17), axis.text.y = element_text(color="black", size=17),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=14, face="bold.italic", hjust=0.5))+
  labs(x = 'Temperature (ºC)',
       y = 'Incubation period in snails (days)') 
dev.off() 




###################
