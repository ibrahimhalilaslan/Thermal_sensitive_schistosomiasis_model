# load packages
library(boot)
library(car)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)
library(dplyr)

#(Prah and James 1977), Bulinus (Physopsis) globosus with S. haematobium.
prah_77_temp <- c(       5,        12,            19,     26)
prah_77_rate <- c(0.001279126, 0.04565688,  0.1084776, 0.1849867)
prah_77_rate <- prah_77_rate/max(prah_77_rate)

#(Shiff 1974), Bulinus globosus infected with S. haematobium.
shiff_74_temp <- c(      25,          20,           18,         16,          13,            11) 
shiff_74_rate <- c(0.03730699,  0.003512327,  0.01022082, 0.001782906,  0.002553961,  0.001281865)
shiff_74_rate <- shiff_74_rate/max(shiff_74_rate)

#Chu et al. 1966), Bulinus trancatus  infected with S. haematobium. 
chu_66_temp <- c(    10,            12,           14,         15,        20,          25,         30,        35,         38)
chu_66_rate <- c(0.002009499,  0.008145561,  0.02029471, 0.02016721, 0.04754069, 0.03577785,  0.06475144, 0.04395093, 0.05660857)
chu_66_rate <- chu_66_rate/max(chu_66_rate)           


#The disease transmission rate in snails (Anderson et al. 1982, experiment with Biomphalaria glabra and schistosoma mansoni)
temp <- c(prah_77_temp, shiff_74_temp, chu_66_temp)
rate <- c(prah_77_rate, shiff_74_rate, chu_66_rate)




# keep just a single curve
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
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))),
         # predict over that data,
         preds =  map2(flinn_1991, new_data, ~augment(.x, newdata = .y)))

# unnest predictions of original model fit
d_preds <- select(d_fit, preds) %>%
  unnest(preds)

# refit model using nlsLM
fit_nlsLM2 <- nlsLM(rate~flinn_1991(temp = temp, a, b, c),
                    data = d,
                    start = coef(d_fit$flinn_1991[[1]]),
                    lower = get_lower_lims(d$temp, d$rate, model_name = 'flinn_1991'),
                    upper = get_upper_lims(d$temp, d$rate, model_name = 'flinn_1991'),
                    control = nls.lm.control(maxiter=500),
                    weights = rep(1, times = nrow(d)))



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
  mutate(pred = flinn_1991(temp = temp, a, b, c))

#calculate the 95 percent credible interval for the optimal temperature 
opt_temp <- boot3_preds %>% group_by(iter) %>% summarise(opt_temp = temp[which.max(pred)])
quantile(opt_temp$opt_temp, c(0.025, 0.975),na.rm = TRUE)

# calculate bootstrapped confidence intervals
boot3_conf_preds <- group_by(boot3_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()

#####
# plot bootstrapped CIs
pdf(file = "beta_s_trans_rate.pdf", width = 5, height = 5)
ggplot() +
  geom_line(aes(temp, .fitted), d_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot3_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d, size = 2) +
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 17),
        axis.text.x = element_text(color="black", size=17), axis.text.y = element_text(color="black", size=17),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=14, face="bold.italic", hjust=0.5))+
  labs(x = 'Temperature (ºC)', 
       y = 'Transmission rate in snails')
dev.off()      



########################

