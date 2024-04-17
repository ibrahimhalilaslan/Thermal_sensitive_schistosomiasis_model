# load packages
library(boot)
library(car)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)


# load in data
#Miracidia hatching rate  (Nguyen et al. 2021, Schistosoma mansoni)
temp <- c(5,  9,  13, 17,  21, 22, 25,  29,  33, 37)
rate <- c(70, 60, 90, 100, 76, 114, 137, 105, 150, 135)/300


# keep just a single curve
d <- data.frame(temp, rate)



# fit Sharpe-Schoolfield model
d_fit <- nest(d, data = c(temp, rate)) %>%
  mutate(flinn_1991 = map(data, ~nls_multstart(rate~flinn_1991(temp = temp, a, b, c),
                                                data = .x,
                                                iter = c(5,5,5),
                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') - 10,
                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') + 10,
                                                lower = get_lower_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                                                upper = get_upper_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                                                supp_errors = 'Y')),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))),
         # predict over that data,
         preds =  map2(flinn_1991, new_data, ~augment(.x, newdata = .y)))

# unnest predictions
d_preds <- select(d_fit, preds) %>%
  unnest(preds)

# refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(rate~flinn_1991(temp = temp, a, b, c), 
                               data = d,
                               start = coef(d_fit$flinn_1991[[1]]),
                               lower = get_lower_lims(d$temp, d$rate, model_name = 'flinn_1991'),
                               upper = get_upper_lims(d$temp, d$rate, model_name = 'flinn_1991'),
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
  mutate(pred = flinn_1991(temp = temp, a, b, c))

# calculate bootstrapped confidence intervals
boot1_conf_preds <- group_by(boot1_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()

#########

# plot bootstrapped CIs
pdf(file = "mir_hatch_rate.pdf", width = 5, height = 5)
ggplot() +
  geom_line(aes(temp, .fitted), d_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d, size = 2, alpha = 0.5) + 
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 17),
        axis.text.x = element_text(color="black", size=17), axis.text.y = element_text(color="black", size=17),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=14, face="bold.italic", hjust=0.5))+
  labs(x = 'Temperature (ºC)',
       y = 'Prob. of miracidia hatching success') 
dev.off() 



