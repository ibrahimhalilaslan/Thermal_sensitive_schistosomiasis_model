# load packages
library(boot)
library(car)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)

#load data 
#S. K. Prah and C. James, 1977 S. haematobium

temp_prah_77_haem <- c(7.5,              20,     27.5,       36.5)
rate_prah_77_haem <- c(2.275331,      2.210028,  2.24239,  4.354702)


temp <- temp_prah_77_haem
rate <- rate_prah_77_haem


d <- data.frame(temp, rate)

# fit Sharpe-Schoolfield model to raw data
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
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))),
         # predict over that data,
         preds =  map2(spain_1982, new_data, ~augment(.x, newdata = .y)))

# refit model using nlsLM
fit_nlsLM2 <- nlsLM(rate~spain_1982(temp = temp, a,b,c,r0),
                    data = d,
                    start = coef(d_fit$spain_1982[[1]]),
                    lower = get_lower_lims(d$temp, d$rate, model_name = 'spain_1982'),
                    upper = get_upper_lims(d$temp, d$rate, model_name = 'spain_1982'),
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
  do(data.frame(temp = seq(min(d$temp), max(d$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = spain_1982(temp = temp, a,b,c,r0))



# plot bootstrapped CIs
pdf(file = "mortality_rate_mir.pdf", width = 5, height = 5)
ggplot() +
  geom_line(aes(temp, .fitted), d_preds, col = 'blue') +
  geom_point(aes(temp, rate), d, size = 2, alpha = 0.5) + 
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 17),
        axis.text.x = element_text(color="black", size=17), axis.text.y = element_text(color="black", size=17),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=14, face="bold.italic", hjust=0.5))+ 
  labs(x = 'Temperature (ºC)',
       y = 'Mortality rate of miracidia per day') 
dev.off() 

