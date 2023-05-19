# load packages
library(boot)
library(car)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)


#R.M. Anderson et al. 1982, experiment with S. mansoni
temp_anderson_82 <- c(5,            10,        15,       20,        25,        30,        35,      40)
rate_anderson_82 <- c(4.957653,  2.141901,  1.50839,  1.974984,  2.513615,  4.322767,  4.490178,  5.114)


#S. K. Prah and C. James, 1977 S. mansoni 
temp_prah_77_man <- c(7.5,            20,      27.5,       36.5)
rate_prah_77_man <- c(2.56078,    0.6536434, 1.499468,   3.940166)


#R.E. Purnell, 1966, experiment with S. mansoni and two time range
temp_purnell_66 <- c(12,                       14,                           16,                 18.5,                        21.5,                     24.8,                  28.6,             32.7)
rate_purnell_66 <- c((3.261705+1.637893)/2,   (2.396054+1.785148)/2,     (2.904859+1.957561)/2,     (1.520372+3.054279)/2,        (2.337589+2.812790)/2,       (2.752958+3.891444)/2,  (3.842463+8.757026)/2, (4.877587+9.721674)/2)

#load data 
temp <- c(temp_anderson_82, temp_prah_77_man, temp_purnell_66)
rate <- c(rate_anderson_82, rate_prah_77_man, rate_purnell_66)


d <- data.frame(temp, rate)


# fit Sharpe-Schoolfield model
d_fit <- nest(d, data = c(temp, rate)) %>%
  mutate(thomas_2017 = map(data, ~nls_multstart(rate~thomas_2017(temp = temp, a,b,c,d,e),
                                               data = .x,
                                               iter = c(3,3,3,3,3),
                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') - 10,
                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') + 10,
                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'thomas_2017'),
                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'thomas_2017'),
                                               supp_errors = 'Y',
                                               convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))),
         # predict over that data,
         preds =  map2(thomas_2017, new_data, ~augment(.x, newdata = .y)))

# unnest predictions
d_preds <- select(d_fit, preds) %>%
  unnest(preds)

# refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(rate~thomas_2017(temp = temp, a,b,c,d,e),
                               data = d,
                               start = coef(d_fit$thomas_2017[[1]]),
                               lower = get_lower_lims(d$temp, d$rate, model_name = 'thomas_2017'),
                               upper = get_upper_lims(d$temp, d$rate, model_name = 'thomas_2017'),
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
  mutate(pred = thomas_2017(temp = temp, a,b,c,d,e))

# calculate bootstrapped confidence intervals
boot1_conf_preds <- group_by(boot1_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()

#########


# plot bootstrapped CIs
pdf(file = "mortality_rate_mir.pdf", width = 5, height = 5)
ggplot() +
  geom_line(aes(temp, .fitted), d_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d, size = 2, alpha = 0.5) + 
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 17),
        axis.text.x = element_text(color="black", size=17), axis.text.y = element_text(color="black", size=17),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=14, face="bold.italic", hjust=0.5))+
   labs(x = 'Temperature (ºC)',
       y = 'Mortality rate of miracidia per day') 
dev.off() 


###################
