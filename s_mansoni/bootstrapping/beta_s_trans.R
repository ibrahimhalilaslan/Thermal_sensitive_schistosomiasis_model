# load packages
library(boot)
library(car)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)


# (Stirewalt 1954) Australorbis glabratus with S. mansoni.
stirewalt_54_temp <- c(    24,         27,        32)
stirewalt_54_rate <- c(0.07515584, 0.4808139, 0.3922359)
stirewalt_54_rate <- stirewalt_54_rate/max(stirewalt_54_rate)

#(Foster 1964), Biomphalaria pfeifferi with S. mansoni.
foster_64_temp <- c(    22.85,      24.01,       26.26,      28.07,       30.04,    31.75)
foster_64_rate <- c(0.04079299,   0.04333207,  0.2203394, 0.08210066,  0.05756868, 0.07958718)
foster_64_rate <- foster_64_rate/max(foster_64_rate)

#(Prah and James 1977),  Biomphalaria pfeifferi with S. mansoni. 
prah_77_temp <- c(    5,             15,         20,         22.5)
prah_77_rate <- c(0.0006344901, 0.05059652, 0.1629059,  0.2249733)
prah_77_rate <- prah_77_rate/max(prah_77_rate)

#(Anderson et al. 1982), Biomphalaria glabrata with S. mansoni 
anderson_82_temp <- c(15,               20,          25,           30,         35)
anderson_82_rate <- c(0.3658224,    0.3865833,    0.9074138,     0.487644, 0.4107836)
anderson_82_rate <-  anderson_82_rate/max(anderson_82_rate)

#(DeWitt 1955), B. glabrata with S. mansoni. 
dewitt_55_temp <-   c(10,      25,       35,    40)
dewitt_55_rate_1 <- c(0,   0.0078523, 0.0283053, 0)
dewitt_55_rate_1 <- dewitt_55_rate_1/max(dewitt_55_rate_1)
dewitt_55_rate_2 <- c(0,   0.032993,  0.0241807, 0)
dewitt_55_rate_2 <- dewitt_55_rate_2/max(dewitt_55_rate_2)
dewitt_55_rate_3 <- c(0,   0.01388562,  0.03692433, 0.01246124)
dewitt_55_rate_3 <- dewitt_55_rate_3/max(dewitt_55_rate_3)
dewitt_55_rate_4 <- c(0,   0.07405795,  0.09039597, 0.00401955)
dewitt_55_rate_4 <- dewitt_55_rate_4/max(dewitt_55_rate_4)

#Coelho and Bezerra 2006), Biomphalari glabrata.
coelho_temp_06 <- c(15, 20, 30)
coelho_rate_06 <- c(0.0001342012, 0.007285127, 0.04153248)
coelho_rate_06 <- coelho_rate_06/max(coelho_rate_06)

#(Upatham 1973), B. glabrata with S. mansoni. 
upatham_temp_73 <- c(10,      13,     16,           19,           22,           25,        28,          31,          34,          37,          40)
upatham_rate_73 <- c(0,        0,  0.002895729, 0.004360468,  0.009807099, 0.01299096, 0.01373957,  0.01528639,  0.02722317, 0.002743025,  0.00199986)
upatham_rate_73 <- upatham_rate_73/max(upatham_rate_73)

#The disease transmission rate in snails (Anderson et al. 1982, experiment with Biomphalaria glabra and schistosoma mansoni)
temp <- c(stirewalt_54_temp, foster_64_temp, prah_77_temp,  anderson_82_temp, dewitt_55_temp,   dewitt_55_temp,   dewitt_55_temp,   dewitt_55_temp,  coelho_temp_06, upatham_temp_73)
rate <- c(stirewalt_54_rate, foster_64_rate, prah_77_rate,  anderson_82_rate, dewitt_55_rate_1, dewitt_55_rate_2, dewitt_55_rate_3, dewitt_55_rate_4,  coelho_rate_06, upatham_rate_73)

d <- data.frame(temp, rate)

# fit Sharpe-Schoolfield model
d_fit <- nest(d, data = c(temp, rate)) %>%
  mutate(spain_1982 = map(data, ~nls_multstart(rate~spain_1982(temp = temp, a,b,c,r0),
                                                  data = d,
                                                  iter = c(4,4,4,4),
                                                  start_lower = get_start_vals(d$temp, d$rate, model_name = 'spain_1982') - 1,
                                                  start_upper = get_start_vals(d$temp, d$rate, model_name = 'spain_1982') + 1,
                                                  lower = get_lower_lims(d$temp, d$rate, model_name = 'spain_1982'),
                                                  upper = get_upper_lims(d$temp, d$rate, model_name = 'spain_1982'),
                                                  supp_errors = 'Y',
                                                  convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))),
         # predict over that data,
         preds =  map2(spain_1982, new_data, ~augment(.x, newdata = .y)))

# unnest predictions
d_preds <- select(d_fit, preds) %>%
  unnest(preds)

# refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(rate~spain_1982(temp = temp, a,b,c,r0),
                               data = d,
                               start = coef(d_fit$spain_1982[[1]]),
                               lower = get_lower_lims(d$temp, d$rate, model_name = 'spain_1982'),
                               upper = get_upper_lims(d$temp, d$rate, model_name = 'spain_1982'),
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
  mutate(pred = spain_1982(temp = temp, a,b,c,r0))

# calculate bootstrapped confidence intervals
boot1_conf_preds <- group_by(boot1_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()

# plot bootstrapped CIs
pdf(file = "beta_s_trans_rate.pdf", width = 5, height = 5)
ggplot() +
  geom_line(aes(temp, .fitted), d_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d, size = 2, alpha = 0.5) + 
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 17),
        axis.text.x = element_text(color="black", size=17), axis.text.y = element_text(color="black", size=17),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=14, face="bold.italic", hjust=0.5))+
  labs(x = 'Temperature (ºC)',
       y = 'Transmission rate in snails') 
dev.off() 


###################
