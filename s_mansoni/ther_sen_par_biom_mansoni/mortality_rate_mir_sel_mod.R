# This script perform the best curve for mu_m
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)


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



# keep just a single curve
d <- data.frame(temp, rate)
.x <- d
# fit model

fit <- nls_multstart(rate~thomas_2017(temp = temp, a,b,c,d,e),
                     data = .x,
                     iter = c(3,3,3,3,3),
                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') - 10,
                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') + 10,
                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'thomas_2017'),
                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'thomas_2017'),
                     supp_errors = 'Y')
fit


# calculate additional traits
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data

new_data <- data.frame(temp = seq(min(d$temp), max(d$temp), 0.5))
preds <- augment(fit, newdata = new_data)


pdf(file = "mortality_rate_mir.pdf", width = 5, height = 5)
par(mar=c(4, 4, 5, 1), xpd = TRUE)
plot(preds$temp, preds$.fitted, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l",
     xlim = c(min(c(preds$temp, d$temp)), max(c(preds$temp, d$temp))), ylim = c(min(c(preds$.fitted, d$rate)), max(c(preds$.fitted, d$rate))))


points(temp_anderson_82,  rate_anderson_82, col = 2, pch = 19)
points(temp_prah_77_man,  rate_prah_77_man, col = 3, pch = 19)
points(temp_purnell_66, rate_purnell_66, col = 4, pch = 19)

legend("bottomleft", inset=c(0, 1.01), legend=c(expression(paste(italic('S. mansoni,'), "  R.M. Anderson et al. 1982")), expression(paste(italic('S. mansoni,'), "  S. K. Prah and C. James, 1977")),
                                              expression(paste(italic('S. mansoni,'), "  R.E. Purnell, 1966"))), 
       pch = 19, col = c(2, 3, 4, 5), title="Data", cex = .8)



text(preds$temp[which.min(preds$.fitted)],   0.7,     labels = "^", cex = 1.5, col = 1)
text(preds$temp[which.min(preds$.fitted)],   1.2,     labels = round(preds$temp[which.min(preds$.fitted)],3), cex = 1.5, col = 1)


mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 2.5, cex = 1.5)
mtext(text = "Mortality rate of miracidia per day", side = 2, line = 2.5, cex = 1.5)

#main = "Death rate of miracidia", cex.main=1.7
dev.off()















pdf(file = "paper_mortality_rate_mir.pdf", width = 5, height = 5)
par(mar=c(4, 4, 1, 1), xpd = TRUE)
plot(preds$temp, preds$.fitted, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l", cex.axis=1.3,
     xlim = c(min(c(preds$temp, d$temp)), max(c(preds$temp, d$temp))), ylim = c(min(c(preds$.fitted, d$rate)), max(c(preds$.fitted, d$rate))))


points(temp_anderson_82,  rate_anderson_82, col = 2, pch = 19)
points(temp_prah_77_man,  rate_prah_77_man, col = 3, pch = 19)
points(temp_purnell_66, rate_purnell_66, col = 4, pch = 19)


text(preds$temp[which.min(preds$.fitted)],   0.7,     labels = "^", cex = 1.5, col = 1)
text(preds$temp[which.min(preds$.fitted)],   1,     labels = round(preds$temp[which.min(preds$.fitted)],3), cex = 1.5, col = 1)


mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 2.5, cex = 1.5)
mtext(text = "Mortality rate of miracidia per day", side = 2, line = 2.5, cex = 1.5)

#main = "Death rate of miracidia", cex.main=1.7
dev.off()










