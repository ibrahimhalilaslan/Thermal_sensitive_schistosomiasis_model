#This script perform the selected model for beta_h
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)

#(Purnell 1966),, S. mansoni
purnell_66_temp <- c(       12,             15,           18,             21,        24,          27,           30,       33)
purnell_66_rate_1 <- c(  0.9246078,      0.810981,    1.320736,       1.092269,   0.6998873,   1.231836,    0.7901292,  0.6292689)
purnell_66_rate_1 <- purnell_66_rate_1/max(purnell_66_rate_1)


purnell_66_temp <- c(       12,             15,           18,             21,        24,          27,           30,       33)
purnell_66_rate_2 <- c(  1.0535,         1.002849,     0.9284408,       1.481275,   0.7225294,   0.6554545,    1.251547,  0.7475422)
purnell_66_rate_2 <- purnell_66_rate_2/max(purnell_66_rate_2)


purnell_66_temp <- c(       12,             15,           18,             21,        24,          27,           30,       33)
purnell_66_rate_3 <- c(  2.182979,      1.523512,    1.777458,       1.840471,   0.8305635,   0.7054314,      1.339291,  0)
purnell_66_rate_3 <- purnell_66_rate_3/max(purnell_66_rate_3)



#purnell_66_rate <- purnell_66_rate/max(purnell_66_rate)


#(Stirewalt 1954), S. mansoni 
stirewalt_54_temp <- c(    24,         27)
stirewalt_54_rate <- c(0.09590499,  0.1552625)
stirewalt_54_rate <- stirewalt_54_rate/max(stirewalt_54_rate)


#(DeWitt 1965) with S. mansoni
dewitt_55_temp_1 <- c(     15,      20,          25,      30,      35)
dewitt_55_rate_1 <- c(0.4618475, 0.8464159,  1.00836, 1.237522, 0.6192221)
dewitt_55_rate_1 <- dewitt_55_rate_1/max(dewitt_55_rate_1)

dewitt_55_temp_2 <- c(     15,      20,          25,      30,         35,        40)
dewitt_55_rate_2 <- c(0.385854, 0.6910937,    0.635236,  0.9783329, 1.0652, 0.2764123)
dewitt_55_rate_2 <- dewitt_55_rate_2/max(dewitt_55_rate_2)

dewitt_55_temp_3 <- c(   10,        15,      20,         25,      30,       35,        40,   45)
dewitt_55_rate_3 <- c(0.7439515,  1.76077,  2.484149,  3.50177,  3.9655, 3.421058,  0.7532138, 0)
dewitt_55_rate_3 <- dewitt_55_rate_3/max(dewitt_55_rate_3)

dewitt_55_temp_4 <- c(0, 5,        10,        15,      20,        25)
dewitt_55_rate_4 <- c(0, 0,    0.1190698, 0.6479505, 1.39834, 1.50036)
dewitt_55_rate_4 <- dewitt_55_rate_4/max(dewitt_55_rate_4) 


# load data
#The transmission rate of schisto in human
temp <- c( purnell_66_temp,   purnell_66_temp,   purnell_66_temp,  stirewalt_54_temp, dewitt_55_temp_1, dewitt_55_temp_2, dewitt_55_temp_3, dewitt_55_temp_4)
rate <- c(purnell_66_rate_1, purnell_66_rate_2, purnell_66_rate_3, stirewalt_54_rate, dewitt_55_rate_1, dewitt_55_rate_2, dewitt_55_rate_3, dewitt_55_rate_4)


# keep just a single curve
.x <- data.frame(temp, rate)

# fit model

fit <- nls_multstart(rate~briere2_1999(temp = temp, tmin, tmax, a,b),
                     data = .x,
                     iter = c(4,4,4,4),
                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') - 10,
                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') + 10,
                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                     supp_errors = 'Y')

fit


# calculate additional traits
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data

new_data <- data.frame(temp = seq(min(.x$temp), max(.x$temp), 0.2))
preds <- augment(fit, newdata = new_data)


pdf(file = "beta_h_trans_rate.pdf", width = 5, height = 5)
par(mar=c(4, 4, 5, 1), xpd = TRUE)
plot(preds$temp, preds$.fitted, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l", cex.axis=1.3,
     xlim = c(min(c(preds$temp, .x$temp)), max(c(preds$temp, .x$temp))), ylim = c(min(c(preds$.fitted, .x$rate)), max(c(preds$.fitted, .x$rate))))

points(c(purnell_66_temp,purnell_66_temp, purnell_66_temp),  c(purnell_66_rate_1,  purnell_66_rate_2,  purnell_66_rate_3), col = 2, pch = 19)
points(stirewalt_54_temp, stirewalt_54_rate, col = 3, pch = 19) 
points(c(dewitt_55_temp_1, dewitt_55_temp_2, dewitt_55_temp_3, dewitt_55_temp_4),  c(dewitt_55_rate_1, dewitt_55_rate_2, dewitt_55_rate_3, dewitt_55_rate_4), col = 4, pch = 19)


legend("bottomleft", inset=c(0, 1), legend=c(expression(paste(italic('S. mansoni,'), "  Purnell, 1966")),
                                             expression(paste(italic('S. mansoni,'), "  Stirewalt, 1954")),
                                             expression(paste(italic('S. mansoni,'), "   DeWitt. 1965"))),
                                             pch = 19, col = c(2, 3, 4), title="Data", cex = .8)

text(preds$temp[which.max(preds$.fitted)],   0,     labels = "^", cex = 1.5, col = 1)
text(preds$temp[which.max(preds$.fitted)],   0.07,     labels = round(preds$temp[which.max(preds$.fitted)],3), cex = 1.5, col = 1)


mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 2.5, cex = 1.5)
mtext(text = "Transmission rate in humans", side = 2, line = 2.5, cex = 1.5)

#main = "Death rate of miracidia", cex.main=1.7
dev.off()




pdf(file = "paper_beta_h_trans_rate.pdf", width = 5, height = 5)
par(mar=c(6, 6, 4, 1), xpd = TRUE)
plot(preds$temp, preds$.fitted, col = 1, las = 1, pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2,
     xlim = c(min(c(preds$temp, .x$temp)), max(c(preds$temp, .x$temp))), ylim = c(min(c(preds$.fitted, .x$rate)), max(c(preds$.fitted, .x$rate))))


points(c(purnell_66_temp,purnell_66_temp, purnell_66_temp),  c(purnell_66_rate_1,  purnell_66_rate_2,  purnell_66_rate_3), col = 2, pch = 19)
points(stirewalt_54_temp, stirewalt_54_rate, col = 3, pch = 19) 
points(c(dewitt_55_temp_1, dewitt_55_temp_2, dewitt_55_temp_3, dewitt_55_temp_4),  c(dewitt_55_rate_1, dewitt_55_rate_2, dewitt_55_rate_3, dewitt_55_rate_4), col = 4, pch = 19)


text(preds$temp[which.max(preds$.fitted)],   0,     labels = "^", cex = 1.5, col = 1)
text(preds$temp[which.max(preds$.fitted)],   0.05,     labels = round(preds$temp[which.max(preds$.fitted)],3), cex = 1.5, col = 1)


mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 4, cex = 2)
mtext(text = expression(beta[h]), side = 2, line = 3.6, cex = 2.5, las = 1)
mtext(text = expression(paste(italic('S. mansoni'))), side = 3, line = 2, cex = 2)


dev.off()





