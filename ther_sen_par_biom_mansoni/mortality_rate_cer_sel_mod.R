#This script perform the selected model for mu_c
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)


#J.R. Lawson and R.A.Wilson 1980, experiment with S. mansoni 
temp_lawson_80 <- c(10,          15,           20,      25,        30,        35,       40)
rate_lawson_80 <- c(0.6844903, 0.347094,   0.9361688, 0.9384906, 1.767248, 3.79802, 7.210983)



#R.E. Purnell, 1966, experiment with S. mansoni
temp_purnell_66 <- c(12,             15,        18,            21,     24,           27,       30)
rate_purnell_66 <- c(0.9958894, 2.150529,  1.729903,       1.026659, 1.085069,  1.230383,  1.718108)



# load in data
temp <- c(temp_lawson_80, temp_purnell_66)
rate <- c(rate_lawson_80, rate_purnell_66)

# keep just a single curve
.x <- data.frame(temp, rate)

# fit model

fit <- nls_multstart(rate~spain_1982(temp = temp, a,b,c,r0),
                     data = .x,
                     iter = c(4,4,4,4),
                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') - 1,
                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') + 1,
                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                     supp_errors = 'Y',
                     convergence_count = FALSE)


fit

d <- .x


# calculate additional traits
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data

new_data <- data.frame(temp = seq(min(d$temp-5), max(d$temp), 0.5))
preds <- augment(fit, newdata = new_data)





pdf(file = "mortality_rate_cer.pdf", width = 5, height = 5)
par(mar=c(4, 4, 5, 1), xpd = TRUE)
plot(preds$temp, preds$.fitted, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l", cex.axis=1.3,
     xlim = c(min(c(preds$temp, d$temp)), max(c(preds$temp, d$temp))), ylim = c(min(c(preds$.fitted, d$rate))-1.2, max(c(preds$.fitted, d$rate))))


points(temp_lawson_80, rate_lawson_80, col = 2, pch = 19)
points(temp_purnell_66, rate_purnell_66, col = 3, pch = 19)


legend("bottomleft", inset=c(0, 1.01), legend=c(expression(paste(italic('S. mansoni,'), "  J.R. Lawson and R.A.Wilson 1980")),expression(paste(italic('S. mansoni,'), "  R.E. Purnell, 1966"))), 
       pch = 19, col = c(2, 3), title="Data", cex = .8)


text(preds$temp[which.min(preds$.fitted)],   -0.8,     labels = "^", cex = 1.5, col = 1)
text(preds$temp[which.min(preds$.fitted)],   -0.2,     labels = round(preds$temp[which.min(preds$.fitted)],3), cex = 1.5, col = 1)

mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1,line = 2.5, cex = 1.5)
mtext(text = "Mortality rate of cercariae per day", side = 2, line = 2.5, cex = 1.5)
# main = "The death rate of cercariae", cex.main=1.7
dev.off()


pdf(file = "paper_mortality_rate_cer.pdf", width = 5, height = 5)
par(mar=c(4, 4, 1, 1), xpd = TRUE)
plot(preds$temp, preds$.fitted, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l", cex.axis=1.3,
     xlim = c(min(c(preds$temp, d$temp)), max(c(preds$temp, d$temp))), ylim = c(min(c(preds$.fitted, d$rate))-1, max(c(preds$.fitted, d$rate))))


points(temp_lawson_80, rate_lawson_80, col = 2, pch = 19)
points(temp_purnell_66, rate_purnell_66, col = 3, pch = 19)


text(preds$temp[which.min(preds$.fitted)],   -.6,     labels = "^", cex = 1.5, col = 1)
text(preds$temp[which.min(preds$.fitted)],   -0.2,     labels = round(preds$temp[which.min(preds$.fitted)],3), cex = 1.5, col = 1)

mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1,line = 2.5, cex = 1.5)
mtext(text = "Mortality rate of cercariae per day", side = 2, line = 2.5, cex = 1.5)
# main = "The death rate of cercariae", cex.main=1.7
dev.off()