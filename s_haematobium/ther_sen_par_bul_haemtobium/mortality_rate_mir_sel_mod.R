# This script perform the best curve for mu_m
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)


#load data 
#S. K. Prah and C. James, 1977 S. haematobium

temp_prah_77_haem <- c(7.5,            20,      27.5,       36.5)
rate_prah_77_haem <- c(2.275331,    2.210028,  2.24239,  4.354702)


temp <- temp_prah_77_haem
rate <- rate_prah_77_haem
  
# keep just a single curve
.x  <- data.frame(temp, rate)

# fit model

fit <- nls_multstart(rate~spain_1982(temp = temp, a,b,c,r0),
                     data = .x,
                     iter = c(4,4,4,4),
                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') - 1,
                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') + 1,
                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
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
par(mar=c(4, 4, 4, 1), xpd = TRUE)
plot(preds$temp, preds$.fitted, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l", cex.axis=1.3,
     xlim = c(min(c(preds$temp, .x$temp)), max(c(preds$temp, .x$temp))), ylim = c(min(c(preds$.fitted, .x$rate))-.4, max(c(preds$.fitted, .x$rate))))


points(temp_prah_77_haem,  rate_prah_77_haem, col = 2, pch = 19)

legend("bottomleft", inset=c(0, 1), legend=c(expression(paste(italic('S. haematobium,'), "  S. K. Prah and C. James, 1977"))), 
       pch = 19, col = 2, title="Data", cex = .8)

text(preds$temp[which.min(preds$.fitted)],   1.8,     labels = "^", cex = 1.5, col = 1)
text(preds$temp[which.min(preds$.fitted)],   1.95,     labels = round(preds$temp[which.min(preds$.fitted)],3), cex = 1.5, col = 1)

mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 2.5, cex = 1.5)
mtext(text = "Mortality rate of miracidia per day", side = 2, line = 2.5, cex = 1.5)

#main = "Death rate of miracidia", cex.main=1.7
dev.off()




pdf(file = "paper_mortality_rate_mir.pdf", width = 5, height = 5)
par(mar=c(4, 4, 1, 1), xpd = TRUE)
plot(preds$temp, preds$.fitted, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l", cex.axis=1.3,
     xlim = c(min(c(preds$temp, .x$temp)), max(c(preds$temp, .x$temp))), ylim = c(min(c(preds$.fitted, .x$rate))-.4, max(c(preds$.fitted, .x$rate))))


points(temp_prah_77_haem,  rate_prah_77_haem, col = 2, pch = 19)


text(preds$temp[which.min(preds$.fitted)],   1.8,     labels = "^", cex = 1.5, col = 1)
text(preds$temp[which.min(preds$.fitted)],   1.95,     labels = round(preds$temp[which.min(preds$.fitted)],3), cex = 1.5, col = 1)

mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 2.5, cex = 1.5)
mtext(text = "Mortality rate of miracidia per day", side = 2, line = 2.5, cex = 1.5)

#main = "Death rate of miracidia", cex.main=1.7
dev.off()






