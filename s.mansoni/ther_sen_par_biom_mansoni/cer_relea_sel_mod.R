# This script perform the best curve for nu_c
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)

#Nguyen et al. 2021, S. mansoni with Biomphalaria glabrata
temp_nguyen_21 <- c(5,  9,     13,     17,     21,     25,     29,     33,     37)
rate_nguyen_21 <- c(0, 15,     260,   1350,   1970,   2900,   1970,   1770,   690)

#(Upatham 1973), B. glabrata with S. mansoni  
temp_upatham_73 <- c(10, 13, 16,   19, 22,    25,   28,  31,   34,   37, 40)
rate_upatham_73 <- c( 0,  0, 206, 114, 442, 1200, 1400, 1233, 1241, 243, 121)



temp <- c(temp_nguyen_21, temp_upatham_73)
rate <- c(rate_nguyen_21, rate_upatham_73)



# keep just a single curve
.x <- data.frame(temp, rate)

# fit model

fit <- nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                     data = .x,
                     iter = c(4,4,4),
                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                     supp_errors = 'Y',
                     convergence_count = FALSE)
fit


# calculate additional traits
calc_params(fit) %>%
  # round for easy viewing
mutate_all(round, 2)

# predict new data

new_data <- data.frame(temp = seq(min(.x$temp), max(.x$temp)+5, 0.5))
preds <- augment(fit, newdata = new_data)


pdf(file = "cercariae_release_rate.pdf", width = 5, height = 5)
par(mar=c(4, 4, 5, 1), xpd = TRUE)
plot(preds$temp, preds$.fitted, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l",  cex.axis=1.3,
     xlim = c(min(c(preds$temp, .x$temp)), max(c(preds$temp, .x$temp))), ylim = c(min(c(preds$.fitted, .x$rate)), max(c(preds$.fitted, .x$rate))))

points(temp_nguyen_21, rate_nguyen_21, col = 2, pch = 19)
points(temp_upatham_73, rate_upatham_73, col = 3, pch = 19)

legend("bottomleft", inset=c(0, 1.01), legend=c(expression(paste(italic('B. glabrata with S. mansoni,'), "  Nguyen et al., 2021")),
                                             expression(paste(italic('B. glabrata with S. mansoni,'), "  Upatham, 1973"))), 
       pch = 19, col = c(2, 3), title="Data", cex = .8)

text(preds$temp[which.max(preds$.fitted)],  0,     labels = "^", cex = 1.5, col = 1)
text(preds$temp[which.max(preds$.fitted)],  200,     labels = round(preds$temp[which.max(preds$.fitted)],3), cex = 1.5, col = 1)


mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 2.5, cex = 1.5)
mtext(text = "Cercarial/sail/day", side = 2, line = 2.5, cex = 1.5)
#main = "The num. of cercariae realease", cex.main=2

dev.off()


pdf(file = "paper_cercariae_release_rate.pdf", width = 5, height = 5)
par(mar=c(6, 7, 4, 1), xpd = TRUE)
plot(preds$temp, preds$.fitted, col = 1, las = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 2, type = "l",  cex.axis=2,
     xlim = c(min(c(preds$temp, .x$temp)), max(c(preds$temp, .x$temp))), ylim = c(min(c(preds$.fitted, .x$rate)), max(c(preds$.fitted, .x$rate))))

points(temp_nguyen_21, rate_nguyen_21, col = 2, pch = 19)
points(temp_upatham_73, rate_upatham_73, col = 3, pch = 19)

text(preds$temp[which.max(preds$.fitted)],  0,     labels = "^", cex = 1.5, col = 1)
text(preds$temp[which.max(preds$.fitted)],  200,     labels = round(preds$temp[which.max(preds$.fitted)],3), cex = 1.5, col = 1)


mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 4, cex = 2)
mtext(text = expression(nu[c]), side = 2, line = 5, cex = 2.5, las = 1)
mtext(text = expression(paste(italic('S. mansoni'))), side = 3, line = 2, cex = 2)

dev.off()