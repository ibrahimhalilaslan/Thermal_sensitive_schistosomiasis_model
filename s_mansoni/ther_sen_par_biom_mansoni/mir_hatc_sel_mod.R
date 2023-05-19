# The perform the best model for delta_e 
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)

# load in data
#Miracidia hatching rate  (Nguyen et al. 2021, Schistosoma mansoni)
temp <- c(5,  9,  13, 17,  21, 22, 25,  29,  33, 37)
rate <- c(70, 60, 90, 100, 76, 114, 137, 105, 150, 135)/300


# keep just a single curve
.x <- data.frame(temp, rate)

# fit model

fit <-nls_multstart(rate~flinn_1991(temp = temp, a, b, c),
                    data = .x,
                    iter = c(5,5,5),
                    start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') - 10,
                    start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') + 10,
                    lower = get_lower_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                    upper = get_upper_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                    supp_errors = 'Y')


fit


# calculate additional traits
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data

new_data <- data.frame(temp = seq(min(.x$temp), max(.x$temp), 0.5))
preds <- augment(fit, newdata = new_data)


pdf(file = "mir_hatch_rate.pdf", width = 5, height = 5)
par(mar=c(4, 4, 4, 1), xpd = TRUE)
plot(preds$temp, preds$.fitted, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.axis=1.3, cex.lab = 1.3, type = "l",
     xlim = c(min(c(preds$temp, .x$temp)), max(c(preds$temp, .x$temp))), ylim = c(min(c(preds$.fitted, .x$rate)), max(c(preds$.fitted, .x$rate))))

points(.x$temp, .x$rate, col = 2, pch = 19)

legend("bottomleft", inset=c(0, 1.01), legend=c(expression(paste(italic('S. mansoni,'), "  Nguyen et al., 2021"))),  pch = 19, col = 2, title="Data",  cex = .8)


mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 2.5, cex = 1.5)
mtext(text = "Prob. of miracidia hatching success", side = 2, line = 2.5, cex = 1.5)

#main = "The probab. of miracidia hatching", cex.main=1.7

dev.off()


pdf(file = "paper_mir_hatch_rate.pdf", width = 5, height = 5)
par(mar=c(4, 4, 1, 1), xpd = TRUE)
plot(preds$temp, preds$.fitted, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l",cex.axis=1.3,
     xlim = c(min(c(preds$temp, .x$temp)), max(c(preds$temp, .x$temp))), ylim = c(min(c(preds$.fitted, .x$rate)), max(c(preds$.fitted, .x$rate))))

points(d$temp, d$rate, col = 2, pch = 19)


mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 2.5,  cex = 1.5)
mtext(text = "Prob. of miracidia hatching success", side = 2, line = 2.5,  cex = 1.5)

#main = "The probab. of miracidia hatching", cex.main=1.7

dev.off()


