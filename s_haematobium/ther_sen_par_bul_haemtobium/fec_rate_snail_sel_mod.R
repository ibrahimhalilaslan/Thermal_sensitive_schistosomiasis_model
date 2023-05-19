# This script perform the best model for nu_c
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)



#Chester Kalinda. 2017, experiment with Bulinus globosus. 
#(Note Age of maturity and hatching rate getting from El-Hassan 1974 for buliunus)
temp_ches_17 <- c(15.5,    21.2,   25.8,    31, 35.5)
rate_ches_17 <- c(  0,    0.027,  0.037,  0.019, 0)


#El- Hassan 1974, Bulinus truncatus
temp_has_74_t <- c(15,        20,       25,       30,    35)
rate_has_74_t <- c(0.0103,   0.025,   0.0348,  0.0199,    0)


#(Shiff 1967), Bulinus globosus
temp_shiff_67 <- c(  18,      22,      25,      27)
rate_shiff_67 <- c(0.0137,  0.0447,  0.0497,  0.0391)


#(Kubirizajournal et al. 2010), Bulinus nyassanus
temp_kubir_10 <- c(  22,     25,    28,    31)
rate_kubir_10 <- c(0.0297, 0.0427, 0.0463, 0.0451)

# We assume there is no egg production at 32, 33, 34 


temp <-c(temp_ches_17, temp_has_74_t,  temp_shiff_67, temp_kubir_10,  32, 33, 34)
rate <-c(rate_ches_17, rate_has_74_t,  rate_shiff_67, rate_kubir_10,   0, 0, 0)

# keep just a single curve
.x <- data.frame(temp, rate)

# fit model
fit <- nls_multstart(rate~johnsonlewin_1946(temp = temp, r0, e, eh, topt),
                     data = .x,
                     iter = c(4,4,4,4),
                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'johnsonlewin_1946') - 1,
                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'johnsonlewin_1946') + 1,
                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'johnsonlewin_1946'),
                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'johnsonlewin_1946'),
                     supp_errors = 'Y',  
                     convergence_count = FALSE)
fit


# calculate additional traits
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data

new_data <- data.frame(temp = seq(min(.x$temp)-5, max(.x$temp), 0.1))
preds <- augment(fit, newdata = new_data)



pdf(file = "fecundity_rate_snail.pdf", width = 5, height = 5)
par(mar=c(4, 4, 6, 1), xpd = TRUE)
plot(preds$temp, preds$.fitted, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l",
      xlim = c(min(preds$temp), max(preds$temp)), ylim = c(0, 0.05))

points(temp_ches_17, rate_ches_17, col = 2, pch = 19)
points(temp_has_74_t, rate_has_74_t, col = 3, pch = 19)
points(temp_shiff_67, rate_shiff_67, col = 4, pch = 19)
points(temp_kubir_10, rate_kubir_10, col = 5, pch = 19)



legend("bottomleft", inset=c(0, 1.01), legend=c(expression(paste(italic('B. Globosus,'),"  Kalinda et al., 2017")), 
                                              expression(paste(italic('B. Truncatus,'),"  A. A. EL-HASSAN, 1974")),
                                             expression(paste(italic('B. globosus,'),"  Shiff,  1967")), 
                                             expression(paste(italic('B. nyassanus,'),"  Kubirizajournal et al., 2010"))),
       pch = 19, col = c(2, 3, 4, 5), title = "Data", cex = .8)

text(preds$temp[which.max(preds$.fitted)],    0,     labels = "^", cex = 1.5, col = 1)
text(preds$temp[which.max(preds$.fitted)],    0.003,     labels = round(preds$temp[which.max(preds$.fitted)],3), cex = 1.5, col = 1)


mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 2.5, cex = 1.5)
mtext(text = "Fecundity rate of snails per day", side = 2, line = 2.5, cex = 1.5)

#main = "The fecundity rate of snail per day", cex.main=1.7
dev.off()




pdf(file = "paper_fecundity_rate_snail.pdf", width = 5, height = 5)
par(mar=c(4, 4, 1, 1), xpd = TRUE)
plot(preds$temp, preds$.fitted, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l", cex.axis=1.3,
     xlim = c(min(preds$temp), max(preds$temp)), ylim = c(0, 0.05))

points(temp_ches_17, rate_ches_17, col = 2, pch = 19)
points(temp_has_74_t, rate_has_74_t, col = 3, pch = 19)
points(temp_shiff_67, rate_shiff_67, col = 4, pch = 19)
points(temp_kubir_10, rate_kubir_10, col = 5, pch = 19)



text(preds$temp[which.max(preds$.fitted)],    0,     labels = "^", cex = 1.5, col = 1)
text(preds$temp[which.max(preds$.fitted)],    0.003,     labels = round(preds$temp[which.max(preds$.fitted)],3), cex = 1.5, col = 1)


mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 2.5, cex = 1.5)
mtext(text = "Fecundity rate of snails per day", side = 2, line = 2.5, cex = 1.5)

#main = "The fecundity rate of snail per day", cex.main=1.7
dev.off()


