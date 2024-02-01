# This script perform the best curve for sigma_s
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)


#W. Pflüger 1981, Biomphalaria glabrata with S. mansoni, diurnally changing water temperature 
temp_pfluger_81 <- c(24.3,     21.9,   20.2,  19.8,  17.5,   17.9,   19.2,  21.1, 18.9,    25,   29.9,  33)
rate_pfluger_81 <- c(1/26,    1/35.5,  1/44, 1/50.5, 1/68,  1/57.5,  1/52,  1/38, 1/57.5, 1/17, 1/16.5, 1/16.5)

#W. Pflüger 1980, Biomphalaria glabrata with S. mansoni, constant temperature 
temp_pfluger_80 <- c(  17,      19,     22,   25,      28,     30,      31,       32,     33,     34,    35)
rate_pfluger_80 <- c(1/92.5,  1/55.5,  1/34, 1/24.8, 1/19.2,  1/16.3,  1/15.9,  1/14.9, 1/15.3, 1/15.9, 1/16)


#Foster 1964, Biomphalaria pfeifleri with s. mansoni 
temp_foster_64 <- c(18,     21,    22.85,   24.01,   26.26,  28.07,  30.04, 31.75)
rate_foster_64 <- c(1/57,  1/37,    1/32,   1/30,     1/23,   1/19,   1/18,   1/16)


# R. M. Gordon, 1934, Planorbis pfeifferi with S. mansoni
temp_gordon_34_pfei <- c(21.7,    22,       26.3,      26.6,       26.9,     27.3,    27.7,    31.8,    32.1,    32.7,   32.8,   35)
rate_gordon_34_pfei <- c(1/36,   1/33,     1/23.5,     1/23,       1/22.2,   1/22,    1/19,   1/15.4, 1/15.8,  1/17.2,  1/14.7, 1/16.4)


#load data, Latent period of snail
temp <- c(temp_pfluger_81,  temp_pfluger_80, temp_foster_64, temp_gordon_34_pfei)
rate <- c(rate_pfluger_81,  rate_pfluger_80, rate_foster_64, rate_gordon_34_pfei)




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

new_data <- data.frame(temp = seq(min(.x$temp)-2, max(.x$temp), 0.5))
preds <- augment(fit, newdata = new_data)





pdf(file = "lat_per_snail.pdf", width = 5, height = 5)
par(mar=c(4, 4, 6, 1), xpd = TRUE)
plot(preds$temp, 1/preds$.fitted, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 2, type = "l",
     xlim = c(min(c(preds$temp, .x$temp)), max(c(preds$temp, .x$temp))), ylim = c(min(c(1/preds$.fitted, 1/.x$rate)), max(c(1/preds$.fitted, 1/.x$rate))))

points(temp_pfluger_81, 1/rate_pfluger_81, col = 2, pch = 19)
points(temp_pfluger_80, 1/rate_pfluger_80, col = 3, pch = 19)
points(temp_foster_64, 1/rate_foster_64, col = 4, pch = 19)
points(temp_gordon_34_pfei, 1/rate_gordon_34_pfei, col = 5, pch = 19)


legend("bottomleft", inset=c(0, 1.01), legend=c(expression(paste(italic('B. glabrata with S. mansoni,'), "  W. Pflüger, 1981")), expression(paste(italic('B. glabrata with S. mansoni,'), "  W. Pflüger, 1980")),
                                             expression(paste(italic('B. pfeifferi with S. mansoni,'), "  R. Foster, 1964")), expression(paste(italic('B. pfeifferi with S. mansoni,'), "  R. M. Gordon et al., 1934"))), 
       pch = 19, col = c(2, 3, 4, 5), title="Data", cex = .8)


mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 2.5, cex = 1.5)
mtext(text = "Incubation period in snails (days)", side = 2, line = 2.5, cex = 1.5)

dev.off()




pdf(file = "paper_lat_per_snail.pdf", width = 5, height = 5)
par(mar=c(6, 6, 4, 1), xpd = TRUE)
plot(preds$temp, 1/preds$.fitted, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 2, type = "l",  cex.axis=2,
     xlim = c(min(c(preds$temp, .x$temp)), max(c(preds$temp, .x$temp))), ylim = c(min(c(1/preds$.fitted, 1/.x$rate)), max(c(1/preds$.fitted, 1/.x$rate))), las = 1)

points(temp_pfluger_81, 1/rate_pfluger_81, col = 2, pch = 19)
points(temp_pfluger_80, 1/rate_pfluger_80, col = 3, pch = 19)
points(temp_foster_64, 1/rate_foster_64, col = 4, pch = 19)
points(temp_gordon_34_pfei, 1/rate_gordon_34_pfei, col = 5, pch = 19)

mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 4, cex = 2)
mtext(expression(frac(1, sigma[s])), side = 2, line = 3.5, cex = 2.5, las = 1)
mtext(text = expression(paste(italic('Biomphalaria'), " spp.")), side = 3, line = 2, cex = 2)


dev.off()