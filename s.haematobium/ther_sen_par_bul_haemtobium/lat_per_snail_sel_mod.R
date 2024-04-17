# This script perform the best curve for sigma_s
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)

#W. Pflüger 1984, Bulinus truncatus with S. haematobium 
temp_pfluger_84 <- c(18,       19,     20,    21,    22,     23,    25,  28,   30,    31,  32)
rate_pfluger_84 <- c(1/117,   1/110,  1/79, 1/60,   1/57,  1/56,  1/44, 1/32, 1/30, 1/31, 1/30)




#R. M. Gordon, 1934, Physopsis (bolinus) globosa  with haematobium 
temp_gordon_34_glo <- c(22,     26.3,      26.4,     26.8,       31.9,    32.1,     33,     35.2)
rate_gordon_34_glo <- c(1/67,  1/38.5,    1/39.3,   1/34.8,     1/27.5,   1/23,   1/25.5,   1/26.3)


#load data, Latent period of snail
temp <- c(temp_pfluger_84, temp_gordon_34_glo)
rate <- c(rate_pfluger_84, rate_gordon_34_glo)

# keep just a single curve
.x <- data.frame(temp, rate)

# fit model

fit <- nls_multstart(rate~lrf_1991(temp = temp, rmax, topt, tmin, tmax),
                     data = .x,
                     iter = c(3,3,3,3),
                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'lrf_1991') - 10,
                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'lrf_1991') + 10,
                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'lrf_1991'),
                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'lrf_1991'),
                     supp_errors = 'Y',
                     convergence_count = FALSE)


fit



# calculate additional traits
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data

new_data <- data.frame(temp = seq(min(.x$temp), max(.x$temp), 0.5))
preds <- augment(fit, newdata = new_data)



pdf(file = "lat_per_snail.pdf", width = 5, height = 5)
par(mar=c(4, 4, 4, 1), xpd = TRUE)
plot(preds$temp, 1/preds$.fitted, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l", 
     xlim = c(min(c(preds$temp, .x$temp)), max(c(preds$temp, .x$temp))), ylim = c(min(c(1/preds$.fitted, 1/.x$rate)), max(c(1/preds$.fitted, 1/.x$rate))))

points(temp_pfluger_84, 1/rate_pfluger_84, col = 2, pch = 19)
points(temp_gordon_34_glo, 1/rate_gordon_34_glo, col = 3, pch = 19)



legend("bottomleft", inset=c(0, 1.01), legend=c(expression(paste(italic('B. truncatus with S. haematobium,'), "  W. Pflüger, 1984")), 
                                             expression(paste(italic('B. globosus with S. haematobium,'), " R. M. Gordon et al., 1934"))),
              pch = 19, col = c(2, 3), title="Data", cex = .8)


mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 2.5, cex = 2)
mtext(text = "Incubation period in snails (days)", side = 2, line = 2, cex = 2)

dev.off()


pdf(file = "paper_lat_per_snail.pdf", width = 5, height = 5)
par(mar=c(6, 6, 4, 1), xpd = TRUE)
plot(preds$temp, 1/preds$.fitted, col = 1, las = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2,
     xlim = c(min(c(preds$temp, .x$temp)), max(c(preds$temp, .x$temp))), ylim = c(min(c(1/preds$.fitted, 1/.x$rate)), max(c(1/preds$.fitted, 1/.x$rate))))

points(temp_pfluger_84, 1/rate_pfluger_84, col = 2, pch = 19)
points(temp_gordon_34_glo, 1/rate_gordon_34_glo, col = 3, pch = 19)


mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 4, cex = 2)
mtext(expression(frac(1, sigma[s])), side = 2, line = 3.5, cex = 2.5, las = 1)
mtext(text = expression(paste(italic('Bulinus'), " spp.")), side = 3, line = 2, cex = 2)

dev.off()