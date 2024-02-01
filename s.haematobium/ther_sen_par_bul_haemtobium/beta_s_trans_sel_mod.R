#This script perform the selected model for beta_s
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)


#(Prah and James 1977), Bulinus (Physopsis) globosus with S. haematobium.
prah_77_temp <- c(       5,        12,            19,     26)
prah_77_rate <- c(0.001279126, 0.04565688,  0.1084776, 0.1849867)
prah_77_rate <- prah_77_rate/max(prah_77_rate)

#(Shiff 1974), Bulinus globosus infected with S. haematobium.
shiff_74_temp <- c(      25,          20,           18,         16,          13,            11) 
shiff_74_rate <- c(0.03730699,  0.003512327,  0.01022082, 0.001782906,  0.002553961,  0.001281865)
shiff_74_rate <- shiff_74_rate/max(shiff_74_rate)

#Chu et al. 1966), Bulinus trancatus  infected with S. haematobium. 
chu_66_temp <- c(    10,            12,           14,         15,        20,          25,         30,        35,         38)
chu_66_rate <- c(0.002009499,  0.008145561,  0.02029471, 0.02016721, 0.04754069, 0.03577785,  0.06475144, 0.04395093, 0.05660857)
chu_66_rate <- chu_66_rate/max(chu_66_rate)           


#The disease transmission rate in snails (Anderson et al. 1982, experiment with Biomphalaria glabra and schistosoma mansoni)
temp <- c(prah_77_temp, shiff_74_temp, chu_66_temp)
rate <- c(prah_77_rate, shiff_74_rate, chu_66_rate)




# keep just a single curve
.x <- data.frame(temp, rate)

# fit model

fit <- nls_multstart(rate~flinn_1991(temp = temp, a, b, c),
                     data = .x,
                     iter = c(5,5,5),
                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') - 10,
                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') + 10,
                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                     supp_errors = 'Y',)


fit


# calculate additional traits
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data

new_data <- data.frame(temp = seq(min(.x$temp), max(.x$temp)+10, 0.5))
preds <- augment(fit, newdata = new_data)


pdf(file = "beta_s_trans_rate.pdf", width = 5, height = 5)
par(mar=c(4, 4, 5, 1), xpd = TRUE)
plot(preds$temp, preds$.fitted, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l", cex.axis=1.3,
     xlim = c(min(c(preds$temp, .x$temp)), max(c(preds$temp, .x$temp))), ylim = c(min(c(preds$.fitted, .x$rate)), max(c(preds$.fitted, .x$rate))))

points(prah_77_temp, prah_77_rate, col = 2, pch = 19)
points(shiff_74_temp, shiff_74_rate, col = 3, pch = 19)
points(chu_66_temp, chu_66_rate, col = 4, pch = 19)



legend("bottomleft", inset=c(0, 1.01), legend=c(expression(paste(italic('B. globosus with S. haematobium,'), "  Prah and James, 1977")), expression(paste(italic('B. globosus with S. haematobium,'), "  Shiff, 1974")),
                                              expression(paste(italic('B. trancatus with S. haematobium,'), "  Chu et al., 1966"))), 
       pch = 19, col = c(2, 3, 4), title="Data", cex = .8)  

text(preds$temp[which.max(preds$.fitted)],   0.03,     labels = "^", cex = 1.5, col = 1)
text(preds$temp[which.max(preds$.fitted)],   0.1,     labels = round(preds$temp[which.max(preds$.fitted)],3), cex = 1.5, col = 1)


mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 2.5, cex = 1.5)
mtext(text = "Transmission rate in snails", side = 2, line = 2.5, cex = 1.5)


#main = "Death rate of miracidia", cex.main=1.7
dev.off()





pdf(file = "paper_beta_s_trans_rate.pdf", width = 5, height = 5)
par(mar=c(6, 6, 4, 1), xpd = TRUE)
plot(preds$temp, preds$.fitted, col = 1, las = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2,
     xlim = c(min(c(preds$temp, .x$temp)), max(c(preds$temp, .x$temp))), ylim = c(min(c(preds$.fitted, .x$rate)), max(c(preds$.fitted, .x$rate))))

points(prah_77_temp, prah_77_rate, col = 2, pch = 19)
points(shiff_74_temp, shiff_74_rate, col = 3, pch = 19)
points(chu_66_temp, chu_66_rate, col = 4, pch = 19)


text(preds$temp[which.max(preds$.fitted)],   0.03,     labels = "^", cex = 1.5, col = 1)
text(preds$temp[which.max(preds$.fitted)],   0.1,     labels = round(preds$temp[which.max(preds$.fitted)],3), cex = 1.5, col = 1)


mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 4, cex = 2)
mtext(text = expression(beta[s]), side = 2, line = 3.5, cex = 2.5, las = 1)
mtext(text = expression(paste(italic('S. haematobium'))), side = 3, line = 2, cex = 2)


dev.off()








