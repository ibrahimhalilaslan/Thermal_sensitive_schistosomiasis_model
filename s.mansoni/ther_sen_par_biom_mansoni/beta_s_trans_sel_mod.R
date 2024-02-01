#This script perform the selected model for beta_s
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)



# (Stirewalt 1954) Australorbis glabratus with S. mansoni.
stirewalt_54_temp <- c(    24,         27,        32)
stirewalt_54_rate <- c(0.07515584, 0.4808139, 0.3922359)
stirewalt_54_rate <- stirewalt_54_rate/max(stirewalt_54_rate)

#(Foster 1964), Biomphalaria pfeifferi with S. mansoni.
foster_64_temp <- c(    22.85,      24.01,       26.26,      28.07,       30.04,    31.75)
foster_64_rate <- c(0.04079299,   0.04333207,  0.2203394, 0.08210066,  0.05756868, 0.07958718)
foster_64_rate <- foster_64_rate/max(foster_64_rate)

#(Prah and James 1977),  Biomphalaria pfeifferi with S. mansoni. 
prah_77_temp <- c(    5,             15,         20,         22.5)
prah_77_rate <- c(0.0006344901, 0.05059652, 0.1629059,  0.2249733)
prah_77_rate <- prah_77_rate/max(prah_77_rate)

#(Anderson et al. 1982), Biomphalaria glabrata with S. mansoni 
anderson_82_temp <- c(15,               20,          25,           30,         35)
anderson_82_rate <- c(0.3658224,    0.3865833,    0.9074138,     0.487644, 0.4107836)
anderson_82_rate <-  anderson_82_rate/max(anderson_82_rate)

#(DeWitt 1955), B. glabrata with S. mansoni. 
dewitt_55_temp <-   c(10,      25,       35,    40)
dewitt_55_rate_1 <- c(0,   0.0078523, 0.0283053, 0)
dewitt_55_rate_1 <- dewitt_55_rate_1/max(dewitt_55_rate_1)
dewitt_55_rate_2 <- c(0,   0.032993,  0.0241807, 0)
dewitt_55_rate_2 <- dewitt_55_rate_2/max(dewitt_55_rate_2)
dewitt_55_rate_3 <- c(0,   0.01388562,  0.03692433, 0.01246124)
dewitt_55_rate_3 <- dewitt_55_rate_3/max(dewitt_55_rate_3)
dewitt_55_rate_4 <- c(0,   0.07405795,  0.09039597, 0.00401955)
dewitt_55_rate_4 <- dewitt_55_rate_4/max(dewitt_55_rate_4)

#Coelho and Bezerra 2006), Biomphalari glabrata.
coelho_temp_06 <- c(15, 20, 30)
coelho_rate_06 <- c(0.0001342012, 0.007285127, 0.04153248)
coelho_rate_06 <- coelho_rate_06/max(coelho_rate_06)

#(Upatham 1973), B. glabrata with S. mansoni. 
upatham_temp_73 <- c(10,      13,     16,           19,           22,           25,        28,          31,          34,          37,          40)
upatham_rate_73 <- c(0,        0,  0.002895729, 0.004360468,  0.009807099, 0.01299096, 0.01373957,  0.01528639,  0.02722317, 0.002743025,  0.00199986)
upatham_rate_73 <- upatham_rate_73/max(upatham_rate_73)

#The disease transmission rate in snails (Anderson et al. 1982, experiment with Biomphalaria glabra and schistosoma mansoni)
temp <- c(stirewalt_54_temp, foster_64_temp, prah_77_temp,  anderson_82_temp, dewitt_55_temp,   dewitt_55_temp,   dewitt_55_temp,   dewitt_55_temp,  coelho_temp_06, upatham_temp_73)
rate <- c(stirewalt_54_rate, foster_64_rate, prah_77_rate,  anderson_82_rate, dewitt_55_rate_1, dewitt_55_rate_2, dewitt_55_rate_3, dewitt_55_rate_4,  coelho_rate_06, upatham_rate_73)

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
                     supp_errors = 'Y')


fit


# calculate additional traits
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data

new_data <- data.frame(temp = seq(min(.x$temp), max(.x$temp), 0.5))
preds <- augment(fit, newdata = new_data)


pdf(file = "beta_s_trans_rate.pdf", width = 5, height = 5)
par(mar=c(4, 4, 8, 1), xpd = TRUE)
plot(preds$temp, preds$.fitted, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l", cex.axis=1.3,
     xlim = c(min(c(preds$temp, .x$temp)), max(c(preds$temp, .x$temp))), ylim = c(min(c(preds$.fitted, .x$rate)), max(c(preds$.fitted, .x$rate))))

points(stirewalt_54_temp,  stirewalt_54_rate, col = 2, pch = 19)
points(foster_64_temp,  foster_64_rate, col = 3, pch = 19)
points(prah_77_temp, prah_77_rate, col = 4, pch = 19)
points(anderson_82_temp, anderson_82_rate, col = 5, pch = 19)
points(c(dewitt_55_temp, dewitt_55_temp,dewitt_55_temp,dewitt_55_temp), c(dewitt_55_rate_1, dewitt_55_rate_2, dewitt_55_rate_3, dewitt_55_rate_4), col = 6, pch = 19)
points(coelho_temp_06, coelho_rate_06, col = 7, pch = 19)
points(upatham_temp_73, upatham_rate_73, col = 8, pch = 19)

legend("bottomleft", inset=c(0, 1.01), legend=c(expression(paste(italic('B. glabrata with S. mansoni,'), "  Stirewalt, 1954")), expression(paste(italic('B. pfeifferi with S. mansoni,'), "  Foster, 1964")),
                                             expression(paste(italic('B. pfeifferi with S. mansoni,'), "  Prah and James, 1977")), expression(paste(italic('B. glabrata with S. mansoni,'), "  Anderson et al., 1982")),
                                             expression(paste(italic('B. glabrata with S. mansoni,'), "  DeWitt, 1955")),expression(paste(italic('B. glabrata. with S. mansoni,'), "  Coelho and Bezerra, 2006")),
                                             expression(paste(italic('B. glabrata with S. mansoni,'), "  Upatham, 1973"))),  
       pch = 19, col = c(2, 3, 4, 5, 6, 7), title="Data", cex = .8)

text(preds$temp[which.max(preds$.fitted)],    0.02,     labels = "^", cex = 1.5, col = 1)
text(preds$temp[which.max(preds$.fitted)],    0.1,     labels = round(preds$temp[which.max(preds$.fitted)],3), cex = 1.5, col = 1)


mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 2.5, cex = 1.5)
mtext(text = "Transmission rate in snails", side = 2, line = 2.5, cex = 1.5)


#main = "Death rate of miracidia", cex.main=1.7
dev.off()





pdf(file = "paper_beta_s_trans_rate.pdf", width = 5, height = 5)
par(mar=c(6, 6, 4, 1), xpd = TRUE)
plot(preds$temp, preds$.fitted, col = 1, las = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2,
     xlim = c(min(c(preds$temp, .x$temp)), max(c(preds$temp, .x$temp))), ylim = c(min(c(preds$.fitted, .x$rate)), max(c(preds$.fitted, .x$rate))))

points(stirewalt_54_temp,  stirewalt_54_rate, col = 2, pch = 19)
points(foster_64_temp,  foster_64_rate, col = 3, pch = 19)
points(prah_77_temp, prah_77_rate, col = 4, pch = 19)
points(anderson_82_temp, anderson_82_rate, col = 5, pch = 19)
points(c(dewitt_55_temp, dewitt_55_temp,dewitt_55_temp,dewitt_55_temp), c(dewitt_55_rate_1, dewitt_55_rate_2, dewitt_55_rate_3, dewitt_55_rate_4), col = 6, pch = 19)
points(coelho_temp_06, coelho_rate_06, col = 7, pch = 19)
points(upatham_temp_73, upatham_rate_73, col = 8, pch = 19)


text(preds$temp[which.max(preds$.fitted)],   0.02,     labels = "^", cex = 1.5, col = 1)
text(preds$temp[which.max(preds$.fitted)],   0.1,     labels = round(preds$temp[which.max(preds$.fitted)],3), cex = 1.5, col = 1)



mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 4, cex = 2)
mtext(text = expression(beta[s]), side = 2, line = 3.7, cex = 2.5, las = 1)
mtext(text = expression(paste(italic('S. mansoni'))), side = 3, line = 2, cex = 2)



dev.off()








