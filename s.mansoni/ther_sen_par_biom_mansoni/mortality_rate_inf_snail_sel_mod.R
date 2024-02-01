# This perform the best curve for mu
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)


#W. Pflüger 1981, Biomphalaria glabrata with S. mansoni
temp_pfluger_81 <- c(24.3,            21.9,          20.2,       19.7)
rate_pfluger_81 <- c(0.02478296,   0.02157024,   0.02852877,  0.02473841)

#W. Pflüger 1981, Biomphalaria glabrata with S. mansoni with sinual function
temp_pfluger_81_sin <- c(17.5,           17.9,         19.2,            21.1,        18.9,      25,           29.9,     33.1)
rate_pfluger_81_sin <- c(0.00524522,  0.01841793,  0.00982357,   0.008281862,  0.02503713,  0.01368767, 0.02098088, 0.03903461)


#W. Pflüger 1980, Biomphalaria glabrata with S. mansoni, constant temperature (This result from prepatent snails)
temp_pfluger_80 <- c(16,               17,            18,              19,       22,         25,            28,      30,            31,           32,         33,        34,       35)
rate_pfluger_80 <- c(0.02354884,   0.01456296,  0.01980421,      0.01161004, 0.02098088,  0.01863046,  0.01136229, 0.01292767,  0.01194658,  0.009346447, 0.01923994, 0.07677861, 0.04147427)


#R. Foster, 1964 experiment with Biomphalaria pfeifferi invasion by S. mansoni
temp_foster_64 <- c(22.85,           24.01,        26.26,      28.07)
rate_foster_64 <- c(0.006674626,  0.01155245,    0.0260108,   0.03678792)

#(Upatham 1973), B. glabrata with S. mansoni  
temp_upatham_73 <- c(      10,          13,          16,        19,          22,          25,           28,          31,          34,          37,      40)
rate_upatham_73 <- c(0.006316809, 0.009430725, 0.006047558, 0.02049171, 0.009058285, 0.003967988, 0.0004177118, 0.006226962, 0.006226962, 0.09977736,  0.0916344)

#Data
temp <- c(temp_pfluger_81, temp_pfluger_81_sin, temp_pfluger_80,  temp_foster_64, temp_upatham_73)
rate <- c(rate_pfluger_81, rate_pfluger_81_sin, rate_pfluger_80,  rate_foster_64, rate_upatham_73)

#The first model is quadratic 
# keep just a single curve
.x <- data.frame(temp, rate)

# fit model

fit_1 <- nls_multstart(rate~flinn_1991(temp = temp, a, b, c),
                       data = .x,
                       iter = c(5,5,5),
                       start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') - 10,
                       start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') + 10,
                       lower = get_lower_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                       upper = get_upper_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                       supp_errors = 'Y',
                       convergence_count = FALSE)


# calculate additional traits
calc_params(fit_1) %>%
  # round for easy viewing
  mutate_all(round, 2)


#P. H. Joubert, 1986, Biomphalaria pfeifferi 
temp_joubert_86_p <- c(     34,         36,       38,       40)
rate_joubert_86_p <- c(0.1457006, 0.6661235, 1.649401,  5.821733)


#El- Hassan 1974, Biomphalaria Alexandrina
temp_has_74_a_h <- c(35,    37)
rate_has_74_a_h <- 1/c(0.03891, 0.3)


# Data 
temp <- c(temp_joubert_86_p,  temp_has_74_a_h) 
rate <- c(rate_joubert_86_p,  rate_has_74_a_h)



.x <- data.frame(temp, rate)

fit_2 <- nls_multstart(rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 15),
                       data = .x,
                       iter = c(4,4,4,4),
                       start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981') - 10,
                       start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981') + 10,
                       lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981'),
                       upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981'),
                       supp_errors = 'Y',
                       convergence_count = FALSE)

# calculate additional traits
calc_params(fit_2) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data
temp <- c(temp_pfluger_81, temp_pfluger_81_sin, temp_pfluger_80,  temp_foster_64, temp_upatham_73)
rate <- c(rate_pfluger_81, rate_pfluger_81_sin, rate_pfluger_80,  rate_foster_64, rate_upatham_73)

# keep just a single curve
.x_1 <- data.frame(temp, rate)

new_data_1 <- data.frame(temp = seq(min(.x_1$temp), max(.x_1$temp)+10, 0.2))
preds_1 <- augment(fit_1, newdata = new_data_1)

# predict new data for other model
temp <- c(temp_joubert_86_p,  temp_has_74_a_h) 
rate <- c(rate_joubert_86_p,  rate_has_74_a_h)

# keep just a single curve
.x_2 <- data.frame(temp, rate)

new_data_2 <- data.frame(temp = seq(min(.x_2$temp)-4, max(.x_2$temp), 0.2))
preds_2 <- augment(fit_2, newdata = new_data_2)


preds_1 <- preds_1$.fitted[which(preds_1$temp <= 34)]
preds_2 <- preds_2$.fitted[which(preds_2$temp > 34)]

preds <- c(preds_1, preds_2)  
temperature <- seq(min(.x_1$temp), max(.x_2$temp), 0.2)



library(TeachingDemos)
pdf(file = "mortality_rate_inf_snail.pdf", width = 5, height = 5)
par(mar=c(4, 4, 6, 1), xpd = TRUE)
plot(temperature, preds, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l", cex.axis=1.3,
     xlim = c(min(temperature), max(temperature)), ylim = c(-1.2, max(preds)+1))

#text(temperature[which.min(preds)],   -1,     labels = "^", cex = 1.5, col = 1)
#text(temperature[which.min(preds)],    -0.5,     labels = round(temperature[which.min(preds)],3), cex = 1.5, col = 1)


points(temp_pfluger_81, rate_pfluger_81, col = 2, pch = 19)
points(temp_pfluger_81_sin, rate_pfluger_81_sin, col = 3, pch = 19)
points(temp_pfluger_80, rate_pfluger_80, col = 4, pch = 19)
points(temp_foster_64, rate_foster_64, col = 5, pch = 19)
points(temp_upatham_73, rate_upatham_73, col = 6, pch = 19)



legend("bottomleft", inset=c(0, 1.01), legend=c(expression(paste(italic('B. glabrata with S. mansoni,'), "  W. Pflüger, 1981")),  expression(paste(italic('B. glabrata with S. mansoni,'), "  W. Pflüger, 1981")),
                                             expression(paste(italic('B. glabrata with S. mansoni,'), "  W. Pflüger, 1980")), expression(paste(italic('B. pfeifferi with S. mansoni,'), "  R. Foster, 1964")),
                                             expression(paste(italic('B. glabrata with S. mansoni,'), "  Upatham, 1973"))), 
       pch = 19, col = c(2, 3, 4, 5, 6), title="Data",  cex = .8)

v <- c(which(temperature == 16):which(temperature == 32))
subplot(
  plot(temperature[v], preds[v],   
       lwd = 3, type = "l", cex.axis = .7, ylim = c(0, 0.05), xlim = c(16, 35)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
  
)

subplot(
  plot(temp_pfluger_81,  rate_pfluger_81, col = 2, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(16, 35)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_pfluger_81_sin,  rate_pfluger_81_sin, col = 3, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(16, 35)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)


subplot(
  plot(temp_pfluger_80,  rate_pfluger_80, col = 4, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(16, 35)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_foster_64,  rate_foster_64, col = 5, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(16, 35)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_upatham_73,  rate_upatham_73, col = 6, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(16, 35)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 2.5, cex = 1.5)
mtext(text = "Mortality rate of infected snails per day", side = 2, line = 2.5, cex = 1.5)
#main = "Death rate of snail per day", cex.main=1.7

dev.off()




pdf(file = "paper_mortality_rate_inf_snail.pdf", width = 5, height = 5)
par(mar=c(6, 5, 4, 1), xpd = TRUE)

plot(temperature, preds, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2,
     xlim = c(min(temperature), max(temperature)), ylim = c(-1, max(preds)+1), las = 1)

points(temp_pfluger_81, rate_pfluger_81, col = 2, pch = 19)
points(temp_pfluger_81_sin, rate_pfluger_81_sin, col = 3, pch = 19)
points(temp_pfluger_80, rate_pfluger_80, col = 4, pch = 19)
points(temp_foster_64, rate_foster_64, col = 5, pch = 19)
points(temp_upatham_73, rate_upatham_73, col = 6, pch = 19)


v <- c(which(temperature == 16):which(temperature == 32))
subplot(
  plot(temperature[v], preds[v],   
       lwd = 3, type = "l", cex.axis = .7, ylim = c(0, 0.05), xlim = c(16, 33), las = 1),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
  
)

subplot(
  plot(temp_pfluger_81,  rate_pfluger_81, col = 2, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(16, 33)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_pfluger_81_sin,  rate_pfluger_81_sin, col = 3, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(16, 33)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)


subplot(
  plot(temp_pfluger_80,  rate_pfluger_80, col = 4, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(16, 33)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_foster_64,  rate_foster_64, col = 5, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(16, 33)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)
subplot(
  plot(temp_upatham_73,  rate_upatham_73, col = 6, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(16, 33)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 4, cex = 2, at = 30)
mtext(text = expression(mu[i]), side = 2, line = 3, cex = 2.5, las = 1)
mtext(text = expression(paste(italic('Biomphalaria'), " spp.")), side = 3, line = 2, cex = 2, at = 30)

dev.off()