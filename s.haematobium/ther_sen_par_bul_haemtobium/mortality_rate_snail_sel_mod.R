# This perform the best curve for mu
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)

# This analysis are up to 31 degree temperature. 

# chester_kalinda_2017_Bulinus_globosus 
temp_ches_17 <- c(15.5,           21.2,        25.8,      31)
rate_ches_17 <- c(0.00244925, 0.006901087, 0.01151302, 0.02073716)


#El- Hassan 1974, Bulinus truncatus
temp_has_74_t <- c(   10,        15,     20,       25,      30)
rate_has_74_t <- c(0.00442,   0.00218, 0.00144, 0.00518, 0.03077)

#(El-Emam and Madsen 1982), Bulinus truncatus
temp_eleman_82 <- c(    10,         18,          26,        28 )
rate_eleman_82 <- c(0.00256859, 0.004560943, 0.003184598, 0.004560943)

#(Shiff 1964), Bulinus globosus
temp_shiff_64 <- c(   18,           22,            25,        27)
rate_shiff_64 <- c(0.005314696, 0.006322683, 0.01045628, 0.01057945)

#(Kubirizajournal et al. 2010),  Bulinus nyassanus
temp_kubirizajournal_10 <- c(22,              25,          28,          31)
rate_kubirizajournal_10 <- c(0.006849573,  0.00513718, 0.02054872, 0.01027436)



temp <- c(temp_ches_17, temp_has_74_t,  temp_eleman_82, temp_shiff_64,  temp_kubirizajournal_10)
rate <- c(rate_ches_17, rate_has_74_t,  rate_eleman_82, rate_shiff_64,  rate_kubirizajournal_10)



#The first model is quadratic 
# keep just a single curve
.x <- data.frame(temp, rate)

# fit model

fit_1 <- nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                     data = .x,
                     iter = c(4,4,4),
                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 0.5,
                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 0.5,
                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                     supp_errors = 'Y',
                     convergence_count = FALSE)


# calculate additional traits
calc_params(fit_1) %>%
  # round for easy viewing
  mutate_all(round, 2)



# chester_kalinda_2017_Bulinus_globosus 
temp_ches_17_2 <-      35.5
rate_ches_17_2 <-  0.04642702


#p_h_joubert_1986_Bul_globosus
temp_joubert_86_g <- c(     34,         36,       38,       40)
rate_joubert_86_g <- c(0.02924029, 0.1490773, 0.8022027, 3.540546)


#p_h_joubert_1986_Bul_africanus
temp_joubert_86_a <- c(     34,         36,       38,       40)
rate_joubert_86_a <- c(0.1044946, 0.5664857, 1.618022,  4.344277)


temp_has_74_t_2 <-     35
rate_has_74_t_2 <-  0.09622


temp <- c(temp_ches_17_2,temp_joubert_86_g, temp_joubert_86_a, temp_has_74_t_2)
rate <- c(rate_ches_17_2, rate_joubert_86_g, rate_joubert_86_a, rate_has_74_t_2)

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

temp <- c(temp_ches_17, temp_has_74_t,  temp_eleman_82, temp_shiff_64,  temp_kubirizajournal_10)
rate <- c(rate_ches_17, rate_has_74_t,  rate_eleman_82, rate_shiff_64,  rate_kubirizajournal_10)


# keep just a single curve
.x_1 <- data.frame(temp, rate)

new_data_1 <- data.frame(temp = seq(min(.x_1$temp), max(.x_1$temp)+10, 0.2))
preds_1 <- augment(fit_1, newdata = new_data_1)


temp <- c(temp_ches_17_2,temp_joubert_86_g, temp_joubert_86_a, temp_has_74_t_2)
rate <- c(rate_ches_17_2, rate_joubert_86_g, rate_joubert_86_a, rate_has_74_t_2)


# keep just a single curve
.x_2 <- data.frame(temp, rate)

new_data_2 <- data.frame(temp = seq(min(.x_2$temp)-5, max(.x_2$temp), 0.2))
preds_2 <- augment(fit_2, newdata = new_data_2)


preds_1 <- preds_1$.fitted[which(preds_1$temp <= 33)]
preds_2 <- preds_2$.fitted[which(preds_2$temp > 33)]

preds <- c(preds_1, preds_2)  
temperature <- seq(min(.x_1$temp), max(.x_2$temp), 0.2)



######################

#####################
pdf(file = "log_mortality_rate_snail.pdf", width = 5, height = 5)
par(mar=c(6, 8, 4, 1), xpd = TRUE)
plot(temperature, log(preds), col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l", cex.axis=1.3,
     xlim = c(min(temperature), max(temperature)), ylim = c(-8.5, 4))

text(temperature[which.min(preds)],   -8.5,     labels = "^", cex = 1.5, col = 1)
text(temperature[which.min(preds)],    -7.5,     labels = round(temperature[which.min(preds)],3), cex = 1.5, col = 1)

points(c(temp_ches_17, temp_ches_17_2), log(c(rate_ches_17,rate_ches_17_2)), col = 2, pch = 19)
points(c(temp_has_74_t, temp_has_74_t_2), log(c(rate_has_74_t, rate_has_74_t_2)), col = 3, pch = 19)
points(temp_joubert_86_a, log(rate_joubert_86_a), col = 4, pch = 19)
points(temp_joubert_86_g, log(rate_joubert_86_g), col = 5, pch = 19)
points(temp_eleman_82, log(rate_eleman_82), col = 6, pch = 19)
points(temp_shiff_64, log(rate_shiff_64), col = 7, pch = 19)
points(temp_kubirizajournal_10, log(rate_kubirizajournal_10), col = 8, pch = 19)



mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 4, cex = 2, at = 25)
mtext(text = expression(log(mu)), side = 2, line = 3, cex = 2, las = 1)
mtext(text = expression(paste(italic('Bulinus'), " spp.")), side = 3, line = 2, cex = 2, at = 25)

dev.off()
#######################






library(TeachingDemos)
pdf(file = "mortality_rate_snail.pdf", width = 5, height = 5)
par(mar=c(4, 4, 8, 1), xpd = TRUE)
plot(temperature, preds, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l", cex.axis=1.3,
     xlim = c(min(temperature), max(temperature)), ylim = c(-.6, max(preds)+1))

text(temperature[which.min(preds)],   -0.6,     labels = "^", cex = 1.5, col = 1)
text(temperature[which.min(preds)],    -0.3,     labels = round(temperature[which.min(preds)],3), cex = 1.5, col = 1)

points(c(temp_ches_17, temp_ches_17_2), c(rate_ches_17,rate_ches_17_2), col = 2, pch = 19)
points(c(temp_has_74_t, rate_has_74_t_2), c(rate_has_74_t, rate_has_74_t_2), col = 3, pch = 19)
points(temp_joubert_86_a, rate_joubert_86_a, col = 4, pch = 19)
points(temp_joubert_86_g, rate_joubert_86_g, col = 5, pch = 19)
points(temp_eleman_82, rate_eleman_82, col = 6, pch = 19)
points(temp_shiff_64, rate_shiff_64, col = 7, pch = 19)
points(temp_kubirizajournal_10, rate_kubirizajournal_10, col = 8, pch = 19)




legend("bottomleft", inset=c(0, 1.01), legend=c(expression(paste(italic('B. globosus,'), "  Kalinda et al., 2017")), expression(paste(italic('B. truncatus,'), "  A. A. EL- HASSAN, 1974")),
                                             expression(paste(italic('B. africanus '), " Joubert et al., 1986")), expression(paste(italic('B. globosus '), " Joubert et al., 1986")),
                                             expression(paste(italic('B. truncatus '), " El-Emam and Madsen, 1982")),  expression(paste(italic('B. globosus '), " Shiff,  1964")),
                                             expression(paste(italic('B. nyassanus '), " Kubirizajournal et al., 2010"))), 
       pch = 19, col = c(2, 3, 4, 5, 6, 7, 8), title="Data", cex = .8)
 
 
v <- c(which(temperature == 10):which(temperature == 32))
subplot(
  plot(temperature[v], preds[v],   
       lwd = 3, type = "l", cex.axis = .7, ylim = c(0, 0.035), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
  
)

subplot(
  plot(temp_ches_17,  rate_ches_17, col = 2, pch = 19, axes = FALSE, ylim = c(0, 0.035), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_has_74_t,  rate_has_74_t, col = 3, pch = 19, axes = FALSE, ylim = c(0, 0.035), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_eleman_82,  rate_eleman_82, col = 6, pch = 19, axes = FALSE, ylim = c(0, 0.035), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_shiff_64,  rate_shiff_64, col = 7, pch = 19, axes = FALSE, ylim = c(0, 0.035), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_kubirizajournal_10,  rate_kubirizajournal_10, col = 8, pch = 19, axes = FALSE, ylim = c(0, 0.035), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)


mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 2.5, cex = 1.5)
mtext(text = "Mortality rate of snails per day", side = 2, line = 2.5, cex = 1.5)
#main = "Death rate of snail per day", cex.main=1.7

dev.off()



pdf(file = "paper_mortality_rate_snail.pdf", width = 5, height = 5)
par(mar=c(6, 6, 4, 1), xpd = TRUE)
plot(temperature, preds, col = 1, las = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2,
     xlim = c(min(temperature), max(temperature)), ylim = c(-.6, max(preds)+1))

text(temperature[which.min(preds)],   -0.6,     labels = "^", cex = 1.5, col = 1)
text(temperature[which.min(preds)],    -0.3,     labels = round(temperature[which.min(preds)],3), cex = 1.5, col = 1)


points(c(temp_ches_17, temp_ches_17_2), c(rate_ches_17,rate_ches_17_2), col = 2, pch = 19)
points(c(temp_has_74_t, temp_has_74_t_2), c(rate_has_74_t, rate_has_74_t_2), col = 3, pch = 19)
points(temp_joubert_86_a, rate_joubert_86_a, col = 4, pch = 19)
points(temp_joubert_86_g, rate_joubert_86_g, col = 5, pch = 19)
points(temp_eleman_82, rate_eleman_82, col = 6, pch = 19)
points(temp_shiff_64, rate_shiff_64, col = 7, pch = 19)
points(temp_kubirizajournal_10, rate_kubirizajournal_10, col = 8, pch = 19)



v <- c(which(temperature == 10):which(temperature == 32))
subplot(
  plot(temperature[v], preds[v],   
       lwd = 3, type = "l", cex.axis = .7, ylim = c(0, 0.03), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
  
)

subplot(
  plot(temp_ches_17,  rate_ches_17, col = 2, pch = 19, axes = FALSE, ylim = c(0, 0.035), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_has_74_t,  rate_has_74_t, col = 3, pch = 19, axes = FALSE, ylim = c(0, 0.035), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_eleman_82,  rate_eleman_82, col = 6, pch = 19, axes = FALSE, ylim = c(0, 0.035), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_shiff_64,  rate_shiff_64, col = 7, pch = 19, axes = FALSE, ylim = c(0, 0.035), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_kubirizajournal_10,  rate_kubirizajournal_10, col = 8, pch = 19, axes = FALSE, ylim = c(0, 0.035), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)


mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 4, cex = 2)
mtext(text = expression(mu), side = 2, line = 4, cex = 2.5, las = 1)
mtext(text = expression(paste(italic('Bulinus'), " spp.")), side = 3, line = 2, cex = 2)

dev.off()