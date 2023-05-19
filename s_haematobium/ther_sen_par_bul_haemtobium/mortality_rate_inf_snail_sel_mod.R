#The script perform the best curve for mu_i
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)



# Chester Kalinda 2017, Bulinus globosu with S. haematobium
temp_chester_17 <- c(15.5,            21.2,        25.8,           31,         36)
rate_chester_17 <- c(0.00244925,  0.006901087,   0.01151302,   0.02073716,  0.04642702)

# W. Pflüger 1984, Bulinus truncatus with S. haematobium
temp_pfluger_84 <- c(18,              19,         20,          21,                 22,      23,             25,            28,        30,          31,    32)
rate_pfluger_84 <- c(0.0131776,  0.01141752,     0.01963277, 0.02130685,     0.02204905, 0.02636316,    0.02185846,   0.01590461, 0.01570362,  0.02473268, 0.0262486)

#(Chu et al. 1966) Bulinus trancatus with S. haematobium  
temp_chu_66 <- c(   9.967,          13,          16.067,         19,         22.2,         25.067,          28,        31.133,        34.2,         37.067,      39.967)
rate_chu_66 <- c(0.0058551278, 0.0085350331, 0.0061192732, 0.0195844770, 0.0085350331, 0.0032587267, 0.0004931228, 0.0061192732, 0.0058551278, 0.0968812480, 0.0915518574)


temp <- c(temp_chester_17, temp_pfluger_84, temp_chu_66)    
rate <- c(rate_chester_17, rate_pfluger_84, rate_chu_66)



# keep just a single curve
.x <- data.frame(temp, rate)

# fit model

fit_1 <- nls_multstart(rate~spain_1982(temp = temp, a,b,c,r0),
                     data = .x,
                     iter = c(4,4,4,4),
                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') - 1,
                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') + 1,
                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                     supp_errors = 'Y')

fit_1

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


temp <- c(temp_chester_17, temp_pfluger_84, temp_chu_66)    
rate <- c(rate_chester_17, rate_pfluger_84, rate_chu_66)



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

preds_1 <- preds_1$.fitted[which(preds_1$temp <= 33.5)]
preds_2 <- preds_2$.fitted[which(preds_2$temp > 33.5)]

preds <- c(preds_1, preds_2)  
temperature <- seq(min(.x_1$temp), max(.x_2$temp), 0.2)





library(TeachingDemos)
pdf(file = "mortality_rate_inf_snail.pdf", width = 5, height = 5)
par(mar=c(4, 4, 5, 1), xpd = TRUE)
plot(temperature, preds, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l", cex.axis=1.3,
     xlim = c(min(temperature), max(temperature)), ylim = c(-.6, max(preds)+1))

#text(temperature[which.min(preds)],   -0.6,     labels = "^", cex = 1.5, col = 1)
#text(temperature[which.min(preds)],    -0.3,     labels = round(temperature[which.min(preds)],3), cex = 1.5, col = 1)

points(temp_chester_17, rate_chester_17, col = 2, pch = 19)
points(temp_pfluger_84, rate_pfluger_84, col = 3, pch = 19)
points(temp_chu_66, rate_chu_66, col = 4, pch = 19)


legend("bottomleft", inset=c(0, 1.01), legend=c(expression(paste(italic('B. globosus with S. haematobium,'), "  Kalinda et al., 2017")), 
                                             expression(paste(italic('B. truncatus with S. haematobium,'), "  W. Pflüger, 1984")),
                                             expression(paste(italic('B. truncatus with S. haematobium,'), "  Chu et al., 1966"))),
                                             pch = 19, col = c(2, 3, 4), title="Data", cex = .8)

 

v <- c(which(34 >temperature & temperature >= 15))
subplot(
  plot(temperature[v], preds[v],   
       lwd = 3, type = "l", cex.axis = .7,  ylim = c(0, 0.05), xlim = c(15, 36)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
  
)

subplot(
  plot(temp_chester_17,  rate_chester_17, col = 2, pch = 19, axes = FALSE,  ylim = c(0, 0.05), xlim = c(15, 36)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_pfluger_84,  rate_pfluger_84, col = 3, pch = 19, axes = FALSE,  ylim = c(0, 0.05), xlim = c(15, 36)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)


subplot(
  plot(temp_chu_66,  rate_chu_66, col = 4, pch = 19, axes = FALSE,  ylim = c(0, 0.05), xlim = c(15, 36)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)



mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 2.5, cex = 1.5)
mtext(text = "Mortality rate of infected snails per day", side = 2, line = 2.5, cex = 1.5)
#main = "Death rate of snail per day", cex.main=1.7

dev.off()


library(TeachingDemos)
pdf(file = "paper_mortality_rate_inf_snail.pdf", width = 5, height = 5)
par(mar=c(4, 4, 1, 1), xpd = TRUE)
plot(temperature, preds, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l", cex.axis=1.3,
     xlim = c(min(temperature), max(temperature)), ylim = c(-.5, max(preds)+1))

#text(temperature[which.min(preds)],   -0.6,     labels = "^", cex = 1.5, col = 1)
#text(temperature[which.min(preds)],    -0.3,     labels = round(temperature[which.min(preds)],3), cex = 1.5, col = 1)

points(temp_chester_17, rate_chester_17, col = 2, pch = 19)
points(temp_pfluger_84, rate_pfluger_84, col = 3, pch = 19)
points(temp_chu_66, rate_chu_66, col = 4, pch = 19)


v <- c(which(34 >temperature & temperature >= 15))
subplot(
  plot(temperature[v], preds[v],   
       lwd = 3, type = "l", cex.axis = .7, ylim = c(0, 0.05), xlim = c(15, 36)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
  
)

subplot(
  plot(temp_chester_17,  rate_chester_17, col = 2, pch = 19, axes = FALSE,  ylim = c(0, 0.05), xlim = c(15, 36)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_pfluger_84,  rate_pfluger_84, col = 3, pch = 19, axes = FALSE,  ylim = c(0, 0.05), xlim = c(15, 36)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_chu_66,  rate_chu_66, col = 4, pch = 19, axes = FALSE,  ylim = c(0, 0.05), xlim = c(15, 36)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)



mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 2.5, cex = 1.5)
mtext(text = "Mortality rate of infected snails per day", side = 2, line = 2.5, cex = 1.5)
#main = "Death rate of snail per day", cex.main=1.7

dev.off()

