# This perform the best curve for mu
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)



#(McCreesh et al. 2014), Biomphalaria sudanica
temp_mccreesh_14 <- c(    13.4,        15.7,       16.7,        18.9,        20.9,        22.8,        26.7,        28.3,        29.5,      32.0)
rate_mccreesh_14 <- c(0.009404046, 0.01682591, 0.005224019, 0.007002545, 0.001578888, 0.008951098, 0.007300832,  0.02229403, 0.03898983, 0.05906629)


#C. C. Appleton, 1977, experiment with Biomphalaria pfeifferi
temp_applet_77 <- c(  25,             27,       29)
rate_applet_77 <- c(0.002832593, 0.01014083, 0.02164846)


#El- Hassan 1974, Biomphalaria Alexandrina
temp_has_74_a <- c(   10,      12.5,    15,      20,      25,      30)
rate_has_74_a <- c(0.01161, 0.01119, 0.01077, 0.00596, 0.00674, 0.01505)


#R. F. Sturrock, 1966, experiment with Biomphalaria pfeifferi
temp_sturr_66 <- c(  19,            25,        30)
rate_sturr_66 <- c(0.01665,      0.01348,   0.03842)


#Foster 1964, Biomphalaria pfeifleri
temp_foster_64 <- c(  22.85,        24.01,      26.26,    28.07)
rate_foster_64 <- c(0.006674626, 0.01155245, 0.0260108, 0.03678792)


#(El-Emam and Madsen 1982),  Biomphalaria alexandrina
temp_eleman_82 <- c(   18,         26,          28)
rate_eleman_82 <- c(0.006188814, 0.006188814, 0.009373412)


#(Shiff and Garnett 1963), Biomphalaria pfeifferi. 
temp_shiff_63 <- c(      18,                 22,            25,                27)
rate_shiff_63 <- c( 0.001650303,     0.002659904,      0.002373008,     0.004170385)



# Data 

temp <- c(temp_mccreesh_14, temp_applet_77,  temp_has_74_a,  temp_sturr_66, temp_foster_64, temp_eleman_82, temp_shiff_63) 
rate <- c(rate_mccreesh_14, rate_applet_77,  rate_has_74_a,  rate_sturr_66, rate_foster_64, rate_eleman_82, rate_shiff_63)


#The first model is quadratic 
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
                       supp_errors = 'Y',
                       convergence_count = FALSE)


# calculate additional traits
calc_params(fit_1) %>%
  # round for easy viewing
  mutate_all(round, 2)


#P. H. Joubert, 1986, Biomphalaria pfeifferi 
#p_h_joubert_1986_Biom_pfeifferi
temp_joubert_86_p <- c(     34,         36,       38,       40)
rate_joubert_86_p <- c(0.1457006, 0.6661235, 1.649401,  5.821733)


#El- Hassan 1974, Biomphalaria Alexandrina
temp_has_74_a_h <- c(35,    37)
rate_has_74_a_h <- c(0.03891, 0.3)

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
temp <- c(temp_mccreesh_14, temp_applet_77,  temp_has_74_a,  temp_sturr_66, temp_foster_64, temp_eleman_82, temp_shiff_63) 
rate <- c(rate_mccreesh_14, rate_applet_77,  rate_has_74_a,  rate_sturr_66, rate_foster_64, rate_eleman_82, rate_shiff_63)

# keep just a single curve
.x_1 <- data.frame(temp, rate)

new_data_1 <- data.frame(temp = seq(min(.x_1$temp), max(.x_1$temp)+10, 0.2))
preds_1 <- augment(fit_1, newdata = new_data_1)

# predict new data for other model
temp <- c(temp_joubert_86_p,  temp_has_74_a_h) 
rate <- c(rate_joubert_86_p,  rate_has_74_a_h)

# keep just a single curve
.x_2 <- data.frame(temp, rate)

new_data_2 <- data.frame(temp = seq(min(.x_2$temp)-5, max(.x_2$temp), 0.2))
preds_2 <- augment(fit_2, newdata = new_data_2)


preds_1 <- preds_1$.fitted[which(preds_1$temp <= 35.5)]
preds_2 <- preds_2$.fitted[which(preds_2$temp > 35.5)]

preds <- c(preds_1, preds_2)  
temperature <- seq(min(.x_1$temp), max(.x_2$temp), 0.2)


#####################
pdf(file = "log_mortality_rate_snail.pdf", width = 5, height = 5)
par(mar=c(6, 8, 4, 1), xpd = TRUE)
plot(temperature, log(preds), col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l", cex.axis=1.3,
     xlim = c(min(temperature), max(temperature)), ylim = c(-8.5, 4))

text(temperature[which.min(preds)],   -8.5,     labels = "^", cex = 1.5, col = 1)
text(temperature[which.min(preds)],    -7.5,     labels = round(temperature[which.min(preds)],3), cex = 1.5, col = 1)


points(c(temp_ches_17, temp_ches_17_2), c(rate_ches_17,rate_ches_17_2), col = 2, pch = 19)
points(c(temp_has_74_t, temp_has_74_t_2), c(rate_has_74_t, rate_has_74_t_2), col = 3, pch = 19)
points(temp_joubert_86_a, rate_joubert_86_a, col = 4, pch = 19)
points(temp_joubert_86_g, rate_joubert_86_g, col = 5, pch = 19)
points(temp_eleman_82, rate_eleman_82, col = 6, pch = 19)
points(temp_shiff_64, rate_shiff_64, col = 7, pch = 19)
points(temp_kubirizajournal_10, rate_kubirizajournal_10, col = 8, pch = 19)



mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 4, cex = 2, at = 25)
mtext(text = expression(log(mu)), side = 2, line = 3, cex = 2, las = 1)
mtext(text = expression(paste(italic('Biomphalaria'), " spp.")), side = 3, line = 2, cex = 2, at = 25)

dev.off()


#########################

library(TeachingDemos)
#pdf(file = "mortality_rate_snail.pdf", width = 5, height = 5)
par(mar=c(4, 4, 8, 1), xpd = TRUE)
plot(temperature, preds, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l", cex.axis=1.3,
     xlim = c(min(temperature), max(temperature)), ylim = c(-1, max(preds)+1))

text(temperature[which.min(preds)],   -0.8,     labels = "^", cex = 1.5, col = 1)
text(temperature[which.min(preds)],    -0.4,     labels = round(temperature[which.min(preds)],3), cex = 1.5, col = 1)


points(temp_applet_77, rate_applet_77, col = 2, pch = 19)
points(c(temp_has_74_a, temp_has_74_a_h), c(rate_has_74_a, temp_has_74_a_h), col = 3, pch = 19)
points(temp_sturr_66, rate_sturr_66, col = 4, pch = 19)
points(temp_foster_64, rate_foster_64, col = 5, pch = 19)
points(temp_joubert_86_p, rate_joubert_86_p, col = 6, pch = 19)
points(temp_mccreesh_14, rate_mccreesh_14, col = 7, pch = 19)
points(temp_eleman_82, rate_eleman_82, col = 8, pch = 19)
points(temp_shiff_63, rate_shiff_63, col = 9, pch = 19)





legend("bottomleft", inset=c(0, 1.01), legend=c(expression(paste(italic('B. pfeifferi,'), "  C. C. APPLETON, 1977")), expression(paste(italic('B. alexandrina,'), "  A. A. EL-HASSAN, 1974")),
                                             expression(paste(italic('B. pfeifferi,'), "  R. F. Sturrock, 1966")), expression(paste(italic('B. pfeifferi,'), "  R. Foster, 1964")), 
                                             expression(paste(italic('B. pfeifferi,'), "  Shiff and Garnett 1963")),
                                             expression(paste(italic('B. sudanica,'), "  McCreesh et al., 2014")),
                                             expression(paste(italic('B. alexandrina,'), "  El-Emam and Madsen, 1982"))),                                           expression(paste(italic('B. pfeifferi'), "  Shiff and Garnett,  1963")),
       pch = 19, col = c(2,  3,  4,  5, 6, 7, 8, 9), title="Data", cex = .8) 


v <- c(which(temperature == 10):which(temperature == 30))
subplot(
  plot(temperature[v], preds[v],   
       lwd = 3, type = "l", cex.axis = .7, ylim = c(0, 0.035), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
  
)

subplot(
  plot(temp_applet_77,  rate_applet_77, col = 2, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_has_74_a,  rate_has_74_a, col = 3, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)


subplot(
  plot(temp_sturr_66,  rate_sturr_66, col = 4, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_foster_64,  rate_foster_64, col = 5, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_mccreesh_14,  rate_mccreesh_14, col = 6, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)


subplot(
  plot(temp_eleman_82,  rate_eleman_82, col = 7, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_shiff_63,  rate_shiff_63, col = 8, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)


mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 2.5, cex = 1.5)
mtext(text = "Mortality rate of snails per day", side = 2, line = 2.5, cex = 1.5)
#main = "Death rate of snail per day", cex.main=1.7

#dev.off()



#pdf(file = "paper_mortality_rate_snail.pdf", width = 5, height = 5)
par(mar=c(6, 5, 4, 1), xpd = TRUE)
plot(temperature, preds, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 2, type = "l", cex.axis=2,
     xlim = c(min(temperature), max(temperature)), ylim = c(-1, max(preds)+1), las =1)

text(temperature[which.min(preds)],   -0.8,     labels = "^", cex = 1.5, col = 1)
text(temperature[which.min(preds)],    -0.4,     labels = round(temperature[which.min(preds)],3), cex = 1.5, col = 1)


points(temp_applet_77, rate_applet_77, col = 2, pch = 19)
points(c(temp_has_74_a, temp_has_74_a_h), c(rate_has_74_a, temp_has_74_a_h), col = 3, pch = 19)
points(temp_sturr_66, rate_sturr_66, col = 4, pch = 19)
points(temp_foster_64, rate_foster_64, col = 5, pch = 19)
points(temp_joubert_86_p, rate_joubert_86_p, col = 6, pch = 19)
points(temp_mccreesh_14, rate_mccreesh_14, col = 7, pch = 19)
points(temp_eleman_82, rate_eleman_82, col = 8, pch = 19)
points(temp_shiff_63, rate_shiff_63, col = 9, pch = 19)


v <- c(which(temperature == 10):which(temperature == 30))
subplot(
  plot(temperature[v], preds[v],   
       lwd = 3, type = "l", cex.axis = .7, ylim = c(0, 0.035), xlim = c(9, 32), las =1),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
  
)

subplot(
  plot(temp_applet_77,  rate_applet_77, col = 2, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_has_74_a,  rate_has_74_a, col = 3, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)


subplot(
  plot(temp_sturr_66,  rate_sturr_66, col = 4, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_foster_64,  rate_foster_64, col = 5, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_mccreesh_14,  rate_mccreesh_14, col = 6, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)


subplot(
  plot(temp_eleman_82,  rate_eleman_82, col = 7, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)

subplot(
  plot(temp_shiff_63,  rate_shiff_63, col = 8, pch = 19, axes = FALSE, ylim = c(0, 0.03), xlim = c(9, 32)),
  x  = grconvertX(c(.25,1), from='npc'),
  y = grconvertY(c(.25, 1), from='npc'),
  type='fig', pars=list(mar=c(3,0,.5,3)+.2),
)


mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 4, cex = 2, at = 31)
mtext(text = expression(mu), side = 2, line = 3, cex = 2.5, las = 1)
mtext(text = expression(paste(italic('Biomphalaria'), " spp.")), side = 3, line = 2, cex = 2, at = 31)


#dev.off()