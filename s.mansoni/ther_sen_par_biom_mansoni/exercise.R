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


#p_h_joubert_1986_Biom_pfeifferi
temp_joubert_86_p <- c(     34,         36,       38,       40)
rate_joubert_86_p <- c(0.1457006, 0.6661235, 1.649401,  5.821733)


#El- Hassan 1974, Biomphalaria Alexandrina
temp_has_74_a_h <- c(35,    37)
rate_has_74_a_h <- c(0.03891, 0.3)

# Data 

temp <- c(temp_mccreesh_14, temp_applet_77,  temp_has_74_a,  temp_sturr_66, temp_foster_64, temp_eleman_82, temp_shiff_63, temp_joubert_86_p,  temp_has_74_a_h) 
rate <- c(rate_mccreesh_14, rate_applet_77,  rate_has_74_a,  rate_sturr_66, rate_foster_64, rate_eleman_82, rate_shiff_63, rate_joubert_86_p,  rate_has_74_a_h)

plot(temp, log(rate))

#The first model is quadratic 
# keep just a single curve
.x <- data.frame(temp, log(rate))

# fit model

fit <- nls_multstart(rate~spain_1982(temp = temp, a,b,c,r0),
                       data = .x,
                       iter = c(4,4,4,4),
                       start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') - 1,
                       start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') + 1,
                       
                       supp_errors = 'Y',
                       convergence_count = FALSE)


# calculate additional traits
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)


new_data <- data.frame(temp = seq(min(.x$temp), max(.x$temp), 0.5))
preds <- augment(fit, newdata = new_data)


#####################

plot(new_data$temp, preds$.fitted, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l", cex.axis=1.3,
     xlim = c(min(new_data$temp), max(new_data$temp)))

#text(temperature[which.min(preds)],   -8.5,     labels = "^", cex = 1.5, col = 1)
#text(temperature[which.min(preds)],    -7.5,     labels = round(temperature[which.min(preds)],3), cex = 1.5, col = 1)


points(temp, rate, col = 2, pch = 19)

mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 4, cex = 2, at = 25)
mtext(text = expression(log(mu)), side = 2, line = 3, cex = 2, las = 1)
mtext(text = expression(paste(italic('Biomphalaria'), " spp.")), side = 3, line = 2, cex = 2, at = 25)



#########################
