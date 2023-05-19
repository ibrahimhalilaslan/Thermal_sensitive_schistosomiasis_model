#This script perform multiple different thermal sensitive models to find the best fit 
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)

# write function to label ggplot2 panels
label_facets_num <- function(string){
  len <- length(string)
  string = paste('(', 1:len, ') ', string, sep = '')
  return(string)
}

# McCreesh et al. 2014, Biomphalaria sudanica.
temp_mccreesh_14 <- c(13.4, 15.7,  16.7,    18.9,    20.9,    22.8,   26.7,  28.3, 29.5, 32.0)
rate_mccreesh_14 <- c(  0,    0,  0.0349,  0.0369,  0.0395,  0.0386, 0.0499,  0,     0,   0)



#C. C. Appleton, 1977, experiment with Biomphalaria pfeifferi
temp_applet_77 <- c(25,       27,      29)
rate_applet_77 <- c(0.07058,  0.06824,    0)


#El- Hassan 1974, Biomphalaria Alexandrina
temp_has_74_a <- c(15,        20,      25,   30,    35)
rate_has_74_a <- c(0.00683, 0.0229, 0.0309, 0.00824, 0)



#R. F. Sturrock, 1966, experiment with Biomphalaria pfeifferi
temp_sturr_66 <- c(19,        25,        30,    35)
rate_sturr_66 <- c(0.02720, 0.07378,    0.06568, 0)

# Michelson 1961, Biomphalaria glabrata.
temp_michel_61 <- c(5, 15,   20,     25,     30,  35)
rate_michel_61 <- c(0, 0, 0.0132, 0.0349, 0.00969, 0)


# Shiff and Garnett 1963, Biomphalaria pfeifferi
temp_shiff_63 <- c(   18,     22,   25,     27)
rate_shiff_63 <- c(0.0138, 0.0507, 0.0349, 0.036)



temp <- c(temp_mccreesh_14, temp_applet_77, temp_has_74_a, temp_sturr_66, temp_michel_61, temp_shiff_63)
rate <- c(rate_mccreesh_14, rate_applet_77, rate_has_74_a, rate_sturr_66, rate_michel_61, rate_shiff_63)




# keep just a single curve
d <- data.frame(temp, rate)

# fit model

fit <- nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                     data = d,
                     iter = c(4,4,4),
                     start_lower = get_start_vals(d$temp, d$rate, model_name = 'gaussian_1987') - 10,
                     start_upper = get_start_vals(d$temp, d$rate, model_name = 'gaussian_1987') + 10,
                     lower = get_lower_lims(d$temp, d$rate, model_name = 'gaussian_1987'),
                     upper = get_upper_lims(d$temp, d$rate, model_name = 'gaussian_1987'),
                     supp_errors = 'Y',
                     convergence_count = FALSE)
fit


# calculate additional traits
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

# predict new data

new_data <- data.frame(temp = seq(min(d$temp), max(d$temp), 0.1))
preds <- augment(fit, newdata = new_data)


pdf(file = "fecundity_rate_snail.pdf", width = 5, height = 5)
par(mar=c(4, 4, 7, 1), xpd = TRUE)
plot(preds$temp, preds$.fitted, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l",  cex.axis=1.3,
     xlim = c(min(preds$temp), max(preds$temp)), ylim = c(0, 0.08))

points(temp_applet_77, rate_applet_77, col = 2, pch = 19)
points(temp_has_74_a, rate_has_74_a, col = 3, pch = 19)
points(temp_sturr_66, rate_sturr_66, col = 4, pch = 19)
points(temp_mccreesh_14, rate_mccreesh_14, col = 5, pch = 19)
points(temp_michel_61, rate_michel_61, col = 6, pch = 19)
points(temp_shiff_63, rate_shiff_63, col = 7, pch = 19)


legend("bottomleft", inset=c(0, 1.01), legend=c(expression(paste(italic('B. pfeifferi,'), "  C. C. APPLETON, 1977")),
                                              expression(paste(italic('B. alexandrina,'), "  A. A. EL-HASSAN, 1974")), expression(paste(italic('B. pfeifferi,'), "  R. F. Sturrock, 1966")),
                                              expression(paste(italic('B. sudanica,'), "  McCreesh et al., 2014")),  expression(paste(italic('B. glabrata,'), "  Michelson, 1961")),
                                               expression(paste(italic('B. pfeifferi,'), "  Shiff and Garnett 1963"))), 
       pch = 19, col = c(2, 3, 4, 5, 6, 7), title = "Data", cex = .8)    

text(preds$temp[which.max(preds$.fitted)],   0,     labels = "^", cex = 1.5, col = 1)
text(preds$temp[which.max(preds$.fitted)],  0.006,     labels = round(preds$temp[which.max(preds$.fitted)],3), cex = 1.5, col = 1)

mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 2.5, cex = 1.5)
mtext(text = "Fecundity rate of snails per day", side = 2, line = 2.5, cex = 1.5)

#main = "The fecundity rate of snail per day", cex.main=1.7
dev.off()



pdf(file = "paper_fecundity_rate_snail.pdf", width = 5, height = 5)
par(mar=c(4, 4, 1, 1), xpd = TRUE)
plot(preds$temp, preds$.fitted, col = 1,  pch = 21, lwd = 3, ylab = "",  xlab = "", cex.lab = 1.3, type = "l",  cex.axis=1.3,
     xlim = c(min(preds$temp), max(preds$temp)), ylim = c(0, 0.08))

points(temp_applet_77, rate_applet_77, col = 2, pch = 19)
points(temp_has_74_a, rate_has_74_a, col = 3, pch = 19)
points(temp_sturr_66, rate_sturr_66, col = 4, pch = 19)
points(temp_mccreesh_14, rate_mccreesh_14, col = 5, pch = 19)
points(temp_michel_61, rate_michel_61, col = 6, pch = 19)
points(temp_shiff_63, rate_shiff_63, col = 7, pch = 19)


text(preds$temp[which.max(preds$.fitted)],   0,     labels = "^", cex = 1.5, col = 1)
text(preds$temp[which.max(preds$.fitted)],  0.005,     labels = round(preds$temp[which.max(preds$.fitted)],3), cex = 1.5, col = 1)

mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 2.5, cex = 1.5)
mtext(text = "Fecundity rate of snails per day", side = 2, line = 2.5, cex = 1.5)

#main = "The fecundity rate of snail per day", cex.main=1.7
dev.off()

