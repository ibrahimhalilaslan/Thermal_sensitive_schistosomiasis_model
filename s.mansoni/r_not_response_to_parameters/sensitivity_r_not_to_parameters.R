#load the parameters 
source("par_set_mansoni.R")
########################## CALCULATE R_NOT ##################


min_temp <- 12
max_temp <- 37

# Generate a sequence of temperature    
temperature <- data.frame(temp = seq(min_temp, max_temp, by = 0.1))

# combine two piece wise function of mu_i
preds_1 <- fn_mu_i_1(temperature)$.fitted[which(temperature$temp <= 37)]
preds_2 <- fn_mu_i_2(temperature)$.fitted[which(temperature$temp > 37)]

# This is mu_i function 
preds_mu_i <- c(preds_1, preds_2)  


# This is analytic representation of R not value without control.
r_not <- (abs(lambda * fn_beta_s(temperature)$.fitted * h * fn_delta_e(temperature)$.fitted * nu_e 
       * fn_beta_h(temperature)$.fitted * fn_nu_c(temperature)$.fitted * fn_sigma_s(temperature)$.fitted/
         (fn_mu_m(temperature)$.fitted * (mu_h + mu_p) * fn_mu_c(temperature)$.fitted * preds_mu_i * 
            (fn_sigma_s(temperature)$.fitted + preds_mu_i))))^(1/2)  
  
# Assume the mu_i is constant 
r_not_mu_i <- (abs(lambda * fn_beta_s(temperature)$.fitted * h * fn_delta_e(temperature)$.fitted * nu_e 
                   * fn_beta_h(temperature)$.fitted * fn_nu_c(temperature)$.fitted * fn_sigma_s(temperature)$.fitted/
                     (fn_mu_m(temperature)$.fitted * (mu_h + mu_p) * fn_mu_c(temperature)$.fitted * mean(preds_mu_i) * 
                        (fn_sigma_s(temperature)$.fitted + mean(preds_mu_i)))))^(1/2)  

# Assume the mu_m is constant 
r_not_mu_m <- (abs(lambda * fn_beta_s(temperature)$.fitted * h * fn_delta_e(temperature)$.fitted * nu_e 
              * fn_beta_h(temperature)$.fitted * fn_nu_c(temperature)$.fitted * fn_sigma_s(temperature)$.fitted/
                (mean(fn_mu_m(temperature)$.fitted) * (mu_h + mu_p) * fn_mu_c(temperature)$.fitted * preds_mu_i * 
                   (fn_sigma_s(temperature)$.fitted + preds_mu_i))))^(1/2)  

# Assume the mu_c is constant 
r_not_mu_c <- (abs(lambda * fn_beta_s(temperature)$.fitted * h * fn_delta_e(temperature)$.fitted * nu_e 
                   * fn_beta_h(temperature)$.fitted * fn_nu_c(temperature)$.fitted * fn_sigma_s(temperature)$.fitted/
                     (fn_mu_m(temperature)$.fitted * (mu_h + mu_p) * mean(fn_mu_c(temperature)$.fitted) * preds_mu_i * 
                        (fn_sigma_s(temperature)$.fitted + preds_mu_i))))^(1/2)  

# Assume the delta_e is constant 
r_not_delta_e <- (abs(lambda * fn_beta_s(temperature)$.fitted * h * mean(fn_delta_e(temperature)$.fitted) * nu_e 
              * fn_beta_h(temperature)$.fitted * fn_nu_c(temperature)$.fitted * fn_sigma_s(temperature)$.fitted/
                (fn_mu_m(temperature)$.fitted * (mu_h + mu_p) * fn_mu_c(temperature)$.fitted * preds_mu_i * 
                   (fn_sigma_s(temperature)$.fitted + preds_mu_i))))^(1/2)  


# Assume the sigma_s is constant 
r_not_sigma_s <- (abs(lambda * fn_beta_s(temperature)$.fitted * h * fn_delta_e(temperature)$.fitted * nu_e 
              * fn_beta_h(temperature)$.fitted * fn_nu_c(temperature)$.fitted * mean(fn_sigma_s(temperature)$.fitted)/
                (fn_mu_m(temperature)$.fitted * (mu_h + mu_p) * fn_mu_c(temperature)$.fitted * preds_mu_i * 
                   (mean(fn_sigma_s(temperature)$.fitted) + preds_mu_i))))^(1/2)  

# Assume the nu_c is constant 
r_not_nu_c <- (abs(lambda * fn_beta_s(temperature)$.fitted * h * fn_delta_e(temperature)$.fitted * nu_e 
              * fn_beta_h(temperature)$.fitted * mean(fn_nu_c(temperature)$.fitted) * fn_sigma_s(temperature)$.fitted/
                (fn_mu_m(temperature)$.fitted * (mu_h + mu_p) * fn_mu_c(temperature)$.fitted * preds_mu_i * 
                   (fn_sigma_s(temperature)$.fitted + preds_mu_i))))^(1/2)  


# This is analytic representation of R not value without control.
r_not_beta_s <- (abs(lambda * mean(fn_beta_s(temperature)$.fitted) * h * fn_delta_e(temperature)$.fitted * nu_e 
              * fn_beta_h(temperature)$.fitted * fn_nu_c(temperature)$.fitted * fn_sigma_s(temperature)$.fitted/
                (fn_mu_m(temperature)$.fitted * (mu_h + mu_p) * fn_mu_c(temperature)$.fitted * preds_mu_i * 
                   (fn_sigma_s(temperature)$.fitted + preds_mu_i))))^(1/2)  


# This is analytic representation of R not value without control.
r_not_beta_h <- (abs(lambda * fn_beta_s(temperature)$.fitted * h * fn_delta_e(temperature)$.fitted * nu_e 
              * mean(fn_beta_h(temperature)$.fitted) * fn_nu_c(temperature)$.fitted * fn_sigma_s(temperature)$.fitted/
                (fn_mu_m(temperature)$.fitted * (mu_h + mu_p) * fn_mu_c(temperature)$.fitted * preds_mu_i * 
                   (fn_sigma_s(temperature)$.fitted + preds_mu_i))))^(1/2)  

# Convert the 'parameters' vector to character
parameters <- c(expression(mu[i]), expression(mu[m]), 
                expression(mu[c]), expression(delta[e]), expression(sigma[s]), 
                expression(nu[c]), expression(beta[s]), expression(beta[h]), "Overall")


df_1 <- data.frame(x = temperature$temp, 
                 y1 = r_not_mu_i/max(r_not_mu_i),
                 y2 = r_not_mu_m/max(r_not_mu_m),
                 y3 = r_not_mu_c/max(r_not_mu_c),
                 y4 = r_not_delta_e/max(r_not_delta_e),
                 y5 = r_not_sigma_s/max(r_not_sigma_s),
                 y6 = r_not_nu_c/max(r_not_nu_c),
                 y7 = r_not_beta_s/max(r_not_beta_s),
                 y8 = r_not_beta_h/max(r_not_beta_h),
                 y9 = r_not/max(r_not))

# Define custom colors for each parameter
custom_colors <- c("blue", "green", "red", "purple", "orange", "cyan", "pink", "brown", "black")

mydata_long_1 <- gather(df_1, key = "values", value = "value", -x)


library(tidyr)

pdf(file = "mansoni_r_not_constant_parameters.pdf", width = 8, height = 5)
ggplot(mydata_long_1, aes(x = x, y = value, color = values)) +
  geom_line(size = 1.2) + scale_color_manual(values = custom_colors, name = "Parameters",
                                             labels = parameters) +
  labs(title = "S. mansoni", x = 'Temperature (ºC)',y = expression(frac("R"[0], max(("R"[0])))))+
  theme(axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5), legend.position="right", legend.key = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 20),
          panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.x = element_text(color="black", size=25), axis.text.y = element_text(color="black", size=25),
        plot.title = element_text(size=20, face="bold.italic", hjust=0.5))
dev.off() 



### The impact of fecundity rate, gaussian_1987 was used

summary(fit_nu_s)
rmax <- 0.050809
topt <- 23.900624
a <- 3.961268

rate_of_change_fecundity <- (-temperature$temp+topt)/(2*a^2)*0  

### The impact of mortality rate, spain_1982 was used

#summary(fit_mu_1)
a <- -9.455e-02
b <- -3.875e-05
c <-  4.269e-01
r0 <- 3.576e-02

rate_of_change_mortality <- (a-a*b*exp(c*temperature$temp)-b*c*exp(c*temperature$temp))/(2*(1-b*exp(c*temperature$temp)))*0


### The impact of mortality rate infected, flinn_1991 was used

#summary(fit_mu_i)
a <- 100.00000
b <- -2.24358
c <-  0

rate_of_change_mortality_of_infected <- (b+2*c*temperature$temp)*(fn_sigma_s(temperature)$.fitted + 2*preds_mu_i)/(2*preds_mu_i*(fn_sigma_s(temperature)$.fitted + preds_mu_i)*(1+a+b*temperature$temp+c*temperature$temp^2)^2)


### The impact of mortality of miracidia, thomas_2017 was used

#summary(fit_mu_m)
a <- 1.666e+01
b <- -2.283e-01
c <- 1.000e+02
d <- -9.867e+01
e <- 1.744e-03

rate_of_change_mortality_of_miracidia <- (-a*b*exp(b*temperature$temp)+d*e*exp(e*temperature$temp))/(2*a*exp(b*temperature$temp)-(c+d*exp(e*temperature$temp)))

### The impact of mortality of cercaria, spain_1982 was used

#summary(fit_mu_c)
a <- 0.022342
b <- -0.004764
c <- 0.196863
r0 <- 1.317761

rate_of_change_mortality_of_cercaria <- -(a-a*b*exp(c*temperature$temp)-b*c*exp(c*temperature$temp))/(2*(1-b*exp(c*temperature$temp)))


### The impact of miracidia hatching rate, flinn_1991 was used

#summary(fit_delta_e)
a <- 4.495569
b <- -0.167971
c <- 0.002090


rate_of_change_miracidia_hatching <- -(b+2*c*temperature$temp)/(2*(1+a+b*temperature$temp+c*temperature$temp^2))

### The impact of prepatent period rate, gaussian_1987 was used

#summary(fit_sigma_s)
rmax <- 0.065309
topt <- 35
a <- 10.251833

rate_of_change_prepatent_period <-  (0.5)*(1/(fn_sigma_s(temperature)$.fitted) -1/(fn_sigma_s(temperature)$.fitted + preds_mu_i))*(-rmax*(temperature$temp-topt)/a^2)*exp(-0.5*((temperature$temp-topt)/a)^2)

### The impact of cercaria release rate, gaussian_1987 was used

#summary(fit_nu_c)
rmax <- 1870.888
topt <- 27.289
a <- 6.461

rate_of_change_cercaria_release_rate <- (-temperature$temp+topt)/(2*a^2) 

### The impact of transmission rate in snail, spain_1982 was used

#summary(fit_beta_s)
a <- 0.11782
b <- 0.12969
c <- 0.05072
r0 <- 0.05104

rate_of_change_snail_transmis_rate <- (a-a*b*exp(c*temperature$temp)-b*c*exp(c*temperature$temp))/(2*(1-b*exp(c*temperature$temp)))

### The impact of transmission rate in human, briere2_1999 was used

#summary(fit_beta_h)
tmin <- 0
tmax <- 4.500e+01
a <- 2.160e-05
b <- 7.388e-01

rate_of_change_human_transmis_rate <- 0.5*(1/temperature$temp+1/(temperature$temp-tmin)-1/(b*(tmax-temperature$temp)))


df_2 <- data.frame(x = temperature$temp, 
                 y1 = rate_of_change_mortality_of_infected,
                 y2 = rate_of_change_mortality_of_miracidia,
                 y3 = rate_of_change_mortality_of_cercaria,
                 y4 = rate_of_change_miracidia_hatching,
                 y5 = rate_of_change_prepatent_period,
                 y6 = rate_of_change_cercaria_release_rate,
                 y7 = rate_of_change_snail_transmis_rate,
                 y8 = rate_of_change_human_transmis_rate,
                 y9 = rate_of_change_mortality_of_infected + rate_of_change_mortality_of_miracidia +
                   rate_of_change_mortality_of_cercaria + rate_of_change_miracidia_hatching +
                   rate_of_change_prepatent_period + rate_of_change_cercaria_release_rate +
                   rate_of_change_snail_transmis_rate + rate_of_change_human_transmis_rate)

mydata_long_2 <- gather(df_2[20:200,], key = "values", value = "value", -x)


pdf(file = "mansoni_r_not_response_to_parameters.pdf", width = 8, height = 5)
ggplot(mydata_long_2, aes(x = x, y = value, color = values)) +
  geom_line(size = 1.2) + scale_color_manual(values = custom_colors, name = "Parameters",
                                             labels = parameters) + 
  labs(title = 'S. mansoni', x = 'Temperature (ºC)',y = expression(frac("1 dR"[0], "R"[0]*"dT")))+
  theme(axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5), legend.position="right", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 20),
        axis.text.x = element_text(color="black", size=25), axis.text.y = element_text(color="black", size=25),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=20, face="bold.italic", hjust=0.5))

dev.off() 


#my_colors <- c("blue", "green", "purple", "orange", "gray", "brown", "yellow", "black", "red")



mydata_long_plot_1 <- ggplot(mydata_long_1, aes(x = x, y = value, color =  values)) +
  geom_line(size = 1.5) + scale_color_manual(values = custom_colors, name = "Parameters",
                                             labels = parameters) + 
  labs(title = 'S. mansoni', x = 'Temperature (ºC)',y = expression("R"[0]))+
  theme(plot.margin = margin(.25, 0, .25, 0, "cm"), axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5), legend.position="right", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 25),
        axis.text.x = element_text(color="black", size=25), axis.text.y = element_text(color="black", size=20),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=25, face="bold.italic", hjust=0.5))



mydata_long_plot_2 <- ggplot(mydata_long_2, aes(x = x, y = value, color = values)) +
  geom_line(size = 1.5) + scale_color_manual(values = custom_colors, name = "Parameters",
                                             labels = parameters) + 
  labs(title = 'S. mansoni', x = 'Temperature (ºC)',y = expression(frac("1 dR"[0], "R"[0]*"dT")))+
  theme(plot.margin = margin(.25, 2, .25, 1, "cm"), axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5), legend.position="none", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 25),
        axis.text.x = element_text(color="black", size=25), axis.text.y = element_text(color="black", size=25),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=20, face="bold.italic", hjust=0.5))


library(gridExtra)

pdf(file = "combine_sensitivity_parameter_to_r_not.pdf", width = 12, height = 4)
grid.arrange(mydata_long_plot_2, mydata_long_plot_1, ncol = 2, nrow = 1)
dev.off() 


