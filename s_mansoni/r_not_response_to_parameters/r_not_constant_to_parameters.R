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
              * fn_beta_h(temperature)$.fitted * fn_nu_c(temperature)$.fitted * fn_sigma_s(temperature)$.fitted/
                (fn_mu_m(temperature)$.fitted * (mu_h + mu_p) * fn_mu_c(temperature)$.fitted * preds_mu_i * 
                   (mean(fn_sigma_s(temperature)$.fitted) + preds_mu_i))))^(1/2)  

# Assume the nu_c is constant 
r_not_nu_c <- (abs(lambda * fn_beta_s(temperature)$.fitted * h * fn_delta_e(temperature)$.fitted * nu_e 
              * fn_beta_h(temperature)$.fitted * fn_nu_c(temperature)$.fitted * fn_sigma_s(temperature)$.fitted/
                (fn_mu_m(temperature)$.fitted * (mu_h + mu_p) * mean(fn_mu_c(temperature)$.fitted) * preds_mu_i * 
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

parameters <- c(expression(mu[i]), expression(mu[m]), 
                expression(mu[c]), expression(delta[e]), expression(sigma[s]), expression(nu[c]), 
                expression(beta[s]), expression(beta[h]), "Overall")

df <- data.frame(x = temperature$temp, 
                 y1 = r_not_mu_i/max(r_not_mu_i),
                 y2 = r_not_mu_m/max(r_not_mu_m),
                 y3 = r_not_mu_c/max(r_not_mu_c),
                 y4 = r_not_delta_e/max(r_not_delta_e),
                 y5 = r_not_sigma_s/max(r_not_sigma_s),
                 y6 = r_not_nu_c/max(r_not_nu_c),
                 y7 = r_not_beta_s/max(r_not_beta_s),
                 y8 = r_not_beta_h/max(r_not_beta_h),
                 y9 = r_not/max(r_not))


library(tidyr)
mydata_long <- gather(df, key = "values", value = "value", -x)

pdf(file = "mansoni_r_not_constant_parameters.pdf", width = 8, height = 5)
ggplot(mydata_long, aes(x = x, y = value, color = values)) +
  geom_line(size = 1) + scale_color_discrete(name = "Parameters",
                                     labels = parameters) + 
  labs(title = "S.mansoni", x = 'Temperature (ºC)',y = expression("R"[0]/max(("R"[0]))))+
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5), legend.position="right", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 17),
        axis.text.x = element_text(color="black", size=17), axis.text.y = element_text(color="black", size=17),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=14, face="bold.italic", hjust=0.5))
dev.off() 

