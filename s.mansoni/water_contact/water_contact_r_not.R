
rm(list = ls(all = TRUE)) 

source("par_set_mansoni.R")

ht <- 35   # higher mean temperature in the hot dry quarter from Sow et al. 2013 paper
lt <- 25   # lower mean temperature in the cold dry quarter from Sow et al. 2013 paper

c.ht <-  1.26  # I assume that contact rate is nerly 90% in the hot dry quarter 
c.lt <-  c.ht * 2/3  # Sow et al. shows that it is about 2/3 lower

c <- 0.28

gamma <- log(1/(c.ht-c) - 1)/log(1/(c.lt-c) - 1)

Tmed_1 <-  (gamma * lt - ht)/(gamma -1)  # first parameter

a_1 <- -log(1/(c.ht-c) - 1)/(ht - Tmed_1)   



# relative contact rate with water as a function of temperature

rwct_1 <-  function(T){
  return(1/(1 + exp(-a_1 * (T - Tmed_1)))+c)
}

c.ht <-  1.1  # I assume that contact rate is nerly 90% in the hot dry quarter 
c.lt <-  c.ht * 2/3  # Sow et al. shows that it is about 2/3 lower

gamma <- log(1/(c.ht-c) - 1)/log(1/(c.lt-c) - 1)

Tmed_2 <-  (gamma * lt - ht)/(gamma -1)  # first parameter

a_2 <- -log(1/(c.ht-c) - 1)/(ht - Tmed_2)   
# second parameter

# relative contact rate with water as a function of temperature

rwct_2 <-  function(T){
  return(1/(1 + exp(-a_2 * (T - Tmed_2)))+c)
}



#temperature boundary value
min_temp <- 12
max_temp <- 37

# Generate a sequence of temperature    
temperature <- data.frame(temp = seq(min_temp, max_temp, by = 0.1))


# combine two piece wise function of mu_i
preds_1 <- fn_mu_i_1(temperature)$.fitted[which(temperature$temp <= 38)]
preds_2 <- fn_mu_i_2(temperature)$.fitted[which(temperature$temp > 38)]

# This is mu_i function 
preds_mu_i <- c(preds_1, preds_2)  

# combine two piece wise function of mu
preds_1 <- fn_mu_1(temperature)$.fitted[which(temperature$temp <= 38)]
preds_2 <- fn_mu_2(temperature)$.fitted[which(temperature$temp > 38)]

# This is mu function 
preds_mu <- c(preds_1, preds_2) 


# This is analytic representation of R not value without control. 
r_not <- (abs(lambda * fn_beta_s(temperature)$.fitted * h * fn_delta_e(temperature)$.fitted * nu_e * fn_beta_h(temperature)$.fitted * fn_nu_c(temperature)$.fitted * fn_sigma_s(temperature)$.fitted/
                        (fn_mu_m(temperature)$.fitted * (mu_h + mu_p) * fn_mu_c(temperature)$.fitted * preds_mu_i * 
                           (fn_sigma_s(temperature)$.fitted + preds_mu_i))))^(1/2)      

# This is analytic representation of R not value without control. 
r_not_water_contact_1 <- (abs(lambda * fn_beta_s(temperature)$.fitted * h * fn_delta_e(temperature)$.fitted * nu_e * fn_beta_h(temperature)$.fitted * rwct_1(temperature$temp) * fn_nu_c(temperature)$.fitted * fn_sigma_s(temperature)$.fitted/
                        (fn_mu_m(temperature)$.fitted * (mu_h + mu_p) * fn_mu_c(temperature)$.fitted * preds_mu_i * 
                           (fn_sigma_s(temperature)$.fitted + preds_mu_i))))^(1/2)      


# This is analytic representation of R not value without control. 
r_not_water_contact_2 <- (abs(lambda * fn_beta_s(temperature)$.fitted * h * fn_delta_e(temperature)$.fitted * nu_e * fn_beta_h(temperature)$.fitted * rwct_2(temperature$temp) * fn_nu_c(temperature)$.fitted * fn_sigma_s(temperature)$.fitted/
                                (fn_mu_m(temperature)$.fitted * (mu_h + mu_p) * fn_mu_c(temperature)$.fitted * preds_mu_i * 
                                   (fn_sigma_s(temperature)$.fitted + preds_mu_i))))^(1/2)      


df <- data.frame(x = temperature$temp, y1 = r_not/max(r_not), y2 = r_not_water_contact_1/max(r_not_water_contact_1), y3 = r_not_water_contact_2/max(r_not_water_contact_2),
               y4 = rep(rwct_1(25), length(temperature$temp)), y5 = rwct_1(temperature$temp), y6 = rwct_2(temperature$temp))

r_not_plot <- ggplot(data = df, aes(x = x)) +
  geom_line(aes(y = y1, colour = "No water contact"), size = 1) +
  geom_line(aes(y = y2, colour = "Mild water contact"), size = 1) +
  geom_line(aes(y = y3, colour = "Severe water contact"), size = 1) +
  scale_colour_manual("", 
                      breaks = c("", "" , ""),
                      values = c("", "", "")) +
  geom_line(aes(x, r_not/max(r_not)), df, size = 1, colour = "red") + 
  geom_line(aes(x, r_not_water_contact_1/max(r_not_water_contact_1)), df, size = 1, colour = "green") + 
  geom_line(aes(x, r_not_water_contact_2/max(r_not_water_contact_2)), df, size = 1, colour = "orange") +
  labs(title = expression(italic("S.mansoni")), x = '',y = expression("R"[0]/max(("R"[0]))))+
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 12),
        axis.text.x = element_text(color="black", size=12), axis.text.y = element_text(color="black", size=12),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=14, face="bold.italic", hjust=0.5)) +
  
  annotate("text", x = temperature$temp[which.max(r_not)], y = 0,
           label = '^', colour = "red", size = 2) +
  annotate("text", x = temperature$temp[which.max(r_not)], y = 0.03,
           label = temperature$temp[which.max(r_not)], colour = "red", size = 2)+


annotate("text", x = temperature$temp[which.max(r_not_water_contact_1)], y = 0,
         label = '^', colour = "green", size = 2) +
  annotate("text", x = temperature$temp[which.max(r_not_water_contact_1)], y = 0.03,
           label = temperature$temp[which.max(r_not_water_contact_1)], colour = "green", size = 2)

contact_function <- ggplot(data = df, aes(x = x)) +
  geom_line(aes(y = y4, colour = "No water contact"), size = 1) +
  geom_line(aes(y = y5, colour = "Mild water contact"), size = 1) +
  geom_line(aes(y = y6, colour = "Severe water contact"), size = 1) +
  scale_colour_manual("", 
                      breaks = c("No water contact", "Mild water contact", "Severe water contact"),
                      values = c("red", "green", "orange")) +
  labs(title = '', x = 'Temperature (ºC)',y = expression("R"[0]/max(("R"[0]))))+
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 12),
        axis.text.x = element_text(color="black", size=12), axis.text.y = element_text(color="black", size=12),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=14, face="bold.italic", hjust=0.5))
  
library(gridExtra)

pdf(file = "compare_r_not_water_cont_rate_mansoni.pdf", width = 5, height = 10)
grid.arrange(r_not_plot, contact_function,  ncol=1, nrow =2)
dev.off()


