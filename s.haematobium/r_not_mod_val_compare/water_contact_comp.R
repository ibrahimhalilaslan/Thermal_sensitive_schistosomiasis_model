# This script compare S.mansoni, S. heamatobium and previous results 
# Clear the global environment 

rm(list = ls(all = TRUE)) 

source("par_set_haematobium.R")

ht <- 29   # higher mean temperature in the hot dry quarter from Sow et al. 2013 paper
lt <- 24   # lower mean temperature in the cold dry quarter from Sow et al. 2013 paper

c.ht <-  1.4  # I assume that contact rate is nerly 90% in the hot dry quarter 
c.lt <-  c.ht * 2/3  # Sow et al. shows that it is about 2/3 lower

k <- .6

gamma <- log(1/(c.ht-k) - 1)/log(1/(c.lt-k) - 1)

Tmed <-  (gamma * lt - ht)/(gamma -1)  # first parameter

a <- log(1/(c.ht-k) - 1)/(ht - Tmed)   
# second parameter


# relative contact rate with water as a function of temperature

rwct <-  function(T){
  return(1/(1 + exp(a * (T - Tmed))) + k)
}

#temperature boundary value
min_temp <- 12
max_temp <- 38

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
r_not_water_contact <- (abs(lambda * fn_beta_s(temperature)$.fitted * h * fn_delta_e(temperature)$.fitted * nu_e * fn_beta_h(temperature)$.fitted * rwct(temperature$temp) * fn_nu_c(temperature)$.fitted * fn_sigma_s(temperature)$.fitted/
                        (fn_mu_m(temperature)$.fitted * (mu_h + mu_p) * fn_mu_c(temperature)$.fitted * preds_mu_i * 
                           (fn_sigma_s(temperature)$.fitted + preds_mu_i))))^(1/2)      


df <- data.frame(x = temperature$temp, y1 = r_not/max(r_not), y2 = r_not_water_contact/max(r_not_water_contact))

pdf(file = "compare_r_not_water_cont_rate_haematobium.pdf", width = 5, height = 5)
ggplot(data = df, aes(x = x)) +
  geom_line(aes(y = y1, colour = "No water contact"), size = 1) +
  geom_line(aes(y = y2, colour = "Water contact"), size = 1) +
  scale_colour_manual("", 
                      breaks = c("No water contact", "Water contact"),
                      values = c("blue", "green")) +
  geom_line(aes(x, r_not/max(r_not)), df, size = 1, colour = "blue") + 
  geom_line(aes(x, r_not_water_contact/max(r_not_water_contact)), df, size = 1, colour = "green") + 
  labs(title = expression(italic("S. haematobium")), x = 'Temperature (ºC)',y = expression("R"[0]/max(("R"[0]))))+
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 17),
        axis.text.x = element_text(color="black", size=17), axis.text.y = element_text(color="black", size=17),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=14, face="bold.italic", hjust=0.5)) +
  #annotate("text", x = temperature$temp[which.max(r_not)], y = 0.05,
  #     label = '^', colour = "red", size = 3) +
  #  annotate("text", x = temperature$temp[which.max(r_not)], y = 0.08,
  #           label = round(temperature$temp[which.max(r_not)]), colour = "red", size = 5) +
annotate("text", x = temperature$temp[which.max(r_not_water_contact)], y = 0.05,
         label = '^', colour = "green", size = 5) +
  annotate("text", x = temperature$temp[which.max(r_not_water_contact)], y = 0.1,
           label = round(temperature$temp[which.max(r_not_water_contact)]), colour = "green", size = 5)
dev.off() 
