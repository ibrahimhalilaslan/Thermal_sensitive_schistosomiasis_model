
rm(list = ls(all = TRUE)) 

source("par_set_mansoni.R")

ht <- 29   # higher mean temperature in the hot dry quarter from Sow et al. 2013 paper
lt <- 24   # lower mean temperature in the cold dry quarter from Sow et al. 2013 paper

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


# The number of year to run the model. 
year <- 30

# set a time span run
run_time <- seq(from = 0, to = 366*year-1, by = 1)

# Set a sample space for temperature 
sample_parameters <- seq(from = 0, to = 366*year, by = 1)

# Min temperature
min_tem <- 12
# Max temperature
max_tem <- 37

# set temperature sequence 
temperature <- seq(min_tem, max_tem, by = 0.1);

# Set a matrix to record outputs 
out_puts <- matrix(0, nrow = 6, ncol = length(temperature))

# run the model for each temperature
for (i in 1:length(temperature)){
  
  #Set a sequence for the temperature 
  temp_data <- data.frame(temp = rep(temperature[i], length(sample_parameters)))
  
  
  # combine two piece wise function of mu_i
  preds_1 <- fn_mu_i_1(temp_data)$.fitted[which(temp_data$temp <= 37)]
  preds_2 <- fn_mu_i_2(temp_data)$.fitted[which(temp_data$temp > 37)]
  
  # This is mu_i function 
  preds_mu_i <- c(preds_1, preds_2)  
  
  # combine two piece wise function of mu
  preds_1 <- fn_mu_1(temp_data)$.fitted[which(temp_data$temp <= 37)]
  preds_2 <- fn_mu_2(temp_data)$.fitted[which(temp_data$temp > 37)]
  
  # This is mu function 
  preds_mu <- c(preds_1, preds_2) 
  
  
  # Generate linearly interpolate point with temperature dependent parameter function 
  nu_s_afun <- approxfun(x = sample_parameters, y = fn_nu_s(temp_data)$.fitted)
  mu_m_afun <- approxfun(x = sample_parameters, y = fn_mu_m(temp_data)$.fitted)
  mu_afun <- approxfun(x = sample_parameters, y = preds_mu)
  sigma_s_afun <- approxfun(x = sample_parameters, y = fn_sigma_s(temp_data)$.fitted)
  mu_i_afun <- approxfun(x = sample_parameters, y = preds_mu_i)       
  nu_c_afun <- approxfun(x = sample_parameters, y = fn_nu_c(temp_data)$.fitted)
  mu_c_afun <- approxfun(x = sample_parameters, y = fn_mu_c(temp_data)$.fitted)           
  delta_e_afun <- approxfun(x = sample_parameters, y = fn_delta_e(temp_data)$.fitted)
  beta_s_afun <- approxfun(x = sample_parameters, y = fn_beta_s(temp_data)$.fitted)
  beta_h_afun <- approxfun(x = sample_parameters, y = fn_beta_h(temp_data)$.fitted)
  
  
  #Call the library 
  library(deSolve)
  
  # the diff equation solver 
  thermal_sensitive_model <- function(t, y, parms){
    with(as.list(c(y, parms)),
         {
           nu_s <- nu_s_afun(t)
           mu_m <- mu_m_afun(t)
           mu <- mu_afun(t) 
           sigma_s <- sigma_s_afun(t)
           mu_i <- mu_i_afun(t)    
           nu_c <- nu_c_afun(t)             
           mu_c <- mu_c_afun(t)
           delta_e <- delta_e_afun(t)
           beta_s <- beta_s_afun(t)
           beta_h <- beta_h_afun(t)
           
           dS <- (nu_s - (S + E + I) * nu) * (S + r * E) -  lambda * (1-exp(-(beta_s * delta_e * nu_e * h * P_m)/(mu_m * (S + E + I)))) * S - mu * S
           dE <- lambda * (1-exp(-(beta_s * delta_e * nu_e *h * P_m)/(mu_m * (S + E + I)))) * S - (sigma_s + mu_i) * E
           dI <- sigma_s * E - mu_i * I
           dP <-  beta_h * (nu_c/mu_c) * I - sigma_p * P
           dP_m <- sigma_p * P - (mu_h + mu_p) * P_m
           
           return(list(c(S = dS, E = dE, I = dI, P = dP, P_m = dP_m)))
         })
  }
  
  # Specified the parameter value. 
  parms0 <- c(nu, lambda)
  
  # Set the initial conditions 
  y0 <- c(S = (60434 + 33232)/2, E = 1285, I = 2340, P = 2, P_m = 105)
  
  ## solve the system 
  model_outputs <- ode(y = y0, times = run_time, func = thermal_sensitive_model, parms = parms0)
  
  # record the outcome 
  out_puts[1:5,i] <- model_outputs[length(model_outputs[, 1]), 2:6]
  out_puts[6,i] <-  model_outputs[length(model_outputs[, 1]),4]/sum(model_outputs[length(model_outputs[, 1]),2:4])
}
# recod mpb
mpb_out_put <- out_puts[5,]


# run the model for each temperature
for (i in 1:length(temperature)){
  
  #Set a sequence for the temperature 
  temp_data <- data.frame(temp = rep(temperature[i], length(sample_parameters)))
  
  
  # combine two piece wise function of mu_i
  preds_1 <- fn_mu_i_1(temp_data)$.fitted[which(temp_data$temp <= 37)]
  preds_2 <- fn_mu_i_2(temp_data)$.fitted[which(temp_data$temp > 37)]
  
  # This is mu_i function 
  preds_mu_i <- c(preds_1, preds_2)  
  
  # combine two piece wise function of mu
  preds_1 <- fn_mu_1(temp_data)$.fitted[which(temp_data$temp <= 37)]
  preds_2 <- fn_mu_2(temp_data)$.fitted[which(temp_data$temp > 37)]
  
  # This is mu function 
  preds_mu <- c(preds_1, preds_2) 
  
  
  # Generate linearly interpolate point with temperature dependent parameter function 
  nu_s_afun <- approxfun(x = sample_parameters, y = fn_nu_s(temp_data)$.fitted)
  mu_m_afun <- approxfun(x = sample_parameters, y = fn_mu_m(temp_data)$.fitted)
  mu_afun <- approxfun(x = sample_parameters, y = preds_mu)
  sigma_s_afun <- approxfun(x = sample_parameters, y = fn_sigma_s(temp_data)$.fitted)
  mu_i_afun <- approxfun(x = sample_parameters, y = preds_mu_i)       
  nu_c_afun <- approxfun(x = sample_parameters, y = fn_nu_c(temp_data)$.fitted)
  mu_c_afun <- approxfun(x = sample_parameters, y = fn_mu_c(temp_data)$.fitted)           
  delta_e_afun <- approxfun(x = sample_parameters, y = fn_delta_e(temp_data)$.fitted)
  beta_s_afun <- approxfun(x = sample_parameters, y = fn_beta_s(temp_data)$.fitted)
  beta_h_afun <- approxfun(x = sample_parameters, y = fn_beta_h(temp_data)$.fitted*rwct_1(temp_data$temp))
  
  
  #Call the library 
  library(deSolve)
  
  # the diff equation solver 
  thermal_sensitive_model <- function(t, y, parms){
    with(as.list(c(y, parms)),
         {
           nu_s <- nu_s_afun(t)
           mu_m <- mu_m_afun(t)
           mu <- mu_afun(t) 
           sigma_s <- sigma_s_afun(t)
           mu_i <- mu_i_afun(t)    
           nu_c <- nu_c_afun(t)             
           mu_c <- mu_c_afun(t)
           delta_e <- delta_e_afun(t)
           beta_s <- beta_s_afun(t)
           beta_h <- beta_h_afun(t)
           
           dS <- (nu_s - (S + E + I) * nu) * (S + r * E) -  lambda * (1-exp(-(beta_s * delta_e * nu_e * h * P_m)/(mu_m * (S + E + I)))) * S - mu * S
           dE <- lambda * (1-exp(-(beta_s * delta_e * nu_e *h * P_m)/(mu_m * (S + E + I)))) * S - (sigma_s + mu_i) * E
           dI <- sigma_s * E - mu_i * I
           dP <-  beta_h * (nu_c/mu_c) * I - sigma_p * P
           dP_m <- sigma_p * P - (mu_h + mu_p) * P_m
           
           return(list(c(S = dS, E = dE, I = dI, P = dP, P_m = dP_m)))
         })
  }
  
  # Specified the parameter value. 
  parms0 <- c(nu, lambda)
  
  # Set the initial conditions 
  y0 <- c(S = (60434 + 33232)/2, E = 1285, I = 2340, P = 2, P_m = 105)
  
  ## solve the system 
  model_outputs <- ode(y = y0, times = run_time, func = thermal_sensitive_model, parms = parms0)
  
  # record the outcome 
  out_puts[1:5,i] <- model_outputs[length(model_outputs[, 1]), 2:6]
  out_puts[6,i] <-  model_outputs[length(model_outputs[, 1]),4]/sum(model_outputs[length(model_outputs[, 1]),2:4])
}
# recod mpb
mpb_out_put_1 <- out_puts[5,]



# run the model for each temperature
for (i in 1:length(temperature)){
  
  #Set a sequence for the temperature 
  temp_data <- data.frame(temp = rep(temperature[i], length(sample_parameters)))
  
  
  # combine two piece wise function of mu_i
  preds_1 <- fn_mu_i_1(temp_data)$.fitted[which(temp_data$temp <= 37)]
  preds_2 <- fn_mu_i_2(temp_data)$.fitted[which(temp_data$temp > 37)]
  
  # This is mu_i function 
  preds_mu_i <- c(preds_1, preds_2)  
  
  # combine two piece wise function of mu
  preds_1 <- fn_mu_1(temp_data)$.fitted[which(temp_data$temp <= 37)]
  preds_2 <- fn_mu_2(temp_data)$.fitted[which(temp_data$temp > 37)]
  
  # This is mu function 
  preds_mu <- c(preds_1, preds_2) 
  
  
  # Generate linearly interpolate point with temperature dependent parameter function 
  nu_s_afun <- approxfun(x = sample_parameters, y = fn_nu_s(temp_data)$.fitted)
  mu_m_afun <- approxfun(x = sample_parameters, y = fn_mu_m(temp_data)$.fitted)
  mu_afun <- approxfun(x = sample_parameters, y = preds_mu)
  sigma_s_afun <- approxfun(x = sample_parameters, y = fn_sigma_s(temp_data)$.fitted)
  mu_i_afun <- approxfun(x = sample_parameters, y = preds_mu_i)       
  nu_c_afun <- approxfun(x = sample_parameters, y = fn_nu_c(temp_data)$.fitted)
  mu_c_afun <- approxfun(x = sample_parameters, y = fn_mu_c(temp_data)$.fitted)           
  delta_e_afun <- approxfun(x = sample_parameters, y = fn_delta_e(temp_data)$.fitted)
  beta_s_afun <- approxfun(x = sample_parameters, y = fn_beta_s(temp_data)$.fitted)
  beta_h_afun <- approxfun(x = sample_parameters, y = fn_beta_h(temp_data)$.fitted*rwct_2(temp_data$temp))
  
  
  #Call the library 
  library(deSolve)
  
  # the diff equation solver 
  thermal_sensitive_model <- function(t, y, parms){
    with(as.list(c(y, parms)),
         {
           nu_s <- nu_s_afun(t)
           mu_m <- mu_m_afun(t)
           mu <- mu_afun(t) 
           sigma_s <- sigma_s_afun(t)
           mu_i <- mu_i_afun(t)    
           nu_c <- nu_c_afun(t)             
           mu_c <- mu_c_afun(t)
           delta_e <- delta_e_afun(t)
           beta_s <- beta_s_afun(t)
           beta_h <- beta_h_afun(t)
           
           dS <- (nu_s - (S + E + I) * nu) * (S + r * E) -  lambda * (1-exp(-(beta_s * delta_e * nu_e * h * P_m)/(mu_m * (S + E + I)))) * S - mu * S
           dE <- lambda * (1-exp(-(beta_s * delta_e * nu_e *h * P_m)/(mu_m * (S + E + I)))) * S - (sigma_s + mu_i) * E
           dI <- sigma_s * E - mu_i * I
           dP <-  beta_h * (nu_c/mu_c) * I - sigma_p * P
           dP_m <- sigma_p * P - (mu_h + mu_p) * P_m
           
           return(list(c(S = dS, E = dE, I = dI, P = dP, P_m = dP_m)))
         })
  }
  
  # Specified the parameter value. 
  parms0 <- c(nu, lambda)
  
  # Set the initial conditions 
  y0 <- c(S = (60434 + 33232)/2, E = 1285, I = 2340, P = 2, P_m = 105)
  
  ## solve the system 
  model_outputs <- ode(y = y0, times = run_time, func = thermal_sensitive_model, parms = parms0)
  
  # record the outcome 
  out_puts[1:5,i] <- model_outputs[length(model_outputs[, 1]), 2:6]
  out_puts[6,i] <-  model_outputs[length(model_outputs[, 1]),4]/sum(model_outputs[length(model_outputs[, 1]),2:4])
}
# recod mpb
mpb_out_put_2 <- out_puts[5,]

# now we use the mpb to calculate the prevalence of diseases 

p <- c()

PairPrevalence <- function(M, k) {
  p = 1 - 2* (1+M/(2*k))^(-k) + (1+M/k)^(-k)    # fraction of humans with at least 2 parasites
  return(p)
}

# Run the function above for the best R_nots  
probability<- PairPrevalence(k = 0.0105*mpb_out_put^0.52, M = mpb_out_put)


# Run the function above for the best R_nots  
probability_1<- PairPrevalence(k = 0.0105*mpb_out_put_1^0.52, M = mpb_out_put_1)


# Run the function above for the best R_nots  
probability_2<- PairPrevalence(k = 0.0105*mpb_out_put_2^0.52, M = mpb_out_put_1)


df <- data.frame(x = temperature, y1 = r_not/max(r_not), y2 = r_not_water_contact_1/max(r_not_water_contact_1), y3 = r_not_water_contact_2/max(r_not_water_contact_2),
                 y4 = probability, y5 = probability_1, y6 = probability_2, y7 = rep(rwct_1(25), length(temperature)), y8 = rwct_1(temperature), y9 = rwct_2(temperature))

r_not_plot <- ggplot(data = df, aes(x = x)) +
  geom_line(aes(y = y1, colour = "No water contact"), size = 1.5) +
  geom_line(aes(y = y2, colour = "Mild water contact"), size = 1.5) +
  geom_line(aes(y = y3, colour = "Severe water contact"), size = 1.5) +
  scale_colour_manual("", 
                      breaks = c("", "" , ""),
                      values = c("", "", "")) +
  geom_line(aes(x, r_not/max(r_not)), df, size = 1, colour = "black") + 
  geom_line(aes(x, r_not_water_contact_1/max(r_not_water_contact_1)), df, size = 1, colour = "red") + 
  geom_line(aes(x, r_not_water_contact_2/max(r_not_water_contact_2)), df, size = 1, colour = "blue") +
  labs(title = 'S. mansoni', x ='Temperature (ºC)', y = expression("R"[0]/max(("R"[0]))))+
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 20),
        axis.text.x = element_text(color="black", size=20), axis.text.y = element_text(color="black", size=20),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=20, face="bold.italic", hjust=0.5)) +
  
  
annotate("text", x = temperature[which.max(r_not_water_contact_1)], y = 0,
         label = '^', colour = "red", size = 7) +
  annotate("text", x = temperature[which.max(r_not_water_contact_1)], y = 0.1,
           label = temperature[which.max(r_not_water_contact_1)], colour = "red", size = 7)

pdf(file = "r_not_plot.pdf", width = 4, height = 4)
r_not_plot
dev.off()


prevalence_plot <- ggplot(data = df, aes(x = x)) +
  geom_line(aes(y = 100*y4, colour = "No water contact"), size = 1.5) +
  geom_line(aes(y = 100*y5, colour = "Mild water contact"), size = 1.5) +
  geom_line(aes(y = 100*y6, colour = "Severe water contact"), size = 1.5) +
  scale_colour_manual("", 
                      breaks = c("", "" , ""),
                      values = c("", "", "")) +
  geom_line(aes(x, 100*probability), df, size = 1, colour = "black") + 
  geom_line(aes(x, 100*probability_1), df, size = 1, colour = "red") + 
  geom_line(aes(x, 100*probability_2), df, size = 1, colour = "blue") +
  labs(title = 'S. mansoni', x = 'Temperature (ºC)', y = 'Prevalence in human (%)')+
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 20),
        axis.text.x = element_text(color="black", size=20), axis.text.y = element_text(color="black", size=20),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=20, face="bold.italic", hjust=0.5)) +

  annotate("text", x = temperature[which.max(probability_1)], y = 1,
           label = '^', colour = "red", size = 7) +
  annotate("text", x = temperature[which.max(probability_1)], y = 6,
           label = paste(temperature[which.max(probability_1)], ".0", sep = ''), colour = "red", size = 7)

pdf(file = "prevalence_plot.pdf", width = 4, height = 4)
prevalence_plot
dev.off()


contact_function <- ggplot(data = df, aes(x = x)) +
  geom_line(aes(y = y7, colour = "No water contact"), size = 1.5) +
  geom_line(aes(y = y8, colour = "Mild water contact"), size = 1.5) +
  geom_line(aes(y = y9, colour = "Severe water contact"), size = 1.5) +
  scale_colour_manual("", 
                      breaks = c("", "", ""),
                      values = c("", "", "")) + 
  geom_line(aes(x, rep(rwct_1(25), length(temperature))), df, size = 1, colour = "black") + 
  geom_line(aes(x, rwct_1(temperature)), df, size = 1, colour = "red") + 
  geom_line(aes(x, rwct_2(temperature)), df, size = 1, colour = "blue") +
  labs(title = '', x = 'Temperature (ºC)',y = "Water contact rate")+
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 20),
        axis.text.x = element_text(color="black", size=20), axis.text.y = element_text(color="black", size=20),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size= 20, face="bold.italic", hjust=0.5) + ylim(0.2, 1.3)) +
    
   annotate("text", x = 13.5, y = 1.35, label = "a = 0.73", colour = "red", size=3.5)+
   annotate("text", x = 13.5, y = 1.29, label = "b = 0.28", colour = "red", size=3.5)+    
   annotate("text", x = 14.5, y = 1.23, label = expression("T"[med]== 23.66), colour = "red", size=3.5)+
  annotate("text", x = 21.5, y = 1.35, label = "a = 0.34", colour = "blue", size=3.5)+
  annotate("text", x = 21.5, y = 1.29, label = "b = 0.28", colour = "blue", size=3.5)+    
  annotate("text", x = 22.5, y = 1.23, label = expression("T"[med]== 24.54), colour = "blue", size=3.5)

contact_function

pdf(file = "contact_function.pdf", width = 4, height = 4)
contact_function
dev.off()










library(gridExtra)







pdf(file = "compare_r_not_prevalence_water_cont_rate_mansoni.pdf", width = 12, height = 4)
grid.arrange(contact_function, r_not_plot, prevalence_plot,  ncol = 3, nrow = 1)
dev.off()


