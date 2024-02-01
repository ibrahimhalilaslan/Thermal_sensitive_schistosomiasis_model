# This script generate the plot of mpb with bootstraps and the best fit dara 

# run the parameter script
source("par_set_mansoni.R")

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

ht <- 29   # higher mean temperature in the hot dry quarter from Sow et al. 2013 paper
lt <- 24   # lower mean temperature in the cold dry quarter from Sow et al. 2013 paper

c.ht <-  1.26  # I assume that contact rate is nerly 90% in the hot dry quarter 
c.lt <-  c.ht * 2/3  # Sow et al. shows that it is about 2/3 lower

c <- 0.28

gamma <- log(1/(c.ht-c) - 1)/log(1/(c.lt-c) - 1)

Tmed <-  (gamma * lt - ht)/(gamma -1)  # first parameter

a <- -log(1/(c.ht-c) - 1)/(ht - Tmed)   


# relative contact rate with water as a function of temperature

rwct <-  function(T){
  return(1/(1 + exp(a * (T - Tmed))) + k)
}


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
  beta_h_afun <- approxfun(x = sample_parameters, y = fn_beta_h(temp_data)$.fitted*rwct(temp_data$temp))
  
  
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


#The best curve for mpb 
best_mpb <- data.frame(temperature, out_puts[5,])




# plot bootstrapped CIs
pdf(file = "water_contact_mansoni_mpb.pdf", width = 5, height = 5)
ggplot() +
  geom_line(aes(temperature,  out_puts[5,]), best_mpb, col = 'red', size = 1) +
  labs(title = "S.mansoni", x = 'Temperature (ºC)',y = "Mean parasite burden")+
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 17),
        axis.text.x = element_text(color="black", size=17), axis.text.y = element_text(color="black", size=17),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=14, face="bold.italic", hjust=0.5)) 
dev.off() 



# now we use the mpb to calculate the prevalence of diseases 

p <- c()

PairPrevalence <- function(M, k) {
  p = 1 - 2* (1+M/(2*k))^(-k) + (1+M/k)^(-k)    # fraction of humans with at least 2 parasites
  return(p)
}




# Run the function above for the best R_nots  
percent_positive<- PairPrevalence(k = 0.0105*best_mpb$out_puts.5...^0.52, M = best_mpb$out_puts.5...)
percent_positive_best_mpb <- data_frame(temperature, percent_positive)



# plot bootstrapped CIs with data 
pdf(file = "water_contact_mansoni_prev.pdf", width = 5, height = 5)
ggplot() +
  geom_line(aes(temperature,  percent_positive), percent_positive_best_mpb, col = 'red', size = 1) +
  labs(title = "S.mansoni", x = 'Temperature (ºC)',y = "Percent of positive in human")+
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 17),
        axis.text.x = element_text(color="black", size=17), axis.text.y = element_text(color="black", size=17),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=14, face="bold.italic", hjust=0.5)) + 
  # Add mark segment
  annotate("text", x = percent_positive_best_mpb$temperature[which.max(percent_positive_best_mpb$percent_positive)], y = 0.05, label = '^', colour = "red", size = 5) + 
  annotate("text", x = percent_positive_best_mpb$temperature[which.max(percent_positive_best_mpb$percent_positive)], y = 0.08,  label = round(percent_positive_best_mpb$temperature[which.max(percent_positive_best_mpb$percent_positive)],1), colour = "red", size = 5)

dev.off() 


