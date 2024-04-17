# This script generate the plot of mpb with bootstraps and the best fit dara 

# run the parameter script
source("par_set_mansoni.R")

# The number of year to run the model. 
year <- 30

# set a time span run
run_time <- seq(from = 0, to = 366*year-1, by = 0.1)

# Set a sample space for temperature 
sample_parameters <- seq(from = 0, to = 366*year, by = 0.1)

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
  seasonal_temperature <- data.frame(temp = rep(temperature[i], length(sample_parameters)))
  
  
  
  # Generate linearly interpolate point with temperature dependent parameter function   
  nu_s_afun <- approxfun(x = sample_parameters, y = fn_nu_s(seasonal_temperature)$.fitted)
  mu_m_afun <- approxfun(x = sample_parameters, y = fn_mu_m(seasonal_temperature)$.fitted)
  mu_afun <- approxfun(x = sample_parameters, y = ifelse(seasonal_temperature$temp <= 36, fn_mu_1(seasonal_temperature)$.fitted,
                                                         fn_mu_2(seasonal_temperature)$.fitted))
  sigma_s_afun <- approxfun(x = sample_parameters, y = fn_sigma_s(seasonal_temperature)$.fitted)
  mu_i_afun <- approxfun(x = sample_parameters, y = ifelse(seasonal_temperature$temp <= 36, fn_mu_i_1(seasonal_temperature)$.fitted,
                                                           fn_mu_i_2(seasonal_temperature)$.fitted))       
  nu_c_afun <- approxfun(x = sample_parameters, y = fn_nu_c(seasonal_temperature)$.fitted)
  mu_c_afun <- approxfun(x = sample_parameters, y = fn_mu_c(seasonal_temperature)$.fitted)           
  delta_e_afun <- approxfun(x = sample_parameters, y = fn_delta_e(seasonal_temperature)$.fitted)
  beta_s_afun <- approxfun(x = sample_parameters, y = ifelse(fn_beta_s(seasonal_temperature)$.fitted > 0,
                                                             fn_beta_s(seasonal_temperature)$.fitted, 0))
  beta_h_afun <- approxfun(x = sample_parameters, y = fn_beta_h(seasonal_temperature)$.fitted)
  
  
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


# now we use the mpb to calculate the prevalence of diseases 

p <- c()

PairPrevalence <- function(M, k) {
  p = 1 - 2* (1+M/(2*k))^(-k) + (1+M/k)^(-k)    # fraction of humans with at least 2 parasites
  return(p)
}

percent_positive <- PairPrevalence(k = 0.0105*out_puts[5,]^0.52, M = out_puts[5,])


#The best curve for mpb 
model_precent_positive_case <- data.frame(temperature, percent_positive*100)

# Upload the GNTD data and choose prevalence with Bio01 
ntd_data <- read.csv(file = 'gntd_vars_all.csv')


sch_mansoni_ind <- which(ntd_data$parasite_s == "S. mansoni")

sch_mansoni <- ntd_data[sch_mansoni_ind, ]
sch_mansoni <- sch_mansoni[-which(sch_mansoni$percent_pos == 0),]
sch_mansoni <- sch_mansoni[-c(which(is.na(sch_mansoni$bio01))),]


ntd_data_frame <- data.frame(sch_mansoni$percent_pos, sch_mansoni$bio01)



for (i in 1:length(sch_mansoni$bio01)){
  
}

# Save as CSV with specified column names
write.csv(ntd_data_frame,"S_mansoni_ntd.csv", row.names = FALSE)
# Save as CSV with specified column names
write.csv(model_precent_positive_case,"S_mansoni_model.csv", row.names = FALSE)
