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



#The best curve for mpb 
best_mpb <- data.frame(temperature, out_puts[5,])

# upload bootstrap outputs from csv file 
mpb_boot_outs <- read.csv(file = 'mansoni_mpb_boot_out_puts.csv')

#generate empty vectors for lower bound of bootstrap
lower_bound <- c()
upper_bound <- c()

### find the confident interval for bootstrap 
for (k in 1:length(temperature)){
  lower_bound[k] <- quantile(mpb_boot_outs[k,], 0.025, na.rm = TRUE)
  upper_bound[k] <- quantile(mpb_boot_outs[k,], 0.975, na.rm = TRUE)
}

# combine the lower and upper bound 
boot_conf_preds <- data.frame(temperature, lower_bound, upper_bound)


# find the 95 percent confident interval for the optimal R_not temperature  
optim_temps <- c()

for (i in 1:length(mpb_boot_outs[1, ])){
  index <-  which.max(mpb_boot_outs[, i])
  optim_temps[i] <- temperature[index] 
}

# calculate the lower and the upper bounds for the optimal R_not temperature
lower_bound_opt_temp <- quantile(optim_temps, 0.025, na.rm = TRUE)
upper_bound_opt_temp <- quantile(optim_temps, 0.975, na.rm = TRUE)

lower_bound_opt_temp
upper_bound_opt_temp

# plot bootstrapped CIs
pdf(file = "mansoni_boot_mpb.pdf", width = 5, height = 5)
ggplot() +
  geom_line(aes(temperature,  out_puts[5,]), best_mpb, col = 'red', size = 1) +
  geom_ribbon(aes(temperature, ymin = lower_bound, ymax = upper_bound), boot_conf_preds, fill = 'red', alpha = 0.1) +
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


# Run the function above for bootstrap 
percent_positive_lower_bound <- PairPrevalence(k = 0.0105*boot_conf_preds$lower_bound^0.52, M = boot_conf_preds$lower_bound)
percent_positive_upper_bound <- PairPrevalence(k = 0.0105*boot_conf_preds$upper_bound^0.52, M = boot_conf_preds$upper_bound)
percent_positive_boot_conf_preds <- data_frame(temperature, percent_positive_lower_bound, percent_positive_upper_bound)



# Run the function above for the best R_nots  
percent_positive<- PairPrevalence(k = 0.0105*best_mpb$out_puts.5...^0.52, M = best_mpb$out_puts.5...)
percent_positive_best_mpb <- data_frame(temperature, percent_positive)


# Upload the GNTD data and choose prevalence with Bio01 
prev_percent <- read.csv(file = 'gntd_vars_all.csv')


sch_mansoni_ind <- which(prev_percent$parasite_s == "S. mansoni")


sch_mansoni <- prev_percent[sch_mansoni_ind, ]
sch_mansoni <- sch_mansoni[-which(sch_mansoni$percent_pos == 0),]
sch_mansoni <- sch_mansoni[-c(which(is.na(sch_mansoni$bio01))),]
sch_mansoni_sort <- sch_mansoni[order(sch_mansoni$bio01), ] 

# Function to generate the bins  

bins_maker <- function(percentile, lenth_of_bin){
  
  out_put <- matrix(0, nrow = 2, ncol = round(length(sch_mansoni_sort$bio01)/lenth_of_bin))
  
  for (i in 1:round(length(sch_mansoni_sort$bio01)/lenth_of_bin)){
    index_of_bin <- (1+(i-1)*lenth_of_bin):(i*lenth_of_bin)
    out_put[1, i] <- sch_mansoni_sort$bio01[index_of_bin[lenth_of_bin/2]]
    out_put[2, i] <- quantile(sch_mansoni_sort$percent_pos[index_of_bin], na.rm = TRUE, probs = percentile)
  }
  return(out_put)
}

#calculate the 95th percentile of the data in each bin

bins_data <-bins_maker(percentile = .98, lenth_of_bin = 400)
d <- data.frame(bins_data[1,], bins_data[2,])
colnames(d) <- c("temp", "prev") 


# plot bootstrapped CIs with data 
pdf(file = "conf_mansoni_boot_mpb.pdf", width = 5, height = 5)
ggplot() +
  geom_line(aes(temperature,  100*percent_positive), percent_positive_best_mpb, col = 'red', size = 1) +
  geom_ribbon(aes(temperature, ymin = 100*percent_positive_lower_bound, ymax = 100*percent_positive_upper_bound), percent_positive_boot_conf_preds, fill = 'red', alpha = 0.1) +
  geom_point(aes(temp, 100*prev*max(percent_positive_best_mpb$percent_positive)/max(prev)), d, size = 2, alpha = 0.5, colour = "red") + 
  labs(title = "S.mansoni", x = 'Temperature (ºC)',y = "Prevalence in human (%)")+
  theme(legend.position="bottom", legend.key = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour = "black"), text = element_text(size = 17),
        axis.text.x = element_text(color="black", size=17), axis.text.y = element_text(color="black", size=17),  panel.border = element_rect(colour = "black", fill=NA, size=1), 
        plot.title = element_text(size=14, face="bold.italic", hjust=0.5)) + 
  # Add mark segment
  annotate("text", x = percent_positive_best_mpb$temperature[which.max(percent_positive_best_mpb$percent_positive)], y = 1, label = '^', colour = "red", size = 5) + 
  annotate("text", x = percent_positive_best_mpb$temperature[which.max(percent_positive_best_mpb$percent_positive)], y = 4,  label = paste(percent_positive_best_mpb$temperature[which.max(percent_positive_best_mpb$percent_positive)],".0", sep = ""), colour = "red", size = 5) +
  # Add horizontal line segment
  geom_segment(aes(x = lower_bound_opt_temp, y = 0.05, xend = upper_bound_opt_temp, yend = 0.05), colour = "red", size = 1)

dev.off() 

d_shade <- data.frame(sch_mansoni$bio01, sch_mansoni$percent_pos*max(percent_positive_best_mpb$percent_positive)/max(sch_mansoni$percent_pos))
colnames(d_shade) <- c("shade_temp", "shade_prev")


# Create a new column in your data frame d for transformed values (example transformation)
d$transformed_prev <- d$prev*max(percent_positive_best_mpb$percent_positive)/max(d$prev)  # Replace with your desired transformation

pdf(file = "conf_mansoni_boot_mpb_shade.pdf", width = 5, height = 5)
ggplot() +
  geom_line(aes(temperature,  100*percent_positive), percent_positive_best_mpb, col = 'red', size = 1) +
  geom_ribbon(aes(temperature, ymin = 100*percent_positive_lower_bound, ymax = 100*percent_positive_upper_bound), percent_positive_boot_conf_preds, fill = 'red', alpha = 0.1) +
  geom_point(aes(shade_temp, 100*shade_prev*max(percent_positive_best_mpb$percent_positive)/max(shade_prev)), d_shade, size = .1, alpha = 0.1, col = "red")+
  geom_point(aes(temp, 100*transformed_prev), d, size = 2, alpha = 0.5, colour = "red") + 
  
  scale_y_continuous(name = expression("Prevalence in humans (%)"), 
                     sec.axis = sec_axis(~ . / max(100 * d$transformed_prev) * max(d$prev), name = "Temperature (ºC)")) +
  labs(title = "", x = 'Temperature (ºC)') +
  theme(legend.position = "bottom", 
        legend.key = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", colour = "black"), 
        text = element_text(size = 17),
        axis.text.x = element_text(color = "black", size = 17), 
        axis.text.y = element_text(colour = "black", size = 17),
        #axis.text.y = element_blank(),  # Hide left y-axis numerical values
        axis.title.y.right  = element_blank(),
        axis.text.y.right = element_text(color = "black", size = 17),  # Keep secondary y-axis text
        panel.border = element_rect(colour = "black", fill = NA, size = 1), 
        plot.title = element_text(size = 14, face = "bold.italic", hjust = 0.5)) +
  annotate("text", x = percent_positive_best_mpb$temperature[which.max(percent_positive_best_mpb$percent_positive)], y = 1, label = '^', colour = "red", size = 5) + 
  annotate("text", x = percent_positive_best_mpb$temperature[which.max(percent_positive_best_mpb$percent_positive)], y = 4,  label = paste(percent_positive_best_mpb$temperature[which.max(percent_positive_best_mpb$percent_positive)],".0", sep = ""), colour = "red", size = 5) +
  # Add horizontal line segment
  geom_segment(aes(x = lower_bound_opt_temp, y = 0.05, xend = upper_bound_opt_temp, yend = 0.05), colour = "red", size = 1)


dev.off() 


