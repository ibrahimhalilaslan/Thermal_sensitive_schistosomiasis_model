#This script is to run thermal sensitive model of Schistosomiasis, June 22, 2022. 

############################# RUN MODEL  ####################################### 

# Set the working directory 
#setwd("~/Documents/research/HMS/r_scripts/june_22_version_R_script")

# run the parameter script
# source("par_set.R")


# Set the number of year to run the model. 
year <- 1

# set time span 
run_time <- seq(from = 0, to = 366*year-1, by = 1)

#we run the model with constant temperature first. 25 degree, this is approximately annual water temperature 
# Set constant temperature

sample_parameters <- seq(from = 0, to = 366*year, by = 1)
temp_data <- data.frame(temp = rep(25, length(sample_parameters)))

# combine two piece wise function of mu_i
preds_1 <- fn_mu_i_1(temp_data)$.fitted[which(temp_data$temp <= 36)]
preds_2 <- fn_mu_i_2(temp_data)$.fitted[which(temp_data$temp > 36)]

# This is mu_i function 
preds_mu_i <- c(preds_1, preds_2)  

# combine two piece wise function of mu
preds_1 <- fn_mu_1(temp_data)$.fitted[which(temp_data$temp <= 34)]
preds_2 <- fn_mu_2(temp_data)$.fitted[which(temp_data$temp > 34)]

# This is mu function 
preds_mu <- c(preds_1, preds_2)  

# Generate linearly interpolate point with temperature dependent parameter function to use them in the solution of diff equations.  
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


#Call the necessary library 
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
           
           dS <- (nu_s - (S + E + I) * nu) * (S + r * E) - lambda * (1-exp(-(beta_s * delta_e * nu_e * h * P_m)/((mu_m * (S + E + I))))) * S - mu * S
           dE <- lambda * (1-exp(-(beta_s * delta_e * nu_e * h * P_m)/(mu_m * (S + E + I)))) * S - (sigma_s + mu_i) * E
           dI <- sigma_s * E -  mu_i * I
           dP <-  beta_h * (nu_c/mu_c) * I - sigma_p * P
           dP_m <- sigma_p * P - (mu_h + mu_p) * P_m
           return(list(c(S = dS, E = dE, I = dI, P = dP, P_m = dP_m)))
         })
}

# Set the initial conditions 
y0 <- c(S = (60434 + 33232)/2, E = 1260, I = 2340, P = 2, P_m = 105)

parms0 <- c(nu, lambda)  
## solve the system 
model_outputs <- ode(y = y0, times = run_time, func = thermal_sensitive_model, parms = parms0)

#plot the simulations and print 
pdf(file="simulation_with_constant_temp.pdf")
par(mfrow=c(2,3))
plot(model_outputs[,1], model_outputs[,2], type = "l",  cex.axis = 2,  lwd = 2, col = 2, xlab = "", ylab = "",  main = "", ylim = c(min(model_outputs[,2]-10), max(model_outputs[,2]+10)))
mtext(text = "Sus. snail", side = 1, line = -19, cex = 1.5)
plot(model_outputs[,1], model_outputs[,3], type = "l", cex.axis = 2, cex.lab = 1.7, lwd = 2, col = 2, xlab = "", ylab = "", main = "", ylim = c(min(model_outputs[,3]-5), max(model_outputs[,3]+5)))
mtext(text = "Exp. snail", side = 1, line = -19, cex = 1.5)
mtext(text = "Day", side = 1, line = 4, cex = 1.5)
plot(model_outputs[,1], model_outputs[,4], type = "l", cex.axis = 2, cex.lab = 1.7, lwd = 2, col = 2, xlab = "", ylab = "", main = "", ylim = c(min(model_outputs[,4]-10), max(model_outputs[,4]+10)))
mtext(text = "Infec. snail", side = 1, line = -19, cex = 1.5)
plot(model_outputs[,1], model_outputs[,5], type = "l", cex.axis = 2, cex.lab = 1.7, lwd = 2, col = 2, xlab = "", ylab = "", main = "", ylim = c(min(model_outputs[,5]-2), max(model_outputs[,5]+2)))
mtext(text = "Im-mat. par.", side = 1, line = -19, cex = 1.5)
plot(model_outputs[,1], model_outputs[,6], type = "l", cex.axis = 2, cex.lab = 1.7, lwd = 2, col = 2, xlab = "", ylab = "", main = "", ylim = c(min(model_outputs[,6]-10), max(model_outputs[,6]+10)))
mtext(text = "Mat. par.", side = 1, line = -19, cex = 1.5)
mtext(text = "Day", side = 1, line = 4, cex = 1.5)
plot(model_outputs[,1], model_outputs[, 4]/sum(model_outputs[length(model_outputs[, 1]),2:4],2), cex.axis = 2, ylab = "", type = "l", cex.lab = 1.7, lwd = 2, col = "darkred", xlab = "", main = "", ylim = c(0.01, max(.1)))
mtext(text = "Prev. in snail", side = 1, line = -19, cex = 1.5)
#mtext("The simulations with constant temperature for two years", side = 3, line = -2, outer = TRUE)
dev.off()



# print the equilibrium point and prevalence 
equilibrium <- round(model_outputs[length(model_outputs[, 1]), 2:6], 0)
prevalence <- round(model_outputs[length(model_outputs[, 1]),4]/sum(model_outputs[length(model_outputs[, 1]),2:4]),2)
print(equilibrium)
print(prevalence)




# Now we want to run the model with seasonal temperature

#set medium temperature 
t_med <- 25


#Set the seasonality 
epsilon <- 0.1


#set seasonal temperature
seasonal_temperature <- data.frame(temp = t_med * (1 + epsilon * sin(2 * pi * sample_parameters/366)))

# combine two piece wise function of mu_i
preds_1 <- fn_mu_i_1(temp_data)$.fitted[which(seasonal_temperature$temp <= 36)]
preds_2 <- fn_mu_i_2(temp_data)$.fitted[which(seasonal_temperature$temp > 36)]

# This is mu_i function 
preds_mu_i <- c(preds_1, preds_2)  

# combine two piece wise function of mu
preds_1 <- fn_mu_1(temp_data)$.fitted[which(seasonal_temperature$temp <= 34)]
preds_2 <- fn_mu_2(temp_data)$.fitted[which(seasonal_temperature$temp > 34)]

# This is mu function 
preds_mu <- c(preds_1, preds_2)  

# Generate linearly interpolate point with temperature dependent parameter function   
nu_s_afun <- approxfun(x = sample_parameters, y = fn_nu_s(seasonal_temperature)$.fitted)
mu_m_afun <- approxfun(x = sample_parameters, y = fn_mu_m(seasonal_temperature)$.fitted)
mu_afun <- approxfun(x = sample_parameters, y = preds_mu)
sigma_s_afun <- approxfun(x = sample_parameters, y = fn_sigma_s(seasonal_temperature)$.fitted)
mu_i_afun <- approxfun(x = sample_parameters, y = preds_mu_i)       
nu_c_afun <- approxfun(x = sample_parameters, y = fn_nu_c(seasonal_temperature)$.fitted)
mu_c_afun <- approxfun(x = sample_parameters, y = fn_mu_c(seasonal_temperature)$.fitted)           
delta_e_afun <- approxfun(x = sample_parameters, y = fn_delta_e(seasonal_temperature)$.fitted)
beta_s_afun <- approxfun(x = sample_parameters, y = fn_beta_s(seasonal_temperature)$.fitted)
beta_h_afun <- approxfun(x = sample_parameters, y = fn_beta_h(seasonal_temperature)$.fitted)



## solve the system 
model_outputs <- ode(y = y0, times = run_time, func = thermal_sensitive_model, parms = parms0)


#Plot the system 
pdf(file="simulation_with_seasonal_temp.pdf")
par(mfrow=c(2,3))
plot(model_outputs[,1]/366, model_outputs[,2], type = "l",  cex.axis = 2,  lwd = 2, col = 2, xlab = "", ylab = "",  main = "", ylim = c(min(model_outputs[,2]-10), max(model_outputs[,2]+10)))
mtext(text = "Sus. snail", side = 1, line = -19, cex = 1.5)
plot(model_outputs[,1]/366, model_outputs[,3], type = "l", cex.axis = 2, cex.lab = 1.7, lwd = 2, col = 2, xlab = "", ylab = "", main = "", ylim = c(min(model_outputs[,3]-5), max(model_outputs[,3]+5)))
mtext(text = "Exp. snail", side = 1, line = -19, cex = 1.5)
mtext(text = "Year", side = 1, line = 4, cex = 1.5)
plot(model_outputs[,1]/366, model_outputs[,4], type = "l", cex.axis = 2, cex.lab = 1.7, lwd = 2, col = 2, xlab = "", ylab = "", main = "", ylim = c(min(model_outputs[,4]-10), max(model_outputs[,4]+10)))
mtext(text = "Infec. snail", side = 1, line = -19, cex = 1.5)
plot(model_outputs[,1]/366, model_outputs[,5], type = "l", cex.axis = 2, cex.lab = 1.7, lwd = 2, col = 2, xlab = "", ylab = "", main = "", ylim = c(min(model_outputs[,5]-2), max(model_outputs[,5]+2)))
mtext(text = "Im-mat. par.", side = 1, line = -19, cex = 1.5)
plot(model_outputs[,1]/366, model_outputs[,6], type = "l", cex.axis = 2, cex.lab = 1.7, lwd = 2, col = 2, xlab = "", ylab = "", main = "", ylim = c(min(model_outputs[,6]-10), max(model_outputs[,6]+10)))
mtext(text = "Mat. par.", side = 1, line = -19, cex = 1.5)
mtext(text = "Year", side = 1, line = 4, cex = 1.5)
plot(model_outputs[,1]/366, model_outputs[, 4]/sum(model_outputs[length(model_outputs[, 1]),2:4],2), cex.axis = 2, ylab = "", type = "l", cex.lab = 1.7, lwd = 2, col = "darkred", xlab = "", main = "", ylim = c(0.01, max(.2)))
mtext(text = "Prev. in snail", side = 1, line = -19, cex = 1.5)
#mtext("The simulations with constant temperature for two years", side = 3, line = -2, outer = TRUE)
dev.off()


# upload prevalence data 
lac_de_guiers_data <- read.csv(file = 'lac_de_guiers.csv')
head(lac_de_guiers_data)
# Note this is monthly temperature 

# Set the number of year to run the model. 
year <- 4

# set time span 
run_time <- seq(from = 0, to = 366*year-1, by = 1)

#we run the model with constant temperature first. 25 degree, this is approximately annual water temperature 
# Set constant temperature

sample_parameters <- seq(from = 0, to = 366*year, by = 1)

inter_time_data <- lac_de_guiers_data$mean


#set seasonal temperature
lac_de_guiers_temperature <- rep(inter_time_data, each = 30)


lac_de_guiers_temperature <- data.frame(temp = lac_de_guiers_temperature[(length(lac_de_guiers_temperature)-
                                                                                 length(sample_parameters)+1):length(lac_de_guiers_temperature)])

# combine two piece wise function of mu_i
preds_1 <- fn_mu_i_1(lac_de_guiers_temperature)$.fitted[which(lac_de_guiers_temperature$temp <= 36)]
preds_2 <- fn_mu_i_2(lac_de_guiers_temperature)$.fitted[which(lac_de_guiers_temperature$temp > 36)]

# This is mu_i function 
preds_mu_i <- c(preds_1, preds_2)  

# combine two piece wise function of mu
preds_1 <- fn_mu_1(lac_de_guiers_temperature)$.fitted[which(lac_de_guiers_temperature$temp <= 34)]
preds_2 <- fn_mu_2(lac_de_guiers_temperature)$.fitted[which(lac_de_guiers_temperature$temp > 34)]

# This is mu function 
preds_mu <- c(preds_1, preds_2)  
# Generate linearly interpolate point with temperature dependent parameter function   
nu_s_afun <- approxfun(x = sample_parameters, y = fn_nu_s(lac_de_guiers_temperature)$.fitted)
mu_m_afun <- approxfun(x = sample_parameters, y = fn_mu_m(lac_de_guiers_temperature)$.fitted)
mu_afun <- approxfun(x = sample_parameters, y = preds_mu)
sigma_s_afun <- approxfun(x = sample_parameters, y = fn_sigma_s(lac_de_guiers_temperature)$.fitted)
mu_i_afun <- approxfun(x = sample_parameters, y = preds_mu_i)       
nu_c_afun <- approxfun(x = sample_parameters, y = fn_nu_c(lac_de_guiers_temperature)$.fitted)
mu_c_afun <- approxfun(x = sample_parameters, y = fn_mu_c(lac_de_guiers_temperature)$.fitted)           
delta_e_afun <- approxfun(x = sample_parameters, y = fn_delta_e(lac_de_guiers_temperature)$.fitted)
beta_s_afun <- approxfun(x = sample_parameters, y = fn_beta_s(lac_de_guiers_temperature)$.fitted)
beta_h_afun <- approxfun(x = sample_parameters, y = fn_beta_h(lac_de_guiers_temperature)$.fitted)

## solve the system 
model_outputs <- ode(y = y0, times = run_time, func = thermal_sensitive_model, parms = parms0)

#Plot the system 
pdf(file="simulation_with_lac_de_guiers_data.pdf")
par(mfrow=c(2,3))
plot(model_outputs[,1]/366, model_outputs[,2], type = "l",  cex.axis = 2,  lwd = 2, col = 2, xlab = "", ylab = "",  main = "", ylim = c(min(model_outputs[,2]-10), max(model_outputs[,2]+10)))
mtext(text = "Sus. snail", side = 1, line = -19, cex = 1.5)
plot(model_outputs[,1]/366, model_outputs[,3], type = "l", cex.axis = 2, cex.lab = 1.7, lwd = 2, col = 2, xlab = "", ylab = "", main = "", ylim = c(min(model_outputs[,3]-5), max(model_outputs[,3]+5)))
mtext(text = "Exp. snail", side = 1, line = -19, cex = 1.5)
mtext(text = "Year", side = 1, line = 4, cex = 1.5)
plot(model_outputs[,1]/366, model_outputs[,4], type = "l", cex.axis = 2, cex.lab = 1.7, lwd = 2, col = 2, xlab = "", ylab = "", main = "", ylim = c(min(model_outputs[,4]-10), max(model_outputs[,4]+10)))
mtext(text = "Infec. snail", side = 1, line = -19, cex = 1.5)
plot(model_outputs[,1]/366, model_outputs[,5], type = "l", cex.axis = 2, cex.lab = 1.7, lwd = 2, col = 2, xlab = "", ylab = "", main = "", ylim = c(min(model_outputs[,5]-2), max(model_outputs[,5]+2)))
mtext(text = "Im-mat. par.", side = 1, line = -19, cex = 1.5)
plot(model_outputs[,1]/366, model_outputs[,6], type = "l", cex.axis = 2, cex.lab = 1.7, lwd = 2, col = 2, xlab = "", ylab = "", main = "", ylim = c(min(model_outputs[,6]-10), max(model_outputs[,6]+10)))
mtext(text = "Mat. par.", side = 1, line = -19, cex = 1.5)
mtext(text = "Year", side = 1, line = 4, cex = 1.5)
plot(model_outputs[,1]/366, model_outputs[, 4]/sum(model_outputs[length(model_outputs[, 1]),2:4],2), cex.axis = 2, ylab = "", type = "l", cex.lab = 1.7, lwd = 2, col = "darkred", xlab = "", main = "", ylim = c(0.01, max(.5)))
mtext(text = "Prev. in snail", side = 1, line = -19, cex = 1.5)
#mtext("The simulations with constant temperature for two years", side = 3, line = -2, outer = TRUE)
dev.off()



