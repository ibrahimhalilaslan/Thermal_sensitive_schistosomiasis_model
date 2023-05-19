############ This script generating thermal function for water contact rate ###########

################### COLD DRY SEASON ############################
# from November to April 

# cold dry season boys,  (0-9) ages
c_d_b_0_9 <- c(0, 0.003, 0.004, 0.005, 0.01, 0.013, 0.018, 0.016, 0.035, 0.028, 0.023, 0.019, 0.015)

# We sum all the frequency to get daily frequency 
sum_c_d_b_0_9 <- sum(c_d_b_0_9)

# cold dry season girls,  (0-9) ages
c_d_g_0_9 <- c(0.002, 0.008, 0.014, 0.022, 0.023, 0.018, 0.023, 0.015, 0.044, 0.036, 0.031, 0.026, 0.016)

# We sum all the frequency to get daily frequency 
sum_c_d_g_0_9 <- sum(c_d_g_0_9)

# Since we do not do sex discrimination we average this 
#c_d_b_g_0_9 <- (sum(c_d_b_0_9)+sum(c_d_g_0_9))/2

# Same is down below 

# cold dry season boys,  (10-19) ages
c_d_b_10_19 <- c(0.01, 0.012, 0.012, 0.012, 0.012, 0.019, 0.034, 0.039, 0.042, 0.03, 0.025,  0.028, 0.065)

sum_c_d_b_10_19 <- sum(c_d_b_10_19)

# cold dry season girls,  (10-19) ages
c_d_g_10_19 <- c(0.007, 0.05, 0.094, 0.104, 0.078, 0.053, 0.061, 0.042, 0.061, 0.07, 0.084, 0.095, 0.079)

sum_c_d_g_10_19 <- sum(c_d_g_10_19)

#c_d_b_g_10_19 <- (sum(c_d_b_10_19)+sum(c_d_g_10_19))/2

# cold dry season boys,  (20+) ages
c_d_b_20_ <- c(0.012, 0.013, 0.011, 0.009,  0.012, 0.016, 0.024, 0.022, 0.017, 0.01, 0.017, 0.019, 0.03)

sum_c_d_b_20_ <- sum(c_d_b_20_)

# cold dry season girls,  (20+) ages
c_d_g_20_ <- c(0.012, 0.033, 0.033, 0.037, 0.03, 0.031, 0.043, 0.031, 0.03, 0.03, 0.036, 0.034, 0.028)

sum_c_d_g_20_ <- sum(c_d_g_20_)

#c_d_b_g_20_ <-  (sum(c_d_b_20_)+sum(c_d_g_20_))/2

# we avarage all them 
#c_d_contact_rate <- mean(c(c_d_b_g_0_9, c_d_b_g_10_19, c_d_b_g_20_))

min(c(c(c_d_b_0_9, c_d_g_0_9, c_d_b_10_19, c_d_g_10_19, c_d_b_20_, c_d_g_20_)))

# air temperature were drawn from andy app https://andychamberlin.users.earthengine.app/view/monthlytempapp1
# from Nowember to April 
c_d_temp <- mean(c(28.991, 26.656, 24.303, 24.896, 26.773))

################### HOT DRY SEASON ############################
# Up 45 degree from April to June 


# hot dry season boys,  (0-9) ages
h_d_b_0_9 <- c(0.013, 0.013, 0.016, 0.021, 0.021, 0.028, 0.028, 0.028, 0.077, 0.037, 0.03, 0.02, 0.023)

sum_h_d_b_0_9 <- sum(h_d_b_0_9)

# hot dry season girls,  (0-9) ages
h_d_g_0_9 <- c(0.022, 0.024, 0.026, 0.033, 0.03, 0.024, 0.027, 0.022, 0.062, 0.041, 0.045, 0.036, 0.026)

sum_h_d_g_0_9 <- sum(h_d_g_0_9)

#h_d_b_g_0_9 <- (sum(h_d_b_0_9)+sum(h_d_g_0_9))/2

# hot dry season boys,  (10-19) ages
h_d_b_10_19 <- c(0.03, 0.025, 0.018, 0.018, 0.025,  0.039, 0.061, 0.045, 0.071, 0.045, 0.043, 0.035, 0.051)

sum_h_d_b_10_19 <- sum(h_d_b_10_19)

# hot dry season girls,  (10-19) ages
h_d_g_10_19 <- c(0.069, 0.125, 0.133, 0.095, 0.067, 0.054, 0.07, 0.029, 0.076, 0.062, 0.101, 0.139, 0.094)

sum_h_d_g_10_19 <- sum(h_d_g_10_19)

#h_d_b_g_10_19 <- (sum(h_d_b_10_19)+sum(h_d_g_10_19))/2

# hot dry season boys,  (20+) ages
h_d_b_20_ <- c(0.02, 0.013, 0.011, 0.013, 0.015, 0.021, 0.04, 0.021, 0.01, 0.011, 0.016, 0.022, 0.029)

sum_h_d_b_20_ <- sum(h_d_b_20_)

# hot dry season boys,  (20+) ages
h_d_g_20_ <- c(0.056, 0.06, 0.062, 0.047, 0.042, 0.047, 0.062, 0.022, 0.027, 0.031, 0.052, 0.062, 0.045)

sum_h_d_g_20_ <- sum(h_d_g_20_)

#h_d_b_g_20_ <- (sum(h_d_b_20_)+sum(h_d_g_20_))/2

# we avarage all them 
#h_d_contact_rate <- mean(c(h_d_b_g_0_9, h_d_b_g_10_19, h_d_b_g_20_))

# from April to june 
h_d_temp <- mean(c(28.266, 28.285, 29.012))

min(c(c(h_d_b_0_9, h_d_g_0_9, h_d_b_10_19, h_d_g_10_19, h_d_b_20_, h_d_g_20_)))


################### HOT WET SEASON ############################
# from July to October 

# hot wet season boys,  (0-9) ages
h_w_b_0_9 <- c(0.001, 0.014, 0.013, 0.019, 0.022, 0.031, 0.033, 0.027, 0.05, 0.037, 0.03, 0.03, 0.026)

sum_h_w_b_0_9 <- sum(h_w_b_0_9)

# hot wet season girls,  (0-9) ages
h_w_g_0_9 <- c(0.009, 0.027, 0.03, 0.032, 0.028, 0.023, 0.03, 0.02, 0.053, 0.034, 0.041, 0.039, 0.027)

sum_h_w_g_0_9 <- sum(h_w_g_0_9)

#h_w_b_g_0_9 <- (sum(h_w_b_0_9)+sum(h_w_g_0_9))/2

# hot wet season boys,  (10-19) ages
h_w_b_10_19 <- c(0.011, 0.025, 0.017, 0.02, 0.032, 0.05, 0.061, 0.044, 0.04, 0.03, 0.026, 0.032, 0.046)

sum_h_w_b_10_19 <- sum(h_w_b_10_19)

# hot wet season girls,  (10-19) ages
h_w_g_10_19 <- c(0.024, 0.103, 0.099, 0.1, 0.047, 0.043, 0.053, 0.025, 0.057, 0.054, 0.08, 0.125, 0.065)

sum_h_w_g_10_19 <- sum(h_w_g_10_19)

#h_w_b_g_10_19 <- (sum(h_w_b_10_19)+sum(h_w_g_10_19))/2


# hot wet season boys,  (20+) ages
h_w_b_20_ <- c(0.027, 0.023, 0.017, 0.011, 0.02, 0.038, 0.041, 0.023, 0.01, 0.007, 0.015, 0.018, 0.027)

sum_h_w_b_20_ <- sum(h_w_b_20_)

# hot wet season girl,  (20+) ages
h_w_g_20_ <- c(0.035, 0.044, 0.037, 0.039, 0.03, 0.039, 0.051, 0.019, 0.025, 0.022, 0.038, 0.047, 0.035)

sum_h_w_g_20_ <- sum(h_w_g_20_)

#h_w_b_g_20_ <- (sum(h_w_b_20_)+sum(h_w_g_20_))/2

min(c(c(h_w_b_0_9, h_w_g_0_9, h_w_b_10_19, h_w_g_10_19, h_w_b_20_, h_w_g_20_)))


# we avarage all them 
#h_w_contact_rate <- mean(c(h_w_b_g_0_9, h_w_b_g_10_19, h_w_b_g_20_))

h_w_temp <- mean(c(29.012, 29.575, 29.869, 30.729, 31.799))

# then the temperature 

senegal_temp_data <- c(c_d_temp, c_d_temp, c_d_temp, c_d_temp, c_d_temp, c_d_temp,
                       h_d_temp, h_d_temp, h_d_temp, h_d_temp, h_d_temp, h_d_temp, 
                       h_w_temp, h_w_temp, h_w_temp, h_w_temp, h_w_temp, h_w_temp)
senegal_contact_rate <- c(sum_c_d_b_0_9, sum_c_d_b_10_19, sum_c_d_b_20_,
                          sum_c_d_g_0_9, sum_c_d_g_10_19, sum_c_d_g_20_,
                          sum_h_d_b_0_9, sum_h_d_b_10_19, sum_h_d_b_20_,
                          sum_h_d_g_0_9, sum_h_d_g_10_19, sum_h_d_g_20_,
                          sum_h_w_b_0_9, sum_h_w_b_10_19, sum_h_w_b_20_,
                          sum_h_w_g_0_9, sum_h_w_g_10_19, sum_h_w_g_20_)

# To model the contact rate we will use sigmoid function f(T) = 1/(1+exp(-a*(T-25)))
# To find the best value of "a"  for the sigmoid function we find the minimum values of 
#fucntion below 


fun <-function(a, b) {(1/(1+b*exp(-a*(senegal_temp_data[1] -25)))-senegal_contact_rate[1])^2+
  (1/(1+b*exp(-a*(senegal_temp_data[2] -25)))-senegal_contact_rate[2])^2+
  (1/(1+b*exp(-a*(senegal_temp_data[3] -25)))-senegal_contact_rate[3])^2+
    (1/(1+b*exp(-a*(senegal_temp_data[4] -25)))-senegal_contact_rate[4])^2+
    (1/(1+b*exp(-a*(senegal_temp_data[5] -25)))-senegal_contact_rate[5])^2+
    (1/(1+b*exp(-a*(senegal_temp_data[6] -25)))-senegal_contact_rate[6])^2+
    (1/(1+b*exp(-a*(senegal_temp_data[7] -25)))-senegal_contact_rate[7])^2+
    (1/(1+b*exp(-a*(senegal_temp_data[8] -25)))-senegal_contact_rate[8])^2+
    (1/(1+b*exp(-a*(senegal_temp_data[9] -25)))-senegal_contact_rate[9])^2
}


a <- 0.2
b <- seq(0, 10, by = 0.01)

plot(b, fun(a,b))

b[which.min(fun(a,b))]


# the root is -0.005, then the sigmoid function is f(T) = 1/(1+exp(-0.0198*(T-25))

sel_b = 2.52

T <- seq(12, 50, by = .1)
plot(T, 1/(1+sel_b*exp(-.2*(T-25))))
points(senegal_temp_data, senegal_contact_rate)




