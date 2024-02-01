# We actually do this to test Nguyen results with our model, we do simulation with increasing some of
# our parameter values as disease control strategies. 
prev_percent <- read.csv(file = 'gntd_vars_all.csv')
head(prev_percent)

sch_mansoni_ind <- which(prev_percent$parasite_s == "S. mansoni")
sch_mansoni <- prev_percent[sch_mansoni_ind, ]
sch_mansoni <- sch_mansoni[-c(which(is.na(sch_mansoni$bio01))),]



sch_haematobium_ind <- which(prev_percent$parasite_s == "S. haematobium")
sch_haematobium <- prev_percent[sch_haematobium_ind, ]
sch_haematobium <- sch_haematobium[-c(which(is.na(sch_haematobium$bio01))),]


############# Nguyen's result plot ######################
temp <- c(15.127, 15.763, 16.525, 16.949, 17.458, 17.881, 18.305, 18.814, 19.322, 19.831, 20.508, 21.186, 21.949, 22.797, 23.644, 24.153, 24.746,  25.085, 25.593,  26.102,  26.61,  27.119,  27.542,  28.136,  28.898, 29.576, 30.339,  31.441, 32.458, 33.814, 34.619)
rate <- c(0.11, 0.172, 0.282, 0.373, 0.467, 0.533, 0.621,  0.712, 0.806, 0.884, 0.953, 0.984, 0.994,  0.95, 0.84, 0.765, 0.665, 0.586,  0.502, 0.439, 0.376, 0.301, 0.254, 0.197, 0.144, 0.11, 0.072, 0.038, 0.025, 0.009, 0.009)


# keep just a single curve
.x <- data.frame(temp, rate)

nguyen <- nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                        data = .x,
                        iter = c(4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                        supp_errors = 'Y')

fn_nguyen <- function(x) augment(nguyen, newdata = x)


nguyen_pred_mansoni <- fn_nguyen(data.frame(temp = sch_mansoni$bio01))$.fitted
nguyen_pred_mansoni <-nguyen_pred_mansoni/max(nguyen_pred_mansoni)

nguyen_pred_haematobium <- fn_nguyen(data.frame(temp = sch_haematobium$bio01))$.fitted
nguyen_pred_haematobium <- max(nguyen_pred_haematobium)


source("par_set_mansoni.R") 

# This is analytic representation of R not value without control. 
our_pred_mansoni <- (abs(lambda * fn_beta_s(data.frame(temp = sch_mansoni$bio01))$.fitted * h * fn_delta_e(data.frame(temp = sch_mansoni$bio01))$.fitted * nu_e 
                   * fn_beta_h(data.frame(temp = sch_mansoni$bio01))$.fitted * fn_nu_c(data.frame(temp = sch_mansoni$bio01))$.fitted * fn_sigma_s(data.frame(temp = sch_mansoni$bio01))$.fitted/
                     (fn_mu_m(data.frame(temp = sch_mansoni$bio01))$.fitted * (mu_h + mu_p) * fn_mu_c(data.frame(temp = sch_mansoni$bio01))$.fitted * fn_mu_i_1(data.frame(temp = sch_mansoni$bio01))$.fitted * 
                        (fn_sigma_s(data.frame(temp = sch_mansoni$bio01))$.fitted + fn_mu_i_1(data.frame(temp = sch_mansoni$bio01))$.fitted))))^(1/2)

#normalize the results 
our_pred_mansoni <- our_pred_mansoni/max(our_pred_mansoni)

source("par_set_haematobium.R")

our_pred_haematobium <- (abs(lambda * fn_beta_s(data.frame(temp = sch_haematobium$bio01))$.fitted * h * fn_delta_e(data.frame(temp = sch_haematobium$bio01))$.fitted * nu_e 
                         * fn_beta_h(data.frame(temp = sch_haematobium$bio01))$.fitted * fn_nu_c(data.frame(temp = sch_haematobium$bio01))$.fitted * fn_sigma_s(data.frame(temp = sch_haematobium$bio01))$.fitted/
                           (fn_mu_m(data.frame(temp = sch_haematobium$bio01))$.fitted * (mu_h + mu_p) * fn_mu_c(data.frame(temp = sch_haematobium$bio01))$.fitted * fn_mu_i_1(data.frame(temp = sch_haematobium$bio01))$.fitted * 
                              (fn_sigma_s(data.frame(temp = sch_haematobium$bio01))$.fitted + fn_mu_i_1(data.frame(temp = sch_haematobium$bio01))$.fitted))))^(1/2)


#normalize the results 
our_pred_haematobium <- our_pred_haematobium/max(our_pred_haematobium)


# Create the csv files 
sch_mansoni_results <- data.frame(sch_mansoni$percent_pos/max(sch_mansoni$percent_pos), sch_mansoni$bio01, our_pred_mansoni, nguyen_pred_mansoni)

colnames(sch_mansoni_results) <- c("Prevalence", "Temperature", "Our Prediction", "Nguyen's Prediction")

# Save the data frame to a CSV file
write.csv(sch_mansoni_results, file = "sch_mansoni_results.csv", row.names = FALSE)


sch_haematobium_results <- data.frame(sch_haematobium$percent_pos/max(sch_haematobium$percent_pos), sch_haematobium$bio01, our_pred_haematobium, nguyen_pred_haematobium)

colnames(sch_haematobium_results) <- c("Prevalence", "Temperature", "Our Prediction", "Nguyen's Prediction")

# Save the data frame to a CSV file
write.csv(sch_haematobium_results, file = "sch_haematobium_results.csv", row.names = FALSE)





