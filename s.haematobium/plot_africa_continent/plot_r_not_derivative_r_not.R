#load the parameters 
source("par_set_haematobium.R")

########################## CALCULATE R_NOT ##################


min_temp <- 12
max_temp <- 37

# Generate a sequence of temperature    
temperature <- data.frame(temp = seq(min_temp, max_temp, by = 0.1))

# combine two piece wise function of mu_i
preds_1 <- fn_mu_i_1(temperature)$.fitted[which(temperature$temp <= 33.5)]
preds_2 <- fn_mu_i_2(temperature)$.fitted[which(temperature$temp > 33.5)]

# This is mu_i function 
preds_mu_i <- c(preds_1, preds_2)  


# This is analytic representation of R not value without control.
r_not <- (abs(lambda * fn_beta_s(temperature)$.fitted * h * fn_delta_e(temperature)$.fitted * nu_e 
       * fn_beta_h(temperature)$.fitted * fn_nu_c(temperature)$.fitted * fn_sigma_s(temperature)$.fitted/
         (fn_mu_m(temperature)$.fitted * (mu_h + mu_p) * fn_mu_c(temperature)$.fitted * preds_mu_i * 
            (fn_sigma_s(temperature)$.fitted + preds_mu_i))))^(1/2)  

# Assuming you already have r_not_afun defined
r_not_afun <- approxfun(x = temperature$temp, y = r_not)

# Numerically compute the derivative
derivative_r_not <- function(x, h = 1e-6) {
  (r_not_afun(x + h) - r_not_afun(x - h)) / (2 * h)
}

derivative_r_not_afun <- approxfun(x = temperature$temp, y = derivative_r_not(temperature$temp))


# Load required libraries
library(maps)
library(ggplot2)
library(sf)
library(tidyr)
library(stars)
library(sfheaders)
library(scales)

# Get world map data
world_map <- map_data("world")

# Filter data for Africa
africa_map <- subset(world_map, region %in% c("Algeria", "Angola", "Benin", "Botswana", "Burkina Faso", "Western Sahara",
                                              "Burundi", "Cameroon", "Cape Verde", "Central African Republic",
                                              "Chad", "Comoros", "Democratic Republic of the Congo", "Republic of Congo",
                                              "Djibouti", "Egypt", "Equatorial Guinea", "Eritrea", "Eswatini", 
                                              "Ethiopia", "Gabon", "Gambia", "Ghana", "Guinea", "Guinea-Bissau", 
                                              "Ivory Coast", "Kenya", "Lesotho", "Liberia", "Libya", "Madagascar", 
                                              "Malawi", "Mali", "Mauritania", "Mauritius", "Morocco", "Mozambique", 
                                              "Namibia", "Niger", "Nigeria", "Rwanda", "Sao Tome and Principe", 
                                              "Senegal", "Seychelles", "Sierra Leone", "Somalia", "South Africa", 
                                              "South Sudan", "Sudan", "Tanzania", "Togo", "Tunisia", "Uganda", 
                                              "Zambia", "Zimbabwe"))


africa_poly <- st_union(sfheaders::sf_polygon(africa_map, x = "long", y = "lat", polygon_id = "group"))
africa_poly <- st_set_crs(africa_poly, 4326)


# upload prevalence data 
df_africa <- read.csv(file = 'africa_mpb.csv')
df_africa_mask <- read.csv(file = 'africa_mask_v1.csv')


df_africa_mask <- replace_na(complete(df_africa_mask, lon, lat), list(grid_code = 0))

df_africa_r_not <- data.frame(Longitude = df_africa$Longitude, Latitude = df_africa$Latitude, 
                              MAT = df_africa$MAT, R_not = r_not_afun(df_africa$MAT), Rate_of_increase = derivative_r_not_afun(df_africa$MAT))

rast <- st_set_crs(st_as_stars(df_africa_r_not, coords = c("Longitude", "Latitude")), 4326)
mask_rast <- st_set_crs(st_as_stars(df_africa_mask, coords = c("lat", "lon")), 4326)
mask_rast <- st_warp(mask_rast, rast)
rast$R_not[mask_rast$grid_code == 1] <- NA
rast <- rast[africa_poly]

# Write the data frame to a CSV file
#write.csv(df_africa_r_not,"s_haematobium_r_not_rate_change.csv", row.names = FALSE)





# Plot the map with intensity of color based on value
pdf(file = "s_haematobium_africa_r_not.pdf", width = 8, height = 8)
ggplot() +
  geom_stars(aes(fill = R_not), data = rast) + 
  scale_fill_gradient(low = "white", high = "red", na.value = "white") +
  
  geom_map(data = africa_map, map = africa_map,
           aes(map_id = region),
           fill = NA, color = "black", size = 0.5) + 
  labs(title = expression(paste(italic(S.haematobium) )), 
       fill = expression(paste(R[0]))) +  # Add title and color legend title
  theme_void() +
  coord_sf() +
  theme(  # Adjust plot margins
    plot.title =  element_text(hjust = 0.5, face = "bold", size = 27), 
    legend.position = c(0.3, 0.2)) # Center plot title 
dev.off()


rast$Rate_of_increase[mask_rast$grid_code == 1] <- NA




pdf(file = "s_haematobium_africa_r_not_derivative.pdf", width = 8, height = 8)

ggplot() +
  geom_stars(aes(fill = Rate_of_increase), data = rast) + 
  scale_fill_gradient2(low = "#1E90FF", high = "#FC6A03", mid = "white", na.value = "white") + 
  geom_map(data = africa_map, map = africa_map,
           aes(map_id = region),
           fill = NA, color = "black", size = 0.5) + 
  labs(title = "", 
       fill =  expression(frac("dR"[0], dT))) +  # Add title and color legend title
  theme_void() +
  coord_sf() +
  theme(  # Adjust plot margins
    plot.title = element_text(hjust = 0.5, face = "bold", size = 27), 
    legend.position = c(0.3, 0.2)) # Center plot title 

dev.off()



