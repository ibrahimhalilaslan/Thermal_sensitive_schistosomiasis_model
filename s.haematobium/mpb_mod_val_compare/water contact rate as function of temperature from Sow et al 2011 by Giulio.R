# determine the parameters of a s-shaped curve that modulate contact rate with water 
# as a function of temperature, for instance
# 
# built on data from:
# Sow et al. 2011. The contribution of water contact behavior to the high 
# Schistosoma mansoni Infection rates observed in the Senegal River Basin
#
# This paper shows that contact rate in the hot dry season is about 50% higher than 
# in the lower cold season
# 
# From Google Earth Engine we have identified the following average temperatures
# in the dry season and in the cold season (by eye) 
#   
#   ht <- 29   #  mean temperature in the hot dry season
#   lt <- 24   #  mean temperature in the cold season   
#
#   c.ht <-  0.9  # this is arbitrary  
#   c.lt <-  c.ht * 2/3   # 2/3 is a qualitative estimation from Sow et al 2011
#
# The S-shaped function we will use is this one:
# 
#   c(T) = 1/(1 + exp(a * (T - Tmed)))   (1)
# 
# therefore (by dropping (T) from c(T) an dusing only c):
# 
#   1/c = 1 + exp(a * (T - Tmed))
# 
#   1/c - 1 = exp(a * (T - Tmed))
# 
# 
#   log(1/c - 1) = a * (T - Tmed)
# 
# Therefore: 
#
#   log(1/c.ht - 1) = a * (ht - Tmed)    (2)
#
# and 
# 
#   log(1/c.lt - 1) = a * (lt - Tmed)    (3)
# 
# 
# I can then devided equation (2) by eq (3)
# 
#   log(1/c.ht - 1) = a * (ht - Tmed)    (4)
#
# and 
# 
#   log(1/c.lt - 1) = a * (lt - Tmed)    (5)
# 
# 
#   log(1/c.ht - 1)/log(1/c.lt - 1) = (ht - Tmed)/(lt - Tmed)    (6)
#
#
# Lets rename gamma <- log(1/c.ht - 1)/log(1/c.lt - 1)     (7)
# 
#   gamma = (ht - Tmed)/(lt - Tmed)    (8)
# 
#   gamma * lt - gamma * Tmed = ht - Tmed   (9)
# 
#   (gamma -1) * Tmed = gamma * lt - ht
# 
#   Tmed <-  (gamma * lt - ht)/(gamma -1)
#   
#   So, finally, we can do all the computation:

ht <- 29   # higher mean temperature in the hot dry quarter from Sow et al. 2013 paper
lt <- 24   # lower mean temperature in the cold dry quarter from Sow et al. 2013 paper

c.ht <-  0.8  # I assume that contact rate is nerly 90% in the hot dry quarter 
c.lt <-  c.ht * 2/3  # Sow et al. shows that it is about 2/3 lower


gamma <- log(1/c.ht - 1)/log(1/c.lt - 1)

Tmed <-  (gamma * lt - ht)/(gamma -1)  # first parameter

a <- log(1/c.ht - 1)/(ht - Tmed)       # second parameter

# relative contact rate with water as a function of temperature

rwct <-  function(T){
  return(1/(1 + exp(a * (T - Tmed))))
}

# no lets' plot the result

plot(10:35, rwct(10:35), type = 'l', lwd = 3, col = 'red', ylim = c(0,1), 
    xlab = 'temperature [C]', ylab = 'relative contact rate')
abline(v = c(lt, ht), col = 'gray', lwd =1, lty = 2)
abline(h = c(c.lt, c.ht), col = 'gray', lwd =1, lty = 2)

# just check that it is correct
(rwct(ht)/rwct(lt)-1)*100

rwct(20)

rwct(30)/rwct(20)
