# This script to estimate the beta_h transition rate 

#Data values temperature and cercariae penetrate the skin (From Foster, 1964)



function_beta_h <- function(time, m, c_i){
  beta_h <- - log(c_i)/(m*(time/24))
  return(beta_h)
}


function_beta_h(time = 2, m = 6, c_i = (100-18.5)/100)
