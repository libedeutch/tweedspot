library(ggplot2)
library(reshape2)
library(gridExtra)


simu_data_continous <- function(loc, pattern, sigma, beta1, beta0, phi_mean = 10, p = 100, p_gamma = 10, gamma = gamma, prop_zero = 0.4, seed = NA) {
  if (!is.na(seed)) {
    set.seed(seed)
  }
  
  n <- nrow(loc)
  s <- exp(rnorm(n, mean = 0, sd = 0.2))
  
  # Generate simulated data
  
  epsilon <- matrix(rnorm(n * p, mean = 0, sd = sigma), nrow = n, ncol = p)
  phi <- rexp(p, 1 / phi_mean)
  
  # low region
  loglambda <- beta0 + epsilon
  count_null <- matrix(rnbinom(n * p, 
                               mu = s * exp(loglambda), 
                               size = matrix(phi, nrow = n, ncol = p, byrow = TRUE)), 
                       nrow = n, ncol = p)
  # gamma <- rep(FALSE, p)
  # gamma[sample(1:p, p_gamma)] <- TRUE
  
  # high region
  if (pattern == 'spot') {
    c.x <- mean(loc$x)
    c.y <- mean(loc$y)
    r <- min(diff(range(loc$x)), diff(range(loc$y))) / 3
    # add pattern
    d <- sqrt((loc$x - c.x)^2 + (loc$y - c.y)^2)
  } else if(pattern == 'linear') {
    c.x <- min(loc$x)
    c.y <- min(loc$y)
    r <- max(diff(range(loc$x)), diff(range(loc$y)))
    # add pattern
    d <- loc$x - c.x + loc$y - c.y
  }
  beta.diff <- (beta0 - beta1) / r * d + beta1 - beta0
  loglambda[, which(gamma==1)] <- loglambda[, which(gamma==1)] + sapply(beta.diff, max, 0)
  count <- matrix(rnbinom(n * p, 
                          mu = s * exp(loglambda), 
                          size = matrix(phi, nrow = n, ncol = p, byrow = TRUE)), 
                  nrow = n, ncol = p)
  
  # zero index
  H <- matrix(0, nrow = n, ncol = p)
  if (prop_zero > 0) {
    H <- matrix(rbinom(n * p, 1, prop_zero), nrow = n, ncol = p)
    count_null[which(H == 1)] <- 0
    count[which(H == 1)] <- 0
  }
  
  return(list(count = count, loc = loc, gamma = gamma, count_null = count_null, zero_index = H, phi = phi, s = s))
}


output_path <- "/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/bigSimul";



# generate data with spot pattern
spot_unit = 50
ct =3
pattern_name <- 'spot_pattern'
pattern_loc <- data.frame(x = rep(1:spot_unit, times = spot_unit*ct), y = rep(1:(spot_unit*ct), each = spot_unit))
cell_type = c(rep(1, spot_unit*spot_unit), rep(2, spot_unit*spot_unit),rep(3, spot_unit*spot_unit))#,rep(4, 16*16),rep(5, 16*16)    )
# 30 replicates 
rpls <- 15
p <- 1000
p_gamma <- 300
beta0 <- 2
beta1 <- beta0 + log(6)
phi_mean <- 10
#zero_p <- c(0, 0.1, 0.3, 0.5)
zero_p = 0.6
zero_names <- 'med'
#zero_names <- c('none','low', 'med', 'high')
sigma <- 0.3




gamma_true= matrix(0, nrow = ct, ncol = p)
gamma_true[1,101:300] = 1
gamma_true[2,201:400] = 1
gamma_true[3,301:500] = 1
# gamma_true[4,401:600] = 1
# gamma_true[5,501:700] = 1
#for (z.i in zero_p) {
z.i = zero_p
for(r in 1:rpls) {
  loc = c()
  count = c()
  s =c()
  for(ctn in 1:ct){
    index = ((ctn-1)*spot_unit*spot_unit+1):(spot_unit*spot_unit*ctn)
    gamma = gamma <- as.vector(gamma_true[ctn,])
    tmp <- simu_data_continous(loc = pattern_loc[index,], pattern = 'spot', sigma = sigma, p = p, p_gamma = p_gamma, gamma = gamma,
                               beta1 = beta1, beta0 = beta0, prop_zero = z.i, phi_mean = phi_mean, seed = 123 + r+ctn)
    count <- rbind(count, tmp$count)
    loc <- rbind(loc, as.matrix(tmp$loc))
    #gamma <- tmp$gamma 
    s <- c(s, tmp$s)
    parameters <- list(sigma = sigma, p=p, p_gamma = p_gamma, beta1 = beta1, beta0=beta0, prop_zero = z.i, phi_mean=phi_mean,cell_type = cell_type) 
  }
  save(count, loc, gamma_true, s, parameters, file = paste0(output_path,'/spot/zero_60/', pattern_name, '_zero_', z.i*100, '_replicate_', r, '.RData'))
}
