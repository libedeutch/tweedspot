pattern_name <- 'bc_pattern'
load("/Users/gilly/BOOST-MI/data/simulated data/generator/bc_patterns.RData")
output_path <- "/Users/gilly/Library/CloudStorage/OneDrive-Personal/cornell/simulMICT/"
n <- nrow(pattern_loc)


simu_data_discrete <- function(loc, pattern, sigma, beta1, beta0, phi_mean = 10, p = 100, p_gamma = 10, gamma= gamma, prop_zero = 0.4, seed = NA) {
  if (!is.na(seed)) {
    set.seed(seed)
  }
  
  n <- nrow(loc)
  s <- exp(rnorm(n, mean = 0, sd = 0.2))
  
  # Generate simulated data
  epsilon <- matrix(rnorm(n * p, mean = 0, sd = sigma), nrow = n, ncol = p)
  #epsilon <- matrix(0, nrow = n, ncol = p)
  phi <- rexp(p, 1 / phi_mean)
  
  # low region
  loglambda <- beta0 + epsilon
  count_null <- matrix(rnbinom(n * p, 
                               mu = s * exp(loglambda), 
                               size = matrix(phi, nrow = n, ncol = p, byrow = TRUE)), 
                       nrow = n, ncol = p)
  
  # high region
  loglambda[pattern, which(gamma==1)] <- loglambda[pattern, which(gamma==1)] + (beta1 - beta0)
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
  
  return(list(count = count, loc = loc, gamma = gamma, count_null = count_null, zero_index = H, phi = phi,s=s))
}

pattern_name <- 'bc_pattern'
#pattern_loc <- data.frame(x = rep(1:16, times = 16*5), y = rep(1:(16*5), each = 16))
#cell_type = c(rep(1, 16*16), rep(2, 16*16),rep(3, 16*16),rep(4, 16*16),rep(5, 16*16)    )

rangeY = range(pattern_loc[,2])
ylap = round(rangeY[2] - rangeY[1], 0)
pattern_loc_cur <- data.frame(x = rep(pattern_loc[,1], times = 3), y = c(pattern_loc[,2], pattern_loc[,2]+ylap, pattern_loc[,2]+ylap*2) )
pattern_cur = rep(pattern, 3)
cell_type = c(rep(1, 250), rep(2, 250),rep(3, 250)   )
# 30 replicates 
rpls <- 15
p <- 1000
p_gamma <- 300
beta0 <- 2
beta1 <- beta0 + log(3)
phi_mean <- 10
#zero_p <- c(0, 0.1, 0.3, 0.5)
zero_p = 0.3
zero_names <- 'med'
#zero_names <- c('none','low', 'med', 'high')
sigma <- 0.3

gamma_true= matrix(0, nrow = 3, ncol = p)
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
  for(ct in 1:3){
    gamma <- as.vector(gamma_true[ct,])
    index = ((ct-1)*250+1):(250*ct)
    tmp <- simu_data_discrete(loc = pattern_loc_cur[index,], pattern = pattern_cur[((ct-1)*250+1):(250*ct)], sigma = sigma, 
                               beta1 = beta1, beta0 = beta0, phi_mean = phi_mean, p= p, p_gamma = p_gamma, gamma = gamma,prop_zero = z.i, seed = 123 + r+ct)
    count <- rbind(count, tmp$count)
    loc <- rbind(loc,as.matrix(tmp$loc))
    #gamma <- tmp$gamma 
    s <- c(s, tmp$s)
    parameters <- list(sigma = sigma, p=p, p_gamma = p_gamma, beta1 = beta1, beta0=beta0, prop_zero = z.i, phi_mean=phi_mean,cell_type = cell_type) 
  }
  save(count, loc, gamma_true, s, parameters,pattern,  file = paste0(output_path,'/bc/', pattern_name, '_zero_', z.i*100, '_replicate_', r, '.RData'))
}

ct=3
index = ((ct-1)*250+1):(250*ct)
plot(pattern_loc_cur[index,1], pattern_loc_cur[index,2], col = pattern_cur[((ct-1)*250+1):(250*ct)]+1, pch = 19)
