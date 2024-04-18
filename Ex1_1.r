theta_0 <- 0.5
m  <- 0
eps  <- 10^(-10)
N_H_obs <- 1
n <- 5
l_1 <- function(theta, N_H_obs = 1, n = 5){
    return((N_H_obs/theta) - (n-N_H_obs)/(1-theta))
}
l_2 <- function(theta, N_H_obs = 1, n = 5){
    return(-N_H_obs/theta^2 - (n-N_H_obs)/(1-theta)^2)
}

NR <- function(theta_0 = 0.4){
    m <- 1
    theta <- c(theta_0)
    while (abs(l_1(theta[m]))>eps){
        theta[m+1] <- theta[m] - l_1(theta[m])/l_2(theta[m])
        m <- m+1
    }
    theta_ML  <- theta[m+1]
    return(data.frame(m = seq(0:(m-1)), theta_m = theta))
}

NR()
