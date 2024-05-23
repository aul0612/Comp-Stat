n  <- 100
X <- rexp(n, rate = 3.5)


ll_1 <- function(theta){
    return(n/theta - sum(X))
}
ll_2 <- function(theta){
    return(-n/theta^2)
}
theta <- NULL
h_step <- NULL
iter <- 0
check <- TRUE
NR_exp <- function(start = 1, crit = 10^(-12)){
    theta <- c(theta, start)
    while(check){
        
        iter <- iter + 1
        h_step_new <- -1*ll_1(theta[iter])/ll_2(theta[iter])
        h_step <- c(h_step, h_step_new)
        theta_new <- theta[iter] + h_step_new
        theta <- c(theta, theta_new)

        if (abs(ll_1(theta_new)) < crit){
            check <- FALSE
        }

    }
    results <- cbind(theta, h_step)
    return(results)
}
NR_exp()

# as comparison the analytical solution for the ML-estimator is: 
theta_ML <- 1/mean(X)
theta_ML
