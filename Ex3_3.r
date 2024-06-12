# Part 1: normal distribution with mean = 1 and standard deviation = 2
# Part 2: chi-squared with 1 df

### a
## Part 1
n <- 20
mu <- 1
sd <- 2
alpha <- 0.05
quantile <- qnorm(1-alpha/2)
CI <- matrix(numeric(3*n), nrow = n, ncol = 3)
for(i in 1:n){
    sample <- rnorm(n, mu, sd)
    mu_est <- mean(sample)
    sd_est <- sd(sample)
    CI[i,1:2] <- c(mu_est - quantile * sd_est/sqrt(n), mu_est + quantile * sd_est/sqrt(n))
    CI[i, 3] <- CI[i,1]<= mu & CI[i,2]>=mu
}
coverage_prob <- sum(CI[,3])/n
coverage_prob


#Part 2
CI <- matrix(numeric(3*n), nrow = n, ncol = 3)
for(i in 1:n){
    sample <- rchisq(n, 1)
    mu_est <- mean(sample)
    sd_est <- sd(sample)
    CI[i,1:2] <- c(mu_est - quantile * sd_est/sqrt(n), mu_est + quantile * sd_est/sqrt(n))
    CI[i, 3] <- CI[i,1]<= mu & CI[i,2]>=mu
}
coverage_prob <- sum(CI[,3])/n
coverage_prob

### b
# Part 2 with Standard Bootstrap
n <- 20
df <- 1
mu <- df
M <- 100000
CI <- matrix(numeric(3*n), nrow = n, ncol = 3)
alpha <- 0.05
bootstrap_func <- function(sample, m = M){
    bootstrap_est <- numeric(m)
    for(i in 1:m){
        bootstrap_sample <- sample(sample, n, replace = TRUE)
        bootstrap_est[i] <- mean(bootstrap_sample)
    }
    return(bootstrap_est)
}

for(i in 1:n){
    sample <- rchisq(n, df)
    mu_est <- mean(sample)
    bootstrap_est <- bootstrap_func(sample)
    quants  <- quantile(bootstrap_est, c(1-alpha/2, alpha/2))
    CI[i,1:2] <- c(2*mu_est - quants[1], 2*mu_est - quants[2])
    CI[i, 3] <- CI[i,1]<= mu & CI[i,2]>=mu
}

coverage_prob <- sum(CI[,3])/n
coverage_prob


###c
# Part 2 with Bootstrap-t
# using sample variance as an unbiased estimator of asymptotic variance
start_time <- Sys.time()
n <- 20
df <- 1
mu <- df
M <- 100000
CI <- matrix(numeric(3*n), nrow = n, ncol = 3)
alpha <- 0.05
bootstrap_t_func <- function(sample, sample_est, m = M){
    bootstrap_t_est <- numeric(m)
    for(i in 1:m){
        bootstrap_sample <- sample(sample, n, replace = TRUE)
        bootstrap_est <- mean(bootstrap_sample)
        bootstrap_sd_est <- sd(bootstrap_sample)
        bootstrap_t_est[i] <- sqrt(n)*(bootstrap_est - sample_est)/bootstrap_sd_est
    }
    return(bootstrap_t_est)
}

for(i in 1:n){
    sample <- rchisq(n, df)
    mu_est <- mean(sample)
    sd_est <- sd(sample)
    bootstrap_est <- bootstrap_t_func(sample, mu_est)
    quants  <- quantile(bootstrap_est, c(1-alpha/2, alpha/2))
    CI[i,1:2] <- c(mu_est - quants[1]*sd_est/sqrt(n), mu_est - quants[2]*sd_est/sqrt(n))
    CI[i, 3] <- CI[i,1]<= mu & CI[i,2]>=mu
}
end_time <- Sys.time()
print(end_time-start_time)
coverage_prob <- sum(CI[,3])/n
coverage_prob
# with M = 10,000 and M = 100,000 it is still relatively inaccurate,
# error or are larger M needed?


## Parallelized version:

library(parallel)
start_time <- Sys.time()
n <- 20
df <- 1
M <- 100000
CI <- matrix(numeric(3*n), nrow = n, ncol = 3)
alpha <- 0.05
num_cores <- detectCores() - 1  # Use one less core than available to avoid freezing the system

bootstrap_t_func <- function(sample, sample_est, n, m = M){
    cl <- makeCluster(num_cores)
    clusterExport(cl, c("sample", "sample_est", "n"), envir = environment())
    bootstrap_t_est <- parSapply(cl, 1:m, function(i) {
        bootstrap_sample <- sample(sample, n, replace = TRUE)
        bootstrap_est <- mean(bootstrap_sample)
        bootstrap_sd_est <- sd(bootstrap_sample)
        sqrt(n) * (bootstrap_est - sample_est) / bootstrap_sd_est
    })
    stopCluster(cl)
    return(bootstrap_t_est)
}

for(i in 1:n){
    sample <- rchisq(n, df)
    mu_est <- mean(sample)
    sd_est <- sd(sample)
    bootstrap_est <- bootstrap_t_func(sample, mu_est, n)
    quants <- quantile(bootstrap_est, c(1-alpha/2, alpha/2))
    CI[i, 1:2] <- c(mu_est - quants[1] * sd_est / sqrt(n), mu_est - quants[2] * sd_est / sqrt(n))
    CI[i, 3] <- (CI[i, 1] <= df) & (CI[i, 2] >= df)
    print(i)
}
end_time <- Sys.time()
print(end_time-start_time)

coverage_prob <- sum(CI[, 3]) / n
coverage_prob
