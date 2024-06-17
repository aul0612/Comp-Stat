library(parallel)
library(tictoc)

tic()
n <- 20
df <- 1
M <- 10000
B <- 10000
alpha <- 0.05
num_cores <- detectCores() - 1  # Use one less core than available

bootstrap_t_func <- function(sample, sample_est, n, m = M) {
    bootstrap_t_est <- sapply(1:m, function(i) {
        bootstrap_sample <- sample(sample, n, replace = TRUE)
        bootstrap_est <- mean(bootstrap_sample)
        bootstrap_sd_est <- sd(bootstrap_sample)
        sqrt(n) * (bootstrap_est - sample_est) / bootstrap_sd_est
    })
    return(bootstrap_t_est)
}

bootstrap_CI_func <- function(i, n, df, M, alpha) {
    sample <- rchisq(n, df)
    mu_est <- mean(sample)
    sd_est <- sd(sample)
    bootstrap_est <- bootstrap_t_func(sample, mu_est, n, M)
    quants <- quantile(bootstrap_est, c(1 - alpha / 2, alpha / 2))
    CI <- c(mu_est - quants[1] * sd_est / sqrt(n), mu_est - quants[2] * sd_est / sqrt(n))
    coverage <- (CI[1] <= df) & (CI[2] >= df)
    return(c(CI, coverage))
}

cl <- makeCluster(num_cores)
clusterExport(cl, c("n", "df", "M", "alpha", "rchisq", "mean", "sd", "quantile", "sqrt", "bootstrap_t_func"))

CI <- parLapply(cl, 1:B, bootstrap_CI_func, n, df, M, alpha)
stopCluster(cl)

CI <- do.call(rbind, CI)
coverage_prob <- mean(CI[, 3])

toc()
coverage_prob
