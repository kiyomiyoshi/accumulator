library(tidyverse)
library(data.table)
library(pracma)
library(matrixStats)
library(doParallel)

# for parallel computing
cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)


#'# data loading
files <- dir("data", ".*tsv$", full.names = TRUE) 
dat <- c()
for (i in 1:length(files)) {
    d <- fread(files[i], header = TRUE)
    dat <- rbind(dat, d)
}

dat <- subset(dat, dat$response_time > 0.2 & dat$response_time < 2)
dat <- subset(dat, dat$confidence > 0)
dat <- mutate(dat, choice_accuracy = ifelse(choice_accuracy == "err", 0, 1)) 
tdat <- subset(dat, dat$trial_type == "act" | dat$trial_type == "act-bad") # select target condition
rdat <- mutate(tdat, response_time = ifelse(choice_accuracy == 0, -response_time, response_time)) # sign reversal for incorrect responses
rdat <- mutate(rdat, confidence = ifelse(choice_accuracy == 0, -confidence, confidence)) 
rt <- rdat$response_time


#'# functions
# fixed parameters
N <- 1000 # trial number
K <- 2000 # evidence number
drift_std <- 0

# trial simulator
evacc <- function(v, B, tnd, tnd_var, zvar) {
    simirt <- rep(NA , N)
    simresp <- rep(NA , N)
    for (i in 1:N) {
        # drift diffusion
        vnew <- v + v * drift_std * rnorm(1)
        # define a random starting point between 0 and B * zvar for each accumulator
        z <- B * zvar * runif(2, 0, 1)
        # dv is the cumulative sum of samples
        dv_ <- cbind(cumsum(-vnew + rnorm(K)) + z[1], cumsum(vnew + rnorm(K)) + z[2])
        # take the maximum of each accumulator
        m <- rowMaxs(dv_)
        r <- max.col(dv_)
        # subtract the boundary
        m <- m - B
        # find when the maximum of each accumulator first reaches the bound
        rt_ <- which(m > 0)[1]
        if (isempty(rt_) == FALSE) { # censoring required for simulations? 
            # save RT and response
            simirt[i] <- rt_
            simresp[i] <- r[rt_]
        }
    }
    simrt <- simirt * 1e-3 + tnd + tnd_var * rnorm(length(simirt))
    simresp <- simresp - 1 # make response 0-1
    return(cbind(simrt, simresp))   
}

# ks loss function 
fit_rt_ks <- function(x, parameters) {
    rt <- parameters$rt
    startvar <- parameters$zvar
    
    drift <-      x[1]
    bound <-      x[2]
    ndtime <-     x[3]
    ndtime_std <- x[4]
    
    trials <- evacc(drift, bound, ndtime, ndtime_std, startvar)
    trials <- as.data.frame(trials)
    trials <- mutate(trials, simrt = ifelse(simresp == 0, -simrt, simrt)) # invert rt for errors
    
    ks <- ks.test(trials$simrt, rt)$statistic # Kolmogorov-Smirnov test
    if (is.nan(ks)) {
        ks <- Inf
    }
    return(ks)
}

# model fit
fit_rt <- function(rt, zvar) {
    drift = 0.05
    bound = 10
    ndtime = median(rt)
    ndtime_std = 0.04
    
    guess <- c(drift, bound, ndtime, ndtime_std) # free parameters to be estimated
    params <- list("rt" = rt, "zvar" = zvar) # inputs
    
    fit <- suppressWarnings(optim(par = guess, fn = fit_rt_ks, gr = NULL, method = "L-BFGS-B", parameters = params,
                                         lower = c(0, 0, 0, 0),
                                         upper = c(1, 20, 1, 1),
                                         control = list("maxit" = 100000,
                                                        "parscale" = c(0.005, 1, 0.05, 0.005))))
   
    est <- data.frame(drift = fit$par[1], drbound = fit$par[2], ndtime = fit$par[3], 
                      ndtime_std =  fit$par[4], startvar = params$zvar, ks = fit$value)
    return(est)
}


#'# fitting & grid search for zvar
zvarinit <- seq(0, 1, 0.1) # values for grid search
models <- foreach(i = 1:length(zvarinit), .packages = c("matrixStats", "pracma", "tidyverse")) %dopar% {
    zvar <- zvarinit[i]
    fit_rt(rt, zvar)
}

kss <- c()
for (i in 1:length(zvarinit)) {
    kss <- rbind(kss, models[[i]]$ks)
}


#'# visualization
f1 <- models[[which.min(kss)]]
simdat <- as.data.frame(evacc(f1[1, 1], f1[1, 2], f1[1, 3], f1[1, 4], f1[1, 5]))
simdat <- na.omit(simdat)
ratio = nrow(simdat) / length(rt)

g1 <- ggplot() +
    geom_histogram(tdat, mapping = aes(x = response_time,
                                      n = nrow(tdat), 
                                      y = ..count.., fill = factor(choice_accuracy)),
                   position = "identity", alpha = 0.5) + guides(fill = F) +
    geom_freqpoly(simdat, mapping = aes(x = simrt,
                                        y = ..count../ratio, 
                                        color = factor(simresp)), 
                  size = 0.6) + guides(color = F) + xlim(0, 2) +
    annotate("text", size = 3, x = 1.5, y = 30, label = paste("ks = " , {round(f1[1, 6], 3)}, sep = "")) +
    annotate("text", size = 3, x = 1.5, y = 28, label = paste("zvar = " , {round(f1[1, 5], 3)}, sep = ""))
g1

ggsave(file = "fig_rt.jpg", plot = g1, dpi = 400, width = 3, height = 3) 