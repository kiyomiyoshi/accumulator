library(tidyverse)
library(data.table)
library(pracma)
library(matrixStats)
library(cowplot)


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


#'# functions
N <- 1000 # trial number
K <- 2000 # evidence number
drift_std <- 0
tconf <- 150 # time (tconf x 10 ms) for confidence construction 
    
# trial simulator with confidence variables
confsim <- function(v, B, tnd, tnd_var, zvar) {
    simirt <- rep(NA , N)
    simresp <- rep(NA , N)
    dce <- rep(NA , N)
    be <- rep(NA , N)
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
        if (isempty(rt_) == FALSE) {
            # save RT and response
            simirt[i] <- rt_
            simresp[i] <- r[rt_]
            dce[i] <- dv_[rt_ + tconf, r[rt_]] 
            be[i] <- dv_[rt_ + tconf, r[rt_]] - dv_[rt_ + 0, abs(r[rt_] - 3)]
        }
    }
    simrt <- simirt * 1e-3 + tnd + tnd_var * rnorm(length(simirt))
    # make response 0-1
    simresp <- simresp - 1
    return(cbind(simrt, simresp, dce, be))   
}

# ks loss function
fit_conf_ks <- function(x, parameters) {
    cv <- parameters$cv
    simresp <- parameters$simresp
    conf <- parameters$conf
    
    slope <-     x[1]
    intercept <- x[2]
    
    sigconf <- exp(slope * cv + intercept) / (1 + exp(slope * cv + intercept))
    simconf <- sigconf * (2 * simresp - 1)
    
    ks <- ks.test(simconf, conf)$statistic # Kolmogorov-Smirnov test
    if (is.nan(ks)) {
        ks <- Inf
    }
    return(ks)
}

# model fit
fit_conf <- function(cv, simresp, conf) {
    slope = 0.05
    intercept = 0.5
    
    guess <- c(slope, intercept) # free parameters to be estimated
    params <- list("cv" = cv, "simresp" = simresp, "conf" = conf) # inputs
    
    fit <- suppressWarnings(optim(par = guess, fn = fit_conf_ks, gr = NULL, method = "L-BFGS-B", parameters = params,
                                         lower = c(0, -10),
                                         upper = c(5, 10, 1, 1),
                                         control = list("maxit" = 100000,
                                                        "parscale" = c(1, 10))))
   
    est <- data.frame(slope = fit$par[1], intercept = fit$par[2], ks = fit$value)
    return(est)
}


#'# fitting & visualization (requires f1 object from fitrt.R)
simconfdat <- as.data.frame(confsim(f1[1, 1], f1[1, 2], f1[1, 3], f1[1, 4], f1[1, 5]))
simconfdat <- na.omit(simconfdat)

fd <- fit_conf(cv = simconfdat$dce, simresp = simconfdat$simresp, conf = rdat$confidence)
fb <- fit_conf(cv = simconfdat$be, simresp = simconfdat$simresp, conf = rdat$confidence)

simconfdat <- mutate(simconfdat, dceconf = exp(fd$slope * dce + fd$intercept) / (1 + exp(fd$slope * dce + fd$intercept)))
simconfdat <- mutate(simconfdat, beconf = exp(fb$slope * be + fb$intercept) / (1 + exp(fb$slope * be + fb$intercept)))

ggplot(tdat) + geom_histogram(alpha = 0.4, position = "identity", aes(x = confidence, fill = factor(choice_accuracy))) +
    xlim(0, 1) + guides(fill = F) + ggtitle("data")

ggplot(simconfdat) + geom_histogram(alpha = 0.4, position = "identity", aes(x = dce, fill = factor(simresp))) +
    xlim(-30, 60) + ylim(0, 150) + guides(fill = F) + ggtitle("dce")

ggplot(simconfdat) + geom_histogram(alpha = 0.4, position = "identity", aes(x = dceconf, fill = factor(simresp))) +
    xlim(0, 1) + ylim(0, 100) + guides(fill = F) + ggtitle("dce")

ggplot(simconfdat) + geom_histogram(alpha = 0.4, position = "identity", aes(x = be, fill = factor(simresp))) +
    xlim(-30, 100) + ylim(0, 150) + guides(fill = F) + ggtitle("be")

ggplot(simconfdat) + geom_histogram(alpha = 0.4, position = "identity", aes(x = beconf, fill = factor(simresp))) +
    xlim(0, 1) + ylim(0, 100) + guides(fill = F) + ggtitle("be")


ratio = nrow(simconfdat) / nrow(tdat)

gd <- ggplot() +
    geom_histogram(tdat, mapping = aes(x = confidence,
                                       n = nrow(tdat), 
                                       y = ..count.., fill = factor(choice_accuracy)),
                   position = "identity", alpha = 0.5) + guides(fill = F) +
    geom_freqpoly(simconfdat, mapping = aes(x = dceconf,
                                        y = ..count../ratio, 
                                        color = factor(simresp)), 
                  size = 0.6) + guides(color = F) + xlim(0, 1) + ggtitle("dce") +
    annotate("text", size = 3, x = 0.1, y = 10, label = paste("ks = " , {round(fd[1, 3], 3)}, sep = ""))


gb <- ggplot() +
    geom_histogram(tdat, mapping = aes(x = confidence,
                                       n = nrow(tdat), 
                                       y = ..count.., fill = factor(choice_accuracy)),
                   position = "identity", alpha = 0.5) + guides(fill = F) +
    geom_freqpoly(simconfdat, mapping = aes(x = beconf,
                                        y = ..count../ratio, 
                                        color = factor(simresp)), 
                  size = 0.6) + guides(color = F) + xlim(0, 1) + ggtitle("be") +
    annotate("text", size = 3, x = 0.1, y = 10, label = paste("ks = " , {round(fb[1, 3], 3)}, sep = ""))

fig_conf <- plot_grid(gd, gb, labels = c("a", "b"), align = "h", scale = 0.99, vjust = 0.87)
fig_conf
ggsave(file = "fig_conf.jpg", plot = fig_conf, dpi = 400, width = 6, height = 3)