## Simple wrapper for PROSAIL
# setwd("../src/RTM")
dyn.load("prosail.so")

prosp.def <- c("N" = 1.5, 
               "Cab" = 40, 
               "Car" = 8, 
               "Cbrown" = 0, 
               "Cw" = 0.01, 
               "Cm" = 0.009)
sail.def = c(prosp.def,
             "LIDFa" = -0.35, 
             "LIDFb" = -0.15,
             "LIDFtype" = as.integer(1), 
             "LAI" = 3,
             "q" = 0.01,
             "tts" = 30,
             "tto" = 10,
             "psi" = 0,
             "psoil" = 1)

pro4sail <- function(params, constants){
    param.order <- c("N", "Cab", "Car", "Cbrown", "Cw", "Cm",
                     "LIDFa", "LIDFb", "LIDFtype", "LAI", "q",
                     "tts", "tto", "psi", "psoil")
    pass.in <- c(as.list(params), as.list(constants))
    pass.in <- pass.in[param.order]
    r <- numeric(2101)
    p <- c(list("PRO4SAIL"), pass.in, rep(list(r),4))
    lp <- length(p)-1
    f <- do.call(.Fortran, p)
    #     out <- do.call(cbind, f[(lp-3):lp])
    out <- f[[(lp-3)]]
    return(out)
}

rtnorm <- function(mu, sd, MIN){
    x <- rnorm(1, mu, sd)
    if(x < MIN)
        x <- qnorm(runif(1, pnorm(MIN, mu, sd, 1, 0), 1), mu, sd, 1, 0)
    return(x)
}

dtnorm <- function(x, mu, sd, MIN){
    if(x < MIN)
        return(-1e15)
    else
        return(dnorm(x, mu, sd, 1) - log(1-pnorm(MIN, mu, sd, 1, 0)))
}

invert.sail <- function(observed, inits, constants, ngibbs, prior, pm){
    observed <- as.matrix(observed)
    nspec <- ncol(observed)
    nwl <- nrow(observed)
    npars <- length(inits)
    rp1 <- 0.001 + nspec*nwl/2
    rsd <- 0.5
    PrevSpec <- pro4sail(inits, constants)
    PrevError <- PrevSpec - observed
    Jump <- inits * 0.05
    results <- matrix(NA, nrow=ngibbs, ncol=npars+1)
    ar <- numeric(npars)
    adapt <- 20
    adj_min <- 0.1
    for(ng in 1:ngibbs){
        if(ng %% adapt < 1){
            adj <- ar / adapt / 0.75
            adj[adj < adj_min] <- adj_min
            Jump <- Jump * adj
            ar <- numeric(npars)
        }
        for(p in 1:npars){
            tvec <- inits
            tvec[p] <- rtnorm(inits[p],Jump[p],pm[p])
            TrySpec <- pro4sail(tvec, constants)
            TryError <- TrySpec - observed
            TryPost <- sum(dnorm(TryError,0,rsd,1)) + prior[[p]](tvec[p])
            PrevPost <- sum(dnorm(PrevError,0,rsd,1)) + prior[[p]](inits[p])
            JN <- dtnorm(tvec[p], inits[p], Jump[p], pm[p])
            JD <- dtnorm(inits[p], tvec[p], Jump[p], pm[p])
            a <- exp((TryPost - JN) - (PrevPost - JD))
            if(is.na(a)) a <- -1
            if(a > runif(1)){
                inits[p] <- tvec[p]
                PrevError <- TryError
                ar[p] <- ar[p] + 1
            }
            results[ng,p] <- inits[p]
        }
        rp2 <- 0.001 + sum(PrevError * PrevError)/2
        rinv <- rgamma(1, rp1, 1/rp2)
        rsd <- 1/sqrt(rinv)
        results[ng,npars+1] <- rsd
    }
    return(results)
}
            

invert.pars <- c("LAI", "Cab", "Car")
nip <- names(sail.def) %in% invert.pars
test.constants <- sail.def[!nip]
test.inits <- numeric(3)
test.inits["LAI"] <- 2
test.inits["Cab"] <- 25
test.inits["Car"] <- 5
test_obs <- pro4sail(test.inits, test.constants)
test.prior <- list()
test.prior[[1]] <- function(x) dnorm(log(x), 0, 100, 1)
test.prior[[2]] <- function(x) dnorm(log(x), 0, 100, 1)
test.prior[[3]] <- function(x) dnorm(log(x), 0, 100, 1)
test.pm <- rep(0,3)

tt <- test_inv <- invert.sail(test_obs, sail.def[nip], 
                              sail.def[!nip], 1000,
                              test.prior, test.pm)

par(mfrow=c(2,2))
plot(test_inv[,1], type='l')
plot(test_inv[,2], type='l')
plot(test_inv[,3], type='l')
plot(test_inv[,4], type='l')
