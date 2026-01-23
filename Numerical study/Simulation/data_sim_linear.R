library(data.table)
library(MASS)

muEvent <- function(t,cov){
  return(exp(log(0.3)+0.2*t+0.1*cov[1]+0.3*cov[2]+0.3*cov[3])) 
}

muCov <- function(t,cov){ return(2) }

muCens <- function(t,cov){ return(cens) }

jump_rate <- function(t, cov){ return(muEvent(t,cov)+muCov(t,cov)+muCens(t,cov)) }

rate_bound <- function(cov){ return(jump_rate(tstop,cov)) }

jump_simulator <- function(tn, b, cov){ 
  u <- runif(1)
  e <- rexp(1, rate = b)
  t <- e
  if (t > tstop) { return(Inf) }
  while(u > jump_rate(tn+t, cov)/b){ 
    u <- runif(1)
    e <- rexp(1, rate = b)
    t <- e+t
    if (t > tstop) { return(Inf) }
  }
  return(t)
}

simulator_n <- function(y0, t0, cov0){ 
  n <- nrow(cov0)
  
  current_max_rows <- n * 5 
  
  v_id     <- integer(current_max_rows)
  v_yn     <- numeric(current_max_rows)
  v_ynNext <- numeric(current_max_rows) 
  v_tn     <- numeric(current_max_rows)
  v_tnNext <- numeric(current_max_rows)
  v_X1     <- numeric(current_max_rows)
  v_X2     <- numeric(current_max_rows)
  v_X3     <- numeric(current_max_rows)
  
  k <- 1L # Global Row Index
  
  for(m in 1:n){
    
    y <- y0
    tn <- t0
    covm <- as.numeric(cov0[m,])
    
    bx <- if(y >= 1) 0 else rate_bound(covm)
    
    if(bx == 0) { break } 
    
    repeat{
      tnNext <- tn + jump_simulator(tn, bx, covm)
      
      if(tnNext > tstop){ break }
      
      dy <- sample(x=0:2,size=1,prob=c(muCov(tnNext,covm)/jump_rate(tnNext,covm),
                                             muEvent(tnNext,covm)/jump_rate(tnNext,covm),
                                             muCens(tnNext,covm)/jump_rate(tnNext,covm)))
      
      #CHECK SPACE AND RESIZE IF NEEDED
      if(k > current_max_rows){
        # Geometric expansion: Double the size
        new_size <- current_max_rows * 2
        length(v_id) <- new_size
        length(v_yn) <- new_size
        length(v_ynNext) <- new_size
        length(v_tn) <- new_size
        length(v_tnNext) <- new_size
        length(v_X1) <- new_size
        length(v_X2) <- new_size
        length(v_X3) <- new_size
        current_max_rows <- new_size
      }

      v_id[k]     <- m
      v_yn[k]     <- y
      v_ynNext[k] <- y + dy
      v_tn[k]     <- tn
      v_tnNext[k] <- tnNext
      v_X1[k]     <- covm[1]
      v_X2[k]     <- covm[2]
      v_X3[k]     <- covm[3]
      
      k <- k + 1L
      
      y <- y + dy
      tn <- tnNext
      
      # Update Covariates
      covm[2] <- rbinom(1,1,0.5)
      covm[3] <- covm[3] + rnorm(1,0.5,0.25)
      
      bx <- if(y >= 1) 0 else rate_bound(covm)
      
      if(bx == 0){ break }
    }
    
    if(k > current_max_rows){
      new_size <- current_max_rows * 2
      length(v_id) <- new_size 
      length(v_yn) <- new_size
      length(v_ynNext) <- new_size
      length(v_tn) <- new_size
      length(v_tnNext) <- new_size
      length(v_X1) <- new_size
      length(v_X2) <- new_size
      length(v_X3) <- new_size
      current_max_rows <- new_size
    }
    
    v_id[k]     <- m
    v_yn[k]     <- y
    v_ynNext[k] <- NA  
    v_tn[k]     <- tn
    v_tnNext[k] <- NA
    v_X1[k]     <- covm[1]
    v_X2[k]     <- covm[2]
    v_X3[k]     <- covm[3]
    
    k <- k + 1L
  }
  
  valid_indices <- 1:(k-1)
  
  simDf <- data.table(
    id     = v_id[valid_indices],
    yn     = v_yn[valid_indices],
    ynNext = v_ynNext[valid_indices],
    tn     = v_tn[valid_indices],
    tnNext = v_tnNext[valid_indices],
    X1     = v_X1[valid_indices],
    X2     = v_X2[valid_indices],
    X3     = v_X3[valid_indices]
  )
  
  return(simDf)
}