# functions for loopless projections

# single loop replications -- function to draw demo rates for each rep
draw_parms = function(iter, N = TRUE){
  
  # survival and transition rates
  surv = avg.mat*surv.mask
  surv.sd = para.sd.mat*surv.mask
  
  avg.surv = matrix(get_beta_vals(16, surv, surv.sd),
                    nrow = 4, ncol = 4)

  
  # fecundity
  fec = avg.mat*fec.mask
  fec.sd = para.sd.mat*fec.mask
  
  avg.fec = matrix(rlnorm(16, log(fec), fec.sd),
                   nrow = 4, ncol = 4)
  
  # output matrices
  avg = avg.surv + avg.fec
  if(N) Ninit = round(runif(1, Ninit.min, Ninit.max))
  K = round(runif(1, K.min, K.max))
  
  if(N) {
    ret <- tibble(avg = list(avg),
                Ninit = Ninit,
                K = K)
  } else ret <- tibble(avg = list(avg),
                       K = K)
  
  return(ret)
  
}

# draw demo rates individually instead of matrix elements
draw_parms_2 = function(iter, N = TRUE){
  
  bparms <- get_beta_vals(length(beta.parms), beta.parms, beta.parms*para.cv)
  names(bparms) = names(beta.parms)
  
  rlparms <- rlnorm(length(lognorm.parms), log(lognorm.parms), log(9*0.15))
  names(rlparms) = names(lognorm.parms)
  
  parms <- c(bparms, rlparms)
  
  if(N) Ninit = round(runif(1, Ninit.min, Ninit.max))
  K = round(runif(1, K.min, K.max))
  
  if(N) {
    ret <- tibble(parms = list(parms),
                  Ninit = Ninit,
                  K = K)
  } else ret <- tibble(parms = list(parms),
                       K = K)
  
  return(ret)
  
}






# function to draw average and sd for each iteration ----
iter_matrix = function(iter){
  
  # survival and transition rates
  surv = avg.mat*surv.mask
  surv.sd = para.sd.mat*surv.mask
  
  avg.surv = matrix(get_beta_vals(16, surv, surv.sd),
                    nrow = 4, ncol = 4)
  sd.surv = matrix(rinvgauss(16, surv.sd, 1),
                   nrow = 4, ncol = 4)
  sd.surv[is.na(sd.surv)] <- 0
  
  # fecundity
  fec = avg.mat*fec.mask
  fec.sd = para.sd.mat*fec.mask
  
  avg.fec = matrix(rlnorm(16, log(fec), fec.sd),
                   nrow = 4, ncol = 4)
  sd.fec = matrix(rinvgauss(16, fec.sd, 1),
                  nrow = 4, ncol = 4)
  sd.fec[is.na(sd.fec)] <- 0
  
  # output matrices
  avg = avg.surv + avg.fec
  sd = sd.surv + sd.fec
  
  Ninit = round(runif(1, Ninit.min, Ninit.max))
  K = round(runif(1, K.min, K.max))
  
  ret <- tibble(avg = list(avg),
                sd = list(sd),
                Ninit = Ninit,
                K = K)
  
  return(ret)
}


# function to draw replicate-level average and sd ----
rep_matrix = function(avg, sd, prop.hab, Ninit){
  
  # survival and transition rates
  surv = avg*surv.mask
  surv.sd = sd*surv.mask
  
  avg.surv = matrix(get_beta_vals(16, surv, surv.sd),
                    nrow = 4, ncol = 4)
  
  # fecundity
  fec = avg*fec.mask
  fec.sd = sd*fec.mask
  
  avg.fec = matrix(rlnorm(16, log(fec), log(fec.sd)),
                   nrow = 4, ncol = 4)
  
  # output matrices
  avg = avg.surv + avg.fec
  
  Ninit_hab = round(Ninit*prop.hab)
  
  ret <- tibble(avg = list(avg),
                Ninit_hab = Ninit_hab)
  
  return(ret)
}


# project populations  ----

project_pop = function(avg, prop.hab, rate, Ninit, K){
  
  realK = numeric(30)
  realK[1] <- K
  
  for(t in 2:30){
    realK[t] <- round(realK[t-1] - rate*realK[t-1])
  }


  ### natural habitat
  # temporal variation
  surv = avg*surv.mask
  
  # to make sure all survivals < 1 (urban)
  if(length(which(surv > 1)) > 0){
    surv[which(surv > 1)] <- 1
  }
  surv.sd = surv*temp.cv
  
  fec = avg*fec.mask
  fec.sd = fec*temp.cv
  if(fec.sd[1,3] < 1){fec.sd[1,3] = 9*0.15}
  if(fec.sd[1,4] < 1){fec.sd[1,4] = 9*0.15}
  
  ### urban habitat
  avg.U <- runif(16, urban[1], urban[2])*avg
  
  # temporal variation
  surv.U = avg.U*surv.mask
  
  # to make sure all survivals < 1 (urban)
  if(length(which(surv.U > 1)) > 0){
    surv.U[which(surv.U > 1)] <- 1
  }
  surv.sd.U = surv.U*temp.cv
  
  fec.U = avg.U*fec.mask
  fec.sd.U = fec.U*temp.cv
  if(fec.sd.U[1,3] < 1){fec.sd.U[1,3] = 9*0.15}
  if(fec.sd.U[1,4] < 1){fec.sd.U[1,4] = 9*0.15}
  
  # 1 = natural, 2 = urban
  mat = array(NA, dim = c(2, n.years, 4, 4))
  for(i in 1:n.years){
    yr.surv = matrix(get_beta_vals(16, surv, surv.sd),
                     nrow = 4, ncol = 4)
    yr.fec = matrix(rlnorm(16, log(fec), log(fec.sd)),
                    nrow = 4, ncol = 4)
    if(length(which(yr.fec > max.fec)) > 0) {
      yr.fec[which(yr.fec > max.fec)] <- max.fec
    }
    yr.fec[which(is.na(yr.fec))] <- 0 
    mat[1,i,,] <- yr.surv+yr.fec
    
    
    yr.surv.U = matrix(get_beta_vals(16, surv.U, surv.sd.U),
                     nrow = 4, ncol = 4)
    yr.fec.U = matrix(rlnorm(16, log(fec.U), log(fec.sd.U)),
                    nrow = 4, ncol = 4)
    if(length(which(yr.fec.U > max.fec)) > 0) {
      yr.fec.U[which(yr.fec.U > max.fec)] <- max.fec
    }
    yr.fec.U[which(is.na(yr.fec.U))] <- 0 
    
    mat[2,i,,] <- yr.surv.U+yr.fec.U
  }
  
  # start population at stable stage distribution from average matrix
  dist = eigen.analysis(avg)$stable.stage
  dist.U = eigen.analysis(avg.U)$stable.stage
  
  # population projection matrix
  N = array(NA, dim = c(2, 4, n.years))
  N[1,,1] <- round((Ninit*(1-prop.hab))*dist)
  N[2,,1] <- round((Ninit*prop.hab)*dist.U)
  
  
  # project with temporal variation
  for(i in 2:n.years){
    
    # set fecundity = 0 if pop exceeds K
    if(sum(N[,,i-1]) > realK[i]){
      mat[,i-1,1,3:4] <- 0
    }
    
    N[1,,i] <- round(mat[1,i-1,,] %*% N[1,,i-1])
    N[2,,i] <- round(mat[2,i-1,,] %*% N[2,,i-1])
    
    switch <- round(rate*N[1,,i])
    
    N[1,,i] <- N[1,,i]-switch
    N[2,,i] <- N[2,,i]+switch
  }
  
  Ntot = apply(N, 3, sum) # total pop size
  lambda = c(Ntot[2:n.years]/Ntot[1:(n.years-1)], NA)
  
  out <- tibble(Ntot = Ntot,
                realK = realK,
                year = c(1:n.years),
                lambda = lambda) 
  
  return(out)
}


# project pop with draws on demo rates not matrix elements
project_pop_2 = function(parms, prop.hab, rate, Ninit, K){
  
  # carrying capacity and proportion urban over time
  realK = numeric(30)
  realK[1] <- K
  
  # prop = numeric(30)
  # prop[1] = prop.hab
  prop = rep(prop.hab, n.years)
  
  for(t in 2:30){
    realK[t] <- round(realK[t-1] - rate*realK[t-1])
    #prop[t] <- prop[t-1] + rate*prop[t-1]
  }
  
  
  # temporal variation
  surv = parms[!str_detect(names(parms), "F")]
  surv.sd = surv*temp.cv
  
  fec = parms[str_detect(names(parms), "F")]
  fec.sd = 9*0.15
  
  # urban habitat effect
  UE <- runif(1, urban[1], urban[2])
  
  mat = array(NA, dim = c(n.years, 4, 4))
  for(i in 1:n.years){
    
    # weighted average of natural and urban
    s = (surv*(1-prop[t])) + (surv*UE*prop[t])
    f = (fec*(1-prop[t])) + (fec*UE*prop[t])
    
    # draw year-specific values
    yr.surv = get_beta_vals(length(surv), s, surv.sd)
    names(yr.surv) = names(surv)
    
    yr.fec = rlnorm(length(fec), log(f), log(fec.sd))
    names(yr.fec) = names(fec)
      
    # limit fecundity to a maximum 
    if(length(which(yr.fec > F.max)) > 0) {
      yr.fec[which(yr.fec > F.max)] <- F.max
    }
    
    T.YY = yr.surv["S.Y"] * (1-yr.surv["G.YJ"])
    T.JY = yr.surv["S.Y"] * yr.surv["G.YJ"]
    T.JJ = yr.surv["S.J"] * (1-yr.surv["G.JS"])
    T.SJ = yr.surv["S.J"] * yr.surv["G.JS"]
    T.SS = yr.surv["S.S"] * (1-yr.surv["G.SA"])
    T.AS = yr.surv["S.S"] * yr.surv["G.SA"]
    T.AA = yr.surv["S.A"]
    F.SY = yr.fec["F.S"] * yr.surv["S.Y"]
    F.AY = yr.fec["F.A"] * yr.surv["S.Y"]
    
    mat[i,,] <- matrix(c(T.YY, 0, F.SY, F.AY,
                       T.JY, T.JJ, 0, 0,
                       0, T.SJ, T.SS, 0,
                       0, 0, T.AS, T.AA),
                     nrow = 4, ncol = 4, byrow = T)
  }
  
  # start population at stable stage distribution from 1st year matrix
  dist = eigen.analysis(mat[1,,])$stable.stage

  # population projection matrix
  N = matrix(NA, nrow = 4, ncol = n.years)
  N[,1] <- round(Ninit*dist)

  
  # project with temporal variation
  for(i in 2:n.years){
    
    # set fecundity = 0 if pop exceeds K
    if(sum(N[,i-1]) > realK[i]){
      mat[i-1,1,3:4] <- 0
    }
    
    N[,i] <- round(mat[i-1,,] %*% N[,i-1])
  }
  
  Ntot = apply(N, 2, sum) # total pop size
  lambda = c(Ntot[2:n.years]/Ntot[1:(n.years-1)], NA)
  
  out <- tibble(Ntot = Ntot,
                realK = realK,
                year = c(1:n.years),
                lambda = lambda) 
  
  return(out)
}


# population growth rate
calc_lambda = function(x){
  lam = x$N[2:length(x$N)]/x$N[1:(length(x$N)-1)]
  return(c(lam, NA))
}

geo_mean = function(lambda, N){
  if(length(which(N == 0)) > 0){
    ext = min(which(is.na(lambda)))
    y = lambda[1:(ext-2)]
  } else y = lambda[1:(length(lambda)-1)]
  n = length(y)
  z = prod(y, na.rm = T)
  z^(1/n)
}

calc_p_change = function(N){
  n = length(N)
  ((N[n]-N[1])/N[1])*100
}
