# Matern correlation
Matern <- function(d, rho, nu=0.5) {
  if(nu==0.5){ # nu=0.5 is exponential model
    result <- exp(-d/rho)
  }else{
    d1 <- sqrt(2*nu)*d/rho
    result <- 2^(1-nu)/gamma(nu)*d1^nu*besselK(d1, nu)
  }
  result[d==0] <- 1
  result
}


# Simulate ONE segment
sim.y <- function(theta, S.dist, TT, T.burn=100) {
  phi <- theta[1]; rho <- theta[2]; sigma2 <- theta[3]; mu <- theta[4]
  if(is.na(mu)) { # fixed mean at 0
    mu <- 0
  }
  Sigma <- sigma2*Matern(S.dist, rho=rho)
  y <- matrix(0, nrow=TT+T.burn, ncol(S.dist))
  e <- rmvnorm(n=TT+T.burn, sigma=Sigma)
  y[1,] <- e[1,]
  for (i in 2:(TT+T.burn)) {
    y[i,] <- phi*y[i-1,] + e[i,]
  }
  y[-c(1:T.burn),] + mu
}


# Auxiliary function to compute pairs within s.lag distance
s.dist.pair <- function(S.dist, s.lag) {
  # for time.lag > 0
  tmp <- cbind(c(row(S.dist)*(S.dist<=s.lag)), c(col(S.dist)*(S.dist<=s.lag)), c(S.dist))
  tmp[tmp[,1]!=0,]
}

s.dist.pair.0 <- function(S.dist, s.lag) {
  # for time.lag = 0 (remove repetition)
  S.dist2 <- row(S.dist) - col(S.dist) < 0
  tmp <- cbind(c(row(S.dist)*(S.dist<=s.lag)), c(col(S.dist)*(S.dist<=s.lag)), c(S.dist))[c(S.dist2),]
  tmp[tmp[,1]!=0,]
}


# Input y and distance matrix find the D.pair (total number of pairs), use.obs (no. of time each obs used) and
# a remedy matrix for the edge effect
D.cal <- function(y, S.dist, s.lag=1, t.lag=1) {
  S <- ncol(y); TT <- nrow(y)
  D.pair <- 0
  
  obs <- 0*y
  s.comb <- s.dist.pair(S.dist, s.lag)
  s.comb.0 <- s.dist.pair.0(S.dist, s.lag)
  s.len <- unique(s.comb[,3])[-1]
  
  # Time lag 1 to max t lag
  for (t in 1:t.lag) {
    # Space-lag 0
    for (j in 1:S){
      obs[1:(TT-t), j] <- obs[1:(TT-t), j] + 1
      obs[(t+1):TT, j] <- obs[(t+1):TT, j] + 1
      x.tmp <- cbind(y[1:(TT-t), j], y[(t+1):TT, j])
      D.pair <- D.pair + nrow(x.tmp)
    }
    # Space-lag > 0
    for (s in s.len) {
      coord.tmp <- s.comb[s.comb[,3]==s,]
      for (i in 1:nrow(coord.tmp)) {
        obs[1:(TT-t), coord.tmp[i,1]] <- obs[1:(TT-t), coord.tmp[i,1]] + 1
        obs[(t+1):TT, coord.tmp[i,2]] <- obs[(t+1):TT, coord.tmp[i,2]] + 1 
        x.tmp <- cbind(y[1:(TT-t), coord.tmp[i,1]], y[(t+1):TT, coord.tmp[i,2]])
        D.pair <- D.pair + nrow(x.tmp)
      }
    }
  }
  # Time-lag 0, Space-lag > 0
  for (s in s.len) {
    coord.tmp <- s.comb.0[s.comb.0[,3]==s,]
    for (i in 1:nrow(coord.tmp)) {
      obs[,coord.tmp[i,1]] <- obs[,coord.tmp[i,1]] + 1
      obs[,coord.tmp[i,2]] <- obs[,coord.tmp[i,2]] + 1
      x.tmp <- cbind(y[,coord.tmp[i,1]], y[,coord.tmp[i,2]])
      D.pair <- D.pair + nrow(x.tmp)
    }
  }
  # By symmetric first row of remedy stands for the number of time of marginal for 1 and T, 
  # (if applicable) second row for 2 and T-1 and so on.
  remedy <- matrix(NA, nrow=t.lag, ncol=S)
  for (i in 1:t.lag) {
    remedy[i,] <- obs[t.lag+1,]-obs[i,]
  }
  D.pair <- D.pair + sum(remedy)
  list(D.pair=D.pair, use.obs=obs, remedy=remedy)
}


# Calculate -logPL and D.pair with given distance matrix, lags and remedy matrix
pl <- function(theta, y, S.dist, t.lag=1, remedy, ts=T, mu_zero=F) {
  if(mu_zero){ # mu is fixed at zero
    phi <- theta[1]; rho <- theta[2]; sigma2 <- theta[3]; mu <- 0
  }else{
    if(ts){ # two-step estimation for computational efficiency; first mean function, then covariance function
      phi <- theta[1]; rho <- theta[2]; sigma2 <- theta[3]; mu <- mean(y)
    }else{
      phi <- theta[1]; rho <- theta[2]; sigma2 <- theta[3]; mu <- theta[4]
    }
  }
  
  S <- ncol(y); TT <- nrow(y)
  logPL <- D.pair <- 0

  Sigma <- sigma2*Matern(S.dist, rho=rho)
  y.var <- sigma2/(1-phi^2)
  y <- y - mu # mean function
  
  c1 <- log(2*pi); c2 <- c1/2
  # Time lag 1 to max t lag
  for (t in 1:t.lag) {
  	# Space-lag 0
  	S.tmp <- matrix(y.var*c(1, phi^t, phi^t, 1), ncol=2)
  	for (j in 1:S){
    	x.tmp <- cbind(y[1:max(TT-t,1),j], y[min(TT,(t+1)):TT,j])
    	logPL <- logPL + sum(dmvnorm(x=x.tmp, sigma=S.tmp, log=TRUE)+c1)
  	}
  	# Space-lag > 0
  	for (s in s.len) {
		  coord.tmp <- s.comb[s.comb[,3]==s,]
		  cov.tmp <- phi^t/(1-phi^2)*Sigma[coord.tmp[1,1], coord.tmp[1,2]]
		  S.tmp <- matrix(c(y.var, cov.tmp, cov.tmp, y.var), ncol=2)
		
  		for (i in 1:nrow(coord.tmp)) {
    		x.tmp <- cbind(y[1:max(TT-t,1), coord.tmp[i,1]], y[min(TT,(t+1)):TT, coord.tmp[i,2]])
    		logPL <- logPL + sum(dmvnorm(x=x.tmp, sigma=S.tmp, log=TRUE)+c1)
  		}
  	}
  }
  # Time-lag 0, Space-lag > 0
  for (s in s.len) {
	  coord.tmp <- s.comb.0[s.comb.0[,3]==s,]
	  cov.tmp <- Sigma[coord.tmp[1,1], coord.tmp[1,2]]/(1-phi^2)
	  S.tmp <- matrix(c(y.var, cov.tmp, cov.tmp, y.var), ncol=2)
	
  	for (i in 1:nrow(coord.tmp)) {
    	x.tmp <- cbind(y[,coord.tmp[i,1]], y[,coord.tmp[i,2]])
    	logPL <- logPL + sum(dmvnorm(x=x.tmp, sigma=S.tmp, log=TRUE)+c1)
  	}
  }
  # Correction for edge effect via marginal likelihood
  sd.remedy <- sqrt(y.var)
  for (i in 1:t.lag) {
    logPL <- logPL + sum(dnorm(rep(y[i,], c(remedy[i,])), sd=sd.remedy, log=TRUE)+c2) + # remedy at start of segment
    				 sum(dnorm(rep(y[(TT-i+1),], c(remedy[i,])), sd=sd.remedy, log=TRUE)+c2)	# remedy at end of segment
  }
  return(-logPL)
}


# Cost function for a segment
cost <- function(y.full, start, end, S.dist, t.lag=1, Ck, min.length, remedy, ini, ts=T, mu_zero=F) {
  S <- ncol(y.full)
  num_para <- length(ini)
  if (end - start + 1 < min.length) {
    return(Inf)
  } else {
    y.seg <- y.full[start:end,]
    if(ts | mu_zero){ # whether to use two-stage estimator (mean and covariance) for computational efficiency (or have known mu)
      ini <- ini[1:3]; lb <- c(-0.7,0.1,1e-3); ub <- c(0.7,3,5)
    }else{
      lb <- c(-0.7,0.1,1e-3,-1); ub <- c(0.7,3,5,1)
    }
    pmle <- optim(ini, pl, y=y.seg, S.dist=S.dist, t.lag=t.lag, remedy=remedy, ts=ts, mu_zero=mu_zero,
    		          method="L-BFGS-B", lower=lb, upper=ub)
    lik <- pmle$value
    pen <- Ck*(num_para/2*(log(end-start+1) + log(S)) + log(end-start+1))
  }
  lik + pen
}


# PELT function
pelt <- function(y.full, S.dist, Ck, remedy, K, t.lag=1, min.length, ini, res=1, ts=T, mu_zero=F) {
  TT <- nrow(y.full)
  S <- ncol(y.full)
  
  # For initializing the objective function F, the cp set cp[0] and the candidate set R[0]
  f <- c(0, rep(Inf,TT))
  cp <- 0
  R <- list(0)
  
  # For minimum segment length initialization
  for (i in 1:(2*min.length-1)) {
  	f[i+1] <- cost(y.full=y.full, start=1, end=i, S.dist=S.dist, t.lag=t.lag,
  	               Ck=Ck, min.length=min.length, remedy=remedy, ini=ini, ts=ts, mu_zero=mu_zero)
  	cpi <- list(c(0,i))
  	cp <- c(cp, cpi)
  	R.next <- list(c(0, i-min.length+1))
  	R <- c(R, R.next)
  	cat(i,"\n")
  }
  # Main function
  for(i in (2*min.length):TT) {
    possible.F <- rep(Inf, length(R[[i]]))
    # Further speed up computation by only considering every res time candidate
    # By default res=1 and thus no speed up
    if(i%%res==0){
      cnt <- 1
      for(j in R[[i]]) {
        possible.F[cnt] <- f[j+1] + cost(y.full=y.full, start=j+1, end=i, S.dist=S.dist,
                                         t.lag=t.lag, Ck=Ck, min.length=min.length,
                                         remedy=remedy, ini=ini, ts=ts, mu_zero=mu_zero)
        cnt <- cnt+1
      }
      f[i+1] <- min(possible.F)
      tau1 <- R[[i]][possible.F==f[i+1]][1]
    }
    cpi <- list(c(cp[[tau1+1]], i))
    cp <- c(cp, cpi)
    cat(i, "est cp:", cp[[i+1]], " candidate:", R[[i]], "\n")	# The current "best" solution by considering the first i data-points
    
  	# Pruning step in PELT
  	R.next <- c(R[[i]][possible.F+K<=f[i+1]], i-min.length+1)
  	R.next <- R.next[R.next%%res==0]
  	
  	R <- c(R, list(R.next))
  }
  cp[[TT+1]] 						# Final output is the change-point position by considering the whole dataset.
}


# Given change point location, calculate CLMDL (cost function)
total_cost <- function(cp_vec, y.full, S.dist, t.lag=1, Ck, remedy, min.length, ini, ts=T, mu_zero=F) {
	ans <- 0
	len <- length(cp_vec)
	for (i in 1:(len-1)) {
		ans <- ans + cost(y.full, start=cp_vec[i]+1, end=cp_vec[i+1], 
				              S.dist, t.lag, Ck, min.length, remedy, ini, ts, mu_zero)
	}
	ans
}


# CI function
Q.grid <- function(lambda.hat, TT, S.dist, y.full, theta1, theta2, t.lag, remedy, q.bound, ts=T, mu_zero=F) {
  q.range <- -q.bound:q.bound
  q.range <- q.range[q.range+lambda.hat >=1 & q.range+lambda.hat <= TT]
  cl <- 0*q.range
  q_min <- max(1, min(q.range) + lambda.hat - 5)
  q_max <- min(TT, max(q.range) + lambda.hat + 5)
  
  for (i in 1:(length(cl))) { #-1 to avoid single-row for theta2?
    cp_tmp <- lambda.hat + q.range[i]
    if (cp_tmp==TT) {
      cl[i] <- pl(theta=theta1, y=matrix(y.full[q_min:cp_tmp,],ncol=ncol(y.full)), S.dist, t.lag, remedy, ts, mu_zero)
    } else {
    cl[i] <- pl(theta=theta1, y=matrix(y.full[q_min:cp_tmp,],ncol=ncol(y.full)), S.dist, t.lag, remedy, ts, mu_zero) +
             pl(theta=theta2, y=matrix(y.full[(cp_tmp+1):q_max,],ncol=ncol(y.full)), S.dist, 
                t.lag,# ifelse(nrow(y.full)<=1,0,t.lag), 
                remedy, ts, mu_zero)
    }
  }
  q.range[cl==min(cl)]
}

