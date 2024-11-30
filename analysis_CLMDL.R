rm(list=ls())
source("Util_Functions.R")
library(mvtnorm)
library(foreach)
library(doParallel)

set.seed(1234)
cl<-makeCluster(detectCores()-1)
registerDoParallel(cl)
timStart<-proc.time()

# ### Read-in the simulated spatio-temporal process
# # data_name_set <- c('y0','y1','y2')
# # data_name <- 'y0' # no change-point
# data_name <- 'y1' # one change-point
# # data_name <- 'y2' # two change-point
# 
# 
# data <- readRDS(paste0(data_name, '.RDS'))
# y <- data$y # T times S spatio-temporal data
# S.dist <- data$S.dist # distance matrix of the S spatial locations
Q<-100
smalSamp<-F
savResults<-list()
foreach(q=1:Q,samp=howSamp,.inorder=F,.packages = c('mvtnorm')) %dopar%  {
  
print(paste('iter',q))
# set.seed(500)
samp<-samp
source('spatial_ar1_data.R',local=T)
y<-result$data
S.dist<-result$distance_matrix
S <- ncol(y) # spatial dimension
TT <- nrow(y) # time dimension

### Parameters for implementing CLDML
# Set spatial and temporal lag in computing CL
s.lag <- t.lag <- 1

# Collect all pairs of observations within s.lag distance (will be used in PL)
s.comb <- s.dist.pair(S.dist, s.lag) # time.lag > 0
s.comb.0 <- s.dist.pair.0(S.dist, s.lag) # time.lag = 0
s.len <- unique(s.comb[,3])
s.len <- s.len[s.len>0] # all unique spatial distance

# Calculate average number of times an obs used in CL for Ck
pair_stat <- D.cal(y, S.dist, s.lag, t.lag)
Ck <- mean(pair_stat$use.obs)
remedy <- pair_stat$remedy # number of marginal likelihood needed for correcting edge effect

# initial parameter estimation based on entire data using the 4-parameter spatial auto-regressive model
ini <- optim(c(0,1,1,0), pl, y=y, S.dist=S.dist, t.lag=t.lag, remedy=remedy,
             method="L-BFGS-B", lower=c(-0.7,0.1,1e-3,-1), upper=c(0.7,3,5,1), ts=F)$par

# Calculate K for pruning step
p.min <- p.max <- length(ini)
K <- floor(Ck)*(floor(log(S*TT))*(p.min/2 - p.max) + (2 + p.max)*log(2) - floor(log(TT)))

######### Main function ##########
### Run PELT for minimizing CLMDL (this is the default algorithm in the paper.)
t0 <- proc.time()
pelt_result <- pelt(y.full=y, S.dist=S.dist, Ck=Ck, remedy=remedy, K=K,
                    t.lag=t.lag, min.length=0.1*TT, ini=ini, res=1)
t_compute <- proc.time() - t0

### One can run PELT with res=2 to further speed up computation
# (res=2 means we assume that cp can only happen on t such that t%%2==0, i.e. t is even.)
# Therefore the resolution of this algorithm is 2.
# t0 <- proc.time()
# pelt_result2 <- pelt(y.full=y, S.dist=S.dist, Ck=Ck, remedy=remedy, K=K,
#                      t.lag=t.lag, min.length=0.1*TT, ini=ini, res=2)
# t_compute2 <- proc.time() - t0

### Subsequent analysis after CP detection
num_cp <- length(pelt_result)-2
ci_result <- c()
if(num_cp>0){
  for(cp_index in 1:num_cp){
    # Parameter estimation on the pre-change segment
    y.seg1 <- y[(pelt_result[cp_index]+1):pelt_result[cp_index+1],]
    pmle1 <- optim(ini, pl, y=y.seg1, S.dist=S.dist, t.lag=t.lag, remedy=remedy,
                   method="L-BFGS-B", lower=c(-0.7,0.1,1e-3,-1), upper=c(0.7,3,5,1), ts=F)
    # Parameter estimation on the post-change segment
    y.seg2 <- y[(pelt_result[cp_index+1]+1):pelt_result[cp_index+2],]
    pmle2 <- optim(ini, pl, y=y.seg2, S.dist=S.dist, t.lag=t.lag, remedy=remedy,
                   method="L-BFGS-B", lower=c(-0.7,0.1,1e-3,-1), upper=c(0.7,3,5,1), ts=F)
    
    seg1_tmp <- nrow(y.seg1)
    seg2_tmp <- nrow(y.seg2)
    TT_tmp <- seg1_tmp + seg2_tmp
    
    # CI construction
    B <- 100
    theta1_est <- pmle1$par
    theta2_est <- pmle2$par
    q_emp <- c()
    t0 <- proc.time()
    for(b in 1:B){
      set.seed(b)
      print(b)
      yB <- rbind(sim.y(theta=theta1_est, S.dist=S.dist, TT=seg1_tmp, T.burn=100),
                  sim.y(theta=theta2_est, S.dist=S.dist, TT=seg2_tmp, T.burn=100))
      #Q.grid is where the failure comes, something in pl()
      q_emp[b] <- Q.grid(lambda.hat=seg1_tmp, TT=TT_tmp, S.dist=S.dist, y.full=yB,
                         theta1=theta1_est, theta2=theta2_est,
                         t.lag=t.lag, remedy=remedy, q.bound=round(0.2*TT_tmp)+1, ts=F,mu_zero=T)
    }
    t_compute3 <- proc.time() - t0
    
    ci <- quantile(q_emp, probs=c(0.025,0.975))
    ci_result <- rbind(ci_result, pelt_result[cp_index+1]+ci)
  }
}
num_cp
# save.image(paste0('Result_', data_name, '.RData'))
if (is.null(ci_result)) {
  # savResults[[q]]<-
    return(NA)
} else {
  # savResults[[q]]<-
    return(ci_result)
}
}
timEnd<-proc.time()
save(file='First100CLMDLMed.RData',list=c('savResults','timStart','timEnd'))
