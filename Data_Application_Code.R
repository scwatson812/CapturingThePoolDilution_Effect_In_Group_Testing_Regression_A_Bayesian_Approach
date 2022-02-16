############################
#####  Load Libraries  #####
############################
library(truncnorm)
library(MCMCpack)
library(BayesLogit)
library(R.utils)
library(coga)

#######################
#####  Functions  #####
#######################

###Binary Link Function for Disease Status
g.link <-function(x){
  res = exp(x)/(1 + exp(x))
  return(res)
}

###Biomarker Distribution Function
f <-function(zeta,theta0,theta1,Y){ 
  if(Y==0){
    res = dgamma(zeta,theta0[1],theta0[2])
  }
  if(Y==1){
    res = dgamma(zeta,theta1[1],theta1[2])
  }
  return(res)
}

###Sample Biomarker Measurement Error Variance
sigma.w.samp<-function(){
  res = rinvgamma(1,J/2 + a.sigma.w,sum((W.tilde.g - W)^2)/2 + b.sigma.w)
  return(res)
}

###Sample the Individual Biomarker Levels
zeta.tilde.samp <-function(zeta.ind){
  zeta.tilde.i.prop = rnorm(1,zeta.tilde.g[zeta.ind],sd.zeta.tilde[zeta.ind])
  zeta.tilde.prop = zeta.tilde.g
  zeta.tilde.prop[zeta.ind] = zeta.tilde.i.prop
  W.tilde.prop = W.tilde.g
  for(z.i in pools.containing.ind[[zeta.ind]]){
    W.tilde.prop[z.i] = mean(zeta.tilde.prop[mathcalP[[z.i]]])
  }
  den.g = sum(dnorm(W[pools.containing.ind[[zeta.ind]]], mean = W.tilde.g[pools.containing.ind[[zeta.ind]]],sqrt(sigma.w.g),log = TRUE))
  den.g = den.g + log(f(zeta.tilde.g[zeta.ind],theta.0.g,theta.1.g,Y.tilde.g[zeta.ind]))
  den.g = den.g + log(dnorm(zeta.tilde.i.prop,zeta.tilde.g[zeta.ind],sd.zeta.tilde[zeta.ind]))
  den.prop = sum(dnorm(W[pools.containing.ind[[zeta.ind]]], W.tilde.prop[pools.containing.ind[[zeta.ind]]],sqrt(sigma.w.g),log = TRUE))
  den.prop = den.prop + log(f(zeta.tilde.prop[zeta.ind],theta.0.g,theta.1.g,Y.tilde.g[zeta.ind]))
  den.prop = den.prop + log(dnorm(zeta.tilde.g[zeta.ind],zeta.tilde.i.prop, sd.zeta.tilde[zeta.ind]))
  A = min(1,exp(den.prop-den.g))
  if(is.na(A)){
    A = 0
  }
  r = rbinom(1,1,A)
  zeta.tilde.i.res = r*zeta.tilde.i.prop + (1-r)*zeta.tilde.g[zeta.ind]
  W.tilde.res = r*W.tilde.prop + (1-r)*W.tilde.g
  if(is.na(r)){
    print(den.g)
    print(den.prop)
    print(zeta.tilde.i.prop)
  }
  return(list(zeta.out = zeta.tilde.i.res,r.out = r,W.tilde.out = W.tilde.res))
}

###Sample the Latent true individual statuses
Y.tilde.samp <-function(){
  p0 = f(zeta.tilde.g,theta.0.g,theta.1.g,0)
  p0 = p0*(1-g.link(X%*%beta.g))
  p1 = f(zeta.tilde.g,theta.0.g,theta.1.g,1)
  p1 = p1*g.link(X%*%beta.g)
  y.out = rbinom(N,1,(p1/(p0+p1)))
  return(y.out)
}

###Sample the Poyla Gamma random variable for the population regression model
omega.samp<-function(){
  res = rpg(N,1,X%*%beta.g)
  return(res)
}

###Sample the population regression coefficients
beta.samp<-function(){
  sigma = solve(t(X)%*%Omega.g%*%X + Phi.inv.g)
  mean = sigma%*%t(X)%*%h.g
  res = mvrnorm(1,mean,sigma)
  return(res)
}

###Sample the biomarker distribution parameters for the negative biomarker distribution
theta.0.samp <-function(theta.i){
  theta.prop <-rtruncnorm(1,theta.0.g[theta.i],sd =sd.theta.0[theta.i],a = theta.0.prior.lb[theta.i], b = theta.0.prior.ub[theta.i])
  theta.0.prop <- theta.0.g
  theta.0.prop[theta.i] = theta.prop
  den.g = sum(log(f(zeta.tilde.g[Y.tilde.g ==0],theta.0.g,theta.1.g,0)))
  #den.g = den.g + log(dtruncnorm(theta.0.g[theta.i],mean = theta.0.prior.mean[theta.i], sd = theta.0.prior.sd[theta.i],a = theta.0.prior.lb[theta.i], b = theta.0.prior.ub[theta.i]))
  den.g = den.g + log(dgamma(theta.0.g[theta.i],theta.0.prior.alpha[theta.i], theta.0.prior.beta[theta.i]))
  den.g = den.g + log(dtruncnorm(theta.prop,theta.0.g[theta.i],sd.theta.0[theta.i],a = theta.0.prior.lb[theta.i], b = theta.0.prior.ub[theta.i]))
  den.p = sum(log(f(zeta.tilde.g[Y.tilde.g ==0],theta.0.prop,theta.1.g,0)))
  #den.p = den.p + log(dtruncnorm(theta.0.prop[theta.i],mean = theta.0.prior.mean[theta.i], sd = theta.0.prior.sd[theta.i],a = theta.0.prior.lb[theta.i], b = theta.0.prior.ub[theta.i]))
  den.p = den.p + log(dgamma(theta.0.prop[theta.i],theta.0.prior.alpha[theta.i],theta.0.prior.beta[theta.i]))
  den.p = den.p + log(dtruncnorm(theta.0.g[theta.i],theta.0.prop[theta.i],sd.theta.0[theta.i],a = theta.0.prior.lb[theta.i], b = theta.0.prior.ub[theta.i]))
  A = min(1,exp(den.p-den.g))
  if(is.na(A)){
    A = 0
  }
  r = rbinom(1,1,A)
  theta.out = r*theta.prop + (1-r)*theta.0.g[theta.i]
  return(list(theta.out = theta.out, r.out = r))
}

###Sample the biomarker distribution parameters for the positive biomarker distribution
theta.1.samp <-function(theta.i){
  theta.prop <-rtruncnorm(1,theta.1.g[theta.i],sd =sd.theta.1[theta.i],a = theta.1.prior.lb[theta.i], b = theta.1.prior.ub[theta.i])
  theta.1.prop <- theta.1.g
  theta.1.prop[theta.i] = theta.prop
  den.g = sum(log(f(zeta.tilde.g[Y.tilde.g ==1],theta.0.g,theta.1.g,1)))
  den.g = den.g + log(dgamma(theta.1.g[theta.i],theta.1.prior.alpha[theta.i],theta.1.prior.beta[theta.i]))
  den.g = den.g + log(dtruncnorm(theta.prop,theta.1.g[theta.i],sd.theta.1[theta.i],a = theta.1.prior.lb[theta.i], b = theta.1.prior.ub[theta.i]))
  den.p = sum(log(f(zeta.tilde.g[Y.tilde.g ==1],theta.0.g,theta.1.prop,1)))
  den.p = den.p + log(dgamma(theta.1.prop[theta.i],theta.1.prior.alpha[theta.i],theta.1.prior.beta[theta.i]))
  den.p = den.p + log(dtruncnorm(theta.1.g[theta.i],theta.1.prop[theta.i],sd.theta.1[theta.i],a = theta.1.prior.lb[theta.i], b = theta.1.prior.ub[theta.i]))
  A = min(1,exp(den.p-den.g))
  if(is.na(A)){
    A = 0
  }
  r = rbinom(1,1,A)
  theta.out = r*theta.prop + (1-r)*theta.1.g[theta.i]
  return(list(theta.out = theta.out, r.out = r))
}

###Helper function to calculate the probability that a negative pool tests negative
w.con.den.S_p<-function(x,w,c.s,theta.0.mean,theta.0.sd,sig.w){
  res = pcoga(w,shape = rep(theta.0.mean,c.s),rate = rep(theta.0.sd,c.s)*c.s)
  res = res*dnorm(x,0,sqrt(sig.w))
  return(res)
}

###Helper function to calculate the probability that a positive pool tests positive
w.con.den.S_e<-function(x,w,c.s,theta.0.mean,theta.0.sd, theta.1.mean,theta.1.sd,sig.w){
  res = pcoga(w,shape = c(rep(theta.0.mean,c.s-1),theta.1.mean), rate = c(rep(theta.0.sd,c.s-1),theta.1.sd)*c.s)
  res = res*dnorm(x,0,sqrt(sig.w))
  return(res)
}

###probability that a negative pool tests negative
S_p <- function(w,c.s,theta.0.mean,theta.0.sd,sig.w){
  res = integrate(w.con.den.S_p,lower = -Inf, upper = Inf,c.s = c.s,w = w,theta.0.mean=theta.0.mean,theta.0.sd = theta.0.sd, sig.w = sig.w )
  return(res$value)
}

###Sensitivity of a pool of size c with threshold of w, conditional on the pool containing 1 positive
S_e1 <- function(w,c.s,theta.0.mean,theta.0.sd,theta.1.mean,theta.1.sd,sig.w){
  res = 1-integrate(w.con.den.S_e,lower = -Inf, upper = Inf,c.s = c.s,w = w,theta.0.mean=theta.0.mean,theta.0.sd = theta.0.sd,theta.1.mean = theta.1.mean, theta.1.sd = theta.1.sd, sig.w = sig.w )$value
  return(res)
}

###Youden.index to be optimized to find optimal threshold
Youden.index <-function(w,fixed.par){
  c.s = fixed.par[1]
  theta.0.mean = fixed.par[2]
  theta.0.sd = fixed.par[3]
  theta.1.mean = fixed.par[4]
  theta.1.sd = fixed.par[5]
  sig.w = fixed.par[6]
  res = S_p(w,c.s,theta.0.mean,theta.0.sd,sig.w)
  res = res + S_e1(w,c.s,theta.0.mean,theta.0.sd,theta.1.mean,theta.1.sd,sig.w)
  res = res -1
  return(-res)
}

##############################
#####  Read in the Data  #####
##############################
###Read in tbe Individual level data (unique ids and covariates)
ind.df = read.csv(file.choose())
###Read in the Pool level data (pool responses and unique ids of individuals in each pool)
pool.df = read.csv(file.choose())

N = dim(ind.df)[1]
X = cbind(1,ind.df$Covariate.1,ind.df$Covariate.2)
Q = dim(X)[2]

###Create Pools for Array Testing  #####
pool.length.list  = pool.df$Pool.Size
J = dim(pool.df)[1]
mathcalP = list()
for(j in 1:J){
  mathcalP[[j]] = pool.df[j,4:(dim(pool.df)[2])][!is.na(pool.df[j,4:(dim(pool.df)[2])])]
}

pools.containing.ind = list()
for(i in 1:N){
  pool.list = c()
  for(j in 1:J){
    if(is.element(i,mathcalP[[j]])){
      pool.list = c(pool.list,j)
    }
  }
  pools.containing.ind[[i]] = pool.list
  print(i)
}

W = pool.df$Observed.Biomarker
plot(W)

#####################################
#####  Specify Hyperparameters  #####
#####################################
###Estimate of the mean and variance of the negative biomarker distribution
mu.0 = 5
var.0 = 10
###Estimate of the mean and variance of the positive biomarker distribution
mu.1 = 19
var.1 = 20 

theta.0.prior.mean = c(mu.0^2/var.0,mu.0/var.0)
theta.0.prior.var= 0.1*theta.0.prior.mean
theta.0.prior.alpha = theta.0.prior.mean^2/theta.0.prior.var
theta.0.prior.beta = theta.0.prior.mean/theta.0.prior.var
theta.1.prior.mean = c(mu.1^2/var.1,mu.1/var.1)
theta.1.prior.var= 0.1*theta.1.prior.mean
theta.1.prior.alpha = theta.1.prior.mean^2/theta.1.prior.var
theta.1.prior.beta = theta.1.prior.mean/theta.1.prior.var
mean(rgamma(1000,theta.0.prior.alpha[1],theta.0.prior.beta[1]))
var(rgamma(1000,theta.0.prior.alpha[1],theta.0.prior.beta[1]))
mean(rgamma(1000,theta.0.prior.alpha[2],theta.0.prior.beta[2]))
var(rgamma(1000,theta.0.prior.alpha[2],theta.0.prior.beta[2]))
mean(rgamma(1000,theta.1.prior.alpha[1],theta.1.prior.beta[1]))
var(rgamma(1000,theta.1.prior.alpha[1],theta.1.prior.beta[1]))
mean(rgamma(1000,theta.1.prior.alpha[2],theta.1.prior.beta[2]))
var(rgamma(1000,theta.1.prior.alpha[2],theta.1.prior.beta[2]))
theta.0.prior.lb = c(0,0)
theta.0.prior.ub = c(Inf,Inf)
theta.1.prior.lb = c(0,0)
theta.1.prior.ub = c(Inf,Inf)

###Specify hyperparmeters for the erro variance
a.sigma.w = 3
b.sigma.w = 3

######################
#####  The MCMC  #####
######################
###Total Number of MCMC Interations
G = 3000
###Burn-in MCMC interations
burn = 2000
###Initialize the MCMC
zeta.tilde.g = rep(0.5,N)
W.tilde.g = rep(NA,J)
for(j in 1:J){
  W.tilde.g[j] = mean(zeta.tilde.g[mathcalP[[j]]])
}
sigma.w.g = 0.005^2
Y.tilde.g = rep(0,N)
theta.0.g = c(2.5,0.5)
theta.1.g = c(20,1)
beta.g = rep(0,Q)
Phi.inv.g = solve(100*diag(Q))
omega.g = rep(1,N)
Omega.g = diag(omega.g)
h.g = (Y.tilde.g -1/2)

###Tuning
acc.zeta.tilde = rep(0,N)
sd.zeta.tilde = rep(0.1,N)
acc.theta.0 = rep(0,length(theta.0.g))
sd.theta.0 = rep(0.5,length(theta.0.g))
acc.theta.1 = rep(0,length(theta.1.g))
sd.theta.1 = rep(0.5,length(theta.1.g))

###Storage 
W.tilde.array = array(NA,c(J,G))
sigma.w.array = array(NA,c(1,G))
zeta.tilde.array = array(NA,c(N,G))
Y.tilde.array = array(NA,c(N,G))
theta.0.array = array(NA,c(length(theta.0.g),G))
theta.1.array = array(NA,c(length(theta.1.g),G))
beta.array = array(NA,c(Q,G))
omega.array = array(NA,c(N,G))

for(g in 1:G){
  ###Sample zeta.tilde
  for(i in 1:N){
    zeta.vals = zeta.tilde.samp(i)
    zeta.tilde.g[i] = zeta.vals$zeta.out
    acc.zeta.tilde[i] = acc.zeta.tilde[i] + zeta.vals$r.out
    W.tilde.g = zeta.vals$W.tilde.out
    if(is.na(zeta.tilde.g[i])){
      break
    }
  }
  if(is.na(zeta.tilde.g[i])){
    break
  }
  ###Sample sigma.w
  sigma.w.g = sigma.w.samp()
  ###Sample Y.tilde
  Y.tilde.g = Y.tilde.samp()
  h.g = (Y.tilde.g -1/2)
  ###Sample omega
  omega.g = omega.samp()
  Omega.g = diag(omega.g)
  h.g = (Y.tilde.g -1/2)
  ###Sample beta
  beta.g = beta.samp()
  ###Sample the biomarker distribution parameters
  for(i in 1:length(theta.0.g)){
    theta.vals = theta.0.samp(i)
    theta.0.g[i] = theta.vals$theta.out
    acc.theta.0[i] = acc.theta.0[i] + theta.vals$r.out
  }
  for(i in 1:length(theta.1.g)){
    theta.vals = theta.1.samp(i)
    theta.1.g[i] = theta.vals$theta.out
    acc.theta.1[i] = acc.theta.1[i] + theta.vals$r.out
  }
  W.tilde.array[,g] = W.tilde.g
  sigma.w.array[,g] = sigma.w.g
  zeta.tilde.array[,g] = zeta.tilde.g
  Y.tilde.array[,g] = Y.tilde.g
  theta.0.array[,g] = theta.0.g
  theta.1.array[,g] = theta.1.g
  beta.array[,g] = beta.g
  omega.array[,g] = omega.g
  if(g%%100 == 0 & g < burn){
    acc.zeta.tilde = acc.zeta.tilde/100
    ind = acc.zeta.tilde < 0.3
    sd.zeta.tilde[ind] = sd.zeta.tilde[ind]*0.9
    ind = acc.zeta.tilde > 0.7
    sd.zeta.tilde[ind] = sd.zeta.tilde[ind]*1.1
    acc.zeta.tilde = rep(0,N)
    acc.theta.0 = acc.theta.0/100
    ind = acc.theta.0 < 0.3
    sd.theta.0[ind] = sd.theta.0[ind]*0.9
    ind = acc.theta.0 > 0.7
    sd.theta.0[ind] = sd.theta.0[ind]*1.1
    acc.theta.0 = rep(0,length(theta.0.g))
    acc.theta.1 = acc.theta.1/100
    ind = acc.theta.1 < 0.3
    sd.theta.1[ind] = sd.theta.1[ind]*0.9
    ind = acc.theta.1 > 0.7
    sd.theta.1[ind] = sd.theta.1[ind]*1.1
    acc.theta.1 = rep(0,length(theta.0.g))
  }
  print(g)
}

###Estimate the Optimal Classification Thresholds
theta.0.g = apply(theta.0.array[,burn:G],1,mean)
theta.1.g = apply(theta.1.array[,burn:G],1,mean)
sigma.w.g = mean(sigma.w.array[,burn:G])
pool.sizes = sort(unique(pool.length.list))
threshold.list = c()
for(c.s in pool.sizes){
  fixed.par.input = c(c.s = c.s,theta.0.mean = theta.0.g[1], theta.0.sd = theta.0.g[2], theta.1.mean = theta.1.g[1], theta.1.sd = theta.1.g[2], sig.w = sigma.w.g)
  pool.optim <- optimize(f = Youden.index, lower = 0, upper = 10, fixed.par = fixed.par.input)
  threshold.list = c(threshold.list,pool.optim$minimum)
  plot.seq = seq(0.0001,30,length.out = 1000)
  plot(plot.seq,dcoga(plot.seq, rep(theta.0.g[1],c.s),rep(theta.0.g[2],c.s)*c.s),type = 'l')
  points(plot.seq,dcoga(plot.seq,c(rep(theta.0.g[1],c.s-1),theta.1.g[1]), c(rep(theta.0.g[2],c.s-1),theta.1.g[2])*c.s),type = 'l',col = 'red')
  abline(v = pool.optim$minimum)
}

########################################################
#####  Calculate Point Estimators from First MCMC  #####
########################################################

beta.hat = apply(beta.array[,burn:G],1,mean)
beta.sd = apply(beta.array[,burn:G],1,sd)
beta.lower = apply(beta.array[,burn:G],1,quantile,prob = 0.025)
beta.upper = apply(beta.array[,burn:G],1,quantile,prob = 0.975)

sigma.w.hat = mean(sigma.w.array[,burn:G])
sigma.w.sd = sd(sigma.w.array[,burn:G])
sigma.w.lower = quantile(sigma.w.array[,burn:G],prob = 0.025)
sigma.w.upper = quantile(sigma.w.array[,burn:G],prob = 0.975)

zeta.tilde.hat = apply(zeta.tilde.array[,burn:G],1,mean)
zeta.tilde.sd = apply(zeta.tilde.array[,burn:G],1,sd)
zeta.tilde.lower = apply(zeta.tilde.array[,burn:G],1,quantile,prob = 0.025)
zeta.tilde.upper = apply(zeta.tilde.array[,burn:G],1,quantile,prob = 0.975)

W.tilde.hat = apply(W.tilde.array[,burn:G],1,mean)
W.tilde.sd = apply(W.tilde.array[,burn:G],1,sd)
W.tilde.lower = apply(W.tilde.array[,burn:G],1,quantile,prob = 0.025)
W.tilde.upper = apply(W.tilde.array[,burn:G],1,quantile,prob = 0.975)

Y.tilde.hat = apply(Y.tilde.array[,burn:G],1,mean)
Y.tilde.sd = apply(Y.tilde.array[,burn:G],1,sd)
Y.tilde.lower = apply(Y.tilde.array[,burn:G],1,quantile,prob = 0.025)
Y.tilde.upper = apply(Y.tilde.array[,burn:G],1,quantile,prob = 0.975)

theta.0.hat = apply(theta.0.array[,burn:G],1,mean)
theta.0.sd = apply(theta.0.array[,burn:G],1,sd)
theta.0.lower = apply(theta.0.array[,burn:G],1,quantile,prob = 0.025)
theta.0.upper = apply(theta.0.array[,burn:G],1,quantile,prob = 0.975)

theta.1.hat = apply(theta.1.array[,burn:G],1,mean)
theta.1.sd = apply(theta.1.array[,burn:G],1,sd)
theta.1.lower = apply(theta.1.array[,burn:G],1,quantile,prob = 0.025)
theta.1.upper = apply(theta.1.array[,burn:G],1,quantile,prob = 0.975)

#############################
#####  Trace Plot MCMC  #####
#############################
plot(beta.array[1,])
plot(beta.array[2,])
plot(beta.array[3,])

plot(sigma.w.array[1,])

plot(theta.0.array[1,])
plot(theta.0.array[2,])
plot(theta.0.array[1,]/theta.0.array[2,])
plot(theta.0.array[1,]/theta.0.array[2,]^2)
plot(theta.1.array[1,])
plot(theta.1.array[2,])
plot(theta.1.array[1,]/theta.1.array[2,])
plot(theta.1.array[1,]/theta.1.array[2,]^2)

plot(zeta.tilde.array[1,])
plot(zeta.tilde.array[100,])
plot(zeta.tilde.array[200,])
plot(zeta.tilde.array[300,])

plot(W.tilde.array[1,])
plot(W.tilde.array[2,])
plot(W.tilde.array[3,])
plot(W.tilde.array[4,])

###################################
#####  Identify Positive Pools ####
###################################

###Classify Pools
pool.status = rep(N,J)
ct = 1
for(c.s in pool.sizes){
  pool.status[pool.length.list == c.s] = W[pool.length.list == c.s] > threshold.list[ct]
  ct = ct + 1
}

positive.pools = which(pool.status == 1)
write.csv(positive.pools,"Positive_Pools.csv",row.names = FALSE)

beta.out = rbind(beta.hat,beta.sd,beta.lower,beta.upper)
write.csv(beta.out,'beta.csv',row.names = FALSE)

sigma.w.out = rbind(sigma.w.hat,sigma.w.sd,sigma.w.lower, sigma.w.upper)
write.csv(sigma.w.out,'sigmaw.csv',row.names = FALSE)

zeta.tilde.out = rbind(zeta.tilde.hat,zeta.tilde.sd,zeta.tilde.lower, zeta.tilde.upper)
write.csv(zeta.tilde.out,'Ztilde.csv',row.names = FALSE)

W.tilde.out = rbind(W.tilde.hat,W.tilde.sd,W.tilde.lower,W.tilde.upper)
write.csv(W.tilde.out,'Wtilde1.csv',row.names = FALSE)

Y.tilde.out = rbind(Y.tilde.hat,Y.tilde.sd,Y.tilde.lower,Y.tilde.upper)
write.table(Y.tilde.out,'Ytilde.csv',row.names = FALSE)

theta.out = rbind(c(theta.0.hat,theta.1.hat),c(theta.0.sd,theta.1.sd),c(theta.0.lower,theta.1.lower),c(theta.0.upper,theta.1.upper))
write.table(theta.out,'theta.csv',row.names = FALSE)

write.table(threshold.list,'threshold.csv',row.names = FALSE)

