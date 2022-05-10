####################################################################################################
###Codes for Group testing with introduction of senstivity $ specificity############################
###############i.e Group testing with errors in inspection##########################################
###############MLE of P.

#########code returns a single value of the estimate.
p.hat=function(n, k, p, sens, spec){
  pro= spec - ((1-p)^k)*(sens+spec-1)
  t=c()
  MLE= c()
  g=c();th=c()
  N=10000
  for (i in 1:N) {
    # if(i%%1000==0){print(i)} ####modulo to show wether the code is running & iterating
    t[i] = rnbinom(1, n, pro)+n
    g=n/t[i]
    th= ((spec-n/t[i])/(sens+spec-1))^(1/k)
    
    if(t[i]==n){
      MLE[i]= 1
    }
    else if(spec<=g){
      MLE[i]=1
    }
    else if(th>=1){
      MLE[i]=0
    }
    else{
      MLE[i] = 1-((spec-n/t[i])/(sens+spec-1))^(1/k)
    }
  }
  return(mean(MLE))
}
options(digits = 5);set.seed(12345); p.hat(1,5,0.3,0.95,0.95)

#############################################
##########Evaluating at different values of p to 6d.p)
####################
p=c(0.005,0.01,0.05,0.10,0.20,0.30)
Frap=rep(NA, length(p))

for (i in 1:length(p)){
  set.seed(12345); Frap[i]=p.hat(1,50,p[i],0.90,0.90)
  d=round(Frap, digits=6)
}

##################################################################
####################Plotting the MLE at different proportions#####
p=seq(0, 1,.05)
set.seed(12345);Frap=rep(NA, length(p))
set.seed(12345);Frap1=rep(NA, length(p))
set.seed(12345);Frap2=rep(NA, length(p))
set.seed(12345);Frap3=rep(NA, length(p))

for (i in 1:length(p)){
  Frap[i]=p.hat(15,20,p[i],0.99,0.99)
  Frap1[i]=p.hat(15,20,p[i],0.98,0.98)
  Frap2[i]=p.hat(15,20,p[i],0.95,0.95)
  Frap3[i]=p.hat(15,20,p[i],0.90,0.90)
}
####Plot for fixed k, varying n values
plot(p, Frap, type="o", col="darkblue", xlab="p", ylab="p.estimate", pch=21, lty=1,
     main= "n=15")
lines(p,Frap1, type="o", col="red4", pch=8, lty=1)
lines(p,Frap2, type="o", col="darkgreen", pch=15, lty=5)
lines(p, Frap3, type="o", pch=19, lty=6)
#legend("bottomright", legend=c("k=5","k=15","k=30", "k=50"), col=c('darkblue',"red4","darkgreen", "black" ), 
# lty=c(1,1,5,6), pch = c(21,8,15, 19), cex=0.8)
legend ("bottomright", legend=c(expression(paste(eta,"=",beta,"=",0.99)), expression(paste(eta,"=",beta,"=",0.98)),expression(paste(eta,"=",beta,"=",0.95)),
                                expression(paste(eta,"=",beta,"=",0.90))),col=c('darkblue',"red4","darkgreen", "black" )
        ,lty=c(1,1,5,6),  pch = c(21,8,15, 19))
cbind(d)


###################code for Bias and MLE
##########Mean Squared Error +Bias############
###############using for function#############
frank.Bias=function(toler, n, k, p, sens, spec){
  tstar=c();g=c();Bias.p=c();mse.p=c()
  t=c();pmft=c();phat=c()
  r=sens+spec-1
  theta=sens-r*(1-(1-p)^k)
  tolcheck=1-toler
  th=c()
  tstar=qnbinom(tolcheck,n, theta)
  t=seq(n, tstar+n, by=1)
  diffl=c(); BIAS=c(); MSE=c()
  #print(t)
  for (i in seq_along(t)) {
    pmft[i]=dnbinom((t[i]-n), n, theta)
    g[i]=n/t[i]
    th[i]= ((spec-n/t[i])/(sens+spec-1))^(1/k)
    if(g[i]==1){
      phat[i]=1
    }
    else if (sens<=g[i]){
      phat[i]=1
    }
    else if(th[i]>=1){
      phat[i]=0
    }
    else{
      phat[i]=1-((sens-n/t[i])/(r))^(1/k)
    }
    diffl[i]=abs(phat[i]-p)
    Bias.p[i]=diffl[i]*pmft[i]
    mse.p[i]=(phat[i]-p)^2 *pmft[i]
    BIAS=sum(Bias.p)*10^4
    MSE=sum(mse.p)*10^4
  }
  #return(output)
  #return(data.frame(t, phat, diffl,BIAS, pmft))
  return(BIAS)
}

frank.Bias (0.00001,20,5,0.003, 0.99,0.99)

p=c(0.005,0.01,0.05,0.10,0.20, 0.30,0.40)
pwbp=rep(NA, length(p))

#out=numeric(length(p))
for (i in 1:length(p)) {
  pwbp[i]=frank.Bias (0.00001,15,50,p[i], 0.99,0.99)
}

cbind(pwbp)

############################################################
#########R code for Mean Squared Error of the estimator#####
############################################################
MSE.RET=function(toler, n, k, p, sens, spec){
  tstar=c();g=c();Bias.p=c();mse.p=c()
  t=c();pmft=c();MLE=c()
  r=(spec^2)-(1-sens)^2
  prob= (((1-sens)^2)*(1-p)^k) + (spec^2)*(1-(1-p)^k)
  tolcheck=1-toler
  tstar=qnbinom(tolcheck,n, prob=prob)
  th=c()
  bias.p=c();mse.p=c(); diffl=c()
  t=seq(n, tstar+n, by=1)
  #print(t)
  for (i in seq_along(t)) {
    pmft[i]=dnbinom((t[i]-n), n, prob=prob)
    th[i]=(((spec^2)-n/t[i])/((spec^2)-(1-sens)^2))^(1/k)
    g[i]=n/t[i]
    
    if(t[i]==n){
      MLE[i]= 1
    }
    else if((spec^2)<=g[i]){
      MLE[i]=1
    }
    else if(th[i]>=1){
      MLE[i]=0
    }
    else{
      MLE[i] = 1-(((spec^2)-n/t[i])/((spec^2)-(1-sens)^2))^(1/k)
    }
    diffl[i]=abs(MLE[i]-p)
    Bias.p[i]=diffl[i]*pmft[i]
    mse.p[i]=(MLE[i]-p)^2 *pmft[i]
    BIAS=sum(Bias.p)*10^4
    MSE=sum(mse.p)*10^4
  }
  return(MSE)
}

MSE.RET(0.00001,1,5,0.3, 0.90,0.90)

######################################################################
#############Confidence Interval######################################
###########Simulating Wald Confidence Interval   #####################
########set.seed=12345
#This code computes 10, 000 Montecarlo confidence interval for each simulation, each MLE for fixed n and different t.
##modify to produce mean value of MLE, and the corresponding mean value of the Confidence interval and the distance between them.
###p.est=mean (out$MLE)
###ll=mean(out$ul)
###ul=mean(out$ll)
###dif=ul-ll
####return(cbind(LL=ll, UL=ul, length=dif)

wald.ci=function (n, k, p, gamma=0.05, sens,spec){
  con <- 1-gamma/2
  r=spec+sens-1
  theta=sens-r*(1-(1-p)^k)
  zalph <- qnorm(con,0,1)
  p.hat=c()
  t=c()
  N=10000
  ci=c()
  ul=c(); ll=c(); var.phat=c(); FisherInfo=c()
  th=c()
  out=data.frame()
  for (i in 1:N) {
    t[i]=rnbinom(1,n=n, theta)+n
    th[i]= ((spec-n/t[i])/(sens+spec-1))^(1/k)
    if (n==t[i]){
      p.hat[i]= 1
      ul[i]=1
      ll[i]=0
      var.phat[i]=0
    }
    else if (sens<= n/t[i]){
      p.hat[i]=1
      ul[i]=1
      ll[i]=0
      var.phat[i]=0
    }
    else if(th[i]>=1){
      p.hat[i]=0
      FisherInfo[i]=((r^2)*n*(k^2)*(1-p.hat[i])^(2*k-2))/(((sens-r*(1-p.hat[i])^k)^2)*(1-sens+r*(1-p.hat[i])^k))
      var.phat[i]=(FisherInfo[i])^-1
      ul[i]=p.hat[i] + zalph*sqrt(var.phat[i])
      ll[i]= p.hat[i] - zalph*sqrt(var.phat[i])
    }
    else{
      p.hat[i]=1-((sens-n/t[i])/(r))^(1/k)
      FisherInfo[i]=((r^2)*n*(k^2)*(1-p.hat[i])^(2*k-2))/(((sens-r*(1-p.hat[i])^k)^2)*(1-sens+r*(1-p.hat[i])^k))
      var.phat[i]=(FisherInfo[i])^-1
      ul[i]=p.hat[i] + zalph*sqrt(var.phat[i])
      ll[i]= p.hat[i] - zalph*sqrt(var.phat[i])
    }
  }
  out=as.data.frame(cbind(n=n, T=t, MLE=p.hat, ll=ll, ul=ul, zalph, var.phat))
  ll.conf=mean(out$ll)
  ul.conf=mean(out$ul)
  MLE=mean(p.hat)
  ci=data.frame(cbind(p=p, MLE=MLE, LL=ll.conf, UL=ul.conf, length.dif=ul.conf-ll.conf))
  #out
  ci
}
set.seed(12345);wald.ci(10,15,0.30,gamma=0.05,0.99,0.99)
########################################################
##########################################################################################
#########################Coverage Probabilities###########################################
#########################Wald Coverage probability with testing errors in Inspection######
##########################################################################################
####################      Working Code at 98% sens and spec             ##################
############# To modify for the 90, 95 and 98% ###########################################
##########################################################################################
wald.Ret=function (toler=0.000001,n, k, p, gamma=0.05, sens, spec){
  con <- 1-gamma/2
  r=spec+sens-1
  theta=spec - ((1-p)^k)*(sens+spec-1)
  zalph <- qnorm(con,0,1)
  p.hat=c(); coverage=c()
  ul=c(); ll=c(); var.phat=c(); FisherInfo=c()
  tolcheck=1-toler
  tstar=qnbinom(tolcheck,n, theta)
  #is.covered=numeric(length(tstar))
  tb=c()
  tb= seq(n, tstar+n, by=1)
  #t=seq(n, tstar+n, by=1)
  for (i in seq_along(tb)) {
    if (n==tb[i]){
      p.hat[i]= 1-((sens-n/tb[i])/(r))^(1/k)
      ul[i]=1
      ll[i]=0
      var.phat[i]=0
    }
    else{
      p.hat[i]=1-((sens-n/tb[i])/(r))^(1/k)
      FisherInfo[i]=((r^2)*n*(k^2)*(1-p.hat[i])^(2*k-2))/(((sens-r*(1-p.hat[i])^k)^2)*(1-sens+r*(1-p.hat[i])^k))
      FisherInfo
      #print(FisherInfo)
      var.phat[i]=(FisherInfo[i])^-1
      ul[i]=p.hat[i] + zalph*sqrt(var.phat[i])
      ll[i]= p.hat[i] - zalph*sqrt(var.phat[i])
    }
  }
  is.covered= (ll <= p) & (p <= ul)###0.3, n=1 ###to check which values of t satisfies the inequality 
  b=tb[is.covered]
  coverage=sum(dnbinom((b-n), n, theta))
  return(coverage)
}
dwd=wald.Ret(0.000001,1,5,p=0.0025, 0.05, 0.99, 0.99)
dwd
#####################################################
p=seq(0.00005,0.10, by=0.00005)
pwbp=rep(NA, length(p))

#out=numeric(length(p))
for (i in 1:length(p)) {
  pwbp[i]=wald.Ret(0.000001,20,5,p[i],0.05, 0.99,0.99)
  
}

pwbp

plot(p, pwbp, type="l", ylab = "coverage Probabilities", col="blue4",
     main=" k= 20")
abline(h=0.95, col="black")