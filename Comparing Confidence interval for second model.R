#########################################################################
#######   COMPILED CODES FOR RETESTING 		###############################
#########################################################################
###Used for comparing the confidence intervals for the retesting model###


##########################################################################
######### 	WALD Confidence Interval	##########################
##########################################################################
Wald.confidence.Ret=function(toler=0.000001,n, k, p,gamma=0.05, sens, spec){
  pro= (((1-sens)^2)*(1-p)^k) + (spec^2)*(1-(1-p)^k)
  con = 1-gamma/2
  zalph = qnorm(con,0,1)
  p.hat= c(); t=c()
  g=c()
  tolcheck=1-toler
  tstar=qnbinom(tolcheck,n, pro)
  t=seq(n, tstar+n, by=1)
  var.phat=c()
  FisherInfo=c();b=c()
  is.covered=c()
  ul=c();ll=c();th=c()
  for (i in seq_along(t)) {
    
    th[i]=(((spec^2)-n/t[i])/((spec^2)-(1-sens)^2))^(1/k)
    g[i]=n/t[i]
    if (t[i]==n){
      p.hat[i]= 1
      ul[i]=1
      ll[i]=0
      var.phat[i]=0
    }
    else if((spec^2)<=g[i]){
      p.hat[i]=1
      ul[i]=1
      ll[i]=0
      var.phat[i]=0
    }
    else if(th[i]>=1){
      p.hat[i]=0
      FisherInfo[i]=n*(k^2)*(1-p.hat[i])^(2*k-2)*((spec^2)-(1-sens)^2)^2 * 1/((1-((1-sens)^2) *(1-p.hat[i])^k)-(spec^2)*(1-(1-p.hat[i])^k))*1/((((1-spec)^2) *(1-p.hat[i])^k)+ (spec^2)*(1-(1-p.hat[i])^k))^2
      var.phat[i]=(FisherInfo[i])^-1
      ul[i]=p.hat[i] + zalph*sqrt(var.phat[i])
      # ll[i]=0
      ll[i]= p.hat[i] - zalph*sqrt(var.phat[i])
    }
    else {
      p.hat[i] = 1-(((spec^2)-n/t[i])/((spec^2)-(1-sens)^2))^(1/k)
      FisherInfo[i]=n*(k^2)*(1-p.hat[i])^(2*k-2)*((spec^2)-(1-sens)^2)^2 * 1/((1-((1-sens)^2) *(1-p.hat[i])^k)-(spec^2)*(1-(1-p.hat[i])^k))*1/((((1-spec)^2) *(1-p.hat[i])^k)+ (spec^2)*(1-(1-p.hat[i])^k))^2
      var.phat[i]=(FisherInfo[i])^-1
      ul[i]=p.hat[i] + zalph*sqrt(var.phat[i])
      ll[i]= p.hat[i] - zalph*sqrt(var.phat[i])
    }
  }
  is.covered= (ll <= p) & (p <= ul)###0.3, n=1 ###to check which values of t satisfies the inequality 
  b=t[is.covered]
  coverage=sum(dnbinom((b-n), n, pro))
  return(coverage)
}

Wald.confidence.Ret(toler=0.000001,5,30,0.0025,gamma=0.05,0.99,0.99)

#######################################
####Plotting coverage Probability#####
p=seq(0.00005,0.10, by=0.00005)
pwbp=rep(NA, length(p))

for (i in 1:length(p)) {
  pwbp[i]=Wald.confidence.Ret(toler=0.000001,100,20,p[i],gamma=0.05,0.99,0.99)
}

pwbp
#################Drawing code########
plot(p, pwbp, type="l", col="blue", ylab="coverage probabilities", 
     main="Wald")
abline(h=0.95, col="red")
##################################################################
#########	EXACT INTERVAL 	 #################################
##################################################################
ExactConf.RET=function(n,p,toler=0.000001,k,gamma=0.05, sens, spec){
  con = 1-gamma/2
  pro= (((1-sens)^2)*(1-p)^k) + (spec^2)*(1-(1-p)^k)
  r=((spec^2)-(1-sens)^2)
  zalph = qnorm(con,0,1)
  p.hat=c(); coverage=c()
  pl=c(); pu=c()
  thetal=c(); thetau=c()
  is.covered=c(); coverage=c();b=c()
  th=c()
  tolcheck=1-toler
  tstar=qnbinom(tolcheck,n, pro)
  tb=c()
  tb= seq(n, tstar+n, by=1)
  for (i in seq_along(tb)) {
    th[i]=(((spec^2)-n/tb[i])/((spec^2)-(1-sens)^2))^(1/k)
    if(tb[i]==n){
      thetau[i]= (sens^2)
      thetal[i]= qbeta(gamma/2,n,tb[i]-n+1)
      pl[i] <- 1-((sens-thetal[i])/r)^(1/k)
      pu[i] <- 1-((sens-thetau[i])/r)^(1/k)
    }
    else if((spec^2)<=n/tb[i]){
      thetau[i]=1
      thetal[i]= qbeta(gamma/2,n,tb[i]-n+1)
      pl[i] = 1-(((spec^2)-n/tb[i])/r)^(1/k)
      pu[i] = 1-(((spec^2)-n/tb[i])/r)^(1/k)
    }
    else if (th[i]>=1){
      thetau[i] = qbeta(1-gamma/2,n,tb[i]-n)
      thetal= 0  ######sens
      pl[i]= 1-(((spec^2)-thetal[i])/r)^(1/k)
      pu[i]= 1-(((spec^2)-thetau[i])/r)^(1/k)
    }
    else{
      thetal[i]= qbeta(gamma/2,n,tb[i]-n+1)
      thetau[i]= qbeta(1-gamma/2,n,tb[i]-n)
      pl[i]= 1-(((spec^2)-thetal[i])/r)^(1/k)
      pu[i]= 1-(((spec^2)-thetau[i])/r)^(1/k)
    }
  }
  # pl[i]= 1-((sens-thetal[i])/r)^(1/k)
  # pu[i]= 1-((sens-thetau[i])/r)^(1/k)
  is.covered= (pl <= p) & (p <= pu)
  b=tb[is.covered]
  coverage=sum(dnbinom((b-n), n, pro), na.rm=T)
  return(coverage)
}
#####################################################
p=seq(0.00005,0.10, by=0.00005)
pwbp=rep(NA, length(p))

#out=numeric(length(p))
for (i in 1:length(p)) {
  pwbp[i]=ExactConf.RET(100,p[i],toler=0.00001,20,gamma=0.05,0.99, 0.99)
  
}

pwbp
min(pwbp);max(pwbp)

plot(p, pwbp, type="l", ylab = "coverage Probabilities", col="blue4",
     main=" Exact" )#, ylim =c(0.5, max(pwbp))
abline(h=0.95, col="red")

######################################################
#############	 AGRESTI-CODE		##############
######################################################
Agresti.coull.RET=function(n, k, p, gamma=0.05, toler=0.00001, spec, sens){
  pro= (((1-sens)^2)*(1-p)^k) + (spec^2)*(1-(1-p)^k)
  con = 1-gamma/2
  zalph = qnorm(con,0,1)
  p.hat= c(); t=c()
  g=c()
  tolcheck=1-toler
  tstar=qnbinom(tolcheck,n, pro)
  tb=seq((n+2), (tstar+n+4), by=1)
  var.phat=c()
  FisherInfo=c();b=c()
  is.covered=c()
  ul=c();ll=c();th=c()
  for (i in seq_along(tb)) {
    th[i]=(((spec^2)-n/tb[i])/((spec^2)-(1-sens)^2))^(1/k)
    g[i]=n/tb[i]
    if(tb[i]==n){
      p.hat[i]= 1
      ul[i]=1
      ll[i]=0
      var.phat[i]=0
    }
    else if((spec^2)<=n/tb[i]){
      p.hat[i]=1
      ul[i]=1
      ll[i]=0
      var.phat[i]=0  
    }
    else if(th[i]>=1){
      p.hat[i]=0
      FisherInfo[i]=n*(k^2)*(1-p.hat[i])^(2*k-2)*((spec^2)-(1-sens)^2)^2 * 1/((1-((1-sens)^2) *(1-p.hat[i])^k)-(spec^2)*(1-(1-p.hat[i])^k))*1/((((1-spec)^2) *(1-p.hat[i])^k)+ (spec^2)*(1-(1-p.hat[i])^k))^2
      var.phat[i]=(FisherInfo[i])^-1
      ul[i]=p.hat[i] + zalph*sqrt(var.phat[i])
      # ll[i]=0
      ll[i]= p.hat[i] - zalph*sqrt(var.phat[i])
      
    }
    else{
      p.hat[i] = 1-(((spec^2)-n/tb[i])/((spec^2)-(1-sens)^2))^(1/k)
      FisherInfo[i]=n*(k^2)*(1-p.hat[i])^(2*k-2)*((spec^2)-(1-sens)^2)^2 * 1/((1-((1-sens)^2) *(1-p.hat[i])^k)-(spec^2)*(1-(1-p.hat[i])^k))*1/((((1-spec)^2) *(1-p.hat[i])^k)+ (spec^2)*(1-(1-p.hat[i])^k))^2
      var.phat[i]=(FisherInfo[i])^-1
      ul[i]=p.hat[i] + zalph*sqrt(var.phat[i])
      ll[i]= p.hat[i] - zalph*sqrt(var.phat[i])
    }
    
  }
  is.covered= (ll <= p) & (p <= ul) 
  b=tb[is.covered]
  coverage=sum(dnbinom((b-n), n, pro))
  return(coverage)
}

Agresti.coull.RET(1,20,0.00025,0.05, 0.00001,0.99,0.99)
#####################################################
p=seq(0.00005,0.10, by=0.00005)
pwbp=rep(NA, length(p))

#out=numeric(length(p))
for (i in 1:length(p)) {
  pwbp[i]=Agresti.coull.RET(100,20,p[i],0.05, 0.00001,0.99,0.99)
  
}
#############PLOTTING CODE#####################
pwbp
min(pwbp);max(pwbp)

plot(p, pwbp, type="l", ylab = "coverage Probabilities", col="blue4",
     main=" Agresti-Coull" )#, ylim =c(0.5, max(pwbp))
abline(h=0.95, col="red")

############################################################################
#######		WILSON INTERVAL		####################################
############################################################################
Wilson.Interv.RET=function (n, k, p, gamma=0.05, toler=0.00001, sens, spec){
  con <- 1-gamma/2
  pro= (((1-sens)^2)*(1-p)^k) + (spec^2)*(1-(1-p)^k)
  zalph <- qnorm(con,0,1)
  p.hat=c(); coverage=c()
  thetal=c(); thetau=c()
  g=c()
  tolcheck=1-toler
  tstar=qnbinom(tolcheck,n, pro)
  tb=c()
  p.hat=c()
  ul=c(); ll=c()
  tb= seq(n, (tstar+n), by=1)
  pl=c(); pu=c()
  is.covered=c()
  coverage=c()
  b=c(); th=c()
  for (i in seq_along(tb)) {
    th[i]=(((spec^2)-n/tb[i])/((spec^2)-(1-sens)^2))^(1/k)
    g[i]=n/tb[i]
    if(tb[i]==n){
      p.hat[i]= 1
    }
    else if((spec^2)<=g[i]){
      p.hat[i]=1
    }
    else if(th[i]>=1){
      p.hat[i]=0
    }
    else{
      p.hat[i] = 1-(((spec^2)-n/tb[i])/((spec^2)-(1-sens)^2))^(1/k)
    }
    ul[i]=(2*tb[i]*p.hat[i]+(zalph^2)+(zalph)*sqrt((zalph^2)+4*tb[i]*p.hat[i]*(1-p.hat[i])))/(2*(tb[i]+(zalph^2)))
    ll[i]=(2*tb[i]*p.hat[i]+(zalph^2)-(zalph)*sqrt((zalph^2)+4*tb[i]*p.hat[i]*(1-p.hat[i])))/(2*(tb[i]+(zalph^2)))
  }
  is.covered= (ll <= p) & (p <= ul)###0.3, n=1 ###to check which values of t satisfies the inequality 
  b=tb[is.covered]
  #return(b)
  coverage=sum(dnbinom((b-n), n, pro))
  return(coverage)
}
Wilson.Interv.RET(1,20,0.000025,0.05, 0.00001,1, 1)

#####################################################
p=seq(0.00005,0.10, by=0.00005)
pwbp=rep(NA, length(p))

#out=numeric(length(p))
for (i in 1:length(p)) {
  pwbp[i]=Wilson.Interv.RET(100,20,p[i],0.05, 0.00001,0.99, 0.99)
  
}

pwbp
min(pwbp);max(pwbp)

plot(p, pwbp, type="l", ylab = "coverage Probabilities", col="blue4",
     main=" Wilson")#, ylim =c(0.5, max(pwbp))
abline(h=0.95, col="red")
