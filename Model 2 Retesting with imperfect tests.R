##########################################################################################
##########################################################################################
###########RETESTING CODES################################################################
##########################################################################################
################MLE of P.hat for the retesting model######################################
###############
###############MLE of P.
#########code returns a 10, 000 montecarlos simulation of MLE for each value of t.
###Modify code to produce single value.
######i.e return(mean(out$mle)) 
######i.e return(mean(out$mle)) 
pest.Ret=function(n, k, p, sens, spec){
  #r=sens+spec-1
  pro= (((1-sens)^2)*(1-p)^k) + (spec^2)*(1-(1-p)^k)
  t=c()
  MLE= c()
  g=c(); th=c()
  out=c()
  N=10000
  for (i in 1:N) {
    # if(i%%1000==0){print(i)} ####modulo to show wether the code is running & iterating
    t[i] = rnbinom(1, n, pro)+n
    th[i]=(((spec^2)-n/t[i])/((spec^2)-(1-sens)^2))^(1/k)
    g[i]=n/t[i]
    
    if(t[i]==n){
      MLE[i]= 1
    }
    else if((spec^2)<=g[i]){
      MLE[i]=1
    }
    else if(th[i]>1){
      MLE[i]=0
    }
    else{
      MLE[i] = 1-(((spec^2)-n/t[i])/((spec^2)-(1-sens)^2))^(1/k)
    }
  }
  return(mean(MLE))
  #out=as.data.frame(cbind(mle=MLE))
  #return(out)
}

################Computation of MLE###########
p=c(0.005,0.01,0.05,0.10,0.20,0.30)
Frap=rep(NA, length(p))

for (i in 1:length(p)){
  set.seed(12345); Frap[i]=pest.Ret(1,50,p[i],0.90,0.90)
  d=round(Frap, digits=6)
}
cbind(d)
#########################################################################################
######################Bias and MSE#######################################################
###############using for function
#####################################################################
####################################################################
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
    else if(th[i]>1){
      MLE[i]=0
    }
    else{
      MLE[i] = 1-(((spec^2)-n/t[i])/((spec^2)-(1-sens)^2))^(1/k)
    }
    diffl[i]=abs(MLE[i]-p)
    Bias.p[i]=diffl[i]*pmft[i]
    mse.p[i]=(MLE[i]-p)^2 *(pmft[i])
    BIAS=sum(Bias.p)*10^4
    MSE=sum(mse.p)*10^4
  }
  return(MSE)
}

MSE.RET(0.00001,1,5,0.40, 0.99,0.99)
MSE.RET(0.00001,1,5,0.005, 0.99,0.99)



##############
p=c(0.005,0.01,0.05,0.10,0.20,0.30,0.40)
#p=seq(0,1, by=0.1)
p
#p=c(0,0.1,0.2)
Frap=rep(NA, length(p))
#Frap1=rep(NA, length(p))
#Frap2=rep(NA, length(p))
#Frap3=rep(NA, length(p))

for (i in seq_along(p)) {
  Frap[i]=MSE.RET(0.00001, 1, 50, p[i], 0.90, 0.90)
  
}
#########################################################################
#########################################################################
##########Alternative code that generates all the Bias for p	#########
#########################################################################
###############using for function
MSE.RET=function(toler, n, k, p, sens, spec){
  tstar=c();g=c();Bias.p=c();mse.p=c()
  t=c();pmft=c();MLE=c()
  r=(spec^2)-(1-sens)^2
  prob= (((1-sens)^2)*(1-p)^k) + (spec^2)*(1-(1-p)^k)
  tolcheck=1-toler
  tstar=qnbinom(tolcheck,n, prob=prob)
  t=seq(n, tstar+n, by=1)
  #print(t)
  for (i in seq_along(t)) {
    pmft[i]=dnbinom((t[i]-n), n, prob=prob)
    g[i]=n/t[i]
    if(g[i]==1){
      MLE[i]=1
    }
    else if ((spec^2)<=g[i]){
      MLE[i]=1
    }
    else{
      MLE[i] = 1-(((spec^2)-n/t[i])/(r))^(1/k)
    }
    #print(t)
    Bias.p=(MLE[i]-p)*pmft
    
    mse.p=(MLE[i]-p)^2 *pmft
    Bias=sum(Bias.p)
    #print(Bias)
    MSE=sum(mse.p)
    #print(MSE)
    out=list(BIAS=Bias*10^4)
    #output=list(BIAS=Bias, MSE=MSE)
    #output=as.data.frame(cbind(BIAS=Bias, MSE=MSE))
    ##output=as.data.frame(cbind(t=t,MLE=MLE, Bias.p=Bias.p, mse.p=mse.p))
    #output=as.data.frame(cbind(t=t,MLE=MLE,pmft=pmft, Bias.p=Bias.p, Bias=Bias, mse.p=mse.p, MSE=MSE))
  }
  return(out)
}

MSE.RET(0.00001,1,5,0.05, 0.99,0.99)

#####################################################################
p=c(0.005, 0.01, 0.05, 0.10, 0.20, 0.30)
Frap=rep(NA, length(p))
for (i in 1:length(p)){
  set.seed(12345);Frap[i]=MSE.RET(0.00001,1,5,p[i], 0.99,0.99)
}
Frap


####################################################################

########################################################################################################
########################################################################################################
########################Confidence Interval#############################################################

########################Confidence Interval#############################################################
confidence.Ret=function(n, k, p, sens, spec, gamma=0.05){
  pro= (((1-sens)^2)*(1-p)^k) + (spec^2)*(1-(1-p)^k)
  con <- 1-gamma/2
  zalph <- qnorm(con,0,1)
  t=c()
  p.hat= c()
  g=c()
  N=10000
  var.phat=c() #output=as.data.frame()
  FisherInfo=c()
  ul=c();ll=c();th=c()
  for (i in 1:N) {
    t[i] = rnbinom(1, n, pro)+n
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
      MLE[i]=0
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
  ##out=as.data.frame(cbind(t=t, MLE=p.hat))
  ##out
  out=as.data.frame(cbind(T=t, MLE=p.hat, ll=ll, ul=ul, zalph, var.phat))
  #out
  ll.conf=mean(out$ll)
  ul.conf=mean(out$ul)
  ci=cbind(LL=ll.conf, UL=ul.conf, length.dif=ul.conf-ll.conf)
  ci
}

options(digits = 5);set.seed(12345); d=confidence.Ret(1,5,0.005,0.99,0.99, 0.05)
d
View(d)





#########################################################################################


#########################################################
#                     MODEL COMPARISON                  #
#########################################################
##########################################################################################
###############ASYMPTOTIC RELATIVE EFFECIENCY#############################################
########################    ARE  ########################################################
##################################################
###############PERFECT CODE FOR ARE ##############
Asympt.ARE=function(n, k, p, sens, spec){
  r= 2*(sens)-1
  FisherWithout=((r^2)*n*(k^2)*(1-p)^(2*k-2))/(((sens-r*(1-p)^k)^2)*(1-sens+r*(1-p)^k))
  varWithout=1/FisherWithout
  #varWithout
  FisherInfoRET=(n*(k^2)*(1-p)^(2*k-2)*((spec^2)-(1-sens)^2)^2) *(1/((1-((1-sens)^2) *(1-p)^k)-(spec^2)*
                                                                       (1-(1-p)^k)))*(1/((((1-spec)^2) *(1-p)^k)+ (spec^2)*(1-(1-p)^k))^2)
  var.RET=(FisherInfoRET)^-1
  ARE=varWithout/var.RET
  #out=(list(ARE=ARE, VAR.Prit=varWithout, var.RET=var.RET))
  #return(out)
  return(ARE)
}
Asympt.ARE(10, 50,0.0025,0.99,0.99)

#############
p=c(0.005, 0.01,0.05,0.10,0.20,0.30)
p
are=rep(NA, length(p))
for (i in 1:length(p)) {
  are[i]=Asympt.ARE(1, 15, p[i], 0.99, 0.99)
  #are
}
cbind(are)

##################################################
###########initial code with error################
##################################################
#			Relative Mean Squared Error (ARE-CODES)    #
##################################################
Asympt.ARE=function(n, k, p, sens, spec){
  r= 2*(sens)-1
  FisherWithout=((r^2)*n*(k^2)*(1-p)^(2*k-2))/(((sens-r*(1-p)^k)^2)*(1-sens+r*(1-p)^k))
  varWithout=1/FisherWithout
  #varWithout
  FisherInfo=n*(k^2)*(1-p)^(2*k-2)*((spec^2)-(1-sens)^2)^2 * 1/((1-((1-sens)^2) *(1-p)^k)-(spec^2)*
                                                                  (1-(1-p)^k))*1/((((1-spec)^2) *(1-p)^k)+ (spec^2)*(1-(1-p)^k))^2
  var.RET=(FisherInfo)^-1
  ARE=varWithout/var.RET
  out=(list(ARE=ARE, VAR.P=varWithout, var.RET=var.RET))
  return(out)
}
Asympt.ARE(10, 50,0.04,0.99,0.99)


###########################################################################################
#			Relative Mean Squared Error (ARE-CODES)                                       #
###########################################################################################
