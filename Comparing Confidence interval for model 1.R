############### Coverage probability #########
##############################################
Pool.conf=function(n, k, p,gamma=0.05, toler=0.000001, sens, spec){
  Agesti.Int=Agresti.Coull(n=n, k=k, p=p, gamma = 0.05, toler=0.00001, sens = sens, spec=spec)
  Exact.Int=ExactConf(n=n, p=p, toler=0.00001, k=k, gamma = 0.05, sens = sens, spec=spec)
  Wilson.Int=Wilson.Interv(n=n,k=k, p=p, gamma = 0.05, toler=0.000001, sens=sens, spec=spec)
  wald.int=wald.conf(toler=0.000001, n=n, k=k, p=p, gamma=0.05, sens = sens, spec=spec)
  out=list(Exact.Int=Exact.Int,wald.int=wald.int,  Wilson.Int=Wilson.Int,  Agesti.Int=Agesti.Int )
  out
}

Pool.conf(n=20, k=20, p=0.1, toler=0.000001, sens=0.99, spec=0.99)

Pool.conf(n=100, k=20, p=0.99, toler=0.000001, sens=0.99, spec=0.99)


round(d, digits=7)


############Code########################
wald.conf=function (toler=0.000001,n, k, p, gamma=0.05, sens, spec){
  theta=spec - ((1-p)^k)*(sens+spec-1)
  r=sens+spec-1
  con <- 1-gamma/2
  zalph <- qnorm(con,0,1)
  p.hat=c(); coverage=c()
  ul=c(); ll=c(); var.phat=c(); FisherInfo=c()
  tolcheck=1-toler
  th=c()
  tstar=qnbinom(tolcheck,n, theta)
  #is.covered=numeric(length(tstar))
  t=seq(n, tstar+n, by=1)
  for (i in seq_along(t)) {
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
      p.hat[i]=1-((sens-n/t[i])/(sens+spec-1))^(1/k)
      FisherInfo[i]=((r^2)*n*(k^2)*(1-p.hat[i])^(2*k-2))/(((sens-r*(1-p.hat[i])^k)^2)*(1-sens+r*(1-p.hat[i])^k))
      var.phat[i]=(FisherInfo[i])^-1
      ul[i]=p.hat[i] + zalph*sqrt(var.phat[i])
      ll[i]= p.hat[i] - zalph*sqrt(var.phat[i])
    }
  }
  is.covered= (ll <= p) & (p <= ul)###0.3, n=1 ###to check which values of t satisfies the inequality 
  b=t[is.covered]
  coverage=sum(dnbinom((b-n), n, theta))
  return(coverage)
}
#############################
ExactConf=function(n,p,toler=0.000001,k,gamma=0.05, sens, spec){
  con = 1-gamma/2
  r=sens+spec-1
  theta= spec-r*((1-p)^k)
  zalph = qnorm(con,0,1)
  p.hat=c(); coverage=c()
  pl=c(); pu=c()
  thetal=c(); thetau=c()
  is.covered=c(); coverage=c();b=c()
  th=c()
  tolcheck=1-toler
  tstar=qnbinom(tolcheck,n, theta)
  tb=c()
  tb= seq(n, tstar+n, by=1)
  for (i in seq_along(tb)) {
    th[i]= ((spec-n/tb[i])/(sens+spec-1))^(1/k)
    if(tb[i]==n){
      thetau[i]= sens
      thetal[i]= qbeta(gamma/2,n,tb[i]-n+1)
      pl[i] <- 1-((sens-thetal[i])/r)^(1/k)
      pu[i] <- 1-((sens-thetau[i])/r)^(1/k)
    }
    else if(spec<=n/tb[i]){
      thetau[i]=sens
      thetal[i]= qbeta(gamma/2,n,tb[i]-n+1)
      pl[i] = 1-((sens-thetal[i])/r)^(1/k)
      pu[i] = 1-((sens-thetau[i])/r)^(1/k)
    }
    else if (th[i]>=1){
      thetau[i] = qbeta(1-gamma/2,n,tb[i]-n)
      thetal= 0  ######sens
      pl[i]= 1-((sens-thetal[i])/r)^(1/k)
      pu[i]= 1-((sens-thetau[i])/r)^(1/k)
    }
    else{
      thetal[i]= qbeta(gamma/2,n,tb[i]-n+1)
      thetau[i]= qbeta(1-gamma/2,n,tb[i]-n)
      pl[i]= 1-((sens-thetal[i])/r)^(1/k)
      pu[i]= 1-((sens-thetau[i])/r)^(1/k)
    }
  }
  # pl[i]= 1-((sens-thetal[i])/r)^(1/k)
  # pu[i]= 1-((sens-thetau[i])/r)^(1/k)
  is.covered= (pl <= p) & (p <= pu)
  b=tb[is.covered]
  coverage=sum(dnbinom((b-n), n, theta), na.rm=T)
  return(coverage)
}

##############################
Agresti.Coull=function (n, k, p, gamma=0.05, toler=0.000001, spec, sens){
  con <- 1-gamma/2
  theta=spec - ((1-p)^k)*(sens+spec-1)
  r=sens+spec-1
  zalph <- qnorm(con,0,1)
  p.hat=c(); coverage=c()
  thetal=c(); thetau=c()
  g=c()
  tolcheck=1-toler
  tstar=qnbinom(tolcheck,n, theta)
  tb=c()
  p.hat=c()
  ul=c(); ll=c()
  tb= seq((n+2), (tstar+n+4), by=1)
  pl=c(); pu=c()
  is.covered=c()
  FisherInfo=c()
  var.phat=c()
  th=c()
  h=c()
  for (i in seq_along(tb)) {
    th[i]= ((spec-n/tb[i])/(sens+spec-1))^(1/k)
    g[i]=n/tb[i]
    if (tb[i]==n){
      p.hat[i]=1
      ul[i]=1
      ll[i]=0
      var.phat[i]=0
    }
    else if (sens<= n/tb[i]){
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
      p.hat[i]=1-(1-n/(tb[i]))^(1/k)
      FisherInfo[i]=(n*(k^2)*(1-p.hat[i])^(k-2))/((1-(1-p.hat[i])^k)^2)
      #var.phat=1/FisherInfo
      var.phat[i]=(FisherInfo[i])^-1
      ul[i]=p.hat[i]+zalph*sqrt(var.phat[i])
      ll[i]=p.hat[i]-zalph*sqrt(var.phat[i])
    }
    
  }
  is.covered= (ll <= p) & (p <= ul) 
  h=tb[is.covered]
  #coverage=sum(dnbinom(b-n, n, theta))
  coverage=sum(dnbinom((h-n), n, theta), na.rm=T)
  return(coverage)
}

###############################################
Wilson.Interv=function (n, k, p, gamma=0.05, toler=0.000001, sens, spec){
  con <- 1-gamma/2
  theta= spec - ((1-p)^k)*(sens+spec-1)
  r= sens+spec-1
  zalph = qnorm(con,0,1)
  p.hat=c(); coverage=c()
  thetal=c(); thetau=c()
  g=c()
  tolcheck=1-toler
  tstar=qnbinom(tolcheck,n, theta)
  tb=c()
  p.hat=c()
  ul=c(); ll=c()
  tb= seq(n, (tstar+n), by=1)
  pl=c(); pu=c()
  is.covered=c()
  coverage=c()
  d=c(); th=c()
  for (i in seq_along(tb)) {
    th[i]= ((spec-n/tb[i])/(sens+spec-1))^(1/k)
    if(tb[i]==n){
      p.hat[i]= 1
    }
    else if(spec<=n/tb[i]){
      p.hat[i]=1
    }
    else if(th[i]>=1){
      p.hat[i]=0
    }
    else{
      p.hat[i]=1-((spec-n/tb[i])/(sens+spec-1))^(1/k)
    }
    ul[i]=(2*tb[i]*p.hat[i]+(zalph^2)+(zalph)*sqrt((zalph^2)+4*tb[i]*p.hat[i]*(1-p.hat[i])))/(2*(tb[i]+(zalph^2)))
    ll[i]=(2*tb[i]*p.hat[i]+(zalph^2)-(zalph)*sqrt((zalph^2)+4*tb[i]*p.hat[i]*(1-p.hat[i])))/(2*(tb[i]+(zalph^2)))
  }
  is.covered= (ll <= p) & (p <= ul)###0.3, n=1 ###to check which values of t satisfies the inequality 
  d=tb[is.covered]
  d
  #return(b)
  coverage=sum(dnbinom((d-n), n, theta))
  return(coverage)
}
