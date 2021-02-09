# Nour Hawila
# Deep Learning HW1
# Implementing PPR for natural cubic spline

library(MASS)
library(splines)


#function to compute basis for natural cubic spline
spline <- function(v) {
  # creating knots
  if(n<50) knot = sort(v,dec=F) else knot = seq(min(v),max(v),length.out=30)
  
  N = matrix(NA,nrow=n,ncol=length(knot)-2)
  for(i in  1:n) {
    for(k in 1:(length(knot)-2)) {
      t_k = ifelse(v[i]>knot[k],(v[i]-knot[k])^3,0)
      t_K = ifelse(v[i]>knot[length(knot)],(v[i]-knot[length(knot)])^3,0)
      t_K1 = ifelse(v[i]>knot[length(knot)-1],(v[i]-knot[length(knot)-1])^3,0)
      N[i,k]=(t_k-t_K)/(knot[k]-knot[length(knot)])-(t_K1-t_K)/(knot[length(knot)-1]-knot[length(knot)])
    }
  }
  N
}

#function to calculate derivative of g(v)
derivative <- function(v) {
  # creating knots
  if(n<50) knot = sort(v,dec=F) else knot = seq(min(v),max(v),length.out=30)
  
  N = matrix(NA,nrow=n,ncol=length(knot)-2)
  for(i in  1:n) {
    for(k in 1:(length(knot)-2)) {
      t_k = ifelse(v[i]>knot[k],(v[i]-knot[k])^2,0)
      t_K = ifelse(v[i]>knot[length(knot)],(v[i]-knot[length(knot)])^2,0)
      t_K1 = ifelse(v[i]>knot[length(knot)-1],(v[i]-knot[length(knot)-1])^2,0)
      N[i,k]=(3*t_k-3*t_K)/(knot[k]-knot[length(knot)])-(3*t_K1-3*t_K)/(knot[length(knot)-1]-knot[length(knot)])
    }
  }
  N
}


# calculating spline function g(v_i)
g = function(omega, x, y) {
  v=t(omega)%*%t(x)
  Z = cbind(data.frame("v"=t(v),spline(v)))
  predictors = apply(Z,2,function(x) {x})
  reg = lm(y~predictors)
  coeff = as.matrix(as.numeric(reg$coefficients))
  coeff[is.na(coeff)]=0
  as.matrix(cbind(rep(1,n),Z))%*%coeff
}

#calculating derivative of spline function g'(v_i)
g.der = function(omega, x, y) {
  v=t(omega)%*%t(x)
  Z = cbind(data.frame("v"=t(v),spline(v)))
  predictors = apply(Z,2,function(x) {x})
  reg = lm(y~predictors)
  coeff = as.matrix(as.numeric(reg$coefficients))
  coeff[is.na(coeff)]=0
  as.matrix(cbind(rep(1,n),derivative(v)))%*%coeff[-1,]
}

#finding new omega from old omega
omega.new = function(omega, x, y) {
  omega.upd = rep(NA,length(omega))
  v=t(omega)%*%t(x)
  for(j in 1:length(omega)){
    num = g.der(omega, x, y)^2 * (t(v)+(y-g(omega,x,y))/g.der(omega,x,y))
    den = g.der(omega,x,y)^2
    omega.upd[j] = sum(num*x[,j])/sum(den*x[,j]^2)
  }
  omega.upd/(sqrt(sum(omega.upd^2)))
}


#iterate until convergence
find.solution <- function(omega.start, x, y, tol=1e-3) {
  
  g.old = g(omega.start,x,y)
  omega.old = omega.start
  
  omega.upd = omega.new(omega.old,x,y)
  g.new = g(omega.upd,x,y)
  counter=1
  while(mean(abs(omega.upd-omega.old))>tol | mean(abs(g.new-g.old))>tol) {
    print(c(counter, omega.upd))
    flush.console()
    
    omega.old = omega.upd
    g.old = g(omega.old,x,y)
    omega.upd = omega.new(omega.old,x,y)
    g.new = g(omega.upd,x,y)
    counter=counter+1
  }

list("omega" = omega.upd, "g" = g.new)
}


#####################################
######   RUNNING THE PROGRAM   ######
#####################################
set.seed(123)
n = 50
x1 = rnorm(n, 0, 1)
x2 = rnorm(n, 0, 1)
x=cbind(x1,x2)
y = exp(x1)/exp(x2+2) + rnorm(n, 0, 0.5)
data = data.frame(y,x)
omega.start = rbind(1/sqrt(2),1/sqrt(2))
tol=1e-3

ex1 = find.solution(omega.start,x,y)
check =   ppr(y ~ x1 + x2, data = data, nterms = 1, max.terms = 5, sm.method='spline', df = 3)

summary(check)
data.frame("y"=y, "Estimate" = ex1$g)

