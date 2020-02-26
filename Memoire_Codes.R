############ MÉMOIRE: Maitrise en Mathématiques Finacières et Actuarielle  #################

library(ggplot2)

# -SECTION 1: Concavification -#

payoff_call<-function(X_T,b,a,K){
  # X_T: valeur du stock à maturité
  # b: Strike value
  # K: valeur initiale du Payoff
  
  # Call de la forme:  a*max(X_T-b,0)+K
  payoff<-a*max(X_T-b,0)+K
  return(payoff)
}
payoff_call(10,5,0.5,8)

P_Utility<-function(X_Tu,gamma){
  res<-(X_Tu)^(1-gamma)/(1-gamma)
  return(res)
}

Find_x_theta_PU<-function(b_o,a_o,K_o,gamma_o){
  f <- function(x) payoff_call(X_T=x,b=b_o,a=a_o,K=K_o)^(-gamma_o)*x*a_o-P_Utility(payoff_call(x,b=b_o,a=a_o,K=K_o),gamma_o)+P_Utility(K_o,gamma_o)        #-Finding the fee
  result<-uniroot(f, c(b_o,b_o+100))$root
  return(result)
}
Find_x_theta_PU(b_o=30,a_o=0.5,K_o=15,gamma_o=4)

P_Utility_con<-function(X_T_con,gamma_c,b_c,a_c,K_c){
  x_theta<-Find_x_theta_PU(b_o=b_c,a_o=a_c,K_o=K_c,gamma_o=gamma_c)
  slope<-(P_Utility(payoff_call(x_theta,b_c,a_c,K_c),gamma_c)-P_Utility(K_c,gamma_c))/(x_theta)
    
  if(X_T_con<=x_theta){Chord<-P_Utility(K_c,gamma_c)+X_T_con*slope}
  else{Chord<-P_Utility(payoff_call(X_T=X_T_con,b=b_c,a=a_c,K=K_c),gamma_c)}
  
  return(Chord)
}


Derivee_P_Utility_con<-function(X_T_con,gamma_c,b_c,a_c,K_c){
  x_theta<-Find_x_theta_PU(b_o=b_c,a_o=a_c,K_o=K_c,gamma_o=gamma_c)
  slope<-(P_Utility(payoff_call(x_theta,b_c,a_c,K_c),gamma_c)-P_Utility(K_c,gamma_c))/(x_theta)
  
  if(X_T_con<=x_theta){Chord<-slope}
  else{Chord<-payoff_call(X_T=X_T_con,b=b_c,a=a_c,K=K_c)^(-gamma_c)*a_c}
  
  return(Chord)
}
  
Convexe_Fun<-function(element){ #Fonction convexe
  res<-(element)^5
  return(res)
}

droite<-function(pente, origine,val_x){
  return(origine+val_x*pente)
}

find_roots<-function(pente1,origine1,beg,end){
  f <- function(x) droite(pente1,origine1,x)-Convexe_Fun(x)
  result1<-uniroot(f, c(beg,end))$root
  return(result1)
}
find_roots(15000,-20000,0,4)

X_opt_given_S<-function(S_T,gamma,sigma,alpha,r){
  result<-(S_T*gamma*sigma^2)/(alpha-r)
  return(result)
}

CRRA<-function(x,gamma){
  result<-(x^(1-gamma)-1)/(1-gamma)
  return(result)
}
####### -SECTION 2: Graphes Concavification -#######
# -Paramètres- #
par_b<-30
par_a<-300
par_K<-15
par_Gamma<-0.75


#-Call Payoff-#
axex<-seq(0,100,1/10)
axey<-lapply(axex,function(x) payoff_call(X_T=x,b=par_b,a=par_a,K=par_K))

#-Utility of a Temrinal Value-#
axey_u<-lapply(axex,function(x) P_Utility(X_Tu=x,gamma=par_Gamma))

#-Utility of a Call Payoff-#
axey_ut<-lapply(axey,function(x) P_Utility(X_Tu=x,gamma=par_Gamma))

#-Concavified Utility of a Call Payoff-#
axey_concav<-lapply(axex,function(x) P_Utility_con(X_T_con=x,gamma_c=par_Gamma,b_c=par_b,a_c=par_a,K_c=par_K))

# -Fonction convexe: Théorique - #
axex2<-seq(0,12,1/10)
axey_convexe<-lapply(axex2,function(x) Convexe_Fun(x))
axey_droite<-lapply(axex2,function(x) droite(11000,-14000,x))

### -GRAPHES- ###
par(mfrow=c(1,1))

# Option d'achat#
plot(axex,axey,ylab = "Bénéfice résultant d'une option d'achat",type='l',lty=1,xlab=expression(paste("Valeur de l'actif à maturité (X"[T],")" )),yaxt = "n" ,xaxt = "n",xaxs="i" )
xtick<-c(0, b=par_b)
axis(side=1, at=xtick,labels=c(0,expression("D"[T])))
axis(side=2, at=par_K,labels=c(expression("K")),las=2)

#Utilité#
axex_Utility<-seq(1,2.5,1/100)
plot(axex_Utility,lapply(axex_Utility,function(x)P_Utility(x,2)),ylab = "U(x)",type='l',lty=1,xlab="x",ylim = c(-1.05,0.1),xlim=c(0.8,2.5))#,yaxt = "n" ,xaxt = "n"
lines(axex_Utility,lapply(axex_Utility,function(x)P_Utility(x,3)),col='saddlebrown',lty=2)
lines(axex_Utility,lapply(axex_Utility,function(x)P_Utility(x,4)),col='tomato3',lty=5)
lines(axex_Utility,lapply(axex_Utility,function(x)P_Utility(x,8)),col='red',lty=4)
text(0.93, -1, labels =  expression(paste(gamma,"=2")))
text(0.93, -0.52, labels =  expression(paste(gamma,"=3")))
text(0.93, -0.35, labels =  expression(paste(gamma,"=4")))
text(0.93, -0.14, labels = expression(paste(gamma,"=8")))

axex_Utility<-seq(0.5,2.5,1/100)
plot(axex_Utility,lapply(axex_Utility,function(x)P_Utility(x,2)),ylab = "U(x)",type='l',lty=1,xlab="x",ylim = c(-1.1,0.03),xlim=c(0.75,2.7))#,yaxt = "n" ,xaxt = "n"
lines(axex_Utility,lapply(axex_Utility,function(x)P_Utility(x,3)),col='saddlebrown',lty=2)
lines(axex_Utility,lapply(axex_Utility,function(x)P_Utility(x,4)),col='tomato3',lty=5)
lines(axex_Utility,lapply(axex_Utility,function(x)P_Utility(x,8)),col='red',lty=4)
text(2.6, -0.39, labels =  expression(paste(gamma,"=2")))
text(2.6, -0.09, labels =  expression(paste(gamma,"=3")))
text(2.6, -0.03, labels =  expression(paste(gamma,"=4")))
text(2.6, 0.02, labels = expression(paste(gamma,"=8")))


#Concavification#
plot(axex,axey_ut,ylab = "Utilité",xlab=expression(paste("Valeur de l'actif à maturité (x)" )),type='l',yaxt = "n" ,xaxt = "n",lty=2,xaxs="i" )
xtick<-c( 0,b=par_b,round(Find_x_theta_PU(b_o=par_b,a_o=par_a,K_o=par_K,gamma_o=par_Gamma),0))
axis(1, at=xtick, labels=c(0 ,"b",expression(tilde(x))))
lines(axex, axey_concav, col = "red",lty=4)
legend(par_b+25, par_K+10, legend=c( expression(paste("Fonction concavifiée  ",tilde(u)(x))), "Fonction originale u(x)"),
       col=c("red", "black"), lty=c(4,2), cex=0.8)

#Fonction convexes#
plot(axex2,axey_convexe,ylab = "",type='l',xlab="", xlim = c(-1,13),ylim = c(-50000,250000),yaxt = "n" ,xaxt = "n")
#abline(-14000,11000,type='l',col='red')
lines(axex2,lapply(axex2,function(x) droite(15000,-20000,x)),col='red')
xtick<-c(find_roots(15000,-20000,0,4),find_roots(15000,-20000,8,12))
axis(1, at=xtick, labels=c(expression(x[1]),expression(x[2])))
ytick<-c(Convexe_Fun(find_roots(15000,-20000,0,4)),Convexe_Fun(find_roots(15000,-20000,8,12)))
axis(2, angle= 45, at=ytick, labels=c(expression(paste("g(x"[1],")")),expression(paste("g(x"[2],")"))))

#X^*_T en fonction de S_T#
par(mfrow=c(1,1))

axex3<-seq(0,40,1)
plot(axex3,lapply(axex3,function(x) X_opt_given_S(x,gamma=2,sigma=0.5,alpha=0.07,r=0.03)),xlab=expression(S[T]),ylab=expression(paste("X","*"[T])),type="l",xaxs="i",yaxs="i",ylim = c(0,2500))
lines(axex3,lapply(axex3,function(x) X_opt_given_S(x,gamma=3,sigma=0.5,alpha=0.07,r=0.03)),col='gray')
lines(axex3,lapply(axex3,function(x) X_opt_given_S(x,gamma=4,sigma=0.5,alpha=0.07,r=0.03)),col='orange')
lines(axex3,lapply(axex3,function(x) X_opt_given_S(x,gamma=5,sigma=0.5,alpha=0.07,r=0.03)),col='red')
legend(5, 2400, legend=c( expression(paste(gamma,"=2")), expression(paste(gamma,"=3")),expression(paste(gamma,"=4")),expression(paste(gamma,"=5"))),
       col=c("black", "gray","orange","red"), lty=c(1,1,1,1), cex=0.8)

plot(axex3,lapply(axex3,function(x) P_Utility(x,gamma=2)),xlab=expression(S[T]),ylab=expression(Utilité),type="l",xaxs="i",yaxs="i",ylim = c(-1,0.2))
lines(axex3,lapply(axex3,function(x) P_Utility(x,gamma=3)),col='gray')
lines(axex3,lapply(axex3,function(x) P_Utility(x,gamma=4)),col='orange')
lines(axex3,lapply(axex3,function(x) P_Utility(x,gamma=5)),col='red')
legend(10, -0.2, legend=c( expression(paste(gamma,"=2")), expression(paste(gamma,"=3")),expression(paste(gamma,"=4")),expression(paste(gamma,"=5"))),
       col=c("black", "gray","orange","red"), lty=c(1,1,1,1), cex=0.8)


#Dérivée de la fonction concavifiée#
plot(axex,lapply(axex,function(x)Derivee_P_Utility_con(X_T_con=x,gamma_c=par_Gamma,b_c=par_b,a_c=par_a,K_c=par_K)),ylab = c( expression(paste(tilde(u)," '",(x)))),xlab=expression(paste("Valeur de l'actif à maturité (x)" )),type='l',col="red" ,xaxt = "n",yaxt = "n",lty=1,xaxs="i",ylim=c(0,1)  )
xtick<-c( 0,round(Find_x_theta_PU(b_o=par_b,a_o=par_a,K_o=par_K,gamma_o=par_Gamma),0))
ytick<-Derivee_P_Utility_con(0,gamma_c=par_Gamma,b_c=par_b,a_c=par_a,K_c=par_K)
axis(1, at=xtick, labels=c(0 ,expression(tilde(x))))
axis(2, at=ytick, labels=c( expression(paste(tilde(u)," '",(tilde(x))))))


##### -TEST SECTION- #####

f <- function(x) x^2-x-4         #-Finding the fee
r<-uniroot(f, c(0,100))

plot(1:10, 1:10, type="l", lty=2, lwd=3)
plot(1:10, 1:10, type="l", lty=1)
plot(axex,axey,ylab = "Résultat de l'option",type='l',lty=3,xlab="Valeur terminale de X",yaxt = "n" ,xaxt = "n")
xtick<-c(0, b=30)

pnorm(0)