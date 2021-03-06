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

lambda_cas_de_base<-function(xini,gamma, Maturite,alpha, r, sigma){
  theta<-(alpha-r)/sigma
  lambda_opti<-(1/xini)^(gamma)*exp((1-gamma)*(r+0.5/gamma*theta^2)*Maturite)
  return(lambda_opti)
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


par(mfrow=c(1,1))
axex_Utility<-seq(0.5,2.5,1/100)
plot(axex_Utility,lapply(axex_Utility,function(x)P_Utility(x,2)),ylab = "U(x)",type='l',lty=1,xlab="x",ylim = c(-1.1,0.03),xlim=c(0.75,2.7))#,yaxt = "n" ,xaxt = "n"
lines(axex_Utility,lapply(axex_Utility,function(x)P_Utility(x,3)),col='saddlebrown',lty=2)
lines(axex_Utility,lapply(axex_Utility,function(x)P_Utility(x,4)),col='tomato3',lty=5)
lines(axex_Utility,lapply(axex_Utility,function(x)P_Utility(x,8)),col='red',lty=4)
text(2.6, -0.39, labels =  expression(paste(gamma,"=2")),cex=0.9)
text(2.6, -0.09, labels =  expression(paste(gamma,"=3")),cex=0.9)
text(2.6, -0.03, labels =  expression(paste(gamma,"=4")),cex=0.9)
text(2.6, 0.02, labels = expression(paste(gamma,"=8")),cex=0.9)


### Dérivée des fonctions d'utilités précédentes ###
#Mettre par(mfrow=c(2,1)) pour mettre les graphiques côte-à-côte.
plot(axex_Utility,lapply(axex_Utility,function(x)derivee_Power_Utility(x,2)),ylab = "U'(x)",type='l',lty=1,xlab="x",ylim = c(-0.1,1.5),xlim=c(0.75,2.7))#,yaxt = "n" ,xaxt = "n"
lines(axex_Utility,lapply(axex_Utility,function(x)derivee_Power_Utility(x,3)),col='saddlebrown',lty=2)
lines(axex_Utility,lapply(axex_Utility,function(x)derivee_Power_Utility(x,4)),col='tomato3',lty=5)
lines(axex_Utility,lapply(axex_Utility,function(x)derivee_Power_Utility(x,8)),col='red',lty=4)
text(2.6, 0.18, labels =  expression(paste(gamma,"=2")),cex=0.9)
text(2.6, 0.095, labels =  expression(paste(gamma,"=3")),cex=0.9)
text(2.6, 0.025, labels =  expression(paste(gamma,"=4")),cex=0.9)
text(2.6, -0.04, labels = expression(paste(gamma,"=8")),cex=0.9)



#Concavification#
plot(axex,axey_ut,ylab = "Utilité",xlab=expression(paste("Valeur de l'actif à maturité (x)" )),type='l',yaxt = "n" ,xaxt = "n",lty=2,xaxs="i" )
xtick<-c( 0,b=par_b,round(Find_x_theta_PU(b_o=par_b,a_o=par_a,K_o=par_K,gamma_o=par_Gamma),0))
axis(1, at=xtick, labels=c(0 ,expression(d),expression(paste(tilde(x),'(d)'))))
axis(2, at=axey_ut[1], labels=c(expression(paste(U,'(K)'))))
lines(axex, axey_concav, col = "red",lty=4)
legend(par_b+25, par_K+10, legend=c( expression(paste("Fonction concavifiée  ",tilde(u)(x,d))), "Fonction originale u(x,d)"),
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
plot(axex,lapply(axex,function(x)Derivee_P_Utility_con(X_T_con=x,gamma_c=par_Gamma,b_c=par_b,a_c=par_a,K_c=par_K)),ylab = c( expression(paste(tilde(u)," '","(x,d)"))),xlab=expression(paste(x)),type='l',col="red" ,xaxt = "n",yaxt = "n",lty=1,xaxs="i",ylim=c(0,1)  )
xtick<-c( 0,round(Find_x_theta_PU(b_o=par_b,a_o=par_a,K_o=par_K,gamma_o=par_Gamma),0))
ytick<-Derivee_P_Utility_con(0,gamma_c=par_Gamma,b_c=par_b,a_c=par_a,K_c=par_K)
axis(1, at=xtick, labels=c(0 ,expression(paste(hat(x),(d)))))
axis(2, at=ytick, labels=expression(paste(tilde(u),"'( ",hat(x),(d),",d )")))

#Graphe: Fonction concave#
x_1<-1.5
x_2<-3.95
psi<-0.25
x_mel<-psi*x_1+(1-psi)*x_2
y_mel<-psi*(x_1^5.5)+(1-psi)*(x_2^5.5)
  
axe_d_x<-seq(1,4,0.05)
plot(axe_d_x,lapply(axe_d_x,function(x)x^5.5), type='l',ylim=c(-200,2048),xaxt = "n",yaxt = "n",xlab = '',ylab='')
xtick<-c(x_1,x_2,x_mel)
ytick<-c(lapply(x_1,function(x)x^5.5),lapply(x_2,function(x)x^5.5),lapply(x_mel,function(x)x^5.5),y_mel)
axis(1, at=xtick,cex=0.1, labels=c(expression(x[1]),expression(x[2]),expression(paste(psi,x[1],' + (1 -',psi,')',x[2]))),las=0,col='azure4')
axis(2, at=ytick, labels=c( expression(paste('f(',x[1],')')),expression(paste('f(',x[2],')')),'',''),las=1,col='azure4',cex=0.5)

segments(x0=x_1, y0=-200, x1 =x_1, y1 =as.numeric(ytick[1]), lty=3, col='azure4')
segments(x0=x_2, y0=-200, x1 =x_2, y1 =as.numeric(ytick[2]), lty=3, col='azure4')
segments(x0=x_mel, y0=-200, x1 =x_mel, y1 =y_mel, lty=3, col='dimgrey')

segments(x0=0, y0=as.numeric(ytick[1]), x1 =x_1, y1 =as.numeric(ytick[1]), lty=3, col='azure4')
segments(x0=0, y0=as.numeric(ytick[2]), x1 =x_2, y1 =as.numeric(ytick[2]), lty=3, col='azure4')
segments(x0=0, y0=as.numeric(ytick[3]), x1 =x_mel, y1 =as.numeric(ytick[3]), lty=3, col='dimgrey')
segments(x0=0, y0=y_mel, x1 =x_mel, y1 =y_mel, lty=3, col='dimgrey')

segments(x0=x_1, y0=as.numeric(ytick[1]), x1 =x_2, y1 =as.numeric(ytick[2]), lty=4, col='red')

text(1.47,as.numeric(ytick[3])+75 , labels = expression(paste(f,'( ',psi,x[1],'+ (1 -',psi,')',x[2],' )')),cex=1,col='dimgrey')
text(1.49,y_mel+75 , labels = expression(paste(psi,f,'(',x[1],') + (1 -',psi,')',f,'(',x[2],')')),cex=1,col='dimgrey')

points(x_mel,as.numeric(ytick[3]),pch = 19,col='grey',cex=0.7)
points(x_mel,y_mel,pch = 19,col='grey',cex=0.7)



#Graphe: fonction convexe#
#Note: Effacer l'historique de la section '#Graphe: Fonction concave#'
x_1<-1.05
x_2<-3.7
psi<-0.6
x_mel<-psi*x_1+(1-psi)*x_2
y_mel<-psi*(-x_1^-1.5)+(1-psi)*(-x_2^-1.5)

axe_d_x<-seq(1,4,0.05)
plot(axe_d_x,lapply(axe_d_x,function(x)-x^-1.5), type='l',xaxt = "n",yaxt = "n",xlab = '',ylab='',xlim=c(0.3,4))
xtick<-c(x_1,x_2,x_mel)
ytick<-c(lapply(x_1,function(x)-x^-1.5),lapply(x_2,function(x)-x^-1.5),lapply(x_mel,function(x)-x^-1.5),y_mel)
axis(1, at=xtick,cex=0.1, labels=c(expression(x[1]),expression(x[2]),expression(paste(psi,x[1],' + (1 -',psi,')',x[2]))),las=0,col='azure4')
axis(2, at=ytick, labels=c( expression(paste('f(',x[1],')')),expression(paste('f(',x[2],')')),'',''),las=1,col='azure4',cex=0.5)

segments(x0=x_1, y0=-525, x1 =x_1, y1 =as.numeric(ytick[1]), lty=3, col='azure4')
segments(x0=x_2, y0=-525, x1 =x_2, y1 =as.numeric(ytick[2]), lty=3, col='azure4')
segments(x0=x_mel, y0=-525, x1 =x_mel, y1 =as.numeric(ytick[3]), lty=3, col='dimgrey')

segments(x0=0, y0=as.numeric(ytick[1]), x1 =x_1, y1 =as.numeric(ytick[1]), lty=3, col='azure4')
segments(x0=0, y0=as.numeric(ytick[2]), x1 =x_2, y1 =as.numeric(ytick[2]), lty=3, col='azure4')
segments(x0=0, y0=as.numeric(ytick[3]), x1 =x_mel, y1 =as.numeric(ytick[3]), lty=3, col='dimgrey')
segments(x0=0, y0=y_mel, x1 =x_mel, y1 =y_mel, lty=3, col='dimgrey')

segments(x0=x_1, y0=as.numeric(ytick[1]), x1 =x_2, y1 =as.numeric(ytick[2]), lty=4, col='red')

text(0.75,as.numeric(ytick[3])+0.025, labels = expression(paste(f,'( ',psi,x[1],'+ (1 -',psi,')',x[2],' )')),cex=1,col='dimgrey')
text(0.75,y_mel+0.025 , labels = expression(paste(psi,f,'(',x[1],') + (1 -',psi,')',f,'(',x[2],')')),cex=1,col='dimgrey')

points(x_mel,as.numeric(ytick[3]),pch = 19,col='grey',cex=0.7)
points(x_mel,y_mel,pch = 19,col='grey',cex=0.7)



### Graphes: optimisation convexe, relaxation lagrangienne ###
# Fonction de contrainte convexe #

axe_d_x<-seq(-2,2,0.01)
plot(axe_d_x,lapply(axe_d_x,function(x)-0.1*(x-5)^(2)-1),type='l',xlab =expression(x),ylab='',ylim=c(-6,2))
#segments(x0=-2, y0=0, x1 =2, y1 =0, lty=3, col='dimgrey')
lines(axe_d_x,lapply(axe_d_x,function(x)2*x^2-1),col='gray0',lty=2)
lines(axe_d_x,lapply(axe_d_x,function(x)-0.1*(x-5)^(2)-1-(0.1)*(2*x^2-1)),col='gray0',lty=3)
lines(axe_d_x,lapply(axe_d_x,function(x)-0.1*(x-5)^(2)-1-(0.2)*(2*x^2-1)),col='gray0',lty=3)
lines(axe_d_x,lapply(axe_d_x,function(x)-0.1*(x-5)^(2)-1-(0.3)*(2*x^2-1)),col='gray0',lty=3)
lines(axe_d_x,lapply(axe_d_x,function(x)-0.1*(x-5)^(2)-1-(0.4)*(2*x^2-1)),col='gray0',lty=3)
lines(axe_d_x,lapply(axe_d_x,function(x)-0.1*(x-5)^(2)-1-(0.5)*(2*x^2-1)),col='gray0',lty=3)
lines(axe_d_x,lapply(axe_d_x,function(x)-0.1*(x-5)^(2)-1-(0.6)*(2*x^2-1)),col='gray0',lty=3)
lines(axe_d_x,lapply(axe_d_x,function(x)-0.1*(x-5)^(2)-1-(0.7)*(2*x^2-1)),col='gray0',lty=3)
lines(axe_d_x,lapply(axe_d_x,function(x)-0.1*(x-5)^(2)-1-(0.8)*(2*x^2-1)),col='gray0',lty=3)
lines(axe_d_x,lapply(axe_d_x,function(x)-0.1*(x-5)^(2)-1-(1)*(2*x^2-1)),col='gray0',lty=3)

points(sqrt(0.5),lapply(sqrt(0.5),function(x)-0.1*(x-5)^(2)-1),pch = 19,col='red',cex=0.7)
abline(v=sqrt(0.5),col='darkgrey',lty=3)
abline(v=-sqrt(0.5),col='darkgrey',lty=3)

plot(seq(0.1,1,0.01),lapply(seq(0.1,1,0.01),function(y)max(as.numeric(lapply(axe_d_x,function(x)-0.1*(x-5)^(2)-1-(y)*(2*x^2-1))))),type='l',xlab = expression(lambda),ylab=expression(paste('J(',lambda,')')))
abline(h=lapply(sqrt(0.5),function(x)-0.1*(x-5)^(2)-1),col='red',lty=2)
points(0.3,-2.842,pch = 19,col='azure4',cex=0.7)
text(0.33,-2.82,labels=expression(j^'*'),col='azure4')
text(1,-2.82,labels=expression(c^'*'),col='red')


# Fonction de contrainte non-convexe #
axe_d_x<-seq(-2,4.5,0.01)
plot(axe_d_x,lapply(axe_d_x,function(x)-0.01*(3*x^4-17*x^3+17*x^2+48)-2),type='l',xlab =expression(x),ylab='',ylim=c(-4,1))
#segments(x0=-2, y0=0, x1 =2, y1 =0, lty=3, col='dimgrey')
lines(axe_d_x,lapply(axe_d_x,function(x)0.5*(x-0.5)^2-0.5),col='gray0',lty=2)
lines(axe_d_x,lapply(axe_d_x,function(x)-0.01*(3*x^4-17*x^3+17*x^2+48)-2-(0.1)*(0.5*(x-0.5)^2-0.5)),col='gray0',lty=3)
lines(axe_d_x,lapply(axe_d_x,function(x)-0.01*(3*x^4-17*x^3+17*x^2+48)-2-(0.2)*(0.5*(x-0.5)^2-0.5)),col='gray0',lty=3)
lines(axe_d_x,lapply(axe_d_x,function(x)-0.01*(3*x^4-17*x^3+17*x^2+48)-2-(0.3)*(0.5*(x-0.5)^2-0.5)),col='gray0',lty=3)
lines(axe_d_x,lapply(axe_d_x,function(x)-0.01*(3*x^4-17*x^3+17*x^2+48)-2-(0.4)*(0.5*(x-0.5)^2-0.5)),col='gray0',lty=3)
lines(axe_d_x,lapply(axe_d_x,function(x)-0.01*(3*x^4-17*x^3+17*x^2+48)-2-(0.5)*(0.5*(x-0.5)^2-0.5)),col='gray0',lty=3)
lines(axe_d_x,lapply(axe_d_x,function(x)-0.01*(3*x^4-17*x^3+17*x^2+48)-2-(0.6)*(0.5*(x-0.5)^2-0.5)),col='gray0',lty=3)
lines(axe_d_x,lapply(axe_d_x,function(x)-0.01*(3*x^4-17*x^3+17*x^2+48)-2-(0.7)*(0.5*(x-0.5)^2-0.5)),col='gray0',lty=3)
lines(axe_d_x,lapply(axe_d_x,function(x)-0.01*(3*x^4-17*x^3+17*x^2+48)-2-(0.8)*(0.5*(x-0.5)^2-0.5)),col='gray0',lty=3)
lines(axe_d_x,lapply(axe_d_x,function(x)-0.01*(3*x^4-17*x^3+17*x^2+48)-2-(0.9)*(0.5*(x-0.5)^2-0.5)),col='gray0',lty=3)

points(1.499,lapply(1.499,function(y)-0.01*(3*y^4-17*y^3+17*y^2+48)-2),pch = 19,col='red',cex=0.7)
abline(v=1.4999,col='darkgrey',lty=3)
abline(v=-0.4999,col='darkgrey',lty=3)

plot(seq(0.1,1,0.01),lapply(seq(0.1,1,0.01),function(y)max(as.numeric(lapply(axe_d_x,function(x)-0.01*(3*x^4-17*x^3+17*x^2+48)-2-(y)*(0.5*(x-0.5)^2-0.5))))),type='l',xlab = expression(lambda),ylab=expression(paste('J(',lambda,')')),ylim=c(-2.5,-2))
abline(h=lapply(1.499,function(x)-0.01*(3*x^4-17*x^3+17*x^2+48)-2),col='red',lty=2)
points(0.22,-2.3849,pch = 19,col='azure4',cex=0.7) #-2.384961
text(0.25,-2.405,labels=expression(j^'*'),col='azure4')
text(1,-2.46,labels=expression(c^'*'),col='red')



#Pour le mémoire: Relation entre lambda optimal et valeur initial du ptf sans call
axe_x3_1<-seq(0.0005,0.2,0.0001)
par(mfrow=c(1,1))
plot(lapply(axe_x3_1,function(x) lambda_cas_de_base(xini=x,gamma=2, Maturite=10,alpha=0.04, r=0.02, sigma=0.2)
),axe_x3_1,xlab=expression(paste(x["ini"])),ylab=expression(paste(lambda,"*")),type="l",xlim =c(0,50000),las=1 )
lines(lapply(axe_x3_1,function(x) lambda_cas_de_base(xini=x,gamma=3, Maturite=10,alpha=0.04, r=0.02, sigma=0.2)
),axe_x3_1,col='red',lty=2)
lines(lapply(axe_x3_1,function(x) lambda_cas_de_base(xini=x,gamma=4, Maturite=10,alpha=0.04, r=0.02, sigma=0.2)
),axe_x3_1,col='blueviolet',lty=3)
lines(lapply(axe_x3_1,function(x) lambda_cas_de_base(xini=x,gamma=5, Maturite=10,alpha=0.04, r=0.02, sigma=0.2)
),axe_x3_1,col='brown4',lty=4)
lines(lapply(axe_x3_1,function(x) lambda_cas_de_base(xini=x,gamma=6, Maturite=10,alpha=0.04, r=0.02, sigma=0.2)
),axe_x3_1,col='blue',lty=5)
legend(2000000,0.17, legend=c( expression(paste(gamma,"=2")), expression(paste(gamma,"=3")),expression(paste(gamma,"=4")),expression(paste(gamma,"=5")),expression(paste(gamma,"=6"))),
       col=c("black","red","blueviolet",'brown4','blue'), lty=c(1,2,3,4,5), cex=0.8)


##### -TEST SECTION- #####
plot(seq(-2,5,0.01),lapply(seq(-2,5,0.01),function(x)0.01*(3*x^4-17*x^3+20*x^2+48)),col='gray0',lty=2,type='l')
f <- function(y) max(as.numeric(lapply(axe_d_x,function(x)-0.01*(3*x^4-17*x^3+17*x^2+48)-2-(y)*(0.5*(x-0.5)^2-0.5))))+2.84      #-Finding the fee
uniroot(f, c(0.1,3))$root

max(as.numeric(lapply(axe_d_x,function(x)-0.1*(x-5)^(2)-1-(0.1)*(2*x^2-1))))

plot(1:10, 1:10, type="l", lty=2, lwd=3)
plot(1:10, 1:10, type="l", lty=1)
plot(axex,axey,ylab = "Résultat de l'option",type='l',lty=3,xlab="Valeur terminale de X",yaxt = "n" ,xaxt = "n")
xtick<-c(0, b=30)

pnorm(0)