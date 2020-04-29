############ MÉMOIRE: Maitrise en Mathématiques Finacières et Actuarielle  #################



# -SECTION 2: Power Utility -#
#-Note: toutes les formules sont UNIQUEMENT applicables à la fonction de Power Utility
#ex: le Lagrangien est "hard codé" pour Power Utility
library("ggplot2")

########## PARAMÈTRES ##########
Maturi<-10        #Time 'til maturity
r<-0.02      #Risk free rate
alpha<-0.04  #Risky rate
sigma<-0.2   #Volatility
gamma<-2     #Parameter of Utility function
S_0<-1       #Initial value of the asset (S_0>0)
budget<-1  #Initial Budget amount
################################

Ptf_optimal<-function(S_t_par,S_0_par,budget_par,alpha_par,r_par,sigma_par,T_par,gamma_par){
  theta<-(alpha_par-r_par)/sigma_par
  lambda_par<-(1/budget_par)^(gamma_par)*exp((-r_par*gamma_par+r_par-(0.5*theta^2)+theta^2/(2*gamma_par))*T_par)
  resultat<-(lambda_par*(S_t_par/S_0_par)^(-theta/sigma_par)*exp(T_par*(theta*alpha_par/sigma_par-theta*sigma_par/2-r_par-0.5*theta^2)))^(-1/gamma_par)
  return(resultat)
}


Ptf_optimal_fee<-function(S_t_par,S_0_par,budget_par,alpha_par,r_par,sigma_par,T_par,gamma_par,c_f_par,c_s_par){
  alpha_par_1<-alpha_par-c_f_par-c_s_par
  r_par_1<-r_par-c_f_par
  theta<-(alpha_par_1-r_par_1)/sigma_par
  S_t_tilde<-S_t_par*exp(-(c_f_par+c_s_par)*T_par)
  
  lambda_par<-(1/budget_par)^(gamma_par)*exp((-r_par_1*gamma_par+r_par_1-(0.5*theta^2)+theta^2/(2*gamma_par))*T_par)
  resultat<-(lambda_par*(S_t_tilde/S_0_par)^(-theta/sigma_par)*exp(T_par*(theta*alpha_par_1/sigma_par-theta*sigma_par/2-r_par_1-0.5*theta^2)))^(-1/gamma_par)
  return(resultat)
}


# 
# f <- function(x) Ptf_optimal(S_t_par=x,S_0_par=S_0,budget_par=budget,alpha_par=alpha,r_par=r,sigma_par=sigma,T_par=T,gamma_par=2)-Ptf_optimal(S_t_par=x,S_0_par=S_0,budget_par=budget,alpha_par=alpha,r_par=r,sigma_par=sigma,T_par=T,gamma_par=4)
# S_Tcommun<-uniroot(f, c(0.5,1.5))$root 



Shape_lambda<-function(budget_1,gamma_1,r_1,alpha_1,sigma_1,T_1){
  theta1<-(alpha_1-r_1)/sigma_1
  lambda_1<-(1/budget_1)^(gamma_1)*exp((-r_1*gamma_1+r-(0.5*theta1^2)+theta1^2/(2*gamma_1))*T_1)
  return(lambda_1)
}
Shape_lambda(budget_1=budget,gamma_1=3,r_1=r,alpha_1=alpha,sigma_1=sigma,T_1=Maturi)


Ptf_optimal_wrt_qsi<-function(qsi,lambda_par,gamma_par){
  resultat<-(qsi*lambda_par)^(-1/gamma_par)
return(resultat)
}

Ptf_optimalCall_wrt_qsi<-function(lambda_par,a_par,gamma_par,K_par,b_par,qsi){
  
  resultat<-((lambda_par*qsi/a_par)^(-1/gamma_par)-K_par)/alpha_par+b_par
  if(resultat<b_par){resultat1<-0}
  else{resultat1<-resultat}
  return(resultat1)
}



Ptf_optimalCall_wrt_qsi2<-function(alpha_par,r_par,sigma_par,gamma_par,K_par,b_par,a_par,T_par,X_0_par,qsi)
{
  lambda_par<-lambda_optimal(r_FL=r_par,alpha_FL=alpha_par,sigma_FL=sigma_par,gamma_FL=gamma_par,K_FL=K_par,b_FL=b_par,a_FL=a_par,X_0_FL=X_0_par,T_FL=T_par)
  x_concavification<-Find_x_theta_PU(b_par,a_par,K_par,gamma_par)
    
  resultat<-((lambda_par*qsi/a_par)^(-1/gamma_par)-K_par)/a_par+b_par
  
  if(resultat<x_concavification){resultat1<-0}
  else{resultat1<-resultat}
  
  print(c(x_concavification,lambda_par))
  
  return(resultat1)
}


Ptf_optimalCall_wrt_qsi3<-function(alpha_par,r_par,sigma_par,gamma_par,K_par,b_par,a_par,T_par,qsi,lambda_par,x_conca_par)
{
  resultat<-(lambda_par*qsi/a_par)^(-1/gamma_par)
  z_conca_par<-(x_conca_par)^(-gamma_par)
  
  if(qsi*lambda_par<z_conca_par){resultat1<-resultat}
  else{resultat1<-0}
  
  #print(c(x_conca_par,lambda_par))
  
  return(resultat1)
}



Find_x_theta_PU<-function(b_o,a_o,K_o,gamma_o){
  f <- function(x) payoff_call(X_T=x,b=b_o,a=a_o,K=K_o)^(-gamma_o)*x*a_o-P_Utility(payoff_call(x,b=b_o,a=a_o,K=K_o),gamma_o)+P_Utility(K_o,gamma_o)        #-Finding the fee
  result<-uniroot(f, c(b_o,b_o+100))$root
  return(result)
}

Find_x_theta_PU(b_o=1,a_o=1,K_o=1,gamma_o=2)
#lambda_optimal(r_FL=0.03,alpha_FL=0.07,sigma_FL=0.3,gamma_FL=2,K_FL=1,b_FL=1,a_FL=1,X_0_FL=1,T_FL=10)
  
Ptf_optimalCall_wrt_ST<-function(alpha_par,r_par,sigma_par,gamma_par,K_par,b_par,a_par,S_T_par,T_par,X_0_par){
  #Pour S_0=1
  theta<-(alpha_par-r_par)/sigma_par
  lambda_par<-lambda_optimal(r_FL=r_par,alpha_FL=alpha_par,sigma_FL=sigma_par,gamma_FL=gamma_par,K_FL=K_par,b_FL=b_par,a_FL=a_par,X_0_FL=X_0_par,T_FL=T_par)

  
  qsi<-(S_T_par)^(-theta/sigma_par)*exp((theta*alpha_par/sigma_par-0.5*theta*sigma_par-r_par-0.5*theta^2)*T_par)
  indicatrice<-((lambda_par*qsi/a_par)^(-1/gamma_par)-K_par)/a_par+b_par
  x_concavification<-Find_x_theta_PU(b_par,a_par,K_par,gamma_par)
  
  if(indicatrice>x_concavification){X_opt<-indicatrice}
  else{X_opt<-0}
  
  #print(c(x_concavification,lambda_par))
  
  return(X_opt)
}

#Ptf_optimalCall_wrt_ST(alpha_par=alpha,r_par=r,sigma_par=sigma,gamma_par=2,K_par=1,b_par=1,a_par=2,S_T_par=x,T_par=Maturi,X_0_par=budget)
#(r_FL=0.02,alpha_FL=0.04,sigma_FL=0.02,gamma_FL=3,K_FL=1,b_FL=1,a_FL=1,X_0_FL=0.1,T_FL=10)

## Evaluating the term in the indicator function ##
constraint<-function(lambda_p,x_tilde,b_p,a_p,K_p,gamma_p,theta_p,T_p){
  y<--(log(a_p/lambda_p*((x_tilde-b_p)*a_p+K_p)^(-gamma_p)-(r-(theta_p^2)/2)*T_p/theta_p))/theta_p
  return(y)
}
constraint(lambda_p=0.1,x_tilde=Find_x_theta_PU(b_o=1,a_o=1,K_o=1,gamma_o=2)
           ,b_p=1,a_p=1,K_p=1,gamma_p=2,theta_p=1/7.5,T_p=10)
#(r_FL=0.03,alpha_FL=0.07,sigma_FL=0.3,gamma_FL=2,K_FL=1,b_FL=1,a_FL=1,T_FL=10,lambda_FL =1)



## Finding the Langrangian when optimizing the payoff of a call ##
Find_lambda<-function(r_FL,alpha_FL,sigma_FL,gamma_FL,K_FL,b_FL,a_FL,T_FL,lambda_FL){
  theta_FL<-(alpha_FL-r_FL)/sigma_FL
  x_tilde_FL<-Find_x_theta_PU(b_FL,a_FL,K_FL,gamma_FL)
  indi<--log(a_FL/lambda_FL*(a_FL*(x_tilde_FL-b_FL)+K_FL)^(-gamma_FL))
  
  tempo<-(r_FL+0.5*theta_FL^2)*T_FL
  pro_1<-pnorm((tempo-T_FL*(1-1/gamma_FL)*theta_FL^2-indi)/(theta_FL*sqrt(T_FL)))
  pro_2<-pnorm((tempo-T_FL*theta_FL^2-indi)/(theta_FL*sqrt(T_FL)))
    
    FL_1<-lambda_FL^(-1/gamma_FL)/(a_FL^(1-1/gamma_FL))*pro_1*exp(-tempo*(1-1/gamma_FL)+0.5*T_FL*((1-1/gamma_FL)*theta_FL)^2)
    FL_2<-(b_FL-K_FL/a_FL)*pro_2*exp(-r_FL*T_FL)
    
    X_0<-FL_1+FL_2 
   return(X_0)
}
Find_lambda(r_FL=0.02,alpha_FL=0.04,sigma_FL=0.2,gamma_FL=4,K_FL=1,b_FL=1,a_FL=1,T_FL=10,lambda_FL =0.1538843)



lambda_optimal<-function(r_FLo,alpha_FLo,sigma_FLo,gamma_FLo,K_FLo,b_FLo,a_FLo,X_0_FLo,T_FLo){
  
  f <- function(x) Find_lambda(r_FL=r_FLo,alpha_FL=alpha_FLo,sigma_FL=sigma_FLo,gamma_FL=gamma_FLo,K_FL=K_FLo,b_FL=b_FLo,a_FL=a_FLo,T_FL=T_FLo,lambda_FL = x)-X_0_FLo

  lambda_result<-uniroot(f, c(0.0000000000000001,100))$root  
  return(lambda_result)
}
Ptf_optimalCall_wrt_ST(alpha_par=0.04,r_par=0.02,sigma_par=0.2,gamma_par=4,K_par=1,b_par=1,a_par=1,S_T_par=1,T_par=10,X_0_par=1)
lambda_optimal(r_FL=0.02,alpha_FL=0.04,sigma_FL=0.2,gamma_FL=4,K_FL=1,b_FL=1,a_FL=1,X_0_FL=1,T_FL=10)
#lambda=0.0015,S_t=x,S_0=10,a=300,sigma=0.35,T_t=10,alpha=0.07,gamma=3,K=10,r=0.03,b=20
#plot(axe_x4,lapply(axe_x4, function(x)Test_Valeur_Optimale(lambda=0.001,S_t=x,S_0=1,a=1,sigma=0.3,T_t=10,alpha=0.07,gamma=2,K=1,r=0.03,b=1)


Test_Valeur_Optimale<-function(lambda,S_t,S_0,a,sigma,T_t,alpha,gamma,K,r,b){
  theta<-(alpha-r)/sigma
  premier_morc<-(lambda*(S_t/S_0)^(-theta/sigma)*exp(T_t*(theta*alpha/sigma-theta*sigma/2-r-0.5*theta^2)))^(-1/gamma)
  x_op<-Find_x_theta_PU(b_o=b,a_o=a,K_o=K,gamma_o=gamma)
    
  ind<-(premier_morc/(a^(-1/gamma))-K)/a+b
  
  if(ind>x_op){X_fin<-ind}
  else{X_fin<-0}
  
  return(X_fin)
}

Test_Valeur_Optimale(lambda=0.001,S_t=25,S_0=1,a=300,sigma=0.35,T_t=10,alpha=0.07,gamma=3,K=10,r=0.03,b=20)


Xi_shape<-function(S_t_xi,S_0_xi,alpha_xi,r_xi,sigma_xi,T_xi){
  theta<-(alpha_xi-r_xi)/sigma_xi
  xi_final<-(S_t_xi/S_0_xi)^(-theta/sigma_xi)*exp(T_xi*(theta*alpha_xi/sigma_xi-theta*sigma_xi/2-r_xi-0.5*theta^2))
  return(xi_final)
}
Xi_shape(S_t_xi=100,S_0=50,alpha_xi=0.07,r_xi=0.03,sigma_xi=0.35,T_xi=10)

log_graphic<-function(r_FL,alpha_FL,sigma_FL,gamma_FL,K_FL,b_FL,a_FL,T_FL,lambda_FL){
  return(log(Find_lambda(r_FL,alpha_FL,sigma_FL,gamma_FL,K_FL,b_FL,a_FL,T_FL,lambda_FL)))
}


############### Variable Annuities #################

## Finding the Langrangian when optimizing the payoff of a call ##
Ptf_optimalCall_wrt_ST_Fee<-function(alpha_par_ini,r_par_ini,sigma_par,gamma_par,K_par,b_par,a_par,S_T_par,T_par,X_0_par,c_f,c_s){
  #Pour S_0=1
  r_par<-r_par_ini-c_f
  alpha_par<-alpha_par_ini-c_f-c_s
  S_T_par_2<-S_T_par*exp((-c_f-c_s)*T_par) #S_T avec frais à partir de S_T sans les frais.
  
  theta<-(alpha_par-r_par)/sigma_par
  lambda_par<-lambda_optimal(r_FL=r_par,alpha_FL=alpha_par,sigma_FL=sigma_par,gamma_FL=gamma_par,K_FL=K_par,b_FL=b_par,a_FL=a_par,X_0_FL=X_0_par,T_FL=T_par)
  
  
  qsi<-(S_T_par_2)^(-(alpha_par-r_par)/(sigma_par^2))*exp((theta*alpha_par/sigma_par-0.5*theta*sigma_par-r_par-0.5*theta^2)*T_par)
  indicatrice<-((lambda_par*qsi/a_par)^(-1/gamma_par)-K_par)/a_par+b_par
  x_concavification<-Find_x_theta_PU(b_par,a_par,K_par,gamma_par)
  
  if(indicatrice>x_concavification){X_opt<-indicatrice}
  else{X_opt<-0}
  
  #print(c(x_concavification,lambda_par))
  
  return(X_opt)
}
Ptf_optimalCall_wrt_ST_Fee(alpha_par=0.04,r_par=0.02,sigma_par=0.2,gamma_par=2,K_par=1,b_par=1,a_par=1,S_T_par=1,T_par=10,X_0_par=1,c_f=0.005,c_s=0.005)










#Calculating the expected utility on a Call Option
Expected_Utility_C<-function(X_T_par,gamma_par,N_Simul,B_par){
  temp1<-max(X_T_par,B_par)
  temp2<-1/(1-gamma_par)*(temp1)^(1-gamma_par)
  EU<-temp2/N_Simul
  return(EU)
}


################### -GRAPHS- ######################
axe_x<-seq(1,40,0.5)
plot(axe_x,lapply(axe_x,function(x) Ptf_optimal(S_t_par=x,S_0_par=S_0,budget_par=budget,alpha_par=alpha,r_par=r,sigma_par=sigma,T_par=T,gamma_par=9)),xlab=expression(S[T]),ylab=expression(paste("X","*"[T])),type="l",xaxs="i",yaxt = "n" ,xaxt = "n",ylim =c(budget,budget*2.2))
lines(axe_x,lapply(axe_x,function(x) Ptf_optimal(S_t_par=x,S_0_par=S_0,budget_par=budget,alpha_par=alpha,r_par=r,sigma_par=sigma,T_par=T,gamma_par=3)),col='gray')
lines(axe_x,lapply(axe_x,function(x) Ptf_optimal(S_t_par=x,S_0_par=S_0,budget_par=budget,alpha_par=alpha,r_par=r,sigma_par=sigma,T_par=T,gamma_par=5)),col='orange')
lines(axe_x,lapply(axe_x,function(x) Ptf_optimal(S_t_par=x,S_0_par=S_0,budget_par=budget,alpha_par=alpha,r_par=r,sigma_par=sigma,T_par=T,gamma_par=7)),col='red')
legend(27, budget+3, legend=c( expression(paste(gamma,"=3")), expression(paste(gamma,"=5")),expression(paste(gamma,"=7")),expression(paste(gamma,"=9"))),
       col=c("gray","orange", "red","black"), lty=c(1,1,1,1), cex=0.8)

#Graphique-retouches
axe_x<-seq(0.1,2,0.001)
plot(axe_x,lapply(axe_x,function(x) Ptf_optimal(S_t_par=x,S_0_par=S_0,budget_par=budget,alpha_par=alpha,r_par=r,sigma_par=sigma,T_par=T,gamma_par=3)),xlab=expression(S[T]),ylab=expression(paste("X","*"[T])),type="l",xaxs="i",yaxt = "n",xaxt = "n")
lines(axe_x,lapply(axe_x,function(x) Ptf_optimal(S_t_par=x,S_0_par=S_0,budget_par=budget,alpha_par=alpha,r_par=r,sigma_par=sigma,T_par=T,gamma_par=3.5)),col='gray',lty=2)
lines(axe_x,lapply(axe_x,function(x) Ptf_optimal(S_t_par=x,S_0_par=S_0,budget_par=budget,alpha_par=alpha,r_par=r,sigma_par=sigma,T_par=T,gamma_par=4)),col='orange',lty=4)
lines(axe_x,lapply(axe_x,function(x) Ptf_optimal(S_t_par=x,S_0_par=S_0,budget_par=budget,alpha_par=alpha,r_par=r,sigma_par=sigma,T_par=T,gamma_par=4.5)),col='red',lty=3)
legend(1.5, 1, legend=c( expression(paste(gamma,"=3")), expression(paste(gamma,"=3.5")),expression(paste(gamma,"=4")),expression(paste(gamma,"=4.5"))),
       col=c("black","grey", "orange","red"), lty=c(1,1,1,1), cex=0.8)


par(mfrow=c(1,2))
axe_x1<-seq(0.5,100,1)
plot(axe_x1,lapply(axe_x1,function(x) Ptf_optimal_wrt_qsi(qsi=x,lambda_par=0.10,gamma_par=3)),ylab=expression(paste("X","*"[T])),xlab=expression(xi[T]),type="l",xaxs="i",col='gray',main='Optimisation de X')
lines(axe_x1,lapply(axe_x1,function(x) Ptf_optimal_wrt_qsi(qsi=x,lambda_par=0.10,gamma_par=7)))
lines(axe_x1,lapply(axe_x1,function(x) Ptf_optimal_wrt_qsi(qsi=x,lambda_par=0.10,gamma_par=4)),col='orange')
lines(axe_x1,lapply(axe_x1,function(x) Ptf_optimal_wrt_qsi(qsi=x,lambda_par=0.10,gamma_par=5)),col='red')
legend(60, 2, legend=c( expression(paste(gamma,"=2")), expression(paste(gamma,"=3")),expression(paste(gamma,"=4")),expression(paste(gamma,"=5"))),
       col=c("black", "gray","orange","red"), lty=c(1,1,1,1), cex=0.8)

plot(axe_x1,lapply(axe_x1,function(x)Ptf_optimalCall_wrt_qsi(lambda_par=0.10,alpha_par=5,gamma_par=3,K_par=1,b_par=2,qsi=x)),ylab=expression(paste("X","*"[T])),xlab=expression(xi[T]),type="l",col='gray',xaxs="i",main='Optimisation d une option d achat')
lines(axe_x1,lapply(axe_x1,function(x) Ptf_optimalCall_wrt_qsi(lambda_par=0.10,alpha_par=5,gamma_par=7,K_par=1,b_par=2,qsi=x)))
lines(axe_x1,lapply(axe_x1,function(x) Ptf_optimalCall_wrt_qsi(lambda_par=0.10,alpha_par=5,gamma_par=4,K_par=1,b_par=2,qsi=x)),col='orange')
lines(axe_x1,lapply(axe_x1,function(x) Ptf_optimalCall_wrt_qsi(lambda_par=0.10,alpha_par=5,gamma_par=5,K_par=1,b_par=2,qsi=x)),col='red')
legend(60, 2, legend=c( expression(paste(gamma,"=2")), expression(paste(gamma,"=3")),expression(paste(gamma,"=4")),expression(paste(gamma,"=5"))),
       col=c("black", "gray","orange","red"), lty=c(1,1,1,1), cex=0.8)

axe_x11<-seq(0.1,2,0.1)
plot(axe_x11,lapply(axe_x11,function(x)Ptf_optimalCall_wrt_qsi(lambda_par=x,alpha_par=5,gamma_par=3,K_par=1,b_par=2,qsi=0.2)),ylab=expression(paste("X","*"[T])),xlab=expression(xi[T]),type="l",col='gray',xaxs="i",main='Optimisation d une option d achat')


#- Comparison Optimal Ptf wrt S_t CALL VS ORIGINAL-#
par(mfrow=c(1,3))
axe_x2<-seq(0.01,0.5,0.01)
plot(axe_x2,lapply(axe_x2,function(x) Ptf_optimalCall_wrt_ST(alpha_par=alpha,r_par=r,sigma_par=sigma,gamma_par=2,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=Maturi,X_0_par=budget)),ylab=expression(paste("X","*"[T])),xlab=expression(S[T]),type="l",xaxs="i",main='Optimisation Resultat Call',ylim=c(0,40))
lines(axe_x2,lapply(axe_x2,function(x) Ptf_optimalCall_wrt_ST(alpha_par=alpha,r_par=r,sigma_par=sigma,gamma_par=2.5,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=Maturi,X_0_par=budget)),col='grey')
lines(axe_x2,lapply(axe_x2,function(x) Ptf_optimalCall_wrt_ST(alpha_par=alpha,r_par=r,sigma_par=sigma,gamma_par=3,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=Maturi,X_0_par=budget)),col='orange')
legend(0.5,40, legend=c( expression(paste(gamma,"=2")), expression(paste(gamma,"=2.5")),expression(paste(gamma,"=3"))),
       col=c("black", "grey","orange"), lty=c(1,1,1), cex=0.8)
abline(h=budget, col='blue',lty=2)
axis(1, at=budget, labels=c(expression(X[0])))

plot(axe_x2,lapply(axe_x2,function(x) Ptf_optimal(S_t_par=x,S_0_par=S_0,budget_par=budget,alpha_par=alpha,r_par=r,sigma_par=sigma,T_par=Maturi,gamma_par=2)),xlab=expression(S[T]),ylab=expression(paste("X","*"[T])),type="l",xaxs="i",ylim =c(0,40),main='Optimisation de X',lty=4)
lines(axe_x2,lapply(axe_x2,function(x) Ptf_optimal(S_t_par=x,S_0_par=S_0,budget_par=budget,alpha_par=alpha,r_par=r,sigma_par=sigma,T_par=Maturi,gamma_par=2.5)),col='gray',lty=4)
lines(axe_x2,lapply(axe_x2,function(x) Ptf_optimal(S_t_par=x,S_0_par=S_0,budget_par=budget,alpha_par=alpha,r_par=r,sigma_par=sigma,T_par=Maturi,gamma_par=3)),col='orange',lty=4)
legend(0.5, 40, legend=c( expression(paste(gamma,"=2")), expression(paste(gamma,"=2.5")),expression(paste(gamma,"=3"))),
       col=c("black","grey", "orange"), lty=c(1,1,1), cex=0.8)
abline(h=budget, col='blue',lty=2)
axis(1, at=budget, labels=c(expression(X[0])))

plot(axe_x2,lapply(axe_x2,function(x) Ptf_optimalCall_wrt_ST(alpha_par=alpha,r_par=r,sigma_par=sigma,gamma_par=2,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=Maturi,X_0_par=budget)),ylab=expression(paste("X","*"[T])),xlab=expression(S[T]),type="l",xaxs="i",main='Optimisation Resultat Call',ylim = c(0,40))
lines(axe_x2,lapply(axe_x2,function(x) Ptf_optimalCall_wrt_ST(alpha_par=alpha,r_par=r,sigma_par=sigma,gamma_par=2.5,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=Maturi,X_0_par=budget)),col='grey')
lines(axe_x2,lapply(axe_x2,function(x) Ptf_optimalCall_wrt_ST(alpha_par=alpha,r_par=r,sigma_par=sigma,gamma_par=3,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=Maturi,X_0_par=budget)),col='orange')
legend(0.5,40, legend=c( expression(paste(gamma,"=2")), expression(paste(gamma,"=2.5")),expression(paste(gamma,"=3"))),
       col=c("black", "grey","orange"), lty=c(1,1,1), cex=0.8)
lines(axe_x2,lapply(axe_x2,function(x) Ptf_optimal(S_t_par=x,S_0_par=S_0,budget_par=budget,alpha_par=alpha,r_par=r,sigma_par=sigma,T_par=Maturi,gamma_par=2)),lty=4)
lines(axe_x2,lapply(axe_x2,function(x) Ptf_optimal(S_t_par=x,S_0_par=S_0,budget_par=budget,alpha_par=alpha,r_par=r,sigma_par=sigma,T_par=Maturi,gamma_par=2.5)),col='gray',lty=4)
lines(axe_x2,lapply(axe_x2,function(x) Ptf_optimal(S_t_par=x,S_0_par=S_0,budget_par=budget,alpha_par=alpha,r_par=r,sigma_par=sigma,T_par=Maturi,gamma_par=3)),col='orange',lty=4)
abline(h=budget, col='blue',lty=2)
axis(1, at=budget, labels=c(expression(X[0])))

# GRAPHIC: Lambda according to the initial value X_0
par(mfrow=c(1,2))
axe_x3<-seq(0.001,0.2,0.001)
plot(axe_x3,lapply(axe_x3,function(x) Find_lambda(r_FL=0.02,alpha_FL=0.04,sigma_FL=0.2,gamma_FL=2,K_FL=1,b_FL=1,a_FL=1,T_FL=10,lambda_FL =x)
),ylab=expression(paste(X[0])),xlab=expression(paste(lambda,"*")),type="l",ylim=c(0,300),main=expression(paste("Value of the initial investment in the portfolio given ",lambda,"*")))
lines(axe_x3,lapply(axe_x3,function(x) Find_lambda(r_FL=0.02,alpha_FL=0.04,sigma_FL=0.2,gamma_FL=3,K_FL=1,b_FL=1,a_FL=1,T_FL=10,lambda_FL =x)
),col='red')
lines(axe_x3,lapply(axe_x3,function(x) Find_lambda(r_FL=0.02,alpha_FL=0.04,sigma_FL=0.2,gamma_FL=4,K_FL=1,b_FL=1,a_FL=1,T_FL=10,lambda_FL =x)
),col='orange')
lines(axe_x3,lapply(axe_x3,function(x) Find_lambda(r_FL=0.02,alpha_FL=0.04,sigma_FL=0.2,gamma_FL=5,K_FL=1,b_FL=1,a_FL=1,T_FL=10,lambda_FL =x)
),col='grey')
lines(axe_x3,lapply(axe_x3,function(x) Find_lambda(r_FL=0.02,alpha_FL=0.04,sigma_FL=0.2,gamma_FL=6,K_FL=1,b_FL=1,a_FL=1,T_FL=10,lambda_FL =x)
),col='blue')
legend( 0.07,170, legend=c( expression(paste(gamma,"=2")), expression(paste(gamma,"=3")),expression(paste(gamma,"=4")),expression(paste(gamma,"=5")),expression(paste(gamma,"=6"))),
       col=c("black","red","orange",'grey','blue'), lty=c(1,1,1,1,1), cex=0.8)

plot(axe_x3,lapply(axe_x3,function(x) log_graphic(r_FL=0.02,alpha_FL=0.04,sigma_FL=0.2,gamma_FL=2,K_FL=1,b_FL=1,a_FL=1,T_FL=10,lambda_FL =x)
),ylab=expression(paste("log(",X[0],")")),xlab=expression(paste(lambda,"*")),type="l",ylim=c(-30,10))
lines(axe_x3,lapply(axe_x3,function(x) log_graphic(r_FL=0.02,alpha_FL=0.04,sigma_FL=0.2,gamma_FL=3,K_FL=1,b_FL=1,a_FL=1,T_FL=10,lambda_FL =x)
),col='red')
lines(axe_x3,lapply(axe_x3,function(x) log_graphic(r_FL=0.02,alpha_FL=0.04,sigma_FL=0.2,gamma_FL=4,K_FL=1,b_FL=1,a_FL=1,T_FL=10,lambda_FL =x)
),col='orange')
lines(axe_x3,lapply(axe_x3,function(x) log_graphic(r_FL=0.02,alpha_FL=0.04,sigma_FL=0.2,gamma_FL=5,K_FL=1,b_FL=1,a_FL=1,T_FL=10,lambda_FL =x)
),col='grey')
lines(axe_x3,lapply(axe_x3,function(x) log_graphic(r_FL=0.02,alpha_FL=0.04,sigma_FL=0.2,gamma_FL=6,K_FL=1,b_FL=1,a_FL=1,T_FL=10,lambda_FL =x)
),col='blue')
legend( 0.0007,2, legend=c( expression(paste(gamma,"=2")), expression(paste(gamma,"=3")),expression(paste(gamma,"=4")),expression(paste(gamma,"=5")),expression(paste(gamma,"=6"))),
        col=c("black","red","orange",'grey','blue'), lty=c(1,1,1,1,1), cex=0.8)
#r_FL=0.03,alpha_FL=0.07,sigma_FL=0.05,gamma_FL=2,K_FL=15,b_FL=30,a_FL=300,X_0_FL=1,T_FL=10,lambda_FL =x

#Pour le mémoire: Relation entre lambda optimal et valeur initial du ptf avec call
axe_x3_1<-seq(0.00001,0.2,0.0001)
par(mfrow=c(1,1))
plot(lapply(axe_x3_1,function(x) Find_lambda(r_FL=0.02,alpha_FL=0.04,sigma_FL=0.2,gamma_FL=2,K_FL=1,b_FL=1,a_FL=1,T_FL=10,lambda_FL =x)
),axe_x3_1,xlab=expression(paste(x)),ylab=expression(paste(lambda,"*")),type="l",xlim=c(0,30),las=1)
lines(lapply(axe_x3_1,function(x) Find_lambda(r_FL=0.02,alpha_FL=0.04,sigma_FL=0.2,gamma_FL=3,K_FL=1,b_FL=1,a_FL=1,T_FL=10,lambda_FL =x)
),axe_x3_1,col='red',lty=2)
lines(lapply(axe_x3_1,function(x) Find_lambda(r_FL=0.02,alpha_FL=0.04,sigma_FL=0.2,gamma_FL=4,K_FL=1,b_FL=1,a_FL=1,T_FL=10,lambda_FL =x)
),axe_x3_1,col='blueviolet',lty=3)
lines(lapply(axe_x3_1,function(x) Find_lambda(r_FL=0.02,alpha_FL=0.04,sigma_FL=0.2,gamma_FL=5,K_FL=1,b_FL=1,a_FL=1,T_FL=10,lambda_FL =x)
),axe_x3_1,col='brown4',lty=4)
lines(lapply(axe_x3_1,function(x) Find_lambda(r_FL=0.02,alpha_FL=0.04,sigma_FL=0.2,gamma_FL=6,K_FL=1,b_FL=1,a_FL=1,T_FL=10,lambda_FL =x)
),axe_x3_1,col='blue',lty=5)
legend(20,0.17, legend=c( expression(paste(gamma,"=2")), expression(paste(gamma,"=3")),expression(paste(gamma,"=4")),expression(paste(gamma,"=5")),expression(paste(gamma,"=6"))),
        col=c("black","red","blueviolet",'brown4','blue'), lty=c(1,2,3,4,5), cex=0.8)

#Valeur Optimale ptf VS S_T
par(mfrow=c(1,1))
axe_x4<-seq(0.1,4,0.05)
plot(axe_x4,lapply(axe_x4, function(x)Test_Valeur_Optimale(lambda=0.001,S_t=x,S_0=1,a=1,sigma=0.3,T_t=10,alpha=0.07,gamma=2,K=1,r=0.03,b=1)
),type='l',xlab=expression(s[T]),ylab=expression(paste("X","*"[T])),ylim=c(0,60))
lines(axe_x4,lapply(axe_x4,function(x)Test_Valeur_Optimale(lambda=0.001,S_t=x,S_0=1,a=1,sigma=0.3,T_t=10,alpha=0.07,gamma=2.2,K=1,r=0.03,b=1)),col='grey')
lines(axe_x4,lapply(axe_x4,function(x)Test_Valeur_Optimale(lambda=0.001,S_t=x,S_0=1,a=1,sigma=0.3,T_t=10,alpha=0.07,gamma=2.5,K=1,r=0.03,b=1)),col='blue')
lines(axe_x4,lapply(axe_x4,function(x)Test_Valeur_Optimale(lambda=0.001,S_t=x,S_0=1,a=1,sigma=0.3,T_t=10,alpha=0.07,gamma=3,K=1,r=0.03,b=1)),col='darkblue')
lines(axe_x4,lapply(axe_x4,function(x)Test_Valeur_Optimale(lambda=0.001,S_t=x,S_0=1,a=1,sigma=0.3,T_t=10,alpha=0.07,gamma=3.2,K=1,r=0.03,b=1)),col='coral')

legend( 0.3,50, legend=c( expression(paste(gamma,"=2")), expression(paste(gamma,"=2.2")),expression(paste(gamma,"=2.3")),expression(paste(gamma,"=3")),expression(paste(gamma,"=3.2"))),
        col=c("black",'grey','blue','darkblue','coral'), lty=c(1,1,1,1,1), cex=0.8)



# Xi with respect to S_t #
par(mfrow=c(1,1))
plot(axe_x4,lapply(axe_x4,function(x)Xi_shape(S_t_xi=x,S_0=1,alpha_xi=0.04,r_xi=0.02,sigma_xi=0.2,T_xi=10)),type='l',xlab=expression(S[t]),ylab=expression(xi[t]))
lines(axe_x4,lapply(axe_x4,function(x)Xi_shape(S_t_xi=x,S_0=50,alpha_xi=0.07,r_xi=0.03,sigma_xi=0.35,T_xi=10)),col='grey')
lines(axe_x4,lapply(axe_x4,function(x)Xi_shape(S_t_xi=x,S_0=50,alpha_xi=0.07,r_xi=0.03,sigma_xi=0.55,T_xi=10)),col='orange')
legend( 34,8, legend=c( expression(paste(sigma,"=0.25")), expression(paste(sigma,"=0.35")),expression(paste(sigma,"=0.55"))),
        col=c("black",'grey','orange'), lty=c(1,1,1), cex=0.8)

plot(axe_x4,lapply(axe_x4,function(x)Xi_shape(S_t_xi=x,S_0=50,alpha_xi=0.05,r_xi=0.03,sigma_xi=0.35,T_xi=10)),type='l',xlab=expression(S[t]),ylab=expression(xi[t]))
lines(axe_x4,lapply(axe_x4,function(x)Xi_shape(S_t_xi=x,S_0=50,alpha_xi=0.07,r_xi=0.03,sigma_xi=0.35,T_xi=10)),col='grey')
lines(axe_x4,lapply(axe_x4,function(x)Xi_shape(S_t_xi=x,S_0=50,alpha_xi=0.09,r_xi=0.03,sigma_xi=0.35,T_xi=10)),col='orange')
legend( 34,1.1, legend=c( expression(paste(alpha,"=0.05")), expression(paste(alpha,"=0.07")),expression(paste(alpha,"=0.09"))),
        col=c("black",'grey','orange'), lty=c(1,1,1), cex=0.8)


####################### TESTS SECTION #######################
par(mfrow=c(1,1))
axe_test<-seq(0.1,3,0.001)

###### Graphiques X*_t differents gamma ######
#-Préparer les données-#
VO_gamma2<-lapply(axe_test, function (x)Ptf_optimalCall_wrt_ST(alpha_par=0.04,r_par=0.02,sigma_par=0.2,gamma_par=3,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=10,X_0_par=1))
VO_gamma3<-lapply(axe_test, function (x)Ptf_optimalCall_wrt_ST(alpha_par=0.04,r_par=0.02,sigma_par=0.2,gamma_par=3,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=10,X_0_par=1))
VO_gamma3_5<-lapply(axe_test, function (x)Ptf_optimalCall_wrt_ST(alpha_par=0.04,r_par=0.02,sigma_par=0.2,gamma_par=3.5,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=10,X_0_par=1))
VO_gamma4<-lapply(axe_test, function (x)Ptf_optimalCall_wrt_ST(alpha_par=0.04,r_par=0.02,sigma_par=0.2,gamma_par=4,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=10,X_0_par=1))
VO_gamma4_5<-lapply(axe_test, function (x)Ptf_optimalCall_wrt_ST(alpha_par=0.04,r_par=0.02,sigma_par=0.2,gamma_par=4.5,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=10,X_0_par=1))

F_gamma2<-lapply(axe_test, function (x)Ptf_optimalCall_wrt_ST_Fee(alpha_par=0.04,r_par=0.02,sigma_par=0.2,gamma_par=2,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=10,X_0_par=1,c_f=0.01224,c_s=0.01224))
F_gamma3<-lapply(axe_test, function (x)Ptf_optimalCall_wrt_ST_Fee(alpha_par=0.04,r_par=0.02,sigma_par=0.2,gamma_par=3,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=10,X_0_par=1,c_f=0.01224,c_s=0.01224))
F_gamma4<-lapply(axe_test, function (x)Ptf_optimalCall_wrt_ST_Fee(alpha_par=0.04,r_par=0.02,sigma_par=0.2,gamma_par=4,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=10,X_0_par=1,c_f=0.01224,c_s=0.01224))
F_gamma5<-lapply(axe_test, function (x)Ptf_optimalCall_wrt_ST_Fee(alpha_par=0.04,r_par=0.02,sigma_par=0.2,gamma_par=5,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=10,X_0_par=1,c_f=0.01224,c_s=0.01224))
F_gamma6<-lapply(axe_test, function (x)Ptf_optimalCall_wrt_ST_Fee(alpha_par=0.04,r_par=0.02,sigma_par=0.2,gamma_par=6,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=10,X_0_par=1,c_f=0.01224,c_s=0.01224))
F_gamma7<-lapply(axe_test, function (x)Ptf_optimalCall_wrt_ST_Fee(alpha_par=0.04,r_par=0.02,sigma_par=0.2,gamma_par=7,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=10,X_0_par=1,c_f=0.01224,c_s=0.01224))

ggp_gamma2<-data.frame(results=as.numeric(F_gamma2))
ggp_gamma3<-data.frame(results=as.numeric(F_gamma3))
ggp_gamma4<-data.frame(results=as.numeric(F_gamma4))
ggp_gamma5<-data.frame(results=as.numeric(F_gamma5))
ggp_gamma6<-data.frame(results=as.numeric(F_gamma6))
ggp_gamma7<-data.frame(results=as.numeric(F_gamma7))

ggp_gamma2$gamma<-'2'
ggp_gamma3$gamma<-'3'
ggp_gamma4$gamma<-'4'
ggp_gamma5$gamma<-'5'
ggp_gamma6$gamma<-'6'
ggp_gamma7$gamma<-'7'

F_ggp_total<-rbind(ggp_gamma2,ggp_gamma3,ggp_gamma4,ggp_gamma5,ggp_gamma6,ggp_gamma7)
F_ggp_total$S_T<-axe_test

# 1) Avec plot #
plot(axe_test,VO_gamma3,type='l',xlab=expression(s[T]),ylab=expression(paste(F[T],"*")))
lines(axe_test,VO_gamma3_5,col='grey',lty=2)
lines(axe_test,VO_gamma4,col='orange',lty=4)
lines(axe_test, VO_gamma4_5,col='red',lty=3)
legend( 0.2,1.7, legend=c( expression(paste(gamma,"=3")), expression(paste(gamma,"=3.5")),expression(paste(gamma,"=4")),expression(paste(gamma,"=4.5"))),
        col=c("black","grey",'orange','red'), lty=c(1,2,4,3), cex=0.8)

# 2) Avec ggplot #
ggplot(data=F_ggp_total,aes(x=F_ggp_total$S_T,y=F_ggp_total$results,group=gamma,color=gamma))+
  geom_line(aes(linetype=gamma))+
 scale_linetype_discrete(name="",
                          breaks=c("2", "3", "4","5","6","7"),
                          labels=c(expression(paste(gamma,"=2")), expression(paste(gamma,"=3")), expression(paste(gamma,"=4")),expression(paste(gamma,"=5")),expression(paste(gamma,"=6")),expression(paste(gamma,"=7"))))+
  scale_color_discrete(name="",
                          breaks=c("2", "3", "4","5","6","7"),
                          labels=c(expression(paste(gamma,"=2")), expression(paste(gamma,"=3")), expression(paste(gamma,"=4")),expression(paste(gamma,"=5")),expression(paste(gamma,"=6")),expression(paste(gamma,"=7"))))+
  theme(legend.position = "bottom")+
  labs(x=expression(paste(S[T])), y=expression(paste(F[T]^"*")))









###### Graphiques X*_t differents frais ######

#-Préparer les données-#
XT_gamma4_fees00et0<-lapply(axe_test, function (x)Ptf_optimalCall_wrt_ST_Fee(alpha_par=0.04,r_par=0.02,sigma_par=0.2,gamma_par=3,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=10,X_0_par=1,c_f=0,c_s=0)
)
XT_gamma4_fees18et0648<-lapply(axe_test, function (x)Ptf_optimalCall_wrt_ST_Fee(alpha_par=0.04,r_par=0.02,sigma_par=0.2,gamma_par=3,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=10,X_0_par=1,c_f=0.00648,c_s=0.018)
)
XT_gamma4_fees1224et1224<-lapply(axe_test, function (x)Ptf_optimalCall_wrt_ST_Fee(alpha_par=0.04,r_par=0.02,sigma_par=0.2,gamma_par=3,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=10,X_0_par=1,c_f=0.01224,c_s=0.01224)
)
XT_gamma4_fees0et2448<-lapply(axe_test, function (x)Ptf_optimalCall_wrt_ST_Fee(alpha_par=0.04,r_par=0.02,sigma_par=0.2,gamma_par=3,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=10,X_0_par=1,c_f=0.02448,c_s=0.0)
)
ggp_fees00et0<-data.frame(results=as.numeric(XT_gamma4_fees00et0))
ggp_fees18et0648<-data.frame(results=as.numeric(XT_gamma4_fees18et0648))
ggp_fees1224et1224<-data.frame(results=as.numeric(XT_gamma4_fees1224et1224))
ggp_fees0et2448<-data.frame(results=as.numeric(XT_gamma4_fees0et2448)) 

ggp_fees00et0$fee<-'c_s=0 et c_f=0'
ggp_fees18et0648$fee<-'c_s=1.8 et c_f=0.648'
ggp_fees1224et1224$fee<-'c_s=1.224 et c_f=1.224'
ggp_fees0et2448$fee<-'c_s=0 et c_f=2.448'

Compa_ggp_total<-rbind(ggp_fees00et0,ggp_fees18et0648,ggp_fees1224et1224,ggp_fees0et2448)
Compa_ggp_total$c_s<-axe_test #Ici, c_s représente plutot S_T
Compa_ggp_total$fee<-factor(Compa_ggp_total$fee , levels = c('c_s=0 et c_f=0', 'c_s=1.8 et c_f=0.648','c_s=1.224 et c_f=1.224','c_s=0 et c_f=2.448'))

#-Graphiques-#
# 1) Avec plot #
plot(axe_test,XT_gamma4_fees00et0,col='purple',lty=1,type='l',xlab=expression(s[T]),ylab=expression(paste(F[T],"*")),xlim =c(0,2),ylim=c(0,2) )
lines(axe_test,XT_gamma4_fees18et0648,col='navy',lty=3)
lines(axe_test,XT_gamma4_fees1224et1224,col='blue',lty=5)
lines(axe_test,XT_gamma4_fees0et2448,col='firebrick3',lty=2)
legend( 1.3,1, legend=c(expression(paste("Aucun frais")), expression(paste(c[s],"=1.800%"," et ",c[f],"=0.648%")),expression(paste(c[s],"=1.224%"," et ",c[f],"=1.224%")),expression(paste(c[s],"=0.000%"," et ",c[f],"=2.448%"))),
        col=c('purple','navy','blue','firebrick3'), lty=c(1,3,5,2), cex=0.8)

# 2) Avec ggplot2, pour le mémoire #
ggplot(data=Compa_ggp_total,aes(x=Compa_ggp_total$c_s,y=Compa_ggp_total$results,group=fee,color=fee))+
  geom_line(aes(linetype=fee))+
scale_color_discrete(name="",
                       breaks=c('c_s=0 et c_f=0','c_s=0 et c_f=2.448','c_s=1.224 et c_f=1.224','c_s=1.8 et c_f=0.648'),
                       labels=c(expression(paste("Aucun frais")), expression(atop(paste(c[s],"=0.000%"),paste(c[f],"=2.448%"))),expression(atop(paste(c[s],"=1.224%"),paste(c[f],"=1.224%"))),expression(atop(paste(c[s],"=1.800%"),paste(c[f],"=0.648%")))))+
  scale_linetype_discrete(name="",
                       breaks=c('c_s=0 et c_f=0','c_s=0 et c_f=2.448','c_s=1.224 et c_f=1.224','c_s=1.8 et c_f=0.648'),
                       labels=c(expression(paste("Aucun frais")), expression(atop(paste(c[s],"=0.000%"),paste(c[f],"=2.448%"))),expression(atop(paste(c[s],"=1.224%"),paste(c[f],"=1.224%"))),expression(atop(paste(c[s],"=1.800%"),paste(c[f],"=0.648%")))))+
  theme(legend.position = "bottom")+
  #theme_classic()+
  labs(x=expression(paste(S[T])), y=expression(paste(F[T]^"*")))
#theme(plot.title = element_text(family = "Helvetica", face = "bold", size = (15),hjust = 0.5))
#guides(colour = guide_legend(nrow = 2))
















#############################################################
################# SIMULATION SECTION ########################

################## -Paramètres- ###############
N_simulations<-1000000
alpha_sim<-0.04
r_sim<-0.02
T_sim<-10
sigma_sim<-0.2
gamma_sim<-7
x_conca<-Find_x_theta_PU(1,1,1,gamma_sim)
c_s_sim<-0.0 
c_f_sim<-0.02448-c_s_sim
alpha_sim_frais<-alpha_sim-c_s_sim-c_f_sim
r_sim_frais<-r_sim-c_f_sim
lbd_sim_call<-lambda_optimal(r_FLo=r_sim,alpha_FLo=alpha_sim,sigma_FLo=sigma_sim,gamma_FLo=gamma_sim,K_FLo=1,b_FLo=1,a_FLo=1,X_0_FLo=1,T_FLo=T_sim)
lbd_sim_call_f<-lambda_optimal(r_FLo=r_sim_frais,alpha_FLo=alpha_sim_frais,sigma_FLo=sigma_sim,gamma_FLo=gamma_sim,K_FLo=1,b_FLo=1,a_FLo=1,X_0_FLo=1,T_FLo=T_sim)
#lbd_sim<-  

library("ggplot2")
################################################


rand_vector<-sqrt(T_sim)*rnorm(N_simulations)
S_T_simulated<-lapply(rand_vector,function (x)1*exp(alpha_sim*T_sim-0.5*T_sim*sigma_sim^2+sigma_sim*x))
X_Optimal_simulated<-lapply(S_T_simulated, function(x) Ptf_optimalCall_wrt_ST(alpha_par=alpha_sim,r_par=r_sim,sigma_par=0.2,gamma_par=gamma_sim,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=T_sim,X_0_par=1)
)

hist(as.numeric(X_Optimal_simulated),main='',xlab=expression(paste("X","*"[T])),ylab='# Realisations',col="azure3",xlim=c(0,2),breaks = 100)
abline(v =x_concavification, col = "black", lwd = 2,lty=2)

sum(X_Optimal_simulated==0)  








#Simualtions of S_T_avec frais#

rand_vector2<-rand_vector
S_T_simulated2<-lapply(rand_vector2,function (x)1*exp((alpha_sim-c_s_sim-c_f_sim)*T_sim-0.5*T_sim*sigma_sim^2+sigma_sim*x))
X_Optimal_simulated_Fee<-lapply(S_T_simulated2, function(x) Ptf_optimalCall_wrt_ST_Fee(alpha_par=0.04,r_par=0.02,sigma_par=0.2,gamma_par=4,K_par=1,b_par=1,a_par=1,S_T_par=x,T_par=10,X_0_par=1,c_f=c_f_sim,c_s=c_s_sim)
)

par(new = TRUE)
hist(as.numeric(X_Optimal_simulated_Fee),main='',xlab=expression(paste("X","*"[T])),ylab='',yaxt = "n",col=rgb(0,0,1,alpha=0.2),xlim=c(0,2),breaks = 100)
abline(v =x_concavification, col = "black", lwd = 2,lty=2)
legend( 0.2,19000, legend=c( expression(paste("X","*"[T]," sans frais")),expression(paste("X","*"[T]," avec frais"),'Point de concavification')),
        col=c("grey",rgb(0,0,1,alpha=0.2),"black"), lty=c(1,1,2),lwd = c(4,4,2), cex=0.8)

sum(X_Optimal_simulated_Fee==0) 


Ptf_optimalCall_wrt_ST_Fee(alpha_par=0.04,r_par=0.02,sigma_par=0.2,gamma_par=4,K_par=1,b_par=1,a_par=1,S_T_par=1,T_par=10,X_0_par=1,c_f=c_f_sim,c_s=c_s_sim)
Ptf_optimalCall_wrt_ST(alpha_par=0.04,r_par=0.02,sigma_par=0.2,gamma_par=4,K_par=1,b_par=1,a_par=1,S_T_par=1,T_par=10,X_0_par=1)

1-pnorm((-(((log(1.5874032*0.1538849^(1/4)))*(-4)+(0.02+0.5*0.1^2)*10)/0.1)/sqrt(10)))


##### New Histograms #####


#### SIMULATIONS 2 ########
#Simulations Xi#
rand_vector3<-sqrt(T_sim)*rnorm(N_simulations)

theta_sim<-(alpha_sim-0.02)/0.2
xi_T_simulated<-lapply(rand_vector3,function (x)exp(-(0.02+0.5*theta_sim^2)*10-theta_sim*x))
X_Optimal_simulated_xi<-lapply(xi_T_simulated, function(x) Ptf_optimalCall_wrt_qsi2(alpha_par=alpha_sim,r_par=0.02,sigma_par=0.2,gamma_par=4,K_par=1,b_par=1,a_par=1,T_par=10,X_0_par=1,qsi=x)
)

hist(as.numeric(X_Optimal_simulated_xi),main='',xlab=expression(paste("X","*"[T])),ylab='# Realisations',col="azure3",xlim=c(0,2),ylim=c(0,25000),breaks = 100)
abline(v =x_concavification, col = "black", lwd = 2,lty=2)

100000-sum(X_Optimal_simulated_xi==0)




#Simulations Xi avec frais#
alpha_tilde<-alpha_sim-c_f_sim-c_s_sim
r_tilde<-r_sim-c_f_sim
theta_tilde_sim<-(alpha_tilde-r_tilde)/sigma_sim

xi_T_simulated_Fee<-lapply(rand_vector3,function (x)exp(-(r_tilde+0.5*theta_tilde_sim^2)*T_sim-theta_tilde_sim*x))
X_Optimal_simulated_xi_Fee<-lapply(xi_T_simulated_Fee, function(x) Ptf_optimalCall_wrt_qsi2(alpha_par=alpha_tilde,r_par=r_tilde,sigma_par=sigma_sim,gamma_par=4,K_par=1,b_par=1,a_par=1,T_par=T_sim,X_0_par=1,qsi=x)
)

par(new = TRUE)
hist(as.numeric(X_Optimal_simulated_xi_Fee),main='',xlab=expression(paste("X","*"[T])),ylab='',yaxt = "n",col=rgb(0.9,0.3,1,alpha=0.2),ylim=c(0,25000),xlim=c(0,2),breaks = 100)
abline(v =x_concavification, col = "black", lwd = 2,lty=2)
legend( 0.2,19000, legend=c( expression(paste("X","*"[T]," sans frais")),expression(paste("X","*"[T]," avec frais"),'Point de concavification')),
        col=c("grey",rgb(0.9,0.3,1,alpha=0.2),"black"), lty=c(1,1,2),lwd = c(4,4,2), cex=0.8)

N_simulations-sum(X_Optimal_simulated_xi_Fee==0) 





######### Simulations Ptf Optimal, Cas de base ##########
rand_vector_base<-sqrt(T_sim)*rnorm(N_simulations)
S_T_simulated_base<-lapply(rand_vector_base,function (x)1*exp(alpha_sim*T_sim-0.5*T_sim*sigma_sim^2+sigma_sim*x))

#Sans frais#
X_Optimal_base_simulated<-lapply(S_T_simulated_base, function(x) Ptf_optimal(S_t_par=x,S_0_par=1,budget_par=1,alpha_par=alpha_sim,r_par=r_sim,sigma_par=sigma_sim,T_par=T_sim,gamma_par=gamma_sim)
)

#Avec frais#
X_Optimal_basefees_simulated<-lapply(S_T_simulated_base, function(x) Ptf_optimal_fee(S_t_par=x,S_0_par=1,budget_par=1,alpha_par=alpha_sim,r_par=r_sim,sigma_par=sigma_sim,T_par=T_sim,gamma_par=gamma_sim,c_f_par=c_f_sim,c_s_par=c_s_sim)
)


### Histogramme combiné ###
#qplot(as.numeric(X_Optimal_base_simulated),geom="histogram",binwidth = 0.005,fill=I("blue"), 
#      col=I("black"),alpha=0.2,xlab=expression(paste("X","*"[T])),ylab='# Realisations') 

No_fee<-data.frame(Optimal_portfolio=as.numeric(X_Optimal_base_simulated))
Fee<-data.frame(Optimal_portfolio=as.numeric(X_Optimal_basefees_simulated))

No_fee$simul<-'No_fee'
Fee$simul<-'Fee'

simul_combin<-rbind(No_fee,Fee)


ggplot(data=simul_combin,aes(simul_combin$Optimal_portfolio,group=simul,fill=simul))+
  geom_histogram(colour='black',binwidth = 0.009,alpha=0.5,position = "identity")+
  scale_fill_manual(name="",values=c("red","darkgray"),labels=c( expression(paste("X","*"[T]," sans frais")),expression(paste("X","*"[T]," avec frais (",c[f],"=0.005"," et ",c[s],"=0.005)"))))+
  theme(legend.position="None")+
  labs(x=expression(paste("X","*"[T])), y="# Réalisations")+
  theme_classic()








############## Demande Anne, 12 décembre 2019 ##############
########### Espérance Utilité, Call, avec frais ##############

#1) Simulation xi_T avec frais;
#2) Calcul du portfeuille optimal avec Call;
#3) Calcul de l'espérance de l'utilité: E[  U( max(X_T,1) ) ]
#4) Faire ça pour c_s=0,0.25, 0.5, 0.75,1,1.25, 1.5, 1.75, 2 (%)
# c_f= 0.02448-c_s


  #1.1)
  rand_vector4<-sqrt(T_sim)*rnorm(N_simulations)
  theta_frais<-(alpha_sim_frais-r_sim_frais)/sigma_sim
  xi_T_simulated_4<-lapply(rand_vector4,function (x)exp(-(r_sim_frais+0.5*theta_frais^2)*T_sim-theta_frais*x))
  
  #2)Ptf_optimalCall_wrt_qsi
  Ptf_call_sim<-as.numeric(lapply(xi_T_simulated_4,function(x)Ptf_optimalCall_wrt_qsi3(alpha_par=alpha_sim_frais,r_par=r_sim_frais,sigma_par=sigma_sim,gamma_par=gamma_sim,K_par=1,b_par=1,a_par=1,T_par=T_sim,qsi=x,lambda_par=lbd_sim_call_f,x_conca_par=x_conca)))
  
  #3) Calcul de l'espérance de l'utilité: E[  U( max(X_T,1) ) ]
  EU_Call_sim<-mean(1/(1-gamma_sim)*(as.numeric(lapply(Ptf_call_sim,function(x)max(x,1))))^(1-gamma_sim))

 result_g7<-data.frame(results=c(-0.05218744,-0.05192372,-0.05082127,-0.04896291))
 result_g6<-data.frame(results=c(-0.06958732, -0.06941325,-0.06857159,-0.06692696))
 result_g5<-data.frame(results=c(-0.09835473,-0.09831733,-0.09749544,-0.09620082))
 result_g4<-data.frame(results=c(-0.1514656,-0.1522911,-0.1523024,-0.1513727))
 result_g3<-data.frame(results=c(-0.2722601,-0.2745076,-0.2756578,-0.2762874))
 result_g2<-data.frame(results=c(-0.6920812,-0.6977341,-0.702836, -0.7063659))
 
 result_g7$gamma<-'g7'   
 result_g6$gamma<-'g6'
 result_g5$gamma<-'g5' 
 result_g4$gamma<-'g4'
 result_g3$gamma<-'g3'
 result_g2$gamma<-'g2'
 
 Compa_EU_gamma<-rbind( result_g7, result_g6, result_g5, result_g4, result_g3, result_g2)
 Compa_EU_gamma$gamma <- factor( Compa_EU_gamma$gamma , levels = c("g2", "g3","g4","g5", "g6","g7"))
 Compa_EU_gamma$c_s<-c(0,0.5,1,1.5)
 
 ggplot(data= Compa_EU_gamma,aes(x=Compa_EU_gamma$c_s,y=Compa_EU_gamma$results,group=gamma,color=gamma))+
    geom_line(aes(linetype=gamma))+
   geom_point(aes(shape=gamma))+
   scale_shape_discrete(name="",
                        breaks=c("g2", "g3", "g4","g5","g6","g7"),
                        labels=c(expression(paste(gamma,"=2")), expression(paste(gamma,"=3")), expression(paste(gamma,"=4")),expression(paste(gamma,"=5")),expression(paste(gamma,"=6")),expression(paste(gamma,"=7"))))+
   scale_color_discrete(name="",
                        breaks=c("g2", "g3", "g4","g5","g6","g7"),
                        labels=c(expression(paste(gamma,"=2")), expression(paste(gamma,"=3")), expression(paste(gamma,"=4")),expression(paste(gamma,"=5")),expression(paste(gamma,"=6")),expression(paste(gamma,"=7"))))+
   scale_linetype_discrete(name="",
                        breaks=c("g2", "g3", "g4","g5","g6","g7"),
                        labels=c(expression(paste(gamma,"=2")), expression(paste(gamma,"=3")), expression(paste(gamma,"=4")),expression(paste(gamma,"=5")),expression(paste(gamma,"=6")),expression(paste(gamma,"=7"))))+
   theme(legend.position = "bottom")+
   #theme_classic()+
   labs(x=expression(paste(c[s],"  (%)")), y=expression(paste("E"^P,"[  U( (",F[T]^"*","-1)"^"+","+1)"," ]")))
   #theme(plot.title = element_text(family = "Helvetica", face = "bold", size = (15),hjust = 0.5))

 

 
 
 
 
 