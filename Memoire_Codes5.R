############ MÉMOIRE: Maitrise en Mathématiques Finacières et Actuarielle  #################

# -SECTION 5: Power Utility - Calculs reliés à la section "Résultats du mémoire" #
#-Note: toutes les formules sont UNIQUEMENT applicables à la fonction de Power Utility

# Notes: ce ficher comporte 2 sections: 
#        -la première: -RÉSULTATS- calculs des résultats pour le mémoire
#        - la deuxième: -FRAIS ÉQUITABLE- contiens les fonctions pour trouver le frais équitable.


############################## PACKAGES ################################################
library(ggplot2)
library(devtools) # installer si ce n'est pas déjà fait.
devtools::install_github("belangerjulie15/Optimisation.Power.Utility")
library(Optimisation.Power.Utility)
########################################################################################




################################# -RÉSULTATS- ##########################################
########################################################################################
#Montant équivalent sûr de l'espérance du bénéfice résultant du fonds distinct
round(Inverse_P_Utility(-0.272,3),3)
round(Inverse_P_Utility(-0.320,3),3)
round(Inverse_P_Utility(-0.358,3),3)
round(Inverse_P_Utility(-0.375,3),3)
round(Inverse_P_Utility(-0.399,3),3)
round(Inverse_P_Utility(-0.423,3),3)
round(Inverse_P_Utility(-0.462,3),3)
round(Inverse_P_Utility(-0.470,3),3)

round(Inverse_P_Utility(-0.276,3),3)
round(Inverse_P_Utility(-0.319,3),3)
round(Inverse_P_Utility(-0.348,3),3)
round(Inverse_P_Utility(-0.375,3),3)
round(Inverse_P_Utility(-0.383,3),3)
round(Inverse_P_Utility(-0.395,3),3)
round(Inverse_P_Utility(-0.414,3),3)
round(Inverse_P_Utility(-0.426,3),3)

round(Inverse_P_Utility(-0.277,3),3)
round(Inverse_P_Utility(-0.315,3),3)
round(Inverse_P_Utility(-0.340,3),3)
round(Inverse_P_Utility(-0.373,3),3)
round(Inverse_P_Utility(-0.375,3),3)
round(Inverse_P_Utility(-0.380,3),3)
round(Inverse_P_Utility(-0.386,3),3)
round(Inverse_P_Utility(-0.381,3),3)





##################################  -FRAIS ÉQUITABLE-  ############################################  
###################################################################################################

########## PARAMÈTRES ##########################################

Maturi<-15         #Time until maturity
r_no_risk<-0.02    #Risk free rate
alpha<-0.04        #Risky rate
sigma<-0.2         #Volatility
gamma<-3           #Parameter of Utility function
S_0<-1             #Initial value of the asset (S_0>0)
B_0<-1             #Initial value of the bank account
budget<-1          #Initial Budget amount
N_Simulations<-100000 #Number of Simulations

Frequ<-52          #Frequency of rebalancing the portfolio

a_call_sim<-1      #Multiplicator of the variable annuity
b_call_sim<-1      #Strike price of the variable annuity
K_call_sim<-1      #Constant of the variable annuity

################################################################


######### Directions ############################################

#(1) Il faut rouler toute la section PARAMÈTRES

#(2)Il faut rouler toute la section FONCTIONS pour activer les fonctions 
#qu'on veut utiliser OU Installer* le package Optimisation.Power.Utility

#(3) Il y a plusieurs sections pour différentes proportions (cte ou dynamiques)
#investies dans l'actif risqué.

#*Pour installer le package, il faut ouvrir le package library(devtools) et rouler:
#build()  install() 
################################################################


####### Installer le package #####################################
library(Optimisation.Power.Utility)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
##################################################################


####### FONCTONS (à rouler)######################################

#1) Payoff d'un call
payoff_call<-function(X_T,b,a,K){
  payoff<-a*max(X_T-b,0)+K
  return(payoff)
}

#2) Fonction de Power Utility
P_Utility<-function(X_Tu,gamma){
  res<-(X_Tu)^(1-gamma)/(1-gamma)
  return(res)
}

#3) Trouve x concavifié
Find_x_theta_PU<-function(b_o,a_o,K_o,gamma_o){
  f <- function(x) payoff_call(X_T=x,b=b_o,a=a_o,K=K_o)^(-gamma_o)*x*a_o-P_Utility(payoff_call(x,b=b_o,a=a_o,K=K_o),gamma_o)+P_Utility(K_o,gamma_o)        #-Finding the fee
  result<-uniroot(f, c(b_o,b_o+100))$root
  return(result)
}

#4) Trouver le langrangien optimal (pour call bien sûr)
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

lambda_optimal<-function(r_FLo,alpha_FLo,sigma_FLo,gamma_FLo,K_FLo,b_FLo,a_FLo,X_0_FLo,T_FLo){
  
  f <- function(x) Find_lambda(r_FL=r_FLo,alpha_FL=alpha_FLo,sigma_FL=sigma_FLo,gamma_FL=gamma_FLo,K_FL=K_FLo,b_FL=b_FLo,a_FL=a_FLo,T_FL=T_FLo,lambda_FL = x)-X_0_FLo
  
  lambda_result<-uniroot(f, c(0.0000000000000001,100))$root  
  return(lambda_result)
}

#5) Fonction de tradding à chaque période

# Inputs:
#invest_risq_tm1-> nombres d'unités dans actif risqué avant le rebalancement au temps t
#invest_srisq_tm1-> nombre d'unités dans actif sans riqueavant le rebalancement au temps t
#Outputs:
#nombres d'unités dans actif risqué/sans risque après le rebalancement.
Investment_Fonct<-function(S_t,B_t,invest_risq_tm1,invest_srisq_tm1,propor){
  Ptf_ini<-S_t*invest_risq_tm1+invest_srisq_tm1*B_t
  
  investment<-rbind(propor*Ptf_ini/S_t,Ptf_ini*(1-propor)/B_t)
  
  return(investment)
}
################################################################



########### SECTION 1-) "Proportion optimale", avec la m?thode martingale ###########

frais_eq_martingale<-function(para_c_s,para_c_f){
  library(Optimisation.Power.Utility)
  #initialisation
  alpha_tilde<-alpha-para_c_s-para_c_f
  r_no_risk_tilde<-r_no_risk-para_c_f
  theta_sim<-(alpha-r_no_risk)/sigma
  theta_sim_tilde<-(alpha_tilde-r_no_risk_tilde)/sigma
  
  cte_Merton<-theta_sim_tilde/(sigma*gamma)
  
  
  lambda_opt<-lambda_optimal(r_FL=r_no_risk_tilde,alpha_FL=alpha_tilde,sigma_FL=sigma,gamma_FL=gamma,K_FL=K_call_sim,b_FL=b_call_sim,a_FL=a_call_sim,X_0_FL=budget,T_FL=Maturi)
  x_concavi<-Find_x_theta_PU(b_call_sim,a_call_sim,K_call_sim,gamma)
  y_concavi<-x_concavi^(-gamma)
  
  random_var<-matrix(rnorm(Maturi*Frequ*N_Simulations),N_Simulations,Maturi*Frequ) #Matrice actif sans risque
  pre1_S_tilde_t<-S_0*exp((alpha_tilde-0.5*sigma^2)/Frequ+sigma*sqrt(1/Frequ)*random_var)
  pre2_S_tilde_t<-cbind((rep(S_0,N_Simulations)),pre1_S_tilde_t)
  
  pre1_xi_tilde_t<-exp(-(r_no_risk_tilde+0.5*theta_sim_tilde^2)/Frequ-theta_sim_tilde*sqrt(1/Frequ)*random_var)
  pre2_xi_tilde_t<-cbind((rep(S_0,N_Simulations)),pre1_xi_tilde_t)
  
  pre1_xi_t<-exp(-(r_no_risk+0.5*theta_sim^2)/Frequ-theta_sim*sqrt(1/Frequ)*random_var)
  pre2_xi_t<-cbind((rep(S_0,N_Simulations)),pre1_xi_t)
  
  
  pre1_B_tilde_t<-rep(exp(r_no_risk_tilde/Frequ),Frequ*Maturi) #vecteur actif sans risque
  B_tilde_t<-B_0*cumprod(c(B_0,pre1_B_tilde_t))
  
  
  matrice_pre2_S<-pre2_S_tilde_t
  matrice_pre2_xi_tilde_t<-pre2_xi_tilde_t
  matrice_pre2_xi<-pre2_xi_t
  vecteur_B<-B_tilde_t
    
  matrice_S<-matrix(0,N_Simulations,Frequ*Maturi+1)
  nu_Act_R<-matrix(0,N_Simulations,Frequ*Maturi+1)
  interior<-matrix(0,N_Simulations,Frequ*Maturi+1)
  processus_ptf<- matrix(0,N_Simulations,Frequ*Maturi+1)
  unit_invest<-c(0,0)
  unit_risque<-matrix(0,N_Simulations,Frequ*Maturi+1)
  unit_n_risque<-matrix(0,N_Simulations,Frequ*Maturi+1)
  funds_d<-rep(0,N_Simulations)
  Uprocessus_ptf<-matrix(0,N_Simulations,Frequ*Maturi+1)
  matrice_xi_tilde<-matrix(0,N_Simulations,Frequ*Maturi+1)
  matrice_xi<-matrix(0,N_Simulations,Frequ*Maturi+1)
  
  for (n in 1:(N_Simulations)){ #Loop sur les simulations
    matrice_S[n,]<-cumprod(matrice_pre2_S[n,])
    matrice_xi_tilde[n,]<-cumprod(matrice_pre2_xi_tilde_t[n,])
    matrice_xi[n,]<-cumprod(matrice_pre2_xi[n,])
    
    for(m in 1:(Frequ*Maturi)){ #Frequ*Maturi
      petit_t<-(m-1)/(Frequ)
      
      interior[n,m]<-((r_no_risk_tilde+(1/gamma-0.5)*theta_sim_tilde^2)*(Maturi-petit_t)+log(y_concavi/(matrice_xi_tilde[n,m]*lambda_opt)))/(theta_sim_tilde*sqrt(Maturi-petit_t))
      if(interior[n,m]<(-30)){interior[n,m]<--30}
      nu_Act_R[n,m]<-max(min((theta_sim_tilde/gamma+dnorm(interior[n,m])/(pnorm(interior[n,m])*sqrt(Maturi-petit_t)))/sigma,2),0)
      
      if(m==1){
        unit_risque[n,m]<-nu_Act_R[n,m]*budget/matrice_S[n,m]
        unit_n_risque[n,m]<-(budget-unit_risque[n,m]*matrice_S[n,m])/B_tilde_t[m]
        
        processus_ptf[n,m]<-budget
      }
      else{
        unit_invest<-Investment_Fonct(matrice_S[n,m],B_tilde_t[m],unit_risque[n,(m-1)], unit_n_risque[n,(m-1)],nu_Act_R[n,m])
        unit_risque[n,m]<-unit_invest[1]
        unit_n_risque[n,m]<-unit_invest[2]
        processus_ptf[n,m]<-unit_risque[n,m]*matrice_S[n,m]+unit_n_risque[n,m]*B_tilde_t[m]
      }
      
    }
    processus_ptf[n,(Frequ*Maturi+1)]<-unit_risque[n,(Frequ*Maturi)]*matrice_S[n,(Frequ*Maturi+1)]+unit_n_risque[n,(Frequ*Maturi)]*B_tilde_t[(Frequ*Maturi+1)]
    funds_d[n]<-a_call_sim*max(0,processus_ptf[n,(Frequ*Maturi+1)]-b_call_sim)+K_call_sim
  }
  
  #CB<-mean(matrice_xi_tilde[,(Frequ*Maturi+1)]*processus_ptf[,(Frequ*Maturi+1)])
  #exercice_guarantie<-sum(processus_ptf[,(Frequ*Maturi+1)]<b_call_sim)/N_Simulations  
  #call_tout_t<-apply(processus_ptf,c(1,2),function(x) a_call_sim*max(0,x-b_call_sim)+K_call_sim)
  #Uty_t<-apply(processus_ptf,c(1,2),function(x) P_Utility(X_Tu=x,gamma=gamma))#utilit? des valeurs du ptf pour tout t
  #Uty_f_t<-apply(call_tout_t,c(1,2),function(x) P_Utility(X_Tu=x,gamma=gamma))#utilit? des valeurs du fonds pour tout t
  
  #Uty<-as.numeric(lapply(funds_d,function(x)P_Utility(X_Tu=x,gamma=gamma)))
  #EU<-mean(Uty)
  #Uprocessus_ptf<-P_Utility(X_Tu=processus_ptf,gamma=gamma)
  #verif<-mean(matrice_xi_tilde[,(Frequ*Maturi+1)]*matrice_S[,(Frequ*Maturi+1)])
  
  budg_frais_equ<-mean(matrice_xi[,(Frequ*Maturi+1)]*funds_d)
    
  return(budg_frais_equ)#processus_ptf[,(Frequ*Maturi+1)]c(EU,CB)colMeans(Uprocessus_ptf)c(verif,CB)exercice_guarantiecolMeans(Uty_t[,])
}

#vérification#
timer<-proc.time()
frais_eq_martingale(para_c_s=0.018,para_c_f=0.02448-0.018)#prend 361.96 sec (pour 1 fois)
proc.time()-timer


#Optimisation#
timer1<-proc.time()
g <- function(x) frais_eq_martingale(para_c_s=0.01224,para_c_f=x)- budget      #-Finding the fee
result<-uniroot(g, c(-0.01,0.01))$root
proc.time()-timer1
#0.04722044
#0.04731921 #Kornos:7297.97 secondes

#Essai parrallélisation#
registerDoParallel(cores=3)
system.time(result_MM<-as.numeric(foreach(i=c(0.0,0.01224,0.018)) %dopar% uniroot(function(x)frais_eq_martingale(para_c_s=i,para_c_f=x)- budget    
                                                                                 , c(0.000001,0.05),tol= 0.000001)$root))

#Kronos: (T=15, cores=3) 17767.19 sec

########### SECTION 2-) "Proportion optimale", avec la m?thode martingale, proportion born?e [0,1] ###########

frais_eq_martingale_Borne<-function(para_c_s,para_c_f){
  library(Optimisation.Power.Utility)
  
  #initialisation
  alpha_tilde<-alpha-para_c_s-para_c_f
  r_no_risk_tilde<-r_no_risk-para_c_f
  theta_sim<-(alpha-r_no_risk)/sigma
  theta_sim_tilde<-(alpha_tilde-r_no_risk_tilde)/sigma
  
  cte_Merton<-theta_sim_tilde/(sigma*gamma)
  
  
  lambda_opt<-lambda_optimal(r_FL=r_no_risk_tilde,alpha_FL=alpha_tilde,sigma_FL=sigma,gamma_FL=gamma,K_FL=K_call_sim,b_FL=b_call_sim,a_FL=a_call_sim,X_0_FL=budget,T_FL=Maturi)
  x_concavi<-Find_x_theta_PU(b_call_sim,a_call_sim,K_call_sim,gamma)
  y_concavi<-x_concavi^(-gamma)
  
  random_var<-matrix(rnorm(Maturi*Frequ*N_Simulations),N_Simulations,Maturi*Frequ) #Matrice actif sans risque
  pre1_S_tilde_t<-S_0*exp((alpha_tilde-0.5*sigma^2)/Frequ+sigma*sqrt(1/Frequ)*random_var)
  pre2_S_tilde_t<-cbind((rep(S_0,N_Simulations)),pre1_S_tilde_t)
  
  pre1_xi_tilde_t<-exp(-(r_no_risk_tilde+0.5*theta_sim_tilde^2)/Frequ-theta_sim_tilde*sqrt(1/Frequ)*random_var)
  pre2_xi_tilde_t<-cbind((rep(S_0,N_Simulations)),pre1_xi_tilde_t)
  
  pre1_xi_t<-exp(-(r_no_risk+0.5*theta_sim^2)/Frequ-theta_sim*sqrt(1/Frequ)*random_var)
  pre2_xi_t<-cbind((rep(S_0,N_Simulations)),pre1_xi_t)
  
  
  pre1_B_tilde_t<-rep(exp(r_no_risk_tilde/Frequ),Frequ*Maturi) #vecteur actif sans risque
  B_tilde_t<-B_0*cumprod(c(B_0,pre1_B_tilde_t))
  
  
  matrice_pre2_S<-pre2_S_tilde_t
  matrice_pre2_xi_tilde_t<-pre2_xi_tilde_t
  matrice_pre2_xi<-pre2_xi_t
  vecteur_B<-B_tilde_t
  
  matrice_S<-matrix(0,N_Simulations,Frequ*Maturi+1)
  nu_Act_R<-matrix(0,N_Simulations,Frequ*Maturi+1)
  interior<-matrix(0,N_Simulations,Frequ*Maturi+1)
  processus_ptf<- matrix(0,N_Simulations,Frequ*Maturi+1)
  unit_invest<-c(0,0)
  unit_risque<-matrix(0,N_Simulations,Frequ*Maturi+1)
  unit_n_risque<-matrix(0,N_Simulations,Frequ*Maturi+1)
  funds_d<-rep(0,N_Simulations)
  Uprocessus_ptf<-matrix(0,N_Simulations,Frequ*Maturi+1)
  matrice_xi_tilde<-matrix(0,N_Simulations,Frequ*Maturi+1)
  matrice_xi<-matrix(0,N_Simulations,Frequ*Maturi+1)
  
  for (n in 1:(N_Simulations)){ #Loop sur les simulations
    matrice_S[n,]<-cumprod(matrice_pre2_S[n,])
    matrice_xi_tilde[n,]<-cumprod(matrice_pre2_xi_tilde_t[n,])
    matrice_xi[n,]<-cumprod(matrice_pre2_xi[n,])
    
    for(m in 1:(Frequ*Maturi)){ #Frequ*Maturi
      petit_t<-(m-1)/(Frequ)
      
      interior[n,m]<-((r_no_risk_tilde+(1/gamma-0.5)*theta_sim_tilde^2)*(Maturi-petit_t)+log(y_concavi/(matrice_xi_tilde[n,m]*lambda_opt)))/(theta_sim_tilde*sqrt(Maturi-petit_t))
      if(interior[n,m]<(-30)){interior[n,m]<--30}
      nu_Act_R[n,m]<-max(min((theta_sim_tilde/gamma+dnorm(interior[n,m])/(pnorm(interior[n,m])*sqrt(Maturi-petit_t)))/sigma,1),0)
      
      if(m==1){
        unit_risque[n,m]<-nu_Act_R[n,m]*budget/matrice_S[n,m]
        unit_n_risque[n,m]<-(budget-unit_risque[n,m]*matrice_S[n,m])/B_tilde_t[m]
        
        processus_ptf[n,m]<-budget
      }
      else{
        unit_invest<-Investment_Fonct(matrice_S[n,m],B_tilde_t[m],unit_risque[n,(m-1)], unit_n_risque[n,(m-1)],nu_Act_R[n,m])
        unit_risque[n,m]<-unit_invest[1]
        unit_n_risque[n,m]<-unit_invest[2]
        processus_ptf[n,m]<-unit_risque[n,m]*matrice_S[n,m]+unit_n_risque[n,m]*B_tilde_t[m]
      }
      
    }
    processus_ptf[n,(Frequ*Maturi+1)]<-unit_risque[n,(Frequ*Maturi)]*matrice_S[n,(Frequ*Maturi+1)]+unit_n_risque[n,(Frequ*Maturi)]*B_tilde_t[(Frequ*Maturi+1)]
    funds_d[n]<-a_call_sim*max(0,processus_ptf[n,(Frequ*Maturi+1)]-b_call_sim)+K_call_sim
  }
  
  #CB<-mean(matrice_xi_tilde[,(Frequ*Maturi+1)]*processus_ptf[,(Frequ*Maturi+1)])
  #exercice_guarantie<-sum(processus_ptf[,(Frequ*Maturi+1)]<b_call_sim)/N_Simulations  
  #call_tout_t<-apply(processus_ptf,c(1,2),function(x) a_call_sim*max(0,x-b_call_sim)+K_call_sim)
  #Uty_t<-apply(processus_ptf,c(1,2),function(x) P_Utility(X_Tu=x,gamma=gamma))#utilit? des valeurs du ptf pour tout t
  #Uty_f_t<-apply(call_tout_t,c(1,2),function(x) P_Utility(X_Tu=x,gamma=gamma))#utilit? des valeurs du fonds pour tout t
  
  #Uty<-as.numeric(lapply(funds_d,function(x)P_Utility(X_Tu=x,gamma=gamma)))
  #EU<-mean(Uty)
  #Uprocessus_ptf<-P_Utility(X_Tu=processus_ptf,gamma=gamma)
  #verif<-mean(matrice_xi_tilde[,(Frequ*Maturi+1)]*matrice_S[,(Frequ*Maturi+1)])
  
  budg_frais_equ<-mean(matrice_xi[,(Frequ*Maturi+1)]*funds_d)
  
  return(budg_frais_equ)#processus_ptf[,(Frequ*Maturi+1)]c(EU,CB)colMeans(Uprocessus_ptf)c(verif,CB)exercice_guarantiecolMeans(Uty_t[,])
}

#vérification#
timer2<-proc.time()
frais_eq_martingale_Borne(para_c_s=0.0,para_c_f=0.02448)#prend 361.96 sec (pour 1 fois)
proc.time()-timer2

#Optimisation#
timer3<-proc.time()
h <- function(x) frais_eq_martingale_Borne(para_c_s=0.01224,para_c_f=x)- budget      #-Finding the fee
result<-uniroot(h, c(-0.01,0.005))$root
proc.time()-timer3

#Essai parallélisation#
registerDoParallel(cores=3)
system.time(result_o<-as.numeric(foreach(i=c(0.0,0.01224,0.018)) %dopar% uniroot(function(x)frais_eq_martingale_Borne(para_c_s=i,para_c_f=x)- budget   
                                                                                   , c(-0.01,0.05),tol= 0.000001)$root))


        #Kronos: cores=3, T=15,   17980.03 sec
########### SECTION 3-) "Proportion" dans l'actif risqué est constante (repré. par un paramètre) ###########

frais_eq_prop_cte<-function(para_c_s,para_c_f,prop_act_r){
  library(Optimisation.Power.Utility)
  
  #initialisation
  alpha_tilde<-alpha-para_c_s-para_c_f
  r_no_risk_tilde<-r_no_risk-para_c_f
  theta_sim<-(alpha-r_no_risk)/sigma
  theta_sim_tilde<-(alpha_tilde-r_no_risk_tilde)/sigma
  
  cte_Merton<-theta_sim_tilde/(sigma*gamma)
  
  lambda_opt<-lambda_optimal(r_FL=r_no_risk_tilde,alpha_FL=alpha_tilde,sigma_FL=sigma,gamma_FL=gamma,K_FL=K_call_sim,b_FL=b_call_sim,a_FL=a_call_sim,X_0_FL=budget,T_FL=Maturi)
  x_concavi<-Find_x_theta_PU(b_call_sim,a_call_sim,K_call_sim,gamma)
  y_concavi<-x_concavi^(-gamma)
  
  random_var<-matrix(rnorm(Maturi*Frequ*N_Simulations),N_Simulations,Maturi*Frequ) #Matrice actif sans risque
  pre1_S_tilde_t<-S_0*exp((alpha_tilde-0.5*sigma^2)/Frequ+sigma*sqrt(1/Frequ)*random_var)
  pre2_S_tilde_t<-cbind((rep(S_0,N_Simulations)),pre1_S_tilde_t)
  
  pre1_xi_tilde_t<-exp(-(r_no_risk_tilde+0.5*theta_sim_tilde^2)/Frequ-theta_sim_tilde*sqrt(1/Frequ)*random_var)
  pre2_xi_tilde_t<-cbind((rep(S_0,N_Simulations)),pre1_xi_tilde_t)
  
  pre1_xi_t<-exp(-(r_no_risk+0.5*theta_sim^2)/Frequ-theta_sim*sqrt(1/Frequ)*random_var)
  pre2_xi_t<-cbind((rep(S_0,N_Simulations)),pre1_xi_t)
  
  pre1_B_tilde_t<-rep(exp(r_no_risk_tilde/Frequ),Frequ*Maturi) #vecteur actif sans risque
  B_tilde_t<-B_0*cumprod(c(B_0,pre1_B_tilde_t))
  
  
  matrice_pre2_S<-pre2_S_tilde_t
  matrice_pre2_xi_tilde_t<-pre2_xi_tilde_t
  matrice_pre2_xi<-pre2_xi_t
  vecteur_B<-B_tilde_t
  
  #initialisation
  
  
  matrice_S<-matrix(0,N_Simulations,Frequ*Maturi+1)
  nu_Act_R<-matrix(0,N_Simulations,Frequ*Maturi+1)
  xi_tilde<-matrix(0,N_Simulations,Frequ*Maturi+1)
  interior<-matrix(0,N_Simulations,Frequ*Maturi+1)
  processus_ptf<- matrix(0,N_Simulations,Frequ*Maturi+1)
  unit_invest<-c(0,0)
  unit_risque<-matrix(0,N_Simulations,Frequ*Maturi+1)
  unit_n_risque<-matrix(0,N_Simulations,Frequ*Maturi+1)
  funds_d<-rep(0,N_Simulations)
  Uprocessus_ptf<-matrix(0,N_Simulations,Frequ*Maturi+1)
  matrice_xi_tilde<-matrix(0,N_Simulations,Frequ*Maturi+1)
  matrice_xi<-matrix(0,N_Simulations,Frequ*Maturi+1)
  
  for (n in 1:(N_Simulations)){ #Loop sur les simulations
    matrice_S[n,]<-cumprod(matrice_pre2_S[n,])
    matrice_xi_tilde[n,]<-cumprod(matrice_pre2_xi_tilde_t[n,])
    matrice_xi[n,]<-cumprod(matrice_pre2_xi[n,])
    
    for(m in 1:(Frequ*Maturi)){ #Frequ*Maturi
      petit_t<-(m-1)/(Frequ)
      
      nu_Act_R[n,m]<-prop_act_r
      
      if(m==1){
        unit_risque[n,m]<-nu_Act_R[n,m]*budget/matrice_S[n,m]
        unit_n_risque[n,m]<-(budget-unit_risque[n,m]*matrice_S[n,m])/B_tilde_t[m]
        
        processus_ptf[n,m]<-budget
      }
      else{
        unit_invest<-Investment_Fonct(matrice_S[n,m],B_tilde_t[m],unit_risque[n,(m-1)], unit_n_risque[n,(m-1)],nu_Act_R[n,m])
        unit_risque[n,m]<-unit_invest[1]
        unit_n_risque[n,m]<-unit_invest[2]
        processus_ptf[n,m]<-unit_risque[n,m]*matrice_S[n,m]+unit_n_risque[n,m]*B_tilde_t[m]
      }
      
    }
    processus_ptf[n,(Frequ*Maturi+1)]<-unit_risque[n,(Frequ*Maturi)]*matrice_S[n,(Frequ*Maturi+1)]+unit_n_risque[n,(Frequ*Maturi)]*B_tilde_t[(Frequ*Maturi+1)]
    funds_d[n]<-a_call_sim*max(0,processus_ptf[n,(Frequ*Maturi+1)]-b_call_sim)+K_call_sim
    
  }
  
  #exercice_guarantie<-sum(processus_ptf[,(Frequ*Maturi+1)]<b_call_sim)/N_Simulations  
  #call_tout_t<-apply(processus_ptf,c(1,2),function(x) a_call_sim*max(0,x-b_call_sim)+K_call_sim)
  #Uty<-as.numeric(lapply(funds_d,function(x)P_Utility(X_Tu=x,gamma=gamma)))
  #Uty_t<-apply(processus_ptf,c(1,2),function(x) P_Utility(X_Tu=x,gamma=gamma))#utilit? des valeurs du ptf pour tout t
  #Uty_f_t<-apply(call_tout_t,c(1,2),function(x) P_Utility(X_Tu=x,gamma=gamma))#utilit? des valeurs du fonds pour tout t
  
  #EU<-mean(Uty)
  
  #xi_tilde[n,(Frequ*Maturi+1)]<-matrice_S[n,(Frequ*Maturi+1)]^(-theta_sim_tilde/sigma)*exp((theta_sim_tilde*alpha_tilde/sigma-theta_sim_tilde*sigma/2-r_no_risk_tilde-0.5*theta_sim_tilde^2)*petit_t)
  
  #contrainte_budget<-mean(matrice_xi_tilde[,(Frequ*Maturi+1)]*processus_ptf[,(Frequ*Maturi+1)])
  #Uprocessus_ptf<-P_Utility(X_Tu=processus_ptf,gamma=gamma)
  #verif<-mean(matrice_xi_tilde[,(Frequ*Maturi+1)]*matrice_S[,(Frequ*Maturi+1)])
  
  budg_frais_equ<-mean(matrice_xi[,(Frequ*Maturi+1)]*funds_d)
  
  return(budg_frais_equ)#c(EU,contrainte_budget)colMeans(Uprocessus_ptf)c(verif,contrainte_budget)exercice_guarantieprocessus_ptf[,(Frequ*Maturi+1)]
}

timer3<-proc.time()
frais_eq_prop_cte(para_c_s=0.00,para_c_f=0.00,prop_act_r=0.2)#prend  168.76  sec (pour 1 fois)  
proc.time()-timer3


timer4<-proc.time()

f <- function(x) frais_eq_prop_cte(para_c_s=0.01224,para_c_f=x,prop_act_r=1)- budget      #-Finding the fee
result<-uniroot(f, c(-0.01,0.02),tol= 0.000001)$root

proc.time()-timer4# Kronos: c_s=1.224 ->3018.17  sec (pas de tol)
                  # Kronos: c_s=0.0%  -> 5796.11 sec (tol=1e-9)
                  # Kronos: cs=0.0% -> 11755.83 
                  # Kronos: cs=0.0% T=20 ->20 10151.64 sec
#0.02438858
#0.02438242
#0.02461108

##### Essai paralléllisation ####
registerDoParallel(cores=2)
system.time(result_p<-as.numeric(foreach(i=c(0.01224,0.018)) %dopar% uniroot(function(x) frais_eq_prop_cte(para_c_s=i,para_c_f=x,prop_act_r=1)- budget      #-Finding the fee
                                                                       , c(-0.01,0.005),tol= 0.000001)$root))
#Kronos: cte=1 (cores=2), T=15 9633.14 sec
#Kronos:  T=15 10102.58 
#Kronos: T=15 (cores=3) 8161.69   
#Kronos: T=15 (cores=2) 8869.61  sec



########### SECTION 4-) Simulations du portefeuille optimal directement ? maturit? ###########

frais_eq_fonds_distinct_maturite<-function(para_c_s,para_c_f){
  library(Optimisation.Power.Utility)
  #initialisation des  param?tres
  alpha_tilde<-alpha-para_c_s-para_c_f
  r_no_risk_tilde<-r_no_risk-para_c_f
  theta_sim<-(alpha-r_no_risk)/sigma
  theta_sim_tilde<-(alpha_tilde-r_no_risk_tilde)/sigma
  
  lambda_opt<-lambda_optimal(r_FL=r_no_risk_tilde,alpha_FL=alpha_tilde,sigma_FL=sigma,gamma_FL=gamma,K_FL=K_call_sim,b_FL=b_call_sim,a_FL=a_call_sim,X_0_FL=budget,T_FL=Maturi)
  x_concavi<-Find_x_theta_PU(b_call_sim,a_call_sim,K_call_sim,gamma)
  y_concavi<-x_concavi^(-gamma)
  
  random_var<-matrix(rnorm(Maturi*Frequ*N_Simulations),N_Simulations,Maturi*Frequ) #Matrice actif sans risque
  pre1_S_tilde_t<-S_0*exp((alpha_tilde-0.5*sigma^2)/Frequ+sigma*sqrt(1/Frequ)*random_var)
  pre2_S_tilde_t<-cbind((rep(S_0,N_Simulations)),pre1_S_tilde_t)
  matrice_pre2_S<-pre2_S_tilde_t
  
  matrice_S<-matrix(0,N_Simulations,Frequ*Maturi+1)
  xi_tilde<-matrix(0,N_Simulations,Frequ*Maturi+1)
  matrice_xi<-matrix(0,N_Simulations,Frequ*Maturi+1)
  funds_d<-rep(0,N_Simulations)
  
  pre1_xi_t<-exp(-(r_no_risk+0.5*theta_sim^2)/Frequ-theta_sim*sqrt(1/Frequ)*random_var)
  pre2_xi_t<-cbind((rep(S_0,N_Simulations)),pre1_xi_t)
  matrice_pre2_xi<-pre2_xi_t
  
  for (n in 1:(N_Simulations)){ #Loop sur les simulations
    matrice_S[n,]<-cumprod(matrice_pre2_S[n,])
    matrice_xi[n,]<-cumprod(matrice_pre2_xi[n,])
    
    for(m in 1:(Frequ*Maturi+1)){ #Frequ*Maturi
      petit_t<-(m-1)/(Frequ) 
      xi_tilde[n,m]<-matrice_S[n,m]^(-theta_sim_tilde/sigma)*exp((theta_sim_tilde*alpha_tilde/sigma-theta_sim_tilde*sigma/2-r_no_risk_tilde-0.5*theta_sim_tilde^2)*petit_t)
    }
  }
  
  pre_ptf_optimal<-as.numeric((lambda_opt* xi_tilde[,(Frequ*Maturi+1)])^(-1/gamma))
  position_ptf_optimal_rej<-pre_ptf_optimal<x_concavi
  pre_ptf_optimal[position_ptf_optimal_rej]<-0
  ptf_optimal<-pre_ptf_optimal
  
  for (n in 1:(N_Simulations)){
    funds_d[n]<-a_call_sim*max(0,ptf_optimal[n]-b_call_sim)+K_call_sim
  }
  #exercice_guarantie<-sum(ptf_optimal[]<b_call_sim)/N_Simulations  
  
  #EU<-mean(as.numeric(lapply(funds_d,function(x)P_Utility(X_Tu=x,gamma=gamma))))
  #CB<-mean(ptf_optimal*xi_tilde[,(Frequ*Maturi+1)])
  #verif<-mean(matrice_S[,(Frequ*Maturi+1)]*xi_tilde[,(Frequ*Maturi+1)])
  budg_frais_equ<-mean(matrice_xi[,(Frequ*Maturi+1)]*funds_d)
  
  return(budg_frais_equ)#c(CB,EU) verifc(CB,EU,verif,exercice_guarantie)
}
#Vérification#
timer5<-proc.time()
 frais_eq_fonds_distinct_maturite(para_c_s=0.018,para_c_f=0.00648)#funds_final
proc.time()-timer5

#Optimisation#
timer6<-proc.time()
f <- function(x) frais_eq_fonds_distinct_maturite(para_c_s=0.0,para_c_f=x)- budget      #-Finding the fee
result<-uniroot(f, c(0.00001,0.04),tol= 0.000001)$root
proc.time()-timer6 #Kornos: c_S=1.8%: 678.66 sec

#Essai parrallélisation#
registerDoParallel(cores=3)
system.time(result_opt<-as.numeric(foreach(i=c(0.0,0.01224,0.018)) %dopar% uniroot(function(x) frais_eq_fonds_distinct_maturite(para_c_s=i,para_c_f=x)- budget
                                                                       , c(-0.1,0.1),tol= 0.000001)$root))
# Kronos, T=15: 2600.06 sec

#0.04026580 0.02265016 0.01612437
proportion_funds_optimale<-function(S_tilde,alpha_prop,r_prop,sigma_prop,c_f_prop,c_s_prop,T_prop1,t_prop2,gamma_prop){
  alpha_tilde<-alpha_prop-c_s_prop-c_f_prop
  r_no_risk_tilde<-r_prop-c_f_prop
  theta_sim_tilde<-(alpha_tilde-r_no_risk_tilde)/sigma_prop
  x_concavi<-Find_x_theta_PU(1,1,1,gamma_prop)
  y_concavi<-x_concavi^(-gamma_prop)
  lambda_prop<-lambda_optimal(r_FL=r_no_risk_tilde,alpha_FL=alpha_tilde,sigma_FL=sigma_prop,gamma_FL=gamma_prop,K_FL=1,b_FL=1,a_FL=1,X_0_FL=1,T_FL=T_prop1)
  
  xi_tilde_prop<-S_tilde^(-theta_sim_tilde/sigma_prop)*exp((theta_sim_tilde*alpha_tilde/sigma_prop-theta_sim_tilde*sigma_prop/2-r_no_risk_tilde-0.5*theta_sim_tilde^2)*t_prop2)

  
  j<-((r_no_risk_tilde+(1/gamma_prop-0.5)*theta_sim_tilde^2)*(T_prop1-t_prop2)+log(y_concavi/(xi_tilde_prop*lambda_prop)))/(theta_sim_tilde*sqrt(T_prop1-t_prop2))
  if(j<(-30)){j<--30}
  proport<-max(min((theta_sim_tilde/gamma_prop+dnorm(j)/(pnorm(j)*sqrt(T_prop1-t_prop2)))/sigma_prop,1),0)
  
  cte_M<-theta_sim_tilde/(sigma_prop*gamma_prop)
  print(cte_M)
  return(proport)
}

axe_d_s<-seq(0.9,5,0.01)
p_g3<-lapply(axe_d_s,function(s)proportion_funds_optimale(S_tilde=s,alpha_prop=0.04,r_prop=0.02,sigma_prop=0.2,c_f_prop=0.01224,c_s_prop=0.01224,T_prop1=10,t_prop2=1,gamma_prop=3))


plot(axe_d_s,p_g3,type='l')
lines(axe_d_s,rep(0.06466667,length(axe_d_s)),col='red',lty=2)



### -Frais équitables (c_s+c_f): Sans simulation- ###
##################### Annexe D ######################

frais_eq_prop_cte_theo<-function(r_s_r,si,m_T,bud,garantie,propo){
  iterat<-function(cf){
    
  r_til<-r_s_r-cf
  k<-(garantie-(1-propo)*exp(r_til*m_T))/propo
  
  d1<-((r_til+0.5*si^2)*m_T-log(k))/(si*sqrt(m_T))
  d2<- (log(k)-(r_til-0.5*si^2)*m_T)/(si*sqrt(m_T))
  budget_tempo<-exp(-cf*m_T)*(propo*pnorm(d1)+(1-propo))+(exp(-r_s_r*m_T)*garantie-(1-propo)*exp(-cf*m_T))*pnorm(d2)-bud
  return(budget_tempo)}
  
  cf_equita<-uniroot(iterat, c(0.000001,0.10),tol= 0.000001)$root
  
  return(cf_equita)
}

frais_eq_prop_cte_theo(r_s_r=0.02,si=0.2,m_T=10,bud=1,garantie=1,propo=1)







