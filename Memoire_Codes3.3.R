############ MÉMOIRE: Maitrise en Mathématiques Finacières et Actuarielle  #################

# -SECTION 4: COMPARAISONS: Power Utility - Simulation de marché avec frais #
#-Note: toutes les formules sont UNIQUEMENT applicables à la fonction de Power Utility


######### DIRECTIONS ##########################################################
#-Note: toutes les formules sont UNIQUEMENT applicables à la fonction Power Utility

#(1) Il faut rouler toute la section PARAMÈTRES

#(2) Il faut rouler toute la section FONCTIONS pour activer les fonctions qu'on
#veut utiliser    OU   installer/activer* le package Optimisation.Power.Utility de Github,
#il suffit de rouler la section PACKAGES.

#(3) Il y a plusieurs sections pour diff?rentes proportions (cte ou dynamiques)
#investies dans l'actif risqu?.

#(4) Pour enregistrer les modifications sur Github: Commit > Push > Pull

#* Pour installer le package, il faut ouvrir le package et rouler:
#build() install()
################################################################################




########## PACKAGES ##########################################
library(ggplot2)

library(devtools) # installer si ce n'est pas déjà fait.
devtools::install_github("belangerjulie15/Optimisation.Power.Utility")
library(Optimisation.Power.Utility)
###################################################################

############## PARAMÈTRES ########################################
Maturi<-10         #Time until maturity
r_no_risk<-0.02    #Risk free rate
alpha<-0.04        #Risky rate
sigma<-0.2         #Volatility
gamma<-3           #Parameter of Utility function
S_0<-1             #Initial value of the asset (S_0>0)
B_0<-1             #Initial value of the bank account
budget<-1          #Initial Budget amount
N_Simulations<-100000 #Number of Simulations
fee_c_s<-0.018    #Fee applied of the risky asset
fee_c_f<-0.00648     #Fee applied of the funds 
Frequ<-52          #Frequency of rebalancing the portfolio

a_call_sim<-1      #Multiplicator of the variable annuity
b_call_sim<-1      #Strike price of the variable annuity
K_call_sim<-1      #Constant of the variable annuity

alpha_tilde<-alpha-fee_c_s-fee_c_f
r_no_risk_tilde<-r_no_risk-fee_c_f
theta_sim<-(alpha-r_no_risk)/sigma
theta_sim_tilde<-(alpha_tilde-r_no_risk_tilde)/sigma

cte_Merton<-theta_sim_tilde/(sigma*gamma)
################################################################

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


########### SECTION 0-) SIMULATIONS CONSTANTES POUR COMPARAISON ###########

lambda_opt<-lambda_optimal(r_FL=r_no_risk_tilde,alpha_FL=alpha_tilde,sigma_FL=sigma,gamma_FL=gamma,K_FL=K_call_sim,b_FL=b_call_sim,a_FL=a_call_sim,X_0_FL=budget,T_FL=Maturi)
x_concavi<-Find_x_theta_PU(b_call_sim,a_call_sim,K_call_sim,gamma)
y_concavi<-x_concavi^(-gamma)

random_var<-matrix(rnorm(Maturi*Frequ*N_Simulations),N_Simulations,Maturi*Frequ) #Matrice actif sans risque
pre1_S_tilde_t<-S_0*exp((alpha_tilde-0.5*sigma^2)/Frequ+sigma*sqrt(1/Frequ)*random_var)
pre2_S_tilde_t<-cbind((rep(S_0,N_Simulations)),pre1_S_tilde_t)

pre1_xi_tilde_t<-exp(-(r_no_risk_tilde+0.5*theta_sim_tilde^2)/Frequ-theta_sim_tilde*sqrt(1/Frequ)*random_var)
pre2_xi_tilde_t<-cbind((rep(S_0,N_Simulations)),pre1_xi_tilde_t)

pre1_B_tilde_t<-rep(exp(r_no_risk_tilde/Frequ),Frequ*Maturi) #vecteur actif sans risque
B_tilde_t<-B_0*cumprod(c(B_0,pre1_B_tilde_t))
#########################################################################


########### SECTION 1-) "Proportion optimale", avec la méthode martingale ###########

E_utility_martingale<-function(matrice_pre2_S,matrice_pre2_xi_tilde_t,vecteur_B){
  #initialisation
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
  
  for (n in 1:(N_Simulations)){ #Loop sur les simulations
    matrice_S[n,]<-cumprod(matrice_pre2_S[n,])
    matrice_xi_tilde[n,]<-cumprod(matrice_pre2_xi_tilde_t[n,])
    
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
  #U_ptf<-mean(P_Utility(X_Tu=processus_ptf[,(Frequ*Maturi+1)],gamma=gamma))
    
  #Uty<-as.numeric(lapply(funds_d,function(x)P_Utility(X_Tu=x,gamma=gamma)))
  #EU<-mean(Uty)
  Utymod<-round(mean(as.numeric(lapply(funds_d,function(x)P_Utility(X_Tu=x,gamma=2)))),3) # l'utlité est mesurée pour gamma=2.
  #Uprocessus_ptf<-P_Utility(X_Tu=processus_ptf,gamma=gamma)
  #verif<-mean(matrice_xi_tilde[,(Frequ*Maturi+1)]*matrice_S[,(Frequ*Maturi+1)])
  
  return(c(Utymod,round(Inverse_P_Utility(Utymod,2),3)))#c(exercice_guarantie,EU,U_ptf,CB)processus_ptf[,(Frequ*Maturi+1)]c(EU,CB)colMeans(Uprocessus_ptf)c(verif,CB)exercice_guarantie
}

timer<-proc.time()

#U5moyenneMM<-
  E_utility_martingale(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t)#ptf_terminal_MMB

proc.time()-timer


######################################################################


##### SECTION 2-) "Proportion optimale", avec la méthode martingale: proportion born?e entre [0,1] #######

E_utility_martingale_Borne<-function(matrice_pre2_S,matrice_pre2_xi_tilde_t,vecteur_B){
  #initialisation
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
  
  for (n in 1:(N_Simulations)){ #Loop sur les simulations
    matrice_S[n,]<-cumprod(matrice_pre2_S[n,])
    matrice_xi_tilde[n,]<-cumprod(matrice_pre2_xi_tilde_t[n,])
    
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
  #U_ptf<-mean(P_Utility(X_Tu=processus_ptf[,(Frequ*Maturi+1)],gamma=gamma))
  
  Utymod<-round(mean(as.numeric(lapply(funds_d,function(x)P_Utility(X_Tu=x,gamma=2)))),3) # l'utlité est mesurée pour gamma=2.
  #Uty<-as.numeric(lapply(funds_d,function(x)P_Utility(X_Tu=x,gamma=gamma)))
  #EU<-mean(Uty)
  #Uprocessus_ptf<-P_Utility(X_Tu=processus_ptf,gamma=gamma)
  #verif<-mean(matrice_xi_tilde[,(Frequ*Maturi+1)]*matrice_S[,(Frequ*Maturi+1)])
  
  return(c(Utymod,round(Inverse_P_Utility(Utymod,2),3)))#c(exercice_guarantie,EU,U_ptf,CB)processus_ptf[,(Frequ*Maturi+1)]c(EU,CB)colMeans(Uprocessus_ptf)c(verif,CB)exercice_guarantie
}

timer<-proc.time()

#U5moyenneMM<-
E_utility_martingale_Borne(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t)#ptf_terminal_MMB

proc.time()-timer
#############################################################################################################


########### SECTION 3-) "Proportion" dans l'actif risqué est constante (repré. par un paramètre) ###########

E_utility_prop_cte<-function(matrice_pre2_S,matrice_pre2_xi_tilde_t,vecteur_B,prop_act_r){
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
  
  for (n in 1:(N_Simulations)){ #Loop sur les simulations
    matrice_S[n,]<-cumprod(matrice_pre2_S[n,])
    matrice_xi_tilde[n,]<-cumprod(matrice_pre2_xi_tilde_t[n,])
    
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
  #U_ptf<-mean(P_Utility(X_Tu=processus_ptf[,(Frequ*Maturi+1)],gamma=gamma))
  Utymod<-round(mean(as.numeric(lapply(funds_d,function(x)P_Utility(X_Tu=x,gamma=2)))),3) # l'utlité est mesurée pour gamma=2.
  
  #EU<-mean(Uty)
  
  #xi_tilde[n,(Frequ*Maturi+1)]<-matrice_S[n,(Frequ*Maturi+1)]^(-theta_sim_tilde/sigma)*exp((theta_sim_tilde*alpha_tilde/sigma-theta_sim_tilde*sigma/2-r_no_risk_tilde-0.5*theta_sim_tilde^2)*petit_t)
  
  #contrainte_budget<-mean(matrice_xi_tilde[,(Frequ*Maturi+1)]*processus_ptf[,(Frequ*Maturi+1)])
  #Uprocessus_ptf<-P_Utility(X_Tu=processus_ptf,gamma=gamma)
  #verif<-mean(matrice_xi_tilde[,(Frequ*Maturi+1)]*matrice_S[,(Frequ*Maturi+1)])
  
  return(c(Utymod,round(Inverse_P_Utility(Utymod,2),3)))#c(exercice_guarantie,EU,U_ptf,contrainte_budget)c(EU,contrainte_budget)colMeans(Uprocessus_ptf)c(verif,contrainte_budget)exercice_guarantieprocessus_ptf[,(Frequ*Maturi+1)]
}

timer2<-proc.time()

#U5moyenne02<-
E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=cte_Merton)

proc.time()-timer2


########### SECTION 4-) Simulations du portefeuille optimal directement ? maturit? ###########

Simulations_fonds_distinct<-function(matrice_pre2_S){
  matrice_S<-matrix(0,N_Simulations,Frequ*Maturi+1)
  xi_tilde<-matrix(0,N_Simulations,Frequ*Maturi+1)
  funds_d<-rep(0,N_Simulations)
  
  for (n in 1:(N_Simulations)){ #Loop sur les simulations
    matrice_S[n,]<-cumprod(matrice_pre2_S[n,])
    
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
  Utymod<-round(mean(as.numeric(lapply(funds_d,function(x)P_Utility(X_Tu=x,gamma=2)))),3) # l'utlité est mesurée pour gamma=2.
  #EU<-mean(as.numeric(lapply(funds_d,function(x)P_Utility(X_Tu=x,gamma=gamma))))
  #CB<-mean(ptf_optimal*xi_tilde[,(Frequ*Maturi+1)])
  #verif<-mean(matrice_S[,(Frequ*Maturi+1)]*xi_tilde[,(Frequ*Maturi+1)])
  return(c(Utymod,round(Inverse_P_Utility(Utymod,2),3)))#c(CB,EU) verifc(CB,EU,verif,exercice_guarantie)
}

timer3<-proc.time()

#ptf_terminal_final<-
  Simulations_fonds_distinct(pre2_S_tilde_t)#funds_final, 37.06 sec ? rouler

proc.time()-timer3








####################### ANALYSE ###########################

          timer4<-proc.time()
ptf_terminal_MMB<-E_utility_martingale_Borne(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t)
          proc.time()-timer4
ptf_terminal_MM<-E_utility_martingale(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t)
          proc.time()-timer4
ptf_terminal_100<-E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=1.0)
            proc.time()-timer4
ptf_terminal_60<-E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=0.6)
            proc.time()-timer4
ptf_terminal_40<-E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=0.4)
            proc.time()-timer4
ptf_terminal_20<-E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=0.2)
            proc.time()-timer4
ptf_terminal_Merton<-E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=cte_Merton)
            proc.time()-timer4
ptf_terminal_final<-Simulations_fonds_distinct(pre2_S_tilde_t)




############## Graphique: densit? de toutes les strat?gies ########################################################
ajust_grap<-2.5 #si ma d'ajustemenbt, alors mettre ex:20. 


funds_rebalancement_MMB<-data.frame(funds_optimal=as.numeric(lapply(ptf_terminal_MMB,function(x)min(x,ajust_grap))))
funds_rebalancement_MM<-data.frame(funds_optimal=as.numeric(lapply(ptf_terminal_MM,function(x)min(x,ajust_grap))))
funds_rebalancement_100<-data.frame(funds_optimal=as.numeric(lapply(ptf_terminal_100,function(x)min(x,ajust_grap))))
funds_rebalancement_60<-data.frame(funds_optimal=as.numeric(lapply(ptf_terminal_60,function(x)min(x,ajust_grap))))
funds_rebalancement_40<-data.frame(funds_optimal=as.numeric(lapply(ptf_terminal_40,function(x)min(x,ajust_grap))))
funds_rebalancement_20<-data.frame(funds_optimal=as.numeric(lapply(ptf_terminal_20,function(x)min(x,ajust_grap))))
funds_rebalancement_Merton<-data.frame(funds_optimal=as.numeric(lapply(ptf_terminal_Merton,function(x)min(x,ajust_grap))))
funds_final<-data.frame(funds_optimal=as.numeric(lapply(ptf_terminal_final,function(x)min(x,ajust_grap))))

funds_rebalancement_100$comp<-'1-100% risque'
funds_rebalancement_MMB$comp<-'6-Optimale, [0,1]'
funds_rebalancement_MM$comp<-'7-Optimale, [0,2]'
funds_rebalancement_60$comp<-'2-60% risque'
funds_rebalancement_40$comp<-'3-40% risque'
funds_rebalancement_20$comp<-'4-20% risque'
funds_rebalancement_Merton$comp<-'5-Constante Merton'
funds_final$comp<-'8-Optimale'


Compa_combin_funds<-rbind(funds_rebalancement_100,funds_rebalancement_MMB,funds_rebalancement_MM,funds_rebalancement_60,funds_rebalancement_40,funds_rebalancement_20,funds_rebalancement_Merton,funds_final)
#max(funds_rebalancement_100$funds_optimal)


ggplot(data=Compa_combin_funds,aes(Compa_combin_funds$funds_optimal,group=comp,fill=comp,color=comp))+
  stat_ecdf(geom='step')+
  #facet_wrap(facets = vars(comp))+
  labs(x=expression(paste("F",""[T])), y="y")+
  theme(legend.position = 'bottom',legend.title = element_blank())

# ggplot(data=Compa_combin_funds,aes(Compa_combin_funds$funds_optimal,group=comp,fill=comp))+
#   geom_histogram(colour='black',alpha=0.5,position = "identity")
####################################################################################







############### Comparaison uniquement du fonds optimal et martingale [0,1] ##########################################
ajust_grap2<-2.5

funds_final_2<-data.frame(funds_optimal=as.numeric(lapply(ptf_terminal_final,function(x)min(x,ajust_grap2))))
funds_rebalancement_MMB_2<-data.frame(funds_optimal=as.numeric(lapply(ptf_terminal_MMB,function(x)min(x,ajust_grap2))))

funds_final_2$comp<-'Optimale'
funds_rebalancement_MMB_2$comp<-'Optimale, [0,1]'

Compa_combin_funds_uniq2<-rbind(funds_rebalancement_MMB_2,funds_final_2)

ggplot(data=Compa_combin_funds_uniq2,aes(Compa_combin_funds_uniq2$funds_optimal,group=comp,fill=comp))+
  geom_histogram(colour='black',bins=100,alpha=0.5,position = "identity")+
  #facet_wrap(facets = vars(comp))+
  labs(x=expression(paste("F",""[T])), y="# R?alisations")+
  theme(legend.position = 'bottom',legend.title = element_blank())
 # theme_classic()
#  scale_fill_manual(name="",values=c("yellow","darkgray"),labels=c( expression(paste("X","*"[T],": Dynamique")),expression(paste("X","*"[T],": Th?orique"))))+


hist(funds_final)
hist(funds_rebalancement)



############## Valeurs moyennes du portefeuille et du fonds #############
par(mfrow=c(1,2))

plot(1:521,moyenneMM,xlab='Ann?es (t)', ylab='Valeur moyenne du portefeuille',type='l',lty=1,col='royalblue1',ylim=c(0.9,1.6),xaxt="n")#,lty='twodash'
xtick<-c(0,52,52*2,52*3,52*4,52*5,52*6,52*7,52*8,52*9,52*10)
axis(side=1, at=xtick, labels =c(0,1,2,3,4,5,6,7,8,9,10))
lines(1:521,moyenneMMB,col='purple',lty=2)#,lty='longdash'
lines(1:521,moyenne02,col='mediumseagreen',lty=6)
lines(1:521,moyenne04,col='darkolivegreen4',lty=5) #lty='dotted'
lines(1:521,moyenne06,col='goldenrod3',lty=4)
lines(1:521,moyenne1,col='red',lty=3)
lines(1:521,moyenneM,col='limegreen',lty=2)
legend(50,1.6, legend=c("Optimale [0,2]", "Optimale [0,1]",'100% actif risqu?','60% actif risqu?','40% actif risqu?','20% actif risqu?','Constante de Merton'),
       col=c('royalblue1', 'purple','red','goldenrod3','darkolivegreen4','mediumseagreen','limegreen'),lty=c(1,2,3,4,5,6,2), cex=0.8)

plot(1:521,fmoyenneMM,xlab='Ann?es (t)', ylab='Valeur moyenne du fonds distinct',type='l',lty=1,col='royalblue1',ylim=c(0.9,1.65),xaxt="n")#,lty='twodash'
xtick<-c(0,52,52*2,52*3,52*4,52*5,52*6,52*7,52*8,52*9,52*10)
axis(side=1, at=xtick, labels =c(0,1,2,3,4,5,6,7,8,9,10))
lines(1:521,fmoyenneMMB,col='purple',lty=2)#,lty='longdash'
lines(1:521,fmoyenne02,col='mediumseagreen',lty=6)
lines(1:521,fmoyenne04,col='darkolivegreen4',lty=5) #lty='dotted'
lines(1:521,fmoyenne06,col='goldenrod3',lty=4)
lines(1:521,fmoyenne1,col='red',lty=3)
lines(1:521,fmoyenneM,col='limegreen',lty=2)
legend(50,1.6, legend=c("Optimale [0,2]", "Optimale [0,1]",'100% actif risqu?','60% actif risqu?','40% actif risqu?','20% actif risqu?','Constante de Merton'),
       col=c('royalblue1', 'purple','red','goldenrod3','darkolivegreen4','mediumseagreen','limegreen'),lty=c(1,2,3,4,5,6,2), cex=0.8)


############## Utilit?s moyennes du portefeuille et du fonds #############
par(mfrow=c(1,2))

plot(1:521,UmoyenneMM,xlab='Ann?es (t)', ylab=expression(paste("Esp?rance de l'utilit? du portefeuille")),type='l',lty=1,col='royalblue1',ylim=c(-2.4,-0.8),xaxt="n")#,lty='twodash'
xtick<-c(0,52,52*2,52*3,52*4,52*5,52*6,52*7,52*8,52*9,52*10)
axis(side=1, at=xtick, labels =c(0,1,2,3,4,5,6,7,8,9,10))
lines(1:521,UmoyenneMMB,col='purple',lty=2)#,lty='longdash'
lines(1:521,Umoyenne02,col='mediumseagreen',lty=6)
lines(1:521,Umoyenne04,col='darkolivegreen4',lty=5) #lty='dotted'
lines(1:521,Umoyenne06,col='goldenrod3',lty=4)
lines(1:521,Umoyenne1,col='red',lty=3)
lines(1:521,UmoyenneM,col='limegreen',lty=2)
legend(50,-1.6, legend=c("Optimale [0,2]", "Optimale [0,1]",'100% actif risqu?','60% actif risqu?','40% actif risqu?','20% actif risqu?','Constante de Merton'),
       col=c('royalblue1', 'purple','red','goldenrod3','darkolivegreen4','mediumseagreen','limegreen'),lty=c(1,2,3,4,5,6,2), cex=0.8)

plot(1:521,UfmoyenneMM,xlab='Ann?es (t)', ylab="Esp?rance de l'utilit? du fonds distinct",type='l',lty=1,col='royalblue1',ylim=c(-1,-0.6),xaxt="n")#,lty='twodash'
xtick<-c(0,52,52*2,52*3,52*4,52*5,52*6,52*7,52*8,52*9,52*10)
axis(side=1, at=xtick, labels =c(0,1,2,3,4,5,6,7,8,9,10))
lines(1:521,UfmoyenneMMB,col='purple',lty=2)#,lty='longdash'
lines(1:521,Ufmoyenne02,col='mediumseagreen',lty=6)
lines(1:521,Ufmoyenne04,col='darkolivegreen4',lty=5) #lty='dotted'
lines(1:521,Ufmoyenne06,col='goldenrod3',lty=4)
lines(1:521,Ufmoyenne1,col='red',lty=3)
lines(1:521,UfmoyenneM,col='limegreen',lty=2)
legend(50,-0.6, legend=c("Optimale [0,2]", "Optimale [0,1]",'100% actif risqu?','60% actif risqu?','40% actif risqu?','20% actif risqu?','Constante de Merton'),
       col=c('royalblue1', 'purple','red','goldenrod3','darkolivegreen4','mediumseagreen','limegreen'),lty=c(1,2,3,4,5,6,2), cex=0.8)

de<-matrix(1:6,3,2)
eee<-apply(de,c(1,2),function(x) max(x,2))
colMeans(eee)
############## Gamma-Modifi?: Utilit?s moyennes du portefeuille et du fonds #############
par(mfrow=c(1,2))

plot(1:521,U5moyenneMM,xlab='Ann?es (t)', ylab=expression(paste("Esp?rance de l'utilit? du portefeuille")),type='l',ylim=c(-3,-0.2),lty=1,col='royalblue1',xaxt="n")#,lty='twodash'
xtick<-c(0,52,52*2,52*3,52*4,52*5,52*6,52*7,52*8,52*9,52*10)
axis(side=1, at=xtick, labels =c(0,1,2,3,4,5,6,7,8,9,10))
lines(1:521,U5moyenneMMB,col='purple',lty=2)#,lty='longdash'
lines(1:521,U5moyenne02,col='mediumseagreen',lty=6)
lines(1:521,U5moyenne04,col='darkolivegreen4',lty=5) #lty='dotted'
lines(1:521,U5moyenne06,col='goldenrod3',lty=4)
lines(1:521,U5moyenne1,col='red',lty=3)
lines(1:521,U5moyenneM,col='limegreen',lty=2)
legend(0,-2.0, legend=c("Optimale [0,2]", "Optimale [0,1]",'100% actif risqu?','60% actif risqu?','40% actif risqu?','20% actif risqu?','Constante de Merton'),
       col=c('royalblue1', 'purple','red','goldenrod3','darkolivegreen4','mediumseagreen','limegreen'),lty=c(1,2,3,4,5,6,2), cex=0.8)

plot(1:521,U5fmoyenneMM,xlab='Ann?es (t)', ylab="Esp?rance de l'utilit? du fonds distinct",type='l',lty=1,col='royalblue1',xaxt="n")#,lty='twodash'
xtick<-c(0,52,52*2,52*3,52*4,52*5,52*6,52*7,52*8,52*9,52*10)
axis(side=1, at=xtick, labels =c(0,1,2,3,4,5,6,7,8,9,10))
lines(1:521,U5fmoyenneMMB,col='purple',lty=2)#,lty='longdash'
lines(1:521,U5fmoyenne02,col='mediumseagreen',lty=6)
lines(1:521,U5fmoyenne04,col='darkolivegreen4',lty=5) #lty='dotted'
lines(1:521,U5fmoyenne06,col='goldenrod3',lty=4)
lines(1:521,U5fmoyenne1,col='red',lty=3)
lines(1:521,U5fmoyenneM,col='limegreen',lty=2)
legend(300,-0.2, legend=c("Optimale [0,2]", "Optimale [0,1]",'100% actif risqu?','60% actif risqu?','40% actif risqu?','20% actif risqu?','Constante de Merton'),
       col=c('royalblue1', 'purple','red','goldenrod3','darkolivegreen4','mediumseagreen','limegreen'),lty=c(1,2,3,4,5,6,2), cex=0.8)



