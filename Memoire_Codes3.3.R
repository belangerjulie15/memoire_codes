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
fee_c_s<-0.0 #Fee applied of the risky asset
fee_c_f<-0.02448 #Fee applied of the funds 
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

pre1_S_t<-S_0*exp((alpha-0.5*sigma^2)/Frequ+sigma*sqrt(1/Frequ)*random_var)
pre2_S_t<-cbind((rep(S_0,N_Simulations)),pre1_S_t)

pre1_xi_t<-exp(-(r_no_risk+0.5*theta_sim^2)/Frequ-theta_sim*sqrt(1/Frequ)*random_var)
pre2_xi_t<-cbind((rep(S_0,N_Simulations)),pre1_xi_t)
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
  cout_guar_ass<-rep(0,N_Simulations)
  Uprocessus_ptf<-matrix(0,N_Simulations,Frequ*Maturi+1)
  matrice_xi_tilde<-matrix(0,N_Simulations,Frequ*Maturi+1)
  matrice_xi<-matrix(0,N_Simulations,Frequ*Maturi+1)
  
  for (n in 1:(N_Simulations)){ #Loop sur les simulations
    matrice_S[n,]<-cumprod(matrice_pre2_S[n,])
    matrice_xi_tilde[n,]<-cumprod(matrice_pre2_xi_tilde_t[n,])
    matrice_xi[n,]<-cumprod(pre2_xi_t[n,])
    
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
    #funds_d[n]<-a_call_sim*max(0,processus_ptf[n,(Frequ*Maturi+1)]-b_call_sim)+K_call_sim
    #cout_guar_ass[n]<-max(0,b_call_sim-processus_ptf[n,(Frequ*Maturi+1)])
  }
  
  #CAss<-mean(matrice_xi[,(Frequ*Maturi+1)]*funds_d)
  #CB<-mean(matrice_xi_tilde[,(Frequ*Maturi+1)]*processus_ptf[,(Frequ*Maturi+1)])
  #exercice_guarantie<-sum(processus_ptf[,(Frequ*Maturi+1)]<b_call_sim)/N_Simulations
  #Esp_cout_garantie<-mean(cout_guar_ass)
  #call_tout_t<-apply(processus_ptf,c(1,2),function(x) a_call_sim*max(0,x-b_call_sim)+K_call_sim)
  Uty_t<-apply(processus_ptf,c(1,2),function(x) P_Utility(X_Tu=x,gamma=gamma))#utilit? des valeurs du ptf pour tout t
  #Uty_f_t<-apply(call_tout_t,c(1,2),function(x) P_Utility(X_Tu=x,gamma=gamma))#utilit? des valeurs du fonds pour tout t
  #U_ptf<-mean(P_Utility(X_Tu=processus_ptf[,(Frequ*Maturi+1)],gamma=gamma))
    
  #Uty<-as.numeric(lapply(funds_d,function(x)P_Utility(X_Tu=x,gamma=gamma)))
  #EU<-mean(Uty)
  #Utymod<-round(mean(as.numeric(lapply(funds_d,function(x)P_Utility(X_Tu=x,gamma=3)))),3) # l'utlité est mesurée pour gamma=3.
  #Uprocessus_ptf<-P_Utility(X_Tu=processus_ptf,gamma=gamma)
  #verif<-mean(matrice_xi_tilde[,(Frequ*Maturi+1)]*matrice_S[,(Frequ*Maturi+1)])
  #ve_a_de<-mean(nu_Act_R>=1)
  #ve_a_2_de<-mean(nu_Act_R>=2)
    
  return(colMeans(Uty_t))#c(round(exercice_guarantie,3),round(Esp_cout_garantie,3))c(Utymod,round(Inverse_P_Utility(Utymod,7),3))#c(exercice_guarantie,EU,U_ptf,CB)processus_ptf[,(Frequ*Maturi+1)]c(EU,CB)colMeans(Uprocessus_ptf)c(verif,CB)exercice_guarantie
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
  cout_guar_ass<-rep(0,N_Simulations)
  Uprocessus_ptf<-matrix(0,N_Simulations,Frequ*Maturi+1)
  matrice_xi_tilde<-matrix(0,N_Simulations,Frequ*Maturi+1)
  matrice_xi<-matrix(0,N_Simulations,Frequ*Maturi+1)
  
  for (n in 1:(N_Simulations)){ #Loop sur les simulations
    matrice_S[n,]<-cumprod(matrice_pre2_S[n,])
    matrice_xi_tilde[n,]<-cumprod(matrice_pre2_xi_tilde_t[n,])
    matrice_xi[n,]<-cumprod(pre2_xi_t[n,])
    
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
    #cout_guar_ass[n]<-max(0,b_call_sim-processus_ptf[n,(Frequ*Maturi+1)])
  }
  
  #CAss<-mean(matrice_xi[,(Frequ*Maturi+1)]*funds_d)
  #CB<-mean(matrice_xi_tilde[,(Frequ*Maturi+1)]*processus_ptf[,(Frequ*Maturi+1)])
  #exercice_guarantie<-sum(processus_ptf[,(Frequ*Maturi+1)]<b_call_sim)/N_Simulations  
  #Esp_cout_garantie<-mean(cout_guar_ass)
  #call_tout_t<-apply(processus_ptf,c(1,2),function(x) a_call_sim*max(0,x-b_call_sim)+K_call_sim)
  #Uty_t<-apply(processus_ptf,c(1,2),function(x) P_Utility(X_Tu=x,gamma=gamma))#utilit? des valeurs du ptf pour tout t
  #Uty_f_t<-apply(call_tout_t,c(1,2),function(x) P_Utility(X_Tu=x,gamma=gamma))#utilit? des valeurs du fonds pour tout t
  #U_ptf<-mean(P_Utility(X_Tu=processus_ptf[,(Frequ*Maturi+1)],gamma=gamma))
  
  #Utymod<-round(mean(as.numeric(lapply(funds_d,function(x)P_Utility(X_Tu=x,gamma=3)))),3) # l'utlité est mesurée pour gamma=7.
  #Uty<-as.numeric(lapply(funds_d,function(x)P_Utility(X_Tu=x,gamma=gamma)))
  #EU<-mean(Uty)
  #Uprocessus_ptf<-P_Utility(X_Tu=processus_ptf,gamma=gamma)
  #verif<-mean(matrice_xi_tilde[,(Frequ*Maturi+1)]*matrice_S[,(Frequ*Maturi+1)])
  
  return(funds_d)#c(round(exercice_guarantie,3),round(Esp_cout_garantie,3))c(round(exercice_guarantie,3),round(Esp_cout_garantie,3))c(exercice_guarantie,EU,U_ptf,CB)processus_ptf[,(Frequ*Maturi+1)]c(EU,CB)colMeans(Uprocessus_ptf)c(verif,CB)exercice_guarantie
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
  cout_guar_ass<-rep(0,N_Simulations)
  Uprocessus_ptf<-matrix(0,N_Simulations,Frequ*Maturi+1)
  matrice_xi_tilde<-matrix(0,N_Simulations,Frequ*Maturi+1)
  matrice_xi<-matrix(0,N_Simulations,Frequ*Maturi+1)
  
  for (n in 1:(N_Simulations)){ #Loop sur les simulations
    matrice_S[n,]<-cumprod(matrice_pre2_S[n,])
    matrice_xi_tilde[n,]<-cumprod(matrice_pre2_xi_tilde_t[n,])
    matrice_xi[n,]<-cumprod(pre2_xi_t[n,])
    
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
    #funds_d[n]<-a_call_sim*max(0,processus_ptf[n,(Frequ*Maturi+1)]-b_call_sim)+K_call_sim
    #cout_guar_ass[n]<-max(0,b_call_sim-processus_ptf[n,(Frequ*Maturi+1)])
  }
  
  #CAss<-mean(matrice_xi[,(Frequ*Maturi+1)]*funds_d)
  #exercice_guarantie<-sum(processus_ptf[,(Frequ*Maturi+1)]<b_call_sim)/N_Simulations
  #Esp_cout_garantie<-mean(cout_guar_ass)
  #call_tout_t<-apply(processus_ptf,c(1,2),function(x) a_call_sim*max(0,x-b_call_sim)+K_call_sim)
  #Uty<-as.numeric(lapply(funds_d,function(x)P_Utility(X_Tu=x,gamma=gamma)))
  Uty_t<-apply(processus_ptf,c(1,2),function(x) P_Utility(X_Tu=x,gamma=gamma))#utilit? des valeurs du ptf pour tout t
  #Uty_f_t<-apply(call_tout_t,c(1,2),function(x) P_Utility(X_Tu=x,gamma=gamma))#utilit? des valeurs du fonds pour tout t
  #U_ptf<-mean(P_Utility(X_Tu=processus_ptf[,(Frequ*Maturi+1)],gamma=gamma))
  #Utymod<-round(mean(as.numeric(lapply(funds_d,function(x)P_Utility(X_Tu=x,gamma=3)))),3) # l'utlité est mesurée pour gamma=7.
  
  #EU<-mean(Uty)
  
  #xi_tilde[n,(Frequ*Maturi+1)]<-matrice_S[n,(Frequ*Maturi+1)]^(-theta_sim_tilde/sigma)*exp((theta_sim_tilde*alpha_tilde/sigma-theta_sim_tilde*sigma/2-r_no_risk_tilde-0.5*theta_sim_tilde^2)*petit_t)
  
  #contrainte_budget<-mean(matrice_xi_tilde[,(Frequ*Maturi+1)]*processus_ptf[,(Frequ*Maturi+1)])
  #Uprocessus_ptf<-P_Utility(X_Tu=processus_ptf,gamma=gamma)
  #verif<-mean(matrice_xi_tilde[,(Frequ*Maturi+1)]*matrice_S[,(Frequ*Maturi+1)])
  
  return(colMeans(Uty_t))#c(exercice_guarantie,EU,U_ptf,contrainte_budget)c(EU,contrainte_budget)colMeans(Uprocessus_ptf)c(verif,contrainte_budget)exercice_guarantieprocessus_ptf[,(Frequ*Maturi+1)]
}#c(round(exercice_guarantie,3),round(Esp_cout_garantie,3))

timer2<-proc.time()
#U5moyenne02<-
E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=0.6)
E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=0.4)
E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=0.2)
E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=cte_Merton)

proc.time()-timer2

##### Essai paralléllisation ####
#registerDoParallel(cores=5)
#system.time(result_p<-as.numeric(foreach(i=c(1,0.6,0.4,0.2,cte_Merton)) %dopar%E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=i)))


########### SECTION 4-) Simulations du portefeuille optimal directement ? maturit? ###########

Simulations_fonds_distinct<-function(matrice_pre2_S){
  matrice_S<-matrix(0,N_Simulations,Frequ*Maturi+1)
  xi_tilde<-matrix(0,N_Simulations,Frequ*Maturi+1)
  funds_d<-rep(0,N_Simulations)
  cout_guar_ass<-rep(0,N_Simulations)
  matrice_xi<-matrix(0,N_Simulations,Frequ*Maturi+1)
  
  for (n in 1:(N_Simulations)){ #Loop sur les simulations
    matrice_S[n,]<-cumprod(matrice_pre2_S[n,])
    matrice_xi[n,]<-cumprod(pre2_xi_t[n,])
    
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
    #cout_guar_ass[n]<-max(0,b_call_sim-ptf_optimal[n])
  }
  #CAss<-mean(matrice_xi[,(Frequ*Maturi+1)]*funds_d)
  #exercice_guarantie<-sum(ptf_optimal[]<b_call_sim)/N_Simulations
  #Esp_cout_garantie<-mean(cout_guar_ass)
  
  #Utymod<-round(mean(as.numeric(lapply(funds_d,function(x)P_Utility(X_Tu=x,gamma=3)))),3) # l'utlité est mesurée pour gamma=7.
  #EU<-mean(as.numeric(lapply(funds_d,function(x)P_Utility(X_Tu=x,gamma=gamma))))
  #CB<-mean(ptf_optimal*xi_tilde[,(Frequ*Maturi+1)])
  #verif<-mean(matrice_S[,(Frequ*Maturi+1)]*xi_tilde[,(Frequ*Maturi+1)])
  return(funds_d)#funds_d,c(CB,EU) verifc(CB,EU,verif,exercice_guarantie)
}#c(round(exercice_guarantie,3),round(Esp_cout_garantie,3))

timer3<-proc.time()

#ptf_terminal_final<-
#ptf_18_6448
#ptf_0_0<-Simulations_fonds_distinct(pre2_S_tilde_t)#funds_final, 37.06 sec ? rouler
#ptf_2448_0
#ptf_1224_1224
#ptf_0_0
Simulations_fonds_distinct(pre2_S_tilde_t)
#fait#cs0<-Simulations_fonds_distinct(pre2_S_tilde_t)
#fait#cs05<-
#fait#cs1<-
#fait#cs15<-
  
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
            proc.time()-timer4

#60 minutes avec Kronos

############## Graphique: densit? de toutes les strat?gies ########################################################
ajust_grap<-2.5 #si pas d'ajustement, alors mettre ex:20. 


funds_rebalancement_MMB<-data.frame(funds_optimal=as.numeric(lapply(ptf_terminal_MMB,function(x)min(x,ajust_grap))))
funds_rebalancement_MM<-data.frame(funds_optimal=as.numeric(lapply(ptf_terminal_MM,function(x)min(x,ajust_grap))))
funds_rebalancement_100<-data.frame(funds_optimal=as.numeric(lapply(ptf_terminal_100,function(x)min(x,ajust_grap))))
funds_rebalancement_60<-data.frame(funds_optimal=as.numeric(lapply(ptf_terminal_60,function(x)min(x,ajust_grap))))
funds_rebalancement_40<-data.frame(funds_optimal=as.numeric(lapply(ptf_terminal_40,function(x)min(x,ajust_grap))))
funds_rebalancement_20<-data.frame(funds_optimal=as.numeric(lapply(ptf_terminal_20,function(x)min(x,ajust_grap))))
funds_rebalancement_Merton<-data.frame(funds_optimal=as.numeric(lapply(ptf_terminal_Merton,function(x)min(x,ajust_grap))))
funds_final<-data.frame(funds_optimal=as.numeric(lapply(ptf_terminal_final,function(x)min(x,ajust_grap))))

funds_rebalancement_100$comp<-'1-100% risque'
funds_rebalancement_MMB$comp<-'6-Optimale [0,1]'
funds_rebalancement_MM$comp<-'7-Optimale [0,2]'
funds_rebalancement_60$comp<-'2-60% risque'
funds_rebalancement_40$comp<-'3-40% risque'
funds_rebalancement_20$comp<-'4-20% risque'
funds_rebalancement_Merton$comp<-'5-Constante Merton'
funds_final$comp<-'8-Optimale'


Compa_combin_funds<-rbind(funds_rebalancement_100,funds_rebalancement_MMB,funds_rebalancement_MM,funds_rebalancement_60,funds_rebalancement_40,funds_rebalancement_20,funds_rebalancement_Merton,funds_final)
#max(funds_rebalancement_100$funds_optimal)


ggplot(data=Compa_combin_funds,aes(Compa_combin_funds$funds_optimal,group=comp,fill=comp,linetype=comp,color=comp))+
  stat_ecdf(geom='step')+
  #facet_wrap(facets = vars(comp))+
  labs(x=expression(paste("F",""[T])), y="Distribution cumulative")+
  theme(legend.position = 'bottom',legend.title = element_blank())

# ggplot(data=Compa_combin_funds,aes(Compa_combin_funds$funds_optimal,group=comp,fill=comp))+
#   geom_histogram(colour='black',alpha=0.5,position = "identity")
####################################################################################




#Ne pas toucher#funds1224_1224_MM<-funds_final$funds_optimal
#Ne pas toucher#funds1224_1224_MM01<-funds_rebalancement_MMB$funds_optimal  


############### Comparaison uniquement du fonds optimal et martingale [0,1] ##########################################
ajust_grap2<-2.5

funds_final_2<-data.frame(funds_optimal=as.numeric(lapply(ptf_terminal_final,function(x)min(x,ajust_grap2))))
funds_rebalancement_MMB_2<-data.frame(funds_optimal=as.numeric(lapply(ptf_terminal_MMB,function(x)min(x,ajust_grap2))))

funds_final_2$comp<-'Optimale'
funds_rebalancement_MMB_2$comp<-'Optimale [0,1]'

Compa_combin_funds_uniq2<-rbind(funds_final_2,funds_rebalancement_MMB_2)

ggplot(data=Compa_combin_funds_uniq2,aes(Compa_combin_funds_uniq2$funds_optimal,group=comp,fill=comp))+
  geom_histogram(colour='black',bins=80,alpha=0.6,position = "identity")+
  #facet_wrap(facets = vars(comp))+
  #scale_colour_gradient(low='#FF61CC',high='blue')+
  scale_fill_grey()+
  labs(x=expression(paste("(F",""[T]," - G)"^"+","+G")), y="Réalisations")+
  theme_classic()+
  theme(legend.position = 'bottom' ,legend.title = element_blank())
 #legend.text=element_text(size=12)
#  scale_fill_manual(name="",values=c("yellow","darkgray"),labels=c( expression(paste("X","*"[T],": Dynamique")),expression(paste("X","*"[T],": Th?orique"))))+


hist(funds_final)
hist(funds_rebalancement)



############## Valeurs moyennes du portefeuille et du fonds #############
########################## À travers le temps ########################
timer10<-proc.time()
moyenneMM<-E_utility_martingale(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t)#ptf_MM
moyenneMMB<-E_utility_martingale_Borne(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t)#ptf_MMB
moyenne02<-E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=0.2)
moyenne04<-E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=0.4)
moyenne06<-E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=0.6)
moyenne1<-E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=1)
moyenneM<-E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=cte_Merton)
proc.time()-timer10 # Kronos:3734.05  sec

fmoyenneMM<-E_utility_martingale(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t)#ptf_MM
fmoyenneMMB<-E_utility_martingale_Borne(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t)#ptf_MMB
fmoyenne02<-E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=0.2)
fmoyenne04<-E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=0.4)
fmoyenne06<-E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=0.6)
fmoyenne1<-E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=1)
fmoyenneM<-E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=cte_Merton)

  
ggp_mMM<-data.frame(ptf_optimal=moyenneMM)
ggp_mMMB<-data.frame(ptf_optimal=moyenneMMB)
ggp_m02<-data.frame(ptf_optimal=moyenne02)
ggp_m04<-data.frame(ptf_optimal=moyenne04)
ggp_m06<-data.frame(ptf_optimal=moyenne06)
ggp_m1<-data.frame(ptf_optimal=moyenne1)
ggp_mM<-data.frame(ptf_optimal=moyenneM)

ggp_mfMM<-data.frame(fonds_optimal=fmoyenneMM)
ggp_mfMMB<-data.frame(fonds_optimal=fmoyenneMMB)
ggp_mf02<-data.frame(fonds_optimal=fmoyenne02)
ggp_mf04<-data.frame(fonds_optimal=fmoyenne04)
ggp_mf06<-data.frame(fonds_optimal=fmoyenne06)
ggp_mf1<-data.frame(fonds_optimal=fmoyenne1)
ggp_mfM<-data.frame(fonds_optimal=fmoyenneM)

ggp_mMM$method<-"7-Optimale [0,2]"
ggp_mMMB$method<-"6-Optimale [0,1]"
ggp_m02$method<-'4-20% risque'
ggp_m04$method<-'3-40% risque'
ggp_m06$method<-'2-60% risque'
ggp_m1$method<-'1-100% risque'
ggp_mM$method<-'5-Constante de Merton'

ggp_mfMM$method<-"7-Optimale [0,2]"
ggp_mfMMB$method<-"6-Optimale [0,1]"
ggp_mf02$method<-'4-20% risque'
ggp_mf04$method<-'3-40% risque'
ggp_mf06$method<-'2-60% risque'
ggp_mf1$method<-'1-100% risque'
ggp_mfM$method<-'5-Constante de Merton'

ggpm_tot<-rbind(ggp_mMM,ggp_mMMB,ggp_m02,ggp_m04,ggp_m06,ggp_m1,ggp_mM)
ggpm_tot$temps<-seq(0,10,1/52)

ggpmf_tot<-rbind(ggp_mfMM,ggp_mfMMB,ggp_mf02,ggp_mf04,ggp_mf06,ggp_mf1,ggp_mfM)
ggpmf_tot$temps<-seq(0,10,1/52)


# -Graphiques- #

#ggplot2: portfeuille moyen à travers le temps #
ggplot(data=ggpm_tot,aes(x=ggpm_tot$temps,y=ggpm_tot$ptf_optimal,group=method,color=method))+
  geom_line(aes(linetype=method))+ 
  scale_color_manual(name="",values=c("#F8766D","#CD9600","#7CAE00","#00BE67","#00BFC4","#00A9FF","#C77CFF"),labels=c('1-100% risque','2-60% risque','3-40% risque','4-20% risque','5-Constante de Merton',"6-Optimale [0,1]","7-Optimale [0,2]"))+
  labs(x=expression(paste(t,"    (années)")), y=expression(paste("E"^P,"[ U( ","F"[t]," ) ]")))+
  scale_linetype_discrete(name="",
                          breaks=c('1-100% risque','2-60% risque','3-40% risque','4-20% risque','5-Constante de Merton',"6-Optimale [0,1]","7-Optimale [0,2]"),
                          labels=c('1-100% risque','2-60% risque','3-40% risque','4-20% risque','5-Constante de Merton',"6-Optimale [0,1]","7-Optimale [0,2]"))+
  theme(legend.position = 'bottom',legend.title = element_blank())##C77CFF
  

#ggplot2: fonds moyen à travers le temps #
ggplot(data=ggpmf_tot,aes(x=ggpmf_tot$temps,y=ggpmf_tot$fonds_optimal,group=method,color=method))+
  geom_line(aes(linetype=method))+
  scale_color_manual(name="",values=c("#F8766D","#CD9600","#7CAE00","#00BE67","#00BFC4","#00A9FF","#C77CFF"),labels=c('1-100% risque','2-60% risque','3-40% risque','4-20% risque','5-Constante de Merton',"6-Optimale [0,1]","7-Optimale [0,2]"))+
  scale_linetype_discrete(name="",
                          breaks=c('1-100% risque','2-60% risque','3-40% risque','4-20% risque','5-Constante de Merton',"6-Optimale [0,1]","7-Optimale [0,2]"),                          labels=c('1-100% risque','2-60% risque','3-40% risque','4-20% risque','5-Constante de Merton',"6-Optimale [0,1]","7-Optimale [0,2]"))+
  theme(legend.position = 'bottom',legend.title = element_blank())+##C77CFF
  labs(x=expression(paste(t,"    (années)")), y=expression(paste("E"^P,"[ U( (","F"[t],"-1)"^"+","+1",") ]")))


#plot
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





######### Graphique de l'espérance de l'utilité du bénéfice résultant du fonds distinct#########
############## en fonction des proportions cte investies dans le fonds #########################

system.time(res<-lapply(seq(0,1,0.1),function(x) E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=x)))
#1852.11 sec sur l'ordi de l'université
#3783.75 sec sur Kronos

plot(seq(0,1,0.1),as.numeric(res),cex.lab=0.85,ylab=expression(paste('E[',(x^{nu[T]}-1)^'+',']')),xlab=expression(paste('Proportion ',(nu[T]))),type='l')#ylab=expression(paste('E[',(x^{nu[T]}-1)^'+',']'))

#fee_c_s<-0.018   
#fee_c_f<-0.006448 
c_s18<-c(-0.3812930, -0.3840199, -0.3859154, -0.3835929, -0.3804927, -0.3777547, -0.3757349, -0.3744792, -0.3739521, -0.3740703, -0.3747655)

#fee_c_s<-0.01224  
#fee_c_f<-0.01224
c_s1224<-c(-0.4281220, -0.4236095, -0.4137572, -0.4035169, -0.3949490, -0.3882057, -0.3830944, -0.3793837, -0.3768577, -0.3753254, -0.3746159)

#fee_c_s<-0.0    
#fee_c_f<-0.02448 
c_s00<-c(-0.5000000, -0.4868001, -0.4615917, -0.4401641, -0.4228823, -0.4091114, -0.3982256, -0.3897042, -0.3831251, -0.3781514, -0.3745170)

#Aucun frais
c_S_c_f_aucun<-c(-0.3351600, -0.3257552, -0.3225388, -0.3200613, -0.3175005,-0.3153406, -0.3138833, -0.3132019, -0.3132425, -0.3139312,-0.3152040)

ggp_c_s18<-data.frame(funds_optimal=c_s18)
ggp_c_s1224<-data.frame(funds_optimal=c_s1224)
ggp_c_s00<-data.frame(funds_optimal=c_s00)
ggp_c_S_c_f_aucun<-data.frame(funds_optimal=c_S_c_f_aucun)

ggp_c_s18$fee<-"cs18"
ggp_c_s1224$fee<-"cs1224"
ggp_c_s00$fee<-"cs0"
ggp_c_S_c_f_aucun$fee<-"aucun"

ggp_cs_tot<-rbind(ggp_c_S_c_f_aucun,ggp_c_s00,ggp_c_s1224,ggp_c_s18)
ggp_cs_tot$prop<-seq(0,1,0.1)
ggp_cs_tot$fee <- factor(ggp_cs_tot$fee, levels = c("aucun", "cs18","cs1224", "cs0"))
#Combin_ptf_frais<-rbind(Ptf_final_2448_0,Ptf_final_1224_1224,Ptf_final_18_6448,Ptf_final_0_0)

# 1) -Graphique avec plot- #
plot(seq(0,1,0.1),c_s00,cex.lab=0.8,ylab=expression(paste('E[',(x[T]^{nu}-1)^'+',']')),xlab=expression(paste("Proportion constante investie dans l'actif risqué ",(nu[t]))),type='l')#ylab=expression(paste('E[',(x^{nu[T]}-1)^'+',']'))
lines(seq(0,1,0.1),c_s1224,col='coral3',lty=2)
lines(seq(0,1,0.1),c_s18,col='darkgoldenrod3',lty=3)
legend(0.6,-0.44, legend=c(expression(paste(c[s],"=0.000%  et  ",c[f],"=2.448%")),expression(paste(c[s],"=1.224%  et  ",c[f],"=1.224%")),expression(paste(c[s],"=1.800%  et  ",c[f],"=0.6448%"))),col=c('black','coral3','darkgoldenrod3'),lty=c(1,2,3))

# 2) -Grpahique avec ggplot2- #
ggplot(data=ggp_cs_tot,aes(x=ggp_cs_tot$prop,y=ggp_cs_tot$funds_optimal,group=fee,color=fee))+
  geom_line(aes(linetype=fee))+
  geom_point(aes(shape=fee))+
  scale_color_discrete(name="",
                       breaks=c("aucun","cs0","cs1224","cs18"),
                       labels=c(expression("Aucun frais"),expression(atop(paste(c[s],"=0.000%"),paste(c[f],"=2.448%"))),expression(atop(paste(c[s],"=1.224%"),paste(c[f],"=1.224%"))),expression(atop(paste(c[s],"=1.800%"),paste(c[f],"=0.648%")))))+
  scale_linetype_discrete(name="",
                       breaks=c("aucun","cs0","cs1224","cs18"),
                       labels=c(expression("Aucun frais"),expression(atop(paste(c[s],"=0.000%"),paste(c[f],"=2.448%"))),expression(atop(paste(c[s],"=1.224%"),paste(c[f],"=1.224%"))),expression(atop(paste(c[s],"=1.800%"),paste(c[f],"=0.648%")))))+
  scale_shape_discrete(name="",
                       breaks=c("aucun","cs0","cs1224","cs18"),
                       labels=c(expression("Aucun frais"),expression(atop(paste(c[s],"=0.000%"),paste(c[f],"=2.448%"))),expression(atop(paste(c[s],"=1.224%"),paste(c[f],"=1.224%"))),expression(atop(paste(c[s],"=1.800%"),paste(c[f],"=0.648%")))))+
  
  #scale_fill_manual(name="",values=c("#7CAE00","yellow","#F8766D","deepskyblue3"),labels=c( expression("Aucun frais"), expression(atop(paste(c[s],"=0.000%"),paste(c[f],"=2.448%"))),expression(atop(paste(c[s],"=1.224%"),paste(c[f],"=1.224%"))),expression(atop(paste(c[s],"=1.800%"),paste(c[f],"=0.648%")))))+
  #theme_classic()+
  #guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = 'bottom',legend.title = element_blank())+##C77CFF
labs(x=expression(paste(nu[t],"  pour tout t compris dans"," [0,T]")), y=expression(paste("E"^P,"[ U( (",F[T],"-G)"^"+","+G",") ]")))


#### -Graphiques: comparaisons, variation des frais- ###
## -Graphique 1: variation de c_f- ##
# (c_s=1.224%) #


#Ne pas toucher#funds_cf0648<-data.frame(funds_optimal=as.numeric(Simulations_fonds_distinct(pre2_S_tilde_t)))
#Ne pas toucher#funds_cf1224<-data.frame(funds_optimal=as.numeric(Simulations_fonds_distinct(pre2_S_tilde_t)))
#Ne pas toucher#funds_cf2448<-data.frame(funds_optimal=as.numeric(Simulations_fonds_distinct(pre2_S_tilde_t)))

funds_cf0648$comp<-'cf0648'
funds_cf1224$comp<-'cf1224'
funds_cf2448$comp<-'cf2448'

Compa_combin_funds_cf<-rbind(funds_cf0648,funds_cf1224,funds_cf2448)

ggplot(data=Compa_combin_funds_cf,aes(Compa_combin_funds_cf$funds_optimal,group=comp,fill=comp))+
  geom_histogram(colour='black',bins=100,alpha=0.5,position = "identity")+
  labs(x=expression(paste(F[T]^"*")), y="Réalisations")+
  #scale_x_continuous(breaks=c(x_concavi),labels=expression(paste(widehat(F),'(',1,')')))+
 theme_classic()+
scale_fill_manual(name="",values=c("yellow","#7CAE00",'#00BFC3'),labels=c( expression(paste(c[f],"=0.648%")),expression(paste(c[f],"=1.224%")),expression(paste(c[f],"=2.448%"))))+
theme(legend.position = 'bottom',legend.title = element_blank())  


## -Graphique 2: variation de c_s- ##
# (c_f=1.224%) #


#Ne pas toucher#funds_cs0<-data.frame(funds_optimal=as.numeric(Simulations_fonds_distinct(pre2_S_tilde_t)))
#Ne pas toucher#funds_cs1224<-data.frame(funds_optimal=as.numeric(Simulations_fonds_distinct(pre2_S_tilde_t)))
#Ne pas toucher#funds_cs18<-data.frame(funds_optimal=as.numeric(Simulations_fonds_distinct(pre2_S_tilde_t)))

funds_cs0$comp<-'cf0'
funds_cs1224$comp<-'cf1224'
funds_cs18$comp<-'cf18'

Compa_combin_funds_cs<-rbind(funds_cs0,funds_cs1224,funds_cs18)

ggplot(data=Compa_combin_funds_cs,aes(Compa_combin_funds_cs$funds_optimal,group=comp,fill=comp))+
  geom_histogram(colour='black',bins=100,alpha=0.5,position = "identity")+
  labs(x=expression(paste(F[T]^"*")), y="Réalisations")+
  theme_classic()+
  #scale_x_continuous(breaks=c(x_concavi),labels=expression(paste(widehat(F),'(',1,')')))+
  scale_fill_manual(name="",values=c('firebrick1','lightsalmon1','darkgoldenrod1'),labels=c( expression(paste(c[s],"=0.000%")),expression(paste(c[s],"=1.224%")),expression(paste(c[s],"=1.800%"))))+
  theme(legend.position = 'bottom',legend.title = element_blank())  


## -Graphique 3: variation de gamma- ##
# (c_f=1.224% et c_s=1.224% ) #

funds_gamma2<-data.frame(funds_optimal=as.numeric(Simulations_fonds_distinct(pre2_S_tilde_t)))
funds_gamma3<-data.frame(funds_optimal=as.numeric(Simulations_fonds_distinct(pre2_S_tilde_t)))
funds_gamma4<-data.frame(funds_optimal=as.numeric(Simulations_fonds_distinct(pre2_S_tilde_t)))
funds_gamma7<-data.frame(funds_optimal=as.numeric(Simulations_fonds_distinct(pre2_S_tilde_t)))

funds_gamma2$comp<-'gamma2'
funds_gamma3$comp<-'gamma3'
funds_gamma4$comp<-'gamma4'
funds_gamma7$comp<-'gamma7'

Compa_combin_funds_gamma<-rbind(funds_gamma2,funds_gamma3,funds_gamma4,funds_gamma7)

ggplot(data=Compa_combin_funds_gamma,aes(Compa_combin_funds_gamma$funds_optimal,group=comp,fill=comp))+
  geom_histogram(colour='black',bins=100,alpha=0.5,position = "identity")+
  labs(x=expression(paste("F",""[T])), y="Réalisations")+
  theme_classic()+
  scale_fill_manual(name="",values=c('gold2','darkorange1',"violet","red"),labels=c( expression(paste(gamma,"=2")),expression(paste(gamma,"=3")),expression(paste(gamma,"=4")),expression(paste(gamma,"=7"))))+
  theme(legend.position = 'bottom',legend.title = element_blank())  




## -Graphique 4: variation de la maturité du fonds distinct- ##
# (c_f=1.224% et c_s=1.224% ) #

#Fait# funds_maturite1<-data.frame(funds_optimal=as.numeric(Simulations_fonds_distinct(pre2_S_tilde_t)))
#Fait# funds_maturite5<-data.frame(funds_optimal=as.numeric(Simulations_fonds_distinct(pre2_S_tilde_t)))
#Fait# funds_maturite10<-data.frame(funds_optimal=as.numeric(Simulations_fonds_distinct(pre2_S_tilde_t)))
#Fait# funds_maturite15<-data.frame(funds_optimal=as.numeric(Simulations_fonds_distinct(pre2_S_tilde_t)))

funds_maturite1$comp<-'T=1'
funds_maturite5$comp<-'T=5'
funds_maturite10$comp<-'T=10'
funds_maturite15$comp<-'T=15'

Compa_combin_funds_maturite<-rbind(funds_maturite1,funds_maturite5,funds_maturite10,funds_maturite15)
Compa_combin_funds_maturite$comp<-factor(Compa_combin_funds_maturite$comp,levels=c('T=1','T=5','T=10','T=15'))


ggplot(data=Compa_combin_funds_maturite,aes(Compa_combin_funds_maturite$funds_optimal,group=comp,fill=comp,order = as.numeric(rev(comp))))+
  geom_histogram(colour='black',bins=150,alpha=0.4,position = "identity")+
  labs(x=expression(paste("F"[T]^"*")), y="Réalisations")+
  theme_classic()+
  #scale_fill_manual(name=""s,values=c('gold2','darkorange1',"violet","red"),labels=c( expression(paste(gamma,"=2")),expression(paste(gamma,"=3")),expression(paste(gamma,"=4")),expression(paste(gamma,"=7"))))+
  theme(legend.position = 'bottom',legend.title = element_blank())  

ggplot(data=Compa_combin_funds_maturite,aes(Compa_combin_funds_maturite$funds_optimal,group=comp,fill=comp,linetype=comp,color=comp,order = as.numeric(rev(comp))))+
  stat_ecdf(geom='step')+
  labs(x=expression(paste("F"[T]^"*")), y="Distribution cumulative")+
  theme(legend.position = 'bottom',legend.title = element_blank())  

## -Graphique 4.1: variation de la maturité du fonds distinct_ MMB [0,1]- ##
# (c_f=1.224% et c_s=1.224% ) #

#Fait# MMB_funds_maturite1<-data.frame(funds_optimal=as.numeric(E_utility_martingale_Borne(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t)))
#Fait# MMB_funds_maturite5<-data.frame(funds_optimal=as.numeric(E_utility_martingale_Borne(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t)))
#Fait# MMB_funds_maturite10<-data.frame(funds_optimal=as.numeric(E_utility_martingale_Borne(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t)))
#Fait# MMB_funds_maturite15<-data.frame(funds_optimal=as.numeric(E_utility_martingale_Borne(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t)))

MMB_funds_maturite1$comp<-'T=1'
MMB_funds_maturite5$comp<-'T=5'
MMB_funds_maturite10$comp<-'T=10'
MMB_funds_maturite15$comp<-'T=15'

MMB_Compa_combin_funds_maturite<-rbind(MMB_funds_maturite1,MMB_funds_maturite5,MMB_funds_maturite10,MMB_funds_maturite15)
MMB_Compa_combin_funds_maturite$comp<-factor(MMB_Compa_combin_funds_maturite$comp,levels=c('T=1','T=5','T=10','T=15'))


ggplot(data=MMB_Compa_combin_funds_maturite,aes(MMB_Compa_combin_funds_maturite$funds_optimal,group=comp,fill=comp,order = as.numeric(rev(comp))))+
  geom_histogram(colour='black',bins=150,alpha=0.4,position = "identity")+
  labs(x=expression(paste("F"[T]^"*")), y="Réalisations")+
  theme_classic()+
  #scale_fill_manual(name=""s,values=c('gold2','darkorange1',"violet","red"),labels=c( expression(paste(gamma,"=2")),expression(paste(gamma,"=3")),expression(paste(gamma,"=4")),expression(paste(gamma,"=7"))))+
  theme(legend.position = 'bottom',legend.title = element_blank())  

ggplot(data=MMB_Compa_combin_funds_maturite,aes(MMB_Compa_combin_funds_maturite$funds_optimal,group=comp,fill=comp,linetype=comp,color=comp,order = as.numeric(rev(comp))))+
  stat_ecdf(geom='step')+
  labs(x=expression(paste("F"[T])), y="Distribution cumulative")+
  theme(legend.position = 'bottom',legend.title = element_blank())  

## -Graphique 5: Pour différentes combinaisons de frais- ##
ajust_graphique5<-2.5
library(scales)

#Ne pas toucher#ptf_2448_0<-Simulations_fonds_distinct(pre2_S_tilde_t)
#Ne pas toucher#ptf_1224_1224<-Simulations_fonds_distinct(pre2_S_tilde_t)
#Ne pas toucher#ptf_18_648<-Simulations_fonds_distinct(pre2_S_tilde_t)
#Ne pas toucher#ptf_0_0<-Simulations_fonds_distinct(pre2_S_tilde_t)
  
Ptf_final_2448_0<-data.frame(valeur_opt=as.numeric(lapply(ptf_2448_0,function(x) min(x,ajust_graphique5))))#Obtenu à l'aide de la méthode 4 MemoireCodes3.3
Ptf_final_1224_1224<-data.frame(valeur_opt=as.numeric(lapply(ptf_1224_1224,function(x)min(x,ajust_graphique5))))#Obtenu à l'aide de la méthode 4 MemoireCodes3.3
Ptf_final_18_648<-data.frame(valeur_opt=as.numeric(lapply(ptf_18_6448,function(x)min(x,ajust_graphique5))))#Obtenu à l'aide de la méthode 4 MemoireCodes3.3
Ptf_final_0_0<-data.frame(valeur_opt=as.numeric(lapply(ptf_0_0,function(x)min(x,ajust_graphique5))))#Obtenu à l'aide de la méthode 4 MemoireCodes3.3

Ptf_final_2448_0$comp<-'2.c_s=0.000% et c_f=2.448%'
Ptf_final_1224_1224$comp<-'3.c_s=1.224% et c_f=1.224%'
Ptf_final_18_648$comp<-'4. c_s=1.800% et c_f=0.648%'
Ptf_final_0_0$comp<-'1. Aucun frais'


Combin_ptf_frais<-rbind(Ptf_final_2448_0,Ptf_final_1224_1224,Ptf_final_18_6448,Ptf_final_0_0)
Combin_ptf_frais$comp <- factor(Combin_ptf_frais$comp, levels = c('1. Aucun frais', '4. c_s=1.800% et c_f=0.648%','3.c_s=1.224% et c_f=1.224%', '2.c_s=0.000% et c_f=2.448%'))


ggplot(data=Combin_ptf_frais,aes(Combin_ptf_frais$valeur_opt,group=comp,fill=comp))+
  geom_histogram(colour='black',binwidth = 0.01,alpha=0.5,position = "identity")+
  labs(x=expression(paste(F[T]^"*")), y="Réalisations")+
  #scale_x_continuous(breaks=c(x_concavi),labels=expression(paste(widehat(x),'(',D[T],')')))+
  scale_fill_manual(name="",values=c("#F8766D","#C77CFF","#00BFC4","#7CAE00"),labels=c( expression("Aucun frais"), expression(atop(paste(c[s],"=0.000%"),paste(c[f],"=2.448%"))),expression(atop(paste(c[s],"=1.224%"),paste(c[f],"=1.224%"))),expression(atop(paste(c[s],"=1.800%"),paste(c[f],"=0.648%")))))+
  theme_classic()+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = 'bottom',legend.title = element_blank())##C77CFF






#### Graphique 6: Valeur présente de la garantie pour diff. combinaison de frais #####
## Pour différentes stratégies, non optimales 

# 1) Changer les formules 1-2-3-4 pour obtenir la valeur présente
# 2) Rouler le code pour les frais c_f=2.448-c_s:

#a) Pour c_s=0.000% et c_f=2.448%#fait
#b) Pour c_s=0.25% et c_f=2.448-0.5%#fait
#c) Pour c_s=0.5% et c_f=2.448-0.5% #fait
#d) Pour c_s=0.75% et c_f=2.448-0.75% #fait
#d) Pour c_s=1.0% et c_f=2.448-0.5% #fait
#e) Pour c_s=1.25% et c_f=2.448-1.25%
#f) Pour c_s=1.5% et c_f=2.448-0.5%#fait


timer5<-proc.time()
ggp_v_ac_MMB_c_s175<-E_utility_martingale_Borne(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t)
proc.time()-timer5
ggp_v_ac_MM_c_s175<-E_utility_martingale(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t)
proc.time()-timer5
ggp_v_ac_100_c_s175<-E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=1.0)
proc.time()-timer5
ggp_v_ac_60_c_s175<-E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=0.6)
proc.time()-timer5
ggp_v_ac_40_c_s175<-E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=0.4)
proc.time()-timer5
ggp_v_ac_20_c_s175<-E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=0.2)
proc.time()-timer5
ggp_v_ac_Merton_c_s175<-E_utility_prop_cte(pre2_S_tilde_t,pre2_xi_tilde_t,B_tilde_t,prop_act_r=cte_Merton)
proc.time()-timer5
ggp_v_ac_terminal_c_s175<-Simulations_fonds_distinct(pre2_S_tilde_t)
proc.time()-timer5 #3529.83 sec avec Kronos






# -Formatage des données- #
ggpVac_tot_c_s0<-c(0.9881922,1.0809471,1.0012364,0.9231385,0.8835673,0.8445091,0.8382593,
                   1.1861791)#c(ggp_v_ac_MMB_c_s0,ggp_v_ac_MM_c_s0,ggp_v_ac_100_c_s0,ggp_v_ac_60_c_s0,ggp_v_ac_40_c_s0,ggp_v_ac_20_c_s0,ggp_v_ac_Merton_c_s0,ggp_v_ac_terminal_c_s0)
ggpVac_tot_c_s025<-c(0.9906008,1.0775105,1.0008413,0.9270760,0.8892291, 0.8510618, 0.8407773,
                     1.1796944)#c(ggp_v_ac_MMB_c_s025,ggp_v_ac_MM_c_s025,ggp_v_ac_100_c_s025,ggp_v_ac_60_c_s025,ggp_v_ac_40_c_s025,ggp_v_ac_20_c_s025,ggp_v_ac_Merton_c_s025,ggp_v_ac_terminal_c_s025)
ggpVac_tot_c_s05<-c(0.9935308, 1.0748678, 1.0000568, 0.9309548, 0.8950383, 0.8583977, 0.8445159,
                   1.1739887)#c(ggp_v_ac_MMB_c_s05,ggp_v_ac_MM_c_s05,ggp_v_ac_100_c_s05,ggp_v_ac_60_c_s05,ggp_v_ac_40_c_s05,ggp_v_ac_20_c_s05,ggp_v_ac_Merton_c_s05,ggp_v_ac_terminal_c_s05)
ggpVac_tot_c_s075<-c(0.9959830, 1.0715459, 1.0003881, 0.9352825, 0.9013353, 0.8669381, 0.8510134,
                     1.1683810)#c(ggp_v_ac_MMB_c_s075,ggp_v_ac_MM_c_s075,ggp_v_ac_100_c_s075,ggp_v_ac_60_c_s075,ggp_v_ac_40_c_s075,ggp_v_ac_20_c_s075,ggp_v_ac_Merton_c_s075,ggp_v_ac_terminal_c_s075)
ggpVac_tot_c_s1<-c(0.9987706,1.0693147, 0.9988078, 0.9387356, 0.9072777, 0.8760920,0.8615410,
                   1.1628909)#c(ggp_v_ac_MMB_c_s1,ggp_v_ac_MM_c_s1,ggp_v_ac_100_c_s1,ggp_v_ac_60_c_s1,ggp_v_ac_40_c_s1,ggp_v_ac_20_c_s1,ggp_v_ac_Merton_c_s1,ggp_v_ac_terminal_c_s1)
ggpVac_tot_c_s125<-c(1.0022422,1.0674588,0.9997839,0.9440816,0.9152262,0.8880798,0.8803042,
                    1.1591854)#c(ggp_v_ac_MMB_c_s125,ggp_v_ac_MM_c_s125,ggp_v_ac_100_c_s125,ggp_v_ac_60_c_s125,ggp_v_ac_40_c_s125,ggp_v_ac_20_c_s125,ggp_v_ac_Merton_c_s125,ggp_v_ac_terminal_c_s125)
ggpVac_tot_c_s15<-c(1.0068096, 1.0675851, 0.9993189, 0.9491258, 0.9234596, 0.9013734, 0.9037023,
                   1.1548775)#c(ggp_v_ac_MMB_c_s15,ggp_v_ac_MM_c_s15,ggp_v_ac_100_c_s15,ggp_v_ac_60_c_s15,ggp_v_ac_40_c_s15,ggp_v_ac_20_c_s15,ggp_v_ac_Merton_c_s15,ggp_v_ac_terminal_c_s15)
ggpVac_tot_c_s175<-c(1.0118382,1.0674703, 1.0008413, 0.9548784, 0.9320741, 0.9152954, 0.9300795,
                   1.1531930)#c(ggp_v_ac_MMB_c_s175,ggp_v_ac_MM_c_s175,ggp_v_ac_100_c_s175,ggp_v_ac_60_c_s175,ggp_v_ac_40_c_s175,ggp_v_ac_20_c_s175,ggp_v_ac_Merton_c_s175,ggp_v_ac_terminal_c_s175)

V_ac_c_s0<-data.frame(value=ggpVac_tot_c_s0)
V_ac_c_s025<-data.frame(value=ggpVac_tot_c_s025)
V_ac_c_s05<-data.frame(value=ggpVac_tot_c_s05)
V_ac_c_s075<-data.frame(value=ggpVac_tot_c_s075)
V_ac_c_s1<-data.frame(value=ggpVac_tot_c_s1)
V_ac_c_s125<-data.frame(value=ggpVac_tot_c_s125)
V_ac_c_s15<-data.frame(value=ggpVac_tot_c_s15)
V_ac_c_s175<-data.frame(value=ggpVac_tot_c_s175)

V_ac_c_s0$c_S<-c(0)
V_ac_c_s025$c_S<-c(0.25)
V_ac_c_s05$c_S<-c(0.5)
V_ac_c_s075$c_S<-c(0.75)
V_ac_c_s1$c_S<-c(1.0)
V_ac_c_s125$c_S<-c(1.25)
V_ac_c_s15$c_S<-c(1.5)
V_ac_c_s175$c_S<-c(1.75)

V_ac_tot<-rbind(V_ac_c_s0,V_ac_c_s025,V_ac_c_s05,V_ac_c_s075,V_ac_c_s1,V_ac_c_s125,V_ac_c_s15,V_ac_c_s175)
V_ac_tot$method<-c('6-Optimale [0,1]','7-Optimale [0,2]','1-100% risque','2-60% risque','3-40% risque','4-20% risque','5-Constante Merton','8-Optimale')
# La ligne d'en haut à vérifier, il faut que la méthode se répète.

ggplot(data=V_ac_tot,aes(x=V_ac_tot$c_S,y=V_ac_tot$value,group=method,color=method))+
  geom_line(aes(linetype=method))+
  geom_point(aes(shape=method))+ 
  theme(legend.position = "bottom",legend.title =element_blank(),legend.text = element_text(size = 5))+
  labs(x=expression(paste(c[s],' (%)')), y=expression(paste("E"^P,"[",xi[T]," ( (","F"[T],"-G)"^"+","+G )"," ]")))


#ggplot_build(g)$data #pour obtenir les couluers,points, etc.


#### Graphique 7: Valeur présente de la garantie pour diff. gamma#####
## Pour la stratégie optimale ##


#a) Pour c_s=0.000% et c_f=2.448 #fait
#b) Pour c_s=0.0025% et c_f=2.448-0.0025% #fait
#c) Pour c_s=0.005% et c_f=2.448-0.005%  #fait
#d) Pour c_s=0.0075% et c_f=2.448-0.0075% #fait
#e) Pour c_s=0.01% et c_f=2.448-0.01%  #fait
#f) Pour c_s=0.0125% et c_f=2.448-0.0125%#fait
#g) Pour c_s=0.015% et c_f=2.448-0.015% #fait
#h) Pour c_s=0.0175% et c_f=2.448-0.0175% 

#-gamma=2-
ggp_v_ac_g2_c_s175<-Simulations_fonds_distinct(pre2_S_tilde_t)

#-gamma=3-
ggp_v_ac_g3_c_s175<-Simulations_fonds_distinct(pre2_S_tilde_t)

#-gamma=4-
ggp_v_ac_g4_c_s175<-Simulations_fonds_distinct(pre2_S_tilde_t)

#-gamma=5-
ggp_v_ac_g5_c_s175<-Simulations_fonds_distinct(pre2_S_tilde_t)

#-gamma=6-
ggp_v_ac_g6_c_s175<-Simulations_fonds_distinct(pre2_S_tilde_t)

#-gamma=7-
ggp_v_ac_g7_c_s175<-Simulations_fonds_distinct(pre2_S_tilde_t)




# -Formatage des données- #
ggpVac_ga_tot_c_s0<-c(1.2553966,1.1858717,1.1418219,1.1073739,1.0814788,1.060112)#c(ggp_v_ac_g2_c_s0,ggp_v_ac_g3_c_s0,ggp_v_ac_g4_c_s0,ggp_v_ac_g5_c_s0,ggp_v_ac_g6_c_s0,ggp_v_ac_g7_c_s0)
ggpVac_ga_tot_c_s025<-c( 1.243003,1.179959,1.136115,1.103264,1.078331,1.057911)#c(ggp_v_ac_g2_c_s025,ggp_v_ac_g3_c_s025,ggp_v_ac_g4_c_s025,ggp_v_ac_g5_c_s025,ggp_v_ac_g6_c_s025,ggp_v_ac_g7_c_s025)
ggpVac_ga_tot_c_s05<-c( 1.234192,1.174178,1.131885,1.100623,1.076360,1.056395)#c(ggp_v_ac_g2_c_s05,ggp_v_ac_g3_c_s05,ggp_v_ac_g4_c_s05,ggp_v_ac_g5_c_s05,ggp_v_ac_g6_c_s05,ggp_v_ac_g7_c_s05)
ggpVac_ga_tot_c_s075<-c( 1.223556,1.167625,1.128256,1.098047,1.074802,1.055590)#c(ggp_v_ac_g2_c_s075,ggp_v_ac_g3_c_s075,ggp_v_ac_g4_c_s075,ggp_v_ac_g5_c_s075,ggp_v_ac_g6_c_s075,ggp_v_ac_g7_c_s075)
ggpVac_ga_tot_c_s1<-c(1.215882,1.163435,1.124941,1.095723,1.073407,1.055695)#c(ggp_v_ac_g2_c_s1,ggp_v_ac_g3_c_s1,ggp_v_ac_g4_c_s1,ggp_v_ac_g5_c_s1,ggp_v_ac_g6_c_s1,ggp_v_ac_g7_c_s1)
ggpVac_ga_tot_c_s125<-c(1.209324,1.158076,1.121104,1.093232,1.070939,1.053376)#c(ggp_v_ac_g2_c_s1255,ggp_v_ac_g3_c_s125,ggp_v_ac_g4_c_s125,ggp_v_ac_g5_c_s125,ggp_v_ac_g6_c_s125,ggp_v_ac_g7_c_s125)
ggpVac_ga_tot_c_s15<-c( 1.203194,1.154905,1.120078,1.093828,1.071933,1.051930)#c(ggp_v_ac_g2_c_s15,ggp_v_ac_g3_c_s15,ggp_v_ac_g4_c_s15,ggp_v_ac_g5_c_s15,ggp_v_ac_g6_c_s15,ggp_v_ac_g7_c_s15)
ggpVac_ga_tot_c_s175<-c(1.196787,1.152575,1.118419,1.092845,1.072702,1.055441)#c(ggp_v_ac_g2_c_s175,ggp_v_ac_g3_c_s175,ggp_v_ac_g4_c_s175,ggp_v_ac_g5_c_s175,ggp_v_ac_g6_c_s175,ggp_v_ac_g7_c_s175)

  
V_ac_ga_c_s0<-data.frame(value=ggpVac_ga_tot_c_s0)
V_ac_ga_c_s025<-data.frame(value=ggpVac_ga_tot_c_s025)
V_ac_ga_c_s05<-data.frame(value=ggpVac_ga_tot_c_s05)
V_ac_ga_c_s075<-data.frame(value=ggpVac_ga_tot_c_s075)
V_ac_ga_c_s1<-data.frame(value=ggpVac_ga_tot_c_s1)
V_ac_ga_c_s125<-data.frame(value=ggpVac_ga_tot_c_s125)
V_ac_ga_c_s15<-data.frame(value=ggpVac_ga_tot_c_s15)
V_ac_ga_c_s175<-data.frame(value=ggpVac_ga_tot_c_s175)

V_ac_ga_c_s0$c_S<-c(0)
V_ac_ga_c_s025$c_S<-c(0.25)
V_ac_ga_c_s05$c_S<-c(0.5)
V_ac_ga_c_s075$c_S<-c(0.75)
V_ac_ga_c_s1$c_S<-c(1.0)
V_ac_ga_c_s125$c_S<-c(1.25)
V_ac_ga_c_s15$c_S<-c(1.5)
V_ac_ga_c_s175$c_S<-c(1.75)

V_ac_ga_tot<-rbind(V_ac_ga_c_s0,V_ac_ga_c_s025,V_ac_ga_c_s05,V_ac_ga_c_s075,V_ac_ga_c_s1,V_ac_ga_c_s125,V_ac_ga_c_s15,V_ac_ga_c_s175)
V_ac_ga_tot$gam<-c("g2","g3","g4","g5","g6","g7")
# La ligne d'en haut à vérifier, il faut que la méthode se répète.

ggplot(data=V_ac_ga_tot,aes(x=V_ac_ga_tot$c_S,y=V_ac_ga_tot$value,group=gam,color=gam))+
  geom_line(aes(linetype=gam))+
  geom_point(aes(shape=gam))+ 
  scale_shape_discrete(name="",
                       breaks=c("g2", "g3", "g4","g5","g6","g7"),
                       labels=c(expression(paste(gamma,"=2")), expression(paste(gamma,"=3")), expression(paste(gamma,"=4")),expression(paste(gamma,"=5")),expression(paste(gamma,"=6")),expression(paste(gamma,"=7"))))+
  scale_color_discrete(name="",
                       breaks=c("g2", "g3", "g4","g5","g6","g7"),
                       labels=c(expression(paste(gamma,"=2")), expression(paste(gamma,"=3")), expression(paste(gamma,"=4")),expression(paste(gamma,"=5")),expression(paste(gamma,"=6")),expression(paste(gamma,"=7"))))+
  scale_linetype_discrete(name="",
                          breaks=c("g2", "g3", "g4","g5","g6","g7"),
                          labels=c(expression(paste(gamma,"=2")), expression(paste(gamma,"=3")), expression(paste(gamma,"=4")),expression(paste(gamma,"=5")),expression(paste(gamma,"=6")),expression(paste(gamma,"=7"))))+
  theme(legend.position = "bottom",legend.text = element_text(size = 10) )+
  labs(x=expression(paste(c[s],' (%)')), y=expression(paste("E"^P,"[",xi[T]," ( (","F"[T]^"*","-G)"^"+","+G )"," ]")))

