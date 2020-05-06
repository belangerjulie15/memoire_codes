############ MÉMOIRE: Maitrise en Mathématiques Finacières et Actuarielle  #################

# -SECTION 4: COMPARAISONS ET PARRALLÉLISATION: Power Utility - Simulation de marché avec frais #
#-Note: toutes les formules sont UNIQUEMENT applicables à la fonction de Power Utility.
#Ce dossier fait exactement les mêmes calculs que Memoire_Codes3, mais où tous les 
#codes ont été parrallélisés pour augmenter la vitesse d'exécution.


########## PARAMÈTRES ##########################################
library(ggplot2)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)



Maturi<-10         #Time until maturity
r_no_risk<-0.02    #Risk free rate 
alpha<-0.04        #Risky rate
sigma<-0.2         #Volatility
gamma<-2           #Parameter of Utility function
S_0<-1             #Initial value of the asset (S_0>0)
B_0<-1             #Initial value of the bank account
budget<-1          #Initial Budget amount
N_Simulations<-1000 #Number of Simulations
fee_c_s<-0.005     #Fee applied of the risky asset
fee_c_f<-0.005     #Fee applied of the funds 
Frequ<-52          #Frequency of rebalancing the portfolio

a_call_sim<-1      #Multiplicator of the variable annuity
b_call_sim<-1      #Strike price of the variable annuity
K_call_sim<-1      #Constant of the variable annuity

alpha_tilde<-alpha-fee_c_s-fee_c_f
r_no_risk_tilde<-r_no_risk-fee_c_f
theta_sim<-(alpha-r_no_risk)/sigma
theta_sim_tilde<-(alpha_tilde-r_no_risk_tilde)/sigma

cte_Merton<-theta_sim_tilde/(sigma*gamma)

numCores_paral <- detectCores() -1
################################################################


######### Directions ############################################
#  Il faut rouler toute la section FONCTIONS pour activer les fonctions 
#qu'on veut utiliser.

#Il y a plusieurs sections pour différentes proportions (cte ou dynamiques) investies 
#dans l'actif risqué.
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

lambda_opt<-lambda_optimal(r_FL=r_no_risk,alpha_FL=alpha,sigma_FL=sigma,gamma_FL=gamma,K_FL=K_call_sim,b_FL=b_call_sim,a_FL=a_call_sim,X_0_FL=budget,T_FL=Maturi)
x_concavi<-Find_x_theta_PU(b_call_sim,a_call_sim,K_call_sim,gamma)
y_concavi<-x_concavi^(-gamma)

random_var<-matrix(rnorm(Maturi*Frequ*N_Simulations),N_Simulations,Maturi*Frequ) #Matrice actif sans risque
pre1_S_tilde_t<-S_0*exp((alpha_tilde-0.5*sigma^2)/Frequ+sigma*sqrt(1/Frequ)*random_var)
pre2_S_tilde_t<-cbind((rep(S_0,N_Simulations)),pre1_S_tilde_t)

pre1_B_tilde_t<-rep(exp(r_no_risk_tilde/Frequ),Frequ*Maturi) #vecteur actif sans risque
B_tilde_t<-B_0*cumprod(c(B_0,pre1_B_tilde_t))
#########################################################################


########### SECTION 1-) "Proportion optimale", avec la méthode martingale ###########
############################# PARRALÉLISATION ###################################### 

E_utility_martingale_paral<-function(matrice_pre2_S,vecteur_B){
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
  
  
  #Loop sur les simulations parrallélisées
  loop_parralelized<-function(n){
    matrice_S[n,]<-cumprod(matrice_pre2_S[n,])
    
    for(m in 1:(Frequ*Maturi)){ #Frequ*Maturi
      petit_t<-(m-1)/(Frequ*Maturi)
      xi_tilde[n,m]<-matrice_S[n,m]^(-theta_sim_tilde/sigma)*exp((theta_sim_tilde*alpha_tilde/sigma-theta_sim_tilde*sigma/2-r_no_risk_tilde-0.5*theta_sim_tilde^2)*petit_t)
      
      interior[n,m]<-((r_no_risk_tilde+(1/gamma-0.5)*theta_sim_tilde^2)*(Maturi-petit_t)+log(y_concavi/(xi_tilde[n,m]*lambda_opt)))/(theta_sim_tilde*sqrt(Maturi-petit_t))
      nu_Act_R[n,m]<-(theta_sim_tilde/gamma+dnorm(interior[n,m])/(pnorm(interior[n,m])*sqrt(Maturi-petit_t)))/sigma
      
      if(m==1){
        unit_risque[n,m]<-nu_Act_R[n,m]*budget/matrice_S[n,m]
        unit_n_risque[n,m]<-(budget-unit_risque[n,m]*matrice_S[n,m])/B_tilde_t[m]
        
        processus_ptf[n,m]<-budget}
      else{
        unit_invest<-Investment_Fonct(matrice_S[n,m],B_tilde_t[m],unit_risque[n,(m-1)], unit_n_risque[n,(m-1)],nu_Act_R[n,m])
        unit_risque[n,m]<-unit_invest[1]
        unit_n_risque[n,m]<-unit_invest[2]
        processus_ptf[n,m]<-unit_risque[n,m]*matrice_S[n,m]+unit_n_risque[n,m]*B_tilde_t[m]}
      
    }
    processus_ptf[n,(Frequ*Maturi+1)]<-unit_risque[n,(Frequ*Maturi)]*matrice_S[n,(Frequ*Maturi+1)]+unit_n_risque[n,(Frequ*Maturi)]*B_tilde_t[(Frequ*Maturi+1)]
    pre_funds_d<-a_call_sim*max(0,processus_ptf[n,(Frequ*Maturi+1)]-b_call_sim)+K_call_sim
    return(pre_funds_d)
  }
  
  
  #registerDoParallel(cores=numCores_paral)
  #pre1_funds_d<-foreach(i=1:N_Simulations) %dopar% loop_parralelized(i)
  pre1_funds_d<-mclapply(1:N_Simulations,loop_parralelized, mc.cores = numCores_paral )#(pour MAC)
  funds_d<-as.numeric(pre1_funds_d)
  
  #xi_tilde[n,(Frequ*Maturi+1)]<-matrice_S[n,(Frequ*Maturi+1)]^(-theta_sim_tilde/sigma)*exp((theta_sim_tilde*alpha_tilde/sigma-theta_sim_tilde*sigma/2-r_no_risk_tilde-0.5*theta_sim_tilde^2)*petit_t)
  #CB<-mean(xi_tilde[n,(Frequ*Maturi+1)]*processus_ptf[n,(Frequ*Maturi+1)])
  
  Uty<-as.numeric(lapply(funds_d,function(x)P_Utility(X_Tu=x,gamma=gamma)))
  EU<-mean(Uty)
  
  return(EU)
}

timer_parrall<-proc.time()
E_utility_martingale_paral(pre2_S_tilde_t,B_tilde_t)
proc.time()-timer_parrall


########################## TEST ##########################






cum_sum <- function(p) {
  retur <- 0
  for (i in 1:p) {
    retur <- retur+ i
  }
  return(retur)
}

n <- 20000
system.time(lapply(1:n, cum_sum))
#system.time( mclapply(n,cum_sum, mc.cores = 6))


registerDoParallel(cores=39)
system.time(as.numeric(foreach(i=1:n) %dopar% cum_sum(i)))

cl <- makeCluster(detectCores())
clusterEvalQ(cl, library(lme4))


cl <- makeForkCluster(detectCores()-1)
clusterSetRNGStream(cl, 1001)

system.time(save3 <- parLapply(cl, 1:n, cum_sum ))
stopCluster(cl)


