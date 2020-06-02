############ MÉMOIRE: Maitrise en Mathématiques Finacières et Actuarielle  #################

# -SECTION 6: Power Utility - Maximisation de l'utilité pour un groupe d'assurés#


######### Directions ############################################
#-Note: toutes les formules sont UNIQUEMENT applicables à la fonction de Power Utility.

#  Il faut installer/ouvrir le package Optimisation.Power.Utility pour activer les fonctions
#qu'on veut utiliser.

#Il y a plusieurs sections pour différentes proportions (cte ou dynamiques) investies
#dans l'actif risqué.
################################################################

########## PACKAGES ##########################################
library(ggplot2)
#library(Optimisation.Power.Utility)
################################################################

########## PARAMÈTRES ##########################################
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
################################################################

########### Fonctions du packages Optimisation.Power.Utility ###########

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

PowerUtility_2<-function(gamma1,gamma2,x){
  sumUtl<-P_Utility(x,gamma1)+P_Utility(x,gamma2)
  return(sumUtl)
}

PowerUtility_3<-function(gamma1,gamma2,gamma3,x){
  sumUtl<-P_Utility(x,gamma1)+P_Utility(x,gamma2)+P_Utility(x,gamma3)
  return(sumUtl)
}

PowerUtility_con_2<-function(gamma1,gamma2,b,a,K,x){
  sumUtl<-P_Utility_con(x, gamma1, b, a, K)+P_Utility_con(x, gamma2, b, a, K)
  return(sumUtl)
}

PowerUtility_con_3<-function(gamma1,gamma2,gamma3,b,a,K,x){
  sumUtl<-P_Utility_con(x, gamma1, b, a, K)+P_Utility_con(x, gamma2, b, a, K)++P_Utility_con(x, gamma3, b, a, K)
  return(sumUtl)
}

Find_x_theta_PU_Numerous<-function(b_o,a_o,K_o,vector_gamma){

  f <- function(x)  sum(P_Utility(payoff_call(x,b_o,a_o,K_o),vector_gamma))-sum(P_Utility(K_o,vector_gamma))-x*sum((x)^(-vector_gamma))
  sln_x<-uniroot(f, c(b_o,b_o+100))$root

  return(sln_x)
}#Verification (dans fichier function, Optimisation.Power.Utility) faite.


#####################################################


########### SECTION 0-) SIMULATIONS CONSTANTES POUR COMPARAISON ###########

lambda_opt<-lambda_optimal(r_FL=r_no_risk,alpha_FL=alpha,sigma_FL=sigma,gamma_FL=gamma,K_FL=K_call_sim,b_FL=b_call_sim,a_FL=a_call_sim,X_0_FL=budget,T_FL=Maturi)
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


######### Somme d'utilités: Comparaison individuel VS en groupe ############

test<-seq(0.7,2.5,0.05)
test_call<-as.numeric(lapply(test,function(x)payoff_call(x,b=1,a=1,K=1)))

plot(test,lapply(test,function(x) PowerUtility_2(2,3,x)),type='l',ylim = c(-2,0),ylab='',col='blue')
lines(test,lapply(test,function(x) P_Utility(x,2)))
lines(test,lapply(test,function(x) P_Utility(x,3)))

plot(test,lapply(test_call,function(x) PowerUtility_2(2,3,x)),type='l',ylim = c(-2,0),ylab='',col='blue')
lines(test,lapply(test_call,function(x) P_Utility(x,2)))
lines(test,lapply(test_call,function(x) P_Utility(x,3)))

plot(test,lapply(test_call,function(x) PowerUtility_3(2,3,5,x)),xaxt = "n",yaxt = "n",type='l',ylim = c(-2,0),ylab = "Utilité",xlab=expression(paste("Valeur de l'actif à maturité (x)" )),col='blue',xlim = c(0.8,2.7))
lines(test,lapply(test_call,function(x) P_Utility(x,2)),lty=2)
lines(test,lapply(test_call,function(x) P_Utility(x,3)),lty=4)
lines(test,lapply(test_call,function(x) P_Utility(x,5)),lty=5)
xtick<-c(b=b_call_sim,round(Find_x_theta_PU(b_o=b_call_sim,a_o=a_call_sim,K_o=K_call_sim,gamma_o=3),4),round(Find_x_theta_PU(b_o=b_call_sim,a_o=a_call_sim,K_o=K_call_sim,gamma_o=2),4),round(Find_x_theta_PU(b_o=b_call_sim,a_o=a_call_sim,K_o=K_call_sim,gamma_o=5),4),round(Find_x_theta_PU_Numerous(b_call_sim,a_call_sim,K_call_sim,c(2,3,5)),4))
axis(1, at=xtick, labels=c(expression(D[T]),expression(paste(tilde(x),"(",gamma[3],")")),expression(paste(tilde(x),"(",gamma[2],")")),expression(paste(tilde(x),"(",gamma[5],")")),''))
text(2.65, -0.38, expression(paste(gamma,'=2')))
text(2.65, -0.1, expression(paste(gamma,'=3')))
text(2.65, -0.001, expression(paste(gamma,'=5')))
text(2.65, -0.49, 'Totale')
abline(v=Find_x_theta_PU_Numerous(b_call_sim,a_call_sim,K_call_sim,c(2,3,5)),  col ='blue',lty=3)




##### Test: Maximisation Lagrangien pour S_{T} #######
#0.9575921 lambda optimal


xi_wrt_s_t<-function(S_repo,al_repo,r_repo,sigma_repo,T_repo){
  theta_repo<-(al_repo-r_repo)/sigma_repo
  
  xi_reponse<-S_repo*exp((theta_repo*al_repo/sigma_repo-theta_repo*sigma_repo/2-r_repo-0.5*theta_repo^2)*T_repo)
  return(xi_reponse)
}

lagrangien<-function(S_test,lambda_test){
  lagr<-P_Utility(S_test,3)-lambda_test*(xi_wrt_s_t(S_test,0.04,0.02,0.2,1)*S_test-1)
  return(lagr)
}
axe_x_test_00<-seq(0.1,2,0.1)


plot(axe_x_test_00,lapply(axe_x_test_00,function(x)P_Utility(x,3)),type='l')
lines(axe_x_test_00,lapply(axe_x_test_00,function(x)lagrangien(x,2)),col='red',lty=2)
lines(axe_x_test_00,lapply(axe_x_test_00,function(x)lagrangien(x,1)),col='blue',lty=2)
lines(axe_x_test_00,lapply(axe_x_test_00,function(x)lagrangien(x,3)),col='blue',lty=2)

l02<-max(as.numeric(lapply(axe_x_test_00,function(x)lagrangien(x,0.2))))
l03<-max(as.numeric(lapply(axe_x_test_00,function(x)lagrangien(x,0.3))))
l05<-max(as.numeric(lapply(axe_x_test_00,function(x)lagrangien(x,0.5))))
l075<-max(as.numeric(lapply(axe_x_test_00,function(x)lagrangien(x,0.75))))
l09<-max(as.numeric(lapply(axe_x_test_00,function(x)lagrangien(x,0.9))))
l1<-max(as.numeric(lapply(axe_x_test_00,function(x)lagrangien(x,1))))
l15<-max(as.numeric(lapply(axe_x_test_00,function(x)lagrangien(x,1.5))))
l2<-max(as.numeric(lapply(axe_x_test_00,function(x)lagrangien(x,2))))

plot(c(0.2,0.3,0.5,0.75,0.9,1,1.5,2),c(l02,l03,l05,l075,l09,l1,l15,l2))
