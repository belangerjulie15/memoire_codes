############ MÉMOIRE: Maitrise en Mathématiques Finacières et Actuarielle  #################

# -SECTION 3: Power Utility - Simulation de marché #
#-Note: toutes les formules sont UNIQUEMENT applicables à la fonction de Power Utility

########## PARAMÈTRES ##########################################
Maturi<-10         #Time until maturity
r_no_risk<-0.02    #Risk free rate
alpha<-0.04        #Risky rate
sigma<-0.2         #Volatility
gamma<-2           #Parameter of Utility function
S_0<-1             #Initial value of the asset (S_0>0)
budget<-1          #Initial Budget amount
N_Simulations<-100000 #Number of Simulations
fee_c_s<-0.00     #Fee applied of the risky asset
fee_c_f<-0.00     #Fee applied of the funds 

Frequ<-52          #Frequency of rebalancing the portfolio
  
a_call_sim<-1
b_call_sim<-1
K_call_sim<-1

library("ggplot2")
alpha_tilde<-alpha-fee_c_s-fee_c_f
r_no_risk_tilde<-r_no_risk-fee_c_f
  
theta_sim<-(alpha-r_no_risk)/sigma
theta_sim_tilde<-(alpha_tilde-r_no_risk_tilde)/sigma
################################################################




########## Graph: xi_t with respect to S_t ##########
ttt<-3
sst<-seq(0.1,3,0.01)
xxi<-lapply(sst,function(x)x^(-theta_sim/sigma)*exp((theta_sim*alpha/sigma-theta_sim*sigma/2-r_no_risk-0.5*theta_sim^2)*ttt))
plot(sst,xxi, type='l',xlab=expression(S[t]),ylab=expression(xi[t]))




########## Graph: xi_t with respect to S_t with fees ##########
ttt<-3
sst<-seq(0.1,3,0.01)
xxi<-lapply(sst,function(x)x^(-theta_sim_tilde/sigma)*exp((theta_sim_tilde*alpha_tilde/sigma-theta_sim_tilde*sigma/2-r_no_risk_tilde-0.5*theta_sim_tilde^2)*ttt))
plot(sst,xxi, type='l',xlab=expression(tilde(S)[t]),ylab=expression(tilde(xi)[t]))
#,yaxt = "n" ,xaxt = "n"

#### "Frequency" - SIMULATIONS - Market ########
# Power Utility, no Call option #

rand_vector4<-rnorm(Maturi*Frequ)
prop_Actif_Risque<-(alpha-r_no_risk)/(gamma*sigma^2)
S_t0<-1
B_t0<-1

S_t_simulated0<-lapply(rand_vector4,function (x)exp((alpha-0.5*sigma^2)/Frequ+sigma*sqrt(1/Frequ)*x))
B_t_simulated0<-rep(exp(r_no_risk/Frequ),Frequ*Maturi)
lambda_0<-(1/budget)^(gamma)*exp((-r_no_risk*gamma+r_no_risk-(0.5*theta_sim^2)+theta_sim^2/(2*gamma))*Maturi)
Ptf_simulated_0<-lambda_0^(-1/gamma)*exp((0.5*theta_sim^2*(1-1/gamma)^2-(r_no_risk+0.5*theta_sim^2)*(1-1/gamma))*Maturi)

Montant_Actif_Risque0<-prop_Actif_Risque*Ptf_simulated_0/S_t0
Montant_Actif_Sans_Risque0<-(Ptf_simulated_0-Montant_Actif_Risque0*S_t0)/B_t0
tool_investment0<-rbind(Montant_Actif_Risque0,Montant_Actif_Sans_Risque0)
  
S_t_simulated1<-cumprod(S_t_simulated0)
B_t_simulated1<-cumprod(B_t_simulated0)

tool_Market0<-rbind(S_t_simulated1[1:(Frequ*Maturi-1)],B_t_simulated1[1:(Frequ*Maturi-1)])
tool_Market1<-rbind(S_t_simulated1,B_t_simulated1)
Ptf_simulated_1<-tool_Market1[,1]%*%tool_investment0
tool_investment1<-rbind(prop_Actif_Risque*Ptf_simulated_1/tool_Market1[1,1],(Ptf_simulated_1-prop_Actif_Risque*Ptf_simulated_1)/tool_Market1[2,1])
#tool_Market1[,1]%*%tool_investment1==tool_Market1[,1]%*%tool_investment0




#Retourne un vecteur contenant les nouveaux montants investis ds l'actif risqué et sans risque, respectivement
Investment_Fonct<-function(S_t,B_t,invest_risq_tm1,invest_srisq_tm1,propor){
  Ptf_ini<-S_t*invest_risq_tm1+invest_srisq_tm1*B_t
  
  investment<-rbind(propor*Ptf_ini/S_t,Ptf_ini*(1-propor)/B_t)

  return(investment)
}
#Investment_Fonct(S_t=1.004773,B_t= 1.000385,invest_risq_tm1=0.2797681,invest_srisq_tm1=0.8393042,propor=0.25)



#Retourne le processus d'investissement de portfeuille
tool_investment[,1]<-tool_investment0
  for (i in 2:(Frequ*Maturi)){
    tool_investment[,i]<-Investment_Fonct(S_t=S_t_simulated1[(i-1)],B_t= B_t_simulated1[(i-1)],invest_risq_tm1=tool_investment[1,(i-1)],invest_srisq_tm1=tool_investment[2,(i-1)],propor=prop_Actif_Risque) 
  }

#Valeur du portfeuille
Ptf_simulated1<-apply(tool_Market1*tool_investment,2,function(x) sum(x))




#### "Frequ" - SIMULATIONS - Market ########
# Power Utility, no Call option #
########### Plusieurs Simulations ###########
tool_investment<-rbind(rep(0,Frequ*Maturi)+1,rep(0,Frequ*Maturi)+1) #Initialisation
tool_invest_risque_multi<-matrix(0,N_Simulations,Frequ*Maturi+1)
tool_invest_srisque_multi<-matrix(0,N_Simulations,Frequ*Maturi+1)
Ptf_simulated1_multi<-matrix(0,N_Simulations,Frequ*Maturi+1)
S_t_simulated1_multi<-matrix(0,N_Simulations,Frequ*Maturi+1)

Montant_Actif_Risque0<-prop_Actif_Risque*Ptf_simulated_0/S_t0
Montant_Actif_Sans_Risque0<-(Ptf_simulated_0-Montant_Actif_Risque0*S_t0)/B_t0
tool_investment0<-rbind(Montant_Actif_Risque0,Montant_Actif_Sans_Risque0)

rand_vector5<-matrix(rnorm(Maturi*Frequ*N_Simulations),N_Simulations,Maturi*Frequ)
S_t_simulated0_multi<-exp((alpha-0.5*sigma^2)/Frequ+sigma*sqrt(1/Frequ)*rand_vector5)
tool_invest_risque_multi[,1]<-tool_investment0[1]
tool_invest_srisque_multi[,1]<-tool_investment0[2]
B_t_simulated0<-rep(exp(r_no_risk/Frequ),Frequ*Maturi)
B_t_simulated0.5<-c(1,B_t_simulated0)
B_t_simulated1<-cumprod(B_t_simulated0.5)


#Retourne l'espérance de l'utilité pour plusieurs simulations
Ptf_optimal_InvStg<-function(props){
  for(j in 1:(N_Simulations)){
  S_t_simulated1_multi[j,]<-c(1,cumprod(S_t_simulated0_multi[j,]))
  
  
  for (i in 2:(Frequ*Maturi+1)){
    tool_each<-Investment_Fonct(S_t=S_t_simulated1_multi[j,i],B_t= B_t_simulated1[i],invest_risq_tm1=tool_invest_risque_multi[j,(i-1)],invest_srisq_tm1=tool_invest_srisque_multi[j,(i-1)],propor=props) 
    tool_invest_risque_multi[j,i]<-tool_each[1]
    tool_invest_srisque_multi[j,i]<-tool_each[2]
  }
  Ptf_simulated1_multi[j,]<-tool_invest_risque_multi[j,]*S_t_simulated1_multi[j,]+tool_invest_srisque_multi[j,]*B_t_simulated1
}
  Utility_T<-1/(1-gamma)*(Ptf_simulated1_multi[,(Frequ*Maturi+1)])^(1-gamma)
  Expected_Utility<-mean(Utility_T)
  return(Expected_Utility)
}

#Ptf_optimal_InvStg(props=0.5) #-0.8181148

#Retourne un vecteur des valeurs terminales 
Ptf_optimal_InvStg_Terminal<-function(props){
  for(j in 1:(N_Simulations)){
    S_t_simulated1_multi[j,]<-c(1,cumprod(S_t_simulated0_multi[j,]))
    
    
    for (i in 2:(Frequ*Maturi+1)){
      tool_each<-Investment_Fonct(S_t=S_t_simulated1_multi[j,i],B_t= B_t_simulated1[i],invest_risq_tm1=tool_invest_risque_multi[j,(i-1)],invest_srisq_tm1=tool_invest_srisque_multi[j,(i-1)],propor=props) 
      tool_invest_risque_multi[j,i]<-tool_each[1]
      tool_invest_srisque_multi[j,i]<-tool_each[2]
    }
    Ptf_simulated1_multi[j,]<-tool_invest_risque_multi[j,]*S_t_simulated1_multi[j,]+tool_invest_srisque_multi[j,]*B_t_simulated1
  }

  return(Ptf_simulated1_multi[,(Frequ*Maturi+1)])
}
############# Vérification de l'utilité à maturité #############











####      "Frequ"     - SIMULATIONS - Market ########
# Power Utility, with Call option & Optimal Proportion #
lambda_opt<-lambda_optimal(r_FL=r_no_risk,alpha_FL=alpha,sigma_FL=sigma,gamma_FL=gamma,K_FL=K_call_sim,b_FL=b_call_sim,a_FL=a_call_sim,X_0_FL=budget,T_FL=Maturi)
x_concavi<-Find_x_theta_PU(b_call_sim,a_call_sim,K_call_sim,gamma)
  
prob_sim<-((r_no_risk+(1/gamma-0.5)*theta_sim^2)*Maturi+log(x_concavi^(-gamma)/lambda_opt))/(theta_sim*sqrt(Maturi))
  
prop_Actif_Risque0<-(theta_sim/gamma+dnorm(prob_sim)/(sqrt(Maturi)*pnorm(prob_sim)))/sigma

#lambda_opt^(-1/gamma)/(a_call_sim^(1-1/gamma))*exp((r_no_risk+0.5*theta_sim^2)*Maturi*(1-1/gamma)+0.5*theta_sim^2*Maturi*(1-1/gamma)^2)*((-1/gamma)*pnorm(((r_no_risk+theta_sim^2)*Maturi-theta_sim^2*Maturi*(1-1/gamma)-indica)/(theta_sim*sqrt(Maturi)))+dnorm(((r_no_risk+theta_sim^2)*Maturi-theta_sim^2*Maturi*(1-1/gamma)-indica)/(theta_sim*sqrt(Maturi))))+tool2
#indica<-((lambda_opt*1/a_call_sim)^(-1/gamma)-K_call_sim)/a_call_sim+b_call_sim 


tool_investment<-rbind(rep(0,Frequ*Maturi)+1,rep(0,Frequ*Maturi)+1) #Initialisation
tool_invest_risque_multi<-matrix(0,N_Simulations,Frequ*Maturi+1)
tool_invest_srisque_multi<-matrix(0,N_Simulations,Frequ*Maturi+1)
Ptf_simulated1_multi<-matrix(0,N_Simulations,Frequ*Maturi+1)
S_t_simulated1_multi<-matrix(0,N_Simulations,Frequ*Maturi+1)
xi_t_simulated1_multi<-matrix(0,N_Simulations,Frequ*Maturi+1)
xi_t_simulated1_multi[,1]<-1

Montant_Actif_Risque0<-prop_Actif_Risque0*Ptf_simulated_0/S_t0
Montant_Actif_Sans_Risque0<-(Ptf_simulated_0-Montant_Actif_Risque0*S_t0)/B_t0
tool_investment0<-rbind(Montant_Actif_Risque0,Montant_Actif_Sans_Risque0)

rand_vector5<-matrix(rnorm(Maturi*Frequ*N_Simulations),N_Simulations,Maturi*Frequ)
S_t_simulated0_multi<-exp((alpha-0.5*sigma^2)/Frequ+sigma*sqrt(1/Frequ)*rand_vector5)
tool_invest_risque_multi[,1]<-tool_investment0[1]
tool_invest_srisque_multi[,1]<-tool_investment0[2]
B_t_simulated0<-rep(exp(r_no_risk/Frequ),Frequ*Maturi)
B_t_simulated0.5<-c(1,B_t_simulated0)
B_t_simulated1<-cumprod(B_t_simulated0.5)
Payoff<-rep(0,N_Simulations)
x_conca<-Find_x_theta_PU(b_o=1,a_o=1,K_o=1,gamma_o=gamma)
lagang<-lambda_optimal(r_FL=r_no_risk,alpha_FL=alpha,sigma_FL=sigma,gamma_FL=gamma,K_FL=1,b_FL=1,a_FL=1,X_0_FL=1,T_FL=Maturi)
vecteur_proportion_opt<-rep(0,(Maturi*Frequ+1)) #Only for the last simulation, so this is why it's a vector

Ptf_optimal_InvStg_Call<-function(a_call,b_call,K_call){ #proportion non-cte and return E(Utility)
  for(j in 1:(N_Simulations)){
    S_t_simulated1_multi[j,]<-c(1,cumprod(S_t_simulated0_multi[j,]))
    
    for (i in 2:(Frequ*Maturi+1)){
      t_sim<-(i-1)/Frequ
      xi_t_simulated1_multi[j,i]<-S_t_simulated1_multi[j,i]^(-theta_sim/sigma)*exp((theta_sim*alpha/sigma-theta_sim*sigma/2-r_no_risk-0.5*theta_sim^2)*t_sim)
      
      if(i!=(Frequ*Maturi+1)){tool3<-((r_no_risk+(1/gamma-0.5)*theta_sim^2)*(Maturi-t_sim)+log(x_conca^(-gamma)/(lagang*xi_t_simulated1_multi[j,i])))/(theta_sim*sqrt(Maturi-t_sim))
    
      props_optimal<-(theta_sim/gamma+dnorm(tool3)/(pnorm(tool3)*sqrt(Maturi-t_sim)))/sigma
      vecteur_proportion_opt[i]<-props_optimal 
      
      tool_each<-Investment_Fonct(S_t=S_t_simulated1_multi[j,i],B_t= B_t_simulated1[i],invest_risq_tm1=tool_invest_risque_multi[j,(i-1)],invest_srisq_tm1=tool_invest_srisque_multi[j,(i-1)],propor=props_optimal) 
      tool_invest_risque_multi[j,i]<-tool_each[1]
      tool_invest_srisque_multi[j,i]<-tool_each[2]
      }
    }
    
    Ptf_simulated1_multi[j,]<-tool_invest_risque_multi[j,]*S_t_simulated1_multi[j,]+tool_invest_srisque_multi[j,]*B_t_simulated1
    Ptf_simulated1_multi[j,(Frequ*Maturi+1)]<-tool_invest_risque_multi[j,(Frequ*Maturi)]*S_t_simulated1_multi[j,(Frequ*Maturi+1)]+tool_invest_srisque_multi[j,(Frequ*Maturi)]*B_t_simulated1[(Frequ*Maturi+1)]
 
     Payoff[j]<-a_call*max(0,Ptf_simulated1_multi[j,(Frequ*Maturi+1)]-b_call)+K_call
      }



  Utility<-1/(1-gamma)*(Payoff)^(1-gamma)
  Expected_Utility<-mean(Utility)
  return(Expected_Utility)
}

#Ptf_optimal_InvStg_Call(a_call=1,b_call=1,K_call=1)# -0.8199844



Ptf_optimal_InvStg_Call_terminal<-function(a_call,b_call,K_call){ #proportion cte and Utility simulated
  for(j in 1:(N_Simulations)){
    S_t_simulated1_multi[j,]<-c(1,cumprod(S_t_simulated0_multi[j,]))
    
    for (i in 2:(Frequ*Maturi+1)){
      t_sim<-(i-1)/Frequ
      xi_t_simulated1_multi[j,i]<-S_t_simulated1_multi[j,i]^(-theta_sim/sigma)*exp((theta_sim*alpha/sigma-theta_sim*sigma/2-r_no_risk-0.5*theta_sim^2)*t_sim)
      
      if(i!=(Frequ*Maturi+1)){tool3<-((r_no_risk+(1/gamma-0.5)*theta_sim^2)*(Maturi-t_sim)+log(x_conca^(-gamma)/(lagang*xi_t_simulated1_multi[j,i])))/(theta_sim*sqrt(Maturi-t_sim))
      
      props_optimal<-(theta_sim/gamma+dnorm(tool3)/(pnorm(tool3)*sqrt(Maturi-t_sim)))/sigma
      vecteur_proportion_opt[i]<-props_optimal 
      
      tool_each<-Investment_Fonct(S_t=S_t_simulated1_multi[j,i],B_t= B_t_simulated1[i],invest_risq_tm1=tool_invest_risque_multi[j,(i-1)],invest_srisq_tm1=tool_invest_srisque_multi[j,(i-1)],propor=props_optimal) 
      tool_invest_risque_multi[j,i]<-tool_each[1]
      tool_invest_srisque_multi[j,i]<-tool_each[2]
      }
    }
    
    Ptf_simulated1_multi[j,]<-tool_invest_risque_multi[j,]*S_t_simulated1_multi[j,]+tool_invest_srisque_multi[j,]*B_t_simulated1
    Ptf_simulated1_multi[j,(Frequ*Maturi+1)]<-tool_invest_risque_multi[j,(Frequ*Maturi)]*S_t_simulated1_multi[j,(Frequ*Maturi+1)]+tool_invest_srisque_multi[j,(Frequ*Maturi)]*B_t_simulated1[(Frequ*Maturi+1)]
    
    Payoff[j]<-a_call*max(0,Ptf_simulated1_multi[j,(Frequ*Maturi+1)]-b_call)+K_call
  }
 
  return(Payoff)
}

#Ptf_optimal_InvStg_Call_terminal(a_call=1,b_call=1,K_call=1)




####################### TESTS SECTION #######################
#Market a Portfolio Overview
par(mfrow=c(3,1))
plot(0:520,c(1,S_t_simulated1),type='l',xlab=' ',ylab=' ',main='Actif risqué')
plot(0:520,c(1,B_t_simulated1),type='l',xlab=' ',ylab=' ',main='Actif sans risque')
plot(0:520,c(Ptf_simulated_0,Ptf_simulated1),type='l',xlab=' ',ylab=' ',main='Portefeuille')

#Utility maximisation, Original Problem 
par(mfrow=c(1,1))
verif<-lapply(seq(0,1,0.05),function(x) Ptf_optimal_InvStg(props=x))
plot(seq(0,1,0.05),verif,type='l',xlab=c( expression(nu[t])),ylab=expression(paste(E,'[ U(',X[T],') ]') ))
abline(v =prop_Actif_Risque, col = "blue", lwd = 2,lty=2)
legend(0.4,-0.80,col='blue', lwd = 2,lty=2,legend=c(expression(paste(nu[t],'*'," = ",frac(alpha-r,gamma*sigma^2)))))

#utility maximisation, On Call Option
verif_Call<-lapply(seq(0,1,0.05),function(x) Ptf_optimal_InvStg_Call(props=x,a_call=1,b_call=1,K_call=1))
plot(seq(0,1,0.05),verif_Call,type='l',xlab=c( expression(nu[t])),ylab=expression(paste(E,'[ U(',X[T],') ]') ))
abline(v =prop_Actif_Risque, col = "blue", lwd = 2,lty=2)
legend(0.4,-0.80,col='blue', lwd = 2,lty=2,legend=c(expression(paste(nu[t],'*'," = ",frac(alpha-r,gamma*sigma^2)))))



########## Histogramme Portefeuille optimal: dynamique et forme fermée ###############
proportion_optimale<-(alpha-r_no_risk)/(gamma*sigma^2)

# Start the clock!
ptm <- proc.time()

Ptf_dyn<-data.frame(valeur_opt=Ptf_optimal_InvStg_Terminal(props=proportion_optimale))
Ptf_theorique<-data.frame(valeur_opt=as.numeric(X_Optimal_base_simulated))

Ptf_dyn$comp<-'Dyn'
Ptf_theorique$comp<-'Theo'

Compa_combin<-rbind(Ptf_dyn,Ptf_theorique)

ggplot(data=Compa_combin,aes(Compa_combin$valeur_opt,group=comp,fill=comp))+
  geom_histogram(colour='black',binwidth = 0.02,alpha=0.5,position = "identity")+
  scale_fill_manual(name="",values=c("yellow","darkgray"),labels=c( expression(paste("X","*"[T],": Dynamique")),expression(paste("X","*"[T],": Théorique"))))+
  labs(x=expression(paste("X","*"[T])), y="# Réalisations")+
  theme(axis.title=element_text(size=16,face="bold"))+
  theme_classic()+
  theme(legend.position=c(0.80,0.65))+
  theme(legend.text=element_text(size=13))+

  geom_vline(aes(xintercept =mean(as.numeric(X_Optimal_base_simulated)), linetype = "Moyenne Théorique: 1.283896"), colour= 'black', data = Compa_combin) +                                                                                                     
  geom_vline(aes(xintercept =mean(as.numeric(Ptf_dyn$valeur_opt)), linetype = "Moyenne Dynamique:1.284377"), colour= 'purple4', data = Compa_combin) +
  scale_linetype_discrete(name = " ")

# Stop the clock
proc.time() - ptm   #136.605 secondes


#  scale_color_manual("Statistics", values = c("Mean" <-"red", "Median" <- "green"))+
#scale_linetype_manual(name = "Thresholds",values = c(2,2))  










########## Histogramme Portefeuille optimal: dynamique et densité ###############
proportion_optimale<-(alpha-r_no_risk)/(gamma*sigma^2)
lbd<-(1/budget)^(gamma)*exp((-r_no_risk*gamma+r_no_risk-(0.5*theta_sim^2)+theta_sim^2/(2*gamma))*Maturi)
mu_mod<-((r_no_risk+0.5*theta_sim^2)*Maturi-log(lbd))/gamma
sigma_mod2<-theta_sim^2*Maturi/gamma^2

mean_theo<-exp(mu_mod+0.5*sigma_mod2)  


# Start the clock!
ptm <- proc.time()

Ptf_dyn<-data.frame(valeur_opt=as.numeric(Ptf_simul))#Ptf_simul is obtained from the Memoire_Codes3.3, 3th method.
x <- seq(0.5, 2, length.out=100)
df <- with(Ptf_dyn, data.frame(x = x, y = dlnorm(x, mu_mod,sqrt(sigma_mod2))))

Ptf_dyn$comp<-'Dyn'


ggplot(data=Ptf_dyn,aes(Ptf_dyn$valeur_opt))+
  geom_histogram(aes(y=..density..),colour='black',binwidth = 0.05,alpha=0.5,position = "identity")+
  geom_line(data = df, aes(x = x, y = y), color = "red")+
  scale_fill_manual(values=c("firebrick2", "grey"))+
  scale_x_continuous(breaks=NULL)+
  scale_y_continuous(breaks=NULL)+
  labs(title='',x=expression(paste("X","*"[T])), y="Densité")+
  theme(axis.title=element_text(size=16,face="bold"))+
  theme_classic()+
  theme(legend.position=c(0.75,0.8))+
  #theme(legend.text=element_text(size=13))+
  #geom_vline(aes(xintercept =mean_theo, linetype = "Moyenne Théorique: 1.268471")) +                                                                                                     
  #geom_vline(aes(xintercept =mean(as.numeric(Ptf_dyn$valeur_opt)), linetype = "Moyenne Dynamique:  1.286424"), data = Ptf_dyn)+                                                                                                     
  #theme(plot.title = element_text(family = "Helvetica", face = "bold", size = (15),hjust = 0.5))+
  scale_linetype_discrete(name = " ")
 
 # Stop the clock
proc.time() - ptm   #136.605 secondes






########## Histogramme Portefeuille optimal Call: dynamique et forme fermée ###############

ptm_call <- proc.time()

Ptf_dyn_call<-data.frame(valeur_opt=as.numeric(Ptf_optimal_InvStg_Call_terminal(a_call=1,b_call=1,K_call=1)))
Ptf_theorique_call<-data.frame(valeur_opt=as.numeric(X_Optimal_simulated)) #(Voir Memoire_Codes1)

Ptf_dyn_call$comp<-'Dyn'
Ptf_theorique_call$comp<-'Theo'

Compa_combin_call<-rbind(Ptf_dyn_call,Ptf_theorique_call)

ggplot(data=Compa_combin_call,aes(Compa_combin_call$valeur_opt,group=comp,fill=comp))+
  geom_histogram(colour='black',binwidth = 0.02,alpha=0.5,position = "identity")+
  scale_fill_manual(name="",values=c("yellow","darkgray"),labels=c( expression(paste("X","*"[T],": Dynamique")),expression(paste("X","*"[T],": Théorique"))))+
  labs(x=expression(paste("X","*"[T])), y="# Réalisations")+
  theme_classic()+
  geom_vline(xintercept = 1.9999793, linetype="dotted")

proc.time() - ptm_call



################# - Comparaison - #################
#### Histogramme Portefeuille optimal Call: dynamique (16-12-2020) ###############

ptm_call <- proc.time()

Ptf_final_2448_0<-data.frame(valeur_opt=as.numeric(ptf_2448_0))#Obtenu à l'aide de la méthode 4 MemoireCodes3.3
Ptf_final_1224_1224<-data.frame(valeur_opt=as.numeric(ptf_1224_1224))#Obtenu à l'aide de la méthode 4 MemoireCodes3.3
Ptf_final_18_6448<-data.frame(valeur_opt=as.numeric(ptf_18_6448))#Obtenu à l'aide de la méthode 4 MemoireCodes3.3
Ptf_final_0_0<-data.frame(valeur_opt=as.numeric(ptf_0_0))#Obtenu à l'aide de la méthode 4 MemoireCodes3.3

Ptf_final_2448_0$comp<-'2.c_s=0.000% et c_f=2.448%'
Ptf_final_1224_1224$comp<-'3.c_s=1.224% et c_f=1.224%'
Ptf_final_18_6448$comp<-'4. c_s=1.800% et c_f=0.6448%'
Ptf_final_0_0$comp<-'1. Aucun frais'


Combin_ptf_frais<-rbind(Ptf_final_2448_0,Ptf_final_1224_1224,Ptf_final_18_6448,Ptf_final_0_0)


ggplot(data=Combin_ptf_frais,aes(Combin_ptf_frais$valeur_opt,group=comp,fill=comp))+
  geom_histogram(colour='black',binwidth = 0.01,alpha=0.7,position = "identity")+
  labs(x=expression(paste("X","*"[T])), y="# Réalisations")+
  scale_x_continuous(breaks=c(x_concavi),labels=expression(paste(widehat(x),'(',D[T],')')))+
  scale_fill_manual(name="",values=c("#F8766D","#00B81F","#00A5FF","#E76BF3"),labels=c( expression(paste("1. Aucun frais")), expression(paste("2. ",c[s],"=0.000% et ",c[f],"=2.448%")),expression(paste("2. ",c[s],"=0.000% et ",c[f],"=2.448%")),expression(paste("2. ",c[s],"=0.000% et ",c[f],"=2.448%"))))+
  theme_classic()+
  theme(legend.position = 'bottom',legend.title = element_blank())

proc.time() - ptm_call















####### Graphique proportion qui varie dans le temps pour le même S_t (20-01-2020)########
S_test<-0.9
xi_tilde<-S_test^(-theta_sim_tilde/sigma)*exp((theta_sim_tilde*alpha_tilde/sigma-theta_sim_tilde*sigma/2-r_no_risk_tilde-0.5*theta_sim_tilde^2)*t_Test)

proportion_test6<-lapply(seq(0,9.9,0.1),function(x){        
  xi_test<-S_test^(-theta_sim_tilde/sigma)*exp((theta_sim_tilde*alpha_tilde/sigma-theta_sim_tilde*sigma/2-r_no_risk_tilde-0.5*theta_sim_tilde^2)*x)
  interior<-((r_no_risk_tilde+(1/gamma-0.5)*theta_sim_tilde^2)*(Maturi- x)+log((x_concavi^(-gamma))/(xi_test*lambda_opt)))/(theta_sim_tilde*sqrt(Maturi-x))  
  return((theta_sim_tilde/gamma+dnorm(interior)/(pnorm(interior)*sqrt(Maturi-x)))/sigma)
})

plot(seq(0,9.9,0.1),proportion_test1,type='l',xlab="Temps (t)",ylab=expression(nu[t]),ylim=c(0,6))
lines(seq(0,9.9,0.1),proportion_test2,col="red")
lines(seq(0,9.9,0.1),proportion_test3,col="blue")
lines(seq(0,9.9,0.1),proportion_test4,col="purple")
lines(seq(0,9.9,0.1),proportion_test5,col="blue4")
lines(seq(0,9.9,0.1),proportion_test6,col="orange")
legend(0, 6, legend=c(expression(paste(S[t],"=0.9")), expression(paste(S[t],"=1")), expression(paste(S[t],"=1.1")),expression(paste(S[t],"=1.2")),expression(paste(S[t],"=1.3")),expression(paste(S[t],"=1.4"))),
       col=c("orange","black", "red","blue","purple","blue4"), lty=c(1,1,1,1,1), cex=0.8)

#### TEST ####
par(mfrow=c(1,1))
tre<-c(1,1/0)
tre[2]<-4




