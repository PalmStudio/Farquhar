
require(stats4)
require(ggplot2)
require(dplyr)

# Functions used for fitting Farquahr parameters from A-Ci and A-PFD curves ------------------------------------------------------------------

# 
# 
# # Parameters --------------------------------------------------------------
# Tref=298.16 # Reference temperature (K)
# Tleaf=28.5+273.16 #Temperature of the leaflet
# abso=0.85 #Absorptance of the leaflet (unitless)
# f=0.15 # Correcting factor for the spectral quality of the light (unitless)
# O=210 #Intracellular concentration of O2 (mmol.mol-1)
# R=8.314 # Ideal gas constant (JK−1 mol−1)
# 
# GstariRef=42.75 # Value of Gstar at the reference temperature (micromol.m-2.s-1)
# GstariHa=37830 # Enthalpy of activation of the temperature response of Gstar (J.mol-1)
# 
# JmaxHa=50300	# Enthalpy of activation of the temperature response of Jmax (J.mol-1)
# JmaxHd=152044 # Enthalpy of desactivation of the temperature response of Jmax (J.mol-1 )
# JmaxS=495 #Entropie of the temperature response of Jmax (J.mol-1.K-1 )
# 
# # cRd=18.714
# HaRd=46390 # Enthalpy of activation of the temperature response of Rd (J.mol-1)
# 
# Kcref=404.90 #Value of kc at Tr (micromol.m-2.s-1 )
# # cKc=32.05
# HaKc=79430 # Enthalpy of activation of the temperature response of kc (J.mol-1)
# 
# Koref=278.4 # Michaelis-Menten constant of Rubisco for O2 (Importance de mettre O en mmol!) (mmol.m-2.s-1 )
# cKo=14.67 
# HaKo=36380 # Enthalpy of activation of the temperature response of ko (J.mol-1 )
# 
# Hav=73637 # Enthalpy of activation of the temperature response of V cmax (J.mol-1)
# Hdv=149252 # Enthalpy of desactivation of the temperature response of V cmax (J.mol-1)
# sv=486 # Entropie of the temperature response of V cmax (J.mol)
# 
# Haj=50300 #Enthalpy of activation of the temperature response of Jmax (J mol-1)
# Hdj=152044 #Enthalpy of desactivation of the temperature response of Jmax (J mol-1)
# sj=495 #Entropie of the temperature response of Jmax (J mol-1 K-1)
# 
# theta=0.853 #Empirical curvature factor for the response of J to PFD
# O=210 #Intracellular concentration of O2 (mmol mol−1)



harleytemp<-function(Pref,Ha,Tleaf){
  P=Pref*exp(Ha/(R*Tref)-Ha/(R*Tleaf))
}

leuningtemp<-function(Pref,Ha,Hd,s,Tleaf){
  P=Pref*(1+exp((s*Tref-Hd)/(R*Tref)))*exp(Ha/(R*Tref)*(1-Tref/Tleaf))/(1+exp((s*Tleaf-Hd)/(R*Tleaf)))
}


Wc<-function(Vcmax=100,Ci=400,Kc=400,O2=400,K0=200){
  Wc<-Vcmax*Ci/(Ci+Kc*(1+O2/K0))
  return(Wc)
}

WJ<-function(Ci=400,G=43,Jmax=120,PFD=1500,theta=theta){
  I2=PFD*abso*(1-f)/2
  J=(I2+Jmax-((I2+Jmax)^2-4*theta*I2*Jmax)^0.5)/(2*theta)
  Wj<-J*Ci/(4*Ci+8*G)
  return(Wj)
}


fA<-function(Ci=400,PFD=1500,O2=210,Vcmax=200,Jmax=100,Rd=1,Tleaf=Tleaf){
  Kc=harleytemp(Kcref,HaKc,Tleaf)    
  K0=harleytemp(Koref,HaKo,Tleaf)
  G=harleytemp(GstariRef,GstariHa,Tleaf)
  Rd=harleytemp(Rd,HaRd,Tleaf)  
  Jmax=leuningtemp(Jmax,Haj,Hdj,sj,Tleaf)
  Vcmax=leuningtemp(Vcmax,Hav,Hdv,sv,Tleaf)
  Ac<-Wc(Vcmax=Vcmax,Ci=Ci,Kc=Kc,O2=O2,K0=K0)*(1-G/Ci)-Rd
  Aj<-WJ(Ci=Ci,PFD=PFD,Jmax=Jmax,G=G,theta=theta)*(1-G/Ci)-Rd
  A<-pmin(Wc(Vcmax=Vcmax,Ci=Ci,Kc=Kc,O2=O2,K0=K0)*(1-G/Ci),WJ(Ci=Ci,Jmax=Jmax,G=G,PFD=PFD,theta=theta)*(1-G/Ci))-Rd
  output=data.frame(A,Aj,Ac)
  
  return(output)
}


#Assimilation formula without Rubisco limitation 
fA_light<-function(theta=0.75,Ci=400,PFD=1500,O2=210,Jmax=100,Rd=1,Tleaf=Tleaf){
  Kc=harleytemp(Kcref,HaKc,Tleaf)    
  Ko=harleytemp(Koref,HaKo,Tleaf)
  G=harleytemp(GstariRef,GstariHa,Tleaf)
  Rd=harleytemp(Rd,HaRd,Tleaf)  
  Jmax=leuningtemp(Jmax,JmaxHa,JmaxHd,JmaxS,Tleaf)
  A<-WJ(Ci=Ci,PFD=PFD,Jmax=Jmax,G=G,theta=theta)*(1-G/Ci)-Rd
  sortie=data.frame(A)
  return(sortie)
}

# Estimate Jmax, Vcmaxref theta and Rd ------------------------------

###fitting from A-PFD curves

SumSq<-function(data,Vcmax=200,Jmax=300,Rd=1){
  y<-data$A-fA(Vcmax=Vcmax,Jmax=Jmax,Rd=Rd,Ci=data$Ci,PFD=data$PFD,Tleaf=data$Tleaf+273.16,O2=210)$A
  return(sum(y^2))
}

SumSq2<-function(data,par){
  return(SumSq(data,Vcmax=par[1],Jmax=par[2],Rd=par[3]))
}

CO2fitting<-function(don,plot=T){
  MoindresCarres<-NULL
  out=NULL
  grph=NULL
  StartMLE<-list(Vcmax=200,Jmax=300,Rd=0.9,sigma=3.31)
  try({
    MoindresCarres<-optim(par=list(Vcmax=200,Jmax=300,Rd=0.9),fn=SumSq2,data=don)
    # print(MoindresCarres)
    StartMLE<-list(Vcmax=MoindresCarres$par[["Vcmax"]],
                   Jmax=MoindresCarres$par[["Jmax"]],
                   Rd=MoindresCarres$par[["Rd"]],
                   sigma=sqrt(MoindresCarres$value/nrow(don)))
  })
  
  
  if (!is.null(MoindresCarres)) { 
    out=don%>%
      mutate(Vcmax=MoindresCarres$par["Vcmax"],
             Jmax=MoindresCarres$par["Jmax"],
             Rd=MoindresCarres$par["Rd"],
             A_sim=fA(Ci=Ci,PFD=PFD,Tleaf=Tleaf+273.16,Jmax=MoindresCarres$par["Jmax"],Vcmax=MoindresCarres$par["Vcmax"],Rd=MoindresCarres$par["Rd"])$A)
    
    
    out=out%>%
      mutate(RMSE=statsimu(X_obs = Ci,Y_obs = A,Y_estim = A_sim)$RMSE,
             RRMSE=statsimu(X_obs = Ci,Y_obs = A,Y_estim = A_sim)$RRMSE,
             Bias=statsimu(X_obs = Ci,Y_obs = A,Y_estim = A_sim)$Bias)
    
    
    
    grph=out%>%
      ggplot()+
      geom_point(aes(x=Ci,y=A,color='Obs'),size=2)+
      geom_line(aes(x=Ci,y=A_sim,color='Sim'))+
      xlab(expression(Ci~(mol~mol^-1)))+
      ylab(expression(A~(mu~mol~m^-2~s^-1)))+
      theme(legend.title=element_blank())+
      annotate(geom = 'text',x =1.2*min(out$Ci) ,y=max(out$A),label=paste('RRMSE=',unique(round(out$RRMSE,3)),'\nBias=',unique(round(out$Bias,3))))+
      annotate(geom = 'text',x =0.8*max(out$Ci) ,y=0.5*max(out$A),label=paste('Vcmax=',unique(round(out$Vcmax,1)),'\nJmax=',unique(round(out$Jmax,1)),'\nRd=',unique(round(out$Rd,2))))+
      ggtitle(paste(unique(out$Progeny),'Tree',unique(out$Tree),'Frond',unique(out$Frond)))
    if (plot==T){
      print(grph)
    }
    
  }
  
  return(list(data=out,graph=grph))
}



###fitting from A-PFD curves
# MinusLogL_light<-function(Jmax=120,theta=0.75,Rd=1,sigma=0.2){
#   data=don
#   y<-dnorm(x=data$A,mean=fA_light(theta=theta,Jmax=Jmax,Rd=Rd,Ci=data$Ci,PFD=data$PFD,Tleaf=data$Tleaf+273.16,O2=210)$A,sd=sigma,log=TRUE)
#   return(-sum(y))
# }

SumSq_light<-function(data,Jmax=150,theta=0.75,Rd=1){
  y<-data$A-fA_light(theta=theta,Jmax=Jmax,Rd=Rd,Ci=data$Ci,PFD=data$PFD,Tleaf=data$Tleaf+273.16,O2=210)$A
  return(sum(y^2))
}
SumSq2_light<-function(data,par){
  return(SumSq_light(data,Jmax=par[1],theta=par[2],Rd=par[3]))
}

Lightfitting=function(don,plot=T){
  out=NULL
  grph_light=NULL
  MoindresCarres<-NULL
  StartMLE<-list(Jmax=240,theta=0.85,Rd=0.9,sigma=3.31)
  try({
    MoindresCarres<-optim(par=list(Jmax=240,theta=0.84,Rd=0.9),fn=SumSq2_light,data=don)
    # print(MoindresCarres)
    StartMLE<-list(Jmax=MoindresCarres$par[["Jmax"]],
                   theta=MoindresCarres$par[["theta"]],
                   Rd=MoindresCarres$par[["Rd"]],
                   sigma=sqrt(MoindresCarres$value/nrow(don)))
  })
  
  
  # Estimation<-NULL
  # try({
  #   Estimation<-mle(minuslog=MinusLogL_light,start=StartMLE,nobs=nrow(don))
  # })
  
  
  # if (!is.null(Estimation)) {
  #   
  #   out=don%>%
  #     mutate(Jmax=coef(Estimation)["Jmax"],
  #            Rd=coef(Estimation)["Rd"],
  #            Theta=coef(Estimation)["theta"],
  #            Sigma=coef(Estimation)["sigma"],
  #            A_sim=fA_light(Ci=Ci,PFD=PFD,Tleaf=Tleaf+273.16,Jmax=coef(Estimation)["Jmax"],theta=coef(Estimation)["theta"],Rd=coef(Estimation)["Rd"])$A)
  
  if (!is.null(MoindresCarres)) {
    
    out=don%>%
      mutate(Jmax=MoindresCarres$par[["Jmax"]],
             theta=MoindresCarres$par[["theta"]],
             Rd=MoindresCarres$par[["Rd"]],
             sigma=sqrt(MoindresCarres$value/nrow(don)),
             A_sim=fA_light(Ci=Ci,PFD=PFD,Tleaf=Tleaf+273.16,Jmax=MoindresCarres$par[["Jmax"]],theta=MoindresCarres$par[["theta"]],Rd=MoindresCarres$par[["Rd"]])$A)
    
    out=out%>%
      mutate(RMSE=statsimu(X_obs = Ci,Y_obs = A,Y_estim = A_sim)$RMSE,
             RRMSE=statsimu(X_obs = Ci,Y_obs = A,Y_estim = A_sim)$RRMSE,
             Bias=statsimu(X_obs = Ci,Y_obs = A,Y_estim = A_sim)$Bias)
    
    
    
    grph_light=out%>%
      ggplot()+
      geom_point(aes(x=PFD,y=A,color='Obs'),size=2)+
      geom_line(aes(x=PFD,y=A_sim,color='Sim'))+
      xlab(expression(PFD~(mu~mol~m^-2~s^-1)))+
      ylab(expression(A~(mu~mol~m^-2~s^-1)))+
      theme(legend.title=element_blank())+
      annotate(geom = 'text',x =1.2*min(out$PFD) ,y=max(out$A),label=paste('RRMSE=',unique(round(out$RRMSE,3)),'\nBias=',unique(round(out$Bias,3))))+
      annotate(geom = 'text',x =0.8*max(out$PFD) ,y=0.5*max(out$A),label=paste('Jmax=',unique(round(out$Jmax,1)),'\nTheta=',unique(round(out$theta,1)),'\nRd=',unique(round(out$Rd,2))))+
      ggtitle(paste(unique(out$Progeny),'Frond',unique(out$Frond)))
    
    if (plot==T){
      print(grph_light)
    }
  }
  return(list(data=out,graph=grph_light))
}



####stats
statsimu=function(X_obs,Y_obs,Y_estim){
  
  nb=length(Y_obs)
  SCR = sum((Y_estim - Y_obs)**2)
  RMSE=sqrt(1/nb * SCR)
  RRMSE=sqrt(1/nb * SCR)/mean(Y_obs,na.rm=T)  
  bias = 1/nb * sum(Y_estim - Y_obs)
  
  return=data.frame(RMSE,RRMSE,bias,row.names='')
  colnames(return)=c("RMSE",'RRMSE','Bias')
  return(return)
}