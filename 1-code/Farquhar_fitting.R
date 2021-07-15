
# Packages ----------------------------------------------------------------

packs <- c( 'stats4','ggplot2','dplyr','cowplot','viridis','lubridate','data.table')
InstIfNec<-function (pack) {
  if (!do.call(require,as.list(pack))) {
    do.call(install.packages,as.list(pack))  }
  do.call(require,as.list(pack)) }
lapply(packs, InstIfNec)


# Theme -------------------------------------------------------------------
myTheme <- theme(
  panel.background=element_rect(fill="transparent", color=NA),  
  plot.background = element_rect(fill = "transparent",colour = NA),
  axis.line=element_line(colour="black"), 
  axis.title=element_text(size=18),
  axis.text.y=element_text(size=16, colour="black"), 
  axis.text.x=element_text(size=16, colour="black", angle=0, hjust=0.5),
  panel.grid.minor = element_blank(), 
  panel.grid.major = element_line(colour="grey90", size=0.2), 
  legend.position=c(0.8,0.1), 
  legend.text=element_text(size=14),
  legend.title=element_text(size=14)
)


# Parameters --------------------------------------------------------------
Tref=298.16 # Reference temperature (K)
Tleaf=28.5+273.16 #Temperature of the leaflet
abso=0.85 #Absorptance of the leaflet (unitless)
f=0.15 # Correcting factor for the spectral quality of the light (unitless)
O=210 #Intracellular concentration of O2 (mmol.mol-1)
R=8.314 # Ideal gas constant (JK−1 mol−1)

GstariRef=42.75 # Value of Gstar at the reference temperature (micromol.m-2.s-1)
GstariHa=37830 # Enthalpy of activation of the temperature response of Gstar (J.mol-1)

JmaxHa=50300	# Enthalpy of activation of the temperature response of Jmax (J.mol-1)
JmaxHd=152044 # Enthalpy of desactivation of the temperature response of Jmax (J.mol-1 )
JmaxS=495 #Entropie of the temperature response of Jmax (J.mol-1.K-1 )

# cRd=18.714
HaRd=46390 # Enthalpy of activation of the temperature response of Rd (J.mol-1)

Kcref=404.90 #Value of kc at Tr (micromol.m-2.s-1 )
# cKc=32.05
HaKc=79430 # Enthalpy of activation of the temperature response of kc (J.mol-1)

Koref=278.4 # Michaelis-Menten constant of Rubisco for O2 (Importance de mettre O en mmol!) (mmol.m-2.s-1 )
cKo=14.67 
HaKo=36380 # Enthalpy of activation of the temperature response of ko (J.mol-1 )

Hav=73637 # Enthalpy of activation of the temperature response of V cmax (J.mol-1)
Hdv=149252 # Enthalpy of desactivation of the temperature response of V cmax (J.mol-1)
sv=486 # Entropie of the temperature response of V cmax (J.mol)

Haj=50300 #Enthalpy of activation of the temperature response of Jmax (J mol-1)
Hdj=152044 #Enthalpy of desactivation of the temperature response of Jmax (J mol-1)
sj=495 #Entropie of the temperature response of Jmax (J mol-1 K-1)

theta=0.853 #Empirical curvature factor for the response of J to PFD
O=210 #Intracellular concentration of O2 (mmol mol−1)

# Import functions --------------------------------------------------------
source('1-code/Functions_Farquhar_Fitting.R')


# Import data -------------------------------------------------------------

###A-Ci curves
datC= 
  data.table::fread("0-data/CO2curveEX.csv", data.table = FALSE)%>%
  mutate(Sample=paste(Progeny,Tree,Frond,sep='_'))

###A-PFD curves
datPFD= data.table::fread("0-data/lightcurveEX.csv", data.table = FALSE)%>%
  mutate(Sample= paste(Progeny,Tree,Frond,Date,sep='_'),
         Date= dmy(Date))

# Fitting CO2 curves  -----------------------------------------------------
          
ResSampleC=NULL

for(i in unique(datC$Sample)){
  # i=unique(datC$Sample)[1]
  sub=datC%>%
    filter(Sample==i)%>%
    arrange(Ci)

  ####!!! data need to be updated each time
  res=CO2fitting(don =sub, plot=F)
  ResSampleC= rbind(ResSampleC,res$data)

}

# ResSampleC%>%
#   group_by(Progeny,Tree,Frond)%>%
#   summarize(Vcmax=mean(Vcmax),Jmax=mean(Jmax),Rd=mean(Rd),RMSE=mean(RMSE),RRMSE=mean(RRMSE))

# Graphics---- ####

###curve fitting per sample
ResSampleC%>%
  ggplot(aes(x=Ci,y=A,group=Sample,col=as.factor(Frond)))+
  geom_point(aes(x=Ci,y=A,group=Sample,col=as.factor(Frond)))+
  geom_line(aes(x=Ci,y=A_sim,group=Sample,col=as.factor(Frond)))+
  scale_color_viridis_d(name='Rank')+
  facet_wrap(Progeny~Tree, ncol = 2)+
  labs(y= "A (µmol m-2 s-1)", x= "Ci (ppm)")

ResSampleC%>%
  # group_by(Progeny,Tree,Frond)%>%
  # summarize(Vcmax=mean(Vcmax),Jmax=mean(Jmax),Rd=mean(Rd))%>%
  ggplot()+
  geom_boxplot(aes(x=Frond,y=Vcmax,fill=as.factor(Frond)))+
  scale_fill_viridis_d(name='Rank')+
  # facet_wrap(Progeny~Tree)+
  myTheme

###value of Farquhar parameters per tree & rank
ResSampleC%>%
  group_by(Progeny,Tree,Frond)%>%
  summarize(Vcmax=mean(Vcmax),Jmax=mean(Jmax),Rd=mean(Rd))%>%
  ggplot()+
  geom_col(aes(x=Frond,y=Vcmax,fill=as.factor(Frond)))+
  scale_fill_viridis_d(name='Rank')+
  facet_wrap(Progeny~Tree)+
  myTheme

ResSampleC%>%
  group_by(Progeny,Tree,Frond)%>%
  summarize(Vcmax=mean(Vcmax),Jmax=mean(Jmax),Rd=mean(Rd))%>%
  ggplot()+
  geom_col(aes(x=Frond,y=Jmax,fill=as.factor(Frond)))+
  scale_fill_viridis_d(name='Rank')+
  facet_wrap(Progeny~Tree)+
  myTheme

ResSampleC%>%
  group_by(Progeny,Tree,Frond)%>%
  summarize(Vcmax=mean(Vcmax),Jmax=mean(Jmax),Rd=mean(Rd))%>%
  ggplot()+
  geom_col(aes(x=Frond,y=Rd,fill=as.factor(Frond)))+
  scale_fill_viridis_d(name='Rank')+
  facet_wrap(Tree~Progeny)+
  myTheme


# Fitting PFD curves  -----------------------------------------------------
###fit per curve

ResSamplePFD=NULL


for(i in 1:length(unique(datPFD$Sample))){
  print(i)
  sub=datPFD%>% 
    filter(Sample==unique(datPFD$Sample)[i])%>%
    arrange(PFD)
  
  ####!!! data need to be updated each time
  res=Lightfitting(don=sub,plot=T)
  ResSamplePFD= rbind(ResSamplePFD,res$data)
}


# graphics---####

###curve fitting per sample
ResSamplePFD%>%
  ggplot(aes(x=PFD,y=A,group=Sample,col=as.factor(Frond)))+
  geom_point(aes(x=PFD,y=A,group=Sample,col=as.factor(Frond)))+
  geom_line(aes(x=PFD,y=A_sim,group=Sample,col=as.factor(Frond)))+
  scale_color_viridis_d(name='Rank')+
  facet_wrap(Progeny~Tree)+
  myTheme

###value of Farquhar parameters per tree & rank
ResSamplePFD%>%
  group_by(Progeny,Tree,Frond)%>%
  summarize(Jmax=mean(Jmax),Rd=mean(Rd),Sigma=mean(sigma),Theta=mean(theta))%>%
  ggplot()+
  geom_col(aes(x=Frond,y=Jmax,fill=as.factor(Frond)))+
  scale_fill_viridis_d(name='Rank')+
  facet_wrap(Tree~Progeny)+
  myTheme

ResSamplePFD%>%
  group_by(Progeny,Tree,Frond)%>%
  summarize(Jmax=mean(Jmax),Rd=mean(Rd),Sigma=mean(sigma),Theta=mean(theta))%>%
  ggplot()+
  geom_col(aes(x=Frond,y=Theta,fill=as.factor(Frond)))+
  scale_fill_viridis_d(name='Rank')+
  facet_wrap(Tree~Progeny)+
  myTheme

ResSamplePFD%>%
  group_by(Progeny,Tree,Frond)%>%
  summarize(Jmax=mean(Jmax),Rd=mean(Rd),Sigma=mean(sigma),Theta=mean(theta))%>%
  ggplot()+
  geom_col(aes(x=Frond,y=Sigma,fill=as.factor(Frond)))+
  scale_fill_viridis_d(name='Rank')+
  facet_wrap(Tree~Progeny)+
  myTheme
