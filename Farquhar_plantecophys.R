library(plantecophys)


# A-Ci curves
datC= 
  data.table::fread("Data/CO2curveEX.csv", data.table = FALSE)%>%
  mutate(Sample=paste(Progeny,Tree,Frond,sep='_'))%>%
  rename(Photo= A, PARi= PFD)

datC%>%
  ggplot(aes(x=Ci,y=Photo))+
  geom_point()
             
datC%>%
  ggplot(aes(x=Ci,y=A,group=Sample,col=as.factor(Frond)))+
  geom_point(aes(x=Ci,y=A,group=Sample,col=as.factor(Frond)))+
  geom_line(aes(x=Ci,y=A_sim,group=Sample,col=as.factor(Frond)))+
  scale_color_viridis_d(name='Rank')+
  facet_wrap(Progeny~Tree)+
  myTheme

# Fit at frond scale for each palm tree: 
fits= 
  datC%>%
  # filter(Sample!="DA8_100_5_9"&Sample!="DA8_101_2_25"&Sample!="DA8_99_3_4")%>%
  fitacis("Sample",id = "Tree") # Nb id is used to colour the plots
fits
coef(fits)
plot(fits, how="manyplots")
plot(fits, how="oneplot")

par(mar= c(5, 4.5, 4, 2) + 0.1)
plot(fits[[5]])

# Get all RMSEs:
rmses <- sapply(fits, "[[", "RMSE")
plot(rmses, type='h', ylab="RMSE", xlab="Curve nr")
# Plotting the fit with the worst RSME:
plot(fits[[which.max(rmses)]])

# Plotting Vcmax by treatment:
boxplot(Vcmax ~ Tree, data=coef(fits), ylim=c(0,130))

# Average parameter values: 
coef(fits)%>%
  summarise_if(is.numeric,mean)

# TPU limitation:
fit_TPU <- fitacis(datC, fitTPU=TRUE)
coef(fit_TPU)


###A-PFD curves
datPFD= data.table::fread("Data/lightcurveEX.csv", data.table = FALSE)%>%
  mutate(Sample= paste(Progeny,Tree,Frond,Date,sep='_'),
         Date= dmy(Date))