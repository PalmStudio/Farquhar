library(data.table)
library(tidyverse)
library(stringr)

# Effect of VPD on photosynthesis ----------------------------------------
# VPD_df= fread("VPDeffect_20190220_R.csv", data.table = F, fill = T)
VPD_df= fread("Data/VPDdatabasedPLPE2018_R.csv", data.table = F, fill = T)

# CAREFULL !!!!
# Transpiration present two weird values, correcting them:
VPD_df$trans[VPD_df$trans>100 & !is.na(VPD_df$trans)]= 
  VPD_df$trans[VPD_df$trans>100 & !is.na(VPD_df$trans)]/100


# Computing other variables: 
VPD_df$Date= lubridate::dmy(VPD_df$Date)
VPD_df$hour= lubridate::hms(VPD_df$HHMMSS)@hour
VPD_df$Rank=as.numeric(str_remove(string = VPD_df$Frond,pattern = 'F'))

###relative position on rachis
VPD_df$PosRel=NA
VPD_df[VPD_df$Position=='A',]$PosRel=1
VPD_df[VPD_df$Position=='1/2_AB',]$PosRel=5/6
VPD_df[VPD_df$Position=='B',]$PosRel=2/3
VPD_df[VPD_df$Position=='1/4_BC',]$PosRel=0.5
VPD_df[VPD_df$Position=='1/2_BC',]$PosRel=1/3


  
str(VPD_df)

summary(VPD_df)

VPD_df= 
  VPD_df%>%
  filter(!is.na(Date))


# graphics ----------------------------------------------------------------
###effect on progeny and season 
VPD_df%>%
  ggplot(aes(y= Photo, x= VpdL))+
  geom_point(aes(color= Season))+
  facet_grid(~Progeny)

###effect of rank
VPD_df%>%
  ggplot(aes(y=Photo , x=VpdL ,col=as.factor(Rank)))+
  facet_grid(Progeny~Season)+
  geom_point()
  

###effect of temperature
VPD_df%>%
  ggplot(aes(y=Tair , x=VpdL ,col=Photo))+
  geom_point()+
  facet_grid(Progeny~Season)

###relationship between Temperature and rank
VPD_df%>%
  ggplot(aes(y=Tair , x=Rank ,col=VpdL))+
  geom_point()+
  facet_grid(Progeny~Season)


###check for confound effect of VPD and Rank du to time of measurments
VPD_df%>%
  ggplot(aes(y=VpdL , x=hour,col=as.factor(Rank)))+
  geom_point()+
  facet_grid(Progeny~Season)

###effect of position on photosynthesis
VPD_df$groups=paste(VPD_df$Season,VPD_df$Palm,VPD_df$Rank,sep='')

VPD_df%>%
  filter(!is.na(Photo))%>%
  ggplot(aes(y=Photo ,x=PosRel,group=groups,col=as.factor(Rank)))+
  geom_point()+
    geom_line()+
  facet_grid(Progeny~Season)

###effect of position on gs
VPD_df%>%
  filter(!is.na(gs))%>%
  ggplot(aes(y=gs ,x=PosRel,group=groups,col=as.factor(Rank)))+
  geom_point()+
  geom_line()+
  facet_grid(Progeny~Season)

VPD_df%>%
  filter(!is.na(gs))%>%
  ggplot(aes(y=gs ,x=VpdL,group=groups,col=as.factor(Rank)))+
  geom_point()+
  # geom_line()+
  facet_grid(Progeny~Season)

VPD_df%>%
  filter(!is.na(gs))%>%
  ggplot(aes(y=gs ,x=VpdL,col=Season))+
  geom_point(alpha=0.2)+
  geom_smooth()+
  # geom_line()+
  facet_grid(~Progeny)



# Statistical analysis ----------------------------------------------------

model_Photo=aov(data=VPD_df,Photo~Season+Progeny+Rank+PosRel+Season:Progeny+Season:PosRel)
summary(model_Photo)

plot(model_Photo)

model_gs=aov(data=VPD_df,gs~Season+Progeny+Rank+PosRel+Season:Rank+Season:PosRel)
summary(model_gs)

plot(model_gs)




# Fitting Medlyn et al. (2011) Gs model ------------------------------------

# Medlyn et al. (2011) model: 
# gs= g0 + (1 + g1/sqrt(VPD)) * (A/Ca)
# Correspondance of the variables: 
# A -> Photo
# gs -> gs
# VPD -> VpdL
# Ca -> 400

# Fitting g0 and g1 on all data:
VPD_df$gs_c= VPD_df$gs/1.57
# VPD_df$CO2S= 400
# VPD_df$RH= plantecophys::VPDtoRH(VPD = VPD_df$VpdL, TdegC = VPD_df$Tair)

plantecophys::fitBB(VPD_df, 
                    varnames = list(ALEAF = "Photo", GS = "gs_c", VPD = "VpdL",
                                    Ca = "CO2S", RH = "RH"),fitg0= TRUE)
Ca= 400 # /!\ Replace by the one measured !!!
Fit_g0_g1_c= nls(gs_c ~ g0 + (1 + g1/sqrt(VpdL)) * (Photo/Ca),
                 data= VPD_df, start = list(g0= 0.0033, g1= 12.5))
VPD_df$gs_medlyn_CO2= coef(Fit_g0_g1_c)[1] + (1 + coef(Fit_g0_g1_c)[2]/sqrt(VPD_df$VpdL)) * (VPD_df$Photo/400)
 
# And using gs for water:
VPD_df$gs_w= VPD_df$gs
Fit_g0_g1_w= nls(gs_w ~ g0 + (1 + g1/sqrt(VpdL)) * (Photo/400), data= VPD_df, start = list(g0= 0.0033, g1= 12.5))
VPD_df$gs_medlyn_H20= coef(Fit_g0_g1_w)[1] + (1 + coef(Fit_g0_g1_w)[2]/sqrt(VPD_df$VpdL)) * (VPD_df$Photo/400)


plot(VPD_df$gs_c,VPD_df$gs_medlyn_CO2)
abline(0,1)

VPD_df%>%
  ggplot(aes(x= gs_c, y= gs_medlyn_CO2, colour= Season))+
  facet_wrap(Progeny~.)+
  geom_point()+
  geom_abline(slope= 1 , intercept = 0)


# Effect of the hour: 
VPD_df%>%
  ggplot(aes(x= gs_c, y= gs_medlyn_CO2, colour= hour))+
  facet_grid(Progeny~Season)+
  geom_point()+
  geom_abline(slope= 1 , intercept = 0)

VPD_df%>%
  ggplot(aes(x= hour, colour= hour))+
  facet_grid(Progeny~Season)+
  geom_point(aes(y= gs_c, colour= "Observed"))+
  geom_point(aes(y= gs_medlyn_CO2, colour= "Medlyn"))

# Effect of the rank:
VPD_df%>%
  ggplot(aes(x= gs_c, y= gs_medlyn_CO2, colour= Rank))+
  facet_grid(Progeny~Season)+
  geom_point()+
  geom_abline(slope= 1 , intercept = 0)

# Effect of the position:
VPD_df%>%
  ggplot(aes(x= gs_c, y= gs_medlyn_CO2, colour= PosRel))+
  facet_grid(Progeny~Season)+
  geom_point()+
  geom_abline(slope= 1 , intercept = 0)


# Fitting according to groups:
VPD_df%>%
  dplyr::group_by(.data$Progeny)%>%
  dplyr::do(mod=nls(gs ~ g0 + (1 + g1/sqrt(VpdL)) * (Photo/400), data= ., start = list(g0= 0.0033, g1= 12.5)))%>%
  dplyr::mutate(g0= stats::coef(.data$mod)['g0'],
                g1= stats::coef(.data$mod)['g1'])

gs_season_progeny= 
  VPD_df%>%
  dplyr::group_by(.data$Season,.data$Progeny)%>%
  dplyr::do(mod=nls(gs ~ g0 + (1 + g1/sqrt(VpdL)) * (Photo/400), data= ., start = list(g0= 0.0033, g1= 12.5)))%>%
  dplyr::mutate(g0= stats::coef(.data$mod)['g0'],
                g1= stats::coef(.data$mod)['g1'])%>%
  dplyr::select(-mod)


VPD_df%>%
  left_join(gs_season_progeny)%>%
  mutate(gs_medlyn_CO2_adjust= g0 + (1 + g1/sqrt(.data$VpdL)) * (.data$Photo/400))%>%
  ggplot(aes(x= gs_c))+
  facet_grid(Progeny~Season)+
  geom_point(aes(y= gs_medlyn_CO2, colour= "Gs medlyn global"))+
  geom_point(aes(y= gs_medlyn_CO2_adjust, colour= "Gs medlyn adjusted"))+
  geom_abline(intercept = 0, slope = 1)


VPD_df%>%
  left_join(gs_season_progeny)%>%
  mutate(gs_medlyn_CO2_adjust= g0+ (1 + g1/sqrt(.data$VpdL)) * (.data$Photo/400))%>%
  ggplot(aes(x= Photo))+
  facet_grid(Progeny~Season)+
  geom_point(aes(y= gs_c, colour= "Observed"))+
  geom_point(aes(y= gs_medlyn_CO2, colour= "Gs medlyn global"))+
  geom_point(aes(y= gs_medlyn_CO2_adjust, colour= "Gs medlyn adjusted"))


# gs~VPD:
VPD_df%>%
  ggplot(aes(x= VpdL, colour= hour))+
  facet_grid(Progeny~Season)+
  geom_point(aes(y= gs, colour= "Observed"))+
  geom_point(aes(y= gs_medlyn_CO2, colour= "Medlyn"))



VPD_df%>%
  ggplot(aes(x= Photo))+
  facet_grid(Progeny~Season)+
  geom_point(aes(y= gs, colour= "Gs observed"))

VPD_df%>%
  ggplot(aes(x= VpdL))+
  facet_grid(Progeny~Season)+
  geom_point(aes(y= Photo, colour= hour))

VPD_df%>%
  ggplot(aes(x= VpdL))+
  facet_grid(Progeny~Season)+
  geom_point(aes(y= trans, colour= gs))

VPD_df%>%
  ggplot(aes(x= trans))+
  facet_grid(Progeny~Season)+
  geom_point(aes(y= gs, colour= hour))


VPD_df%>%
  ggplot(aes(x= VpdL, y= gs_c))+
  facet_grid(Progeny~Season)+
  geom_point()
