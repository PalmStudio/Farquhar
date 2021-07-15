library(data.table)
library(tidyverse)
library(stringr)
require(ggpmisc)


# Import data ----------------------------------------
map= fread("A_mapping_PLPE.csv", data.table = F, fill = T)
map$Progeny=paste0('P',map$Progeny)
###relative position on rachis
map$PosRel=NA
map[map$Position=='a',]$PosRel=1
map[map$Position=='ab',]$PosRel=5/6
map[map$Position=='b',]$PosRel=2/3
map[map$Position=='1/4_bc',]$PosRel=0.5
map[map$Position=='1/2_bc',]$PosRel=1/3

map$SLA_N=map$SLA*map$N_cont

# relation beween N and photsynthesis -------------------------------------

map%>%
  ggplot(aes(x=N_cont,y=Photo,col=PosRel))+
  geom_point()+
  geom_smooth(method='lm',se=F,col=1)+
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc =0.2 , label.y.npc =1 ,
               formula = y~x, parse = TRUE, size = 5)+
  facet_grid(~Progeny)

map%>%
  ggplot(aes(x=N_cont,y=Photo,col=Rank))+
  geom_point()+
  geom_smooth(method='lm',se=F,col=1)+
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc =0.2 , label.y.npc =1 ,
               formula = y~x, parse = TRUE, size = 5)+
  facet_grid(~Progeny)

map%>%
  ggplot(aes(x=Rank,y=N_cont,col=Photo))+
  geom_point()+
  geom_smooth(method='lm',se=F,col=1)+
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc =0.2 , label.y.npc =1 ,
               formula = y~x, parse = TRUE, size = 5)+
  facet_wrap(~Palm)

map%>%
  ggplot(aes(x=N_cont,y=Photo,col=Rank))+
  geom_point()+
  geom_smooth(method='lm',se=F,col=1)+
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc =0.2 , label.y.npc =1 ,
               formula = y~x, parse = TRUE, size = 5)+
  facet_wrap(~Palm)

map%>%
  ggplot(aes(x=SLA,y=Photo,col=PosRel))+
  geom_point()+
  geom_smooth(method='lm',se=F,col=1)+
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc =0.2 , label.y.npc =1 ,
               formula = y~x, parse = TRUE, size = 5)+
  facet_wrap(~Palm)


map%>%
  ggplot(aes(x=Rank,y=SLA,col=PosRel))+
  geom_point()+
  geom_smooth(method='lm',se=F,col=1)+
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc =0.2 , label.y.npc =1 ,
               formula = y~x, parse = TRUE, size = 5)+
  facet_wrap(~Palm)

map%>%
  ggplot(aes(x=Rank,y=SLA,col=PosRel))+
  geom_point()+
  geom_smooth(method='lm',se=F,col=1)+
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc =0.2 , label.y.npc =1 ,
               formula = y~x, parse = TRUE, size = 5)+
  facet_wrap(~Palm)


# Statistics --------------------------------------------------------------
modelN=aov(data=map,N_cont~Progeny*Rank*PosRel)
summary(modelN)
plot(modelN)

modelSLA=aov(data=map,SLA~Progeny*Rank+PosRel+Rank:PosRel)
summary(modelSLA)
plot(modelSLA)
