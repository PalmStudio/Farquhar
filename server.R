
# Load packages -----------------------------------------------------------

packs <- c('shiny','datasets',"lubridate", "stringr", "ggplot2",'dplyr','viridis','plotly')
InstIfNec<-function (pack) {
  if (!do.call(require,as.list(pack))) {
    do.call(install.packages,as.list(pack))  }
  do.call(require,as.list(pack)) }
lapply(packs, InstIfNec)

# Import functions --------------------------------------------------------

library(shiny)

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


shinyServer(function(input, output) {
  
  source('1-code/Functions_Farquhar_Fitting.R')
  
  filedata <- reactive({
    infile <- input$file1
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,header=TRUE,sep=";")
  })
  
  leafInput<- reactive({
    input$leaf
  })
  
  typeInput<- reactive({
    input$type
  })
  
  output$leaf <- renderUI({
    df <- filedata()
    if (is.null(df)) return(NULL)
    df=df%>%
      mutate(Sample=paste(Progeny,Tree,Frond,sep='_'))
    items=unique(df$Sample)
    names(items)=items
    selectInput("leaf","Select the leaf:",items,multiple=TRUE)
  })
  
  
  
  
  # Dataset of sampled leaves with predictions of assimilation-----------------------------------------------------------------
  Dataset<-reactive({
    input$action
    isolate({
      df <- filedata()
      if (is.null(df)) return(NULL)
      
      leaf<-leafInput()
      if (is.null(leaf)) return(NULL)
      
      type<-typeInput()
      if (type=='CO2 curve'){
        print(paste(type))
        
        # Fitting CO2 curves per leaf -----------------------------------------------------
        
        resLeaf=NULL
        
        for(i in leaf){
          
          print(i)
          
          sub=df%>%
            mutate(Sample=paste(Progeny,Tree,Frond,sep='_'),Rank=as.factor(Frond))%>%
            dplyr::filter(Sample %in% i)%>%
            arrange(Ci)
          
          
          res=CO2fitting(don =sub, plot=F)
          resLeaf=rbind(resLeaf,res$data)
        }
        
      }
      
      if (type=='Light curve'){
        print(paste(type))
        
        # Fitting PFD curves  -----------------------------------------------------
        
        resLeaf=NULL
        
        for(i in leaf){
          
          print(i)
          
          sub=df%>%
            mutate(Sample=paste(Progeny,Tree,Frond,sep='_'),Rank=as.factor(Frond))%>%
            dplyr::filter(Sample %in% i)%>%
            arrange(PFD)
          
          
          res=Lightfitting(don =sub, plot=F)
          resLeaf=rbind(resLeaf,res$data)
        }
        
      }
      
      Dataset=resLeaf
      return(Dataset)
      
    })
  })
  
  
  # Parameters set ----------------------------------------------------------
  Paramset<-reactive({
    input$action
    isolate({
      dat <- Dataset()
      type<-typeInput()
      
      if (is.null(dat)) return(NULL)
      
      if (type=='CO2 curve'){
        Paramset=dat%>%
          dplyr::group_by(Progeny,Tree,Rank)%>%
          summarize(Vcmax=mean(Vcmax),Jmax=mean(Jmax), Rd=mean(Rd),RMSE=mean(RMSE),RRMSE=mean(RRMSE))
      }
      
      if (type=='Light curve'){
        Paramset=dat%>%
          dplyr::group_by(Progeny,Tree,Rank)%>%
          summarize(Jmax=mean(Jmax),Theta=mean(theta), Rd=mean(Rd),RMSE=mean(RMSE),RRMSE=mean(RRMSE))
      }
      
      return(Paramset)
    })
  })
  
  
  
  # Plot  -------------------------------------------------------------------
  
  output$plot <- renderPlotly({
    input$action
    isolate({   
      df <- filedata()
      if (is.null(df)) return(NULL)
      
      leaf<-leafInput()
      if (is.null(leaf)) return(NULL)
      
      dat <- Dataset()
      if (is.null(dat)) return(NULL)
      
      type<-typeInput()
      
      if (type=='CO2 curve'){
        
        graph=dat%>%
          ggplot(aes(x=Ci,y=A,col=Sample))+
          geom_point(aes(x=Ci,y=A,col=Sample))+
          geom_line(aes(x=Ci,y=A_sim,col=Sample))+
          scale_color_viridis_d(name='')+
          ylab('A (micromol m-2 s-1)')+
          xlab('Ci (ppm)')
        
      }
      
      if (type=='Light curve'){
        
        graph=dat%>%
          ggplot(aes(x=PFD,y=A,col=Sample))+
          geom_point(aes(x=PFD,y=A,col=Sample))+
          geom_line(aes(x=PFD,y=A_sim,col=Sample))+
          scale_color_viridis_d(name='')+
          ylab('A (micromol per m2 per s))')+
          xlab('PFD (micromol m-2 s-1)')
      }
      print(graph)
      
    })   
  })
  
  
  
  # Table of parameters-------------------------------------------------------------------
  output$tab <- renderTable({Paramset()
  })
  
  
  # export parameters -------------------------------------------------------
  output$downloadData <- downloadHandler(
    filename = function() {
      paste('Parameters',input$file1$name, sep = "")
    },
    content = function(file) {
      write.csv(Paramset(), file, row.names = FALSE)
    }
  )
  
  
  output$ui.action <- renderUI({
    if (is.null(input$file1)) return()
    actionButton("action", "plot and fit")
  })
  
}) 
