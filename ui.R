# Load packages -----------------------------------------------------------

packs <- c('shiny','datasets',"lubridate",'plotly')
InstIfNec<-function (pack) {
  if (!do.call(require,as.list(pack))) {
    do.call(install.packages,as.list(pack))  }
  do.call(require,as.list(pack)) }
lapply(packs, InstIfNec)

library(shiny)

shinyUI(fluidPage(
  fluidRow(
    column(12,wellPanel(
      titlePanel("Fitting Farquhar parameters")
    ))
  ),
  fluidRow(
    column(5, wellPanel(
      fileInput('file1', 'Choose CSV File',
                accept=c('text/csv', 
                         'text/comma-separated-values,text/plain', 
                         '.csv')),
      selectInput("type", "Type of response curve:", choices = c("CO2 curve", "Light curve")),
      tags$hr(),
      uiOutput("leaf"),
      tags$hr(),
      uiOutput('ui.action'))
    ),
    column(7, wellPanel(
      tableOutput("tab"),
      downloadButton("downloadData", "Download")
    ))
  ),
  fluidRow(
    column(12,wellPanel(
      plotlyOutput('plot')
    ))
  )
))

# shinyUI(
#   fluidPage(
#     fluidRow(
#       column(12,wellPanel(
#         titlePanel("Fitting Farquhar parameters")
#       ))
#     ),
#     fluidRow(
#       column(6, wellPanel(
#         selectInput("type", "Choose a type curve:", choices = c("CO2 curve", "Light curve")),
#         fileInput('dataInput', 'Import the data file (.csv)',
#                   accept=c('text/csv',
#                            'text/comma-separated-values,text/plain',
#                            '.csv'))
#       )),
#       column(6, wellPanel(
#       tableOutput("contents"),
#       tableOutput("tab")
#       ))
#     ),
#     fluidRow(
#       column(12,wellPanel(
#         plotOutput("graph1")
#       ))
#     )
#   )
# )
