#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(dplyr)
library(unmarked)
library(knitr)
library(ggplot2)
library(grid)
library(gridExtra)
my.theme <- theme_classic() +
  theme(text=element_text(size=20, family="Arial"),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(margin=margin(20,0,0,0)),
        axis.title.y=element_text(margin=margin(0,20,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1))

##################################################################
# Run a null-occupancy null-detection single-season model on 
# simulated point count data with a habitat predictor and detection
# predictor.
##################################################################


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Calculate Single Season Occupancy and Detection Probability"),

    
    ###Calculate power to detect a proportional change in occupancy between two distributions.
    

    sidebarLayout(
        sidebarPanel(
          sliderInput("S",
                      "#Sites (First Season):",
                      min = 10,
                      max = 1000,
                      value = 100),
          sliderInput("K",
                      "#Visits/Site:",
                      min = 2,
                      max = 20,
                      value = 4),
          sliderInput("hab0",
                      "Occ. Intercept:",
                      min = -4,
                      max = -2,
                      value = -3),
          sliderInput("hab1",
                      "Habitat Effect (1):",
                      min = -4.5,
                      max = 4.5,
                      value = 1.05),
          sliderInput("hab2",
                      "Habitat Effect (2):",
                      min = -4.5,
                      max = 4.5,
                      value = 1.05),
          sliderInput("det0",
                      "Detect. Intercept:",
                      min = -4,
                      max = 1,
                      value = 0.5),
          sliderInput("det1",
                      "Detection Effect (1):",
                      min = -4.5,
                      max = 4.5,
                      value = 1.05),
          sliderInput("det2",
                      "Detection Effect (2):",
                      min = -4.5,
                      max = 4.5,
                      value = -1.05),width=3),

        # Show a plot of the generated distribution
        mainPanel(
          fluidRow(
            splitLayout(cellWidths = c("33%", "67%"), tableOutput("modelResults"), plotOutput("Predictions"))
          ),width=9
        )
    )
)

# Define server logic required to generate tables based on input files and rank thresholds used
server <- function(input, output) {

  nullmodeltable <-function(S,K,hab0,hab1,hab2,det0,det1,det2){
    # Simulate an occupancy dataset
    # Covariates to include in simulation
    forms <- list(state=~habitatvar1+habitatvar2, det=~detectionvar1+detectionvar2)
    
    # Covariate effects and intercept values
    coefs <- list(state=c(intercept=hab0, habitatvar1=hab1, habitatvar2=hab2), 
                  #probability of occupancy ~ 0 when habitatvar=0
                  det=c(intercept=det0, detectionvar1=det1, detectionvar2=det2))

    # Study design
    design <- list(M=S, J=K) # S sites, K occasions per site
    
    #create occupancy data frame
    set.seed(1234)#necessary, to make sure that same dataset is used for both outputs
    occu_umf <- unmarked::simulate("occu", nsim=1, formulas=forms, coefs=coefs, design=design)
    
    #run occupancy model
    occmodel<-occu(formula = ~ detectionvar1+detectionvar2 ~ habitatvar1+habitatvar2, data = occu_umf)
    Values<-data.frame(coef(occmodel))
    Coefficients<-c("Occ-int","Occ-hab1","Occ-hab2","Det-int","Det-hab1","Det-hab2")
    df<-data.frame(Coefficients,Values)
    colnames(df)<-c("Coefficients","Values")
    H<-glimpse(data.frame(df))
    return(H)
  }
  
  allPlots <-function(S,K,hab0,hab1,hab2,det0,det1,det2){
    # Simulate an occupancy dataset
    # Covariates to include in simulation
    forms <- list(state=~habitatvar1+habitatvar2, det=~detectionvar1+detectionvar2)
    
    # Covariate effects and intercept values
    coefs <- list(state=c(intercept=hab0, habitatvar1=hab1, habitatvar2=hab2), 
                  #probability of occupancy ~ 0 when habitatvar=0
                  det=c(intercept=det0, detectionvar1=det1, detectionvar2=det2))
    
    # Study design
    design <- list(M=S, J=K) # S sites, K occasions per site
    
    #create occupancy data frame
    set.seed(1234)#necessary, to make sure that same dataset is used for both outputs
    occu_umf <- unmarked::simulate("occu", nsim=1, formulas=forms, coefs=coefs, design=design)
    
    #run occupancy model
    occmodel<-occu(formula = ~ detectionvar1+detectionvar2 ~ habitatvar1+habitatvar2, data = occu_umf)

    newdata<-data.frame(habitatvar1=occu_umf@siteCovs$habitatvar1,
                        habitatvar2=mean(occu_umf@siteCovs$habitatvar2))
    m1<-predict(occmodel, type="state", newdata)
    newdata$preds<-m1$Predicted
    p1<-ggplot(newdata, aes(habitatvar1, preds))+
      geom_point()+
      geom_line(aes(habitatvar1, preds), col="blue")+
      xlab("hab1")+ylab("Prob.occ.")+my.theme
    
    newdata2<-data.frame(habitatvar1=mean(occu_umf@siteCovs$habitatvar1),
                         habitatvar2=occu_umf@siteCovs$habitatvar2)
    m2<-predict(occmodel, type="state", newdata2)
    newdata2$preds<-m2$Predicted
    p2<-ggplot(newdata2, aes(habitatvar2, preds))+
      geom_point()+
      geom_line(aes(habitatvar2, preds), col="blue")+
      xlab("hab2")+ ylab("Prob.occ.")+my.theme
    
    newdata3<-data.frame(habitatvar1=mean(occu_umf@siteCovs$habitatvar1),
                         habitatvar2=mean(occu_umf@siteCovs$habitatvar2),
                         detectionvar1=occu_umf@obsCovs$detectionvar1,
                         detectionvar2=mean(occu_umf@obsCovs$detectionvar2))
    m3<-predict(occmodel, type="det", newdata3)
    newdata3$preds<-m3$Predicted
    p3<-ggplot(newdata3, aes(detectionvar1, preds))+
      geom_point()+
      geom_line(aes(detectionvar1, preds), col="blue")+
      xlab("det1")+ ylab("Prob.det.")+my.theme
    
    
    newdata4<-data.frame(habitatvar1=mean(occu_umf@siteCovs$habitatvar1),
                         habitatvar2=mean(occu_umf@siteCovs$habitatvar2),
                         detectionvar1=mean(occu_umf@obsCovs$detectionvar1),
                         detectionvar2=occu_umf@obsCovs$detectionvar2)
    m4<-predict(occmodel, type="det", newdata4)
    newdata4$preds<-m4$Predicted
    p4<-ggplot(newdata4, aes(detectionvar2, preds))+
      geom_point()+
      geom_line(aes(detectionvar2, preds), col="blue")+
      xlab("det2")+ ylab("Prob.det.")+my.theme    
    G1<-grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2)
    
    I<-glimpse(G1)
    return(I)
  }
  
    output$modelResults<-renderTable({
      #create table showing mean probability of occupancy and detection per site
      glimpse(nullmodeltable(S=input$S,
                             K=input$K,
                             hab0=input$hab0,
                             hab1=input$hab1,
                             hab2=input$hab2,
                             det0=input$det0,
                             det1=input$det1,
                             det2=input$det2))
    })
    output$Predictions<-renderPlot({
      #create table showing mean probability of occupancy and detection per site
      glimpse(allPlots(S=input$S,
                             K=input$K,
                             hab0=input$hab0,
                             hab1=input$hab1,
                             hab2=input$hab2,
                             det0=input$det0,
                             det1=input$det1,
                             det2=input$det2))
    })

}

# Run the application 
shinyApp(ui = ui, server = server)
