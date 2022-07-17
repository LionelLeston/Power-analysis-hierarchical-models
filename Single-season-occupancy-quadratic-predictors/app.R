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
library(scales)
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
    titlePanel("Simulate Single Season Occupancy and Detection: Quadratic Effects"),

    
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
                      "Habitat Effect (Linear):",
                      min = -4.5,
                      max = 4.5,
                      value = 1.05),
          sliderInput("hab2",
                      "Habitat Effect (Quadratic):",
                      min = -4.5,
                      max = 4.5,
                      value = 1.05),
          sliderInput("det0",
                      "Det. Intercept:",
                      min = -4,
                      max = 1,
                      value = 0.5),
          sliderInput("det1",
                      "Detection Effect (Linear):",
                      min = -4.5,
                      max = 4.5,
                      value = 1.05),
          sliderInput("det2",
                      "Detection Effect (Quadratic):",
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
    
    #make habitatvar2 a quadratic term of habitatvar1, assuming
    #that prior to scaling for model, habitatvar1 ranges from 0 to 1
    habitatvar1.unscaled<-rescale(occu_umf@siteCovs$habitatvar1, to=c(0, 1))
    #now create quadratic term after centering
    habitatvar2.unscaled<-(habitatvar1.unscaled-0.5)^2
    occu_umf@siteCovs$habitatvar2<-scale(habitatvar2.unscaled)
    
    #make detectionvar2 a quadratic term of detectionvar1, assuming
    #that prior to scaling for model, detectionvar1 ranges from 0 to 1
    detectionvar1.unscaled<-rescale(occu_umf@obsCovs$detectionvar1, to=c(0, 1))
    #now create quadratic term after centering
    detectionvar2.unscaled<-(detectionvar1.unscaled-0.5)^2
    occu_umf@obsCovs$detectionvar2<-scale(detectionvar2.unscaled)
    
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
    #make habitatvar2 a quadratic term of habitatvar1
    occu_umf@siteCovs$habitatvar2<-occu_umf@siteCovs$habitatvar1^2
    #make detectionvar2 a quadratic term of detectionvar1
    occu_umf@obsCovs$detectionvar2<-occu_umf@obsCovs$detectionvar1^2
    
    newdata<-data.frame(habitatvar1=occu_umf@siteCovs$habitatvar1,
                        habitatvar2=occu_umf@siteCovs$habitatvar1^2)
    m1<-predict(occmodel, type="state", newdata)
    newdata$preds<-m1$Predicted
    p1<-ggplot(newdata, aes(habitatvar1, preds))+
      geom_point()+
      geom_line(aes(habitatvar1, preds), col="blue")+
      xlab("Habitat Var.")+ylab("Prob.occ.")+my.theme
    
    newdata3<-data.frame(habitatvar1=mean(occu_umf@siteCovs$habitatvar1),
                         habitatvar2=mean(occu_umf@siteCovs$habitatvar1)^2,
                         detectionvar1=occu_umf@obsCovs$detectionvar1,
                         detectionvar2=occu_umf@obsCovs$detectionvar1^2)
    m3<-predict(occmodel, type="det", newdata3)
    newdata3$preds<-m3$Predicted
    p3<-ggplot(newdata3, aes(detectionvar1, preds))+
      geom_point()+
      geom_line(aes(detectionvar1, preds), col="blue")+
      xlab("Detection Var.")+ ylab("Prob.det.")+my.theme
    
    G1<-grid.arrange(p1, p3, nrow=1, ncol=2)
    
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
