#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(pwr)
library(ggplot2)

toLookupPercent <- function(v) {
    as.integer(100 * (v + .Machine$double.eps))
}

data <- readRDS("millionRuns.RDS")
data <- data[order(data$trialN),]
spData <- split(data, data$riskCtrl==data$riskRx)

ctrlCtrl <- subset(spData$`TRUE`, select=c("trialN", "riskCtrl","nonSig"))
convertDta <- 1e-4 # div 1 million but as percentage
ctrlCtrl$falsePos <- 100 - ctrlCtrl$nonSig * convertDta
ctrlCtrl$riskCtrl <- toLookupPercent(ctrlCtrl$riskCtrl)
ctrlCtrl$nonSig <- NULL


ctrlRx <- subset(spData$`FALSE`, select=-c(nonSig))
ctrlRx$RxHarm <- ctrlRx$RxHarm * convertDta
ctrlRx$RxBenefit <- ctrlRx$RxBenefit * convertDta
ctrlRx$riskCtrl <- toLookupPercent(ctrlRx$riskCtrl)
ctrlRx$riskRx <- toLookupPercent(ctrlRx$riskRx)

data <- NULL
spData <- NULL

mmRr <- c(min(ctrlRx$riskRx), max(ctrlCtrl$riskCtrl))
mmN <- c(min(ctrlCtrl$trialN), max(ctrlCtrl$trialN))


# Define UI
ui <- fluidPage(

    # Application title
    titlePanel("Monte Carlo Simulation Data"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("nRng",
                        "Min to Max Trial Participants(N):",
                        min = mmN[1],
                        max = mmN[2],
                        value = c(200, 400)),
            sliderInput("riskRng",
                        "% Risk in treatment (lower) and control (upper) group:",
                        min = mmRr[1],
                        max = mmRr[2],
                        value = c(25,30),
                        step = 1),
            checkboxGroupInput("show", "Variables to show:",
                               c("False Positive" = "falsePos",
                                 "True Positive" = "RxBenefit",
                                 "Incorrect harm" = "RxHarm",
                                 "Power Calculation" = "pwr"),
                               selected = c("RxBenefit")),
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput(outputId = "againstN", width="100%"),
           plotOutput(outputId = "againstXRisk", width="100%")
        )
    )
)

server <- function(input, output){
    reactive_XN <- reactive({
        ctrlNames <- intersect(c("trialN", input$show), names(ctrlCtrl))
        intervNames <- intersect(input$show, names(ctrlRx))
        returnVar <- data.frame(
            ctrlCtrl[input$nRng[1] <= ctrlCtrl$trialN & ctrlCtrl$trialN <= input$nRng[2]
                     & ctrlCtrl$riskCtrl == input$riskRng[2], ctrlNames],
            ctrlRx[input$nRng[1] <= ctrlRx$trialN & ctrlRx$trialN <= input$nRng[2] 
                   & ctrlRx$riskCtrl == input$riskRng[2] & ctrlRx$riskRx == input$riskRng[1], intervNames])
        colnames(returnVar) <- c(ctrlNames, intervNames)
        if ("pwr" %in% input$show) {
            h <- 2 * asin(sqrt(input$riskRng[1] / 100)) - 2 * asin(sqrt(input$riskRng[2] / 100))
            i <- 1
            returnVar$pwr <- rep(0,nrow(returnVar))
            for (n in returnVar[,1]) {
                returnVar$pwr[i] <- unname(pwr.p.test(h=h,n=n,sig.level=0.05,alternative="two.sided")[["power"]]) * 100
                i <- i + 1
            }
        }
        reshapeNm <- colnames(returnVar)[-1]

        return(reshape(data=returnVar, idvar="trialN",
                                 varying = reshapeNm,
                                 v.names=c("value"),
                                 times=reshapeNm,
                                 direction="long"))
    })
    
    reactive_XRisk <- reactive({
        return(ctrlRx[input$nRng[1] <= ctrlRx$trialN & ctrlRx$trialN <= input$nRng[2] 
                            & input$riskRng[1] <= ctrlRx$riskRx & ctrlRx$riskRx < input$riskRng[2]
                            & ctrlRx$riskCtrl == input$riskRng[2],c("riskRx", "RxBenefit", "trialN")])
    })
    
    output$againstN <- renderPlot({
        ggplot(data = reactive_XN(), aes(x=trialN, y=value, group=time)) +
            geom_line(aes(color=time))+
            geom_point(aes(color=time))+
            labs(y="Percentage (%)", x = "Trial Participants (N)", colour = "Variable:",
                 caption=paste("Baseline (control) risk ",input$riskRng[2], "% vs ",input$riskRng[1],"%"))
    })
    output$againstXRisk <- renderPlot({
        ggplot(data = reactive_XRisk(), aes(x=riskRx, y=RxBenefit, group=trialN)) +
            geom_line(aes(color=trialN))+
            geom_point(aes(color=trialN))+
            labs(y="Percentage (%)", x = "Intervention Group Risk", colour = "N:",
                 caption=paste("Baseline (control) risk ",input$riskRng[2], "%"))
    })
}
shinyApp(ui=ui, server=server)
