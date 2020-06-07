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

ctrlCtrl <- spData$`TRUE`[,c("trialN", "nonSig","riskCtrl")]
milToPercent <- 1e-4 # div 1 million but as percentage
ctrlCtrl$trueNeg <- ctrlCtrl$nonSig * milToPercent
ctrlCtrl$falsePos <- 100 - ctrlCtrl$trueNeg
ctrlCtrl$riskCtrl <- toLookupPercent(ctrlCtrl$riskCtrl)
ctrlCtrl$nonSig <- NULL


ctrlRx <- spData$`FALSE`[,c("trialN", "nonSig", "RxBenefit", "riskCtrl", "riskRx")]
ctrlRx$falseNeg <- ctrlRx$nonSig * milToPercent
ctrlRx$truePos <- 100 - ctrlRx$falseNeg
ctrlRx$trueTruePos <- ctrlRx$RxBenefit * milToPercent
ctrlRx$riskCtrl <- toLookupPercent(ctrlRx$riskCtrl)
ctrlRx$riskRx <- toLookupPercent(ctrlRx$riskRx)
ctrlRx$nonSig <- NULL

data <- NULL
spData <- NULL

mmRr <- c(min(ctrlRx$riskRx), max(ctrlCtrl$riskCtrl))
mmN <- c(min(ctrlCtrl$trialN), max(ctrlCtrl$trialN))

displayVariables <- c("false positive = p < 0.05 despite control and treatment equal risk" = "falsePos",
                      "True Positive = p < 0.05 when treatment works, regardless of direction of risk demonstrated" = "truePos",
                      "True True Positive = p < 0.05 and also treatment risk correctly found lower than control risk" = "trueTruePos",
                      "True Negative (1 - false pos) only applies to equal risk in control and treatment" = "trueNeg",
                      "False Negative (1 - true pos)" = "falseNeg",
                      "Power Calculation (Binomial)" = "pwr",
                      "Power Calculation (Chi-squared)" = "pwrChi")
# Define UI
ui <- fluidPage(

    # Application title
    titlePanel("Monte Carlo Simulation Data"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("nRng",
                        "Min to Max Trial Participants(N):",
                        min = 60,
                        max = 6000,
                        value = c(200, 400)),
            sliderInput("riskRng",
                        "% Risk in treatment (lower slider) and control (upper slider) groups:",
                        min = mmRr[1],
                        max = mmRr[2],
                        value = c(25,30),
                        step = 1),
            checkboxGroupInput("show", "Variables to show:",
                               displayVariables,
                               selected = c("truePos")),
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
        for (pwr in intersect(c("pwrChi","pwr"),input$show)) {
            # h <- 2 * asin(sqrt(input$riskRng[1] / 100)) - 2 * asin(sqrt(input$riskRng[2] / 100))
            # same as above
            h <- ES.h(input$riskRng[2] / 100, input$riskRng[1] / 100)
            w <- ES.w1(input$riskRng[2] / 100, input$riskRng[1] / 100)
            i <- 1
            if (pwr == "pwr"){
                getPwr <- function(n) {
                    return(unname(pwr.2p2n.test(h = h , n1 = n/2, n2 = n/2, sig.level = 0.05, alternative = "two.sided")[["power"]])*100)
                }
                returnVar$pwr <- sapply(returnVar[,1], getPwr)
            } else {
                getPwr <- function(n) {
                    return(unname(pwr.chisq.test(N = n, w = w, df = 1, sig.level = 0.05)[["power"]])*100)
                }
                returnVar$pwrChi <- sapply(returnVar[,1], getPwr)
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
                            & ctrlRx$riskCtrl == input$riskRng[2],c("riskRx", "truePos", "trialN")])
    })
    
    output$againstN <- renderPlot({
        nDta <- reactive_XN()
        nDta$time <- factor(nDta$time, levels = displayVariables)
        colScheme <- setNames( c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d')
                  , levels(nDta$time))
        ggplot(nDta, aes(x=trialN, y=value, group=time)) +
            geom_line(size=1, aes(col=time))+
            # geom_point(aes(col=time))+
            scale_color_manual(values = colScheme, aesthetics = "colour")+
            labs(y="Percentage (%)", x = "Trial Participants (N)", colour = "Variable:",
                 caption=paste("Baseline (control) risk ",input$riskRng[2], "% vs ",input$riskRng[1],"%"))
    })
    output$againstXRisk <- renderPlot({
        ggplot(data = reactive_XRisk(), aes(x=riskRx, y=truePos, group=trialN)) +
            geom_line(aes(color=trialN))+
            geom_point(aes(color=trialN))+
            labs(y="Percentage true positives(%)", x = "Intervention Group Risk", colour = "N:",
                 caption=paste("as compared to a baseline (control) risk ",input$riskRng[2], "%"))
    })
}
shinyApp(ui=ui, server=server)
