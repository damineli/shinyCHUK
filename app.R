# Version 2 of the App Layout for the TipFinding Algorythm:

# Libraries used in the app: ----
library("devtools")
library("utils")
library("shiny")
library("fields")
library("colorRamps")
library("waveslim")
library("xlsx")
library("markdown")
library("dygraphs")
source("TipFindingApp_Funcs.R")


# Defining the app UI: ---------------------------------------------------
# JScode <-
#   "$(function() {
#     setTimeout(function(){
#       var vals = [0];
#       var powStart = -2;
#       var powStop = 2;
#       for (i = powStart; i >= powStop; i++) {
#         var val = Math.pow(10, i);
#         val = parseFloat(val.toFixed(8));
#         vals.push(val);
#       }
#       $('#pvalue').data('ionRangeSlider').update({'values':vals})
#     }, 5)})"
# Using navbarPage as the main Layout
ui <- navbarPage(title = "CHUKNORRIS 2.0 !!", id = "navbar.tabs",
  tabPanel(title = "'First Step'",
    fluidRow(
      column(4,
        fluidRow(
          radioButtons(inputId = "Decision_KymoTS", label = "Choose the type of data to analyze:",
            choiceNames = c("Kymographs", "Time Series"), 
            choiceValues = c(1, 2))
        ),
        conditionalPanel(condition = "input.Decision_KymoTS == 1",
          radioButtons(inputId = "Decision_Kymo.type", label = "Choose option:",
            choiceNames = c("1 channel", "multi-channel (with same time 'range')"),
            choiceValues = c(1, 2))
        ),
        conditionalPanel(condition = "input.Decision_KymoTS == 2",
          radioButtons(inputId = "Decision_TS.type", label = "Choose option:",
            choiceNames = c("Single Time Series", "Multiple time Series (simultaneous series)"),
            choiceValues = c(1, 2))
        ),
        offset = 2),
      column(4,
        h5("For now, we're only working with 1 channel Kymographs (the other features are still under construction)"),
        br(),
        br(),
        actionButton(inputId = "Decision_Go.button", label = "Open analysis tab"))
    )),
  # tags$head(tags$script(HTML(JScode))),
  
  # First Panel ("HOME"): ----
  tabPanel(title = "'Home'",
    fluidRow(  ## Start of fluidRow 1.1
      column(6,
        includeMarkdown("www/AppDescription_v1.Rmd")
      ),
      column(6,
        img(src = "DefaultKymoInput.png", height = 400, width = 550), # file source
        br(),
        br(),
        tags$h6(strong("Alogrythm designed by ",
          a(href = "https://www.researchgate.net/profile/Daniel_Damineli",
            "Daniel S C Damineli"),
          ",",
          a(href = "https://academic.oup.com/jxb/article/68/12/3267/3091950",
            " \"Oscillatory signatures underlie growth regimes in Arabidopsis pollen tubes: 
            computational methods to estimate tip location, periodicity, and synchronization 
            in growing cells\""))
        ),
        br(),
        tags$h6(strong("Interface designed by Francisco F. A. Neves 
          (Undergraduate Student under the supervision of",
          a(href = "https://www.researchgate.net/profile/Daniel_Damineli",
            "Daniel Damineli"),
          "and ", 
          a(href = "https://cbmg.umd.edu/faculty/josefeijo",
            "Jose Feijo"), 
          " at ",
          a(href = "https://cbmg.umd.edu", 
            "Cell Biology and Molecular Genetics department at University of Maryland")))
        )
  )## End fluidRow 1.1
    ),# End of First Panel ("HOME")
  # Second Panel ("Data Upload"): ----
  tabPanel(title = "'Data Upload'",
    sidebarLayout( ## Start of sidebarLayout 2.1
      # sidebarPanel for inputs in Tab2
      sidebarPanel( 
        fileInput(inputId = "user.data", label = "Upload Data:", 
          buttonLabel = "Choose file (.csv/.txt)", 
          accept = c(
            'text/csv',
            'text/comma-separated-values',
            'text/tab-separated-values',
            'text/plain',
            '.csv',
            '.tsv')), # upload data 
        hr(),
        ## fluidRow, to be able to create columns with initial parameters (inputs)
        checkboxInput(inputId = "changePixel", label = "Specify custom Pixel Size (default = 1 [pixel])", value = F),
        checkboxInput(inputId = "changeTime", label = "Specify custom Time Step (default = 1 [frame])", value = F),
        fluidRow(
          column(6,
            conditionalPanel(condition = "input.changePixel == true",
              numericInput(inputId = "pixel.sz", label = "Pixel Size:", 
                value = 1, step = 0.01, min = 0, max = 1)
            ),
            conditionalPanel(condition = "input.changeTime == true", 
              numericInput(inputId = "t.step", label = "Time Step:", 
                value = 1, step = 0.01, min = 0, max = 1)
            )),
          column(6,
            conditionalPanel(condition = "input.changePixel == true",
              selectInput(inputId = "u.pixel", label = "Pixel Unit:", 
                choices = c("um", "nm"))  # "pixel" is not an option in this case (it'll be put on the code in the server)
            ),
            conditionalPanel(condition = "input.changeTime == true",
              selectInput(inputId = "u.time", label = "Time Unit:",
                choices = c("minute", "second"))  # "frame" is not an option in this case (it'll be put on the code in the server)
            ))
        ),
        numericInput(inputId = "redBias.tab2", label = "Red Bias ('Graph color range')",
          min = 0, max = 100, step = 0.01, value = 1, width = '60%'), # Might change to a slider later on
        hr(),
        h5("Please confirm this is the data you want analyzed:"),
        actionButton(inputId = "confirm.user.data", label = "Confirm Data"), # user confirms data
        h6("Check the next tab ('Tip Detection') for continuing with the Analysis")#,
        # radioButtons(inputId = "analysisDecision_start", label = "Choose the analysis:", 
        #              choiceNames = c("Ratiometrics", "Multiple-Channels",  "Time Series Analysis"), 
        #              choiceValues = c("ratiometric", "multi.channel", "time.series")),
        # actionButton(inputId = "Decision_Go.button", label = "Open analysis tab")
      ),
      # mainPanel for outputs in Tab2
      mainPanel(
        h4("How the Data should look like:"),
        img(src = "DefaultKymoInput.png", height = 400, width = 550),
        hr(),
        h4("Your Data:"),
        plotOutput(outputId = "initial.plot_tab2", height = 500),
        fluidRow(
          column(3,
            actionButton(inputId = "row.invert", label = "Invert Y-axis"), 
            offset = 2),
          column(3,
            actionButton(inputId = "col.invert", label = "Invert X-axis"),
            offset = 3)
        ),
        fluidRow(
          column(3, 
            actionButton(inputId = "clock.rot", label = "Rotate"),
            offset = 5)
        ),
        br(),
        br()
      )
    )## End of sidebarLayout 2.1
  ),# End of Second Panel ("Seiing the Data")
  # Third Panel ("Tip Detection"): ----
  tabPanel(title = "'Tip Detection'",
    fluidRow(
      column(6, 
        actionButton(inputId = "next.tab34", label = "Show Tab 'Smooth Tip Location'"),
        offset = 3, align = "center")
    ),
    fluidRow(
      column(3,
        checkboxInput(inputId = "advSetting", label = "Advanced Settings", value = FALSE),
        conditionalPanel(condition = "input.advSetting == true",
          wellPanel(
            selectInput(inputId = "algorithmChoice", label = "Algorithm:", 
              choices = c("simple", "oldonline"), selected = "simple"),
            helpText("Option 'oldonline' is under work for now"),
            checkboxInput(inputId = "useCoarse", label = "Coarse Smoothing as ref.",
              value = FALSE),
            checkboxInput(inputId = "useSmooth", label = "Smooth Kymograph", #strong("Smooth Kymograph"),
              value = TRUE),
            checkboxInput(inputId = "rmvMinimum", label = "Remove Minimum", 
              value = TRUE)
          )))
    ),
    fluidRow( 
      column(3,
        h4("Controls:"),
        wellPanel(
          # checkboxInput for useSmooth is in the advanced settings box
          conditionalPanel(condition = "input.useSmooth == true",
            fluidRow(
              column(6,
                numericInput(inputId = "smooth.wnd.sz", label = "Window Size for Smoothing:",
                  value = 7, min = 0, max = 20)),
              column(6,
                selectInput(inputId = "smooth.dg", label = "Degree used for Smoothing:",
                  choices = c(1, 2, 0)))
            )
          ),
          conditionalPanel(condition = "input.useCoarse == true",
            fluidRow(
              column(6,
                numericInput(inputId = "coarse.wnd.sz", label = "Window Size for Coarse Smoothing:",
                  value = 50, min = 0, max = 100)),
              column(6,
                selectInput(inputId = "coarse.dg", label = "Degree used for Coarse Smoothing:",
                  choices = c(2, 1, 0)))
            )
          )
        ),
        wellPanel(
          selectInput(inputId = "tip.estimate", label = "Tip Estimate",
            choices = c("subpixel", "pixel", "pixel.min", "pixel.max", "global.max", "max.peak")),
          fluidRow(
            column(6,
              checkboxInput(inputId = "restrict.X", label = "Restrict 'X'", value = FALSE)),
            column(6,
              conditionalPanel(condition = "input.algorithmChoice == 'oldonline'",
                checkboxInput(inputId = "fix.slope", label = "Fix Slope", value = FALSE)
              ))
          ),
          numericInput(inputId = "redBias.tab3", label = "Red Bias ('Graph color range')",
            min = 0.01, max = 100, step = 0.01, value = 1, width = '80%') # Might change to a slider later on
        )
      ),
      column(9,
        plotOutput(outputId = "main.plot_tab3", height = 600))
    ), 
    hr(), # can accept other colors for the line ( tags$hr(style="border-color: gray;") )
    fluidRow( 
      column(9,
        plotOutput(outputId = "slice.plot_tab3")),
      column(3,
        sliderInput(inputId = "t.slice", label = "Time Slice:", 
          value = 100, min = 1, max = 500, step = 1),
        # numericInput(inputId = "t.slice", label = "Time Slice:",
        #              value = 100, min = 1, max = 50000, step = 1),  # change to slider - update depending on unit !!!
        fluidRow(
          column(7,
            numericInput(inputId = "n.pts", label = "Number of points for fit:", #check the label
              value = 5, min = 3, max  = 20, step = 1)),
          column(5,
            conditionalPanel(condition = "input.algorithmChoice == 'oldonline'",
              sliderInput(inputId = "mad.tol", label = "Tolerance used for MAD:",
                value = 3.5, min = 1, max = 10, step = 0.25)
            ))
        ),
        conditionalPanel(condition = "input.algorithmChoice == 'simple'",
          wellPanel(
            fluidRow(
              column(6, 
                numericInput(inputId = "chunk.size", label = "Chunk Size:",
                  value = 7, min = 7, max = 500, step = 1)),
              column(6,
                sliderInput(inputId = "fluo.frac", label = "Fluorescence Fraction (threshold):",
                  value = 0.25, min = 0, max = 1, step = 0.01))
            )
          ))
      )
    )
  ),# End of Third Tab ("TIP DETECTION")
  # Fourth Panel("SMOOTH TIP LOCATION"): ----
  tabPanel(title = "'Smooth Tip Location'",
    fluidRow(
      column(6, 
        actionButton(inputId = "next.tab45", label = "Show Tab 'Extract Fluorescence Series'"),
        offset = 3, align = "center")
    ),
    fluidRow(
      column(3,
        wellPanel(     # Check Inputs' labels!
          fluidRow(
            column(6,
              numericInput(inputId = "tip.spn", label = "Tip Smooth Window Size:",
                value = 7, min = 0, max = 20)),
            column(6,
              selectInput(inputId = "tip.spn.dg",label = "Tip Smooth Degree:",
                choices = c(1, 2, 0)))
          ),
          checkboxInput(inputId = "removeTipOut", label = "Remove Tip Outlier", value = FALSE),
          conditionalPanel(condition = "input.removeTipOut == true",
            fluidRow(
              column(12, 
                selectInput(inputId = "outlierAlgorithm", label = "Ourlier removal Algorithm:",
                  choices = c("maddiff"))
              )
            ),
            fluidRow(
              column(6,
                numericInput(inputId = "remove.tip.out.spn", label = "Span:",
                  value = 0.4, min = 0.1, max = 20, step = 0.1)),
              column(6,
                selectInput(inputId = "remove.tip.out.dg", label = "Degree:",
                  choices = c(1, 2, 0)))
            ),
            fluidRow(
              column(6,
                sliderInput(inputId = "tol.tab4", label = "Tolerance:", 
                  value = 10, min = 1, max = 20, step = 0.25)),
              column(6,
                sliderInput(inputId = "px.tol", label = "Pixel Tolerance:",
                  value = 0, min = 0, max = 20))
            )),
          numericInput(inputId = "redBias.tab4", label = "Red Bias ('Graph color range')",
            min = 0.01, max = 100, step = 0.01, value = 1, width = '80%'),
          checkboxInput(inputId = "rmv.legend", label = "Remove plot legends", value = FALSE)
        ),
        helpText("There is an interactive version of the growth rate plot below")),
      column(9,
        h5(strong("Don't flip out! We are flipping the axis to have time on the x-axis"), align = "center"),
        plotOutput(outputId = "main.plot_tab4", height = 800))
    ), 
    fluidRow(
      column(3,
        helpText("Click and drag to zoom in (double click to zoom back out).")),
      column(9,
        dygraphOutput(outputId = "dygraphTab4", width = "88%")#,
      ),#offset = 3),
      hr(),
      hr()
    )
  ),# End of Fourth Tab("SMOOTH TIP LOCATION")
  # Fifth Panel ("EXTRACT FLUORESCENCE SERIES"): ----
  tabPanel(title = "'Extract Fluorescence Series'",
    fluidRow(
      column(4, 
        actionButton(inputId = "next.tab56", label = "Show Tab 'Save'"),
        offset = 4, align = "center")#,
      # column(3, 
      #        actionButton(inputId = "tab.save", label = "Show 'Save' Tab"), 
      #        align = "center")
    ),
    fluidRow(
      column(3,
        wellPanel(
          sliderInput(inputId = "tip.margin", label = "Tip Margin:",
            value = 0, min = -100, max = 100),
          sliderInput(inputId = "avrg.width", label = "Average Width:",
            value = 5, min = 1, max = 15),
          checkboxInput(inputId = "useMedian.avrg", label = "Use Median for Average"),
          checkboxInput(inputId = "manualROI", label = "Manual Region of Interest", value = FALSE),
          conditionalPanel(condition = "input.manualROI == true",
            numericInput(inputId = "roi.px", label = "Region of Interest (Pixel):",
              min = 0.01, max = 200, value = 123)
          ),
          numericInput(inputId = "redBias.tab5", label = "Red Bias ('Graph color range')",
            min = 0.01, max = 100, step = 0.01, value = 0.75, width = '60%'),
          textInput(inputId = "y.lab.tab5", label = "Y-axis legend:", value = "Fluorescence (AU)")
        ),
        helpText("There is an interactive version of the fluorescence series plots below"),
        helpText("Click and drag to zoom in (double click to zoom back out) - for the interactive plots.")
      ),
      column(9,
        plotOutput(outputId = "main.plot_tab5", height = 900),
        dygraphOutput(outputId = "dygraph1Tab5", width = "88%"),
        hr(),
        dygraphOutput(outputId = "dygraph2Tab5", width = "88%")
      )
    )
    
  ),# End of Fifth Tab ("EXTARCT FLUORESCENCE SERIES")
  # Sixth Panel ("FILTER KYMOGRAPH"): ----
            # tabPanel(title = "'Filter Kymograph'",
            #   fluidRow(
            #     column(6,
            #       h5("Please decide whether or not to run this tab (if your file is too big, ~ 10MB, the app might crash):"),
            #       fluidRow(
            #         column(6,
            #           actionButton(inputId = "decision.tab6.run", label = "Run")),
            #         column(6,
            #           actionButton(inputId = "decision.tab6", label = "Don't run"))
            #       ),
            #       h6("If you choose to run the tab, it'll still take a few seconds for the graphs to show up (if necessary you can still choose to not run it)"),
            #       h6("If you chose to 'Don't run', but changed your mind, please restart the app (we're still working on this transition)"),
            #       h6("Go to the 'Save' tab to save the rest of your analysis"),
            #       offset = 3, align = "center")
            #   ),
            #   fluidRow(
            #     column(3,
            #       wellPanel(
            #         sliderInput(inputId = "period.slider", label = "Period Range:",
            #           min = 0, max = 800, value = c(0, 800), step = 0.01),
            #         sliderInput(inputId = "trimTime", label = "Trim Time (unit):",
            #           min = 0, max = 800, value = c(1, 654)),
            #         sliderInput(inputId = "trimLength", label = "Trim Length (unit):",
            #           min = 0, max = 876, value = c(1, 543)),
            #         numericInput(inputId = "tip.mrgn.filt", label = "Tip Margin Filter:",
            #           min = 0, max = 10, value = 0),
            #         fluidRow(
            #           column(6,
            #             checkboxInput(inputId = "manualROI_tab6", label = "Manual ROI", value = FALSE)),
            #           column(6,
            #             conditionalPanel(condition = "input.manualROI_tab6 == true",
            #               numericInput(inputId = "roi.px.tab6", label = "Region of Interest (Pixel):",
            #                 min = 0.01, max = 200, value = 123))
            #           )
            #         ),
            #         fluidRow(
            #           column(6,
            #             sliderInput(inputId = "avrg.width_tab6", label = "Average Width:",
            #               value = 5, min = 1, max = 15)),
            #           column(6,
            #             checkboxInput(inputId = "useMedian.avrg_tab6", label = "Use Median for Average"))
            #         ),
            #         fluidRow(
            #           column(6,
            #             numericInput(inputId = "redBias.tab6", label = "Red Bias ('Graph color range')",
            #               min = 0.01, max = 100, step = 0.01, value = 1.65)),
            #           column(6,
            #             numericInput(inputId = "redBias.tab6.trim", label = "Red Bias ('Graph color range') - Second Plot",
            #               min = 0.01, max = 100, step = 0.01, value = 1.65))
            #         ),
            #         textInput(inputId = "y.lab.tab6", label = "Y-axis legend:", value = "Fluorescence (AU)"),
            #         checkboxInput(inputId = "rmv.legend.tab6", label = "Remove plot legends", value = FALSE)#,
            #         #actionButton(inputId = "next.tab6save", label = "Show 'Save' Tab")
            #       )),
            #     column(9,
            #       plotOutput(outputId = "main.plot_tab6", height = 600),
            #       hr(),
            #       plotOutput(outputId = "plot2.tab6", height = 600))
            #   )
            # ),
  # 'Multi-Channel' and 'Time-Series' Tabs: ----
  tabPanel(title = "'Multi-Channel'",
    h4("Tab for the ", strong("Multi-Channel Analysis"), " - we're stil working on it...")
  ),
  tabPanel(title = "'Single Time-Series Analysis'",
    h4("Tab for ", strong("Single Time-Series Analysis"), " still under construction...")
  ),
  tabPanel(title = "'Multiple Time-Series Analysis'",
    h4("Tab for ", strong("Multiple Time-Series Analysis"), " still under construction...")
  ),
  # Last Panel ("SAVE"): ----
  tabPanel(title = "'Save'",
    fluidRow(
      column(3,
        fluidRow(
          radioButtons(inputId = "dataFileType", label = "Choose file type for the data:",
            choices = c("txt", "csv"))
        ),
        conditionalPanel(condition = "input.dataFileType == 'txt'",
          radioButtons(inputId = "txt.sep", label = "Separator for txt files:",
            choiceNames = c("space", "tab"), choiceValues = c(" ", "\t"))
        ),
        conditionalPanel(condition = "input.dataFileType == 'csv'",
          radioButtons(inputId = "csv.sep", label = "Separator for csv files:",
            choiceNames = c("comma", "semicolon"), choiceValues = c(",", ";"))
        ),
        offset = 3),
      column(3,
        textInput(inputId = "out.Data_name", label = "File name:", value = "TipFinding"),
        downloadButton(outputId = "out.Data", label = "Download Data from Analysis")
      )),
    br(),
    hr(),
    br(),
    fluidRow(
      column(3,
        fileInput(inputId = "BatchZip", label = "Zip with files for 'Batch Mode'", multiple = T, 
          accept = c(
            'text/csv',
            'text/comma-separated-values',
            'text/tab-separated-values',
            'text/plain',
            '.csv',
            '.tsv')),
        h5("For 'Batch Mode', the app will re-run the analysis with the same parameters that have been set previously,
          for each file selected, and return all the results of the analysis through a zip folder"),
        offset = 3), 
      column(3,
        downloadButton(outputId = "Batch.out", label = "Download Data from Batch Analysis"))
      )
    )
  )# End of UI




# Defining Server --------------------------------------------------------

# File size limit for uploading is 5MB, by default, but it can be raised using these functions:
# options(shiny.maxRequestSize = 30*1024^2) # Raising limit to 30MB
# options(shiny.maxRequestSize = 9*1024^2)  # Raising limit to 9MB
# options(shiny.maxRequestSize = 50*1024^2) # Here we're trying to raise limit to 50MB
options(shiny.maxRequestSize = 30*1024^2)

server <- function(input, output, clientData, session) {
  #options(shiny.reactlog = TRUE)
  Sys.setenv(R_ZIPCMD="/usr/bin/zip")
  
  # __ # __ # ----
  
  hideTab(inputId = "navbar.tabs", target = "'Home'")
  hideTab(inputId = "navbar.tabs", target = "'Data Upload'")
  hideTab(inputId = "navbar.tabs", target = "'Multi-Channel'")
  hideTab(inputId = "navbar.tabs", target = "'Single Time-Series Analysis'")
  hideTab(inputId = "navbar.tabs", target = "'Multiple Time-Series Analysis'")
  hideTab(inputId = "navbar.tabs", target = "'Tip Detection'")
  hideTab(inputId = "navbar.tabs", target = "'Smooth Tip Location'")
  hideTab(inputId = "navbar.tabs", target = "'Extract Fluorescence Series'")
  hideTab(inputId = "navbar.tabs", target = "'Save'")
  
  Anal.decision_KymoTS <- reactive({ input$Decision_KymoTS })
  Anal.decision_Kymo <- reactive({ input$Decision_Kymo.type })
  Anal.decision_TS <- reactive({ input$Decision_TS.type })
  
  observeEvent(eventExpr = {
    input$Decision_Go.button
  }, handlerExpr = {
    if(Anal.decision_KymoTS() == 1) {
      if(Anal.decision_Kymo() == 1) {
        showTab(inputId = "navbar.tabs", target = "'Home'")
        showTab(inputId = "navbar.tabs", target = "'Data Upload'")
        hideTab(inputId = "navbar.tabs", target = "'Multi-Channel'")
        hideTab(inputId = "navbar.tabs", target = "'Single Time-Series Analysis'")
        hideTab(inputId = "navbar.tabs", target = "'Multiple Time-Series Analysis'")
      }
      if(Anal.decision_Kymo() == 2) {
        showTab(inputId = "navbar.tabs", target = "'Multi-Channel'")
        hideTab(inputId = "navbar.tabs", target = "'Home'")
        hideTab(inputId = "navbar.tabs", target = "'Data Upload'")
        hideTab(inputId = "navbar.tabs", target = "'Single Time-Series Analysis'")
        hideTab(inputId = "navbar.tabs", target = "'Multiple Time-Series Analysis'")
      }
    }
    if(Anal.decision_KymoTS() == 2) {
      if(Anal.decision_TS() == 1) {
        showTab(inputId = "navbar.tabs", target = "'Single Time-Series Analysis'")
        hideTab(inputId = "navbar.tabs", target = "'Home'")
        hideTab(inputId = "navbar.tabs", target = "'Data Upload'")
        hideTab(inputId = "navbar.tabs", target = "'Multi-Channel'")
        hideTab(inputId = "navbar.tabs", target = "'Multiple Time-Series Analysis'")
      }
      if(Anal.decision_TS() == 2) {
        showTab(inputId = "navbar.tabs", target = "'Multiple Time-Series Analysis'")
        hideTab(inputId = "navbar.tabs", target = "'Home'")
        hideTab(inputId = "navbar.tabs", target = "'Data Upload'")
        hideTab(inputId = "navbar.tabs", target = "'Multi-Channel'")
        hideTab(inputId = "navbar.tabs", target = "'Single Time-Series Analysis'")
      }
    }
  })
  
  
  
  observeEvent(eventExpr = input$confirm.user.data, {
    showTab(inputId = "navbar.tabs", target = "'Tip Detection'")
  })
  
 
  observeEvent(eventExpr = input$next.tab34, {
    showTab(inputId = "navbar.tabs", target = "'Smooth Tip Location'")
  })
  
  observeEvent(eventExpr = input$next.tab45, {
    showTab(inputId = "navbar.tabs", target = "'Extract Fluorescence Series'")
  })
  
  observeEvent(eventExpr = input$next.tab56, {
    showTab(inputId = "navbar.tabs", target = "'Save'")
  })
        
  
  Anal.decision_start <- reactive({ input$analysisDecision_start })
  
              
  
  
  # Inputs and calls with User's data (Second Tab): ----
  ## input$file1 will be NULL initially. After the user selects
  ## and uploads a file, it will be a data frame with 'name',
  ## 'size', 'type', and 'datapath' columns. The 'datapath'
  ## column will contain the local filenames where the data can
  ## be found - but the file names are in the 'name' column!
  inFile <- reactive({ input$user.data  })
  fl.nm <- reactive({ inFile()$datapath })
  
  changePxl <- reactive({ input$changePixel })
  changeTm <- reactive({ input$changeTime })
  px.sz <- reactive({
    if(changePxl() == TRUE) { return(input$pixel.sz) }
    return(1) # default = 1pixel
  })
  px.unit <- reactive({ 
    if(changePxl() == TRUE) { return(input$u.pixel) }
    return("pixel")  # default = 1pixel
  })
  time.step <- reactive({ 
    if(changeTm() == TRUE) { return(input$t.step ) }
    return(1) # default = 1frame
  })
  time.unit <- reactive({ 
    if(changeTm() == TRUE) { return(input$u.time) }
    return("frame")
  })
  
  red.bias.tab2 <- reactive({ input$redBias.tab2})
  confirmed <- reactive({ input$confirm.user.data  })
  
  ## Calls:
  kymo.1 <- eventReactive(eventExpr = input$user.data,
    valueExpr = {
      #file.name <- fl.nm()
      ReadKymo(fl.nm())
    })
  clock <- reactive({ input$clock.rot })
  # counter <- reactive({ input$counter.rot })
  inv.col <- reactive({ input$col.invert })
  inv.row <- reactive({ input$row.invert })
  
  kymo_untreated <- reactive({
    kmo <- kymo.1()
    
    for (i in seq(1, clock(), length.out = clock())) {
      kmo <- clock.rotation(kmo)
    }
    # for (i in seq(1, counter(), length.out = counter())) {
    #   kmo <- counter.rotation(kmo)
    # }
    for (i in seq(1, inv.col(), length.out = inv.col())) {
      kmo <- invert.column(kmo)
    }
    for (i in seq(1, inv.row(), length.out = inv.row())) {
      kmo <- invert.row(kmo)
    }
    return(kmo)
  })
  
  kymo <- reactive({ TreatKymo(kymo_untreated()) })
  
  # Inputs and calls for Tip Detection (Third Tab): ----
  qntl <- 0.95
  fit.right <- FALSE
  
  use.smooth <- reactive({ input$useSmooth })
  kymo.span.n <- reactive({ input$smooth.wnd.sz })
  kymo.loess.dg <- reactive({ input$smooth.dg })
  
  use.coarse <- reactive({ input$useCoarse })
  coarse.kymo.span.n <- reactive({ input$coarse.wnd.sz })
  coarse.kymo.loess.dg <- reactive({ input$coarse.dg })
  
  n.pts <- reactive({ input$n.pts })
  mad.tol <- reactive({ input$mad.tol })
  fluo.frac <- reactive({ input$fluo.frac })
  min.chunk.size <- reactive({ input$chunk.size })
  
  fix.slope <- reactive({ input$fix.slope })
  
  algorithm <- reactive({ input$algorithmChoice })
  rmv.min <- reactive({ input$rmvMinimum })
  tip.estimate <- reactive({ input$tip.estimate })
  
  
  
  observe({
    label.stg <- paste("Time Slice (", time.unit(), "):", sep = "")
    if(time.unit() == "frame") {
      updateSliderInput(session, inputId = "t.slice", label = label.stg,
        value = time.step()*100, min = 1, max = dim(kymo())[1]*time.step(), step = time.step())
    }else{
      updateSliderInput(session, inputId = "t.slice", label = label.stg,
        value = time.step()*100, min = 0, max = dim(kymo())[1]*time.step(), step = time.step())
    }
  })
  
  # t.slice.ind <- reactive({ input$t.slice })
  
  t.slice.ind <- reactive({
    if(time.unit() == "frame") {
      return(input$t.slice)
    }else{
      return(input$t.slice / time.step())
    }
  })
  
  
  red.bias.tab3 <- reactive({ input$redBias.tab3 })
  restrict.x <- reactive({ input$restrict.X })
  
  ## Calls:
  tip.lst <- reactive({ FindTip( kymo(), #NEW!!!
    use.smooth(), kymo.span.n(), kymo.loess.dg(), 
    n.pts(), mad.tol(), qntl, fit.right, fix.slope(),     
    fluo.frac(), min.chunk.size(), 
    coarse.kymo.span.n(), coarse.kymo.loess.dg(),
    use.coarse(), rmv.min(), algorithm() ) #TODO Change FindTip function
  })
  kymo.smth <- reactive({ tip.lst()$kymo.smth })
  tip.tbl <- reactive({ tip.lst()$tip.tbl })
  tip.loc <- reactive({ tip.tbl()[, which(colnames(tip.tbl()) == tip.estimate())] })
  
  # Inputs and calls for Tip location smoothing (Fourth Tab): ----
  tip.span.n <- reactive({ input$tip.spn })
  tip.span.dg <- reactive({ input$tip.spn.dg })
  
  out.rm.algorithm <- reactive({ input$outlierAlgorithm })
  
  rm.tip.out <- reactive({ input$removeTipOut })
  rm.tip.out.spn <- reactive({ input$remove.tip.out.spn }) 
  rm.tip.out.dg <- reactive({ input$remove.tip.out.dg }) 
  
  tol <- reactive({ input$tol.tab4 }) 
  px.tol <- reactive({
    if(input$px.tol == 0){
      return(NULL)
    }else{
      return(input$px.tol)
    }
  })
  
  red.bias.tab4 <- reactive({ input$redBias.tab4 })
  rmv.lgnd <- reactive({ input$rmv.legend })
  
  ## Calls:
  tip.loc.smth.lst <- reactive({ RefineTipLoc(tip.loc(), tip.span.n(), tip.span.dg(), rm.tip.out(), rm.tip.out.dg(), 
    rm.tip.out.spn(), tol(), px.tol(), out.rm.algorithm()) })
  
  tip.loc.out <- reactive({ tip.loc.smth.lst()$tip.loc })
  tip.loc.smth <- reactive({ tip.loc.smth.lst()$tip.loc.smth })
  out.indx <- reactive({ tip.loc.smth.lst()$out.indx })
  
  
  # Inputs and calls for Extracting Flourescence Series (Fifth Tab): ----
  tip.mrgn <- reactive({ input$tip.margin })
  avg.width <- reactive({ input$avrg.width })
  use.median.for.avg <- reactive({ input$useMedian.avrg })
  
  red.bias.tab5 <- reactive({ input$redBias.tab5 })
  y.lab.tab5 <- reactive({ input$y.lab.tab5 })
  
  manual.ROI <- reactive({ input$manualROI })
  
  roi.px <- reactive({ input$roi.px })
  
  
  observeEvent(eventExpr = {
    input$manualROI
    input$next.tab45
  },
    handlerExpr = {
      if(!manual.ROI()) {
        updateNumericInput(session, inputId = "roi.px", 
          value = max(0, median((tip.loc.smth() + tip.mrgn()) - tip.tbl()[, 5], na.rm = TRUE)))
      }
    })
  
  ## Calls:
  kymo.align <- reactive({ AlignByTip(tip.loc = tip.loc.smth() + tip.mrgn(), imaj = kymo()) })
  
  tip.fluo.lst <- reactive({ ExtractFluoTimeSeries(kymo.align(), roi.px(), avg.width(), use.median = use.median.for.avg(), max.tip = min(tip.loc.smth(), na.rm = TRUE)) }) # Wait for Dani to finish this and put '()' in all of them
  
  
  
  # Inputs and calls for Filtering the Kymograph (Sixth Tab): ----
                # per.min <- reactive({ input$period.slider[1] })
                # per.max <- reactive({ input$period.slider[2] })
                # 
                # trim.time <- reactive({ input$trimTime })
                # trim.length <- reactive({ input$trimLength })
                # 
                # tip.mrgn.filt <- reactive({ input$tip.mrgn.filt })
                # manual.ROI.tab6 <- reactive({ input$manualROI_tab6 })
                # roi.px.tab6 <- reactive({ input$roi.px.tab6 })
                # 
                # avg.width.tab6 <- reactive({ input$avrg.width_tab6 })
                # use.median.for.avg.tab6 <- reactive({ input$useMedian.avrg_tab6 })
                # 
                # red.bias.tab6 <- reactive({ input$redBias.tab6 })
                # red.bias.tab6.trim <- reactive({ input$redBias.tab6.trim })
                # y.lab.tab6 <- reactive({ input$y.lab.tab6 })
                # rmv.lgnd.tab6 <- reactive({ input$rmv.legend.tab6 })
                # 
                # decision.Tab6 <- reactive({ input$decision.tab6 })
                # decision.Tab6.Run <- reactive({ input$decision.tab6.run })
                # 
                # ## Calls:
                # observeEvent(eventExpr = {
                #   input$manualROI_tab6
                #   input$next.tab56
                # },
                #   handlerExpr = {
                #     if(!manual.ROI.tab6()) {
                #       updateNumericInput(session, inputId = "roi.px.tab6", 
                #         value = max(0, median((tip.loc.smth() + tip.mrgn.filt()) - tip.tbl()[, 5], na.rm = TRUE)))
                #     }
                #   })
                # 
                # observeEvent(eventExpr = {
                #   input$t.step
                #   kymo()
                # },
                #   handlerExpr = {
                #     updateSliderInput(session, inputId = "period.slider", 
                #       value = c((time.step() * 4), ((time.step() * dim(kymo())[1]) / 3)),
                #       min = (time.step() * 2), max = (time.step() * dim(kymo())[1]), step = time.step())
                #   })
                # 
                # 
                # observeEvent(eventExpr = {
                #   kymo()
                # },
                #   handlerExpr = {
                #     updateSliderInput(session, inputId = "trimTime", 
                #       value = c(1, dim(kymo())[1]), max = dim(kymo())[1])
                #   })
                # 
                # observeEvent(eventExpr = {
                #   kymo()
                # },
                #   handlerExpr = {
                #     updateSliderInput(session, inputId = "trimLength", 
                #       value = c(1, dim(kymo())[2]), max = dim(kymo())[2])
                #   })
                # 
                # kymo.filt <- reactive({ FilterKymo(kymo(), time.step(), low.per = per.min(), high.per = per.max()) })
                # kymo.filt.align <- reactive({ AlignByTip(tip.loc = tip.loc.smth() + tip.mrgn.filt(), imaj = kymo.filt()) })
                # kymo.trim <- reactive({ kymo.filt.align()[trim.time()[1]:trim.time()[2], trim.length()[1]:trim.length()[2]] })
                # tip.filt.fluo.lst <- reactive({ ExtractFluoTimeSeries(kymo.trim(), roi.px.tab6(), avg.width.tab6(), use.median = use.median.for.avg.tab6(), max.tip = min(tip.loc.smth(), na.rm = TRUE)) })
                # 
  
  # Calls for Last Tab ("Save"): ----
  
  # File.type <- reactive({ input$dataFileType })  # input$dataFileType is in 'extension.out' below
  
  all.ts <- reactive({
    cbind("time" = 0:(dim(kymo())[1] - 1) * time.step(), 
      "tip.loc.raw" = tip.loc() * px.sz(), 
      "tip.loc.smth" = tip.loc() * px.sz(), 
      "growth" = c(NA, diff(tip.loc.smth() * px.sz())/time.step()), 
      tip.fluo.lst()[[2]])
  })
  
  # all.trim.fluo <- reactive({
  #   cbind("time" = ((trim.time()[1]:trim.time()[2]) - 1) * time.step(),
  #     tip.filt.fluo.lst()[[2]])
  # })
  
  px.tol_forPAr <- reactive({
    if(is.null(px.tol())){
      return(0)
    }else{
      return(px.tol())
    }
  })
  
  all.par.vals <- reactive({
    c(fl.nm(), 
      px.sz(), 
      px.unit(), 
      time.step(), 
      time.unit(), 
      red.bias.tab2(), 
      use.smooth(), 
      kymo.span.n(), 
      kymo.loess.dg(), 
      n.pts(), 
      mad.tol(),
      fix.slope(), 
      qntl,
      fit.right,
      
      fluo.frac(),
      min.chunk.size(),
      coarse.kymo.loess.dg(),
      coarse.kymo.span.n(),
      use.coarse(),
      rmv.min(),
      algorithm(),
      
      tip.estimate(),
      t.slice.ind(),
      red.bias.tab3(), 
      restrict.x(),
      tip.span.n(),
      tip.span.dg(),
      
      out.rm.algorithm(),
      
      rm.tip.out(),
      rm.tip.out.dg(), 
      rm.tip.out.spn(),
      tol(),
      px.tol_forPAr(),
      rmv.lgnd(),
      red.bias.tab4(),
      tip.mrgn(),
      manual.ROI(),
      roi.px(), 
      avg.width(),
      use.median.for.avg(),
      red.bias.tab5(), 
      y.lab.tab5()#,
      # per.min(), 
      # per.max(),
      # trim.time()[1],
      # trim.time()[2],
      # trim.length()[1],
      # trim.length()[2],
      # tip.mrgn.filt(), 
      # manual.ROI.tab6(), 
      # roi.px.tab6(), 
      # avg.width.tab6(),
      # use.median.for.avg.tab6(),
      # red.bias.tab6(),
      # red.bias.tab6.trim(),
      #y.lab.tab6()
      )
  })
  
  all.par.nms <- list( "fl.nm",
    "px.sz",
    "px.unit", 
    "time.step", 
    "time.unit",
    "red.bias.tab2", 
    "use.smooth", 
    "kymo.span.n", 
    "kymo.loess.dg",
    "n.pts", 
    "mad.tol",
    "fix.slope", 
    "qntl", 
    "fit.right",
    
    "fluo.frac",
    "min.chunk.size",
    "coarse.kymo.loess.dg",
    "coarse.kymo.span.n",
    "use.coarse",
    "rmv.min",
    "algorithm",
    
    "tip.est",
    "t.slice.ind",
    "red.bias.tab3", 
    "restrict.x",
    "tip.span.n",
    "tip.span.dg",
    
    "out.rm.algorithm",
    
    "rm.tip.out",
    "rm.tip.out.dg", 
    "rm.tip.out.spn",
    "tol",
    "px.tol",
    "rmv.lgnd",
    "red.bias.tab4",
    "tip.mrgn",
    "manual.ROI",
    "roi.px", 
    "avg.width",
    "use.median.for.avg",
    "red.bias.tab5", 
    "y.lab"#,
    # "per.min", 
    # "per.max",
    # "trim.time.1",
    # "trim.time.2",
    # "trim.length.1",
    # "trim.length.2",
    # "tip.mrgn.filt", 
    # "manual.ROI.tab6", 
    # "roi.px.tab6", 
    # "avg.width.tab6",
    # "use.median.for.avg.tab6",
    # "red.bias.tab6",
    # "red.bias.tab6.trim",
    #"y.lab.tab6"
    )
  all.par.tbl <- reactive({
    cbind("parameter" = all.par.nms, "value" = all.par.vals())
  })
  
  
  extension.out <- reactive({ input$dataFileType })
  sep.out <- reactive({
    if(extension.out() == "txt"){
      return(input$txt.sep)
    }
    if(extension.out() == "csv"){
      return(input$csv.sep)
    }
  })
  
  
  Batch.InFile <- reactive({ input$BatchZip })
  Batch.paths <- reactive({ Batch.InFile()$datapath })  # datapath used to read the data
  Batch.names <- reactive({ Batch.InFile()$name })  # name used for naming the files
  
  Batch.mats <- reactive({ lapply(Batch.paths(), ReadKymo) })  # reading the files and putting them all in one list 
  
 
  # Outputs for the Second Tab ("Data Upload"): ----
  output$initial.plot_tab2 <- renderPlot({
    if(is.null(inFile())) { # no file from user
      return(NULL)
    }
    PlotKymoInput(kymo_untreated(), px.sz(), px.unit(), time.step(), time.unit(), red.bias = red.bias.tab2()) # keep in mind that if necessary 'kymo_v <- kymo()' will also work
  })
  
  # Outputs for the Third Tab ("Tip Detection"): ----
  output$main.plot_tab3 <- renderPlot({
    PlotKymoWithAllTip(kymo.smth(), tip.tbl(), px.sz(), px.unit(), time.step(), 
      time.unit(), brks = NULL, red.bias = red.bias.tab3(), restrict.x = restrict.x())
  })
  output$slice.plot_tab3 <- renderPlot({
    t.slice.raw <- (kymo() - min(c(kymo())))[t.slice.ind(), ]
    t.slice.smth <- kymo.smth()[t.slice.ind(), ]
    tip.ests <- tip.tbl()[t.slice.ind(), ]
    
    PlotTimeSlice(t.slice.raw = t.slice.raw, 
      t.slice.smth = t.slice.smth,
      tip.ests = tip.ests, restrict.x = restrict.x()) # there is no input for 'cx', so I'm leaving it out of the call
  })
  
  # Outputs for the Fourth Tab ("Tip Location Smoothing"): ----
  output$main.plot_tab4 <- renderPlot({
    PlotKymoWithTip(kymo.smth(), tip.loc.raw = tip.loc(), tip.loc.smth = tip.loc.smth(), tip.loc.out = tip.loc.out(), out.indx = out.indx(),
      px.sz(), px.unit(), time.step(), time.unit(), brks = NULL, red.bias = red.bias.tab4(), rmv.lgnd = rmv.lgnd())
  })
  
  output$dygraphTab4 <- renderDygraph({
    PlotKymoWithDygraph(kymo.smth(), tip.loc.raw = tip.loc(), tip.loc.smth = tip.loc.smth(), px.sz(), px.unit(),
      time.step(), time.unit(), brks = NULL, red.bias = red.bias.tab4(), rmv.lgnd = rmv.lgnd())
  })
  # Outputs for the Fifth Tab ("Extract Fluorescence Series"): ----
  output$main.plot_tab5 <- renderPlot({
    PlotKymoAndFluo(kymo.align(), tip.fluo.ts = tip.fluo.lst()[[2]], roi.ini.pxs = tip.fluo.lst()[[1]], avg.width(), max.tip = max(tip.loc.smth(), na.rm = TRUE), 
      px.sz(), px.unit(), time.step(), time.unit(), brks = NULL, red.bias = red.bias.tab5(), y.lab = y.lab.tab5())
  })
  output$dygraph1Tab5 <- renderDygraph({
    PlotKymoAndFluoWithDygraph1(kymo.align(), tip.fluo.ts = tip.fluo.lst()[[2]], roi.ini.pxs = tip.fluo.lst()[[1]], avg.width(), max.tip = max(tip.loc.smth(), na.rm = TRUE), 
      px.sz(), px.unit(), time.step(), time.unit(), brks = NULL, red.bias = red.bias.tab5(), y.lab = y.lab.tab5())
  })
  output$dygraph2Tab5 <- renderDygraph({
    PlotKymoAndFluoWithDygraph2(kymo.align(), tip.fluo.ts = tip.fluo.lst()[[2]], roi.ini.pxs = tip.fluo.lst()[[1]], avg.width(), max.tip = max(tip.loc.smth(), na.rm = TRUE), 
      px.sz(), px.unit(), time.step(), time.unit(), brks = NULL, red.bias = red.bias.tab5(), y.lab = y.lab.tab5())
  })
  
  # Outputs for the Sixth Tab ("Filter Kymograph"): ----
      # output$main.plot_tab6 <- renderPlot({
      #   if(decision.Tab6() > 0 ) {  #| input$tab.save > 0
      #     return(NULL)
      #   } 
      #   if(decision.Tab6() == 0 & decision.Tab6.Run() > 0) {
      #     PlotFilteredKymo(kymo(), kymo.filt(), max.tip = max(tip.loc.smth(), na.rm = TRUE), px.sz(), px.unit(), 
      #       time.step(), time.unit(), brks = NULL, red.bias.raw = red.bias.tab4(), 
      #       red.bias.filt = red.bias.tab6(), rmv.lgnd.tab6(), cex.ax = 1.15)
      #   }
      # })
      # 
      # output$plot2.tab6 <- renderPlot({
      #   if(decision.Tab6() > 0 ) {  #| input$tab.save > 0
      #     return(NULL)
      #   }
      #   if(decision.Tab6() == 0 & decision.Tab6.Run() > 0) {
      #     PlotTrimKymoAndFluo(kymo.trim(), trim.length(), trim.time(), tip.filt.fluo.ts = tip.filt.fluo.lst()[[2]],
      #       roi.ini.pxs = tip.filt.fluo.lst()[[1]], avg.width = avg.width.tab6(), max.tip = max(tip.loc.smth(), na.rm = TRUE), 
      #       px.sz(), px.unit(), time.step(), time.unit(), rmv.lgnd.tab6(), brks = NULL, red.bias = red.bias.tab6.trim(), y.lab = y.lab.tab6())
      #   }
      # })
  
  
  
  # Outputs for the Last Tab ("Save"): ----
  
  output$out.Data <- downloadHandler(
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      paste(input$out.Data_name, "zip", sep = ".")
    },
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(fname) {
      tmpdir <- tempdir()
      setwd(tempdir())
      
      #print(tempdir())
      # fs <- c("rock.csv", "pressure.csv", "cars.csv")
      # write.csv(datasetInput()$rock, file = "rock.csv", sep =",")
      # write.csv(datasetInput()$pressure, file = "pressure.csv", sep =",")
      # write.csv(datasetInput()$cars, file = "cars.csv", sep =",")
      #print (fs)
      #Switched the order of fs and the 'write.'s (to check for batch mode)
      pdf("TipDetection.pdf")
      PlotKymoWithAllTip(kymo.smth(), tip.tbl(), px.sz(), px.unit(), time.step(), 
        time.unit(), brks = NULL, red.bias = red.bias.tab3(), restrict.x = restrict.x())
      dev.off()
      pdf("TipLocation.pdf")
      PlotKymoWithTip(kymo.smth(), tip.loc.raw = tip.loc(), tip.loc.smth = tip.loc.smth(), tip.loc.out = tip.loc.out(), out.indx = out.indx(),
        px.sz(), px.unit(), time.step(), time.unit(), brks = NULL, red.bias = red.bias.tab4(), rmv.lgnd = rmv.lgnd())
      dev.off()
      pdf("FluorescenceSeries.pdf")
      PlotKymoAndFluo(kymo.align(), tip.fluo.ts = tip.fluo.lst()[[2]], roi.ini.pxs = tip.fluo.lst()[[1]], avg.width(), max.tip = max(tip.loc.smth(), na.rm = TRUE), 
        px.sz(), px.unit(), time.step(), time.unit(), brks = NULL, red.bias = red.bias.tab5(), y.lab = y.lab.tab5())
      dev.off()
      # pdf("FilterKymograph.pdf")
      #               # if (decision.Tab6() == 0){  # & input$tab.save == 0
      #               #   PlotFilteredKymo(kymo(), kymo.filt(), max.tip = max(tip.loc.smth(), na.rm = TRUE), px.sz(), px.unit(), 
      #               #     time.step(), time.unit(), brks = NULL, red.bias.raw = red.bias.tab4(), 
      #               #     red.bias.filt = red.bias.tab6(), rmv.lgnd.tab6(), cex.ax = 1.15)
      #               #   PlotTrimKymoAndFluo(kymo.trim(), trim.length(), trim.time(), tip.filt.fluo.ts = tip.filt.fluo.lst()[[2]],
      #               #     roi.ini.pxs = tip.filt.fluo.lst()[[1]], avg.width = avg.width.tab6(), max.tip = max(tip.loc.smth(), na.rm = TRUE), 
      #               #     px.sz(), px.unit(), time.step(), time.unit(), rmv.lgnd.tab6(), brks = NULL, red.bias = red.bias.tab6.trim(), y.lab = y.lab.tab6())
      #               # }else{
      #   wrng.txt <- paste("The Parameters 'per.min', 'per.max', 'trim.time.1', 'trim.time.2', 
      #     'trim.length.1', 'trim.length.2', 'tip.mrgn.filt', 'manual.ROI.tab6', 
      #     'roi.px.tab6', 'avg.width.tab6',' use.median.for.avg.tab6', 
      #     'red.bias.tab6', 'red.bias.tab6.trim' and 'y.lab.tab6' present in 
      #     \"Parameters.txt(.csv)\" all pertain to the sixth tab (\"Filter Kymograph\") 
      #     and since it still under construction, such paramenters do not have any reliable value")
      #   plot(NA, xlim=c(0,50), ylim=c(0,50), bty='n', xaxt='n', yaxt='n', xlab='', ylab='')
      #   text(20, 40, labels = wrng.txt, cex = 0.75)
      #               #}
      # dev.off()
      SaveTable(kymo.align(), fl.nm = "Tip_aligned_kymograph", extension.out(), sep.out())
      # SaveTable(kymo.filt.align(), fl.nm = "Tip_aligned_filtered_kymograph", extension.out(), sep.out())
      # SaveTable(kymo.trim(), fl.nm = "Tip_aligned_filtered_trimmed_kymograph", extension.out(), sep.out())
      SaveTable(all.ts(), fl.nm = "All_time_series", extension.out(), sep.out())
      #SaveTable(all.trim.fluo(), fl.nm = "All_filtered_time_series", extension.out(), sep.out())
      SaveTable(all.par.tbl(), fl.nm = "Parameters", extension.out(), sep.out())
      
      fs <- c("TipDetection.pdf", "TipLocation.pdf", "FluorescenceSeries.pdf", #"FilterKymograph.pdf", 
        paste("Tip_aligned_kymograph.", extension.out(), sep = ""), 
        # paste("Tip_aligned_filtered_kymograph.", extension.out(), sep = ""), paste("Tip_aligned_filtered_trimmed_kymograph.", extension.out(), sep = ""), 
        paste("All_time_series.", extension.out(), sep = ""), paste("Parameters.", extension.out(), sep = "")) #paste("All_filtered_time_series.", extension.out(), sep = ""),
      zip(zipfile=fname, files=fs)
    },
    contentType = "application/zip"
  )
  
  output$Batch.out <- downloadHandler(
    filename = function() {
      paste(input$out.Data_name, "_Batch", ".zip", sep = "")
    },
    content = function(fname) {
      tmpdir <- tempdir()
      setwd(tempdir())
      
      fs <- c()
      for (i in 1:length(Batch.paths())) {
        # running the function for each df/name - it'll create all the files and return the names of such files
        nms <- RunBatch(df = Batch.mats()[[i]], df.name = Batch.names()[i], parameters = all.par.tbl(), 
          extension = extension.out(), sep = sep.out()) #, decision.Tab6 = decision.Tab6())  
        fs <- c(fs, nms) # getting the names of all the files of all dfs
      }
      
      zip(zipfile=fname, files=fs) # creating the zip with all the files generated
      
    },
    contentType = "application/zip"
  )
  
  
    }# End of Server _________________________________________________________


# Calling ShinyApps ------------------------------------------------------
shinyApp(ui = ui, server = server)
