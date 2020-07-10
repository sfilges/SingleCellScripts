suppressMessages(library(tidyverse, quietly = TRUE))
suppressMessages(library(patchwork, quietly = TRUE))
suppressMessages(library(shiny, quietly = TRUE))
suppressMessages(library(shinyFiles, quietly = TRUE))
suppressMessages(library(shinyWidgets, quietly = TRUE))
suppressMessages(library(DT, quietly = TRUE))
suppressMessages(library(shinydashboard, quietly = TRUE))

### ui end, to browse to desired folder
ui = dashboardPage(
  dashboardHeader(
    title = 'Shiny Seurat'
  ),
  
  # Define menu items on the sidebar
  dashboardSidebar(
    sidebarMenu(
      menuItem(
        text = 'Dashboard',
        tabName = 'dashboard',
        icon = icon('dna')
      ),
      menuItem(
        text = 'Github',
        icon = icon('git'),
        href = 'https://github.com/ozimand1as/SingleCellScripts'
      )
    )
  ),
  
  dashboardBody(
    
    fluidRow(
      box(
        title = "Input",
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        width = 12,
        
        textInput("text1", "Sample Name:"),
            
        shinyDirButton('directory', 'Import data', 'Please select a folder'),
            
        textOutput("selected_var",inline = TRUE),
            
        selectInput(
              inputId = 'samples', width = "100%",
              label = 'Select samples to analyse:',
              choices = '',
              multiple = TRUE
        ),
        actionButton(
              inputId = 'import_data',
              label = 'Load data',
              icon = icon('file-import'),
              style = "margin-bottom: 10px;margin-left: 5px;"
        )
      ),
      
      box(
        title = "Quality Control",
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        
        tabBox(
          width = 12,
          type = "tabs",
          tabPanel(
            title = "QC- Plot",
            fluidRow(
              # Option for plot customisation
              dropdownButton(
                tags$h3('Customise plot'),
                sliderInput(
                  inputId = "nFeature",
                  label =  "Number of genes to use:",
                  min = 0, max = 20000,
                  value = c(2000,8000), step = 500,
                  post = " genes", sep = ","
                ),
                sliderInput(
                  inputId = "maxMito",
                  label =  "Maximum mitochondrial percentage:",
                  min = 0, max = 100,
                  value = 10, step = 1,
                  post = " %", sep = ","
                ),
                circle = FALSE,
                status = 'default',
                icon = icon('gear'),
                width = '300px',
                tooltip = tooltipOptions(title = 'Click to customise plot!')
              ),
              plotOutput("qc_scatter"),
              downloadButton(
                outputId = 'download_qc_scatter.pdf',
                label = 'Download figure'
              )
            )
          )
        )
      )
      
    )
  )
)

### extracting the folder path
server = function(input, output, session) {
  
  source("import_functions.R")
  
  volumes <- c(Home = fs::path_home(),
               'R Installation' = R.home(),
               getVolumes()())
  shinyDirChoose(input, 'directory', roots=volumes, session=session)
  
  # Values is a reactive object to which a umiExperiment object is added in
  # the data slot.
  values <- reactiveValues()
  values$list <- data.frame(path = character(0), name = character(0))
  
  path1 <- reactive({
    main <- parseDirPath(volumes, input$directory)
    
    if(! identical(main, character(0))){
      return(main)
    } else {
      return(NULL)
    }
  })
  
  ### constructing the 3 file paths
  list <- observe({
    if(is.null(path1())){
      return(NULL)
    }
      
      
      
      print(path1())
      x <- isolate(c(path = path1(), name = input$text1))
      isolate(values$list <- dplyr::bind_rows(values$list,x))
      
      # Render variant call table in app
      output$varDataTable <- DT::renderDataTable({
        values$list
      })
      
      output$selected_var <- renderText({
        paste('You selected <',input$text1, '> at: ',  path1(), sep = '')
      })
      
      return(values$list)
  })
  
  # make reactive expresion of input values
  sample_settings <- reactive({input$samples})
  
  #------------- Update assays list -------------
  # Update assay and sample list based on initially loaded object, meaning that
  # the lists will be visible even if filter are applied
  observe({
    
    if (is.null(values$list)){
      return(NULL)
    }
    
    data <- values$list
    
    print(data)
    
    updateSelectInput(
      session = session,
      inputId = 'samples',
      choices = data$name
    )
    
  })
  
  seurat <- reactiveValues(object=NULL)
  
  #------ Import data --------
  observeEvent(input$import_data, {
    
    df <- values$list
    samples = sample_settings()
    
    #print(df)
    
    object <- mergeSeuratObjects(df)
    
    print(object)
    
    seurat$object <- object

  })

  
  
  observe({
    
    if(is.null(seurat$object)){
      return(NULL)
    }
    
    object <- seurat$object
    print(object)
    
    x <- input$nFeature[1]
    y <- input$nFeature[2]
    z <- input$maxMito
    
    expr <- FetchData(object = object, vars = var1)
    object[, which(x = expr > low & expr < high)]
    
    output$qc_scatter <- renderPlot({
      qc_scatter(object = object)
    })
    
  })
}

shinyApp(ui = ui, server = server)