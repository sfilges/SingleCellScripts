library(shiny)
library(shinyFiles)

### ui end, to browse to desired folder
ui = fluidPage(
  textInput("text1", "Sample Name:"),
  shinyDirButton('directory', 'Folder select', 'Please select a folder'),
  textOutput("selected_var"),
  DT::dataTableOutput('varDataTable'),
  actionButton(
    inputId = 'import_data',
    label = 'Load data',
    icon = icon('file-import'),
    style = "margin-bottom: 10px;margin-left: 5px;"
  )
)

### extracting the folder path
server = function(input, output, session) {
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
    if(!is.null(path1())){
      print(path1())
      x <- isolate(c(path = path1(), name = input$text1))
      isolate(values$list <- dplyr::bind_rows(values$list,x))
      
      # Render variant call table in app
      output$varDataTable <- DT::renderDataTable({
        values$list
      })
      
      output$selected_var <- renderText({
        paste('You selected',input$text1, ' at: ',  path1(), sep = '')
      })
      
      return(values$list)
    }
  })
  
  
  #------ Import data --------
  
  seurat <- observeEvent(input$import_data, {
    
    df <- values$list
    
    #print(df)
    
    object <- mergeSeuratObjects(df)
    
    print(object)
    
    return(object)

  })
  
  
}

shinyApp(ui = ui, server = server)