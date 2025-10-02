library(shiny)
library(bslib)
library(Seurat)
library(DT)
library(rlang)
library(shinythemes)
library(shinyFeedback)
library(rsconnect)
library(qs)
library(shinyBS)

# source other R scripts
source('helpers.R')
# increase the file upload size allowed
options(shiny.maxRequestSize = 4000 * 1024^2)
# specify an application port
options(shiny.host = "0.0.0.0")
options(shiny.port = 8180)


# define the ui ----
ui <- fluidPage(theme = shinytheme("yeti"),
  # title styling
  div(
    style = "text-align: center; display: flex; flex-direction: column; align-items: center;",
    titlePanel("Goessling Lab Single Cell Browser"),
    img(src = 'logo.png', width = "200px") 
  ),
  
  div(style = "height: 20px;"), # add spacing
  
  tabsetPanel(
    # data import tab ----
    tabPanel("Import Data",
             
       # title and description above everything ----
       h3("Upload or Select Dataset", 
          style = "font-weight: bold; font-size: 20px; color: #1E3A8A; text-align: center; margin-top: 20px; margin-bottom: 15px;"),
       
       p("To begin analysis, upload your annotated single cell dataset or select an existing dataset from the list below. Wait for the data to load and process. Then, set the appropriate cell type column within the metadata for analysis.",
         style = "font-size: 15px; color: #555; text-align: center; margin-bottom: 20px;"),
       
       p("Note: if uploading a file, file extension must be .qs",
         style = "font-weight: bold; font-size: 15px; color: #555; text-align: center; margin-bottom: 30px;"),
       
       # line break here
       hr(style = "border: 1px solid #ccc; margin-top: 20px; margin-bottom: 20px;"), 
             
       # main layout ----
       fluidRow(
         column(6,
                fileInput("file_upload", "Upload dataset:", accept = c(".qs")),
                div(style = "margin-top: -30px"),
                p("Please wait for data to finish upload before processing",
                  style = "font-size: 10px; font-weight: bold; margin-top: 30px;"),
                # buttons in the same row
                fluidRow(
                  column(9,
                         loadingButton("upload_data", "Process Dataset for Analysis",
                                       style = "width: 100%; margin: 8px 0;")
                  ),
                  column(3,
                         div(style = "margin-left: -20px; margin-top: 20px;",
                             loadingButton("cancel_button", "Cancel",
                                           style = "width: 70%; font-size: 11px; padding: 4px 6px; background-color: #f8d7da; color: #721c24; border-color: #f5c6cb;")
                         )
                  )
                ),
                selectInput("file", "OR choose a dataset:",
                            choices = list.files("/mnt/s3-bucket/data", pattern = "\\.qs$", full.names = FALSE),
                            selected = NULL),
                loadingButton("load_data", "Load & Process Dataset",
                              style = "width: 71%; margin: 10px auto;")
         ),
         column(6,
                selectInput("cell.type.meta", "Choose the metadata cell type column:", choices = NULL),
                actionButton("set_cell_type", "Set Cell Type Column"),
                DT::dataTableOutput("cell_type_table")
         )
       ),
    ),
    
    # umap tab: all cell types -----
    tabPanel("UMAP",
       # title and description above everything ----
       h3("Gene Expression Visualization with UMAPs", 
          style = "font-weight: bold; font-size: 20px; color: #1E3A8A; text-align: center; margin-top: 20px; margin-bottom: 15px;"),
       
       p("In this section, you can enter a gene of interest to visualize its expression across different cells using UMAPs and a violin plot. 
          The generated plots will show how the gene is distributed across the dataset, helping you to identify clusters or patterns of gene expression.",
         style = "font-size: 15px; color: #555; text-align: center; margin-bottom: 30px;"),
       
       # line break here
       hr(style = "border: 1px solid #ccc; margin-top: 20px; margin-bottom: 20px;"), 
       
       # main layout - sidebar and main panel ----      
      sidebarLayout(
        sidebarPanel(
          textInput("gene", "Enter a gene of interest:", placeholder = "Type gene name..."),
          actionButton("generate_umaps", "Generate UMAPs") 
          ),
          mainPanel(
            plotOutput("feature_plot"),  
            plotOutput("umap_plot"),
            plotOutput("violin_celltype_plot")
        )
      )
    ),
    
    # differential expression table: all cell types ----
    tabPanel("Differential Expression",
       # title and description above everything ----
       h3("Differential Expression Analysis", 
          style = "font-weight: bold; font-size: 20px; color: #1E3A8A; text-align: center; margin-top: 20px; margin-bottom: 15px;"),
       
       p("In this section, you can compare gene expression in the selected cell type to all other cell types in the dataset. 
       Use the controls to select a cell type and view the corresponding differential expression table for detailed insights.",
         style = "font-size: 15px; color: #555; text-align: center; margin-bottom: 30px;"),
       
       p("Note: You must select the cell type metadata column in the Import Data tab to use this feature",
         style = "font-weight: bold; font-size: 15px; color: #555; text-align: center; margin-bottom: 30px;"),
       
       # line break here
       hr(style = "border: 1px solid #ccc; margin-top: 20px; margin-bottom: 20px;"), 
       
       # main layout - sidebar and main panel ----
       sidebarLayout(
         sidebarPanel(
           selectInput("cell.type.DE", "Choose a cell type:", choices = NULL), # default empty
           actionButton("display.de.table", "Display Table"),
           hr(),  # Horizontal line for separation
           h4("How to interpret this table:", 
              style = "font-weight: bold; margin-bottom: 6px; font-size: 15px;"),
           # how to interpret the table
           tags$ul(
             tags$li(tags$b("p_val:"), " The statistical significance of the gene’s differential expression in the selected cell type compared to all others."),
             tags$li(tags$b("avg_log2FC:"), " The average log2 fold change in gene expression between the selected cell type and other cell types."),
             tags$li(tags$b("pct.1:"), " The percentage of cells in the selected cell type that express the gene."),
             tags$li(tags$b("pct.2:"), " The percentage of cells in all other cell types that express the gene."),
             tags$li(tags$b("p_val_adj:"), " The adjusted p-value, corrected for multiple testing.")
           ),
         ),
         mainPanel(
           DT::dataTableOutput("de_table")
         )
       )
    ),
    
    
    tabPanel("Cell Type Specific",
     # title and description above everything ----
     h3("Cell Type Specific Analysis", 
        style = "font-weight: bold; font-size: 20px; color: #1E3A8A; text-align: center; margin-top: 20px; margin-bottom: 15px;"),
     
     p("This section allows you to explore gene expression patterns for specific cell types. 
        Choose a cell type, enter a gene of interest, and customize the data groupings to visualize the results.",
       style = "font-size: 15px; color: #555; text-align: center; margin-bottom: 20px;"),
     
     p("Note: You must select the cell type metadata column in the Import Data tab to use this feature",
       style = "font-weight: bold; font-size: 15px; color: #555; text-align: center; margin-bottom: 30px;"),
     
     # line break here
     hr(style = "border: 1px solid #ccc; margin-top: 20px; margin-bottom: 20px;"), 
     
     # main layout - sidebar and main panel ----
      sidebarLayout(
        sidebarPanel(
          selectInput("cell.type.subcluster", "Choose a cell type:", choices = NULL), # default empty
          textInput("gene.subcluster", "Enter a gene of interest:", placeholder = "Type gene name..."),
          # Wrap the radio button and info icon together for alignment
          tags$div(
            style = "display: flex; align-items: center;",
            radioButtons("non.zero", 
                         "Display cells with zero expression of selected gene?", 
                         choices = list("Yes" = 1, "No" = 2), 
                         selected = 1),
            tags$i(class = "fa fa-info-circle", 
                   id = "info_non_zero", 
                   style = "margin-left: 5px; cursor: pointer; font-size: 16px; color: #007BFF;")
          ),
          
          # Tooltip for the radio button info icon
          bsTooltip("info_non_zero", 
                    "Choose whether to display only cells with non-zero expression of the selected gene in the violin plot.",
                    placement = "right", trigger = "hover"),
          

          # wrap the dropdown and info icon together for alignment
          tags$div(
            style = "display: flex; align-items: center;",
            selectInput("meta.subcluster", "Choose a metadata column to group by:", choices = NULL), # default empty
            tags$i(class = "fa fa-info-circle", 
                   id = "info_meta", 
                   style = "margin-left: 5px; cursor: pointer; font-size: 16px; color: #007BFF;")
          ),
          actionButton("display.subcluster", "Submit"),
          # define tooltip content
          bsTooltip("info_meta", 
                    "Select a metadata column (e.g., treatment, timepoint) to group/label the data in visualizations.", 
                    placement = "right", trigger = "hover")
        ),
        mainPanel(
          plotOutput("umap.subcluster"),
          plotOutput("umap.gene.subcluster"),
          plotOutput("violin.subcluster")  
        )
      )
    )
  )
)

# define the server logic ----
server <- function(input,output,session){
  # create a var to store the seurat object globally
  seurat_obj <- reactiveValues(data = NULL, cell_type_column = NULL) # initally empty
  
  # load the seurat object ----
  observeEvent(input$load_data, {
    # ensure input is not null
    req(input$file)
    
    # Reset old data and metadata UI selections
    seurat_obj$data <- NULL
    seurat_obj$cell_type_column <- NULL
    updateSelectInput(session, "cell.type.meta", choices = character(0), selected = NULL)
    updateSelectInput(session, "cell.type.DE", choices = character(0), selected = NULL)
    updateSelectInput(session, "cell.type.subcluster", choices = character(0), selected = NULL)
    updateSelectInput(session, "meta.subcluster", choices = character(0), selected = NULL)
    output$cell_type_table <- DT::renderDT({ NULL })
    output$de_table <- DT::renderDT({ NULL })
    output$feature_plot <- renderPlot({ NULL })
    output$umap_plot <- renderPlot({ NULL })
    output$violin_celltype_plot <- renderPlot({ NULL })
    output$violin.subcluster <- renderPlot({ NULL })
    output$umap.subcluster <- renderPlot({ NULL })
    output$umap.gene.subcluster <- renderPlot({ NULL })
    
    # get the file path of the chosen dataset
    dataset_path <- file.path("/mnt/s3-bucket/data", input$file)
    
    # check the file extension
    file_extension <- tools::file_ext(dataset_path)
    
    # load the dataset based on the file extension
    if (file_extension == "qs") {
      seurat_obj$data <- qread(dataset_path)  # Load .qs file
    }
    
    # reset loading button
    resetLoadingButton("load_data")

    # get the metadata cols
    metacols <- colnames(seurat_obj$data@meta.data)
    # update the cell type selection input
    updateSelectInput(session, "cell.type.meta", choices = metacols)
    
    # update the group by col in the subcluster tab
    updateSelectInput(session, "meta.subcluster", choices = metacols)
    
    showNotification("Data Loaded!", type = "message")
  })
  
  # handle cancel button ----
  observeEvent(input$cancel_button, {
    resetLoadingButton("upload_data")
    resetLoadingButton("cancel_button")
  })
  
  # handle file upload ----
  observeEvent(input$upload_data, {
    req(input$file_upload)  # ensure a file is uploaded
    
    # Reset old data and metadata UI selections
    seurat_obj$data <- NULL
    seurat_obj$cell_type_column <- NULL
    updateSelectInput(session, "cell.type.meta", choices = character(0), selected = NULL)
    updateSelectInput(session, "cell.type.DE", choices = character(0), selected = NULL)
    updateSelectInput(session, "cell.type.subcluster", choices = character(0), selected = NULL)
    updateSelectInput(session, "meta.subcluster", choices = character(0), selected = NULL)
    output$cell_type_table <- DT::renderDT({ NULL })
    output$de_table <- DT::renderDT({ NULL })
    output$feature_plot <- renderPlot({ NULL })
    output$umap_plot <- renderPlot({ NULL })
    output$violin_celltype_plot <- renderPlot({ NULL })
    output$violin.subcluster <- renderPlot({ NULL })
    output$umap.subcluster <- renderPlot({ NULL })
    output$umap.gene.subcluster <- renderPlot({ NULL })
    
    # get the file path and name of the uploaded file
    file_path <- input$file_upload$datapath
    file_name <- input$file_upload$name
    
    # save the uploaded file to the directory
    saved_file_path <- file.path("/mnt/s3-bucket/data", file_name)
    file.copy(file_path, saved_file_path, overwrite = TRUE)
    
    # check the file extension
    file_extension <- tools::file_ext(saved_file_path)
    
    # load the dataset based on the file extension
    if (file_extension == "qs") {
      seurat_obj$data <- qread(saved_file_path)  # Load .qs file
    } else if (file_extension == "h5Seurat") {
      seurat_obj$data <- LoadH5Seurat(saved_file_path)  # Load .h5Seurat file
    }
    
    # update the files in the drop down to include the updated file
    files_in_dir <- list.files("/mnt/s3-bucket/data", pattern = "\\.qs$", full.names = FALSE)
    updateSelectInput(session, "file", choices = files_in_dir, selected = file_name)
    
    # reset loading button
    resetLoadingButton("upload_data")
    
    # get the metadata cols
    metacols <- colnames(seurat_obj$data@meta.data)
    # update the cell type selection input
    updateSelectInput(session, "cell.type.meta", choices = metacols)
    
    # update the group by col in the subcluster tab
    updateSelectInput(session, "meta.subcluster", choices = metacols)
  
    # update the selectInput choices to include the newly uploaded file
    showNotification("File uploaded successfully!", type = "message")
    
  })
  
  # store cell type column ----
  observeEvent(input$set_cell_type, {
    req(input$cell.type.meta)
    seurat_obj$cell_type_column <- input$cell.type.meta
    
    # update the cell type drop down for the de table
    cell.col <- input$cell.type.meta
    cell_types <- unique(seurat_obj$data@meta.data[[cell.col]])
    updateSelectInput(session, "cell.type.DE", choices = cell_types)
    
    # update the cell type drop down for the subcluster tab
    updateSelectInput(session, "cell.type.subcluster", choices = cell_types)
    
    # update the identity of the seurat obj
    Idents(seurat_obj$data) <- input$cell.type.meta
    
    showNotification("Cell type column set!", type = "message")
    
    # show a summary of the cell types
    output$cell_type_table <- DT::renderDT({
      df <- cell_type_count(seurat_obj$data, cell.col)
      DT::datatable(df, options = list(
        searching = FALSE,   # enable searching
        paging = FALSE,      # enable pagination
        scrollX = TRUE      # allow horizontal scrolling
      ), rownames = FALSE)
    })
  })
  
  # generate umap plots ----
  observeEvent(input$generate_umaps, {
    req(seurat_obj$data) # ensure data is loaded
    req(input$gene) # ensure a gene has been chosen
    
    # isolate the gene selected
    gene_selected <- isolate(input$gene)
    
    # check if the gene exists in the dataset
    if (!(gene_selected %in% rownames(seurat_obj$data))) {
      showNotification("Error: Gene not found in dataset!", type = "error", duration = 3)
      return(NULL)  # stop execution if gene is not found
    }
    
    # render the feature plot
    output$feature_plot <- renderPlot({
      generate_feature_plot(seurat_obj$data, gene_selected)
    })
    
    # render the umap
    output$umap_plot <- renderPlot({
      generate_umap_plot(seurat_obj$data, seurat_obj$cell_type_column)  
    })
    
    # render the violin plot
    output$violin_celltype_plot <- renderPlot({
      gene_across_cell_type(seurat_obj$data, gene_selected, seurat_obj$cell_type_column)  
    })
  })
  
  # generate de table ----
  observeEvent(input$display.de.table, {
    req(seurat_obj$data) # ensure data is loaded
    
    output$de_table <- DT::renderDT({
      df <- generate_de_table(seurat_obj$data, input$cell.type.DE, seurat_obj$cell_type_column)
      DT::datatable(df, options = list(
        searching = TRUE,   # enable searching
        paging = TRUE,      # enable pagination
        filter = 'top',     # show filter inputs above columns
        scrollX = TRUE      # allow horizontal scrolling
      ))
    })
  })
  
  # generate subcluster plots ----
  observeEvent(input$display.subcluster, {
    req(seurat_obj$data) # ensure data is loaded
    req(input$cell.type.subcluster) # ensure a cell type has been chosen
    req(input$gene.subcluster) # ensure a gene has been chosen
    req(input$meta.subcluster) # ensure a group by column has been chosen
    
    # isolate the gene selected
    gene_selected <- isolate(input$gene.subcluster)
    # isolate the group.by col selected
    group_by <- isolate(input$meta.subcluster)
    
    # check if the gene exists in the dataset
    if (!(gene_selected %in% rownames(seurat_obj$data))) {
      showNotification("Error: Gene not found in dataset!", type = "error", duration = 3)
      return(NULL)  # stop execution if gene is not found
    }
    
    # subset the seurat obj
    subset <- subset(x = seurat_obj$data, idents = input$cell.type.subcluster)
    
    # isolate wether non-zero or not
    non_zero <- isolate(input$non.zero)
    
    # render the violin plot
    output$violin.subcluster <- renderPlot({
      if (non_zero == 1){
        generate_violin_plot(subset, gene_selected, group_by)
      }
      else{
        subset_nonzero <- subset(subset, subset = !!sym(gene_selected) > 0)
        generate_violin_plot(subset_nonzero, gene_selected, group_by)
      }
    })
    
    # render the umap for a subcluster
    output$umap.subcluster <- renderPlot({
      generate_subcluster_umap(subset, group_by)
    })
    
    # render the umap for a gene in the subcluster
    output$umap.gene.subcluster <- renderPlot({
      generate_subcluster_featureplot(subset, gene_selected)
    })
  })
}

# run the app ----
shinyApp(ui = ui, server = server)
