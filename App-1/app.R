library(shiny)
library(bslib)
library(Seurat)
library(SeuratDisk)
library(DT)
library(rlang)
library(shinythemes)
library(shinyFeedback)
library(rsconnect)
library(qs)

# source other R scripts
source('helpers.R')
# define the path to the data directory
data_dir <- "data"
# increase the file upload size allowed
options(shiny.maxRequestSize = 4000 * 1024^2)


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
       fluidRow(
         column(6, # left column for file upload and select input
                fileInput("file_upload", "Upload dataset:", accept = c(".qs")),
                div(style = "margin-top: -30px"), # Adjust space between elements
                loadingButton("upload_data", "Process Dataset for Analysis", style = "width: 71%; margin: 10px auto;"),
                
                selectInput("file", "OR choose a dataset:", choices = list.files(data_dir, pattern = "\\.qs$", full.names = FALSE),
                            selected = NULL), # default empty
                loadingButton("load_data", "Load & Process Dataset", style = "width: 71%; margin: 10px auto;")
         ),
         column(6, # right column for metadata cell type selection
                selectInput("cell.type.meta", "Choose the metadata cell type column:", choices = NULL), # default empty
                actionButton("set_cell_type", "Set Cell Type Column")
         )
       ),
    ),
    
    # umap tab: all cell types -----
    tabPanel("UMAP",
      sidebarLayout(
        sidebarPanel(
          textInput("gene", "Enter a gene of interest:", placeholder = "Type gene name..."),
          actionButton("generate_umaps", "Generate UMAPs") 
          ),
          mainPanel(
            plotOutput("feature_plot"),  
            plotOutput("umap_plot")     
        )
      )
    ),
    
    # differential expression table: all cell types ----
    tabPanel("Differential Expression",
      sidebarLayout(
        sidebarPanel(
          selectInput("cell.type.DE", "Choose a cell type:", choices = NULL), # default empty
          actionButton("display.de.table", "Display Table")
        ),
        mainPanel(
          DT::dataTableOutput("de_table")
        )
      )),
    
    tabPanel("Cell Type Specific",
      sidebarLayout(
        sidebarPanel(
          selectInput("cell.type.subcluster", "Choose a cell type:", choices = NULL), # default empty
          textInput("gene.subcluster", "Enter a gene of interest:", placeholder = "Type gene name..."),
          radioButtons("non.zero", "Display cells with non-zero expression of selected gene?", choices = list("Yes" = 1, "No" = 2),
                       selected = 1),
          selectInput("meta.subcluster", "Choose a metadata column to group by:", choices = NULL), # default empty
          actionButton("display.subcluster", "Submit")
        ),
        mainPanel(
          plotOutput("umap.subcluster"),  
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
    
    if (!dir.exists("data")) {
      dir.create("data", recursive = TRUE)
    }
    
    # get the file path of the chosen dataset
    dataset_path <- file.path("data/", input$file)
    
    # load the seurat object
    seurat_obj$data <- qread(dataset_path)
    resetLoadingButton("load_data")

    # get the metadata cols
    metacols <- colnames(seurat_obj$data@meta.data)
    # update the cell type selection input
    updateSelectInput(session, "cell.type.meta", choices = metacols)
    
    # update the group by col in the subcluster tab
    updateSelectInput(session, "meta.subcluster", choices = metacols)
    
    showNotification("Data Loaded!", type = "message")
  })
  
  # handle file upload ----
  observeEvent(input$upload_data, {
    req(input$file_upload)  # ensure a file is uploaded
    
    if (!dir.exists("data")) {
      dir.create("data", recursive = TRUE)
    }
    
    # get the file path and name of the uploaded file
    file_path <- input$file_upload$datapath
    file_name <- input$file_upload$name
    
    # save the uploaded file to the directory
    saved_file_path <- file.path(data_dir, file_name)
    file.copy(file_path, saved_file_path, overwrite = TRUE)
    
    # load the seurat object
    seurat_obj$data <- qread(saved_file_path)
    resetLoadingButton("upload_data")
    
    # get the metadata cols
    metacols <- colnames(seurat_obj$data@meta.data)
    # update the cell type selection input
    updateSelectInput(session, "cell.type.meta", choices = metacols)
    
    # update the group by col in the subcluster tab
    updateSelectInput(session, "meta.subcluster", choices = metacols)
  
    # Update the selectInput choices to include the newly uploaded file
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
  })
}

# run the app ----
shinyApp(ui = ui, server = server)
