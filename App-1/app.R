library(shiny)
library(bslib)
library(Seurat)
library(SeuratDisk)
library(DT)
# source other R scripts
source('helpers.R')
# define the path to the data directory
data_dir <- "data/"

# define the ui ----
ui <- fluidPage(
  titlePanel("Goessling Lab Single Cell Browser"),
  tabsetPanel(
    
    br(), # page break
    
    # data import tab ----
    tabPanel("Import Data",
      selectInput("file", "Choose a dataset:", choices = list.files(data_dir, pattern = "\\.h5Seurat$", full.names = FALSE),
                  selected = NULL), # default empty
      actionButton("load_data", "Load Dataset"),
      
      br(), br(), # page break
      
      selectInput("cell.type.meta", "Choose cell type column from metadata", choices = NULL), # default empty
      actionButton("set_cell_type", "Set Cell Type Column") 
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
          selectInput("cell.type.DE", "Choose a cell type", choices = NULL), # default empty
          actionButton("display.de.table", "Display Table")
        ),
        mainPanel(
          DT::dataTableOutput("de_table")
        )
      )),
    
    tabPanel("Cell Type Specific UMAP"),
    
    tabPanel("Cell Type Specific Differential Expression")
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
    # get the file path of the chosen dataset
    dataset_path <- file.path("data/", input$file)
    # load the seurat object
    seurat_obj$data <- LoadH5Seurat(dataset_path)
    # get the metadata cols
    metacols <- colnames(seurat_obj$data@meta.data)
    # update the cell type selection input
    updateSelectInput(session, "cell.type.meta", choices = metacols)
    showNotification("Data Loaded!", type = "message")
  })
  
  # store cell type column ----
  observeEvent(input$set_cell_type, {
    req(input$cell.type.meta)
    seurat_obj$cell_type_column <- input$cell.type.meta
    
    # update the cell type drop down for the de table
    cell.col <- input$cell.type.meta
    cell_types <- unique(seurat_obj$data@meta.data[[cell.col]])
    updateSelectInput(session, "cell.type.DE", choices = cell_types)
    
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
}

# run the app ----
shinyApp(ui = ui, server = server)
