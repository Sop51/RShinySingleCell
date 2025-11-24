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
    id = "tabs",
    # how to tab ----
    tabPanel("How To",
             
      # title for how to video ----
      h3("How to Video",
          style = "font-weight: bold; font-size: 20px; color: #1E3A8A; 
            text-align: center; margin-top: 20px; margin-bottom: 15px;"),
      
      # youtube video ----
      div(
        style = "position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden;",
        tags$iframe(
          src = "https://www.youtube.com/embed/xBjEfo02og8?si=-VK404u21xj4-2hW",
          style = "position:absolute; top:0; left:0; width:100%; height:100%;",
          frameborder = "0",
          allow = "accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture",
          allowfullscreen = NA
        )
      ),
      
      # upload data title ----
      h3("Thinking of uploading data? Please review the requirements:",
         style = "font-weight: bold; font-size: 20px; color: #1E3A8A; 
            text-align: center; margin-top: 20px; margin-bottom: 15px;"),
      tags$ul(
        style = "font-size: 15px; color: #555; text-align: left; 
           margin: 0 auto; width: 80%; line-height: 1.4;",
        tags$li("Dataset must be a Seurat object created using Seurat v5.3.0."),
        tags$li("The object must be saved and uploaded as a .qs file."),
        tags$li("Cell type annotations must be present in the meta.data slot and contain valid labels."),
        tags$li("Raw counts or normalized expression values must be included."),
        tags$li("Dimensionality reductions (e.g., PCA, UMAP) must already be computed."),
      ),
      
      # title for contact/github info ----
      h3("Need to get in contact with us?",
         style = "font-weight: bold; font-size: 20px; color: #1E3A8A; 
            text-align: center; margin-top: 20px; margin-bottom: 15px;"),
      
      p("Please reach out to smarcotte2@mgh.harvard.edu",
        style = "font-size: 15px; color: #555; text-align: center; margin-bottom: 20px;"),
      
      div(
        style = "text-align: center; margin-top: 20px;",
        tags$p(
          "Repo link:",
          style = "font-size: 14px; color: #555; margin-bottom: 5px;"
        ),
        tags$a(
          href = "https://github.com/Sop51/RShinySingleCell/",
          target = "_blank",
          tags$i(class = "fa fa-github", style = "margin-right: 6px;"),
          "RShinySingleCell"
        )
      )
      
    ),
    # data import tab ----
    tabPanel("Import Data",
             
       # title and description above everything ----
       h3("Select and Load a Dataset", 
          style = "font-weight: bold; font-size: 20px; color: #1E3A8A; text-align: center; margin-top: 20px; margin-bottom: 15px;"),
       
       p("To begin, select a dataset from the list below. Wait for the data to load and process. Then, set the appropriate cell type annotations within the metadata.",
         style = "font-size: 15px; color: #555; text-align: center; margin-bottom: 20px;"),
       
       # line break here
       hr(style = "border: 1px solid #ccc; margin-top: 20px; margin-bottom: 20px;"), 
       
       # fix wrapping within the file drop down list for long file names ----
       tags$head(
         tags$style(HTML("
           .selectize-dropdown-content .option,
           .selectize-dropdown-content .item {
             white-space: normal !important;
             word-wrap: break-word !important;
           }

           .selectize-input > div {
             white-space: normal !important;
             word-wrap: break-word !important;
             overflow-wrap: break-word !important;
             max-width: 100% !important;   
           }

           .selectize-input {
             height: auto !important;
             min-height: 38px;
             width: 100% !important;  
           }
        "))
       ),
       
       # main layout ----
       fluidRow(
         column(6,
                selectInput("file", "1. Choose a dataset:",
                            choices = list.files("/home/single-cell-shiny-app/data", pattern = "\\.qs$"),
                            selected = NULL),
                loadingButton("load_data", "Load & Process Dataset",
                              style = "width: 71%; margin: 10px auto;")
         ),
         column(6,
                selectInput("cell.type.meta", "2. Choose the cell type annotation label:", choices = NULL),
                actionButton("set_cell_type", "Set Cell Annotations"),
                tags$i(class = "fa fa-info-circle", 
                       id = "info_cell_col", 
                       style = "margin-left: 5px; cursor: pointer; font-size: 16px; color: #007BFF;"),
                DT::dataTableOutput("cell_type_table"),
                
                # Tooltip for the radio button info icon
                bsTooltip("info_cell_col", 
                          "Choose the label in your single cell object that represents the cell type annotations",
                          placement = "right", trigger = "hover")
         )
       ),
    ),
    
    # umap tab: all cell types -----
    tabPanel("UMAP",
       # title and description above everything ----
       h3("Gene Expression Visualization with UMAPs", 
          style = "font-weight: bold; font-size: 20px; color: #1E3A8A; text-align: center; margin-top: 20px; margin-bottom: 15px;"),
       
       p("Enter a gene of interest to visualize its expression across different cell types.",
         style = "font-size: 15px; color: #555; text-align: center; margin-bottom: 30px;"),
       
       # line break here
       hr(style = "border: 1px solid #ccc; margin-top: 20px; margin-bottom: 20px;"), 
       
       # main layout - sidebar and main panel ----      
      sidebarLayout(
        sidebarPanel(
          textInput("gene", "Enter a gene of interest (case sensitive):", placeholder = "Type gene name..."),
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
    tabPanel("Marker Gene Analysis",
       # title and description above everything ----
       h3("Differential Expression Between Cell Types", 
          style = "font-weight: bold; font-size: 20px; color: #1E3A8A; text-align: center; margin-top: 20px; margin-bottom: 15px;"),
       
       p("Compare gene expression in the selected cell type to all other cell types in the dataset.",
         style = "font-size: 15px; color: #555; text-align: center; margin-bottom: 30px;"),
       
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
     
     p("Explore gene expression patterns for specific cell types. 
        Choose a cell type, enter a gene of interest, and customize the data groupings to visualize.",
       style = "font-size: 15px; color: #555; text-align: center; margin-bottom: 20px;"),
     
     # line break here
     hr(style = "border: 1px solid #ccc; margin-top: 20px; margin-bottom: 20px;"), 
     
     # main layout - sidebar and main panel ----
      sidebarLayout(
        sidebarPanel(
          selectInput("cell.type.subcluster", "Choose a cell type:", choices = NULL), # default empty
          textInput("gene.subcluster", "Enter a gene of interest (case sensitive):", placeholder = "Type gene name..."),
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
            selectInput("meta.subcluster", "Choose a metadata label to group by:", choices = NULL), # default empty
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
    dataset_path <- file.path("/home/single-cell-shiny-app/data", input$file)
    
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
  
  # control tab visibility ----
  data_loaded <- reactiveVal(FALSE)
  cell_type_set <- reactiveVal(FALSE)
  
  # check if data is loaded
  observeEvent(input$load_data, {
    data_loaded(TRUE)
  })
  
  # check if data loaded AND cell type col is set
  observeEvent(input$set_cell_type, {
    if (data_loaded()) {
      cell_type_set(TRUE)
    }
  })
  
  # prevent switching tabs until both dataset is loaded and column is set 
  observeEvent(input$tabs, {
    
    # allowed tabs before required data is set
    allowed_tabs_before_ready <- c("How To", "Import Data")
    
    if (input$tabs %in% allowed_tabs_before_ready) return()
    
    # block if BOTH are not set
    if (!data_loaded() && !cell_type_set()){
      showNotification(
        "Please load a dataset and set the cell type column.",
        type = "error", duration = 4
      )
      updateTabsetPanel(session, "tabs", selected = "Import Data")
      return()
    }
    
    # block if ONLY data not loaded
    if (!data_loaded()) {
      showNotification(
        "Please import and load a dataset.",
        type = "error", duration = 4
      )
      updateTabsetPanel(session, "tabs", selected = "Import Data")
      return()
    }
    
    # block if ONLY cell type not set
    if (!cell_type_set()) {
      showNotification(
        "Please set the cell type annotation label.",
        type = "error", duration = 4
      )
      updateTabsetPanel(session, "tabs", selected = "Import Data")
      return()
    }
  })
}


# run the app ----
shinyApp(ui = ui, server = server)
