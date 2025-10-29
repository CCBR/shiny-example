library(shiny)
library(bslib)

# Load the plot scripts
source("R/violin_boxplot.R")

# Load the input data
tutorial.folder <- getwd()
data.folder <- paste0(tutorial.folder, "/test_datasets/Human_Kidney_DSP/")
counts.data <- read.csv(paste0(data.folder, "human_kidney_test_quantile_norm_counts.csv"))
annotation.data <- read.csv(paste0(data.folder, "human_kidney_test_annotation.csv"))

# Define the column names to offer
annotation.cols <- colnames(annotation.data %>% 
                              select_if(is.character))

# Reformat the counts and annotation data for mapping and transforming
rownames(counts.data) <- counts.data$gene
counts.data <- counts.data %>% select(-gene)
colnames(counts.data) <- sub("\\.dcc", "", colnames(counts.data))
colnames(counts.data) <- gsub("\\.", "-", colnames(counts.data))

rownames(annotation.data) <- sub("\\.dcc", "", annotation.data$sample_ID)

# Define UI
ui <- page_sidebar(
  # App title 
  title = "DSP Analysis - QC Plots",
  # Sidebar panel for inputs
  sidebar = sidebar(
    
    selectizeInput("annotation", 
                   "Annotation to display", 
                   choices = NULL, 
                   multiple = FALSE), 
    
    selectizeInput("gene", 
                   "Gene to display", 
                   choices = NULL, 
                   multiple = FALSE), 
    
    checkboxInput("summary_stats", 
                  "Display Summary Stats", 
                  value = FALSE)
    
  ),
  
  # Tabs in the app
  navset_card_underline(
    nav_panel("Boxplot", 
              plotOutput("Boxplot")))
  
  
  
)


# Define server logic required to create the plot
server <- function(input, output, session) {
  
  # Reactive annotation selection
  observe({
    updateSelectizeInput(session, 
                         "annotation", 
                         choices = annotation.cols, 
                         server = TRUE, 
                         selected = "slide_name")
    
    
    
  })
  
  observe({
    
    updateSelectizeInput(session, 
                         "gene", 
                         choices = rownames(counts.data), 
                         server = TRUE)
    
  })
  
  summary_stats <- reactive({input$summary_stats})
  
  output$Boxplot <- renderPlot({
    
    violin_boxplot(counts = counts.data,  
                                   annotation.df = annotation.data, 
                                   gene.list = input$gene, 
                                   annotation.field = input$annotation, 
                                   display.summary.stat = input$summary_stats, 
                                   compare = FALSE, 
                                   compare.groups = NULL)
  
  })
    
    
}

shinyApp(ui = ui, server = server)
