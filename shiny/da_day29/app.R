### EDGER EXPLORER

libDir <- "/cluster/gjb_lab/mgierlinski/R_shiny/library/4.1"
if (dir.exists(libDir)) .libPaths(libDir)

library(shiny)
library(tidyverse)
library(DT)
source("../shiny_func.R")

select <- dplyr::select
filter <- dplyr::filter

css <- "table{font-size: 11px; background-color: #EAF5FF}"

### Read data ###

cat("Reading data")
set <- read_rds("../data.rds")
da <- read_rds("../da_day29.rds") |> 
  mutate(sel = FDR < 0.05)
genes <- read_rds("../genes.rds")
all_terms <- read_rds("../terms.rds")
cat("\n")

all_genes <- da$gene_symbol |> unique()
contrasts <- levels(da$contrast)
max_points <- 500



#######################################################################

ui <- shinyUI(fluidPage(
  
  tags$style(css),
  
  titlePanel("COVID-19 proteomics: day 29 brensocatib vs placebo"),

  fluidRow(
    column(12,
      fluidRow(
        column(4,
          radioButtons("contrast_sel", "Contrast:", choices = contrasts, inline = TRUE),
          radioButtons("plot_type", "Plot type:", choices = c("MA", "Volcano"), inline = TRUE),
          plotOutput("main_plot", height = "480px", width = "100%", brush = "plot_brush", hover = "plot_hover")
        ),
        column(3,
          plotOutput("protein_plot", height = "300px",width = "100%")
        ),
        column(5,
          p("Protein list"),
          div(style = 'height: 200px; overflow-y: scroll', tableOutput("protein_info")),
          br(),
          radioButtons("enrichment_sel", "Enrichment:",
                       choiceNames = c("GO", "Reactome", "KEGG"),
                       choiceValues = c("go", "re", "kg"),
                       inline = TRUE
          ),
          div(style = 'height: 400px; overflow-y: scroll', tableOutput("enrichment")),
        )
      ),
      fluidRow(
        DT::dataTableOutput("all_protein_table")
      )
    )
  )
)
)


########################################################################################


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  get_xy_data <- function() {
    if (input$plot_type == "Volcano") {
      xy_data <- da |> 
        filter(contrast == input$contrast_sel) |> 
        mutate(x = logFC, y = -log10(PValue))
    } else if (input$plot_type == "MA") {
      xy_data <- da |> 
        filter(contrast == input$contrast_sel) |> 
        mutate(x = AveExpr, y = logFC)
    } else {
      stop("Cowabunga!")
    }
    xy_data
  }
  
  select_protein <- function(max_hover=1) {
    xy_data <- get_xy_data()
    sel <- NULL
    tab_idx <- as.numeric(input$all_gene_table_rows_selected)
    if (!is.null(input$plot_brush)) {
      brushed <- na.omit(brushedPoints(xy_data, input$plot_brush))
      sel <- brushed$id
    } else if (!is.null(input$plot_hover)) {
      near <- nearPoints(xy_data, input$plot_hover, threshold = 20, maxpoints = max_hover)
      sel <- near$id
    } else if (length(tab_idx) > 0) {
      sel <- xy_data[tab_idx, ] |> pull(id)
    }
    return(sel)
  }
  
  output$protein_info <- renderTable({
    sel <- select_protein()
    df <- NULL
    if (!is.null(sel) && length(sel) >= 1 && length(sel) <= max_points) {
      df <- set$info |>
        filter(id %in% sel) |> 
        arrange(gene_names) |> 
        select(gene_names, protein_descriptions)
    } else if (length(sel) > max_points) {
      df <- data.frame(Error = paste0('only ',max_points,' points can be selected.'))
    }
    df
  })

  enrichment_table <- function(terms) {
    xy_data <- get_xy_data()
    sel <- NULL
    fe <- NULL
    if (!is.null(input$plot_brush)) {
      brushed <- na.omit(brushedPoints(xy_data, input$plot_brush))
      sel <- brushed$id
      n <- length(sel)
      if (n > 0 && n <= max_points) {
        sel_genes <- tibble(id = sel) |> 
          left_join(xy_data, by = "id") |> 
          pull(gene_symbol) |> 
          unique()
        fe <- sh_enrichment(all_genes, sel_genes, terms)
      } else if (n > 0) {
        fe <- data.frame(Error = paste0('only ',max_points,' points can be selected.'))
      }
    }
    fe
  }
  
  output$enrichment <- renderTable({
    enrichment_table(all_terms[[input$enrichment_sel]])
  })
  
  
  output$protein_plot <- renderPlot({
    sel <- select_protein()
    if (!is.null(sel) && length(sel) == 1) {
      sh_plot_protein(set, sel)
    }
  })
  
  output$main_plot <- renderPlot({
    xy_data <- get_xy_data()
    tab_idx <- as.numeric(input$all_protein_table_rows_selected)

    if (input$plot_type == "Volcano") {
      g <- sh_plot_volcano(xy_data)
    } else if (input$plot_type == "MA") {
      g <- sh_plot_ma(xy_data)
    }
    if (length(tab_idx) >= 1) {
      g <- g + geom_point(data = xy_data[tab_idx, ], colour = "red", size = 2)
    }
    g
  })

  output$all_protein_table <- DT::renderDataTable({
    d <- get_xy_data() |>
        select(id, gene_symbol, logFC, FDR) |> 
        mutate(across(c(logFC, FDR), ~signif(.x, 3))) |> 
        left_join(select(set$info, id, protein_descriptions), by = "id")
    DT::datatable(d, class = 'cell-border strip hover', selection = "single", rownames = FALSE)
  })
}

# Run the application
shinyApp(ui = ui, server = server)




# input <- list(contrast_sel = "treatmentdrug", plot_type = "MA")
