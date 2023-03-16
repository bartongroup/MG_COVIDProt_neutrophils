# GSE table
# Creates a table with GSE results

# ----- UI definitions -----

mod_gse_table_ui <- function(id) {
  ns <- NS(id)
  tagList(
    dropdownButton(
      radioGroupButtons(
        inputId = ns("ontology"),
        label = "Ontology",
        choices = c(
          "GO" = "go",
          "Reactome" = "re",
          "KEGG" = "kg",
          "BioPlanet" = "bp"
        )
      ),
      numericInput(
        inputId = ns("fdr_limit"),
        label = "FDR limit",
        value = 0.05,
        min = 0,
        max = 1
      ),
      circle = FALSE,
      status = "primary",
      icon = icon("gear"),
      size = "xs",
      tooltip = tooltipOptions(title = "Table configuration")
    ),
    
    DT::dataTableOutput(
      outputId = ns("gse_table")
    )
  )
  
}

# ----- Server logic -----

mod_gse_table_server <- function(id, state) {
  
  server <- function(input, output, session) {
    
    # Prepare updated feature info table
    make_gse_table <- reactive({
      ctr <- state$contrast
      req(ctr)
      DATA$gse[[input$ontology]] |>
        filter(contrast == ctr & fdr <= input$fdr_limit) |> 
        select(TermId = term_id, Name = term_name, nes) |> 
        mutate(across(where(is.numeric), ~signif(.x, 2))) |> 
        mutate(TermId = as.character(TermId))
    })
    
    # Wait for row selection in the GSE, pass it to app state
    observeEvent(input$gse_table_rows_selected, ignoreNULL = FALSE, {
      rows_sel <- input$gse_table_rows_selected
      ids <- NULL
      if(!is.null(rows_sel)) {
        gt <- make_gse_table()
        termid <- gt[rows_sel, ]$TermId
        r <- DATA$gse[[input$ontology]] |> 
          filter(term_id == termid)
        if(nrow(r) > 0) {
          ids <- r$leading_edge[[1]]
        }
      }
      state$sel_term <- ids
    })
    
    # feature_info table
    output$gse_table <- DT::renderDataTable({
      gt <- make_gse_table()
      DT::datatable(
        gt,
        options = list(paging = FALSE),
        style = "bootstrap",
        selection = "single",
        rownames = FALSE
      ) |> 
        DT::formatStyle(columns = colnames(gt), fontSize = '80%')
    })
  }
  
  moduleServer(id, server)
}
