# Module for global input for all other modules

require(dqshiny)

# ----- UI definitions -----

mod_global_input_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    selectInput(
      inputId = ns("contrast"),
      label = "Contrast",
      choices = c("")
    ),
    autocomplete_input(
      id = ns("search"),
      label = "Search",
      options = NULL
    ),
    actionButton(
      inputId = ns("clear"),
      label = "Clear"
    )
  )
}

# ----- Server logic -----

mod_global_input_server <- function(id, state) {
  
  server <- function(input, output, session) {
    
    # Update dummy contrast selection
    contrasts <- levels(DATA$de$contrast)
    observe({
      updateSelectInput(session, "contrast", choices = contrasts)
    })
    
    # Update search input with feature names
    all_names <- DATA$features$name |> unique()
    observe({
      update_autocomplete_input(session, "search", options = c("", all_names))
    })
    
    # Observe contrast selection, copy to state
    observeEvent(input$contrast, state$contrast <- input$contrast)
    
    # Observe search input, copy to state, convert "" to NULL
    observeEvent(input$search, {
      state$search <- input$search;
      if(state$search == "") 
        state$search <- NULL
    },
    ignoreNULL = FALSE
    )
    
    # Clear input
    observeEvent(input$clear, update_autocomplete_input(session, "search", value = "", placeholder = ""))
    
  }
  
  moduleServer(id, server)
}
