# Module TIMELINE PLOT
#
# Creates an intensity vs TIME plot for one feature ID
#
# Input:
#    state$sel_feature_plot - selection of feature IDs to plot
# 
# Uses:
#    DATA$data - a tibble with values
#    DATA$features - a tibble with feature ID and feature name
#    DATA$metadata - a tibble with data grouping
#

require(shinyWidgets)
require(dplyr)
require(tibble)
require(ggplot2)
require(ggbeeswarm)

# ----- UI definitions -----

mod_timeline_plot_ui <- function(id) {
  ns <- NS(id)
  
  # Read data columns from config
  choices <- unlist(CONFIG$data_columns)
  col_choice <- rlang::set_names(names(choices), choices)
  
  tagList(
    dropdownButton(

      selectInput(
        inputId = ns("intensity_value"),
        label = "Intensity value",
        choices = col_choice,
        selected = CONFIG$default_data_column
      ),
      circle = FALSE,
      status = "primary",
      icon = icon("gear"),
      size = "xs",
      tooltip = tooltipOptions(title = "Plot configuration")
    ),
    
    plotOutput(
      outputId = ns("timeline_plot"),
      width = "100%"
    )
  )
  
}

# ----- Server logic -----

mod_timeline_plot_server <- function(id, state) {
  
  server <- function(input, output, session) {
    
    output$timeline_plot <- renderPlot({
      ids <- state$sel_feature_plot
      req(ids, length(ids) == 1)
      d <- DATA$data |> 
        dplyr::filter(id %in% ids) |> 
        dplyr::left_join(DATA$features, by = "id") |> 
        dplyr::left_join(DATA$metadata, by = "sample") 
      sh_plot_timeline(d)
    })
  }
  
  moduleServer(id, server)
}


sh_plot_timeline <- function(d, what = "abu_norm", min_n = 3, exclude_days = 4, dodge_width = 1,
                             text_size = 18, point_size = 2) {
  dm <- d |> 
    mutate(val = get(what)) |> 
    filter(!(day %in% exclude_days)) |> 
    select(day, treatment, val) |> 
    mutate(day = as.numeric(as.character(day))) |> 
    group_by(treatment, day) |> 
    summarise(
      m = mean(val),
      s = sd(val),
      n = n(),
      se = s / sqrt(n),
      t_crit = qt(0.975, df = n - 1),
      ci = se * t_crit
    ) |> 
    ungroup() |> 
    filter(n >= min_n)
  days <- unique(dm$day) |> sort()
  pd <- position_dodge(width = dodge_width)
  clrs <- c("#E69F00", "#56B4E9")
  dm |> 
    ggplot(aes(x = day, y = m, ymin = m - ci, ymax = m + ci, colour = treatment, fill = treatment)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "top",
      text = element_text(size = text_size)
    ) +
    geom_line(alpha = 0.5, position = pd) +
    geom_errorbar(position = pd, width = 1) +
    geom_point(position = pd, shape = 21, colour = "grey30", size = point_size) +
    scale_colour_manual(values = clrs, name = "Treatment") +
    scale_fill_manual(values = clrs, name = "Treatment") +
    scale_x_continuous(breaks = days) +
    labs(x = "Day", y = "Normalised abundance")
}

