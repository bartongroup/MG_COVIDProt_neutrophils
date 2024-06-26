# Module FEATURE PLOT
#
# Creates an intensity vs group plot for one feature ID or a entity vs ID vs
# intensity heatmap for multiple feature IDs,
#
# Input:
#    state$sel_feature_plot - selection of feature IDs to plot
# 
# Uses:
#    DATA$data - a tibble with values per feature
#    DATA$features - a tibble with feature ID and feature name
#    DATA$metadata - a tibble with data grouping
#

require(shinyWidgets)
require(dplyr)
require(tibble)
require(ggplot2)
require(ggbeeswarm)

# ----- UI definitions -----

mod_feature_plot_ui <- function(id) {
  ns <- NS(id)
  
  # Read data columns from config
  choices <- unlist(CONFIG$data_columns)
  col_choice <- rlang::set_names(names(choices), choices)
  
  tagList(
    dropdownButton(
      
      checkboxInput(
        inputId = ns("norm_fc"),
        label = "Heatmap normalised per row"
      ),
      
      conditionalPanel(
        condition = 'input.norm_fc == 0',
        ns = ns,
        
        selectInput(
          inputId = ns("intensity_value"),
          label = "Intensity value",
          choices = col_choice,
          selected = CONFIG$default_data_column
        ),
        
        radioButtons(
          inputId = ns("intensity_scale"),
          label = "Intesity scale",
          choices = c("Linear" = "lin", "Logarithmic" = "log"),
          inline = TRUE
        )
      ),
      
      circle = FALSE,
      status = "primary",
      icon = icon("gear"),
      size = "xs",
      tooltip = tooltipOptions(title = "Plot configuration")
    ),
    
    plotOutput(
      outputId = ns("feature_plot"),
      width = "100%"
    )
  )
  
}

# ----- Server logic -----

mod_feature_plot_server <- function(id, state) {
  
  server <- function(input, output, session) {
    
    output$feature_plot <- renderPlot({
      ids <- state$sel_feature_plot
      req(ids)
      d <- DATA$data |> 
        dplyr::filter(id %in% ids) |> 
        dplyr::left_join(DATA$features, by = "id")
      sh_plot_features(d, meta = DATA$metadata, scale = input$intensity_scale,
                    what = input$intensity_value, max_n_lab = 50, norm_fc = input$norm_fc)
    })
  }
  
  moduleServer(id, server)
}



#' Plot one feature
#'
#' @param d Tibble with feature intensities. Columns needed: group, val, replicate
#' @param ylab Label on y axis
#' @param scale Scale of y axis (lin or log)
#' @param text_size Text size
#' @param point_size Point size
#' @param cex Point spread scaling for beeswarm
#'
#' @return ggplot object
sh_plot_one_feature <- function(d, ylab, scale = c("lin", "log"), text_size, point_size, cex) {
  
  d <- d |>
    dplyr::mutate(shape = if_else(val == 0, 24, 21))
  
  ncond <- length(unique(d$group))
  vlines <- tibble::tibble(x = seq(0.5, ncond + 0.5, 1))
  
  g <- ggplot() +
    theme_bw() +
    theme(
      text = element_text(size = text_size),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    scale_shape_identity() +  # necessary for shape mapping
    geom_beeswarm(data = d, aes(x = group, y = val, fill = colour, shape = shape),
                  colour = "grey40", size = point_size, cex = cex) +
    geom_vline(data = vlines, aes(xintercept = x), colour = "grey80", alpha = 0.5) +
    scale_fill_viridis_d(option = "cividis") +
    labs(x = NULL, y = ylab)
  
  #if(scale == "lin")
  #  g <- g + scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA))
  
  return(g)
}

#' Plot a heatmap of multiple features
#'
#' @param d Feature data
#' @param ylab Label on y axis
#' @param text_size Text size
#' @param max_n_lab Limit of features above which feature names are not displayed
#' @param norm_fc Logical, normalise in each feature to its mean and plot logFC
#'
#' @return ggplot object
sh_plot_feature_heatmap <- function(d, lab, text_size, max_n_lab, norm_fc) {
  
  if (norm_fc) {
    d <- d |>
      dplyr::mutate(val = get(CONFIG$default_data_column)) |> 
      dplyr::group_by(id) |> 
      dplyr::mutate(M = mean(val, na.rm = TRUE)) |> 
      dplyr::mutate(val = log2(val / M)) |> 
      dplyr::ungroup()
    lab <- expression(log[2]~FC)
  }
  
  g <- d |> 
    ggplot(aes(x = entity, y = name, fill = val)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      text = element_text(size = text_size),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    geom_tile() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL, fill = lab)
  
  if(norm_fc) {
    g <- g + scale_fill_distiller(type = "div", palette = "RdBu")
  } else {
    g <-  g + scale_fill_viridis_c(option = "cividis")
  }
  
  if(length(unique(d$id)) > max_n_lab)
    g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  return(g)
}

#' Make a feature plot: intensity vs entity or group
#'
#' @param dat Tibble with feature intensities
#' @param meta Metadata with grouping of entities
#' @param scale Intensity scale, "lin" or  "log"
#' @param what Which column to plot
#' @param text_size Text size
#' @param point_size Point size
#' @param cex Point spread scaling for beeswarm
#' @param max_n_lab Limit of features above which feature names are not displayed
#' @param norm_fc Logical, normalise in each feature to its mean and plot logFC
#'
#' @return A ggplot object
sh_plot_features <- function(dat, meta, scale, what = "rpkm", text_size = 14, point_size = 3, cex = 2,
                          max_n_lab = 30, norm_fc = FALSE) {
  
  if(nrow(dat) == 0) return(NULL)
  
  d <- dat |> 
    dplyr::mutate(val = get(what)) |> 
    dplyr::left_join(meta, by = "entity")
  
  lab <- what
  
  if(scale == "log"){
    d$val <- log10(d$val)
    lab <-  str_glue("Log {what}")
  }
  
  n_feat <- length(unique(d$id)) 
  if(n_feat == 1) {
    sh_plot_one_feature(d, lab, scale, text_size, point_size, cex)
  } else {
    sh_plot_feature_heatmap(d, lab, text_size, max_n_lab, norm_fc)
  }
}