#' Launch an Interactive ROI Annotation App for Seurat Spatial Data
#'
#' This function opens a Shiny app to allow users to manually annotate Regions of Interest (ROIs)
#' on spatial transcriptomics or metabolomics data stored in a Seurat object. Users can interactively
#' select regions using lasso selection, assign custom names to each ROI, and save the results as
#' binary (0/1) metadata columns in the Seurat object.
#'
#' @param seurat_obj A \code{Seurat} object containing spatial coordinates and metadata.
#' @param image A character string specifying the spatial image source used by \code{GetTissueCoordinates()}.
#'        Default is \code{"fov"}. Set to \code{"VisiumV2"} for Visium-style spatial data.
#'
#' @return A modified \code{Seurat} object with added metadata columns for each saved ROI (1 = inside ROI, 0 = outside).
#'
#' @details
#' The app includes options to:
#' \itemize{
#'   \item Choose a metadata column to display.
#'   \item Adjust the spot size for display.
#'   \item Use lasso to select spots and name each ROI.
#'   \item Save each ROI to the Seurat object metadata.
#'   \item Export the final Seurat object upon clicking "Finish".
#' }
#'
#' Spatial selection uses \code{sf} geometry tools. The plot is rendered using \code{plotly}.
#'
#' @importFrom Seurat GetTissueCoordinates
#' @importFrom shiny shinyApp fluidPage titlePanel sidebarLayout sidebarPanel mainPanel
#' @importFrom shiny selectInput radioButtons sliderInput textInput actionButton verbatimTextOutput
#' @importFrom shiny plotOutput renderPlot reactiveVal runApp showNotification stopApp
#' @importFrom plotly plotlyOutput renderPlotly plot_ly add_trace layout event_data
#' @importFrom sf st_polygon st_sfc st_make_valid st_as_sf st_within
#' @importFrom dplyr %>%
#' @export
SelectROIs <- function(seurat_obj, image = "fov") {

  if (!inherits(seurat_obj, "Seurat")) stop("seurat_obj must be a Seurat object")

  coords <- GetTissueCoordinates(seurat_obj, image = image)
  if (!all(c("x", "y") %in% colnames(coords))) stop("Spatial coordinates 'x' and 'y' not found")

  meta_cols <- colnames(seurat_obj@meta.data)

  return(runApp(shinyApp(
    ui = fluidPage(
      titlePanel("Select Regions of Interest"),
      sidebarLayout(
        sidebarPanel(
          selectInput("meta_col", "Select Metadata Column to Plot", choices = meta_cols),
          radioButtons("plot_type", "Plot Type",
                       choices = c("Feature (Continuous)" = "feature", "Categorical" = "categorical")),
          sliderInput("pt_size", "Spot Size", min = 1, max = 20, value = 5),
          textInput("roi_name", "Name for ROI", value = "ROI_1"),
          actionButton("save_roi", "Save ROI"),
          actionButton("reset_roi", "Reset ROI Selection"),
          actionButton("done", "Finish & Return Object"),
          verbatimTextOutput("status")
        ),
        mainPanel(
          plotlyOutput("spatial_plot", height = "600px")
        )
      )
    ),
    server = function(input, output, session) {
      rv <- reactiveVal(seurat_obj)
      roi_mask <- reactiveVal(rep(0, nrow(coords)))

      output$spatial_plot <- renderPlotly({
        meta_col <- input$meta_col
        plot_type <- input$plot_type
        meta_data <- rv()@meta.data[[meta_col]]
        spot_size <- input$pt_size

        plot_ly() %>%
          add_trace(
            x = coords$x,
            y = coords$y,
            type = "scattergl",
            mode = "markers",
            marker = if (plot_type == "feature") {
              list(color = meta_data, colorscale = "Viridis", showscale = TRUE, size = spot_size)
            } else {
              list(color = as.factor(meta_data), showscale = FALSE, size = spot_size)
            },
            text = rownames(coords),
            hoverinfo = "text"
          ) %>%
          layout(
            title = paste("Spatial Plot:", meta_col),
            xaxis = list(title = "X Coordinate"),
            yaxis = list(title = "Y Coordinate"),
            dragmode = "lasso"
          )
      })

      observeEvent(event_data("plotly_selected"), {
        sel_data <- event_data("plotly_selected")
        if (!is.null(sel_data)) {
          sel_points <- sel_data$pointNumber + 1
          sel_coords <- coords[sel_points, ]
          if (nrow(sel_coords) >= 3) {
            sel_coords <- sel_coords[chull(sel_coords[, c("x", "y")]), ]
            poly_coords <- as.matrix(rbind(sel_coords[, c("x", "y")], sel_coords[1, c("x", "y")]))
            poly <- st_polygon(list(poly_coords))
            poly_sf <- st_sfc(poly)
            poly_sf <- st_make_valid(poly_sf)

            points_sf <- st_as_sf(coords, coords = c("x", "y"))
            inside <- st_within(points_sf, poly_sf, sparse = FALSE)[, 1]

            new_mask <- roi_mask()
            new_mask[inside] <- 1
            roi_mask(new_mask)
            output$status <- renderText("Points selected for ROI.")
          } else {
            showNotification("Select at least 3 points", type = "warning")
          }
        }
      })

      observeEvent(input$save_roi, {
        name <- input$roi_name
        if (name == "") {
          showNotification("Please enter a name for the ROI.", type = "error")
          return()
        }

        seurat <- rv()
        seurat[[name]] <- roi_mask()
        rv(seurat)
        roi_mask(rep(0, nrow(coords)))
        output$status <- renderText(paste0("Saved ROI to metadata column: ", name))
      })

      observeEvent(input$reset_roi, {
        roi_mask(rep(0, nrow(coords)))
        output$status <- renderText("ROI selection reset.")
      })

      observeEvent(input$done, {
        stopApp(rv())
      })
    }
  )))
}
