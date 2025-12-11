################################################################################
## PURPOSE: Interactive Shiny App to manually select "Endemic" LGAs by clicking
##          on a map, mimicking NCDC Situation Reports.
## OUTPUT:  Saves a CSV of selected LGA Codes/Names.
################################################################################

library(shiny)
library(leaflet)
library(sf)
library(terra)
library(dplyr)
library(here)
library(leaflet.extras)
library(DT) # For the table

# --- 1. Data Setup -----------------------------------------------------------
# Load Nigeria Admin 2 Shapes
gadm_dir <- here("data", "raw", "gadm")
nga_files <- list.files(gadm_dir, pattern = "NGA.*_2_.*\\.rds", full.names = TRUE)

if (length(nga_files) == 0) stop("Nigeria GADM file not found.")

# Prepare LGA Data (Admin 2)
nga_lga <- readRDS(nga_files[1]) |> 
  st_as_sf() |> 
  st_transform(crs = 4326) |> 
  select(GID_2, NAME_1, NAME_2, geometry) |> 
  mutate(
    Label = paste0(NAME_2, " (", NAME_1, ")"),
    layerId = GID_2 
  )

# Prepare State Data (Admin 1) for Visual Overlay
nga_state <- nga_lga |> 
  group_by(NAME_1) |> 
  summarise(geometry = st_union(geometry), .groups = "drop")

print("Data loaded. Launching App...")


# --- 2. UI Definition --------------------------------------------------------
ui <- fluidPage(
  titlePanel("NCDC Validation Tool"),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("Workflow"),
      tags$ol(
        tags$li("Enter the Report Year."),
        tags$li("Find red areas on the NCDC PDF."),
        tags$li("Click corresponding LGAs on the map to turn them RED."),
        tags$li("Click 'Download CSV' when finished."),
        tags$li("Click 'Clear Selection' to start the next year.")
      ),
      hr(),
      
      # Controls
      textInput("year_input", "Report Year:", value = "2023"),
      
      div(style="display:inline-block", 
          downloadButton("downloadData", "Download CSV", class = "btn-primary")),
      div(style="display:inline-block", 
          actionButton("clear_all", "Clear Selection", icon = icon("trash"), class = "btn-danger")),
      
      hr(),
      
      # Live Table
      h4("Selected LGAs"),
      DTOutput("selected_table")
    ),
    
    mainPanel(
      width = 8,
      leafletOutput("map", height = "90vh")
    )
  )
)


# --- 3. Server Logic ---------------------------------------------------------
server <- function(input, output, session) {
  
  # Reactive storage for selected GIDs
  selected_ids <- reactiveVal(character(0))
  
  # A. Render Base Map
  output$map <- renderLeaflet({
    leaflet() |> 
      addProviderTiles(providers$CartoDB.Positron) |> 
      
      # Create a custom pane for States to ensure they are ALWAYS on top
      addMapPane("state_pane", zIndex = 450) |> 
      addMapPane("lga_pane", zIndex = 400) |>
      
      # Layer 1: LGAs (Clickable) - Thinner Borders
      addPolygons(
        data = nga_lga,
        layerId = ~layerId,
        fillColor = "white",
        fillOpacity = 0.1,
        color = "black",       
        weight = 0.5,          # Thinner borders (was 1.5)
        opacity = 1,
        label = ~Label,
        group = "LGAs",
        options = pathOptions(pane = "lga_pane"), # Send to lower pane
        highlightOptions = highlightOptions(
          color = "red", weight = 2, bringToFront = TRUE
        )
      ) |> 
      
      # Layer 2: States (Visual Overlay) - On Top
      addPolygons(
        data = nga_state,
        fill = FALSE,
        color = "purple",      
        weight = 3,            
        opacity = 1,
        options = pathOptions(clickable = FALSE, pane = "state_pane") # Force to top pane
      ) |> 
      
      # Search Bar
      addSearchFeatures(
        targetGroups = "LGAs", 
        options = searchFeaturesOptions(
          zoom = 10, openPopup = TRUE, firstTipSubmit = TRUE,
          textPlaceholder = "Search LGA..."
        )
      )
  })
  
  # B. Handle Clicks (Toggle Logic)
  observeEvent(input$map_shape_click, {
    click <- input$map_shape_click
    
    # Ignore clicks on things without IDs (background)
    if (is.null(click$id)) return()
    
    current_set <- selected_ids()
    clicked_id <- click$id
    
    # Get the geometry of the clicked polygon to redraw it
    # (This is faster than redrawing the whole map)
    poly_to_update <- nga_lga[nga_lga$GID_2 == clicked_id, ]
    
    if (clicked_id %in% current_set) {
      # -- DESELECT --
      new_set <- setdiff(current_set, clicked_id)
      
      # Redraw as WHITE
      leafletProxy("map") |> 
        addPolygons(
          data = poly_to_update,
          layerId = clicked_id,
          fillColor = "white", fillOpacity = 0.1,
          color = "black", weight = 1.5, opacity = 1,
          label = poly_to_update$Label,
          highlightOptions = highlightOptions(color = "red", weight = 3)
        )
      
    } else {
      # -- SELECT --
      new_set <- c(current_set, clicked_id)
      
      # Redraw as RED
      leafletProxy("map") |> 
        addPolygons(
          data = poly_to_update,
          layerId = clicked_id,
          fillColor = "red", fillOpacity = 0.6, # Bright Red
          color = "black", weight = 1.5, opacity = 1,
          label = poly_to_update$Label,
          highlightOptions = highlightOptions(color = "darkred", weight = 3)
        )
    }
    
    # Update the reactive variable
    selected_ids(new_set)
  })
  
  
  # C. Clear All Button
  observeEvent(input$clear_all, {
    # 1. Reset reactive
    selected_ids(character(0))
    
    # 2. Reset Map Visuals (Redraw the base white layer)
    # This is the fastest way to clear many red polygons at once
    leafletProxy("map") |> 
      clearGroup("LGAs") |> # Clear existing
      addPolygons(          # Re-add white
        data = nga_lga,
        layerId = ~layerId,
        fillColor = "white", fillOpacity = 0.1,
        color = "black", weight = 1.5, opacity = 1,
        label = ~Label,
        group = "LGAs",
        highlightOptions = highlightOptions(color = "red", weight = 3, bringToFront = FALSE)
      )
  })
  
  
  # D. Table Output
  output$selected_table <- renderDT({
    req(length(selected_ids()) > 0)
    
    nga_lga |> 
      st_drop_geometry() |> 
      filter(GID_2 %in% selected_ids()) |> 
      select(State = NAME_1, LGA = NAME_2, Code = GID_2) |> 
      arrange(State, LGA)
    
  }, options = list(pageLength = 10, dom = 'tp', ordering=FALSE))
  
  
  # E. Download Handler
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("ncdc_affected_lgas_", input$year_input, ".csv")
    },
    content = function(file) {
      out_df <- nga_lga |> 
        st_drop_geometry() |> 
        filter(GID_2 %in% selected_ids()) |> 
        select(GID_2, NAME_1, NAME_2) |> 
        rename(State = NAME_1, LGA = NAME_2, LGA_Code = GID_2) |> 
        mutate(Year = input$year_input, Status = 1)
      
      write.csv(out_df, file, row.names = FALSE)
    }
  )
}

# Run
shinyApp(ui, server)