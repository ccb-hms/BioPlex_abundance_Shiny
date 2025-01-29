# Load required libraries for Shiny
library(pins)
library(shiny)
library(scater)

# Read the pinned data from Posit Connect
board <- board_connect(
    server = Sys.getenv("POSIT_SERVER_URL"), 
    key = Sys.getenv("POSIT_API_KEY")
)
sce <- pin_read(board, name = "anthony-alexander_christidis@hms.harvard.edu/marrow_single_cell_data")

# UI
ui <- fluidPage(
    titlePanel("Single Cell RNA-seq Analysis"),
    sidebarLayout(
        sidebarPanel(
            radioButtons("dimred", "Choose Dimensionality Reduction:",
                         choices = c("UMAP", "t-SNE"), selected = "UMAP"),
            width = 3  # Reduce sidebar width
        ),
        mainPanel(
            # Wrap the plot in a column to control its width
            column(
                width = 8,  # Adjust this value to change plot width
                plotOutput("scatterPlot")
            )
        )
    )
)

# Server
server <- function(input, output, session) {
    
    output$scatterPlot <- renderPlot({
        if (input$dimred == "UMAP") {
            plotReducedDim(sce, dimred = "UMAP", color_by = "reclustered.broad")
        } else {
            plotReducedDim(sce, dimred = "TSNE", color_by = "reclustered.broad")
        }
    }, height = 500)  # Set a fixed height for the plot
}

# Run the app
shinyApp(ui = ui, server = server)