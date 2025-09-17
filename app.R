library(shiny)
library(bslib)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(DT)
library(pins)

## Sep 12, 2025
#---------------------
# READ DATA FROM PINS
#---------------------

# Read the server URL and API keys
readRenviron(".Renviron")
server_url <- Sys.getenv("POSIT_SERVER_URL")
api_key <- Sys.getenv("POSIT_API_KEY")

# Read the pinned data from Posit Connect
board <- board_connect(
    server = Sys.getenv("POSIT_SERVER_URL"), 
    key = Sys.getenv("POSIT_API_KEY")
)

prot_df <- pin_read(board, 
                    name = "tram_nguyen@hms.harvard.edu/protein_abundance_data")
ptm_df <- pin_read(board, 
                   name = "tram_nguyen@hms.harvard.edu/ptm_exp_2025_data")

#---------------------
# FORMAT SIDEBAR PANEL
#---------------------

# Get all unique proteins from both dataframes for selection
all_proteins <- unique(c(prot_df$Protein, ptm_df$Protein))

ui <- page_sidebar(
    title = "Differential protein abundance and PTM explorer",
    
    ## Sidebar selection options
    sidebar = sidebar(
        selectizeInput("protein_choice", 
               "Select Protein:",
               choices = sort(all_proteins), 
               options = list(maxOptions = NULL)), 
         selectInput("PTM_type",
                     "Select PTM type:",
                     choices = c("Acetylation" = "acetylation",
                                 "Phosphorylation" = "phosphorylation",
                                 "Ubiquitination" = "ubiquitination")),
        selectInput("analysis_type",
                     "Select analysis:",
                     choices = c("ANOVA" = "anova",
                                 "Kruskal-Wallis Rank Sum" = "kruskal.test")),
        selectInput("ref_group", "Select Reference Cell Line:", 
                    choices = c("HEK_293T", "HCT116", "U2OS", "RPE1", "HeLa")),
    ),
    
    ## Format layout
    div(
    style = "height: calc(100vh - 60px);",  
    layout_columns(
        col_widths = c(6, 6),
        
        # PROTEIN BOXPLOT
        card(
            card_header(paste("Protein abundance boxplot")),
            plotOutput("protein_boxplot", height = "300px"),
            textOutput("stats_result")
        ),
        
        # PTM BOXPLOT
        card(
          card_header(paste("PTM boxplot")),
            plotOutput("ptm_boxplot", height = "300px"),
            textOutput("ptm_selection")
        )
    ),
    
    # PROTEIN TABLE
    card(
        card_header("Selected protein abundance"),
        tableOutput("protein_table")
    ),
    
    # PTM TABLE
    card(
        card_header("Selected ptm abundance"),
        DT::dataTableOutput("ptm_table"), rownames = FALSE)
    )
)

# -----------------------------
# ANALYSIS AND PLOT RENDERING
# -----------------------------

server <- function(input, output) {
    
    # Initialize reactive value for selected row - ADD THIS LINE!
    selected_row <- reactiveVal(NULL)
    
    
    # Reactive to check data availability
    data_availability <- reactive({
        protein <- input$protein_choice
        
        prot_exists <- protein %in% prot_df$Protein
        ptm_exists <- protein %in% ptm_df$Protein
        
        list(
            protein = protein,
            prot_exists = prot_exists,
            ptm_exists = ptm_exists
        )
    })
    
    
    
    #------------------------
    # SELECTED PROTEIN TABLE
    #------------------------
    output$protein_table <- renderTable({
        availability <- data_availability()
        
        if (availability$prot_exists) {
            # Filter the dataframe to show only the selected protein
            prot_df[prot_df$Protein == input$protein_choice, ]
        } else {
            # Return empty dataframe with message
            data.frame(Message = paste("No protein data available for", 
                                       input$protein_choice))
        }
    })
    
    #---------------------------
    # PROTEIN ABUNDANCE BOXPLOT
    #---------------------------
    output$protein_boxplot <- renderPlot({
        availability <- data_availability()
        
        if (availability$prot_exists) {
            # Filter for selected protein and reshape data for plotting
            df_filtered <- prot_df[prot_df$Protein == input$protein_choice, ]
            
            df_filtered <- df_filtered |>
                pivot_longer(
                    cols = -Protein,
                    names_to = c("cell_line", "replicate"),
                    names_pattern = "(.*)_([0-9])$",
                    values_to = "relative expression"
                ) |>
                mutate("cell line" = gsub("_$", "", cell_line))
            
            # Change stats based on analysis type
            if (input$analysis_type == "anova") {
                # Perform ANOVA
                ggboxplot(df_filtered, x = "cell line", y = "relative expression",
                          color = "cell line", palette = "jco",
                          add = "jitter") +
                    theme(legend.position = "") +
                    ggtitle(paste("Protein:", input$protein_choice)) +
                    stat_compare_means(method = "anova",
                                       vjust = 15) +
                    stat_compare_means(label = "p.signif", method = "t.test",
                             ref.group = input$ref_group) 
                
            } else if (input$analysis_type == "kruskal.test") {
                # Perform kruskal.test t-test
                ggboxplot(df_filtered, x = "cell line", y = "relative expression",
                          color = "cell line", palette = "jco",
                          add = "jitter") +
                    theme(legend.position = "") +
                    ggtitle(paste("Protein:", input$protein_choice)) +
                    stat_compare_means(method = "kruskal.test",
                                       vjust = 15) +
                    stat_compare_means(label = "p.signif", method = "t.test",
                             ref.group = input$ref_group)
            }
        } else {
            # Create empty plot with message
            ggplot() +
                annotate("text", x = 0.5, y = 0.5, 
                        label = paste("No protein abundance data available for", input$protein_choice),
                        size = 4, hjust = 0.5, vjust = 0.5) +
                xlim(0, 1) + ylim(0, 1) +
                theme_void() +
                theme(plot.background = element_rect(fill = "white", color = NA))
        }
    })
    
    #---------------------------
    # CALCULATE PAIRWISE T-TEST
    #---------------------------
    ## Not rendered currently
    output$stats <- renderTable({
        availability <- data_availability()
        
        if (availability$prot_exists) {
            df_filtered <- prot_df[prot_df$Protein == input$protein_choice, ]
                
            df_filtered <- df_filtered |>
                    pivot_longer(
                        cols = -Protein,
                        names_to = c("cell_line", "replicate"),
                        names_pattern = "(.*)_([0-9])$",
                        values_to = "relative expression"
                    ) |>
                    mutate("cell line" = gsub("_$", "", cell_line))
            
            # Convert 'Group' to a factor (if it's not already)
            df_filtered$cell_line <- factor(df_filtered$cell_line, 
                                            levels = unique(df_filtered$cell_line))
            
            # Set reference group
            df_filtered$cell_line <- relevel(df_filtered$cell_line, 
                                             ref = input$ref_group)
            
            # Now perform the pairwise t-test
            t_result <- pairwise.t.test(df_filtered$relative_expression, 
                           df_filtered$cell_line, p.adjust.method = "bonferroni")
            
            t_matrix <- t_result$p.value
            
            p_df <- as.data.frame(as.table(t_matrix)) |> 
                mutate(pval = as.numeric(Freq)) |> 
                rename(Group1 = Var1, Group2 = Var2, P_Value = pval) |> 
                select(Group1, Group2, P_Value)
            
            p_df$P_Value <- formatC(p_df$P_Value, digits = 3) # no rounding
            
            p_df
        } else {
            # Return empty results
            data.frame(Message = "No data available for statistical analysis")
        }
    }, digits = 3)
    
   #------------
    # PTM TABLE 
    #------------
    
    # Modify the table output to store the selected row
    output$ptm_table <- renderDataTable({
        availability <- data_availability()
        
        if (availability$ptm_exists) {
            # Use filter() to avoid NA rows
            cdf <- ptm_df %>% 
                filter(Protein == input$protein_choice) %>%
                mutate(across(where(is.numeric), ~format(round(., 2), nsmall = 2)))
                   
            if (input$PTM_type == "acetylation") {
                cdf <- cdf %>% filter(PTM == "Acetylation")
            } else if (input$PTM_type == "phosphorylation") {
                cdf <- cdf %>% filter(PTM == "Phosphorylation")
            } else if (input$PTM_type == "ubiquitination") {
                cdf <- cdf %>% filter(PTM == "Ubiquitination")
            }
            
            if (nrow(cdf) > 0) {
                DT::datatable(cdf, 
                          selection = 'single',
                          options = list(pageLength = 5), rownames = FALSE)
            } else {
                empty_df <- data.frame(Message = paste("No", input$PTM_type, "data available for", input$protein_choice))
                DT::datatable(empty_df, 
                          options = list(pageLength = 5), rownames = FALSE)
            }
        } else {
            empty_df <- data.frame(Message = paste("No PTM data available for", input$protein_choice))
            DT::datatable(empty_df, 
                      options = list(pageLength = 5), rownames = FALSE)
        }
    })
    
    # Update the selected row when user clicks
    observeEvent(input$ptm_table_rows_selected, {
        availability <- data_availability()

        # Only proceed if PTM data exists for this protein
        if (availability$ptm_exists) {
            # First get the PTM data for the selected protein (using filter to match table code)
            cdf <- ptm_df %>% 
                filter(Protein == input$protein_choice)
            
            # Then filter by PTM type (using filter to match table code)
            if (input$PTM_type == "acetylation") {
                cdf <- cdf %>% filter(PTM == "Acetylation")
            } else if (input$PTM_type == "phosphorylation") {
                cdf <- cdf %>% filter(PTM == "Phosphorylation")
            } else if (input$PTM_type == "ubiquitination") {
                cdf <- cdf %>% filter(PTM == "Ubiquitination")
            }
            
            # Store the selected row data (only if row was actually selected)
            if (length(input$ptm_table_rows_selected) > 0) {
                selected_row(cdf[input$ptm_table_rows_selected, ])
            }
        } else {
            # Clear selected row if no PTM data exists
            selected_row(NULL)
        }
    })
    
 
    
    #-----------------------
    # REACTIVE PTM BOXPLOT 
    #-----------------------
    
    # Modify the boxplot to use the selected row data
    output$ptm_boxplot <- renderPlot({
    availability <- data_availability()
    
    # Check if PTM data exists and if a row is selected
    if (availability$ptm_exists && !is.null(selected_row())) {
        # Get the selected row data
        selected_data <- selected_row()
        
        # Create long format data from the selected row
        ptm_long <- selected_data |> 
            pivot_longer(
                cols = -c(Protein, Site, PTM),
                names_to = c("cell_line", "replicate"),
                names_pattern = "(.*)_([0-9])$",
                values_to = "relative abundance"
              ) |>
              mutate("cell line" = gsub("_$", "", cell_line))
        
        # Change stats based on analysis type
        if (input$analysis_type == "anova") {
            # Perform ANOVA
            ggboxplot(ptm_long, 
                  x = "cell line", 
                  y = "relative abundance",
                  color = "cell line", 
                  palette = "jco",
                  add = "jitter") +
            theme(legend.position = "") +
            ggtitle(paste("Site:", selected_data$Site)) +
                stat_compare_means(method = "anova",
                                   vjust = 15) +
                stat_compare_means(label = "p.signif", method = "t.test",
                         ref.group = input$ref_group) 
        } else if (input$analysis_type == "kruskal.test") {
            # Perform kruskal.test
            ggboxplot(ptm_long, 
                  x = "cell line", 
                  y = "relative abundance",
                  color = "cell line", 
                  palette = "jco",
                  add = "jitter") +
            theme(legend.position = "") +
            ggtitle(paste("Site:", selected_data$Site)) +
                stat_compare_means(method = "kruskal.test",
                                   vjust = 15) +
                stat_compare_means(label = "p.signif", method = "t.test",
                         ref.group = input$ref_group)
        }
    } else {
        # Create empty plot with appropriate message
        message <- if (!availability$ptm_exists) {
            paste("No PTM data available for", input$protein_choice)
        } else {
            paste("Please select a PTM site from the table below to view the boxplot")
        }
        
        ggplot() +
            annotate("text", x = 0.5, y = 0.5, 
                    label = message,
                    size = 4, hjust = 0.5, vjust = 0.5) +
            xlim(0, 1) + ylim(0, 1) +
            theme_void() +
            theme(plot.background = element_rect(fill = "white", color = NA))
    }
})

}
shinyApp(ui, server)