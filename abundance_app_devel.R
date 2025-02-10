library(shiny)
library(bslib)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(DT)

## Feb 10, 2025
# ---------------------------------------
# LOAD DATA
# ---------------------------------------
df <- readRDS("~/Desktop/R/docker-data/prot_assay.rds")
ptm_df <- readRDS("~/Desktop/R/docker-data/ptm_assay.rds")

# ---------------------------------------
# SIDEBAR PANEL
# ---------------------------------------

ui <- page_sidebar(
    title = "Differential protein abundance and PTM explorer",
    
    # SIDEBAR
    sidebar = sidebar(
        selectInput("protein_choice", 
                    "Select Protein:",
                    choices = unique(df$Protein)),
         selectInput("PTM_type",
                     "Select PTM type:",
                     choices = c("Acetylation" = "acetylation",
                                 "Phosphorylation" = "phosphorylation")),
        selectInput("analysis_type",
                     "Select analysis:",
                     choices = c("ANOVA" = "anova",
                                 "Kruskal-Wallis Rank Sum" = "kruskal.test")),
        selectInput("ref_group", "Select Reference Cell Line:", 
                    choices = c("HEK_293T", "HCT116", "U2OS", "RPE1", "HeLa")),
    ),
    
    div(
    style = "height: calc(100vh - 60px);",  
    # Makes main content scrollable
    layout_columns(
        col_widths = c(6, 6),
        
        # PROTEIN BOXPLOT
        card(
            card_header(paste("Protein abundance boxplot")),
            plotOutput("value_boxplot", height = "300px"),
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
        DT::dataTableOutput("ptm_table"))
    )
)

# ---------------------------------------
# ANALYSIS AND RENDER PLOTTING
# ---------------------------------------

server <- function(input, output) {
    
    # ---------------------------------------
    # PROTEIN TABLE
    # ---------------------------------------
    
    output$protein_table <- renderTable({
        # Filter the dataframe to show only the selected protein
        df[df$Protein == input$protein_choice, ]
    })
    
    
    # Wrangle and boxplot
    output$value_boxplot <- renderPlot({
        
        # Filter for selected protein and reshape data for plotting
        df_filtered <- df[df$Protein == input$protein_choice, ]
        
        df_filtered <- df_filtered |>
            pivot_longer(
                cols = -Protein,
                names_to = c("cell_line", "replicate"),
                names_pattern = "(.*)_([0-9])$",
                values_to = "relative_expression"
            ) |>
            mutate(cell_line = gsub("_$", "", cell_line)) ##|> 
            ##mutate(relative expression = log2(relative expression))


        ggboxplot(df_filtered, x = "cell_line", y = "relative_expression",
                  color = "cell_line", palette = "jco",
                  add = "jitter") +
            theme(legend.position = "")
        
        # Change stats based on analysis type
        if (input$analysis_type == "anova") {
            # Perform ANOVA
            ggboxplot(df_filtered, x = "cell_line", y = "relative_expression",
                      color = "cell_line", palette = "jco",
                      add = "jitter") +
                theme(legend.position = "") +
                stat_compare_means(method = "anova",
                                   vjust = 25) +
                stat_compare_means(label = "p.signif", method = "t.test",
                         ref.group = input$ref_group) 
        } else if (input$analysis_type == "kruskal.test") {
            # Perform kruskal.test t-test
            ggboxplot(df_filtered, x = "cell_line", y = "relative_expression",
                      color = "cell_line", palette = "jco",
                      add = "jitter") +
                theme(legend.position = "") +
                stat_compare_means(method = "kruskal.test",
                                   vjust = 25) +
                stat_compare_means(label = "p.signif", method = "t.test",
                         ref.group = input$ref_group)
        }
    })
    
    # Pairwise t-test table
    output$stats <- renderTable({
    
    df_filtered <- df[df$Protein == input$protein_choice, ]
        
    df_filtered <- df_filtered |>
            pivot_longer(
                cols = -Protein,
                names_to = c("cell_line", "replicate"),
                names_pattern = "(.*)_([0-9])$",
                values_to = "relative_expression"
            ) |>
            mutate(cell_line = gsub("_$", "", cell_line))
    
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
   
    }, digits = 3)
    
    # ---------------------------------------
    # PTM TABLE AND PLOT
    # ---------------------------------------
    
    selected_row <- reactiveVal()
    
    # Modify the table output to store the selected row
    output$ptm_table <- renderDataTable({
        cdf <- ptm_df[ptm_df$Protein == input$protein_choice, ]
        
        if (input$PTM_type == "acetylation") {
            cdf <- cdf[cdf$PTM == "Acetylation", ]
        } else if (input$PTM_type == "phosphorylation") {
            cdf <- cdf[cdf$PTM == "Phosphorylation", ]
        }
        
        DT::datatable(cdf, 
                  selection = 'single',
                  options = list(pageLength = 5))
    })
    
    # Update the selected row when user clicks
    observeEvent(input$ptm_table_rows_selected, {
        cdf <- ptm_df[ptm_df$Protein == input$protein_choice, ]
        
        if (input$PTM_type == "acetylation") {
            cdf <- cdf[cdf$PTM == "Acetylation", ]
        } else if (input$PTM_type == "phosphorylation") {
            cdf <- cdf[cdf$PTM == "Phosphorylation", ]
        }
        
        # Store the selected row data
        selected_row(cdf[input$ptm_table_rows_selected, ])
    })
    
    # Modify the boxplot to use the selected row data
    output$ptm_boxplot <- renderPlot({
        req(selected_row())  # Ensure we have a selected row
        
        # Get the selected row data
        selected_data <- selected_row()
        
        # Create long format data from the selected row
        ptm_long <- selected_data |> 
            pivot_longer(
                cols = -c(Protein, Site, PTM),
                names_to = c("cell_line", "replicate"),
                names_pattern = "(.*)_([0-9])$",
                values_to = "average relative abundance"
              ) |>
              mutate(cell_line = gsub("_$", "", cell_line))
        
        # Create the boxplot
        ggboxplot(ptm_long, 
                  x = "cell_line", 
                  y = "average relative abundance",
                  color = "cell_line", 
                  palette = "jco",
                  add = "jitter") +
            theme(legend.position = "") +
            ggtitle(paste("Site:", selected_data$Site))
    })
}

shinyApp(ui, server)