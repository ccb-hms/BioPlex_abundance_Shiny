#--------------
# PIN DATASETS
#--------------

# library(SummarizedExperiment)
#library(dplyr)

#se <- readRDS("~/Desktop/R/docker-data/BioPlex_data/data/bioplex_data_08_2025/ptm_abundance_data/PTM_Summary_se_10Sep25.rds")

# Does now include ubiquintination 08/2025
#colData(se) 
#ptm_new <- as.data.frame(rowData(se)) |>  select(Symbol, PTM, Site)

# Join the assay to the row data
#ptm_new <- cbind(ptm_new, (assay(se)))

# Remove multiple sites - from
#ptm_new <- ptm_new |> 
#  filter(!grepl(";", Site))

# Rename "Symbol" to "Protein" columns
#names(ptm_new)[1] <- "Protein"

#saveRDS(ptm_new, "~/Desktop/R/docker-data/ptm_assay_2025.rds")


# Navigate to working directory 

# Load in rds
prot_df <- readRDS("~/Desktop/R/docker-data/prot_assay.rds")
ptm_df <- readRDS("~/Desktop/R/docker-data/ptm_assay_2025.rds")

# Read the server URL and API keys
readRenviron(".Renviron")
server_url <- Sys.getenv("POSIT_SERVER_URL")
api_key <- Sys.getenv("POSIT_API_KEY")

# Connect to Posit server using environment variables
board <- board_connect(
    server = server_url,
    key = api_key
)

# Pin the protein and ptm objects
pin_write(board, prot_df, name = "protein_abundance_data", type = "rds")
pin_write(board, ptm_df, name = "ptm_exp_2025_data", type = "rds")
