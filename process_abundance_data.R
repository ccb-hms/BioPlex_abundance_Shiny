#--------------
# PIN DATASETS
#--------------

# Load in rds
prot_df <- readRDS("~/Desktop/R/docker-data/prot_assay.rds")
ptm_df <- readRDS("~/Desktop/R/docker-data/ptm_assay.rds")

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
pin_write(board, ptm_df, name = "ptm_exp_data", type = "rds")