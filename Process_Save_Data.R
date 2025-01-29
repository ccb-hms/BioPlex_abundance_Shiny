# Load necessary libraries
library(scRNAseq)
library(scuttle)
library(scran)
library(scater)

# Fetch the He Organ Atlas Data for the Marrow tissue
sce <- HeOrganAtlasData(tissue = c("Marrow"), ensembl = FALSE)

# Calculate mitochondrial protein-coding feature percentages if needed
is_mito <- grepl("^MT-", rowData(sce)$symbol)
sce <- addPerCellQC(sce, subsets = list(Mito = is_mito))

# Run quality control by flagging low-quality cells based on mitochondrial content
# percent_subsets looks within the data subsets calculated by addPerCellQC
qc_metrics <- quickPerCellQC(colData(sce), 
                             percent_subsets = "subsets_Mito_percent")

# Filter cells using the QC metrics obtained
sce <- sce[, !qc_metrics$discard]

# Normalization using scran
sce <- logNormCounts(sce)

# Selecting highly variable genes
sce <- sce[getTopHVGs(sce, n = 1e3),]

# Compute PCA
sce <- runPCA(sce)

# Compute UMAP
sce <- runUMAP(sce, dimred = "PCA")

# Compute t-SNE
sce <- runTSNE(sce, dimred = "PCA")

# Read the server URL and API keys
readRenviron(".Renviron")
server_url <- Sys.getenv("POSIT_SERVER_URL")
api_key <- Sys.getenv("POSIT_API_KEY")

# Connect to Posit server using environment variables
board <- board_connect(
    server = server_url,
    key = api_key
)

# Pin the SingleCellExperiment object
pin_write(board, sce, name = "marrow_single_cell_data", type = "rds")
