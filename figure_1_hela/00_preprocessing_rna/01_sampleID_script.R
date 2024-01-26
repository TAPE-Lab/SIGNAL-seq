# Load devtools
library(devtools)

# download required packages from bioconductor if needed for first time
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("zellkonverter", "scater", "ShortRead", "DropletUtils"))

# download and install splitRtools from github
devtools::install_github("https://github.com/TAPE-Lab/splitRtools")

# Load splitRtools
library(splitRtools)

# specify raw reads per sublibrary
reads_df =  data.frame(sl_name = c('ex0003_hela'), reads = c(183030567))

# Run the splitRtool pipeline
# Each sublibrary is contained within its own folder in the data_folder folder and must contain zUMIs output, named by sublib name.
run_split_pipe(mode = 'single', # Merge sublibraries or process seperately
               n_sublibs = 1, # How many to sublibraries are present
               data_folder = "~/PATH/sulib_data/", # Location of zUMIs data directory
               output_folder = "~/PATH/pipe_outputs", # Output folder path
               filtering_mode = "manual", # Filter by knee (standard) or manual value (default 1000) transcripts
               filter_value = 500,
               count_reads = FALSE,
               total_reads = reads_df,
               fastq_path = NA, # Path to folder containing subibrary raw FastQ
               rt_bc = "~/PATH/barcode_maps/barcodes_v2_48.csv", # RT barcode map
               lig_bc = "~/PATH//barcode_maps/barcodes_v1.csv", # Ligation barcode map
               sample_map = "~/PATH/barcode_maps/cell_metadata.xlsx" # RT plate layout file
)
