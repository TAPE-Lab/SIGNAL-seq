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
reads_df =  data.frame(sl_name = c('ex0015_rna_76', 'ex0015_rna_77', 'ex0015_rna_78', 'ex0015_rna_79'), 
                       reads = c(860141497, 770196821, 800992534, 913336602))

# Run the splitRtool pipeline
# Each sublibrary is contained within its own folder in the data_folder folder and must contain zUMIs output, named by sublib name.
run_split_pipe(mode = 'single', # Merge sublibraries or process seperately
               n_sublibs = 4, # How many to sublibraries are present
               data_folder = "./data_folder/", # Location of zUMIs data directory
               output_folder = "./splitRtools_outputs", # Output folder path
               filtering_mode = "manual", # Filter by knee (standard) or manual value (default 1000) transcripts
               filter_value = 500,
               count_reads = FALSE,
               total_reads = reads_df,
               gene_names = FALSE,
               fastq_path = NA, # Path to folder containing subibrary raw FastQ
               rt_bc = "./barcode_maps/barcodes_v2_48.csv", # RT barcode map
               lig_bc = "./barcode_maps/barcodes_v1.csv", # Ligation barcode map
               sample_map = "./barcode_maps/cell_metadata.xlsx" # RT plate layout file
)
