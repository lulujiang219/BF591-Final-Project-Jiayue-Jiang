
# Interactive Transcriptomics Web App

## Overview
Welcome to the Interactive Transcriptomics Web App. This project was developed as part of the BF591 course, with the primary goal of enhancing the understanding of transcriptomic changes across various biological conditions. The application was built using R Shiny to provide an interactive web interface that allows for comprehensive genomic data analysis.

## Application Layout
The application features a multi-tab layout that guides the user through several stages of transcriptomic analysis:

1. **Sample Information Exploration**: Load and examine the sample information matrix.
2. **Counts Matrix Exploration**: Visualize gene expression counts and apply filtering thresholds.
3. **Differential Expression Analysis**: Explore differential expression results.
4. **Gene Set Enrichment Analysis**: Perform and visualize gene set enrichment analysis.
5. **Individual Gene Expression Visualization**: Plot gene expression levels for selected genes and categories.

## Installation and Usage
To use the application, follow these steps:

1. Clone the repository to your local machine.
2. Ensure R and the required packages listed in `requirements.txt` are installed.
3. Run the app by opening the `app.R` file in RStudio and clicking 'Run App'.

## Data source: GEO Dataset GSE64810 
mRNA-Seq Expression profiling of human post-mortem BA9 brain tissue for Huntington's Disease and neurologically normal individuals

## Data Preparation
The application accepts CSV formatted data for analysis. Users are expected to process raw data to format it appropriately for the Shiny app. Example datasets and their formats are provided for reference.

## Code Organization
The repository is organized as follows:

- `app.R`: The application file containing server logic.
- `data/`: Directory containing example datasets and processed data files.
- `ui.R`: The application file containing UI logic.
- `main.R`: The main application file to run both R files

