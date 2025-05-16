# 🧬 Single Cell RNA-seq Shiny App

This **RShiny** application was developed by the **Goessling Lab** to support members of the **zebrafish research community** in the analysis and visualization of single-cell RNA sequencing (scRNA-seq) data.

The tool provides an intuitive interface that enables users to:

- 📁 Upload and explore their own scRNA-seq datasets  
- 📊 Generate UMAP visualizations for genes of interest  
- 🧪 Perform differential expression analyses between cell types  
- 🔍 Investigate gene expression patterns within specific cell populations  

Designed to facilitate discovery and data-driven insights, this application empowers researchers to interactively explore cellular and transcriptional dynamics in **zebrafish**.

## Requirements

This application requires the following:

- [Docker](https://docs.docker.com/get-docker/) must be installed on your system.

## Installation

This application can be easily installed and run using **Docker**, ensuring a consistent and dependency-free setup across systems.

1. **Clone the repository**

   Open a terminal and run:

   ```bash
   git clone https://github.com/sop51/RShinySingleCell.git
   cd RShinySingleCell/App-1

3. **Build the Docker image**
   
   Run the following command in the project directory (App-1):

   💡 Note for Mac users: --platform flag will have to be specified in this command (linux/amd64 or linux/x86_64)

   ```bash
   docker build -t scRNAseq-shiny-app .

5. **Run the application**

   ```bash
   docker run -p 3838:3838 scRNAseq-shiny-app
  
