# base R shiny image
FROM rocker/shiny:latest

# install system dependencies
RUN apt-get update && apt-get install -y \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libglpk-dev \
    libwebp-dev \
    libssl-dev \
    libgit2-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# make a directory in the container
RUN mkdir /home/single-cell-shiny-app
RUN mkdir /home/single-cell-shiny-app/data
RUN mkdir /home/single-cell-shiny-app/www

# install R dependencies
RUN R -e "install.packages(c('Seurat', 'bslib', 'DT', 'rlang', 'shinythemes', 'shinyFeedback', 'rsconnect', 'qs', 'shinyBS', 'ggplot2', 'devtools', 'ragg'))"

# install presto from GitHub
RUN R -e "devtools::install_github('immunogenomics/presto')"

# copy the shiny app code
COPY app.R /home/single-cell-shiny-app/app.R
COPY helpers.R /home/single-cell-shiny-app/helpers.R
COPY www/logo.png /home/single-cell-shiny-app/www/logo.png

# expose the application port
EXPOSE 8180

# run the R shiny app
CMD ["R", "-e", "shiny::runApp('/home/single-cell-shiny-app', host = '0.0.0.0', port = 3838)"]

