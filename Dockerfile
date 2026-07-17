FROM bioconductor/bioconductor_docker:latest
MAINTAINER Disease Transcriptomics Lab <diseasetranscriptomicslab@gmail.com>

RUN apt-get update && apt-get -y upgrade && apt-get -y autoremove && rm -rf /var/lib/apt/lists/*

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'https://cloud.r-project.org'; options(repos = r);" > ~/.Rprofile

# CRAN packages used by app.R / helper_functions.R
RUN Rscript -e "install.packages(c( \
      'shiny', 'shinythemes', 'fontawesome', 'DT', 'plotly', 'ggplot2', \
      'ggrepel', 'shinycssloaders', 'dplyr', 'tidyr', 'ggpubr', 'png', \
      'base64enc', 'cowplot', 'msigdbr', 'qs2' \
    ))"

# Bioconductor packages used by app.R
RUN Rscript -e "BiocManager::install(c('limma', 'fgsea'), update = FALSE, ask = FALSE)"

# Copy app source code (the data/ directory is excluded via .dockerignore;
# it is supplied at runtime through a mounted volume, see below)
WORKDIR /home/app
COPY . .

EXPOSE 3838

# Launch the CHARM Shiny app when the container starts, listening on all
# interfaces so it's reachable from outside the container. Shiny sets its
# working directory to the app directory, so app.R's relative
# "data/QS2_Files" / "data/QS_Files" lookups resolve against
# /home/app/data -- no path changes are needed as long as the data volume
# is mounted there:
#
#   docker run -p 3838:3838 -v /path/to/data:/home/app/data <image>
CMD ["R", "-e", "shiny::runApp(appDir = '/home/app', host = '0.0.0.0', port = 3838, launch.browser = FALSE)"]
