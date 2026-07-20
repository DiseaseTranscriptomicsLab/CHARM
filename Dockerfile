FROM bioconductor/bioconductor_docker:latest
MAINTAINER Disease Transcriptomics Lab <diseasetranscriptomicslab@gmail.com>

RUN apt-get update && apt-get -y upgrade && apt-get -y autoremove && rm -rf /var/lib/apt/lists/*

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'https://cloud.r-project.org'; options(repos = r);" > ~/.Rprofile

# CRAN packages used by app.R / helper_functions.R -- includes packages only
# referenced via pkg::fn() (Rtsne, reshape2, qs, purrr, scales, tibble),
# which install.packages() only guarantees for packages actually requested
# here, not arbitrary code elsewhere that happens to call into them.
#
# Interactive plots use highcharter, not plotly/ggplotly -- CHARM was
# migrated off plotly because its large JS bundle was the most likely
# casualty of the httpuv keep-alive bug patched below, and highcharter is
# also what the lab's other ShinyProxy-deployed apps (betAS, voyAGEr) use.
RUN Rscript -e "install.packages(c( \
      'shiny', 'shinythemes', 'fontawesome', 'DT', 'highcharter', 'ggplot2', \
      'ggrepel', 'shinycssloaders', 'dplyr', 'tidyr', 'ggpubr', 'png', \
      'base64enc', 'cowplot', 'msigdbr', 'qs2', 'qs', 'Rtsne', 'reshape2', \
      'purrr', 'scales', 'tibble' \
    ))"

# Bioconductor packages used by app.R
RUN Rscript -e "BiocManager::install(c('limma', 'fgsea'), update = FALSE, ask = FALSE)"

# Rebuild shiny with a raised HTTP header-size limit. httpuv's underlying
# http-parser has a bug that intermittently 503s requests when a keep-alive
# connection is reused for several sequential requests -- exactly what
# happens loading a page's static assets. This is what originally broke
# CHARM's plotly-based plots in production (client-side "Plotly is not
# defined", which aborted Shiny's output-update batch and took down
# whichever other outputs were queued alongside it) before the move to
# highcharter. Kept as general hardening regardless of widget library,
# per https://www.shinyproxy.io/documentation/troubleshooting/ -- also
# applied in voyAGEr's Dockerfile for the same reason.
RUN Rscript -e "install.packages(c('withr'), repos='https://cloud.r-project.org/')"
RUN Rscript -e "withr::with_makevars(c(PKG_CPPFLAGS='-DHTTP_MAX_HEADER_SIZE=0x7fffffff'), {install.packages(c('shiny'), repos='https://cloud.r-project.org/')}, assignment = '+=')"

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
