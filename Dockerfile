FROM rocker/r-ver:latest
MAINTAINER Disease Transcriptomics Lab <diseasetranscriptomicslab@gmail.com>
# rocker/r-ver is the same lineage as bioconductor/bioconductor_docker
# (which this Dockerfile used before) but WITHOUT RStudio Server bolted on
# top - this app is a headless Shiny app and never uses it, so skipping it
# saves close to 1.7GB. The trade-off: bioconductor_docker also came with
# BiocManager and various system libraries pre-installed for building
# Bioconductor/CRAN packages from source; both are added explicitly below
# since rocker/r-ver doesn't include them.
RUN apt-get update && apt-get -y upgrade && apt-get install -y --no-install-recommends \
      build-essential \
      gfortran \
      cmake \
      pkg-config \
      libcurl4-openssl-dev \
      libssl-dev \
      libxml2-dev \
      libgit2-dev \
      libuv1-dev \
      zlib1g-dev \
      libpng-dev \
      liblapack-dev \
      libblas-dev \
      libglpk-dev \
    && apt-get -y autoremove && rm -rf /var/lib/apt/lists/*
# cmake: nloptr (a transitive dependency of ggpubr, via rstatix -> car ->
# pbkrtest -> lme4 -> nloptr) has required cmake to build its bundled
# NLopt C library from source since nloptr 2.0. bioconductor_docker ships
# a much broader build toolchain (cmake included) for compiling
# genomics/Bioconductor packages from source; rocker/r-ver doesn't, so
# this needs to be added back explicitly, same reasoning as the other
# system libraries above.
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'https://cloud.r-project.org'; options(repos = r);" > ~/.Rprofile
# install.packages() does NOT raise an error or give a non-zero exit code
# just because a package failed to install - it only prints a warning, so
# a transient mirror/network hiccup here would previously be cached by
# Docker as a "successful" layer even though BiocManager was never
# actually installed, surfacing as a confusing failure several steps
# later (at BiocManager::install()) instead of right here. This
# explicitly checks that it's actually present afterwards and fails the
# build immediately and clearly if not.
RUN Rscript -e '\
  install.packages("BiocManager"); \
  if (!"BiocManager" %in% rownames(installed.packages())) { \
    stop("Failed to install: BiocManager") \
  }'
# CRAN packages used by app.R / helper_functions.R -- includes packages only
# referenced via pkg::fn() (Rtsne, reshape2, purrr, scales, tibble), which
# install.packages() only guarantees for packages actually requested here,
# not arbitrary code elsewhere that happens to call into them. `qs` (the
# legacy serializer) is deliberately NOT installed anywhere in this
# Dockerfile: app.R's DATA_DIR / DATA_EXT logic picks data/QS2_Files (and
# therefore qs2::qs_read()) for the whole session whenever that directory
# exists, and only falls back to qs::qread() on data/QS_Files if it
# doesn't -- our deployed data volume only has data/QS2_Files, so the old
# `qs` package is never actually called and isn't worth installing (it's
# also been archived on CRAN, since superseded by qs2, which would
# otherwise need a CRAN-Archive workaround for no benefit).
#
# Interactive plots use highcharter, not plotly/ggplotly -- CHARM was
# migrated off plotly because its large JS bundle was the most likely
# casualty of the httpuv keep-alive bug patched below, and highcharter is
# also what the lab's other ShinyProxy-deployed apps (betAS, voyAGEr) use.
#
# Same install.packages()-doesn't-error caveat as the BiocManager step
# above applies here, and it matters more now than it did under
# bioconductor_docker: that image had most of these already compiled and
# pre-installed, so a from-source build failing here (e.g. a system -dev
# library rocker/r-ver doesn't ship that bioconductor_docker did) would
# previously go unnoticed at build time and only surface later as
# "container did not respond in time" -- library(x) erroring the instant
# app.R starts, which kills the R process before Shiny ever binds the
# port, which ShinyProxy can only see as a timeout. Checking here instead
# fails the build immediately and names the actual missing package.
RUN Rscript -e '\
  cran_pkgs <- c( \
      "shiny", "shinythemes", "fontawesome", "DT", "highcharter", "ggplot2", \
      "ggrepel", "shinycssloaders", "dplyr", "tidyr", "ggpubr", "png", \
      "base64enc", "cowplot", "msigdbr", "qs2", "Rtsne", "reshape2", \
      "purrr", "scales", "tibble" \
    ); \
  install.packages(cran_pkgs); \
  missing_pkgs <- setdiff(cran_pkgs, rownames(installed.packages())); \
  if (length(missing_pkgs) > 0) { \
    stop("Failed to install: ", paste(missing_pkgs, collapse = ", ")) \
  }'
# Bioconductor packages used by app.R
RUN Rscript -e '\
  bioc_pkgs <- c("limma", "fgsea"); \
  BiocManager::install(bioc_pkgs, update = FALSE, ask = FALSE); \
  missing_pkgs <- setdiff(bioc_pkgs, rownames(installed.packages())); \
  if (length(missing_pkgs) > 0) { \
    stop("Failed to install: ", paste(missing_pkgs, collapse = ", ")) \
  }'
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
