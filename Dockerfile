# Dockerfile for Reproducible R Environment
# For maximum reproducibility across different systems
# R version 4.5.1 (2025-06-13)

FROM rocker/r-ver:4.5.1

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
# Based on requirements from R/Utilities/Helpers/check_system_dependencies.R
RUN apt-get update && apt-get install -y \
    # Core build tools
    build-essential \
    gfortran \
    cmake \
    # Required system tools (from check_system_dependencies.R)
    ghostscript \
    pandoc \
    imagemagick \
    # XML and networking
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libgit2-dev \
    # Graphics and fonts
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libcairo2-dev \
    libxt-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    # LaTeX/PDF generation dependencies
    wget \
    perl \
    # HDF5 for bioinformatics (MSnbase, xcms)
    libhdf5-dev \
    libnetcdf-dev \
    # Additional dependencies for Bioconductor packages
    libfftw3-dev \
    libgsl-dev \
    libgmp-dev \
    libglpk-dev \
    # GraphViz for network plots (needed for Rgraphviz, RBGL)
    graphviz \
    libgraphviz-dev \
    # Additional system dependencies identified by renv
    gdal-bin \
    git \
    libgdal-dev \
    libmagick++-dev \
    # Additional dependencies for spatial packages (raster, sp)
    libproj-dev \
    libgeos-dev \
    libudunits2-dev \
    && rm -rf /var/lib/apt/lists/*

# Set up renv for exact package restoration
RUN Rscript -e "install.packages('renv', repos='https://cloud.r-project.org')"

# Copy project files needed for renv
WORKDIR /analysis
COPY renv.lock renv.lock
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json

# Copy DESCRIPTION file (needed for renv explicit snapshot type)
COPY DESCRIPTION .

# Copy R utility functions (needed before restore for dependencies)
COPY R/ R/

# Copy analysis scripts
COPY All_Run/ All_Run/

# Copy database files
COPY Databases/ Databases/

# Create Outputs directory structure
RUN mkdir -p Outputs/Balloon_and_Volcano \
    Outputs/Enrichment \
    Outputs/Heatmaps \
    Outputs/Legacy \
    Outputs/LIMMA \
    Outputs/Mummichog\ Outputs \
    Outputs/PCA \
    Outputs/PLSDA \
    Outputs/Supplemental

# Restore R packages from renv.lock (this captures exact versions)
# This step takes ~30-45 minutes on first build but is cached
RUN Rscript -e "renv::restore(prompt = FALSE)"

# Install tinytex for PDF generation (if needed)
RUN Rscript -e "if (requireNamespace('tinytex', quietly = TRUE)) { tinytex::install_tinytex() }"
ENV PATH="${PATH}:/root/bin"

# Verify renv status and package installation
RUN Rscript -e "cat('Verifying renv environment...\n'); renv::status(); cat('\nEnvironment ready!\n')"

# Set working directory
WORKDIR /analysis

# Default command runs the full pipeline
CMD ["Rscript", "All_Run/run.R"]
