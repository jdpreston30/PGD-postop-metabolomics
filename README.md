# Oxidative Stress-Dominant Metabolomic Signatures Characterize Primary Graft Dysfunction After Heart Transplantation

**Reproducible analysis code for academic publication**

## 📖 Citation

This code is associated with the analysis presented in the following manuscript:
> Rust et al. (2025). Oxidative Stress-Dominant Metabolomic Signatures Characterize Primary Graft Dysfunction After Heart Transplantation. *EJCTS*. (Submitted)

## 🚀 Quick Start for Reproduction

### Option 1: Using Docker (Recommended for Exact Reproducibility)

**Prerequisites**: 
- Install [Docker Desktop](https://www.docker.com/products/docker-desktop)
- (Optional) Create free [Docker Hub](https://hub.docker.com) account

#### Method A: Pull Pre-built Image (Fastest - Recommended)

```bash
# 1. Clone the repository
git clone https://github.com/jdpreston30/PGD-postop-metabolomics.git
cd PGD-postop-metabolomics

# 2. Pull the pre-built Docker image (~5-10 minutes)
docker pull jdpreston30/pgd-postop-metabolomics:latest

# 3. Run the complete analysis pipeline
docker run --rm -v $(pwd)/Outputs:/analysis/Outputs jdpreston30/pgd-postop-metabolomics:latest
```

#### Method B: Build Image Locally

```bash
# 1. Clone the repository
git clone https://github.com/jdpreston30/PGD-postop-metabolomics.git
cd PGD-postop-metabolomics

# 2. Build the Docker image from Dockerfile (~30-45 minutes)
export DOCKER_BUILDKIT=1
docker build -t pgd-postop-metabolomics .

# 3. Run the complete analysis pipeline
docker run --rm -v $(pwd)/Outputs:/analysis/Outputs pgd-postop-metabolomics
```

#### What's Included

The Docker container provides a completely isolated, reproducible environment with:
- **R 4.5.1** (2025-06-13) with all required packages at pinned versions
- **Package management**: Uses `renv` to lock exact versions of 57 R packages
- **Bioconductor 3.21**: mixOmics, limma, xcms, CAMERA, fgsea, and more
- **System dependencies**: Pandoc, ImageMagick, GraphViz, GDAL, HDF5, LaTeX
- **Guaranteed identical results** regardless of host system or when run
- **Key packages**: mixOmics 6.32.0, limma, ggplot2 3.5.2, igraph, MetaboAnalystR
- **Automatic loading**: All dependencies from `DESCRIPTION` loaded via helper functions

All outputs (figures, tables, pathway results) will be saved to your local `Outputs/` directory.

#### Testing the Container

To verify the Docker image before running the full analysis:

```bash
# Quick verification (< 1 minute)
docker run --rm pgd-postop-metabolomics Rscript -e "packageVersion('mixOmics'); packageVersion('limma')"
```

This should display package versions. If successful, the container is ready.

#### Troubleshooting

**Build fails or is very slow**: 
- First build takes 30-45 minutes (compiling packages like xcms, mixOmics)
- If interrupted, run `docker build` again - it resumes from cached layers

**"No space left on device"**:
- Final image is ~3-5 GB
- Check space: `docker system df`
- Clean up: `docker image prune`

**No output files created**:
- Verify mount: `-v $(pwd)/Outputs:/analysis/Outputs`
- Check permissions: `ls -la Outputs/`

### Option 2: Manual Installation (Without Docker)

**Prerequisites**: 
- R >= 4.5.1
- Git

**Note**: This project uses `renv` for package management. The `renv.lock` file contains exact versions of all packages used in the manuscript.

```r
# 1. Clone the repository (from terminal)
git clone https://github.com/jdpreston30/PGD-postop-metabolomics.git
cd PGD-postop-metabolomics

# 2. Start R in the project directory
# (renv automatically activates via .Rprofile)

# 3. Restore all packages at exact versions (first time: ~10-20 minutes)
renv::restore()

# 4. Check system dependencies
source("R/Utilities/Helpers/check_system_dependencies.R")
check_system_dependencies()

# 5. Run the complete analysis pipeline
source("All_Run/run.R")
```

**What happens during `renv::restore()`**:
- Installs 57 R packages at exact versions from `renv.lock`
- CRAN packages (ggplot2, dplyr, purrr, tidyr, etc.)
- Bioconductor 3.21 packages (limma, mixOmics, xcms, CAMERA, fgsea, etc.)
- GitHub packages (TernTablesR, MetaboAnalystR)
- Creates isolated project library (doesn't affect system R packages)
- Only needed once per computer
- Packages auto-loaded from `DESCRIPTION` during pipeline execution

## 📁 Project Structure

```
├── DESCRIPTION              # R package dependencies (CRAN, Bioconductor, GitHub)
├── Dockerfile              # Docker container for reproducibility
├── renv.lock               # Exact package versions
├── All_Run/
│   ├── config_dynamic.yaml # Analysis configuration
│   └── run.R              # Main pipeline script
├── R/
│   ├── Scripts/           # Analysis workflow (00a-08)
│   └── Utilities/         # Custom functions
│       ├── Analysis/      # Statistical and pathway analysis
│       ├── Helpers/       # Utility functions
│       ├── Preprocessing/ # Data cleaning
│       └── Visualization/ # Plotting
├── Databases/
│   ├── IDX_IROA/         # IROA quantification standards
│   └── MetaboAnalystR/   # KEGG pathway databases
├── Outputs/              # Generated results
│   ├── Balloon_and_Volcano/
│   ├── Enrichment/       # Pathway analysis
│   ├── Heatmaps/
│   ├── LIMMA/           # Differential analysis
│   ├── PCA/
│   └── PLSDA/
└── Supporting Information/
```

## 🔬 Analysis Workflow

The pipeline executes in sequence:

1. **00a-00d**: Environment setup, clinical metadata, feature tables
2. **01**: Data cleanup and QC filtering
3. **02**: PCA, PLS-DA, heatmap generation
4. **03**: LIMMA differential analysis (group, time, interaction)
5. **04**: Pathway enrichment (mummichog)
6. **05**: Targeted metabolite balloon/volcano plots
7. **06**: Subject-based trajectories
8. **07**: Results number compilation
9. **08**: Supplemental materials

## 💻 System Requirements

### Computational Requirements
- **R**: Version 4.5.1 or higher
- **Platform**: Developed on macOS; compatible with Windows/Linux
- Standard modern computer sufficient

### System Dependencies
- **Pandoc**: R Markdown rendering
- **ImageMagick**: Image processing
- **GraphViz**: Network visualization (Rgraphviz)
- **TinyTeX/LaTeX**: PDF generation

*All system dependencies automatically installed in Docker container. For manual installation, run `check_system_dependencies()` for platform-specific instructions.*

## 📦 Package Dependencies

All R package dependencies specified in `DESCRIPTION`. Key packages:

### CRAN Packages
- **Data manipulation**: tidyverse (dplyr, tidyr, purrr, readr, stringr)
- **Visualization**: ggplot2, ggraph, patchwork, viridis, Cairo
- **Statistical modeling**: mixOmics, caret, randomForest, e1071
- **Reporting**: rmarkdown, knitr, officer, flextable

### Bioconductor 3.21
- **Metabolomics**: xcms, CAMERA, MSnbase
- **Differential analysis**: limma
- **Pathway analysis**: fgsea, globaltest, GlobalAncova
- **Network analysis**: RBGL, Rgraphviz

### GitHub Packages
- `jdpreston30/TernTablesR`: Epidemiologic tables
- `xia-lab/MetaboAnalystR`: Metabolomics and pathway enrichment

*See `DESCRIPTION` for complete list.*

## 🔄 Reproducibility Features

Best practices implemented:

- ✅ **Version Control**: Complete code on GitHub
- ✅ **Package Management**: `renv` with `renv.lock` pinning 57 packages
- ✅ **Dependency Declaration**: All dependencies in `DESCRIPTION` with auto-loading
- ✅ **Containerization**: Docker image with R 4.5.1, Bioconductor 3.21
- ✅ **Conflict Resolution**: `conflicted` package for predictable behavior
- ✅ **Configuration-Driven**: All parameters in `config_dynamic.yaml`
- ✅ **System Validation**: Automated dependency checking
- ✅ **Documentation**: Function documentation and workflow comments
- ✅ **Session Info**: Exact package versions documented

## 🤝 For Reviewers & Collaborators

**Easiest reproduction**: Pull the pre-built Docker image from Docker Hub (Option 1A above). This ensures exact computational environment used for all manuscript results.

### Quick Verification (5 minutes)
```bash
# Clone, pull image, verify
git clone https://github.com/jdpreston30/PGD-postop-metabolomics.git
cd PGD-postop-metabolomics
docker pull jdpreston30/pgd-postop-metabolomics:latest
docker run --rm jdpreston30/pgd-postop-metabolomics:latest Rscript -e "packageVersion('limma')"
# Should output version number
```

### Full Analysis Run (~30-60 minutes)
```bash
docker run --rm -v $(pwd)/Outputs:/analysis/Outputs jdpreston30/pgd-postop-metabolomics:latest
```

If issues occur:
1. Verify Docker Desktop is running
2. Ensure you're in repository root directory
3. Check output directories are writable
4. Review logs: `docker logs <container-id>`

For questions, open a GitHub issue or contact the corresponding author.

## 📧 Contact

**Senior Author**: Joshua Chan, MD, PhD  
- **Email**: joshua.chan@emory.edu  
- **ORCID**: [0000-0001-7220-561X](https://orcid.org/0000-0001-7220-561X)  
- **Institution**: Department of Surgery, Emory University School of Medicine

**Corresponding Author**: Joshua Preston  
- **Email**: joshua.preston@emory.edu  
- **ORCID**: [0000-0001-9834-3017](https://orcid.org/0000-0001-9834-3017)  
- **Institution**: Department of Surgery, Emory University School of Medicine

---

**Repository**: https://github.com/jdpreston30/PGD-postop-metabolomics
