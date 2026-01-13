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
- **Package management**: Uses `renv` to lock exact versions of 246+ R packages
- **Bioconductor 3.21**: 16 packages including mixOmics, limma, xcms, CAMERA, fgsea, RBGL, Rgraphviz
- **GitHub packages**: TernTablesR, MetaboAnalystR (both automatically installed from lockfile)
- **System dependencies**: Pandoc, ImageMagick, GraphViz, GDAL, HDF5, LaTeX
- **Guaranteed identical results** regardless of host system or when run
- **Key packages**: mixOmics 6.32.0, limma 3.64.3, ggplot2 4.0.1, igraph 2.2.1, MetaboAnalystR 4.2.0
- **Automatic loading**: All 79 dependencies from `DESCRIPTION` loaded via helper functions

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
# renv automatically activates and prompts to restore on first launch
renv::restore()

# 4. Check system dependencies
source("R/Utilities/Helpers/check_system_dependencies.R")
check_system_dependencies()

# 5. Run the complete analysis pipeline
source("All_Run/run.R")
```

**What happens during `renv::restore()`**:
- Installs 246+ R packages at exact versions from `renv.lock`
- CRAN packages: ggplot2 4.0.1, dplyr, tidyr, purrr, data.table, Cairo, and 200+ more
- Bioconductor 3.21 packages (16 total): limma 3.64.3, mixOmics 6.32.0, xcms 4.6.4, CAMERA 1.64.0, fgsea, RBGL, Rgraphviz, and 9 more
- GitHub packages (2): jdpreston30/TernTablesR, xia-lab/MetaboAnalystR 4.2.0
- Creates isolated project library (doesn't affect system R packages)
- Only needed once per computer (or after deleting renv/library)
- All 79 packages from `DESCRIPTION` auto-load during pipeline execution

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

The pipeline executes 12 scripts in sequence:

1. **00a_environment_setup**: Load 79 packages, set conflict preferences, verify dependencies
2. **00b_setup**: Load utility functions from R/Utilities/
3. **00c_clinical_metadata**: Process patient demographics and clinical data
4. **01_FTs**: Generate feature tables, apply QC filters (80% uniqueness threshold)
5. **02_PCA_PLSDA_heatmaps**: Multivariate analysis with permutation testing (1000 permutations)
6. **03_limma**: Differential analysis for PGD effect, time effect, and interaction
7. **04_pathway_enrichment**: Mummichog pathway analysis (3 contrasts, 100 permutations each)
8. **05_targeted_volcano_diverge**: Volcano and balloon plots for identified metabolites
9. **06_targeted_subject_based**: Individual subject trajectory plots
10. **07_results_numbers**: Compile statistics for manuscript text
11. **08_assign_figures**: Map plots to final figure panels
12. **09_render_figures**: Generate publication-ready PDFs/PNGs
13. **10_tables**: Create descriptive and results tables (using TernTablesR)
14. **11_supplementary**: Generate supplementary materials
15. **12_session_info**: Document complete session information

**Runtime**: ~5-10 minutes on standard laptop (excluding initial package installation)

**Outputs**: All results saved to `Outputs/` with organized subdirectories for figures, tables, pathway analyses

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
- **� Data Availability

**Note**: The raw metabolomics data underlying this analysis are not included in this repository due to patient privacy protections and institutional data sharing agreements.

**For Researchers**:
- De-identified data may be available upon reasonable request to the corresponding author
- Requests will be evaluated on a case-by-case basis following institutional review and appropriate data use agreements
- This code repository demonstrates the complete analysis pipeline and can be adapted for similar metabolomics datasets

**Repository Purpose**: This repository provides:
- ✅ Complete, reproducible analysis code
- ✅ All custom functions and utilities
- ✅ Exact computational environment (via renv + Docker)
- ✅ Analysis workflow documentation
- ✅ Example database structures (IDX_IROA, KEGG pathways)

## �Differential analysis**: limma
- **Pathway analysis**: fgsea, globaltest, GlobalAncova
- **Network analysis**: RBGL, Rgraphviz

### GitHub Packages
- `jdpreston30/TernTablesR`: Epidemiologic tables
- `xia-lab/MetaboAnalystR`: Metabolomics and pathway enrichment

*See `DESCRIPTION` for complete list.*

## 🔄 Reproducibility Features
 for maximum reproducibility:

- ✅ **Version Control**: Complete code on GitHub with full history
- ✅ **Package Management**: `renv` with `renv.lock` pinning 246+ packages at exact versions
- ✅ **Dependency Declaration**: All 79 dependencies in `DESCRIPTION` with auto-loading via helper functions
- ✅ **Containerization**: Docker image with R 4.5.1, Bioconductor 3.21, all system dependencies
- ✅ **Explicit Snapshots**: renv configured with `snapshot.type: explicit` for complete control
- ✅ **GitHub Packages**: TernTablesR and MetaboAnalystR captured in lockfile with full definitions
- ✅ **Conflict Resolution**: `conflicted` package for predictable function behavior (32 preferences set)
- ✅ **Configuration-Driven**: All analysis parameters in `config_dynamic.yaml`
- ✅ **System Validation**: Automated dependency checking for Pandoc, ImageMagick, LaTeX
- ✅ **Documentation**: Function documentation with roxygen2, workflow comments following style guide
- ✅ **Session Info**: Complete package versions, system info saved in `session_info.txt`
- ✅ **Automated Testing**: Nuclear test validated (rm -rf renv/library → restore → run pipeline)low comments
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

## 🐛 Troubleshooting

### Docker Issues

**Container fails to start**:
```bash
# Check Docker is running
docker info

# View recent logs
docker ps -a
docker logs <container_id>
```

**No output files generated**:
- Verify volume mount is correct: `-v $(pwd)/Outputs:/analysis/Outputs`
- Check directory permissions: `ls -la Outputs/`
- Ensure you're in repository root: `pwd` should show `.../PGD-postop-metabolomics`

### Manual Installation Issues

**renv::restore() fails on Bioconductor packages**:
```r
# Verify Bioconductor version
BiocManager::version()  # Should be 3.21

# Manual install if needed
BiocManager::install(c("limma", "mixOmics", "xcms", "CAMERA"))
```

**GitHub packages fail to install**:
```r
# Install remotes if needed
install.packages("remotes")

# Manual GitHub package installation
remotes::install_github("jdpreston30/TernTablesR")
remotes::install_github("xia-lab/MetaboAnalystR")

# Then update lockfile
renv::snapshot()
```

**Package conflicts or loading errors**:
```r
# Nuclear option: delete library and restore
unlink("renv/library", recursive = TRUE)

# Restart R, then:
renv::restore()
```

**System dependencies missing** (Linux/macOS):
```bash
# Check what's needed
Rscript -e "source('R/Utilities/Helpers/check_system_dependencies.R'); check_system_dependencies()"

# Ubuntu/Debian example
sudo apt-get install pandoc imagemagick graphviz

# macOS example
brew install pandoc imagemagick graphviz
```

If issues persist:
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
