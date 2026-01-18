# RNA-seq Analysis Pipeline with Reproducible Nix Environment

[![Reproducible with Nix](https://img.shields.io/badge/Reproducible-Nix-blue?logo=nixos)](https://nixos.org)
[![GitHub Actions CI](https://github.com/Sowmya-Harsh/RNA-Seq_RAP/actions/workflows/ci.yml/badge.svg)](https://github.com/Sowmya-Harsh/rnaRAP/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains a complete, reproducible RNA-seq differential expression analysis pipeline using DESeq2. The entire analysis environment is managed with **Nix**, ensuring that anyone can run this analysis and get **identical results**, regardless of their operating system or local dependencies.

**Key Features:**
- ✅ Fully reproducible analysis environment with pinned dependencies (managed by `rix`)
- ✅ One-command execution: `nix-shell --run "Rscript run_pipeline.R"`
- ✅ Automated CI/CD with GitHub Actions
- ✅ Complete documentation and methodology

---

## Why Nix for Reproducibility?

### The Problem We're Solving

Imagine you publish a research paper with computational results. A colleague tries to reproduce your analysis 6 months later and gets **different numbers**. Why? 

Common reasons:
- 📦 Package versions changed on their system
- 🐍 Python/R upgraded with breaking changes
- 🔧 System libraries updated (BLAS, LAPACK versions differ)
- 💻 Different operating system (macOS vs Linux results slightly differ)
- 📝 Missing notes about exact versions you used
- 🎲 Floating-point randomness in different environments

**This is a crisis in computational science.** The National Institutes of Health (NIH) found that ~50% of published research cannot be reproduced due to computational irreproducibility.

### What is Nix?

**Nix is a functional package manager** that solves reproducibility by:

1. **Declarative Configuration**: You specify *exactly* what software you need in `default.nix`
2. **Cryptographic Hashing**: Each package and dependency is hashed; identical inputs always produce identical outputs
3. **Immutable Environment**: Your analysis always runs with the same versions, even if the package manager updates globally
4. **Bit-for-Bit Reproducibility**: Multiple runs on different machines produce identical results (important for scientific publishing)

### Key Advantages for Science

| Aspect | Traditional | Nix |
|--------|-----------|-----|
| **Version Control** | ❌ Manual tracking | ✅ Automatic, cryptographic |
| **Cross-platform** | ⚠️ OS-specific issues | ✅ Same on macOS, Linux, Windows |
| **5-year reproducibility** | ❌ Usually broken | ✅ Works identically |
| **Dependency conflicts** | ⚠️ Common ("dependency hell") | ✅ Impossible |
| **Collaboration** | ⚠️ "Works on my machine" | ✅ Works on everyone's machine |
| **Publication** | ⚠️ Results may not replicate | ✅ Mathematically identical results |

---

## Why This Matters in Biomedical Research

### The Reproducibility Crisis

**In biomedicine specifically, computational irreproducibility is a critical issue:**

#### **Scientific Integrity**
- Published computational results must be independently verifiable
- Regulatory agencies (FDA, EMA) increasingly require reproducible analysis
- Clinical decision-making depends on robust, validated statistical methods
- Peer reviewers cannot reliably validate results without reproducibility

#### **Economic Impact**
- A 2015 meta-analysis estimated **$28.4 billion annually** spent on irreproducible research in the US
- Drug development failures often trace back to irreproducible computational findings
- Failed clinical trials from non-reproducible upstream analysis

#### **Clinical Consequences**
- Gene expression signatures used for patient stratification must be reproducible
- Personalized medicine relies on reproducible bioinformatics pipelines
- Biomarker discovery must withstand independent validation

#### **RNA-seq Specific Issues**
RNA-seq analysis has particular reproducibility challenges:
- **Tool versions matter**: STAR, RSEM, DESeq2 version changes can affect results
- **Statistical parameters**: Different R/Bioconductor versions produce different p-values
- **Data processing**: Normalization methods vary between versions
- **System dependencies**: BLAS/LAPACK library versions affect numerical stability
---

## Quick Start

### One-Command Execution

If you have Nix installed, run the entire analysis and generate the report in one command:

```bash
nix-shell --run "Rscript run_pipeline.R"
```

That's it! This will:
1. ✅ Download and configure exact R version (R 4.5.2)
2. ✅ Install all Bioconductor packages (DESeq2, etc.)
3. ✅ Run the complete analysis pipeline
4. ✅ Generate visualizations and tables
5. ✅ Generate the HTML report
6. ✅ Save results to `reports/`

**No manual dependency installation needed.**

### Prerequisites

- **macOS/Linux**: [Install Nix](https://nixos.org/download.html)
- **Windows**: Use [WSL2](https://docs.microsoft.com/en-us/windows/wsl/install) with Linux, or Docker

### Verify Installation

```bash
# Check if Nix is installed
nix --version

# Check if you have default.nix in the project
ls -la default.nix
```

---

## Installation

### Step 1: Install Nix

**macOS:**
```bash
brew install nix
```

**Linux:**
```bash
curl -L https://nixos.org/nix/install | sh
source ~/.profile
```

**Windows (WSL2):**
```bash
# Inside WSL2
curl -L https://nixos.org/nix/install | sh
source ~/.profile
```

### Step 2: Clone This Repository

```bash
git clone https://github.com/Sowmya-Harsh/RNA-Seq_RAP.git
cd rnaRAP
```

---

## Project Structure

```
rnaRAP/
├── .github/
│   └── workflows/
│       └── ci.yml              # GitHub Actions workflow
├── default.nix                 # Nix environment (generated by rix) - THE KEY FILE
├── gen-env.R                   # Script to regenerate default.nix
├── README.md                   # This file
├── run_tests.R                 # Test suite
├── run_pipeline.R              # Main analysis pipeline
├── rnaseq_report.qmd           # Quarto report template
├── analysis_core.R             # Core analysis functions
└── reports/                    # Generated outputs
    ├── rnaseq_report.html      # Main HTML report
    ├── pca_condition.png       # PCA plot
    ├── top_genes.csv           # Top differentially expressed genes
    └── [other outputs...]
```

**Important files:**
- `default.nix`: Automatically generated by `rix`, declares all dependencies with exact versions
- `gen-env.R`: Regenerates `default.nix` if you need to change packages
- `run_pipeline.R`: Orchestrates the entire analysis

---

## Running the Analysis

### Option 1: One-Command Execution (Recommended)

```bash
# Run everything with pinned dependencies
nix-shell --run "Rscript run_pipeline.R"
```

This is the simplest and most reproducible approach.

### Option 2: Interactive Development

```bash
# Drop into the Nix development environment
nix-shell

# Now you can run R directly
Rscript run_pipeline.R

# Or start an R interactive session
R

# Or examine functions before running
Rscript -e "source('analysis_core.R'); head(load_example_data())"

# Exit when done
exit
```

### Option 3: Step-by-Step Execution

```bash
# Enter the Nix environment
nix-shell

# Step 1: Run tests to verify setup
Rscript run_tests.R

# Step 2: Run the pipeline
Rscript run_pipeline.R

# Exit environment
exit
```

---

## Understanding the Pipeline

### The `default.nix` File: Heart of Reproducibility

This project uses **`default.nix`** (managed by the `rix` R package) for reproducible environments. The file contains all pinned dependencies.

**Key points:**
- `pkgs` is pinned to a specific nixpkgs commit (2026-01-14)
- All R packages are at exact versions locked at that date
- The `default.nix` is generated by: `rix(date = "2026-01-14", r_pkgs = c(...))`
- Run from inside `nix-shell -p R rPackages.rix` environment

### The `gen-env.R` File: Generating `default.nix`

This R script generates the `default.nix` file with your exact environment:

```bash
# Generate a new environment spec
nix-shell -p R rPackages.rix

# Inside that shell, run:
Rscript gen-env.R

# This creates/updates default.nix with locked versions
exit

# Now enter the environment with updated dependencies:
nix-shell
```

**Workflow for adding/removing packages:**
1. Edit `gen-env.R` to add/remove R packages
2. Run `Rscript gen-env.R` in a `nix-shell -p R rPackages.rix` environment
3. This updates `default.nix` with pinned versions
4. Exit and run `nix-shell` to enter the new environment
5. Commit both files to git
6. Anyone else can now run `nix-shell` and get identical versions

### Data Flow

```
1. Raw Data (Pasilla data)
        ↓
2. Load data via [load_example_data()]
        ↓
3. DESeq2 Differential Expression Analysis
   - Condition effect analysis
   - Type effect analysis
        ↓
4. Generate Results:
   - Statistics (p-values, log2 fold-changes)
   - Significant genes (padj < 0.05)
        ↓
5. Visualizations [analysis_core.R plotting functions]
        ↓
6. Export Tables [top_genes_table()]
   - Top 20 DE genes and Full results matrix
        ↓
7. Generate Report [rnaseq_report.qmd → quarto_render()]
        ↓
8. Final Output [reports/rnaseq_report.html]
```

---

## GitHub Actions Setup

Automated testing and report generation on every commit!

### CI/CD Pipeline

The `.github/workflows/ci.yml` file automatically:

1. **Installs Nix** on GitHub Actions runner
2. **Creates reproducible environment** using `default.nix`
3. **Runs tests** with `run_tests.R`
4. **Executes pipeline** with `run_pipeline.R`
5. **Verifies outputs** (all expected files exist)
6. **Uploads artifacts** for 30 days
7. **Deploys report** to GitHub Pages

### Workflow Features

CI/CD pipeline:

**On every push to main**: Run full analysis, validate reproducibility  
**On PRs to main**: Test analysis with proposed changes  
**Cross-platform**: Uses standardized Nix environment (runs identically on any CI system)

---

## Troubleshooting

### Common Issues and Solutions

#### Issue: `command not found: nix-shell`

```bash
# Nix not installed
# Solution: Install Nix
curl -L https://nixos.org/nix/install | sh
source ~/.profile
```

#### Issue: `error: file 'default.nix' not found`

```bash
# You're not in the project directory or default.nix is missing
# Solution: Check you're in the right directory
pwd
ls -la default.nix  # Should exist

# If missing, regenerate it:
nix-shell -p R rPackages.rix
Rscript gen-env.R
exit
```

## Modifying the Pipeline

### Adding New R Packages

1. Edit `gen-env.R` - add package name to the list:

```r
# In gen-env.R, find the r_pkgs vector and add your package:
r_pkgs = c("DESeq2", "dplyr", "ggplot2", 
           "YOUR_NEW_PACKAGE",  # ← Add here
           "quarto")
```

2. Regenerate the environment:

```bash
nix-shell -p R rPackages.rix
Rscript gen-env.R
exit

# Now enter the new environment
nix-shell
```

3. Commit the updated `default.nix` to git

### Modifying Analysis Steps

1. Edit `analysis_core.R` for core functions
2. Edit `run_pipeline.R` to change pipeline steps
3. Run locally to test: `nix-shell --run "Rscript run_pipeline.R"`
4. Commit and push to GitHub

---

## Updating Dependencies

To update all packages to the latest compatible versions:

```bash
# Enter the rix environment
nix-shell -p R rPackages.rix

# Edit gen-env.R to change the date:
# Change: date = "2026-01-14"
# To: date = "2026-01-18" (or today's date)

Rscript gen-env.R

exit

# Enter the new environment
nix-shell

# Test the analysis
Rscript run_pipeline.R

# If it works, commit the new default.nix
git add default.nix
git commit -m "Update dependencies"
```

---

## Performance Tips

### Speed Up Initial Build

```bash
# Use the rstats-on-nix binary cache for pre-compiled packages
nix-shell -p nix

# Inside nix-shell:
nix-env -iA cachix -f https://cachix.org/api/v1/install
cachix use rstats-on-nix
exit

# Now subsequent runs are much faster
nix-shell --run "Rscript run_pipeline.R"
```

---

## Advanced: Using with Docker

For even broader compatibility, combine Nix with Docker:

### Build Docker Image

```dockerfile
FROM nixos/nix:latest

RUN mkdir /work
WORKDIR /work

COPY . .

# Build the analysis environment
RUN nix-shell --run "echo 'Environment ready'"

ENTRYPOINT ["nix-shell", "--run", "Rscript run_pipeline.R"]
```

### Run Analysis in Docker

```bash
docker build -t rna-seq-pipeline .
docker run --rm -v $(pwd)/reports:/work/reports rna-seq-pipeline
```

---

## Resources

### Learning Nix and Reproducibility

- [Reproducible Analysis with R in Nix](https://github.com/b-rodrigues/rap4mads_2024) - Full course
- [Nix Manual](https://nixos.org/manual/nix/stable/)
- [Nix Package Manager Guide](https://nixos.wiki/)
- [rstats-on-nix](https://github.com/ropensci/rix) - R + Nix integration

### DESeq2 and RNA-seq

- [DESeq2 Vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [RNA-seq Analysis Course](https://www.bioconductor.org/packages/release/workflows/html/rnaseqGene.html)
- [Modern Statistics for Modern Biology](https://www.huber.embl.de/msmb/)

### Related Tools

- 🐳 [Docker](https://www.docker.com/) - Container alternative to Nix
- [Singularity](https://sylabs.io/) - Another container system
- [Conda](https://conda.io/) - Alternative package manager (less reproducible than Nix)

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Support

Have questions or issues?

-  [Open an issue](https://github.com/Sowmya-Harsh/RNA-Seq_RAP/issues)

## Acknowledgments

- **Reproducible Pipeline Course**: https://github.com/b-rodrigues/rap4mads_2024
- **rstats-on-nix community**: For excellent R + Nix integration tools
- **Bioconductor**: For DESeq2 and RNA-seq analysis tools

---
