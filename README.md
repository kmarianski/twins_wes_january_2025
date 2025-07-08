# Snakemake pipeline for WES of Twins cohort

## Overview
This repository contains a Snakemake workflow for processing whole exome sequencing (WES) data from a twins cohort. The pipeline performs quality control, trimming, alignment, variant calling, and annotation.

---

## Checklist: First-Time Execution on an Interactive SLURM Node

1. **Clone the repository and navigate to the project directory:**
   ```bash
   git clone <repo_url>
   cd twins_wes_january_2025
   ```

2. **Check and populate the `raw/` directory:**
   - Ensure your raw FASTQ files are placed in the `raw/` directory.
   - File naming should follow the pattern: `sampleID_laneInfo_1.fq.gz` and `sampleID_laneInfo_2.fq.gz` (e.g., `s10003_EKDN210017711-1A_HNT5FDSX2_L4_1.fq.gz`).

3. **Check reference and resource files:**
   - Ensure all paths in `config.yaml` (e.g., reference genome, adapter, annovar, known_sites) are correct and files are accessible from the compute node.

4. **Create the conda environment:**
   ```bash
   conda env create -f envs/biotools.yml
   # Or, if already created:
   conda activate biotools
   ```

5. **Check Singularity image (for DeepVariant):**
   - Ensure `singularity/deepvariant_latest.sif` exists and is accessible.

6. **Check SLURM profile configuration:**
   - Review and edit `profile/config.yaml` if needed for your cluster.

7. **Start an interactive SLURM session:**
   ```bash
   srun --pty -p <partition> -c <cpus> --mem=<mem>G -t <time> /bin/bash
   # Example:
   srun --pty -p long -c 8 --mem=32G -t 24:00:00 /bin/bash
   ```

8. **Activate the conda environment:**
   ```bash
   conda activate biotools
   ```

9. **Run Snakemake:**
   ```bash
   snakemake --profile profile --use-conda --cores <n>
   # Example:
   snakemake --profile profile --use-conda --cores 8
   ```

---

## Troubleshooting
- **Missing files or directories:**
  - Ensure all required files and directories exist as specified in `config.yaml`.
- **Conda environment issues:**
  - Recreate the environment: `conda env remove -n biotools && conda env create -f envs/biotools.yml`
- **Singularity errors:**
  - Check that the image path is correct and accessible.
- **Cluster errors:**
  - Adjust resources in `profile/config.yaml` or your SLURM command as needed.

---

## Pipeline Structure
- `Snakefile`: Main workflow definition
- `config.yaml`: Pipeline configuration
- `envs/biotools.yml`: Conda environment definition
- `profile/`: Snakemake profile for SLURM
- `raw/`: Place your raw FASTQ files here
- `results/`: Output directory (created by pipeline)
- `logs/`: Log files (created by pipeline)

---

## Contact
For questions or issues, contact the pipeline maintainer.