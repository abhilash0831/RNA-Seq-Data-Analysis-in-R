# Copper Stress Transcriptomics in Desulfovibrio alaskensis G20

R scripts for transcriptomic data analysis and scientific visualization associated with our published study of copper stress in the sulfate-reducing bacterium *Desulfovibrio alaskensis* G20.

## Publication

Tripathi, A. K., Saxena, P., et al. (2022). **Transcriptomics and Functional Analysis of Copper Stress Response in the Sulfate-Reducing Bacterium Desulfovibrio alaskensis G20.** *International Journal of Molecular Sciences*, 23(3), 1396.

[Read the publication](https://doi.org/10.3390/ijms23031396) · [PubMed](https://pubmed.ncbi.nlm.nih.gov/35163324/)

## Research Overview

Copper is an essential micronutrient, but elevated concentrations can disrupt cellular functions. This study investigated how *D. alaskensis* G20 responds to copper exposure by combining RNA sequencing with growth measurements and RT-PCR analysis.

The transcriptomic analysis compared three conditions:

* Untreated control versus 5 µM Cu(II).
* Untreated control versus 15 µM Cu(II).
* 5 µM versus 15 µM Cu(II).

The published findings describe changes in gene expression associated with ion transport, translation, transcription, and signal transduction, providing a basis for further investigation of bacterial copper-stress responses.

## Code in This Repository

The scripts support downstream analysis of DESeq2 results generated in Galaxy and visualization of transcriptomic and experimental data.

| Script                                     | Purpose                                                                                              |
| ------------------------------------------ | ---------------------------------------------------------------------------------------------------- |
| `Deseq2 Data cleaning and Fig 2, S2, S3.R` | Processing differential-expression tables, comparing gene sets, and visualizing expression patterns. |
| `Bubble_Plot_Figure 3.R`                   | Bubble-plot visualization.                                                                           |
| `ChordPlot_Figure 4.R`                     | Chord-plot visualization.                                                                            |
| `GO_BarPlot_Fig S5.R`                      | Gene Ontology bar-plot visualization.                                                                |
| `Growth Curve _Figure 1A,1B.R`             | Growth-curve visualization.                                                                          |
| `RTPCR_BarGraph_Fig S4.R`                  | RT-PCR result visualization.                                                                         |

Figure references follow the script filenames.

## Tools and Methods

* **R:** Data processing, analysis, and visualization.
* **Galaxy / DESeq2:** Upstream differential-expression analysis.
* **tidyverse and readxl:** Data import and manipulation.
* **ggplot2, ggvenn, GOplot, EnhancedVolcano, and pheatmap:** Scientific visualization.

## Using the Scripts

These are research scripts associated with the publication. Individual scripts list their required R packages and reference input tables by filename.

Input tables are not bundled with this repository. Raw-read processing and the upstream Galaxy workflow are outside the scope of the code provided here.

For experimental design, detailed methods, results, and supplementary materials, consult the [published article](https://doi.org/10.3390/ijms23031396).

## Citation

If this work informs your research, please cite:

Tripathi, A. K., Saxena, P., et al. (2022). Transcriptomics and Functional Analysis of Copper Stress Response in the Sulfate-Reducing Bacterium Desulfovibrio alaskensis G20. *International Journal of Molecular Sciences*, 23(3), 1396. https://doi.org/10.3390/ijms23031396
