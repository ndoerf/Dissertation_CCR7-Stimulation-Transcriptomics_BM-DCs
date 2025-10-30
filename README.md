# Dissertation_CCR7-Stimulation-Transcriptomics_BM-DCs
A reproducible bioinformatics pipeline (using R and Bioconductor) for the transcriptomic analysis of murine BM-DC following CCR7 stimulation.

Project Overview

This repository contains the R analysis scripts and associated files used to generate the figures and statistical results presented in Chapter V of the submitted PhD thesis, "CCR7-Mediated Dendritic Cell Migration: A Multilevel Analysis of Metabolic and Regulatory Mechanisms"

Data Source and Experimental Design

The data analysied here is derived from Next-Generation Sequencing (NGS) experiments on mature bone marrow-derived dendritic cells (BM-DC) stimulated with the CCR7 ligand CCL21 for 24 h or 48 h.

Experimental Groups:

The mature BM-DCs were prepared and treated according to the following 2x2 factorial design:
Stimulus for Maturation	Chemokine Stimulation	Harvest Time (After Chemokine)
CpG		500 ng/mL CCL21 or Unstimulated	24 h or 48 h
IFN-γ		500 ng/mL CCL21 or Unstimulated	24 h or 48 h
CpG+IFN-γ	500 ng/mL CCL21 or Unstimulated	24 h or 48 h

Data Availability:

Due to the substantial size of the NGS data and its continued use for ongoing projects within the lab, the raw NGS data files are not included in this public repository.

The final, processed dataset (the gene expression matrix) necessary to execute the R script and reproduce the figures in Chapter V is available to academic reviewers and researchers upon formal request to the corresponding author, subject to a reasonable Material Transfer Agreement (MTA) or Data Use Agreement (DUA) to protect the integrity of the lab’s separate, unpublished work.

Please contact the corresponding author Nicole Doerffer at nicole.doerffer@uni-bonn.de for access.
Repository Structure and Analysis


Script Usage:

The R scripts were used to perform the DESeq2 analysis and to generate the figures presented in Chapter V of the manuscript.

Contact

    Author: Nicole Doerffer

    Contact: nicole.doerffer@uni-bonn.de
