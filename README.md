This repository contains most of the data and code used in the study “Cost-saving, trading, and internalization of externality: three economic strategies underlying amino acid archetypes of human metabolic enzymes”, including scripts for data analysis and mathematical modeling. All scripts have been implemented and tested under `R 4.4.2`, `Python 3.12.7`, and `MATLAB R2023b`. Solving the linear programming problems requires `Mosek Optimization Tools Version 10.1.12`.

## 1-Definition_of_AArchetype
This section provides code for defining and characterizing AArchetypes of human metabolic enzymes. Relevant datasets are stored in the `data/` directory. Key scripts include:

**1-Calculate_AA_composition.py**: calculates amino acid composition fractions of human metabolic enzymes

**2-Calculate_AA_enrichment_scores.py**: computes amino acid enrichment scores for metabolic pathways

**3-Enrichment_analysis_functional_category.R**: performs enrichment analysis of functional categories for AArchetypes

**4-Classify_proteins.R**: classifies proteins into predefined categories

**5-PDB_to_3Di_sequence.R**: converts PDB structures of metabolic enzymes into 3Di sequences

**6-Corr_AADistanceMatrix~3DiDistanceMatrix.R**: evaluates the correlation between the amino-acid distance matrix (AADistanceMatrix) and the 3Di distance matrix (3DiDistanceMatrix)

## 2-AArchetype_in_other_organisms
This section performs MDS dimensionality reduction on amino acid composition profiles of metabolic pathways across 10 organisms using MDS_AA_comp_organisms.R. Relevant data are stored in the `data/` directory.

## 3-Calculation_of_proteinCost_and_energyCost
This section contains scripts for calculating protein cost and energy cost, located under `1-Protein_cost/` and `2-Energy_cost/`, respectively. Relevant data are stored in the `data/` directory.

## 4-Analysis_about_dietary_patterns
This section compares amino acid composition profiles between human metabolic pathways and dietary patterns:

**1-tSNE_of pathways+diets.R**: performs t-SNE embedding of amino acid compositions for pathways and diets

**2-Dotplot_of_distance_pathways_diets.R**: generates a bubble plot showing the median distance in amino acid composition between human dietary patterns and AArchetypes

Relevant data are stored in the `data/` directory.

## 5-Analysis_about_IDR_and_LLPS
This section includes scripts for GSEA (GSEA_metabolic_enzymes.R, GSEA_proteome.R) and provides the following datasets under the `data/` directory:

**1-human_protein_abundance.csv**

**2-human_proteome_amino_acid_sequences.csv**

**3-protein_cost_and_energy_cost_human_metabolic_enzymes.csv**

**4-protein_cost_and_energy_cost_human_proteome.csv**

**5-IDR_annotations_human_proteome.csv**

**6-subcellular_localization_annotations_human_proteome.csv**

## 6-LLPS_economics_mathematical_model
This section contains the implementation of the LLPS economics mathematical model (**LLPS_economics_mathematical_model.m**). Relevant data are stored in `data/` directory.