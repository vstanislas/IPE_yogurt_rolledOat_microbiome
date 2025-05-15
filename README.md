

Repository for the code and data of the results presented in Thriene et all. *"Impact of Yogurt and Rolled Oats on the Gut Microbiome: A Randomized Crossover Study displaying Individual Responses and General Resilience"* (under revision).


# Structure of the repository
- `DataAnalysis`: contains code and results for all analyses presented in the manuscript
- `DataPreparation`: contains all data files and code used to prepare data for analysis
- `FiguresPaper`: code used to generate the figures shown in the manuscript
- `R functions`: a collection of R scripts with helper functions used throughout the data analysis

# DATA files

## Metagenomic data
Raw metagenomic sequencing data have been deposited in the NCBI Sequence Read Archive (SRA) under accession number [PRJNA125884](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1258884)

Processed data utilized in this study are stored in the `DataPreparation/Microbiome/` directory within the following files:

- Taxonomic abundance table at SGB level: `Taxonomy/merged_abundance_table_ALL_SGB.txt`
- Taxonomic abundance table at species level: `Taxonomy/merged_abundance_table_ALL_species.txt`
- Genes content table: `Functional/genefamilies_exp86_440_cpm_unstratified_10prev_1cpm.tsv`
- Pathways table: `Functional/pathways_all_704_samples_cpm.tsv`  
- Pathways table for participant JH144: `Functional/JH144_pathabundance_cpm.tsv`
- KEGGs table: `Functional/ko_merged_440samples_1cpm_10prev.tsv`
- GO terms table : `Functional/go_merged_440samples_1cpm_10prev.tsv`



## Targeted analysis of SCFA
The processed metabolite concentration table and the associated analysis report have been deposited on Zenodo (https://doi.org/10.5281/zenodo.15363886). 

The table is also available in this repository under `DataPreparation/Metabolomic/FREI-0301-21TASA CDT.xlsx`



## Blood markers measurements

Blood markers measurements are stored in the `DataPreparation/BloodMarkers/` directory within the following files:

- Zonulin: `22-026 Michels 8 serum.xlsx`
- TNFR2, Interleukin-6, sRAGE, 8-OHdG: `22-027 Michels 9 plasma.xlsx`
- Fructosamine, CRP: `JH_Serum_Zentrallabor_2022-08-19.xlsx`



       