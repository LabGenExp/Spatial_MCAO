# Spatiotemporal transcriptomic map of glial cell response in mouse model of acute brain ischemia

Zucha et al. 2024


Hello, 
Welcome to the repository for the codebase used for data processing, analysis and visualization in our manuscript titled "Spatiotemporal transcriptomic map of glial cell response in mouse model of acute brain ischemia". In this study, we used spatial and single-cell transcriptomics to map out processes governing the response to ischemic stroke during the first week after the injury. Focusing on glial cells, we documented their activation, formation of glial scar and changes in cell-cell communication patterns. Here, we list the major sections of the analysis and links to for data download.

## Interactive ST Data Exploration
Our spatial dataset can be interactively explored at [Nygen portal](https://scarfweb.nygen.io/eu-central-1/public/xv2x2szz). It provides a UMAP or Spatial view of the plots, allowing to browse discussed biological themes of the paper - brain and lesion regions, gene ontology term enrichment, cell type content, ligand-receptor expression, gene ontology processes affected downstream of Apoe-Trem2 signaling, and expression signatures for reactive glial populations.

## Data Availability
- Raw sequencing data for both ST and single-nucleus / cell experiments are available for download at [GEO GSE233815](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233815).
- Seurat objects for the ST and single-nucleus / cell experiments are available for download at [Mendeley data](https://data.mendeley.com/preview/gnb2dsjms2?a=1e744314-eb08-4c66-abe5-e3885b8415c7).


## Chapter 1 - Ischemic brain injury severely disrupts cortical gene expression landscape
We supply scripts on the processing and visualization for quality control of the control section, integration with ischemic sections, differential gene expression analysis and gene set enrichment.

- [Processing and brain region QC for the control section.](PNAS_Review/1DP_01_MCAO_Ctrl_spatial.Rmd)
- [Integration of ST data across batches.](PNAS_Review/1DP_02_MCAO_Integration.Rmd)
- [Differential gene expression analysis and gene ontology.](PNAS_Review/1DP_03_MCAO_DEGs_and_GeneOntology.Rmd)
- [Figure 1 Visualization](PNAS_Review/2V_01_Fig1_SpatialOverview.Rmd)

## Chapter 2 - Ischemic lesion perturbs cortical processes and cellular composition
To provide characterization of the perturbed state through the lense of molecular processes and cellular composition. We provide markdowns for processing and visualization for reference-based deconvolution.

- [Reference-based deconvolution of ST data employing RCTD algorithm.](PNAS_Review/1DP_04_MCAO_deconvolution.Rmd)
- [Figure 2 Visualization](PNAS_Review/2V_02_Fig2_gsea_deconvolution.Rmd)

## Chapter 3 - Cell-cell communication analysis identifies increased glia-oriented crosstalk in the lesion periphery
To decompose the interactions between cells in the lesion periphery. We provide markdowns for processing and visualization the cell-cell interactions in the lesion peripheries and bulk RNA-seq validations.

- [Cell-cell interactions in the entire cortex (for ctrl section) or full lesions.](PNAS_Review/1DP_05_MCAO_spatial_cell_cell_interactions_SpaTalk.Rmd)
- [Cell-cell interactions for the lesion periphery.](PNAS_Review/1DP_06_MCAO_spatial_cell_cell_interactions_SpaTalk_periphery.Rmd)
- [Validation using bulk data.](PNAS_Review/1DP_07_MCAO_bulk.Rmd)
- [Figure 3 Visualization](PNAS_Review/2V_03_Fig3_CCI_Apoe_Trem2.Rmd)

## Chapter 4 & 5 - A spectrum of reactive glia emerges after ischemic brain injury
To supplement the findings in ST data for presence of reactive glia, we performed single-nucleus profiling followed by an oligodendrocyte-enriched single-cell profiling. We provide markdowns for single-nucleus and single-cell processing, background RNA clean up, their integration, analysis, immunohistochemistry validation and visualization.

- [SN processing](PNAS_Review/1DP_08_MCAO_snRNA_preprocessing.Rmd)
- [SN background RNA removal (soupX)](PNAS_Review/1DP_09_MCAO_snRNA_soupX.Rmd)
- [SC processing](PNAS_Review/1DP_10_MCAO_scRNA_preprocessing.Rmd)
- [Integration and respective cell type analysis](PNAS_Review/1DP_11_MCAO_scsn.Rmd)
- [Figure 4 & 5 visualization](PNAS_Review/2V_04_Fig04_05_MCAO_single_nucleus.Rmd)

## Chapter 6 - Glial prominence in the ischemic lesion periphery is recapitulated in published ST datasets
To robustify our findings, through analysis as well as number of replicates, we re-analyzed two publicaly available ischemic ST datasets: [Han et al 2024](https://www.science.org/doi/10.1126/scitranslmed.adg1323) and [Scott et al 2024](https://www.nature.com/articles/s41467-024-45821-y). We looked at the lesino separation into core-periphery, temporal and spatial localization of the expression patterns associated with glial accumulation, and reactive glia expression signatures.

- [Spatial meta-analysis processing](PNAS_Review/1DP_12_MCAO_SpatialMetanalysis.Rmd)
- [Figure 6 visualization](PNAS_Review/2V_05_Fig06_MCAO_Metanalysis.Rmd)