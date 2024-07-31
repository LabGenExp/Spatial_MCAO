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
- Processing and brain region QC for the control section.
- Integration of ST data across batches.
- Differential gene expression analysis and gene ontology.
- Figure 1 Visualization

## Chapter 2 - Ischemic lesion perturbs cortical processes and cellular composition
To provide characterization of the perturbed state through the lense of molecular processes and cellular composition. We provide markdowns for processing and visualization for reference-based deconvolution.
- Reference-based deconvolution of ST data employing RCTD algorithm.
- Figure 2 Visualization

## Chapter 3 - Cell-cell communication analysis identifies increased glia-oriented crosstalk in the lesion periphery
To decompose the interactions between cells in the lesion periphery. We provide markdowns for processing and visualization the cell-cell interactions in the lesion peripheries and bulk RNA-seq validations.
- Cell-cell interactions in the entire cortex (for ctrl section) or full lesions.
- Cell-cell interactions for the lesion periphery.
- Validation using bulk data.
- Figure 3 Visualization

## Chapter 4 & 5 - A spectrum of reactive glia emerges after ischemic brain injury
To supplement the findings in ST data for presence of reactive glia, we performed single-nucleus profiling followed by an oligodendrocyte-enriched single-cell profiling. We provide markdowns for single-nucleus and single-cell processing, background RNA clean up, their integration, analysis, immunohistochemistry validation and visualization.
- SN processing
- SN background RNA removal (soupX)
- SC processing
- Integration and respective cell type analysis
- Figure 4 & 5 visualization

## Chapter 6 - Glial prominence in the ischemic lesion periphery is recapitulated in published ST datasets
To robustify our findings, through analysis as well as number of replicates, we re-analyzed two publicaly available ischemic ST datasets: [Han et al 2024](https://www.science.org/doi/10.1126/scitranslmed.adg1323) and [Scott et al 2024](https://www.nature.com/articles/s41467-024-45821-y). We looked at the lesino separation into core-periphery, temporal and spatial localization of the expression patterns associated with glial accumulation, and reactive glia expression signatures.
- Spatial meta-analysis processing
- Figure 6 visualization




