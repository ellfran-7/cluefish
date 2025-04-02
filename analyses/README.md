# 📄Content

This directory contains structured analysis scripts organised in organism and chemical exposure folders, corresponding to different datasets, and contains in this instance three main folders:

-  `drerio-dbp`: A folder containing pipeline scripts associated to the main dataset from Franklin et al. (submitted). This contains:
   - `dromics_dre_transcriptomic_pipeline.R`: Step-by-step `DRomics` pipeline for processing dose-response transcriptomic raw count data. Generates input data for the Cluefish workflow case study.
   - `dromics_dre_anchoring_pipeline`: Step-by-step `DRomics` pipeline for processing dose-response anchoring measurements (body length, eye surface) providing experimental validation for Cluefish workflow results.
   - `standard_approach_dre_pipeline.R`: Standard functional enrichment pipeline for benchmarking against Cluefish workflow outputs, using the deregulated transcript list derived from the `DRomics` pipeline.

*Note: The Cluefish analysis for this main dataset is implemented in the root make.R script rather than as a separate file in this directory. This dataset provided the building blocks for the Cluefish workflow.*

-  `pcanadensis-phe`: A folder containing pipeline scripts associated to an external dataset (**GSE263776**) aimed to provide further validation behind the Cluefish workflow in the Franklin et al. (submitted) paper. This contains:
   - `dromics_pca_transcriptomic_pipeline.R`: Step-by-step `DRomics` pipeline for processing dose-response transcriptomic raw count data. Generates input data for the Cluefish workflow case study.
   - `standard_approach_pca_pipeline.R`: Standard functional enrichment pipeline for benchmarking against Cluefish workflow outputs, using the deregulated transcript list derived from the `DRomics` pipeline.
   - `cluefish_approach_pca_pipeline.R`: Complete Cluefish pipeline implementation for this dataset.
-  `rnorvegicus-pfoa`: A folder containing pipeline scripts associated to an external dataset (**GSE147072**) aimed to provide further validation behind the Cluefish workflow in the Franklin et al. (submitted) paper. This contains:
   - `dromics_rno_transcriptomic_pipeline.R`: Step-by-step `DRomics` pipeline for processing dose-response transcriptomic raw count data. Generates input data for the Cluefish workflow case study.
   - `standard_approach_rno_pipeline.R`: Standard functional enrichment pipeline for benchmarking against Cluefish workflow outputs, using the deregulated transcript list derived from the `DRomics` pipeline.
   - `cluefish_approach_rno_pipeline.R`: Complete Cluefish pipeline implementation for this dataset.

These pipelines are formatted to enable the reproduction of the results required for the Cluefish workflow and those published in Franklin et al. (submitted). However, they can also be used as adaptable templates for your own case studies. 

For help with the `DRomics` R package, visit the [DRomics documentation](https://lbbe-software.github.io/DRomics/).

An additional directory named `scripts-for-paper-content` holds additional file scripts, specifically to generate most of the figures and tables in the paper associated to Cluefish (Franklin et al. submitted). 

## 📍Note

The suggested content and structure here follow best practices for research compendiums. However, feel free to organize and store your project files in a way that best suits your needs and workflow.
