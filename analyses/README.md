# 📄Content

This directory contains structured analysis scripts organised in organism and chemical exposure folders, corresponding to different datasets, and contains in this instance three main folders:

-  `drerio-dbp`: A folder containing pipeline scripts associated to the main dataset from [Franklin *et al.* (2025)](https://doi.org/10.1093/nargab/lqaf103). This contains:
   - `dromics_dre_anchoring_pipeline`: Step-by-step `DRomics` pipeline for processing dose-response anchoring measurements (body length, eye surface) providing experimental validation for Cluefish workflow results.
   - `dromics_dre_transcriptomic_pipeline.R`: Step-by-step `DRomics` pipeline for processing dose-response transcriptomic raw count data. Generates input data for the Cluefish workflow case study.
   - `standard_approach_dre_pipeline.R`: Standard functional enrichment pipeline for benchmarking against Cluefish workflow outputs, using the deregulated transcript list derived from the `DRomics` pipeline.

*Note: The Cluefish analysis for this main dataset is implemented in the root `make.R` script rather than as a separate file in this directory. This dataset provided the building blocks for the Cluefish workflow.*

-  `pcanadensis-phe`: A folder containing pipeline scripts for external dataset validation (**GSE263776**, [Gréau *et al.* 2025](https://doi.org/10.1007/s11356-025-36002-5)). This contains:
   - `dromics_pca_transcriptomic_pipeline.R`: Step-by-step `DRomics` pipeline for processing dose-response transcriptomic raw count data. Generates input data for the Cluefish workflow case study.
   - `standard_approach_pca_pipeline.R`: Standard functional enrichment pipeline for benchmarking against Cluefish workflow outputs, using the deregulated transcript list derived from the `DRomics` pipeline.
   - `cluefish_approach_pca_pipeline.R`: Complete Cluefish pipeline implementation for this dataset.
-  `rnorvegicus-pfoa`: A folder containing pipeline scripts for external dataset validation (**GSE147072**, [Gwinn *et al.* 2020](https://www.ncbi.nlm.nih.gov/pubmed/32492150)). This contains:
   - `dromics_rno_transcriptomic_pipeline.R`: Step-by-step `DRomics` pipeline for processing dose-response transcriptomic raw count data. Generates input data for the Cluefish workflow case study.
   - `standard_approach_rno_pipeline.R`: Standard functional enrichment pipeline for benchmarking against Cluefish workflow outputs, using the deregulated transcript list derived from the `DRomics` pipeline.
   - `cluefish_approach_rno_pipeline.R`: Complete Cluefish pipeline implementation for this dataset.

These pipelines enable reproduction of results published in [Franklin *et al.* (2025)](https://doi.org/10.1093/nargab/lqaf103) and can be used as adaptable templates for your own case studies. 

For help with the `DRomics` R package, visit the [DRomics documentation](https://lbbe-software.github.io/DRomics/).

An additional directory named `scripts-for-paper-content` holds scripts to generate figures and tables from the associated paper.

## 📍Note

The suggested content and structure here follow best practices for research compendiums. However, feel free to organise and store your project files in a way that best suits your needs and workflow.
