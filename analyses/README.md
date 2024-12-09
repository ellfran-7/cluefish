# 📄Content

This directory contains three original R scripts designed to generate derived data outside of the Cluefish workflow:

-  `dromics-transcriptomics-pipeline.R`: A step-by-step `DRomics` pipeline for processing the dose-response transcriptomic raw count dataset published in Franklin et al. (submitted). This script generates the input data for the Cluefish workflow in our case study.
-  `dromics-acnhoring-pipeline.R`: A step-by-step `DRomics` pipeline for analyzing dose-response anchoring (individual) measurements, including body length and eye surface data, as described in Franklin et al. (submitted). The results provide experimental validation for the Cluefish workflow results.
-  `standard-approach-pipeline.R`: A standard functional enrichment pipeline to generate results from the functional enrichment analysis of the deregulated gene list derived from the DRomics pipeline. These results serve as a benchmark for comparison with the Cluefish workflow outputs.

These pipelines are formatted to enable the direct reproduction of the results required for the Cluefish workflow and those published in Franklin et al. (submitted). However, they can also be used as adaptable templates for your own case studies. 

For help with the `DRomics` R package, visit the [DRomics documentation](https://lbbe-software.github.io/DRomics/).

## 📍Note

The suggested content and structure here follow best practices for research compendiums. However, feel free to organize and store your project files in a way that best suits your needs and workflow.
