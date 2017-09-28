
This repository contains all the data and scripts necessary to replicate the statistical analysis in in **"Big Data" for quantifying recreational ecosystem services and characterising people’s connection to nature**.
For a description of the data and the methods see the manuscript (pre-print coming soon).

The script *DistanceMarineSAC.R* contains code to calculate the distance from each cell centre to the boundaries of the nearest marine Special Area of Conservation.

The file *GooglePlaceAPI.R* contains the code necessary to download infrastructure data from Google Places API.

The folder *data* contains all the datasets used in the analysis; the file *CombineData.R* contains the script to combine them. 

To replicate the mixture model to classify observations according to GBIF records, run *MixtureModels.R*.

All the code necessary to replicate the statistical analysis is in the file *Analysis.R*, while the file *AnalysisReport.Rmd* contains a report of the analysis.
The files *Collinearity.R*, *Multiplot.R* and *nseq.R* contain functions used in the script *Analysis.R*.

*Env_ModSel_bash* contains a bash script used to run model selection on the University of Aberdeen HPC Maxwell, while the results of model averaging on the environmetal model are stored in the file *Env_avg_results.txt*.
