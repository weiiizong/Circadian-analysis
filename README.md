## Circadian Analysis: Rhythmicity detection with confidence intervals 
The purpose of this software is to demonstrate rhythmicity detection of circadian patterns in gene expression profiles. Full analysis using this code can be found in the paper: Diurnal rhythms across the human dorsal and ventral striatum (Kyle D. Ketchesin, under review)

# Required R packages: 
- minpack.lm v1.2-1
- Snowfall v1.84-6.1
- doParallel v1.014
- AWfisher https://github.com/Caleb-Huo/AWFisher
*This code has been tested on R version 4.0.2 

# Installation Instructions: 
Download and unzip the folders to the your working directory. 
Save each folder to the desktop
R/: code for analysis demostration
data/: example data used in code
ouput/: results generated from code in R/
Install required packages listed above

# Clinical Data (Example_clinical.csv):
Clinical variables of interest include:
- CorrectedTOD: corrected time of death (TOD) in ZT. 
- pair: sample names matched to data columns.

# Expression data (Example_data.csv)::  
- 15330 genes (all genes that were used in the full analysis of NAc region in the paper)
- Data has already been filtered according to method discussed in the paper 
- Expression units are in log2 CPM 

# Fitted estimates data (Example_observed_NAc.csv and Example_observed_Caudate.csv):
- Observed fitted estimates in NAc and Caudate regions used in paper for phase concordance plots demostration.

# Demo:
- fitSinCurve.R and Curve_Drawing_axis.R are external functions required to run the main analysis in Rhythmicity.R 
- Run all external functions, so that they are in the R environment 
- Run code as is in Rhythmicity.R to derive fitted parameters
- Run code as is in Rhythmicity_peakCI.R to derive confidence interval of peak
- Run code as is in PhaseShift.R to derive phase concordance plot
*NOTE: Permutations and Bootstraps were set to 10 to save time and computation space. Permuted and bootstrap files are included in the output/Rhythmicity/null_control and output/Rhythmicity/bootstrap_control respectively to save reviewers time. If files are recreated by the reviewer, output will be slightly different due to the randomness. 

# Expected output:
- output/Rhythmicity/observed_para_demo.csv: the example output of curve fitting parameters and p values. 
- output/Rhythmicity/TopPlots: PDF plots of top circadian genes 
- output/Rhythmicity/CorePlots: PDF plots of core circadian genes
- output/Rhythmicity/observed_para_peakCI_demo.csv: the example output of curve fitting parameters, p values and confidence interval of peak estimates. 
- output/PhaseShift: PDF plots of phase concordance plots 
