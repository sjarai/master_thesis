# Correcting for Linkage Errors in Contingency Tables: New Methods to Improve the Correction Approach
## Master Thesis by Sjarai Eikenhout | Statistics and Data Science

### Introduction
Record linkage is a solution to the problem of recognizing records in two files that represent identical persons, objects, or events and aims to bring those records together. During the record linkage process, linkage errors can occur. These linkage errors are a particular
type of measurement error, which can lead to biased inference if no appropriate steps are taken to control and/or adjust this bias.

National Statistical Institutes aim to publish high-quality and accurate descriptive statistics of the population, which are used to inform policies and stakeholders. Nowadays, the use of administrative and new forms of data in statistical systems increases and linking data from multiple sources is also becoming more useful. More research is needed to understand the impact of linkage error on categorical data and to continue developing record linkage approaches.

In this project, two existing methods from Scholtus et al. (2021) and three new methods to correct for linkage errors are tested by means of a simulation study. The quality of the new estimators is measured by computing the errors, bias, variance, and
mean square error. By comparing the performances of the different correction methods on simulated data, it is ascertained which correction method is best for a given situation. This repository contains all R-code that was used in this project.

### Simulation study
In section 4.1 in the thesis, the design of the simulation study is described. To reproduce the simulation study, two files are needed, namely:  and . The first file, [Functions.R](Functions.R), contains all the functions that are needed in the simulation study. The second file [Simulation_study.R](Simulation_study.R) is a parallelised code that performs the simulation study. It also loads the functions from the first file. It is therefore recommended to save both R-files in the same directory and set the working directory in the code to this  directory. The simulation study takes approximately three days to run.

### Results
In section 4.2 in the thesis, the results of the simulation study are described. The figures and tables that were presented in this section can be reproduced by using the file [Results.R](Results.R). To be able to reproduce the results, the simulation study has to be performed in the same R-session, or the saved results from the simulation study from a previous R-session have to be loaded.

### Contact

For questions, comments or further information, please contact me: 

Sjarai Eikenhout - sjaraieikenhout@gmail.com
