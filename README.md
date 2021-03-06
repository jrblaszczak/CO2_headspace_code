# Hall lab headspace equilibration correction comparison

## Main folder
1. "Differences between Kos and Hall.Rmd" - Rmarkdown file of compare major differences between our group (Hall lab) and that of Koschorreck et al. for calculating the CO2 concentration in a water sample after headspace equilibration as a part of our comment to their manuscript in review with Biogeosciences. 
2. "Differences-between-Kos-and-Hall.html" - Knit version of Rmarkdown file

## Code folder
1. "2020_10 Data Compile.R" - Code to compile raw data files from the Data folder into appropriate format for code
2. "2020_10 HS correction data test.R" - Initial test of Koschorreck code using formatted Blaszczak data from MT and AZ
3. "Comparing approaches to solve for H.R" - Initial code comparing approaches to solving for H<sup>+</sup>
4. "Hall lab CO2 headspace equilibration.R" - Hall lab code for solving carbonate chemistry
5. "Koschorreck et al. hs_equil_correction.R" - code from https://github.com/icra/headspace included in Koschorreck et al.
6. "Koschorreck_Rheadspace.R" - code from https://github.com/icra/headspace included in Koschorreck et al. with some lines commented out so it can be sourced in "2020_10 HS correction data test.R"
7. Oct2020_CO2HeadspaceEquations.html" - Knit version of Rmarkdown file with equations to solve for carbonate equilibrium
8. Oct2020_CO2HeadspaceEquations.Rmd" - Rmarkdown file with derivation of equations for solving carbonate chemistry
9. "picarro_data_harvest_withQAQC.R" - Code that processes raw Picarro output from gas sample injection


## Data folder
Data files from Montana and Arizona streams including among raw data files:
1. "~Diel_QAQC.csv" - processed Picarro gas analyzer output from headspace equilibrated gas sample injections
2. "2020_10_02_HScorrection_allformatted.csv" - compiled and formatted data for headspace equilibration correction comparison
