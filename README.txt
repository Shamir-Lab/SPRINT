This code was devloped using R version 3.4.1 and runs exclusively on the datasets from the SPRINT trial (NEJM 2015)
This data can be requested through BioLINCC: https://biolincc.nhlbi.nih.gov/studies/sprint_pop/
Needed libraries and their versions (updated versions probebly will work as well): ggplot2_2.2.1,sas7bdat_0.5,survival_2.41-3 impute_1.50.1,PRROC_1.3
This code assumes all datasets and scripts are in the current working directory (i.e. under "./")
in order to reproduce the results correctly please execute the code as follows:
1) run the source() function for the following scripts: 
generate_sprint_raw_data.R
fix_framingham_score.R
generate_accord_raw_data.R
extract_dynamics_SPRINT.R
extract_dynamics_ACCORD.R
predictors.R
run_predictor.R
treatment_reccomendation.R
(or run the first lines of "analyze_sprint.R")
2) reproduce the results on the SPRINT database (including the treatment assignment simulation) by running "analyze_sprint.R" top to bottom
3) reproduce the results on the ACCORD database by running "analyze_accord.R" top to bottom
