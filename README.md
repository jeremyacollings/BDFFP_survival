The code in this repository was used to generate the analyses reported in "Climate change aggravates bird mortality in pristine tropical forests" (Wolfe et al.)

Initial analyses were based off of Wolfe's preliminary analyses in Program MARK. Collings then modified these analyses using Bayesian hierarchical modeling, following Yachulic et al. 2020. 
The "run_mod" scripts were used to fit models and reference the Stan scripts found in the "stan_files" directory. The "CJS_Results" scripts use the RDS files generated in the "run_mod" scripts to summarize the results of the models and generate figures. All estimates generated from the climate models and null model can be found in the "estimates" csvs. 
