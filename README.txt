Supplementary information / reproducible research files for the manuscript 
Title: "A general approach to fitting multistate cure models based on an extended-long-format data structure"

Authors: Yilin Jiang, Harm van Tinteren and Marta Fiocco
Code was written by Yilin Jiang
In case of questions or comments please contact y.jiang@math.leidenuniv.nl.

The code was written/evaluated in R with the following software versions:
R version 4.6.1 (2026-06-24 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=English_United States.utf8 
[2] LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8
[4] LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

time zone: Europe/Amsterdam
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] latex2exp_0.9.8      icmstate_0.2.0       flexsurv_2.3.2      
 [4] mice_3.19.0          MultiCure_0.0.0.9000 GenKern_1.2-60      
 [7] KernSmooth_2.23-26   cowplot_1.2.0        data.table_1.18.4   
[10] dplyr_1.2.1          ggplot2_4.0.3        mstate_0.3.3        
[13] survival_3.8-6      

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3  rstudioapi_0.19.0   shape_1.4.6.1      
  [4] magrittr_2.0.5      rainbow_3.8         jomo_2.7-6         
  [7] farver_2.1.2        nloptr_2.2.1        fs_2.1.0           
 [10] fields_17.3         vctrs_0.7.3         memoise_2.0.1      
 [13] minqa_1.2.8         RCurl_1.98-1.19     JOPS_0.2.0         
 [16] rstatix_0.7.3       usethis_3.2.1       broom_1.0.13       
 [19] deSolve_1.42        Formula_1.2-5       mitml_0.4-5        
 [22] hdrcde_3.5.0        parallelly_1.47.0   pracma_2.4.6       
 [25] cachem_1.1.0        igraph_2.3.3        lifecycle_1.0.5    
 [28] cmprsk_2.2-12       iterators_1.0.14    pkgconfig_2.0.3    
 [31] Matrix_1.7-5        R6_2.6.1            fastmap_1.2.0      
 [34] rbibutils_2.4.1     future_1.70.0       digest_0.6.39      
 [37] numDeriv_2016.8-1.1 colorspace_2.1-2    rprojroot_2.1.1    
 [40] pkgload_1.5.3       ggpubr_0.6.3        labeling_0.4.3     
 [43] abind_1.4-8         compiler_4.6.1      remotes_2.5.0      
 [46] withr_3.0.3         S7_0.2.2            backports_1.5.1    
 [49] carData_3.0-6       pkgbuild_1.4.8      maps_3.4.3         
 [52] ggsignif_0.6.4      pan_1.9             MASS_7.3-65        
 [55] lava_1.9.1          sessioninfo_1.2.4   tools_4.6.1        
 [58] future.apply_1.20.2 nnet_7.3-20         glue_1.8.1         
 [61] quadprog_1.5-8      nlme_3.1-169        grid_4.6.1         
 [64] checkmate_2.3.4     cluster_2.1.8.2     generics_0.1.4     
 [67] gtable_0.3.6        SpATS_1.0-20        tidyr_1.3.2        
 [70] survminer_0.5.2     car_3.1-5           foreach_1.5.2      
 [73] pillar_1.11.1       spam_2.11-4         splines_4.6.1      
 [76] lattice_0.22-9      ks_1.15.2           tidyselect_1.2.1   
 [79] fds_1.9             muhaz_1.2.6.4       reformulas_0.4.4   
 [82] gridExtra_2.3       expm_1.0-0          statmod_1.5.2      
 [85] devtools_2.5.2      boot_1.3-32         codetools_0.2-20   
 [88] msm_1.8.2           tibble_3.3.1        cli_3.6.6          
 [91] rpart_4.1.27        Rdpack_2.6.6        Rcpp_1.1.1-1.1     
 [94] globals_0.19.1      parallel_4.6.1      ellipsis_0.3.3     
 [97] dotCall64_1.2       mclust_6.1.2        bitops_1.0-9       
[100] lme4_2.0-1          listenv_0.10.1      glmnet_5.0         
[103] viridisLite_0.4.3   mvtnorm_1.4-1       scales_1.4.0       
[106] prodlim_2026.03.11  pcaPP_2.0-5         purrr_1.2.2        
[109] rlang_1.2.0        

This folder contains the following files that can be used to reproduce all the analyses, result figures and tables of the manuscript.
It contains four subfolders containing the following files:

./function/:
This folder contains all the general functions required for running our proposed EM algorithm based on the extended long format (ELF).

./ebmt/:
This folder contains the code files and results for the EBMT application example.
	Core EM algorithm.R
	An R script that contains the code for the main EBMT analysis. The application data comes from the R package mstate. Patient 	baseline characteristics table (table 1), point estimates for the regression coefficients (table 2), cumulative hazards by 	transition plot (figure A.1), the relevant tables in Appendix B, Appendix F can be yielded from this code.

	Bootstrap.R
	An R script that contains the code for the bootstrap procedure. Standard errors for the regression coefficients are calculated 	using this code (table 2).

	Prediction.R
	An R script that contains the code for prediction. Figure 6 and 7 are generated from this code.

	./ebmt results/
	A subfolder that contains the results for the EBMT application. The main results are stored in result_elf.rda. cox_elf.csv and 	logit_elf.csv contain the tables of the extended long data format from the last iteration of the proposed EM algorithm. 	Bootstrap_results.csv contains the coefficient estimates from all the bootstrap samples, while bootstrap_failures.csv records the 	2 bootstrap samples with numerical issues and the reasons reported from R. Prob of being cured after transplant.pdf, dynamic 	prediction day0.pdf, dynamic prediction day100.pdf and dynamic prediction day310.pdf are the illustrations for dynamic prediction 	in section 8(figure 6 and 7). 

./simulation 1/:
This folder contains the code files and results for simulation 1 in Section 7.1 and Appendix D.
	sim1_zerotail.R
	An R script that contains the code for running ELF method and MultiCure method under zero-tail constraint in the setting as described in Section 7.1.1. Data used for simulation 1 is generated from the MultiCure package (Taylor and Beesley (2019)).

	sim1_nozerotail.R
	An R script that contains the code for running ELF method and MultiCure method without the zero-tail constraint in the setting as described in Section 7.1.1. Data used for simulation 1 is generated from the MultiCure package (Taylor and Beesley (2019)).

	process_results.R
	An R script that contains the code for calculating the performance measures, drawing the comparison plots per performance measures (including figure 2 and 3) and generating computation time summaries in simulation 1, along with the specific functions needed for result processing.

	./results_zerotail/
	A subfolder that contains the results for simulation 1 under zero-tail constraint. results_zerotail.rds stores the combined results of the ELF and MultiCure methods across the 500 repetitions under zero-tail constraint. performance_results_zerotail.csv contains the calculated performance measures along with their MC errors per parameter per method. runtime_summary_zerotail.csv contains the running time summary for both methods. nested_bias_zerotail.pdf (figure 2), nested_emp_se_zerotail.pdf and nested_rmse_zerotail.pdf compares the performances of both methods under zero-tail constraint.


	./results_nozerotail/
	A subfolder that contains the results for simulation 1 without zero-tail constraint. The structure is identical to that of the results_zerotail folder. results_nozerotail.rds stores the combined results of the ELF and MultiCure methods across the 500 repetitions without zero-tail constraint. performance_results_nozerotail.csv contains the calculated performance measures along with their MC errors per parameter per method. runtime_summary_nozerotail.csv contains the running time summary for both methods. nested_bias_nozerotail.pdf (figure 3), nested_emp_se_nozerotail.pdf and nested_rmse_nozerotail.pdf compares the performances of both methods without zero-tail constraint.

./simulation 2/:
This folder contains the code files and results for simulation 2 in Section 7.2 and Appendix E.
	generate_data.R
	An R script that contains the code to generate simulation data. It stores several functions specific for simulation 2, which are used in Make scenarios.R and run_scenario.R.

	Make scenarios.R
	An R script that contains the code to make the full factorial design for simulation 2. To be specific, it finds the values of alpha_0 at predefined cure proportions and the values of the scale parameter for the censoring distribution at predefined censoring rates.

	scenarios.rda
	This file contains the full factorial design for simulation 2, generated from Make scenarios.R.

	run_scenario.R
	An R script that contains a function to run ELF and MultiCure methods for a pre-specified scenario in Simulation 2. The function is called in run_simulation.R.

	run_simulation.R
	An R script containing the code to execute the simulation across the scenarios specified in Simulation 2. 

	pool simulation results.R
	An R script that contains the code for calculating the performance measures for ELF and MultiCure methods per scenario and generating computation time summaries in simulation 2, along with the specific functions needed for result processing.

	nested loop across scenarios.R
	An R script that contains the code for generating nested loop plots across scenarios per performance metric per parameter (including figure 4 and 5).

	./simulation 2 results/
	A subfolder that contains all Simulation 2 results. For each scenario, you can find the combined results of the ELF and MultiCure methods across the 500 repetitions per scenario (named as results_scenario_x.rds), the performance measures along with their MC errors broken down by parameter and method (named as performance_x.csv) and the nested loop plots comparing performance across scenarios for each metric and parameter (named as plot_<parameter>_<metric>.pdf) and a summary of computation time for both methods (runtime_summary.csv).



 


