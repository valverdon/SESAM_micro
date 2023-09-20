Modelling pipeline from raw to figures 15/02/2023

1) 03_Prep_data.R (Use metagenomic data + Soil data + ENVstack(itself using olivier.B + antoine A. layers)
--> Creates ENVdata, TotSeqSum, OTUdata and seqvec files for each group.


						- On cluster -  

2) 01_Fit.R (needs ENVdata, TotSeqSum and OTUdata + R codes for functions)
Fits and Evaluates models
--> Stores Models, Fit_met, Eval_met, Fit_data, Pred_data, VarSel, UnivExpl, Deviances files

2.1) Gather_VarImp_RFGBM.R (need Varsel)
--> Stores Var imp for Rf and GBM
						- Off Cluster -


3) 04_Gather.R  (needs Eval_met Fit_met folders in output folder)
--> Creates  Eval_met and Fit_mat files in data folder.


				- PA models only -

4) explo_nullmodels.R (needs Eval_met and Fit_mat files)
--> Creates Eval file in data folder

				- end PA only -



5) Table_numberofmodels. (needs Eval file)
--> Creates Fig 1 & 2 about number of models

6) Phylum_goodmodels_better.R
-->Creates fig 3 & 4 about Phlums

(following needs VarSel in outputs folder)
7)  VarImp_tests.R
--> Creates Variable_importance, Variable_ranks and Variable_preselected

8) VarImp_figures.R
--> Creates Varimp figures