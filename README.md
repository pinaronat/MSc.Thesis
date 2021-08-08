# MSc.Thesis
Contains various code files I wrote during my thesis.

*Abeta_PAMclustering.R*
*asyn_PAMClustering.R*
*tau_PAMClustering.R*  -  Contains the separate R code files for PAM Clustering applied to the proteins that have a homologous sequence to the APRs of three causal proteins selected for this study (A-beta, tau and alpha-synuclein)

*aa_freq.py* - calculates the freuqencies of each aminoacid in the human proteome (FASTA file obtained from UniProt)

*take_snippets_repeat.py*  -  python code file which takes 500 three complexity, 500 four comlexity and 500 five complexity sequences from the human proteome randomly, for 100 repeats.

*separate_plots.R*  -  R code file that fits a second degree regression to predict the number of matches of randomly taken sippets in the human proteome.

*SeqComplexity_predict_error.R*  -  R code file that repeats the regression applied in *separate_plots.R* for 100 sets of randomly taken snippets
