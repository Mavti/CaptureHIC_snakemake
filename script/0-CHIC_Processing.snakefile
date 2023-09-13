# Snakemake pipeline
# CHIC processing pipeline
# Author Martina Rimoldi

import pandas as pd

wd = "/hps/software/users/flicek/rimoldi/Scripts/3DGenomeEvo/SM_PCE/"
rd1 = "/nfs/research/flicek/user/rimoldi/Projects/PCHIC/Analyses/Processing/"

chic=pd.read_csv(wd + "config/CHIC_samples_selected.txt",header=0,sep="\t")
species=chic["Species"].unique()
species

do=chic["do"].unique()
do

tissue=chic["Tissue"].unique()
tissue

rep=chic["replicate"].unique()
rep

sample=chic["do"]+"_PCHIC_"+chic["Species"]+"_"+chic["Tissue"]+"_"+chic["replicate"]
sample

cond = chic["Species"]+"_"+chic["Tissue"]
cond

#####################################################################################

rule all:
	input:
		expand("/nfs/research/flicek/user/rimoldi/Projects/PCHIC/Analyses/Processing/peakMatrix_f/{species}_toPeakMatrix.tab",  species=species),
		expand("/nfs/research/flicek/user/rimoldi/Projects/PCHIC/Analyses/Processing/peakMatrix_f/{species}_peakMatrix_avg_cutoff3.Rds", species=species)
	


rule create_matrixInput:
	output:
		otp_fl = "/nfs/research/flicek/user/rimoldi/Projects/PCHIC/Analyses/Processing/peakMatrix_f/{species}_toPeakMatrix.tab"
	params:
		sample_list =  expand("{sample}", sample=sample),
		otpDir = rd1,
		sp = "{species}"
	script:
		'''mktable_peakMatrix.py'''

#################################################################

rule peakMatrix:
	input:
		script = "/hps/software/users/flicek/rimoldi/chicagoTeam-chicago-e5a7d9247f7a/chicagoTools/makePeakMatrix.R",
		tab = "/nfs/research/flicek/user/rimoldi/Projects/PCHIC/Analyses/Processing/peakMatrix_f/{species}_toPeakMatrix.tab"
	output:
		"/nfs/research/flicek/user/rimoldi/Projects/PCHIC/Analyses/Processing/peakMatrix_f/{species}_peakMatrix_avg_cutoff3.Rds"
	params:
		fname = "/nfs/research/flicek/user/rimoldi/Projects/PCHIC/Analyses/Processing/peakMatrix_f/{species}_peakMatrix_avg_cutoff3",
		options = "--cutoff 5"
	conda:
		"R-4.0"
	shell:
		" Rscript {input.script} \
 				 {input.tab} \
 				{params.options} \
 				{params.fname}"

