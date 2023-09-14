# Snakemake pipeline
# CHIC processing pipeline
# Author Martina Rimoldi

import pandas as pd
import glob
configfile: "config/config.yaml"

wd = config["wd"] # snakemake directory
rd1 = config["rd1"] # results directory
rawdata = config["rawdata"]
probesDir = config["probesDir"]
chicagoTools = config["chicagoTools"]

chic=pd.read_csv(wd + "config/CHIC_samples_selected.txt",header=0,sep="\t")
species=sorted(chic["Species"].unique())
species

do=chic["do"].unique()
do

tissue=chic["Tissue"].unique()
tissue

rep=chic["replicate"].unique()
rep

sample=chic["do"]+"_PCHIC_"+chic["Species"]+"_"+chic["Tissue"]+"_"+chic["replicate"]
sample
# sort according to species
sample= sorted(sample, key=lambda x: x.split("_")[2])

cond = chic["Species"]+"_"+chic["Tissue"]
cond


### GET a species list with as many species entry as the number of samples
## for that species
species_rep=[]
for samp in sample:
	sp=samp.split("_")[2]
	species_rep.append(sp)

# file extensions of the chicago pipeline outputs
chicago_ext=[".rds", "_table.txt", ".ibed", "_seqmonk.txt", "_washU_text.txt"]

#print(expand( "{SPECIES_rep}_DesignDir/{SAMPLE}.chicinput", zip, SPECIES_rep=species_rep, SAMPLE=sample))
#####################################################################################

rule all:
	input:
		expand( rawdata + "Hicup/{SAMPLE}_r_1_2.hicup.bam", SAMPLE=sample),
		expand( rd1 + "{SPECIES}_DesignDir/{SPECIES}.rmap", SPECIES=species ),
		expand( rd1 + "{SPECIES}_DesignDir/{SPECIES}.baitmap", SPECIES=species ),
		expand( rd1 + "{SPECIES}_DesignDir/{SPECIES}.baitmap_4col.txt", SPECIES=species ),
		expand( rd1 + "{SPECIES}_DesignDir/{SPECIES}.{ext}", SPECIES=species, ext=["poe","npb","nbpb"]),
		expand( rd1 + "{SPECIES_rep}_DesignDir/{SAMPLE}/{SAMPLE}.chicinput", zip, SPECIES_rep=species_rep, SAMPLE=sample ), # species_rep is a list containg as many repeated species as the number of samples for that species I have
		expand( rd1 + "{SPECIES_rep}_DesignDir/{SAMPLE}/{SAMPLE}_bait2bait.bedpe", zip, SPECIES_rep=species_rep, SAMPLE=sample ),
		expand( rd1 + "{SPECIES_rep}_DesignDir/{SAMPLE}.chicinput", zip, SPECIES_rep=species_rep, SAMPLE=sample ),
		expand( rd1 + "{SPECIES_rep}_DesignDir/{SAMPLE}_bait2bait.bedpe", zip, SPECIES_rep=species_rep, SAMPLE=sample ),
		expand( rd1 + "peakMatrix_f/{SPECIES}_toPeakMatrix.tab",  SPECIES=species),
		expand( rd1 + "peakMatrix_f/{SPECIES}_peakMatrix_avg_cutoff5.Rds", SPECIES=species),
		expand( rd1 + "{SPECIES}_{TISSUE}/{SPECIES}_{TISSUE}{chicago_ext}", SPECIES=species, TISSUE=tissue, chicago_ext=chicago_ext)

# Getting digest file for each genome
def get_digest(wildcards):
	# code that returns a list of fastq files for read for one *do* library
	return glob.glob(config["digest_path"] + "Digest_" + wildcards.SPECIES + "_HindIII_*txt")

### RUN HICUP:
# Takes care of: QC, alignment, filtering etc.. 
rule HICUP:
	input:
		r1 = rawdata + "Raw/{SAMPLE}_r_1.fq.gz",
		r2 = rawdata + "Raw/{SAMPLE}_r_2.fq.gz",
		index = expand( config["genome_path"] +  "{SPECIES}/*_filtered,fa", SPECIES=species)
		
	output:
		rawdata + "Hicup/{SAMPLE}_r_1_2.hicup.bam"
	params:
		digest = get_digest,
		bowtie2 = config["bowtie2"],
		otpDir = rawdata + "Hicup/",
	threads: 16 
	conda:
		"PCHIC"
	shell:
		'''hicup  \
		--digest {params.digest} \
		--format Sanger \
		--index {input.index} \
		--longest 700 \ 
		--zip 1 \
		--keep 1 \
		--bowtie2 {params.bowtie2} \
		--outdir  {params.otpDir} \
		--shortest 50 \
		--threads {threads} \
		{input.r1} {input.r2} '''


#######################################################
################# CHICAGO DESIGN FILES ###############
rule Chicago_RMAP:
	input:
		digest = get_digest
	output:
		rd1 + "{SPECIES}_DesignDir/{SPECIES}.rmap" 
	params:
		DesDir = rd1 + "{SPECIES}_DesignDir/"
	conda:
		"PCHIC"
	shell:
		'''
		mkdir -p {params.DesDir}
		tail -n +3 {input.digest} |
            awk 'BEGIN{OFS=\"\t\";n=0};{n+=1; print \$1, \$2, \$3, n }' > {output}
		'''

rule Chicago_BAITMAP:
	input:
		rmap = rd1 + "{SPECIES}_DesignDir/{SPECIES}.rmap",
		probes = probesDir + "{SPECIES}_filt_finalTargets.bed"
	output:
		rd1 + "{SPECIES}_DesignDir/{SPECIES}.baitmap"
	conda:
		"PCHIC"
	shell:
		'''
		bedtools intersect -a {input.rmap}  -b {input.probes} |
            bedtools groupby -i stdin -g 1-4 -c 9 -o distinct > {output}
		'''

rule Chicago_BAITMAP_4col:
	input:
		rd1 + "{SPECIES}_DesignDir/{SPECIES}.baitmap"
	output:
		rd1 + "{SPECIES}_DesignDir/{SPECIES}.baitmap_4col.txt"
	conda:
		"PCHIC"
	shell:
		'''
		cat {input} | cut -f1-4  > ${output}
		'''

rule Chicago_DesignFiles:
	input:
		rmap = rd1 + "{SPECIES}_DesignDir/{SPECIES}.rmap",
		baitmap = rd1 + "{SPECIES}_DesignDir/{SPECIES}.baitmap"
	output:
		rd1 + "{SPECIES}_DesignDir/{SPECIES}{.poe,.npb,.nbpb}"
	params:
		chicago_script = chicagoTools + "makeDesignFiles.py",
		outPrefix = rd1 + "{SPECIES}_DesignDir/{SPECIES}",
		desDir = rd1 + "{SPECIES}_DesignDir"
	conda:
		"PCHIC"
	shell:
		'''
		python2.7 {params.chicago_script} \
                    --rmapfile = {input.rmap}   \
                    --baitmapfile = {input.baitmap} \
                    --outfilePrefix= {params.outPrefix} \
                    --minFragLen=150 \
                    --maxFragLen=40000 \
                    --designDir= {params.desDir}
		'''

######################################################
# def decompose_sample(wildcards):
# 	sample_comp={}
# 	# SAMPLE = [ SPECIES, TISSUE, REP ]
# 	sample_comp[wildcards.SAMPLE]=[wildcards.SAMPLE.split("_")[2],wildcards.SAMPLE.split("_")[3], wildcards.SAMPLE.split("_")[4] ]
# 	return(sample_comp)



rule Chicago_inputs:
	input:
		bam =  rawdata + "Hicup/{SAMPLE}_r_1_2.hicup.bam",
		rmap = rd1 + "{SPECIES}_DesignDir/{SPECIES}.rmap",
		baitmap = rd1 + "{SPECIES}_DesignDir/{SPECIES}.baitmap"
	params:
		chicago_script = chicagoTools + "bam2chicago_V02.sh",
		outname = rd1 +"{SPECIES}_DesignDir/{SAMPLE}/"
		#get_sp = decompose_sample
	output:
		rd1 + "{SPECIES}_DesignDir/{SAMPLE}/{SAMPLE}.chicinput",
		rd1 + "{SPECIES}_DesignDir/{SAMPLE}/{SAMPLE}_bait2bait.bedpe"
	conda:
		"PCHIC"
	shell:
		'''
		{params.chicago_script} \
            --bamfile {input.bam} \
            --baitmap {input.baitmap} \
            --rmap {input.rmap} \
            --outname {params.outname}
		'''

#HOUSEKEEPING: move chicinputs into new folder
rule HK_moveInputs:
	input:
		rd1 + "{SPECIES}_DesignDir/{SAMPLE}/{SAMPLE}.chicinput",
		rd1 + "{SPECIES}_DesignDir/{SAMPLE}/{SAMPLE}_bait2bait.bedpe"
	output:
		rd1 + "{SPECIES}_DesignDir/{SAMPLE}.chicinput",
		rd1 + "{SPECIES}_DesignDir/{SAMPLE}_bait2bait.bedpe"
	params:
		destination = rd1 + "{SPECIES}_DesignDir/",
		origin = rd1 + "{SPECIES}_DesignDir/{SAMPLE}/"
	shell:
		'''
		mv {input} ../ {params.destination}
        rm -r rd1 + {params.origin}
		'''


#################################################################

rule create_matrixInput:
	output:
		otp_fl = rd1 + "peakMatrix_f/{SPECIES}_toPeakMatrix.tab"
	params:
		sample_list =  expand("{SAMPLE}", SAMPLE=sample),
		otpDir = rd1,
		sp = "{SPECIES}"
	script:
		'''mktable_peakMatrix.py'''


rule peakMatrix:
	input:
		script = chicagoTools + "/makePeakMatrix.R",
		tab = rd1 + "peakMatrix_f/{SPECIES}_toPeakMatrix.tab"
	output:
		rd1 + "peakMatrix_f/{SPECIES}_peakMatrix_avg_cutoff5.Rds"
	params:
		fname = rd1 + "peakMatrix_f/{SPECIES}_peakMatrix_avg_cutoff5",
		options = "--cutoff 5"
	conda:
		"R-4.0"
	shell:
		" Rscript {input.script} \
 				 {input.tab} \
 				{params.options} \
 				{params.fname}"

#################################################################
################ CHICAGO PIPELINE ###############################


rule CHICAGO:
	input:
		script = "script/CHICprocessing.R"
	output:
		rd1 + "{SPECIES}_{TISSUE}/{SPECIES}_{TISSUE}{chicago_ext}"
	params:
		mode= "standard",
		cond = "{SPECIES}_{TISSUE}",
		wdir = rd1
	conda:
		"R-4.0"
	shell:
		'''
		 Rscript {input.script} {params.mode} {params.wdir} {params.wdir} {params.cond}
		'''


# rule Chicago_BGsparsity:
# 	input:

# 	output:

# 	conda:
# 		"PCHIC"

