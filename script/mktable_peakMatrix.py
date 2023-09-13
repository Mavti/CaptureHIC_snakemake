#! /usr/bin/env python

def make_table(sample_list, otp_Dir, species):
    samples_of_species = list(filter(lambda x: species in x, sample_list))
    flname = otp_Dir + "peakMatrix_f/" + species + "_toPeakMatrix.tab"

    for samp in samples_of_species:
        sp = samp.split("_")[2]
        tiss = samp.split("_")[3]
        combo = sp + "_" + tiss
        row = samp + "\t" + otp_Dir + combo + "/" + samp + ".rds"
       
        with open(flname, "a") as otp:
            otp.write(row + "\n")
    return

make_table(snakemake.params["sample_list"], snakemake.params["otpDir"], snakemake.params["sp"] )


