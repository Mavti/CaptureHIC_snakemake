renv::load(project="/hps/software/users/flicek/rimoldi/Scripts/3DGenomeEvo/")

library("Chicago")
library("tidyverse")
args <- commandArgs(trailingOnly=TRUE)

phase=args[1]
wDir=args[2]
desDir=args[3]
bname=args[4]

sp=str_split(bname,"_")[[1]][1]
inputDir=paste0(wDir,"/",bname,"/selected_BR3-BR4")

# get input files
inputs=list.files(path = inputDir, full.names = TRUE, pattern="*chinput")

if( phase == "standard") {
    testDir=paste0(desDir, "/",sp,"_DesignDir")

    # create blank object with design files
    chicagoData <- setExperiment(designDir = testDir)
    # read in input data
    chicagoData <- readAndMerge(files=inputs, cd=chicagoData)

    prefixName=paste0(wDir,"/",bname, "/selected_BR3-BR4/", bname)
    chicagoData <- chicagoPipeline(chicagoData, prefixName)

    exportResults(chicagoData, file.path(inputDir, bname))

    rdsName=paste0(prefixName, ".rds")
    saveRDS(chicagoData, file=rdsName)

    pngDir=paste0(wDir,"/",bname,"/selected_BR3-BR4/baitExamples/", bname)
    for (i in c(1:6)) {
        pngName=paste0(pngDir,"_",i, ".pdf")
        pdf(pngName)
        plotBaits(chicagoData, n=1, plotBprof=T)
        dev.off()
    }


    flname=paste0(prefixName, "_table.txt")
    write.table(x=chicagoData@x, file=flname, sep="\t", quote=F, row.names=F)

} else if( phase == "toPeakMatrix" ){
    #sample=args[3]
    # here bname corresponds to sample (sample example: do26735_PCHIC_Dog_Brain_BR1)

    do=str_split(bname,"_")[[1]][1]
    sp=str_split(bname,"_")[[1]][3]
    tissue=str_split(bname,"_")[[1]][4]

    testDir=paste0(desDir, "/",sp,"_DesignDir")
    inputDir=paste0(wDir,"/",sp,"_",tissue)

    # get input file
    inputFl=paste0(inputDir,"/",bname, ".chinput")
    #input=list.files(path = inputDir, full.names = TRUE, pattern="*chinput")
    # create blank object with design files
    chicagoData <- setExperiment(designDir = testDir)
    # read in input data
    chicagoData <- readAndMerge(files=inputFl, cd=chicagoData)

    prefixName=paste0(inputDir, "/",bname)
    chicagoData <- chicagoPipeline(chicagoData, prefixName)

    rdsName=paste0(prefixName, ".rds")
    saveRDS(chicagoData, file=rdsName)


} else if( phase == "recalibrated" ){
    testDir=paste0(desDir, "/",sp,"_DesignDir")

    # create blank object with design files
    settingFl=paste0(inputDir,"/",bname,"_pvalweights.settings")
    chicagoData <- setExperiment(designDir = testDir,settingsFile=settingFl)
    # read in input data
    chicagoData <- readAndMerge(files=inputs, cd=chicagoData)

    prefixName=paste0(wDir,"/",bname, "_recab/", bname)
    chicagoData <- chicagoPipeline(chicagoData, prefixName)

    otpDir=paste0(wDir,"/",bname, "_recab")
    exportResults(chicagoData, file.path(otpDir, bname))

    rdsName=paste0(prefixName, ".rds")
    saveRDS(chicagoData, file=rdsName)

    pngDir=paste0(wDir,"/",bname,"/baitExamples/", bname, "_recab")
    for (i in c(1:6)) {
        pngName=paste0(pngDir,"_",i, ".png")
        png(pngName)
        plotBaits(chicagoData, n=1, plotBprof=T)
        dev.off()
    }

    flname=paste0(prefixName, "_table.txt")
    write.table(x=chicagoData@x, file=flname, sep="\t", quote=F, row.names=F)

} else if( phase == "testis" ){
    kb=args[4]
    testDir=paste0(desDir, "/",sp,"_DesignDir_", kb)

    # create blank object with design files
    chicagoData <- setExperiment(designDir = testDir)
    # read in input data
    chicagoData <- readAndMerge(files=inputs, cd=chicagoData)

    prefixName=paste0(wDir,"/",bname, "_",kb,"/", bname)
    chicagoData <- chicagoPipeline(chicagoData, prefixName)

    otpDir=paste0(wDir,"/",bname,  "_",kb)
    exportResults(chicagoData, file.path(otpDir, bname))

    rdsName=paste0(prefixName, ".rds")
    saveRDS(chicagoData, file=rdsName)

    pngDir=paste0(wDir,"/",bname, "_",kb,"/baitExamples/", bname, "_",kb)
    for (i in c(1:6)) {
        pngName=paste0(pngDir,"_",i, ".png")
        png(pngName)
        plotBaits(chicagoData, n=1, plotBprof=T)
        dev.off()
    }

    flname=paste0(prefixName, "_table.txt")
    write.table(x=chicagoData@x, file=flname, sep="\t", quote=F, row.names=F)
}
