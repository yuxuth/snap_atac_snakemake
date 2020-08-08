library(ArchR)
set.seed(1)
addArchRThreads(threads = 11) 
addArchRGenome("Mm10")

inputFiles = snakemake@input[[1]]
outdir = snakemake@params[[1]]
sample = snakemake@wildcards[['sample']]

if(!dir.exists(outdir)) dir.create(outdir)
setwd(outdir)
new_file = paste0('../../',inputFiles)

ArrowFiles <- createArrowFiles(
    inputFiles = new_file,
    sampleNames = sample,
    filterTSS = 4, #Dont set this too high because you can always increase later
    filterFrags = 1000, 
    addTileMat = TRUE,
    addGeneScoreMat = T
)

