library(ArchR)
set.seed(1)
addArchRThreads(threads = 8) 
addArchRGenome("Mm10")

inputFiles = snakemake@input[[1]]
outdir = snakemake@output[[1]]
sample = snakemake@wildcards[['sample']]


ArrowFiles <- createArrowFiles(
    inputFiles = inputFiles,
    sampleNames = sample,
    filterTSS = 4, #Dont set this too high because you can always increase later
    filterFrags = 1000, 
    addTileMat = TRUE,
    addGeneScoreMat = T
)

# proj <- ArchRProject(
#     ArrowFiles = ArrowFiles, 
#     outputDirectory = outdir,
#     copyArrows = F #This is recommened so that you maintain an unaltered copy for later usage.
# )
