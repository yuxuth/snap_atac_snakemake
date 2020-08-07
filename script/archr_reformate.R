library(ArchR)
set.seed(1)
addArchRThreads(threads = 5) 
addArchRGenome("Mm10")

frag = snakemake@input[[1]]

reformatFragmentFiles(frag)
