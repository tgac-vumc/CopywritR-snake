
library("CopywriteR")

if(exists('snakemake')) {
    binsize = snakemake@params[['binsize']]
    refgenome = snakemake@params[['ref']]
    chrname = snakemake@params[['chrname']]
    path = snakemake@params[['path']]
} else{
    stop("Currently only support Snakemake")
}

dir.create(path)
preCopywriteR(output.folder = file.path(path),
                bin.size = binsize,
                ref.genome = refgenome,
                prefix = chrname)




