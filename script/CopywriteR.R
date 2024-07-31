#install.packages("Rmpi")
#library(Rmpi)
library(CopywriteR)
library(CGHcall)

if(exists('snakemake')) {
    control = snakemake@input[["control"]]
    tumor = snakemake@input[["tumor"]]
    reference = snakemake@input[['reference']]
    #output = snakemake@output[[1]]
    output = snakemake@params[["output"]]

    manifest = snakemake@params[["manifest"]]
    binsize = snakemake@params[["binsize"]]
    thread = snakemake@threads
} else{
    stop("currently only support Snakemake")
}

print('control')
print(control)
print('tumor')
print(tumor)

## Set number of CPU (max usefull == number of samples)
bp.param <- SnowParam(workers = thread, type = "SOCK")

samples <- c(control, tumor)
controls <- samples[rep(1, length(samples))]
# controls <- samples
sample.control <- data.frame(samples, controls)

# RUN copywriter
# Set reference.folder to folder generated with preCopywriteR!
# create destination folder before running
if (!dir.exists(output)) {
   dir.create(output, recursive = T)
} else{
	system(paste0("mv ", output, " ", output, "_backup"))
	dir.create(output, recursive=T)
}

CopywriteR(sample.control = sample.control, 
    destination.folder = output,
    reference.folder = reference,
    capture.regions.file = manifest,
    keep.intermediary.files=TRUE,
    bp.param = bp.param)


## Read output file from CopywriteR
data <- read.table(file = paste0(output, "/CNAprofiles/log2_read_counts.igv"), sep="\t", header=T, quote="", fill=T, skip=1, as.is=T)

## Run CGHcall for segmentation and calling of tha data (Wiel et al. BioInformatics 2007)
Corrected_logratio<-data[,c(4,1,2,3,5:ncol(data))]

raw <- make_cghRaw(Corrected_logratio)
prep <- preprocess(raw, maxmiss = 0, nchrom = 22) # TODO: from here, also try my QDNAseq pipeline; compare

saveRDS(prep, paste0(output, '/preprocessed_Copywriter.rds'))


nor <-  normalize(prep, method = "median", smoothOutliers = TRUE)  

seg <-  segmentData(nor, method = "DNAcopy",nperm=2000,undo.splits="sdundo", min.width=5,undo.SD=1, clen=25, relSDlong=5)
segnorm <- postsegnormalize(seg,inter=c(-0.4,0.4))

#save(segnorm, file="CopywriteR_100kbp_seg.Rdata")  
saveRDS(segnorm, file=paste0(output, "CopywriteR_seg.rds"))  


## plot outcome
plot.profiles <- function(cgh, directory, byChr=FALSE) {
  tmp <- sampleNames(cgh)
  if (!file.exists(directory))
    dir.create(directory, recursive = T)
  if ('filter' %in% colnames(fData(cgh))) {
    chrs <- unique(chromosomes(cgh)[fData(cgh)$filter])
  } else {
    chrs <- unique(chromosomes(cgh))
  }
  for (i in 1:length(sampleNames(cgh))) {
    if (byChr) {
      for (chr in chrs) {
        png(file.path(directory, paste(tmp[i], '-chr', chr, '.png', sep='')), width=297, height=210, units='mm', res=150)
	print(sum(chromosomes(cgh) == chr))
	Sys.sleep(1)
        plot(cgh[chromosomes(cgh) == chr,i], ylab=expression(normalized~log[2]~read~count), dotres=1)
        dev.off()
      }
    } else {
      png(file.path(directory, paste(tmp[i], '.png', sep='')), width=297, height=210, units='mm', res=150)
      plot(cgh[,i], ylab=expression(normalized~log[2]~read~count), dotres=1)
      dev.off()
    }
  }
}


plot.profiles(segnorm, output)
