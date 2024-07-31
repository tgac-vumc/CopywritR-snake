import pandas as pd
import glob
import os.path as path
configfile: "./config.yaml"

## configurations ###########
PATH_BAM = config['path']['bam']
PATH_LOG = config['path']['log']
PATH_TEMP = config['path']['temp']
PATH_PIPELINE = srcdir('.')
#################################

PLATFORM = config['platform']['SRorPE'] # SE or PE
PREFIX = config['platform']['prefix'] # prefix before "bam"
#

##### Configurations VariantDetection ######
PATH_CNV = config['path']['copywriter']


## obtain sample list ###########
all_Samples=pd.read_csv(config['path']['sampleList'])
print(all_Samples)
Tumor=list(all_Samples['samples']) # these are IDs, and use as the key to fetch files.
print(Tumor)
Normal=list(all_Samples['controls'])
print(Normal)

def getnames(samplelist, platform, pathdata):
    Files = []
    RNAIDs = []
    SAMPLES = dict()
    for sample in samplelist:
        if platform in ['SR', 'sr']:
            for prefix in PREFIX:
                SAMPLES[sample] = glob.glob(path.join(pathdata, sample+'*'+prefix+'.bam'))
                SAMPLES[sample].sort()
        if platform in ['PE', 'pe']:
            SAMPLES[sample] = dict()
            SAMPLES[sample]['R1'] = glob.glob(path.join(pathdata, sample+'*'+PREFIX[0]+'.fastq.gz'))
            SAMPLES[sample]['R1'].sort()
            SAMPLES[sample]['R2'] = glob.glob(path.join(pathdata, sample+'*'+PREFIX[1]+'.fastq.gz'))
            SAMPLES[sample]['R2'].sort()
    return(SAMPLES)

Tumor_samples = getnames(Tumor, PLATFORM, PATH_BAM)
Normal_samples = getnames(Normal, PLATFORM, PATH_BAM)
AllFiles = Tumor_samples
AllFiles.update(Normal_samples)

#print(Tumor_samples)
#print(Normal_samples)

pairs=dict(zip(Tumor, Normal)) # dictionary containing Tumorname as key and normalname as value
print(pairs)

def getNormalSample(wildcards):
    normal = pairs[wildcards.sample]
    print(normal)   
    #name=re.match('[a-zA-Z0-9\-]*', normal).group(0)
    name = normal
    return(name)

# def getTumorSample(wildcards):
#     #name=re.match('[a-zA-Z0-9\-]*', wildcards.sample).group(0)
#     name = wildcards.sample
#     print(name)
#     return(name)


rule all:
    input:
        expand(path.join(PATH_CNV, '{binsize}', '{sample}CopywriteR_seg.rds'), sample=Tumor, binsize=['100kb'])


rule preCopywriteR:
    output:
        'reference/{ref}_{binsize}_{chrname}'
    params:
        binsize = lambda wildcards: int(wildcards.binsize.split('kb')[0])*1000,
        ref = '{ref}',
        chrname = '{chrname}',
        path = 'reference'
    conda: 'envs/copywriter.yaml'
    script:
        'script/preCopywriteR.R'


rule index:
    input:
    	'{sample}.bam'
    output:
    	'{sample}.bam.bai'
    conda: 'envs/samtools.yaml'
    shell:
    	"""
	samtools index {input}
	"""

rule filtering:
    input:
    	bam='{sample}.bam',
	index='{sample}.bam.bai'
    output:
    	'{sample}.filtered.bam'
    conda: 'envs/samtools.yaml'
    shell:
    	"""
	samtools view -b {input.bam} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > {output}
	"""


rule CopywriteR:
    input:
        tumor = path.join(PATH_BAM, '{sample}.filtered.bam'),
        control = lambda wildcards: path.join(PATH_BAM, getNormalSample(wildcards)+'.filtered.bam'),
        reference = 'reference/hg19_{binsize}_chr'
    output:
        path.join(PATH_CNV, '{binsize}', '{sample}CopywriteR_seg.rds')
    params:
        manifest = config['all']['bait_intervals'],
        binsize = lambda wildcards: int(wildcards.binsize.split('kb')[0])*1000,
        pairs = pairs,
        output = path.join(PATH_CNV, '{binsize}', '{sample}')
    conda: 'envs/copywriter.yaml'
    threads: 1 #20
    script:
        'script/CopywriteR.R'
