CCBR ChIP-Seq Pipeline

#Installation

1. install nextflow
	Visit https://www.nextflow.io/

2. clone this git repository


#Runing the pipeline
	nextflow run ChIP-Seq-Pipeline/main.nf --reads='ChIP-Seq-Pipeline/examples/*.fastq' -config ChIP-Seq-Pipeline/config --genome='hg19'
	
   the following additional options can be used 
        -resume
        -profile 'local'
        -with-timeline 'timeline.html'
	-with-dag 'flowchart.dag'


We implemented this pipeline using Nextflow.


# Current implementation includes following tools.

1. Trimgalor
2. BWA mem
3. Picard (MarkDuplicate)
4. FASTQC
5. DeepTools
6. Macs2
7. Sicer
8. MEME-ChIP
9. PhantomPeakqualTool


#Todo list

1. ChipSeeker
2. Homer
3. PeakDiff

#Thanks
I got many nice implementation ideas from Nextflow examples, especially from NGI-ChIP-seq pipeline.
Many thanks to the NGI pipeline developers and our CCBR team mebers.
