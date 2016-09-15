# CCBR ChIP-Seq Pipeline

##Installation

1. install nextflow
	curl -fsSL get.nextflow.io | bash

	For details, visit https://www.nextflow.io/

2. clone this git repository
	git clone https://github.com/CCBR/ChIP-Seq-Pipeline.git


##How to run the pipeline
	./nextflow run ChIP-Seq-Pipeline/main.nf --reads='ChIP-Seq-Pipeline/examples/*.fastq' -config ChIP-Seq-Pipeline/config --genome='hg19'
	
   the following additional options can be used 
        -resume
        -profile 'local' for 
        -with-timeline 'timeline.html'
	-with-dag 'flowchart.dag'


We implemented the pipeline using Nextflow.


## Current implementation includes following tools.
Thanks to the authors of the tools!

1. Trimgalor
2. BWA mem
3. Picard (MarkDuplicate)
4. FASTQC
5. DeepTools
6. Macs2
7. Sicer
8. MEME-ChIP
9. PhantomPeakqualTool


##Todo list

1. ChipSeeker
2. Homer
3. PeakDiff

##Thanks
I got many nice implementation ideas from Nextflow examples, especially from NGI-ChIP-seq pipeline.
Many thanks to the NGI pipeline developers and our CCBR team mebers.
