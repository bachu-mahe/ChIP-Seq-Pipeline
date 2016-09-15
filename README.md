# CCBR ChIP-Seq Pipeline

##Installation

1. install nextflow
	curl -fsSL get.nextflow.io | bash

	For details, visit https://www.nextflow.io/

2. clone this git repository

	git clone https://github.com/CCBR/ChIP-Seq-Pipeline.git


##How to run the pipeline
	./nextflow run ChIP-Seq-Pipeline/main.nf --macsconfig='examples/macs.config' --reads='ChIP-Seq-Pipeline/examples/*.fastq' -config ChIP-Seq-Pipeline/config --genome='hg19'
	
* Note that the reads names are supposed to end with ".fastq" or ".fq" or "*.fastq.gz".
If you experience an error message related to missing file(s), this might be the reason.

* Config file is designed to work at our NIH biowulf2 HPC. If you plan to use on helix, add -profile 'local' option at the end of command line arguments. If you plan to use outside NIH, you need to edit the "genome" section according to your paths information.



###Additional options can be used: 
* -resume                         --> to resume the previous failed run
* -profile 'local'                --> for running tools "locally" not thourgh high performance computer queueing mechanisms.
* -with-timeline 'timeline.html'  --> record the run time of each process.
* -with-dag 'flowchart.png'       --> draw the flowchart

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
10. CEAS


##Todo list

1. ChipSeeker
2. Homer
3. PeakDiff

##Thanks
I got many nice implementation ideas from Nextflow examples, especially from NGI-ChIP-seq pipeline.
Many thanks to the NGI pipeline developers and our CCBR team mebers.


###Questions or Suggestions
Email to Bong-Hyun.Kim at NIH dot GOV.

