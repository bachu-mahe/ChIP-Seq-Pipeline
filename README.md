CCBR ChIP-Seq Pipeline

We implemented this pipeline using Nextflow.

Current implementation includes following tools.

1. Trimgalor
2. BWA mem
3. Picard (MarkDuplicate)
4. FASTQC
5. DeepTools
6. Macs2
7. Sicer
8. MEME-ChIP
9. PhantomPeakqualTool

Todo list

1. ChipSeeker
2. Homer
3. PeakDiff

I got many nice implementation ideas from Nextflow examples, especially from NGI-ChIP-seq pipeline.
Many thanks to the NGI pipeline developers.
