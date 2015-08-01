==================
mtPipeline READ ME
===================

mtPipeline is a bioinformatics pipeline written by Tina Quach (2015). It is aimed at analyzing mitchondrial DNA and was created in order to augment current genome sequence analysis pipelines that exclude mitochondrial DNA analysis. The main goal is to call variants in the mitochondrial genome. 

mtPipeline incorporates the variant calling portion of the mitochondrial analysis pipeline in MToolBox (PMID:25028726)

---------------------
System Requirements
---------------------
- UNIX-based OS
- Python2.7 (www.python.org)
for the following programs, specify path in parameters.sh:
- samtools
- GATK
- BWA
- Picard 
- bedtools
- Java
- R

-------------------------------------------
What's included in the mtPipeline package?
--------------------------------------------
/scripts:               (all the code for the pipeline)
mtPipeline.qsub     	(an example script of how to submit a job/run mtPipeline)
mtPipeline.sh       	(the main script you will want to run)
simplepipe.py       	(python script that takes care of coverage analysis and mt extraction)
myMtoolbox.sh       	(bash script that process the extracted mitchondrial DNA bam file and calls 					    MToolBox functions to analyze it--variant call)
parameters.sh       	(bash script that sets up the environment for the pipeline) 
myAssembleMTgenome.py   (augmented MToolBox script for assembling the MTgenome)
mtVariantCaller.py  	(MToolbox script for variant calling)
VCFoutput.py 			(MToolBox script for creating the VCF)
pipeFunctions.py 		(contains the sam/bam analysis functions called by simplepipe.py)
histDepth.R 			(creates a histogram of the depth of coverage for each; called by simplepipe.						 py)
parameters.sh 			(sourced in mtPipeline, denotes the paths to various tools and executables)


/ref:					(reference data)
chrRCRS.fa
chrRCRS.fa.fai
chrRCRS.fa.sa
chrRCRS.fa.pac
chrRCRS.fa.bwt
chrRCRS.fa.ann
chrRCRS.fa.amb
chrRCRS.dict

/externalTools
MToolBox 				


----------------------------------------------------------
What happens in mtPipeline? (Outline + Specifications)
----------------------------------------------------------
INPUTS:
The pipeline takes in two paths as inputs: 
-path to the mtPipeline folder 
-path to the folder that contains a subdirectory for each sample. 
Each  subdirectory is expected to contain part.bam files to be merged to create a single bam file for the sample.

For each sample directory, 
1) merge the part.bam files > exome.bam file
2) get and analyze depth of coverage of the chrM region of exome.bam file
3) extract chrM portion of exome.bam > mtExtract.bam
4) get and analyze the depth of coverage of the mtExtract.
5) process the resulting mtExtract file 
	-index (samtools index)
	-sort (samtools sort)
	-add read groups (Picard AddOrReplaceReadGroups)
	-execute indel realignment (GATK IndelRealigner)
	-mark duplicates (Picard MarkDuplicates)
6) assemble the mtGenome 
7) create pileup (samtools mpileup)
8) generate VCF 

OUTPUTS:
The pipeline results in: 
two histogram plots of the coverage before 
[c/v the stuff from mtoolbox regarding the stuff they wrote]

INSTALLATION
============
git clone [url to my git project]

or 

download the folder off sourceforge?

or

Download MToolBox.tar.gz. Decompress the file and copy the package folder in a folder of your choice. Add this folder to your system PATH with the following command:

export PATH=/dir/MToolBox/:$PATH


==========TODO========= 

--------------------
Running mtPipeline
-------------------

Basic execution of mtPipeline can be run as follows:

	mtPipeline.sh -i <pathToSampleDirs> -m <pathTomtPipelineFolder>

Both -i and -m options must be specified.
The command must be executed inside the folder containing your input files. Please note that only one of the supported formats can be used within a single MToolBox run.
For a complete list of MToolBox.sh options please run as follows:
	
	MToolBox.sh -h

IMPORTANT: if you wish to run again MToolBox in the same folder, it is preferable that you delete all the files produced during the previous execution.