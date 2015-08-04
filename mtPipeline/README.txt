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
mtPipeline.qsub     	(an example script that you may adapt or use an an example to run mtPipeline)
mtPipeline.sh       	(the main script you will want to run)
simplepipe.py       	(python script that takes care of coverage analysis and mt extraction)
myMtoolbox.sh       	(bash script that process the extracted mitchondrial DNA bam file and calls MToolBox functions to analyze it--variant call)
parameters.sh       	(bash script that sets up the environment for the pipeline) 
myAssembleMTgenome.py   (augmented MToolBox script for assembling the MTgenome)
mtVariantCaller.py  	(MToolbox script for variant calling)
VCFoutput.py 			(MToolBox script for creating the VCF)
combineVCF.sh 			(combines the VCFs of a cohort into a single VCF for haplogroup analysis)
pipeFunctions.py 		(contains the sam/bam analysis functions called by simplepipe.py)
histDepth.R 			(creates a histogram of the depth of coverage for each; called by simplepipe.py)
parameters.sh 			(sourced in mtPipeline, denotes the paths to various tools and executables)


/ref:					(reference data = Revised Cambridge Reference Sequence)
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
The pipeline needs two things: 
-path to the parameters.sh file (should be in mtPipeline folder; set the paths as needed)
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

After that, a new VCF folder for the cohort is created and the VCF for the cohort are merged into a single "combined.vcf".
OUTPUTS:
IN THE COHORT DIRECTORY (e.g. ID18)
- directory for each SAMPLE (e.g. ID18_Father, ID18_Mother, ID18_Proband)
- "VCF" directory that contains a VCF for each sample & the combined VCF (Proband, Mother, Father order)
- log.txt (the outputs of each script called in the pipeline is written to the log for easier debugging)

IN EACH SAMPLE DIRECTORY (e.g. ID18_Proband)
- "chrM" Directory:
- "myMtoolbox_out" Directory:

- "partbam" Directory: contains the initial input files (part.bam), and those files merged and indexed.
	- *part.bam
	- *exome.bam
	- *exome.bam.bai
-"results" Directory:
	- 
-"tmp" Directory: used during analysis to store temporary files. Empty when processes are complete.


INSTALLATION
============
git clone [url to my git project]

or 

Download the mtPipeline folder off of SourceForge

Edit the variables in mtPipeline/scripts/parameters.sh that DO NOT start with an underscore in order to make the your paths are correct.

==========TODO========= 

--------------------
Running mtPipeline
-------------------

Basic execution of mtPipeline can be run as follows:

	bash mtPipeline.sh -i <pathToSampleDirs> -p <pathToParametersFile>

Both -i and -p options must be specified.
For a complete list of mtPipeline.sh options please run as follows:
	
	bash mtPipeline.sh -h

IMPORTANT: if you wish to run again mtPipeline in the same folder, you should delete all the files produced during the previous execution.