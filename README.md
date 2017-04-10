# mtPipeline

mtPipeline is a bioinformatics pipeline written by Tina Quach (2015). 
It is aimed at analyzing mitchondrial DNA and was created in order to augment 
current genome sequence analysis pipelines that exclude mitochondrial DNA analysis. 
The main goal is to call variants in the mitochondrial genome. 

mtPipeline incorporates the variant calling portion of 
the mitochondrial analysis pipeline in MToolBox (PMID:25028726)

##System Requirements
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


##What's included in the mtPipeline package?
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
completeVCF.py   		(generates completed.vcf from combined.vcf using info extracted from the samples' assembly-table.txt to fill missing info)
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


##What happens in mtPipeline? (Outline + Specifications)
###INPUTS
The pipeline needs two things: 
-path to the parameters.sh file (should be in mtPipeline folder; set the paths as needed)
-path to the folder that contains a subdirectory for each sample. 
Each  subdirectory is expected to contain part.bam files to be merged to create a single bam file for the sample.

###PROCESSING
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

###OUTPUTS
IN THE COHORT DIRECTORY (e.g. ID18)
- directory for each SAMPLE (e.g. ID18_Father, ID18_Mother, ID18_Proband)
- "VCF" directory 
	- "tables" directory (holds copy of each sample's assembly-table.txt)
	- a VCF for each sample (<sampleName>_VCF.vcf)
	- combined.vcf (output of combineVCF.sh; Father, Mother, Proband order)
	- completed.vcf (output of completeVCF.py; Father, Mother, Proband order; no missing info)
	- ordered.vcf (FINAL VCF OUTPUT of completeVCF.py; Proband, Mother, Father order.)
- log.txt (the outputs of each script called in the pipeline is written to the log for easier debugging)

IN EACH SAMPLE DIRECTORY (e.g. ID18_Proband)
- "bam" Directory:
	- "chrM" Directory: (files created during mtDNA extraction and analysis by simplepipe.py)
		*mtExtract.bam
		*mtExtract.bam.bai
		*mtExtract.qsort.bam (sorted by query name)
		*mtExtract_1.fq & *mtExtract_2.fq (created for remapping of mtDNA to chrRCRS reference sequence)
		*mtExtractremap.bam
		*mtExtactremap.sam
		*remap.csort.bam
	- "myMtoolbox_out" Directory:
		(KEY: 
			rg = read groups added (Picard)
			ra = indel realigned (GATK)
			marked = duplicates marked (GATK)
			norg = header removed
			sorted = sorted by coordinate (samtools))
		*rg.bam 
		*rg.bam bai
		*rg.ra.bam
		*rg.ra.bam.bai
		*rg.ra.marked.bam
		*rg.ra.marked.bam.sam
		*rg.ra.marked.norg.sam
		*rg.ra.marked-sorted.bam (USED IN MYASSEMBLEMTGENOME)
		*rg.ra.marked-sorted.bam.bai		
	- "partbam" Directory (contains the initial input files indexed (part.bam) and those files merged and indexed):
		 *part.bam
		 *part.bam.bai
		 *exome.bam
		 *exome.bam.bai
-"results" Directory: 
NOTE: results are only organized here for reference--not in the actual folder.
	General:
		-log.txt (output of all the processes except for assembly (see logassemble.txt for that))
	Coverage:
		Before Remap
			- <sampleName>.depthAnalysis_beforeExtract.txt (coverage statistical analysis)
			- hist_chrM_<samplName>.pdf  (depth of coverage histogram)
		After Remap
			- hist_chrM_<sampleName>_exome_mtExtractremap.pdf (depth of coverage histogram after remapping chrM)
			- <sampleName>.depthAnalysis_afterExtract.txt (coverage statistical analysis)
	PreAssembly Processing:
		<sampleName>-metrics.txt
		<sampleName>.rg.ra.marked.pileup
	Assembly:
		logassemble.txt (the log file of the myAssembleMTgenome.py script)
		<sampleName>-mtDNAassembly-contigs.fasta (a fasta file including all reconstructed contigs)
		<sampleName>-mtDNAassembly-coverage.txt (a text file including the coverage per contig and per mitochondrial known annotation)
		<sampleName>-mtDNAassembly-table.txt (the main table describing the assembly position by position)
	Variant Calling:
		<sampleName>_VCF.vcf (contains all the mitochondrial variant positions against RSRS and other meta-information)
-"tmp" Directory: used during analysis to store temporary files. Empty when processes are complete.

###TO DO

##Prepping Inputs
"-i" PREPPING THE PATH TO THE SAMPLE DIRECTORIES (aka cohort)
	This is the path to your cohort directory (e.g. ID18)
	It should hold ONLY the sample directories (e.g. ID18_Father, ID18_Mother, ID18_Proband)
	Directories must include slash at the end.
	These sample directories should ONLY contain  all the exome part.bam files for the particular sample along with their indexes.

"-p" PREPPING PARAMETERS.SH
	Make sure that the path to the various executables are correct.
	UPDATE the mtPipelineFolder variable (different for every user)
 	when inputting in the path to parameters.sh, DO NOT include the slash for parameters.sh

-Check EXAMPLE FILE NAMES (for part.bam files)
	ID18_Father_exome_L001_001.part.bam
	ID18_Mother_exome_L002_001.part.bam



##Running mtPipeline (Overview)

NOTE: Make sure you read the TO DO section below before running mtPipeline.

Basic execution of mtPipeline can be run as follows:
	bash mtPipeline.sh -i <pathToSampleDirs> -p <pathToParametersFile>
Both -i and -p options must be specified

For a complete list of mtPipeline.sh options please run as follows:
	bash mtPipeline.sh -h
IMPORTANT: if you wish to run again mtPipeline in the same folder, you should delete all the files produced during the previous execution.

TWO WAYS: 
- run the script directly on the command line as shown above 
- use the mtPipeline_qsub.sh as a template to submit a job
	- edit the "param" and "pathToSampleDirs" variables. 


# INSTALLATION

git clone https://github.com/quachtina96/mtPipeline.git
- make sure you know where your copy of the mtPipeline folder is
- edit the contents of mtPipeline/scripts/parameters.sh to match your environment


