#!/bin/bash
#export PATH=/Volumes/RAID5/UTENTI_FTP/1000genomes/MToolBox:$PATH


check_exit_status()
{
rc=$?
if [[ $rc != 0 ]]
then
	echo ""
	echo "The last process reported an error. Exit."
	exit $rc
else
	echo "Success."
	echo ""
fi
}

usage()
{
	USAGE="""
	MToolBox: a tool for heteroplasmy annotation and accurate functional analysis of mitochondrial variants from high throughput sequencing data.
	Written by Domenico Simone, Claudia Calabrese and Maria Angela Diroma 2013-2014.
	https://sourceforge.net/projects/mtoolbox/

	You must run the MToolBox command on only one of the supported input file formats (bam, sam, fastq, fastq.gz, fasta).

	Input & workflow execution options (must include -i):

		-p	path to input folder.
		-i	input file format. Mandatory argument [bam|sam|fastq|fasta]. FASTQ files may be compressed or uncompressed.
		-o	path to output folder.
		-l	list of samples to be analyzed. [comma separated sample names]
		-X	extraction of mitochondrial reads from bam file, avoiding realignment of all bam file input
		-m	options for mapExome script [see mapExome.py -h for details]
		-M	remove duplicate reads with PicardTools MarkDuplicates after mapExome [default: no]
		-I	perform local realignment of reads on known indels with GATK IndelRealigner [default: no]
		-a	options for assembleMTgenome script [see assembleMTgenome.py -h for details]
		-c	options for mt-classifier script [see mt-classifier.py -h for details]
		-r	reference sequence to use for read mapping (VCF output will use the same reference) [RSRS|rCRS; default: RSRS]

	Help options:

		-h	show this help message
		-v	show version

	"""
	echo "$USAGE"
}

version()
{
	VERSION=$(echo "MToolBox v0.3")
	echo $VERSION
}

# Default command lines and behaviours for scripts and programs used in the workflow
#assembleMTgenome_OPTS=""
#mt_classifier_OPTS=""
#mapExome_OPTS=""
UseMarkDuplicates=false
UseIndelRealigner=false
MitoExtraction=false
# export folder where MToolBox.sh is placed, it is the same folder of PicardTools and GATK jars
me=`basename $0`
export mtoolbox_folder=$(which $me | sed "s/$me//g")
export externaltoolsfolder=${mtoolbox_folder}ext_tools/


# Environment variables for executables and files required by MToolBox
export ref=RSRS
export fasta_path=/usr/local/share/genomes/
export mtdb_fasta=chrRSRS.fa
export hg19_fasta=hg19RSRS.fa
export gsnapexe=/usr/local/bin/gsnap
export gsnapdb=/usr/local/share/gmapdb/
export mtdb=chrRSRS
export humandb=hg19RSRS
export samtoolsexe=/usr/local/bin/samtools
export muscleexe=/usr/local/bin/muscle


while getopts ":hva:c:f:p:o:i:l:m:r:j:MIX" opt; do
	case $opt in
		h)
			usage
			exit 1
			;;
		v)
			version
			exit 1
			;;
		a)
			assembleMTgenome_OPTS=$OPTARG
			;;
		c)
			mt_classifier_OPTS=$OPTARG
			;;
		f)
			variants_functional_annotation_OPTS=$OPTARG
			;;
		p)
			input_path=$OPTARG
			;;
		o)
			output_name=$OPTARG
			;;		
		i)
			input_type=$OPTARG
			;;			
		l)
			list=$OPTARG
			;;
		m)
			mapExome_OPTS=$OPTARG
			;;
		r)
			ref=$(echo $OPTARG | tr '[:lower:]' '[:upper:]')
			;;
		M)
			UseMarkDuplicates=true
			;;
		I)
			UseIndelRealigner=true
			;;
		X)
			MitoExtraction=true
			;;						
		\?)
			echo "Invalid option: -$OPTARG" >&2
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			exit 1
			;;
	esac
done

# define reference
if [[ $ref == 'RCRS' ]]
then 
	export mtdb_fasta=chrRCRS.fa
	export hg19_fasta=hg19RCRS.fa
	export mtdb=chrRCRS
	export humandb=hg19RCRS
elif [[ $ref != 'RSRS' ]]
then
	echo "Reference name not valid. Abort."
	exit 1
fi


# The following lines are commented since the involved parameters
# are specified elsewhere
#
# Set thresholds for hf and tail
#export hfthreshold=0.8
#export taillength=7
#

#if (( $taillength < 5 ))
#then
#	echo "Minumum tail length required >= 5. Tail length will be set to 5."
#	export taillength=5
#i

# Check python version (2.7 required)
echo ""
echo "Check python version... (2.7 required)"
min=$(python -c "import sys; print (sys.version_info[:])[1]")
maj=$(python -c "import sys; print (sys.version_info[:])[0]")
if [[ $maj != 2 ]] || [[ $min != 7 ]]
then
echo "You need Python2.7 in order to run MToolBox. Abort."
exit 1
else
echo "OK."
echo ""
fi


# Check existence of files to be used in the pipeline; if any of them does not exist, the pipeline will be aborted.
echo "Checking files to be used in MToolBox execution..."

#-t ${mt_classifier_OPTS} \

check_files.py \
--assembleMTgenome_OPTS="${assembleMTgenome_OPTS}" \
--mapExome_OPTS="${mapExome_OPTS}" \
--mt_classifier_OPTS="${mt_classifier_OPTS}"
# Check exit status of check_files.py
rc=$?
if [[ $rc != 0 ]] ; then
	exit $rc
fi

# Function definition
in-out_folders()
{ # create output folder and enter input folder
	if [[ "${input_path}" ]]
	then
		cd ${input_path}
	fi
	if [[ "${output_name}" ]]
	then
		mkdir ${output_name}
	fi
}

fastq_input()
{ # run mapExome directly.
	# get unique list of sample IDs
	# sampleIDs=$(ls *fastq* | awk 'BEGIN{FS="."}{count[$1]++}END{for (j in count) print j}')

	# map against mt genome and human genome
	# for i in $sampleIDs; do datasets=$(echo $i.*fastq*); mapExome_RSRS_SamHeader.py -g ${gsnapexe} -D ${gsnapdb} -M ${mtdb} -H ${humandb} -a "${datasets}" -o ${output_name}/OUT_${i}; done &> log_mapexome.txt
	echo ""
	echo "##### EXECUTING READ MAPPING WITH MAPEXOME..."
	echo ""
	if [[ ! "${list}" ]]
	#all the input files
	then
		if [ $input_type = 'bam' ]
		then
			sampleIDs=$(ls *.bam | grep -v ".MT.bam" | grep -v ".sorted.bam" | awk 'BEGIN{FS="."}{print $1}')
		elif [ $input_type = 'sam' ]
		then
			sampleIDs=$(ls *.sam | awk 'BEGIN{FS="."}{print $1}')			
		else			
			sampleIDs=$(ls *fastq* | awk 'BEGIN{FS="."}{count[$1]++}END{for (j in count) print j}')
		fi
	elif [ -s list.txt ]
	#samples in list.txt
	then
		sampleIDs=$(cat list.txt | awk 'BEGIN{FS="."}{count[$1]++}END{for (j in count) print j}')	
	else
	#list of files defined
		sampleIDs=$(echo "${list}" | sed 'y/,/\n/' | awk 'BEGIN{FS="."}{count[$1]++}END{for (j in count) print j}')	
	fi
	echo ""
		
	for i in $sampleIDs; do
		if [[ "${output_name}" ]]
		#output folder defined
		then
			if [ $input_type = 'bam' -o $input_type = 'sam' ]
			then
				#fastq files are in output folder
				echo "mapExome for sample" ${i}", files found:" $(ls ${output_name}/$i.*fastq*)					
				cd ${output_name}
			else
				#fastq files are in input folder
				echo "mapExome for sample" ${i}", files found:" $(ls $i.*fastq*)
			fi
			if (( $(ls $i.*fastq* | wc -l) == 1 ))
			then
				#echo $i is 1
				mapExome.py -g ${gsnapexe} -D ${gsnapdb} -M ${mtdb} -H ${humandb} -a $i.fastq* -o ${output_name}/OUT_${i} ${mapExome_OPTS}
			elif (( $(ls $i.*fastq* | wc -l) == 2 ))
			then
				#echo $i is 2
				mapExome.py -g ${gsnapexe} -D ${gsnapdb} -M ${mtdb} -H ${humandb} -a $i.R1.fastq* -b $i.R2.fastq* -o ${output_name}/OUT_${i} ${mapExome_OPTS}
			elif (( $(ls $i.*fastq* | wc -l) == 3 ))
			then
				if [ -s $i.fastq* ]
				then 
					mapExome.py -g ${gsnapexe} -D ${gsnapdb} -M ${mtdb} -H ${humandb} -a $i.R1.fastq* -b $i.R2.fastq* -c $i.fastq* -o ${output_name}/OUT_${i} ${mapExome_OPTS}
				else
					rm $i.fastq*
					echo "$i.fastq is an empty unpaired fastq. File has been removed."
					mapExome.py -g ${gsnapexe} -D ${gsnapdb} -M ${mtdb} -H ${humandb} -a $i.R1.fastq* -b $i.R2.fastq* -o ${output_name}/OUT_${i} ${mapExome_OPTS}
				fi			
			#then
				#echo $i is 3
				#mapExome.py -g ${gsnapexe} -D ${gsnapdb} -M ${mtdb} -H ${humandb} -a $i.R1.fastq* -b $i.R2.fastq* -c $i.fastq* -o ${output_name}/OUT_${i} ${mapExome_OPTS}
			else (( $(ls $i.*fastq* | wc -l) > 3 ))
				echo "$i not processed. Too many files."
				:
			fi
		else
		#no output folder defined
			echo "mapExome for sample" ${i}", files found:" $(ls $i.*fastq*)
			if (( $(ls $i.*fastq* | wc -l) == 1 ))
			then
				#echo $i is 1
				mapExome.py -g ${gsnapexe} -D ${gsnapdb} -M ${mtdb} -H ${humandb} -a $i.fastq* -o OUT_${i} ${mapExome_OPTS}
			elif (( $(ls $i.*fastq* | wc -l) == 2 ))
			then
				#echo $i is 2
				mapExome.py -g ${gsnapexe} -D ${gsnapdb} -M ${mtdb} -H ${humandb} -a $i.R1.fastq* -b $i.R2.fastq* -o OUT_${i} ${mapExome_OPTS}
			elif (( $(ls $i.*fastq* | wc -l) == 3 ))
			then
				if [ -s $i.fastq* ]
				then 
					mapExome.py -g ${gsnapexe} -D ${gsnapdb} -M ${mtdb} -H ${humandb} -a $i.R1.fastq* -b $i.R2.fastq* -c $i.fastq* -o OUT_${i} ${mapExome_OPTS}
				else
					rm $i.fastq*
					echo "$i.fastq is an empty unpaired fastq. File has been removed."
					mapExome.py -g ${gsnapexe} -D ${gsnapdb} -M ${mtdb} -H ${humandb} -a $i.R1.fastq* -b $i.R2.fastq* -o OUT_${i} ${mapExome_OPTS}
				fi	
			#then
				#echo $i is 3
				#mapExome.py -g ${gsnapexe} -D ${gsnapdb} -M ${mtdb} -H ${humandb} -a $i.R1.fastq -b $i.R2.fastq -c $i.fastq -o OUT_${i} ${mapExome_OPTS}
			else (( $(ls $i.*fastq* | wc -l) > 3 ))
				echo "$i not processed. Too many files."
				:
			fi
		fi				
	done

	echo ""
	

	if [ $input_type = 'bam' -o $input_type = 'sam' ]
	then
		echo "Compression of fastq files from bam/sam input files..."
		if [[ "${output_name}" ]]
		#output folder defined
		then
			mkdir ${output_name}/processed_fastq
			for i in $sampleIDs; do mv $i*fastq ${output_name}/processed_fastq; done
			cd ${output_name}							
		else
		#no output folder defined
			mkdir processed_fastq
			for i in $sampleIDs; do mv $i*fastq processed_fastq; done
		fi	
		tar czf processed_fastq.tar.gz processed_fastq
		rm -r processed_fastq 
		echo "Done."
	fi

	if [[ "${output_name}" ]]
	then
		cd ${output_name}
	fi
	
	echo ""
	echo "SAM files post-processing..."
	echo ""
	# SORT SAM WITH PICARD TOOLS
	echo "##### SORTING OUT.sam FILES WITH PICARDTOOLS..."
	echo ""
	for i in $(ls -d OUT_*); do cd ${i}; java -Xmx4g \
	-Djava.io.tmpdir=`pwd`/tmp \
	-jar ${externaltoolsfolder}SortSam.jar \
	SORT_ORDER=coordinate \
	INPUT=OUT.sam \
	OUTPUT=OUT.sam.bam \
	TMP_DIR=`pwd`/tmp; cd ..; done
	check_exit_status
	# INDEXING BAM FILES WITH SAMTOOLS
	for i in $(ls -d OUT_*); do cd ${i}; ${samtoolsexe} index OUT.sam.bam; cd ..; done
	
	# REALIGN KNOWN INDELS WITH GATK
	if $UseIndelRealigner
	then
		echo ""
		echo "##### REALIGNING KNOWN INDELS WITH GATK INDELREALIGNER..."
		echo ""
		for i in $(ls -d OUT_*); do cd ${i}; \
		echo "Realigning known indels for file" ${i}"/OUT.sam.bam using" ${mtoolbox_folder}"data/MITOMAP_HMTDB_known_indels.vcf as reference..."
		java -Xmx4g \
		-Djava.io.tmpdir=`pwd`/tmp \
		-jar ${externaltoolsfolder}GenomeAnalysisTK.jar \
		-T IndelRealigner \
		-R ${mtoolbox_folder}/data/chr${ref}.fa \
		-I OUT.sam.bam \
		-o OUT.realigned.bam \
		-targetIntervals ${mtoolbox_folder}/data/intervals_file_${ref}.list  \
		-known ${mtoolbox_folder}/data/MITOMAP_HMTDB_known_indels_${ref}.vcf \
		-compress 0;
		check_exit_status; cd ..; done
	else
		for i in $(ls -d OUT_*); do cd ${i}; cat OUT.sam.bam > OUT.realigned.bam; cd ..; done
	fi
	# MARK DUPLICATES WITH PICARD TOOLS
	if $UseMarkDuplicates
	then
		echo ""
		echo "##### ELIMINATING PCR DUPLICATES WITH PICARDTOOLS MARKDUPLICATES..."
		echo ""
		for i in $(ls -d OUT_*); do cd ${i}; java -Xmx4g \
		-Djava.io.tmpdir=`pwd`/tmp \
		-jar ${externaltoolsfolder}MarkDuplicates.jar \
		INPUT=OUT.realigned.bam \
		OUTPUT=OUT.sam.bam.marked.bam \
		METRICS_FILE=OUT.sam.bam.metrics.txt \
		ASSUME_SORTED=true \
		REMOVE_DUPLICATES=true \
		TMP_DIR=`pwd`/tmp; cd ..; done
	else
		for i in $(ls -d OUT_*); do cd ${i}; cat OUT.realigned.bam > OUT.sam.bam.marked.bam; cd ..; done
	fi
	# RE-CONVERT BAM OUTPUT FROM MARKDUPLICATES IN SAM.
	for i in $(ls -d OUT_*); do cd ${i}; java -Xmx4g -Djava.io.tmpdir=`pwd`/tmp -jar ${externaltoolsfolder}SamFormatConverter.jar INPUT=OUT.sam.bam.marked.bam OUTPUT=OUT.sam.bam.marked.bam.marked.sam TMP_DIR=`pwd`/tmp; cd ..; done

	for i in $(ls -d OUT_*); do cd ${i}; grep -v "^@" *marked.sam > OUT2.sam; mkdir MarkTmp; mv OUT.sam.bam MarkTmp; mv OUT.sam.bam.marked.bam MarkTmp; mv OUT.sam.bam.marked.bam.marked.sam MarkTmp; tar -czf MarkTmp.tar.gz MarkTmp; rm -R MarkTmp/; cd ..; done

	# ASSEMBLE CONTIGS, GET MT-TABLES...
	echo ""
	echo "##### ASSEMBLING MT GENOMES WITH ASSEMBLEMTGENOME..."
	echo ""
	echo "WARNING: values of tail < 5 are deprecated and will be replaced with 5"
	echo ""	
	for i in $(ls -d OUT_*); do outhandle=$(echo ${i} | sed 's/OUT_//g'); cd ${i}; assembleMTgenome.py -i OUT2.sam -o ${outhandle} -r ${fasta_path} -f ${mtdb_fasta} -a ${hg19_fasta} -s ${samtoolsexe} -FCP ${assembleMTgenome_OPTS}; cd ..; done > logassemble.txt
	echo ""
	echo "##### GENERATING VCF OUTPUT..."
	# ... AND VCF OUTPUT
	VCFoutput.py -r ${ref}
}

fasta_input()
{
	if [[ $input_type = 'fasta' ]]
	then
		echo ""
		echo "##### PRE-PROCESSING OF FASTA INPUT FILES..."
		echo ""
		echo "Files to be analyzed:"
		if [[ "${list}" ]]
		then
			if [ -s list.txt ]
			#samples in list.txt
			then
				filelist=$(cat list.txt | tr '\n' '\t')
			else
			#list of input files defined
				filelist=$(echo "${list}" | tr ',' '\t')
			fi
		fi	
		
		
		if [[ "${output_name}" ]]
		#output folder defined
		then
			if [[ ! "${list}" ]]
			#all the input files
			then
				for i in $(test_fasta.py); do bname=$(echo ${i} | awk 'BEGIN{FS="."}{print $1}'); bname_dir=OUT_${bname}; mkdir ${output_name}/${bname_dir}; cp ${i} ${output_name}/${bname_dir}/${bname}-contigs.fasta; echo ${bname}; done
				#for i in $(ls); do bname=$(echo ${i} | awk 'BEGIN{FS="."}{print $1}'); bname_dir=OUT_${bname}; mkdir ${bname_dir}; cp ${i} ${bname_dir}/${bname}-contigs.fasta; done
			else
				for i in $filelist; do bname=$(echo ${i} | awk 'BEGIN{FS="."}{print $1}'); bname_dir=OUT_${bname}; mkdir ${output_name}/${bname_dir}; cp ${i} ${output_name}/${bname_dir}/${bname}-contigs.fasta; echo ${bname}; done
				#for i in $(ls); do bname=$(echo ${i} | awk 'BEGIN{FS="."}{print $1}'); bname_dir=OUT_${bname}; mkdir ${bname_dir}; cp ${i} ${bname_dir}/${bname}-contigs.fasta; done
			fi
		else
		#no output folder defined
			if [[ ! "${list}" ]]
			#all the input files
			then
				for i in $(test_fasta.py); do bname=$(echo ${i} | awk 'BEGIN{FS="."}{print $1}'); bname_dir=OUT_${bname}; mkdir "${bname_dir}"; cp ${i} ${bname_dir}/${bname}-contigs.fasta; echo ${bname}; done
				#for i in $(ls); do bname=$(echo ${i} | awk 'BEGIN{FS="."}{print $1}'); bname_dir=OUT_${bname}; mkdir ${bname_dir}; cp ${i} ${bname_dir}/${bname}-contigs.fasta; done
			else
			#list of input files defined	
				for i in $filelist; do bname=$(echo ${i} | awk 'BEGIN{FS="."}{print $1}'); bname_dir=OUT_${bname}; mkdir "${bname_dir}"; cp ${i} ${bname_dir}/${bname}-contigs.fasta; echo ${bname}; done
				#for i in $(ls); do bname=$(echo ${i} | awk 'BEGIN{FS="."}{print $1}'); bname_dir=OUT_${bname}; mkdir ${bname_dir}; cp ${i} ${bname_dir}/${bname}-contigs.fasta; done
			fi				
		fi
		check_exit_status
	fi	
	echo ""
	echo "##### PREDICTING HAPLOGROUPS AND ANNOTATING/PRIORITIZING VARIANTS..."
	echo ""
	if [[ "${output_name}" ]]
	then
		cd "${output_name}"
	fi
	#### Haplogroup prediction and functional annotation
	# Brand new haplogroup prediction best file
	hpbest="mt_classification_best_results.csv" # change just this name for changing filename with most reliable haplogroup predictions
	echo "Haplogroup predictions based on RSRS Phylotree build 16"
	echo "SampleID,Best predicted haplogroup(s)" > ${hpbest}
	for i in $(ls -d OUT_*); do inhandle=$(echo ${i} | sed 's/OUT_//g'); cd ${i}; mt-classifier.py -i ${inhandle}-contigs.fasta -s ${hpbest} -b ${inhandle} -m ${muscleexe} ${mt_classifier_OPTS}; cd ..; done

	# Functional annotation of variants
	#for i in $(ls -d OUT_*); do cd $i; variants_functional_annotation.py $hpbest ; cd ..; done
	variants_functional_annotation.py #${hpbest}
	# Collect all prioritized variants from all the samples
	for i in $(ls -d OUT_*/*annotation.csv); do tail -n+2 $i | awk 'BEGIN {FS="\t"}; {if ($5 == "yes" && $6 == "yes" && $7 == "yes") {print $1"\t"$2"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$28"\t"$29"\t"$30"\t"$31"\t"$32"\t"$33"\t"$34"\t"$35"\t"$36"\t"$37"\t"$38"\t"$39}}' >> priority_tmp.txt; done
	for i in $(ls -d OUT_*/*annotation.csv); do tail -n+2 $i | awk 'BEGIN {FS="\t"}; {if ($5 == "yes" && $6 == "yes" && $7 == "yes") count++} END {print $1"\t"NR"\t"count}' >> variant_number.txt; done
	prioritization.py priority_tmp.txt
	rm priority_tmp.txt
	echo ""
	echo "Prioritization analysis done."
	echo ""
	if [[ $input_type = 'fasta' ]]
	then
		summary.py
	else
		for i in $(ls -d OUT_*); do name=$(echo $i | sed 's/OUT_//g'); cd $i; coverage=$(cat *coverage.txt | grep "Assemble"); cd ..; echo "Sample:" "$name" "$coverage"; done >> coverage_tmp.txt
		if [[ "${assembleMTgenome_OPTS}" ]]
		then
			HFthreshold=$(echo "$assembleMTgenome_OPTS" | grep -oh "\w*-z[[:space:]][0-9]\.[0-9]\w*" | tr '\ ' '\n' | awk 'NR==2')
			REdistance=$(echo "$assembleMTgenome_OPTS" | grep -oh "\w*-t[[:space:]][0-300]\w*" | tr '\ ' '\n' | awk 'NR==2')
			if [ -z "$HFthreshold" ]
			then
				HFthreshold=$(echo "0.8")
			elif [ -z "$REdistance" ]
			then
				REdistance=$(echo "5")
			fi
		else
			HFthreshold=$(echo "0.8")
			REdistance=$(echo "5")
		fi	
		#if [[ "${assembleMTgenome_OPTS}" ]]
		#then
		for i in $(ls -d OUT_*); do name=$(echo $i | sed 's/OUT_//g'); cd $i; heteroplasmy=$(echo "$HFthreshold"); homo_variants=$(awk 'BEGIN {FS="\t"}; {if ($3 == "1.0") count++} END {print count}' *annotation.csv); above_threshold=$(awk -v thrsld=$heteroplasmy 'BEGIN {FS="\t"};{if ( $3 >= thrsld && $3 < "1.0" ) count++} END {print count}' *annotation.csv); under_threshold=$(awk -v thrsld=$heteroplasmy 'BEGIN {FS="\t"};{if ( $3 < thrsld && $3 > "0" ) count++} END {print count}' *annotation.csv); cd ..; echo "$name" "$homo_variants" "$above_threshold" "$under_threshold"; done >> heteroplasmy_count.txt
		#else
		#	for i in $(ls -d OUT_*); do name=$(echo $i | sed 's/OUT_//g'); cd $i; homo_variants=$(awk 'BEGIN {FS="\t"}; {if ($3 == "1.0") count++} END {print count}' *annotation.csv); above_threshold=$(awk 'BEGIN {FS="\t"};{if ( $3 >= "0.8" && $3 < "1.0" ) count++} END {print count}' *annotation.csv); under_threshold=$(awk 'BEGIN {FS="\t"};{if ( $3 < "0.8" && $3 > "0" ) count++} END {print count}' *annotation.csv); cd ..; echo "$name" "$homo_variants" "$above_threshold" "$under_threshold"; done >> heteroplasmy_count.txt
		#fi		
		summary.py coverage_tmp.txt heteroplasmy_count.txt
		rm coverage_tmp.txt
		rm heteroplasmy_count.txt
	fi
	rm variant_number.txt
	if [[ $input_type = 'fasta' ]]
	then
		echo -e "Selected input format\t$(echo "$input_type")\nReference sequence used for haplogroup prediction\tRSRS\n\n==============================\n\n$(cat summary_tmp.txt)\n\n==============================\n\nTotal number of prioritized variants\t$(awk 'END{print NR-1}' prioritized_variants.txt)" > summary_`date +%Y%m%d_%H%M%S`.txt
	else		
		echo -e "Selected input format\t$(echo "$input_type")\nReference sequence chosen for mtDNA read mapping\t$(echo "$ref")\nReference sequence used for haplogroup prediction\tRSRS\nDuplicate read removal?\t$(echo "$UseMarkDuplicates")\nLocal realignment around known indels?\t$(echo "$UseIndelRealigner")\nMinimum distance of indels from read end\t$(echo "$REdistance")\nHeteroplasmy threshold for FASTA consensus sequence\t$(echo "$HFthreshold")\n\nWARNING: If minimum distance of indels from read end is < 5, it is deprecated and replaced with 5\n\n==============================\n\n$(cat summary_tmp.txt | sed "s/thrsld/$HFthreshold/g")\n\n==============================\n\nTotal number of prioritized variants\t$(awk 'END{print NR-1}' prioritized_variants.txt)"  >  summary_`date +%Y%m%d_%H%M%S`.txt
	fi	
	rm summary_tmp.txt
	echo ""
	echo "Analysis completed!"

}

sam_input()
{ # convert sam to fastq.
	echo ""
	if [[ ! "${list}" ]]
	#all the input files
	then
		sam_samples=$(ls *.sam | awk 'BEGIN{FS="."}{print $1}')
	elif [ -s list.txt ]
	#samples in list.txt
	then
		sam_samples=$(cat list.txt | awk 'BEGIN{FS="."}{print $1}')
	else
	#list of input files defined
		sam_samples=$(echo "${list}" | sed 's/,/\n/g' | awk 'BEGIN{FS="."}{print $1}')
	fi	
		
	if [[ "${output_name}" ]]
	#output folder defined
	then
		for i in ${sam_samples}; do echo "Converting sam to fastq..." ${i}.sam; java -Xmx4g \
		-Djava.io.tmpdir=`pwd`/tmp \
		-jar ${externaltoolsfolder}SamToFastq.jar \
		INPUT=${i}.sam \
		FASTQ=${output_name}/${i}.R1.fastq \
		SECOND_END_FASTQ=${output_name}/${i}.R2.fastq \
		UNPAIRED_FASTQ=${output_name}/${i}.fastq \
		VALIDATION_STRINGENCY=SILENT \
		TMP_DIR=${output_name}/tmp; echo "Done."; done
	else
	#no output folder defined
		for i in ${sam_samples}; do echo "Converting sam to fastq..." ${i}.sam; java -Xmx4g \
		-Djava.io.tmpdir=`pwd`/tmp \
		-jar ${externaltoolsfolder}SamToFastq.jar \
		INPUT=${i}.sam \
		FASTQ=${i}.R1.fastq \
		SECOND_END_FASTQ=${i}.R2.fastq \
		UNPAIRED_FASTQ=${i}.fastq \
		VALIDATION_STRINGENCY=SILENT \
		TMP_DIR=`pwd`/tmp; echo "Done."; done
	fi
#echo "Converting sam input(s) to fastq"
#echo "Done."
}

bam_input()
{ # convert bam to fastq or extract mitochondrial reads from bam and then convert MT.bam file in fastq.
	echo ""
	if [[ ! "${list}" ]]
	#all the input files
	then
		bam_samples=$(ls *.bam | awk 'BEGIN{FS="."}{print $1}')
	elif [ -s list.txt ]
	#samples in list.txt
	then
		bam_samples=$(cat list.txt | awk 'BEGIN{FS="."}{print $1}')
	else
	#list of input files defined
		bam_samples=$(echo "${list}" | sed 's/,/\n/g' | awk 'BEGIN{FS="."}{print $1}')
	fi	
		
	if $MitoExtraction
	#extract mitochondrial reads from bam input files and then convert in fastq files
	then
		if [[ "${output_name}" ]]
		#output folder defined
		then
			for i in ${bam_samples}; do echo "Sorting, indexing and extraction of mitochondrial reads from bam file..." ${i}.bam; samtools sort $i.bam ${output_name}/$i.sorted; samtools index ${output_name}/$i.sorted.bam; samtools view -b ${output_name}/$i.sorted.bam MT M chrMT chrM > ${output_name}/$i.MT.bam; echo "Done."; done
			echo ""
			for i in $(ls ${output_name}/*.MT.bam); do echo "Converting bam to fastq..." ${i}; n=$(echo $i | awk 'BEGIN{FS="."}{print $1}'); java -Xmx4g \
				-Djava.io.tmpdir=`pwd`/tmp \
				-jar ${externaltoolsfolder}SamToFastq.jar \
				INPUT=${n}.MT.bam \
				FASTQ=${n}.R1.fastq \
				SECOND_END_FASTQ=${n}.R2.fastq \
				UNPAIRED_FASTQ=${n}.fastq \
				VALIDATION_STRINGENCY=SILENT \
				TMP_DIR=${output_name}/tmp; echo "Done."; done
			mkdir ${output_name}/processed_bam
			mv ${output_name}/*MT.bam ${output_name}/*bai ${output_name}/*sorted.bam ${output_name}/processed_bam
			echo ""
			echo "Compression of processed bam files..."
			tar -czf ${output_name}/processed_bam.tar.gz ${output_name}/processed_bam
			rm -r ${output_name}/processed_bam
		else
		#no output folder defined
			for i in ${bam_samples}; do echo "Sorting, indexing and extraction of mitochondrial reads from bam file..." ${i}.bam; samtools sort $i.bam $i.sorted; samtools index $i.sorted.bam; samtools view -b $i.sorted.bam MT M chrMT chrM > $i.MT.bam; echo "Done."; done
			echo ""	
			for i in $(ls *.MT.bam); do echo "Converting bam to fastq..." ${i}; n=$(echo $i | awk 'BEGIN{FS="."}{print $1}'); java -Xmx4g \
				-Djava.io.tmpdir=`pwd`/tmp \
				-jar ${externaltoolsfolder}SamToFastq.jar \
				INPUT=${n}.MT.bam \
				FASTQ=${n}.R1.fastq \
				SECOND_END_FASTQ=${n}.R2.fastq \
				UNPAIRED_FASTQ=${n}.fastq \
				VALIDATION_STRINGENCY=SILENT \
				TMP_DIR=`pwd`/tmp; echo "Done."; done
			mkdir processed_bam
			mv *MT.bam *bai *sorted.bam processed_bam
			echo ""
			echo "Compression of processed bam files..."
			tar -czf processed_bam.tar.gz processed_bam	
			rm -r processed_bam			
		fi	
	else	
	#convert all bam input files in fastq files 
		if [[ "${output_name}" ]]
		#output folder defined
		then
			for i in ${bam_samples}; do echo "Converting bam to fastq..." ${i}.bam; java -Xmx4g \
			-Djava.io.tmpdir=`pwd`/tmp \
			-jar ${externaltoolsfolder}SamToFastq.jar \
			INPUT=${i}.bam \
			FASTQ=${output_name}/${i}.R1.fastq \
			SECOND_END_FASTQ=${output_name}/${i}.R2.fastq \
			UNPAIRED_FASTQ=${output_name}/${i}.fastq \
			VALIDATION_STRINGENCY=SILENT \
			TMP_DIR=${output_name}/tmp; echo "Done."; done
		else
		#no output folder defined
			for i in ${bam_samples}; do echo "Converting bam to fastq..." ${i}.bam; java -Xmx4g \
			-Djava.io.tmpdir=`pwd`/tmp \
			-jar ${externaltoolsfolder}SamToFastq.jar \
			INPUT=${i}.bam \
			FASTQ=${i}.R1.fastq \
			SECOND_END_FASTQ=${i}.R2.fastq \
			UNPAIRED_FASTQ=${i}.fastq \
			VALIDATION_STRINGENCY=SILENT \
			TMP_DIR=`pwd`/tmp; echo "Done."; done
		fi
	fi
#echo "Converting bam input(s) to fastq"
#echo "Done."
}

	
if (( $# >= 1 ))
then
	if [[ $input_type = 'fasta' ]]
	then
		echo "Input type is fasta."
		in-out_folders
		fasta_input
	elif [[ $input_type = 'fastq' ]]
	then
		echo "Input type is fastq."
		in-out_folders
		fastq_input
		fasta_input
	elif [[ $input_type = 'sam' ]]
	then
		echo "Input type is sam."
		in-out_folders
		sam_input
		fastq_input
		fasta_input
	elif [[ $input_type = 'bam' ]]
	then
		echo "Input type is bam."
		in-out_folders
		bam_input
		fastq_input
		fasta_input
	else
		echo "Input format not recognized."
		exit 1
	fi
else
	echo "Input type not specified."
	exit 1
fi
