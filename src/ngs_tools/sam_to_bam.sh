#!/bin/bash
# $1 tracks_dir
# $2 genome_fasta
# [$3] sam file

set -e

basedr=$1
cd $basedr

fastaname=$2

#samtoolsdr="/fml/ag-raetsch/share/software/samtools/"
samtoolsdr="/fml/ag-raetsch/share/software/samtools.svn_up/"


if [ -e ${fastaname}.fai ]
then
    echo "${fastaname}.fai already exists"
else
    echo "${fastaname}.fai is generated"
    $samtoolsdr./samtools faidx $fastaname
fi

if [ $# -lt 3 ]; then 
	# apply to all sam files in basedr
	echo "Compressing SAM files..."
	# do gzip in a for loop 
	# then it will not fail if 
	# sam files are already zipped
	for f in `ls $basedr/*.sam`
	do
		gzip -9 $f
	done
	echo " done."

	sam_files=`ls ${1}/*.sam.gz`
else
	sam_files=$3
	if [[ ${sam_files} != *.gz ]]
	then
		echo gzip file $sam_files
		gzip -9 $sam_files
		sam_files=$3.gz
	fi
fi
    
echo $sam_files

for n in $sam_files
do
	echo "starting with file $n"
    m=`echo $n | sed "s/.sam.gz/.bam/g"`
	echo "creating file $m"
    if [ ! -f $m ]
    then
		echo $n
		echo "creating file $m"
		$samtoolsdr./samtools import ${fastaname}.fai $n $m
		# BAM must be sorted by start position to use random access
		$samtoolsdr./samtools sort $m ${m}_sorted
		mv ${m}_sorted.bam $m
		# index for BAM file
		$samtoolsdr./samtools index $m
    fi
done
