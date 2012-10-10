
sam=/tmp/all_sam_file.sam

for f in $*
do
	echo $f
	samtools view $f >> $sam
	echo -e "\n" >> $sam
done
gzip $sam
sh ./sam_to_bam.sh /tmp/ /fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana/genome/A_thaliana/A_thaliana.fasta $sam.gz


#sh ./merge_bam_files.sh /fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana/reads/SRX00*.new.bam 
