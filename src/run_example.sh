

gio_file=/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana_magic/genomes/A_thaliana_magic/Col_0/Col_0.gio/genome.config
bam_file=/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana_magic/reads/Col_0_Chr_nss.bam
gff_file=genes.gff3
reg_file=regions.txt
options="-maxel 7000 -minel 10 -reglen 0.25 -maxic 10000"
echo "./infer_genes -gio $gio_file -bam $bam_file -gff $gff_file -reg $reg_file $options"
valgrind --leak-check=full ./infer_genes -gio $gio_file -bam $bam_file -gff $gff_file -reg $reg_file $options  
