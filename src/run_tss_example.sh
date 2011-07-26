

gio_file=/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana_magic/genomes/A_thaliana_magic/Col_0/Col_0.gio/genome.config
bam_file=/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana_magic/reads/Col_0_Chr_nss.bam
gff_file=~/tmp/genes.gff3
reg_file=~/tmp/regions.txt
tss_file=~/tmp/tss
tis_file=~/tmp/tis
options="-maxel 7000 -minel 10 -reglen 0.25 -maxic 10000"
echo "./infer_genes -gio $gio_file -bam $bam_file -gff $gff_file -reg $reg_file $options"
./infer_genes -tss $tss_file -tis $tis_file -gio $gio_file -bam $bam_file -gff $gff_file -reg $reg_file $options  
