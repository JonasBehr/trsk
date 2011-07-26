

gio_file=/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana_magic/genomes/A_thaliana_magic/Col_0/Col_0.gio/genome.config
bam_file=/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana_magic/reads/Col_0_Chr_nss.bam
gff_file=genes_noorf.gff3
options="-noorf"
./infer_genes -gio $gio_file -bam $bam_file -gff $gff_file -reg $reg_file $options  
