
#include "region.h"
#include "gene.h"

int main(int argc, char* argv[])
{   
	int start = 100000;//zero based
	int stop  = 101000;//zero based
	int chr_num = 1;
	char strand = '+';
	const char* gio_file = "/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana_magic/genomes/A_thaliana_magic/Col_0/Col_0.gio/genome.config";
	
	printf("test initialization of Region\n");

	Region* reg = new Region(start, stop, chr_num, strand, gio_file);

	reg->print(stdout);
	printf("contig_len: %i\n", reg->gio->contig_len(1));
	reg->load_genomic_sequence();
	//fprintf(stdout, "seq: %s\n", reg->seq);

	/*
	fprintf(stdout, "seq: ");
	char* seq = reg->gio->read_flat_file(chr_num);
	for (int i=start; i<=stop; i++)
	{
		fprintf(stdout, "%c", seq[i]);
	} 
	fprintf(stdout, "\n");
	*/

	char* bam_file = (char*) "/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana_magic/reads/Col_0_Chr_nss.bam";

	reg->get_reads(&bam_file, 1);

	Gene* gene = new Gene(start, stop, chr_num, strand, gio_file); 

	gene->get_reads(&bam_file, 1);

	exon exon1(100050, 100060); 
	exon exon2(100070, 100090); 
	gene->exons.push_back(exon1);
	gene->exons.push_back(exon2);

	gene->print(stdout);
}
