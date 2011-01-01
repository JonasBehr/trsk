#include <mex.h>

extern char *get_string(const mxArray *prhs);

/*
 * prhs[0] gio_file: Genome Information Object file
 * prhs[1] bam_file: Alignment file in bam format
 * prhs[2] gff_file: genes in gff3 format will be written to that file
 * prhs[3] reg_file: list of intergenic regions will be written to that file
 *
 * return:
 * plhs[0] 
 *
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	char* gio_file = get_string(prhs[0]);
	char* bam_file = get_string(prhs[1]);
	char* gff_file = get_string(prhs[2]);
	char* reg_file = get_string(prhs[3]);

	vector<Region*> regions = GeneTools::init_regions(gio_file);
	printf("regions.size(): %i\n", (int) regions.size());

	regions[0]->stop = 1e6;
	//regions[1]->stop = 1e6;
	vector<Gene*> genes;
	//for (int r=0; r<regions.size(); r++)
	for (int r=0; r<1; r++)
	{
		printf("Starting with contig %i, strand %c\n", regions[r]->chr_num, regions[r]->strand);
		regions[r]->get_reads(&bam_file, 1);
		regions[r]->compute_coverage();
		regions[r]->compute_intron_list();
		regions[r]->compute_intron_coverage();

		infer_genes(regions[r], &genes);
		delete regions[r];
	}
	printf("found %i genes\n", (int) genes.size()); 

	regions.clear();

	FILE* gff_fd = fopen(gff_file, "w"); 
	FILE* reg_fd = fopen(reg_file, "w"); 

	for (int r=0; r<genes.size(); r++)
	{
		genes[r]->find_orf(300, 0.7);
		genes[r]->print_gff3(gff_fd, r+1);
		genes[r]->print_region(reg_fd);
		delete genes[r];
	}
	fclose(gff_fd);
	fclose(reg_fd);
	genes.clear();
}
}
