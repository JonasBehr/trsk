#include <string>
	using std::string;
#include <vector>
	using std::vector;
#include "region.h"
#include "genome.h"
#include "gene.h"
#include "gene_tools.h"
#include "infer_genes.h"


int main(int argc, char* argv[])
{   

	InferGenes* ig = new InferGenes(argc, argv);

	ig->run_infer_genes();
}
