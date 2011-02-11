
#ifndef __GENE_TOOLS_
#define __GENE_TOOLS_

#include <string>
	using std::string;
#include <vector>
	using std::vector;
#include "region.h"
#include "genome.h"
#include "gene.h"

class GeneTools
{
	public:
		static bool check_consensus(int pos, const char* seq, int len, vector<string*> motifs);
		static void find_max_orf(char* seq, int len, int* tis, int* stop, int* second_best);
		static void find_all_orfs(char* seq, int len, vector<int>* tis, vector<int>* stop, int min_len);
		static vector<Region*> init_regions(const char* gio_fname);
};
#endif
