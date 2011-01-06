#ifndef _INFER_GENES_H__
#define _INFER_GENES_H__

#include <string>
	using std::string;
#include <vector>
	using std::vector;
#include "region.h"
#include "genome.h"
#include "gene.h"
#include "gene_tools.h"

class InferGenes
{
	public: 
	//	InferGenes(){conf = new Config{};}; 
		static void infer_genes(Region* region, vector<Gene*>* genes);
		static int run_infer_genes(char* gio_file, char* bam_file, char* gff_file, char* reg_gile);

	private:
		static int find_next_intron(int idx, Region* region);
		static int find_previous_intron(int idx, Region* region);
		static int score_cand_intron(int idx, Region* region);
		static int score_cand_exon(segment exon, Region* region);
		static segment find_terminal_exon(segment exon, Region* region);
		static segment find_initial_exon(segment exon, Region* region);
		static void find_intergenic_region(Region* region, Gene* gene);
		static vector<segment>* greedy_extend(int intron_idx, Region* region, bool* intron_used);
		static float mean(uint32_t* arr, int from, int to);
		
};
#endif
