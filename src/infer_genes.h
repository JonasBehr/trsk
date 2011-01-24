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
#include "Config.h"

class InferGenes
{
	public: 
		InferGenes(){conf = new Config();}; 
		InferGenes(int argc, char *argv[]){conf = new Config(argc, argv);}; 
		void infer_genes(Region* region, vector<Gene*>* genes);
		int run_infer_genes();
		Config* conf; 

	private:
		int find_next_intron(int idx, Region* region);
		int find_previous_intron(int idx, Region* region, bool verb);
		int score_cand_intron(int idx, Region* region);
		int score_cand_exon(segment exon, Region* region);
		segment find_terminal_exon(segment exon, Region* region);
		segment find_initial_exon(segment exon, Region* region);
		void find_intergenic_region(Region* region, Gene* gene);
		vector<segment>* greedy_extend(int intron_idx, Region* region, bool* intron_used);
		float mean(uint32_t* arr, int from, int to);
		bool check_segment(segment seg, Region* region);
		
		//statistics: 
		int no_upstream_intron; 
		int no_upstream_intron_reject; 
		int no_downstream_intron; 
		int no_downstream_intron_reject; 
		int init_reject_intron; 
		int term_reject_intron;
		int init_reject_cov;
		int term_reject_cov;
};
#endif
