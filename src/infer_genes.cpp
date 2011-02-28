#include <string>
	using std::string;
#include <vector>
	using std::vector;
#include <math.h>
#include <algorithm>
#include <assert.h>
#include "region.h"
#include "genome.h"
#include "gene.h"
#include "gene_tools.h"
#include "infer_genes.h"
#include "Config.h"

//int median(vector<int>* vec)
//{
//	vector<int> vec2 ///copy
//	sort
//	take elem
//}

int InferGenes::find_next_intron(int idx, Region* region)
{
	segment start_intron = region->unique_introns[idx];
	// unique introns sorted by start
	int cidx = idx+1;
	segment cintron = region->unique_introns[cidx];
	while(cidx<region->unique_introns.size())
	{
		//skip over introns that have the same start
		cintron = region->unique_introns[cidx];
		if (cintron.first!=start_intron.first)
			break;
		cidx++; 
	}
	int next_start = cintron.first;
	int max_cov = 0;
	int max_idx = -1;
	while(cidx<region->unique_introns.size())
	{
		// find next intron with max coverage
		// starting at the next downstream position
		cintron = region->unique_introns[cidx];
		if (cintron.first!=next_start)
			break;
		if (region->intron_counts[cidx]>max_cov)
		{
			max_cov = region->intron_counts[cidx];
			max_idx = cidx; 
		}
		cidx++;
	}
	while(cidx<region->unique_introns.size())
	{
		// find another intron with max coverage
		// starting at a nearby position
		cintron = region->unique_introns[cidx];
		if (cintron.first>next_start+conf->intron_dist)
			break;
		if (region->intron_counts[cidx]>max_cov)
		{
			max_cov = region->intron_counts[cidx];
			max_idx = cidx; 
		}
		cidx++;
	}

	return max_idx;
}
/** find intron with maximal end position to the left of the 
 *  intron handed over 
 */
int InferGenes::find_previous_intron(int idx, Region* region, bool verb)
{
	if (idx==0)
		return -1;
	segment start_intron = region->unique_introns[idx];
	if (verb)
		printf("start: [%i, %i]\n", start_intron.first, start_intron.second);
	// unique introns sorted by start
	int cidx = idx-1;
	int max_len = conf->max_exon_len+conf->max_intron_len; 
	int max_second = 0;
	int max_cov = 0;
	int best_idx = -1;
	segment cintron = region->unique_introns[cidx];
	while(cidx>0 && cintron.first>start_intron.first-max_len)
	{
		cintron = region->unique_introns[cidx];
		if (verb)
			printf("cintron: [%i, %i]\n", cintron.first, cintron.second);
		if (cintron.second>=start_intron.first)
		{
			//skip overlapping intron
			if (verb)
				printf("skip: [%i, %i]\n", cintron.first, cintron.second);
		}
		else if (cintron.second>=max_second+conf->intron_dist || 
				(cintron.second>=max_second-conf->intron_dist && region->intron_counts[cidx]>max_cov))
		{
			max_second = cintron.second;
			max_cov = region->intron_counts[cidx];
			best_idx = cidx;
			if (verb)
				printf("better: [%i, %i]\n", cintron.first, cintron.second);
		}
		else if (verb)
		{
				printf("ignore: [%i, %i]\n", cintron.first, cintron.second);
		}
		cidx--; 
	}
	return best_idx;
}
int InferGenes::score_cand_intron(int idx, Region* region)
{
	segment intron = region->unique_introns[idx];

	if (region->intron_counts[idx]<conf->intron_conf)
		return -1;

	if (conf->reject_retained_introns)
	{
		float mean_cov = mean(region->coverage, intron.first, intron.second); 
		if (mean_cov>region->intron_counts[idx])
			return -1;//intron seems to be minor isoform
	}
	return 1;
}
int InferGenes::score_cand_exon(segment exon, Region* region)
{
	int exon_len = exon.second-exon.first+1; 
	if (exon_len<conf->min_exon_len)
		return 0; // discard transcript
	if (exon_len>conf->max_exon_len)
		return -1; // search for end of gene
	int mean_cov = 0;
	int min_cov  = 1000;// any large number... this does not matter 
	int zero_cov = 0;
	int low_cov  = 0;
	for (int i=exon.first; i<exon.second; i++)
	{
		mean_cov+=region->coverage[i];
		if (region->coverage[i]<min_cov)
			min_cov = region->coverage[i];
		if (region->coverage[i]==0)
			zero_cov++;
		if (region->coverage[i]<conf->exon_cut)
			low_cov++;
	}
	mean_cov = mean_cov/exon_len;

	//printf("mean_cov: %i, min_cov: %i zero_frac: %f\n", mean_cov, min_cov, ((double)zero_cov)/exon_len);
	if (((double)zero_cov)/exon_len>0.2)
	{
		//looks like an intergenic region 
		return -1;
	}
	if (mean_cov>conf->exon_mean&&((double)low_cov)/exon_len<0.01)
	{
		// the more than 99% of all nucleotides are covered 
		// higher than some threshold (exon_cut)
		//looks like a true exon
		return 1;
	}
	// undecided
	return 0;
}
float InferGenes::mean(uint32_t* arr, int from, int to)
{
	//printf("initial exon[%i, %i], cov: %i\n", from, to, (int) arr[to]);
	float mean_val = 0; 
	for (int i=from; i<=to; i++)
		mean_val+=(float) arr[i];
	//printf("initial exon[%i, %i], sum: %f\n", from, to, mean_val);
	mean_val = mean_val/(to-from+1);
	return mean_val;
}
segment InferGenes::find_terminal_exon(segment exon, Region* region)
{
	int orig_stop = exon.second;
	int threshold = conf->exon_term_thresh;
    int win = 20;
	int j;
	float start_cov = mean(region->coverage, exon.first, exon.first+win);
	for (j=exon.first; j<exon.second-win; j+=10)//jump to speed up computation
	{
        float win_cov = mean(region->coverage, j, j+win);
        if (win_cov<threshold)
            break;

        if (win_cov<start_cov/conf->exon_drop) //drop in coverage
            break;
	}
	if (j >= exon.second-win-1) // stop conditions never met
	{
		exon.second = -1;
		return exon;
	}
	else
		exon.second = std::min(j+conf->term_offset, exon.second);

	float mean_cov = mean(region->coverage, exon.first, j);
	//printf("terminal exon: mean cov: %f\n", mean_cov);
	float mean_intron_cov = mean(region->intron_coverage, exon.second, std::min(exon.second+100, orig_stop));
	if (exon.second-exon.first<conf->min_exon_len || mean_cov<conf->term_filter*threshold)
	{
		exon.second = -1;
		term_reject_cov++;
	}
	if (mean_intron_cov>conf->intron_cut)//gene is likely to be cut
	{
		exon.second = -1;
		term_reject_intron++;
	}

	return exon;
}
segment InferGenes::find_initial_exon(segment exon, Region* region)
{
	//printf("initial exon[%i, %i], cov: %i\n", exon.first, exon.second, (int) region->coverage[exon.second]);
	assert(exon.first>=0);
	assert(exon.second<region->stop);
	int orig_start = exon.first;
	int threshold = conf->exon_term_thresh;
    int win = 20;
	int j;
	float start_cov = mean(region->coverage, exon.second-win, exon.second);
	for (j=exon.second-win; j>exon.first; j-=10)
	{
        float win_cov = mean(region->coverage, j, j+win);
        if (win_cov<threshold)
            break;

        if (win_cov<start_cov/conf->exon_drop) //drop in coverage
            break;
	}
	if (j <= exon.first+win) // stop conditions never met
	{
		exon.first = -1;
		return exon;
	}
	else
		exon.first = std::max(j-conf->term_offset, exon.first);

	float mean_cov = mean(region->coverage, j, exon.second);
	float mean_intron_cov = mean(region->intron_coverage, std::max(0, exon.first-100), exon.first);
	//printf("initial exon[%i, %i]: mean cov: %f\n\n", exon.first, exon.second, mean_cov);
	if (exon.second-exon.first<conf->min_exon_len || mean_cov<conf->term_filter*threshold)
	{
		exon.first = -1;
		init_reject_cov++;
	}
	if (mean_intron_cov>conf->intron_cut)//gene is likely to be cut
	{
		exon.first = -1;
		init_reject_intron++;
	}
	return exon;
}
bool InferGenes::check_segment(segment seg, Region* region)
{
	if (seg.first<region->start || seg.first>region->stop)
		return false;
	if (seg.second<region->start || seg.second>region->stop)
		return false;
	return true;
}
vector<segment>* InferGenes::greedy_extend(int intron_idx, Region* region, bool* intron_used)
{
	int cur_idx = intron_idx; 
	segment cur_intron = region->unique_introns[cur_idx];

	vector<segment>* exons = new vector<segment>(); 

	if (score_cand_intron(intron_idx, region)==-1)
		return exons;

	// extent upstream
	while (true)
	{
		int prev_idx = find_previous_intron(cur_idx, region, false);
		//printf("extend upstream: cur_idx: %i, prev_idx: %i\n", cur_idx, prev_idx);
		segment prev_intron;
		if (prev_idx==-1)
		{
			no_upstream_intron++;	
			//printf("no upstream intron found: %i [%i %i]\n", no_upstream_intron, cur_intron.first, cur_intron.second);
			// cannot find previous intron => search for initial exon
			segment initial_exon(std::max(0, cur_intron.first-conf->max_exon_len), cur_intron.first-1);
			initial_exon = find_initial_exon(initial_exon, region);
			if (initial_exon.first==-1)
			{
				no_upstream_intron_reject++;	
				exons->clear();
				return exons; 
			}
			exons->push_back(initial_exon);
			break;
		}
		else
		{
			prev_intron = region->unique_introns[prev_idx];
		}
		intron_used[prev_idx] = true;
		segment exon(prev_intron.second+1, cur_intron.first-1);
		if (score_cand_intron(prev_idx, region)==-1 && exon.second-exon.first<conf->max_exon_len)
		{
			//todo: differentiate between introns with low coverage and introns 
			//that are often skipped
			exons->clear();
			return exons;
		}
		int exon_score = score_cand_exon(exon, region);
		if (exon_score==1)
		{
			// accept current exon and move on to next intron
			exons->push_back(exon); 
			cur_idx = prev_idx;
			cur_intron = region->unique_introns[cur_idx];
		}
		else if (exon_score==0)
		{
			// unclear => reject transcript
			exons->clear();
			return exons;
		}
		else
		{
			exon = find_initial_exon(exon, region); 
			if (exon.first==-1)
			{
				exons->clear();
				return exons;
			}
			exons->push_back(exon);
			break; 
		}
	}	

	// extend transcript downstream
	cur_idx = intron_idx; 
	cur_intron = region->unique_introns[cur_idx];
	while (true)
	{
		int next_idx = find_next_intron(cur_idx, region);
		segment next_intron;
		if (next_idx>-1)
			next_intron = region->unique_introns[next_idx];
		else
		{
			no_downstream_intron++;
			//printf("no downstream intron found: %i [%i %i]\n", no_downstream_intron, cur_intron.first, cur_intron.second);
			segment exon(cur_intron.second+1, std::min(region->stop, cur_intron.second+conf->max_exon_len));
			exon = find_terminal_exon(exon, region);
			if (exon.second==-1)
			{
				no_downstream_intron_reject++;
				exons->clear();
				return exons;
			}
			exons->push_back(exon);
			break;
		}
		intron_used[next_idx] = true;
		segment exon(cur_intron.second+1, next_intron.first-1);

		if (score_cand_intron(next_idx, region)==-1 && exon.second-exon.first<conf->max_exon_len)
		{
			// take next intron only into account if it 
			// is not to far away. Otherwise it is more likely 
			// that we are looking at the first intron of the next gene. 
			// Still one could try to recover some of the cases where 
			// another downstream intron would fit the criteria

			//if (score_cand_intron(next_idx, region)==-1)
			//	printf("intron score -1\n");
			exons->clear();
			return exons;
		}
		int exon_score = score_cand_exon(exon, region);
		if (exon_score==1)
		{
			// accept current exon and move on to next intron
			exons->push_back(exon); 
			cur_idx = next_idx;
			cur_intron = region->unique_introns[cur_idx];
		}
		else if (exon_score==0)
		{
			// unclear => reject transcript
			exons->clear();
			return exons;
		}
		else
		{
			exon = find_terminal_exon(exon, region); 
			if (exon.second==-1)
			{
				exons->clear();
				return exons;
			}
			exons->push_back(exon);
			//printf("terminal exon\n");
			break; 
		}	
	}
	return exons;
}

void InferGenes::find_intergenic_region(Region* region, Gene* gene)
{
	// find start
	float relative_length = (float) conf->region_rel_length;
	int threshold = 3;
	int min_intergenic_len = conf->min_intergenic_len; 
    int win = conf->intergenic_win;
	bool strand_specific = conf->strand_specific;
	int j;
	for (j=gene->start-win-min_intergenic_len; j>0; j-=10)
	{
        float win_cov = mean(region->intron_coverage, j, j+win);
        if (strand_specific)
			win_cov += mean(region->coverage, j, j+win);
        
		if (win_cov>threshold)
            break;
	}
	int region_start = j+win;
	gene->intergenic_region_start = (int) ceil((1-relative_length)*gene->start+relative_length*region_start);
	if (gene->intergenic_region_start<gene->start-conf->max_intergenic_len)
		gene->intergenic_region_start = gene->start-conf->max_intergenic_len;
	// find end
	for (j=gene->stop+min_intergenic_len; j<region->stop-win-1; j+=10)
	{
        float win_cov = mean(region->intron_coverage, j, j+win);
        if (strand_specific)
			win_cov += mean(region->coverage, j, j+win);
        
		if (win_cov>threshold)
            break;
	}
	gene->intergenic_region_stop = (int) ceil((1-relative_length)*gene->stop+relative_length*j);
	if (gene->intergenic_region_stop>gene->stop+conf->max_intergenic_len)
		gene->intergenic_region_stop = gene->stop+conf->max_intergenic_len;
		
	//printf("reg: %i int_start: %i, start:%i, stop: %i, int_stop: %i\n", region_start, gene->intergenic_region_start, gene->start, gene->stop, gene->intergenic_region_stop);
}

void InferGenes::infer_genes(Region* region, vector<Gene*>* genes)
{
	bool intron_used[region->unique_introns.size()];
	for (int i=0; i<region->unique_introns.size(); i++)
		intron_used[i] = false;

	//int median_cnt = median(region->intron_counts);
	//int median_cnt = 3;

	for (int i=0; i<region->unique_introns.size(); i++)
	{
		if (intron_used[i])
			continue;
		
		if (region->intron_counts[i]<conf->intron_seed_conf)
			continue;

		intron_used[i] = true;
		
		//printf("starting greedy extend \n");
		vector<segment>* exons = greedy_extend(i, region, intron_used);
		//printf("greedy extend found %i exons \n", (int) exons->size());

		if (exons->size()>0)
		{
			// sort exons by start position
			sort(exons->begin(), exons->end());
			Gene* g = new Gene(exons, region->chr_num, region->strand, region->gio);
			find_intergenic_region(region, g);
			genes->push_back(g);
		}
	}
	if (conf->find_single_exon_genes_orf)
	{
		segment gseq(region->start, region->stop);
		vector<segment> exons;
		exons.push_back(gseq);
		Gene* dummy_gene = new Gene(&exons, region->chr_num, region->strand, region->gio);
		char* seq;
		int len;
		dummy_gene->get_mRNA_seq(&seq, &len);
		vector<int> tis; 
		vector<int> stop;
		GeneTools::find_all_orfs(seq, len,  &tis, &stop, 1000);
		for (int i=0; i<tis.size(); i++)
		{
			tis[i] = dummy_gene->map_rna_to_dna(tis[i]);
			stop[i] = dummy_gene->map_rna_to_dna(stop[i]);
		}
		printf("found %i orfs longer than 1000\n", (int) tis.size()); 

		for (int i=0; i<tis.size(); i++)
		{
			// define to be strand independent
			int cds_start = std::min(tis[i], stop[i]);
			int cds_stop = std::max(tis[i], stop[i]);

			// make sure no introns are in that regions
			int from = std::max(region->start, cds_start-1000);
			int to = std::min(region->stop, cds_stop+1000);
			float mean_cov = mean(region->coverage, cds_start, cds_stop);
			float mean_intron_cov = mean(region->intron_coverage, from, to);
			if (mean_intron_cov<1e-5 && mean_cov>10)
			{
				//find tss and cleave
				segment exon(std::max(0, cds_start-conf->max_exon_len), cds_start-1);
				exon = find_initial_exon(exon, region);
				if (exon.first==-1)
					continue;
				int first = exon.first;
				exon.first = cds_stop+1;
				exon.second = cds_start+conf->max_exon_len;
				exon = find_terminal_exon(exon, region);
				if (exon.second==-1)
					continue;
				exon.first = first;

				exons.clear();
				exons.push_back(exon);
				Gene* g = new Gene(&exons, region->chr_num, region->strand, region->gio);
				find_intergenic_region(region, g);
				genes->push_back(g);
			}
				
		}
	}
}

void InferGenes::process_gene(Gene* gene)
{
	// set the exon boundaries to the original place
	// (these have been moved before to find more ORFs)
	int minlen = conf->min_exon_len;

	int end = std::max(gene->exons.back().first+minlen, gene->exons.back().second-conf->term_offset);
	if (gene->is_coding())
		end = std::max(end, gene->cds_exons.back().second+minlen);
	gene->exons.back().second = end;

	int start = std::min(gene->exons.front().second-minlen, gene->exons.front().first+conf->term_offset);
	if (gene->is_coding())
		start = std::min(start, gene->cds_exons.front().first-minlen);
	gene->exons.front().first = start;

	//update UTR exons
	if (gene->is_coding() && gene->strand=='+')
		gene->split_exons(gene->cds_exons.front().first-1, gene->cds_exons.back().second-4);
	else if (gene->is_coding())
		gene->split_exons(gene->cds_exons.back().second-4, gene->cds_exons.front().first-1);
}

int InferGenes::run_infer_genes()
{   
	if (!conf->have_bam_file)
	{
		fprintf(stderr, "Error: no bam file given (-bam)\n");
		return -1;
	}
	if (!conf->have_gff_file)
	{
		fprintf(stderr, "Error: no gff output file given (-gff)\n");
		return -1;
	}
	if (!conf->have_gio_file)
	{
		fprintf(stderr, "Error: no genome informatio object given (-gio)\n");
		return -1;
	}
	// statistics
	no_upstream_intron = 0;	
	no_upstream_intron_reject = 0;	
	no_downstream_intron = 0;	
	no_downstream_intron_reject = 0;	
	init_reject_intron = 0;
	term_reject_intron = 0;
	init_reject_cov = 0;
	term_reject_cov = 0;
	//

	vector<Region*> regions = GeneTools::init_regions(conf->gio_file);
	printf("regions.size(): %i\n", (int) regions.size());

	conf->print(stdout);

	vector<Gene*> genes;
	//for (int r=0; r<2; r++)
	for (int r=0; r<regions.size(); r++)
	{
		printf("Starting with contig %s, strand %c\n", regions[r]->get_region_str(), regions[r]->strand);
		regions[r]->get_reads(&(conf->bam_files[0]), conf->bam_files.size(), conf->max_intron_len, conf->mm_filter, conf->el_filter);
		regions[r]->compute_coverage();
		regions[r]->compute_intron_list();
		regions[r]->compute_intron_coverage();

		infer_genes(regions[r], &genes);
		printf("found %i genes\n", (int) genes.size());

        // write seqs
        //...

		delete regions[r];
	}

	regions.clear();

	if (conf->have_reg_file)
	{
		FILE* reg_fd = fopen(conf->reg_file, "w"); 
		for (int r=0; r<genes.size(); r++)
			genes[r]->print_region(reg_fd);
		fclose(reg_fd);
	}

	FILE* gff_fd = fopen(conf->gff_file, "w"); 
    //TODO set from config file
	FILE* tis_fd = fopen("/fml/ag-raetsch/home/cwidmer/svn/projects/transcript_skimmer/src/tis.flat", "w"); 
	int coding_cnt = 0;
	int non_coding_cnt = 0;
	int single = 0;
    int num_tis_labels = 0;

	for (int r=0; r<genes.size(); r++)
	{
		if (conf->find_orf)
			genes[r]->find_orf(300, 0.7);

		if (false)//move transcript start and stop back to the original place
			process_gene(genes[r]);

		if (genes[r]->is_coding())
			coding_cnt++;
		else
			non_coding_cnt++;
		if (genes[r]->exons.size()==1)
			single++;
		genes[r]->print_gff3(gff_fd, r+1);

        // create tis label, increment counter if successful
        if (genes[r]->generate_tis_labels(tis_fd))
            num_tis_labels++;

		delete genes[r];
	}
	fclose(gff_fd);
    
    /*
	if (true)
	{
		//FILE* reg_fd = fopen(conf->reg_file, "w"); 
		for (int r=0; r<genes.size(); r++)
			genes[r]->generate_tis_labels(stdin);
		//fclose(reg_fd);
	}
    */
	printf("statistics:\n"); 
	printf("\tfound %i genes (%i coding, %i non-coding, %i single-exon)\n\n", (int) genes.size(), coding_cnt, non_coding_cnt, single);
	printf("\tno upstream intron: %i\n", no_upstream_intron);
	printf("\tno upstream intron reject: %i\n", no_upstream_intron_reject);
	printf("\tno downstream intron: %i\n", no_downstream_intron);
	printf("\tno downstream intron reject: %i\n", no_downstream_intron_reject);
	printf("\tinitial exon reject (cov): %i\n", init_reject_cov);
	printf("\tinitial exon reject (intron): %i\n", init_reject_intron);
	printf("\tterminal exon reject (cov): %i\n", term_reject_cov);
	printf("\tterminal exon reject (intron): %i\n", term_reject_intron);
	printf("\tno generated tis labels: %i\n", num_tis_labels);
	genes.clear();

}
