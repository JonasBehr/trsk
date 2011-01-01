#include <string>
	using std::string;
#include <vector>
	using std::vector;
#include <math.h>
#include <assert.h>
#include "region.h"
#include "genome.h"
#include "gene.h"
#include "gene_tools.h"
#include "infer_genes.h"



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
	return max_idx;
}
int InferGenes::score_cand_intron(int idx, Region* region)
{
	segment intron = region->unique_introns[idx];

	if (region->intron_counts[idx]<3)
		return -1;

	int mean_cov = 0; 
	for (int i=intron.first; i<intron.second; i++)
		mean_cov+=region->coverage[i];
	mean_cov = mean_cov/(intron.second-intron.first+1);

	if (mean_cov>region->intron_counts[idx])
		return -1;//intron seems to be minor isoform

	return 1;
}
int InferGenes::score_cand_exon(segment exon, Region* region)
{
	int exon_len = exon.second-exon.first+1; 
	int min_len = 10;
	if (exon_len<min_len)
		return 0;
	int mean_cov = 0;
	int min_cov  = 1000; 
	int zero_cov = 0;
	int exon_cut = 5;
	int low_cov  = 0;
	for (int i=exon.first; i<exon.second; i++)
	{
		mean_cov+=region->coverage[i];
		if (region->coverage[i]<min_cov)
			min_cov = region->coverage[i];
		if (region->coverage[i]==0)
			zero_cov++;
		if (region->coverage[i]<exon_cut)
			low_cov++;
	}
	mean_cov = mean_cov/exon_len;

	//printf("mean_cov: %i, min_cov: %i zero_frac: %f\n", mean_cov, min_cov, ((double)zero_cov)/exon_len);
	if (((double)zero_cov)/exon_len>0.5)
	{
		//looks like an intergenic region 
		return -1;
	}
	if (mean_cov>5&&((double)low_cov)/exon_len<0.01)
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
	int threshold = 3;
	int min_len = 30;
    int win = 20;
	int j;
	float start_cov = mean(region->coverage, exon.first, exon.first+win);
	for (j=exon.first; j<exon.second-win; j+=10)//jump to speed up computation
	{
        float win_cov = mean(region->coverage, j, j+win);
        if (win_cov<threshold)
            break;

        if (win_cov<start_cov/5) //drop in coverage
            break;
	}
	if (j >= exon.second-win-1) // stop conditions never met
		exon.second = -1;
	else
		exon.second = std::min(j+3*win, exon.second);

	float mean_cov = mean(region->coverage, exon.first, exon.second);
	//printf("terminal exon: mean cov: %f\n", mean_cov);
	if (exon.second-exon.first<min_len || mean_cov<2*threshold)
		exon.second = -1;
	return exon;
}
segment InferGenes::find_initial_exon(segment exon, Region* region)
{
	//printf("initial exon[%i, %i], cov: %i\n", exon.first, exon.second, (int) region->coverage[exon.second]);
	assert(exon.first>=0);
	assert(exon.second<region->stop);
	int orig_start = exon.first;
	int threshold = 3;
	int min_len = 30; 
    int win = 20;
	int j;
	float start_cov = mean(region->coverage, exon.second-win, exon.second);
	for (j=exon.second-win; j>exon.first; j-=10)
	{
        float win_cov = mean(region->coverage, j, j+win);
        if (win_cov<threshold)
            break;

        if (win_cov<start_cov/5) //drop in coverage
            break;
	}
	if (j <= exon.first+win) // stop conditions never met
		exon.first = -1;
	else
		exon.first = std::max(j-3*win, exon.first);

	float mean_cov = mean(region->coverage, exon.first, exon.second);
	float mean_intron_cov = mean(region->intron_coverage, orig_start, exon.second);
	//printf("initial exon[%i, %i]: mean cov: %f\n\n", exon.first, exon.second, mean_cov);
	if (exon.second-exon.first<min_len || mean_cov<2*threshold)
		exon.first = -1;
	if (mean_intron_cov>2)//gene is likely to be cut
		exon.first = -1;
	return exon;
}

vector<segment>* InferGenes::greedy_extend(int intron_idx, Region* region, bool* intron_used)
{
	int cur_idx = intron_idx; 
	int max_exon_len = 8000; 
	segment cur_intron = region->unique_introns[cur_idx];

	vector<segment>* exons = new vector<segment>(); 

	segment initial_exon(std::max(0, cur_intron.first-max_exon_len), cur_intron.first-1);
	initial_exon = find_initial_exon(initial_exon, region);

	if (initial_exon.first==-1)
		return exons; 
	exons->push_back(initial_exon);
	// extend transcript downstream
	while (true)
	{
		int next_idx = find_next_intron(cur_idx, region);
		segment next_intron;
		if (next_idx>-1)
			next_intron = region->unique_introns[next_idx];
		else
		{
	//		printf("no intron found\n");
	//		segment exon(cur_intron.second+1, std::min(region->stop, cur_intron.second+max_exon_len));
	//		exon = find_terminal_exon(exon, region);
	//		if (exon.second==-1)
	//		{
	//			exons->clear();
	//			return exons;
	//		}
	//		exons->push_back(exon);
	//		break;

			exons->clear();
			return exons;
		}
		intron_used[next_idx] = true;
		segment exon(cur_intron.second+1, next_intron.first-1);

		if (score_cand_intron(next_idx, region)==-1 && exon.second-exon.first<max_exon_len)
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
		if (exon.second-exon.first>max_exon_len)
		{
			exon = find_terminal_exon(exon, region);
			if (exon.second==-1)
			{
				exons->clear();
				return exons;
			}
			exons->push_back(exon);
			break;
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
	float relative_length = 0.5;
	int threshold = 3;
	int min_len = 50; 
    int win = 100;
	bool strand_specific = false;
	int j;
	for (j=gene->start-win-min_len; j>0; j-=10)
	{
        float win_cov = mean(region->intron_coverage, j, j+win);
        if (strand_specific)
			win_cov += mean(region->coverage, j, j+win);
        
		if (win_cov>threshold)
            break;
	}
	int region_start = j+win;
	gene->intergenic_region_start = (int) ceil((1-relative_length)*gene->start+relative_length*region_start);
	// find end
	for (j=gene->stop+min_len; j<region->stop-win-1; j+=10)
	{
        float win_cov = mean(region->intron_coverage, j, j+win);
        if (strand_specific)
			win_cov += mean(region->coverage, j, j+win);
        
		if (win_cov>threshold)
            break;
	}
	gene->intergenic_region_stop = (int) ceil((1-relative_length)*gene->stop+relative_length*j);

	//printf("reg: %i int_start: %i, start:%i, stop: %i, int_stop: %i\n", region_start, gene->intergenic_region_start, gene->start, gene->stop, gene->intergenic_region_stop);
}

void InferGenes::infer_genes(Region* region, vector<Gene*>* genes)
{
	bool intron_used[region->unique_introns.size()];
	for (int i=0; i<region->unique_introns.size(); i++)
		intron_used[i] = false;

	//int median_cnt = median(region->intron_counts);
	int median_cnt = 3;

	for (int i=0; i<region->unique_introns.size(); i++)
	{
		if (intron_used[i])
			continue;
		
		if (region->intron_counts[i]<median_cnt)
			continue;

		intron_used[i] = true;
		
		//printf("starting greedy extend \n");
		vector<segment>* exons = greedy_extend(i, region, intron_used);
		//printf("greedy extend found %i exons \n", (int) exons->size());

		if (exons->size()>0)
		{
			Gene* g = new Gene(exons, region->chr_num, region->strand, region->gio);
			find_intergenic_region(region, g);
			genes->push_back(g);
		}
	}
}

int InferGenes::run_infer_genes(char* gio_file, char* bam_file, char* gff_file, char* reg_file)
{   
	vector<Region*> regions = GeneTools::init_regions(gio_file);
	printf("regions.size(): %i\n", (int) regions.size());

	//regions[0]->stop = 1e6;
	//regions[1]->stop = 1e6;
	vector<Gene*> genes;
	for (int r=0; r<regions.size(); r++)
	//for (int r=1; r<2; r++)
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

	//genes[0]->find_orf();
	//genes[0]->print_gff3(stdout, 0);
	//printf("\n\n\n"); 
	//genes[0]->strand='+';
	//genes[0]->find_orf();
	//genes[0]->print_gff3(stdout, 0);
	//printf("\n\n\n"); 
	//genes[1]->find_orf();
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
