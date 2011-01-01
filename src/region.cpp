#include <algorithm>
#include "assert.h"
#include "region.h"
#include "genome.h"
#include "get_reads_direct.h"
#include "read.h"
#include <vector>
  using std::vector;

bool compare_second(segment* intr1, segment* intr2);

/** default constructor*/
Region::Region()
{
	start = NULL; 
	stop = NULL;
	strand = NULL;
	chr_num = NULL;
	seq = NULL;
	coverage = NULL;
	intron_coverage = NULL;
	
}

/** constructor*/
Region::Region(int pstart, int pstop, int pchr_num, char pstrand, const char* gio_fname)
{
	start = pstart; 
	stop = pstop;
	strand = pstrand;
	chr_num = pchr_num;
	seq = NULL;
	coverage = NULL;
	intron_coverage = NULL;

	// initialize genome information object
	gio = new Genome(); 
	int ret = gio->init_genome((char*) gio_fname);
	if (ret<0)
	{   
		fprintf(stderr, "error reading genome info file: %s\n", gio_fname);
		return;
	}
}
/** constructor*/
Region::Region(int pstart, int pstop, int pchr_num, char pstrand)
{
	start = pstart; 
	stop = pstop;
	strand = pstrand;
	chr_num = pchr_num;
	seq = NULL;
	coverage = NULL;
	intron_coverage = NULL;
}

/** destructor*/
Region::~Region()
{
	delete[] seq;	
	for (int i=0; i<all_reads.size(); i++)
		delete all_reads[i];
	all_reads.clear();
	
	// reads is a subset of all_reads, therefore the 
	// destructor for each read has already been called
	// only the pointers have to be removed
	reads.clear(); 

	delete[] coverage;
	delete[] intron_coverage;
	intron_list.clear();
	unique_introns.clear();
	intron_counts.clear();
}

void Region::set_gio(Genome* pgio)
{
	gio = pgio;
}

void Region::load_genomic_sequence()
{
	seq = gio->read_flat_file(chr_num, start, stop);
}

void Region::print(_IO_FILE*& fd)
{
	fprintf(fd, "region start:\t%i\n", start);
	fprintf(fd, "region stop:\t%i\n", stop);
	fprintf(fd, "region strand:\t%c\n", strand);
	fprintf(fd, "region chr_num:\t%i\n", chr_num);
	if (gio)
		fprintf(fd, "region chr:\t%s\n", gio->get_contig_name(chr_num));
}

char* Region::get_region_str()
{
	char* reg_str = new char[1000];
	if (!gio)
	{
		fprintf(stderr, "genome information object not set");
		exit(-1);
	}
	sprintf(reg_str, "%s:%i-%i", gio->get_contig_name(chr_num), start, stop);
	return reg_str;
}

void Region::get_reads(char** bam_files, int num_bam_files)
{
	char* reg_str = get_region_str();
	int subsample = 1000;
	for (int i=0; i<num_bam_files; i++)
	{
    	char* fn_bam = bam_files[i];
	    fprintf(stdout, "getting reads from file: %s\n", fn_bam);
		get_reads_from_bam(fn_bam, reg_str, &all_reads, strand, subsample);
	}
	delete[] reg_str; 
	fprintf(stdout, "number of reads (not filtered): %d\n", (int) all_reads.size());
	/* filter reads
	* **************/
	int exon_len_filter = 8;
	int filter_mismatch = 1;
	int intron_len_filter = 100000;
	for (int i=0; i<all_reads.size(); i++)
	{
	    if (all_reads[i]->max_intron_len()<intron_len_filter && all_reads[i]->min_exon_len()>exon_len_filter && all_reads[i]->get_mismatches()<=filter_mismatch && all_reads[i]->multiple_alignment_index==0)
	        reads.push_back(all_reads[i]);
	}
	fprintf(stdout, "number of reads: %d\n", (int) reads.size());
}

void Region::compute_coverage() 
{
	int num_pos = stop-start+1;
	if (!coverage)
		coverage = new uint32_t[num_pos];

	for (int i=0; i<num_pos; i++)
		coverage[i] = 0;
	for (int i=0; i<reads.size(); i++)
	{
		// add read contribution to coverage
		reads[i]->get_coverage(start, stop, coverage);
	}
}
void Region::compute_intron_coverage() 
{
	int num_pos = stop-start+1;
	if (!intron_coverage)
		intron_coverage = new uint32_t[num_pos];

	for (int i=0; i<num_pos; i++)
		intron_coverage[i] = 0;
	for (int i=1; i<unique_introns.size(); i++)
	{
		int from_pos = std::max(start, unique_introns[i].first);
        int to_pos = std::min(stop, unique_introns[i].second);
		for (int j=from_pos; j<to_pos; j++)
		{
			intron_coverage[j] += intron_counts[i];
		}
	}
}

void Region::compute_intron_list()
{
	vector<int> introns;
	for (int i=0; i<reads.size(); i++)
	{
		reads[i]->get_introns(&introns);
    }
	for (int i=0; i<introns.size(); i+=2)
	{
		segment intr(introns[i], introns[i+1]);
		intron_list.push_back(intr);
	}
	printf("found %i introns\n", (int) intron_list.size());
	// sort by intron start
	sort(intron_list.begin(), intron_list.end());

	vector<segment*> same_start;
	for (int i=0; i<intron_list.size(); i++)
	{
		if (i<intron_list.size()-1 && intron_list[i].first==intron_list[i+1].first)
		{
			// collect introns with same start
			same_start.push_back(&intron_list[i]);
		}
		else
		{
			same_start.push_back(&intron_list[i]);
			sort(same_start.begin(), same_start.end(), compare_second);
			unique_introns.push_back(*(same_start[0]));
			intron_counts.push_back(1);
			for (int j=1; j<same_start.size(); j++)
			{
				if (same_start[j]->second!=same_start[j-1]->second)
				{
					unique_introns.push_back(*(same_start[j]));
					intron_counts.push_back(1);
				}
				else
				{
					//printf("intron_counts before: %i ", (*(intron_counts.back())));
					(intron_counts.back())++; 
					//printf("after: %i \n", (*intron_counts.back()));
				}
			}
			same_start.clear();
		}
	}
	printf("found %i unique introns\n", (int) unique_introns.size());

	int sum = 0;
	for (int i=0; i<intron_counts.size(); i++)
		sum+=intron_counts[i]; 
	assert(sum==intron_list.size());
}

bool compare_second(segment* intr1, segment* intr2)
{
	 return (intr1->second<intr2->second);
}
