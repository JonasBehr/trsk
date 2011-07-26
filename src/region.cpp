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
	chr = NULL;
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
	chr = NULL;
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
	chr = NULL;
	seq = NULL;
	coverage = NULL;
	intron_coverage = NULL;
}
/** constructor*/
Region::Region(int pstart, int pstop, char* pchr, char pstrand)
{
	start = pstart; 
	stop = pstop;
	strand = pstrand;
	chr_num = NULL;
	chr = pchr;
	seq = NULL;
	coverage = NULL;
	intron_coverage = NULL;
	fd_out = stdout;
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

	delete chr;
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
	if (!check_region())
	{
		printf("load_genomic_sequence: check_region failed\n");
		print(stderr);
		exit(-1);
	}
	seq = gio->read_flat_file(chr_num, start, stop);
}

void Region::print(_IO_FILE*& fd)
{
	fprintf(fd, "region %s\n", get_region_str());
	fprintf(fd, "region start:\t%i\n", start);
	fprintf(fd, "region stop:\t%i\n", stop);
	fprintf(fd, "region strand:\t%c\n", strand);
	if (chr_num)
		fprintf(fd, "region chr_num:\t%i\n", chr_num);
	if (gio && chr_num)
		fprintf(fd, "region chr:\t%s\n", gio->get_contig_name(chr_num));
	else if (chr)
		fprintf(fd, "region chr:\t%s\n", chr);
}

char* Region::get_region_str()
{
	char* reg_str = new char[1000];
	if (!gio && !chr)
	{
		fprintf(stderr, "genome information object not set");
		exit(-1);
	}
	if (gio && chr_num)
		sprintf(reg_str, "%s:%i-%i", gio->get_contig_name(chr_num), start, stop);
	else if (chr)
		sprintf(reg_str, "%s:%i-%i", chr, start, stop);
	return reg_str;
}

void Region::get_reads(char** bam_files, int num_bam_files, int intron_len_filter, int filter_mismatch, int exon_len_filter)
{
	char* reg_str = get_region_str();
	int subsample = 1000;
	for (int i=0; i<num_bam_files; i++)
	{
    	char* fn_bam = bam_files[i];
	    fprintf(fd_out, "getting reads from file: %s\n", fn_bam);
		get_reads_from_bam(fn_bam, reg_str, &all_reads, strand, subsample);
	}
	delete[] reg_str; 
	fprintf(fd_out, "number of reads (not filtered): %d\n", (int) all_reads.size());
	/* filter reads
	* **************/
	//int exon_len_filter = 8;
	//int filter_mismatch = 1;
	//int intron_len_filter = 100000;
	for (int i=0; i<all_reads.size(); i++)
	{
	    if (all_reads[i]->max_intron_len()<intron_len_filter && all_reads[i]->min_exon_len()>exon_len_filter && all_reads[i]->get_mismatches()<=filter_mismatch && all_reads[i]->multiple_alignment_index==0)
	        reads.push_back(all_reads[i]);
	}
	fprintf(fd_out, "number of reads: %d\n", (int) reads.size());
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

float Region::get_coverage_global(int pstart, int pstop)
{
	if (!coverage)
		compute_coverage();
	
	assert(pstart<pstop);

	int sum = 0;
	for (int i = pstart; i<pstop; i++)
	{
		assert(i>=start && i<stop);
		sum+=coverage[i-start];
	}
	return ((float) sum)/(pstop-pstart+1);
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

int Region::get_intron_conf(int intron_start, int intron_stop)
{
	//printf("get_intron_conf for %i -> %i\n", intron_start, intron_stop);
	for (int i=1; i<unique_introns.size(); i++)
	{
		//if (unique_introns[i].first==intron_start || unique_introns[i].second==intron_stop)
		//printf("found intron %i -> %i\n", unique_introns[i].first, unique_introns[i].second);
		if (unique_introns[i].first==intron_start && unique_introns[i].second==intron_stop)
			return intron_counts[i];
	}
	return 0;
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
	fprintf(fd_out, "found %i introns\n", (int) intron_list.size());
	
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
	fprintf(fd_out, "found %i unique introns\n", (int) unique_introns.size());

#ifdef READ_DEBUG
	//check unique list
	printf("DEBUG: Check 'unique introns'-list contains all introns\n");
	int idx = 0;
	for (int i=0; i<intron_list.size(); i++)
	{
		bool match = false;
		while (unique_introns[idx].first<intron_list[i].first && idx<unique_introns.size())
			idx++;
		int idx2 = idx;
		while (unique_introns[idx2].first==intron_list[i].first && idx2<unique_introns.size())
		{
			if (unique_introns[idx2].second==intron_list[i].second)
				match = true;
			idx2++;
		}
		assert(match);
	}
#endif

	int sum = 0;
	for (int i=0; i<intron_counts.size(); i++)
		sum+=intron_counts[i]; 
	assert(sum==intron_list.size());
}

bool compare_second(segment* intr1, segment* intr2)
{
	 return (intr1->second<intr2->second);
}
