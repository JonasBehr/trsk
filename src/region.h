
#ifndef _REGION_H__
#define _REGION_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "genome.h"
#include "read.h"
#include <utility>
	using std::pair;
#include <vector>
	using std::vector;

typedef pair<int, int> segment;

class Region
{
	public:
		int start; 
		int stop;
		int chr_num;
		char strand; 
		char* seq;
		vector<CRead*> all_reads; 
		vector<CRead*> reads; 
		uint32_t* coverage;
		uint32_t* intron_coverage;
		vector<segment> intron_list;
		vector<segment> unique_introns;
		vector<int> intron_counts;

		Genome* gio; 
    	
			
		/** default constructor*/	
		Region();

		/** constructor*/
		Region(int pstart, int pstop, int pchr_num, char pstrand, const char* gio_fname);
		Region(int pstart, int pstop, int pchr_num, char pstrand);

		/** destructor*/
		~Region();

		void set_gio(Genome* pgio);
	
		void load_genomic_sequence(); 			

		void get_reads(char** bam_files, int num_bam_files);

		void compute_coverage();

		void compute_intron_coverage();

		void compute_intron_list();

		char* get_region_str();

		virtual void print(_IO_FILE*& fd);

};

#endif
