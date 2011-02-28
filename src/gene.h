
#ifndef _GENE_H__
#define _GENE_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "genome.h"
#include "region.h"
#include "read.h"
#include <utility>
	using std::pair;
#include <vector>
	using std::vector;


class Gene: public Region
{
	public:
		Gene(int pstart, int pstop, int pchr_num, char pstrand, const char* gio_fname):Region(pstart, pstop, pchr_num, pstrand, gio_fname){};
		Gene(int pstart, int pstop, int pchr_num, char pstrand):Region(pstart, pstop, pchr_num, pstrand){};
		Gene(vector<segment>* pexons, int pchr_num, char pstrand, Genome* pgio):Region()
		{
			// assert that exons are sorted
			for (int i=1; i<exons.size(); i++)
			{
				assert(pexons->at(i-1).first <pexons->at(i-1).second);
				assert(pexons->at(i-1).second<pexons->at(i).first);
				assert(pexons->at(i).first   <pexons->at(i).second);
			}
			exons = *pexons; 
			start = pexons->at(0).first; 
			stop = pexons->back().second;
			chr_num = pchr_num; 
			strand = pstrand;
			gio = pgio; 
		};
		~Gene()
		{
			exons.clear();
			cds_exons.clear();	
			utr5exons.clear();	
			utr3exons.clear();	
		};
	
		vector<segment> exons;
		vector<segment> cds_exons;
		vector<segment> utr3exons;
		vector<segment> utr5exons;
		int intergenic_region_start;
		int intergenic_region_stop;

		bool is_coding(){return cds_exons.size()>0;};
	
		void find_orf(int min_len, float separation);
		void get_mRNA_seq(char** seq, int* len);
		int map_rna_to_dna(int rna_pos);
		void split_exons(int tis, int stop);

		void print_region(_IO_FILE*& fd)
		{
			char* chr = gio->get_contig_name(chr_num);
			fprintf(fd,"%s\t%i\t%i\t%c\n", chr, intergenic_region_start, intergenic_region_stop, strand);
		}

		void print(_IO_FILE*& fd)
		{
			Region::print(fd);
			fprintf(fd, "gene exons:\t%i\n", (int) exons.size());
			for (int i=0; i<exons.size(); i++)
				fprintf(fd, "\t\t%i\t%i\n", exons[i].first, exons[i].second);
			fprintf(fd, "gene cds_exons:\t%i\n", (int) cds_exons.size());
			for (int i=0; i<cds_exons.size(); i++)
				fprintf(fd, "\t\t%i\t%i\n", cds_exons[i].first, cds_exons[i].second);
			fprintf(fd, "gene utr5exons:\t%i\n", (int) utr5exons.size());
			for (int i=0; i<utr5exons.size(); i++)
				fprintf(fd, "\t\t%i\t%i\n", utr5exons[i].first, utr5exons[i].second);
			fprintf(fd, "gene utr3exons:\t%i\n", (int) utr3exons.size());
			for (int i=0; i<utr3exons.size(); i++)
				fprintf(fd, "\t\t%i\t%i\n", utr3exons[i].first, utr3exons[i].second);

		};

		void print_gff3(_IO_FILE*& fd, int gene_no)
		{
			const char* source = "TranscriptSkimmer";
			char name[1000];
			char* chr = gio->get_contig_name(chr_num);
			sprintf(name, "TrSk%i", gene_no);
			const char* score = "."; 
			const char* phase = ".";

			char attr_str[1000];
			const char* type = "gene";
			sprintf(attr_str, "ID=Gene:%s",name);
			fprintf(fd,"%s\t%s\t%s\t%i\t%i\t%s\t%c\t%s\t%s\n", chr, source, type, start, stop, score, strand, phase,attr_str);

			type = "mRNA";
			sprintf(attr_str, "ID=Transcript:%s.1;Parent=Gene:%s",name,name);
			fprintf(fd,"%s\t%s\t%s\t%i\t%i\t%s\t%c\t%s\t%s\n", chr, source, type, start, stop, score, strand, phase,attr_str);
			
			type = "exon";
			for (int i=0; i<exons.size(); i++)
			{
				sprintf(attr_str, "ID=%s:%s.1.%i;Parent=Transcript:%s.1", type, name, i+1, name);
				fprintf(fd,"%s\t%s\t%s\t%i\t%i\t%s\t%c\t%s\t%s\n", chr, source, type, exons[i].first, exons[i].second, score, strand, phase,attr_str);
			}
			type = "CDS";
			for (int i=0; i<cds_exons.size(); i++)
			{
				sprintf(attr_str, "ID=%s:%s.1.%i;Parent=Transcript:%s.1", type, name, i+1, name);
				fprintf(fd,"%s\t%s\t%s\t%i\t%i\t%s\t%c\t%s\t%s\n", chr, source, type, cds_exons[i].first, cds_exons[i].second, score, strand, phase,attr_str);
			}
			type = "five_prime_UTR";
			for (int i=0; i<utr5exons.size(); i++)
			{
				sprintf(attr_str, "ID=%s:%s.1.%i;Parent=Transcript:%s.1", type, name, i+1, name);
				fprintf(fd,"%s\t%s\t%s\t%i\t%i\t%s\t%c\t%s\t%s\n", chr, source, type, utr5exons[i].first, utr5exons[i].second, score, strand, phase,attr_str);
			}
			type = "three_prime_UTR";
			for (int i=0; i<utr3exons.size(); i++)
			{
				sprintf(attr_str, "ID=%s:%s.1.%i;Parent=Transcript:%s.1", type, name, i+1, name);
				fprintf(fd,"%s\t%s\t%s\t%i\t%i\t%s\t%c\t%s\t%s\n", chr, source, type, utr3exons[i].first, utr3exons[i].second, score, strand, phase,attr_str);
			}

		}

        bool write_window(_IO_FILE*& fd, string tmp_seq, int center_pos, int left_offset, int right_offset, int label);

		bool generate_tis_labels(_IO_FILE*& fd);

		bool generate_tss_labels(_IO_FILE*& fd);

        int get_length()
        {
            return stop-start+1;
        }


};

#endif
