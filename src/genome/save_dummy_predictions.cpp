#include <stdio.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <iostream>
#include <fstream>
using namespace std;
//#include "common.h"
#include "genome.h"
#define MAXLINE 1000

extern int* find_all_occurences(char* seq, int seq_len, char** consensus, int cons_num, int cons_len, int* num_pos, int offset);
extern void save_score_pos(char* path_prefix, char** score_names, float** scores, int num_scores, int* pos, int num_pos);
extern void complement(char* str, const int len);

/* main
 * **********/
int main(int argc, char** args)
{
	Genome* g = new Genome();

	fprintf(stdout, "expected genome config file as arg1\n");
	fprintf(stdout, "argc: %i\n", argc);
	
	char* fname = NULL;
	char* out_dir = NULL;
	if (argc>1)
	{
		fname = args[1];
		fprintf(stdout, "found arg1: %s\n", fname);
	}
	else
		fprintf(stdout, "expected genome config file as arg1\n");
	if (argc>2)
	{
		out_dir = args[2];
		fprintf(stdout, "found arg2: %s\n", out_dir);
	}
	g->init_genome(fname);

	for (int sig_type=0; sig_type<2; sig_type++)
	{
		char** consensus;	
		int cons_len;
		int offset;
		int cons_num;
		char* sig_name;
		if (sig_type==0)
		{
			sig_name = (char*) "acc";
			cons_num = 1;
			consensus = new char*[cons_num];
			consensus[0] = (char*) "ag";
			cons_len = 2;
			offset = 2;
		}
		else
		{
			sig_name = (char*) "don";
			cons_num = 2;
			consensus = new char*[cons_num];
			consensus[0] = (char*) "gc";
			consensus[1] = (char*) "gt";
			cons_len = 2;
			offset = -1;
		}

		for (int chr_num=0; chr_num<g->num_contigs; chr_num++)
		{
			for (int strand=0; strand<2; strand++)
			{
				char* DNA_seq = g->read_flat_file(chr_num);
				int seq_len = g->contig_len(chr_num);
				if (strand==1)
					g->complement(DNA_seq, seq_len);
	
	
				int* pos = NULL; 
				int num_pos = 0;
				pos = find_all_occurences(DNA_seq, seq_len, consensus, cons_num, cons_len, &num_pos, offset);
	
				char path_prefix[MAXLINE];
				sprintf(path_prefix, "%s/%s/", out_dir, sig_name);
				struct stat st;
				if(stat(path_prefix,&st) != 0)
					if (mkdir(path_prefix, 01777))
					{
						fprintf(stderr, "cannot create dir: %s\n", path_prefix );
						exit(-1);
					}
				if (strand==0)
					sprintf(path_prefix, "%scontig_%d+", path_prefix, chr_num+1);
				else
					sprintf(path_prefix, "%scontig_%d-", path_prefix, chr_num+1);
				char** score_names = new char*[1];
				score_names[0] = (char*) "Conf_cum";
				float** scores = new float*[1]; 
				scores[0] = new float[num_pos];
				for (int i=0; i<num_pos; i++)
					scores[0][i] = 0.0;
				int num_scores = 1;
	
				save_score_pos(path_prefix, score_names, scores, num_scores, pos, num_pos);
			
				delete[] pos;
				delete[] DNA_seq;
				delete[] scores[0];
				delete[] scores;
				delete[] score_names;
			}
		}
		delete[] consensus;
	}
	delete g;	

	return 0;
}

int* find_all_occurences(char* seq, int seq_len, char**  consensus, int cons_num, int cons_len, int* num_pos, int offset)
{
	fprintf(stdout, "find_all_occurences: seq_len:%i cons_len:%i\n", seq_len, cons_len);	

	for (int i=0; i<seq_len-cons_len+1; i++)
	{
		bool match = false;
		for (int c=0; c<cons_num; c++)
		{
			bool cmatch = true;
			for (int j=0; j<cons_len; j++)
				cmatch = cmatch && seq[i+j]==consensus[c][j];
			match = match || cmatch;
		}

		if (match)
			(*num_pos)++;
	}

	fprintf(stdout, "find_all_occurences: seq_len:%i num_pos:%i\n", seq_len, *num_pos);	

	int* pos = new int[*num_pos]; 

	*num_pos = 0;
	for (int i=0; i<seq_len-cons_len+1; i++)
	{
		bool match = false;
		for (int c=0; c<cons_num; c++)
		{
			bool cmatch = true;
			for (int j=0; j<cons_len; j++)
				cmatch = cmatch && seq[i+j]==consensus[c][j];
			match = match || cmatch;
		}

		if (match)
		{
			pos[*num_pos] = i+offset+1;
			(*num_pos)++;
		}
	}
	return pos; 
}


