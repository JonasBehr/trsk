#ifndef _SPLICE_LABELS_FROM_RNA_SEQ_H__
#define _SPLICE_LABELS_FROM_RNA_SEQ_H__

#include <stdio.h>
#include <assert.h>
#include <sys/stat.h>
#include <stdint.h>
#include <algorithm>
#include "get_reads_direct.h"
#include "genome.h"

#include <vector>
  using std::vector;
#include <string>
  using std::string;
#include "splice_labels_from_RNA_seq.h"


#define MAXLINE 1000
typedef struct {
	char** consensus;	 
	int cons_len; 
	int offset; 
	int cons_num; 
	char* sig_name; 
} signal_t;
typedef struct {
	int pos;
	int label;
} example_t;
extern void save_score_pos(char* path_prefix, char** score_names, float** scores, int num_scores, int* pos, int num_pos);

void vec_unique(vector<int>* vec);
void find_all_occurences(vector<int>* pos, char* seq, int seq_len, signal_t sig, int start_pos);
void efficient_clear(vector<int>* vec, uint32_t pos);
void remove_intersection(vector<int>* vec1, vector<int>* vec2);
bool exp_compare(example_t exp1, example_t exp2);
void combine_with_label(vector<example_t>* all_exp, vector<int> pos, int label);
void create_dir(char* path_prefix, char* dn_output, int chr_num, signal_t sig, char strand);
bool check_consensus(int pos, signal_t sig, char* seq);
void cons_hist(vector<int> pos, signal_t sig, char* DNA_seq);
int get_splice_labels(int argc, char* argv[]);

#endif
