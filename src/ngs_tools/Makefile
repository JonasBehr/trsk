CXX = g++
CXXFLAGS = -O9 -g -Wall -fmessage-length=0
MATLAB_DIR = /fml/ag-raetsch/share/software
MEX_BIN ?= $(MATLAB_DIR)/matlab-7.6/bin/mex
MATLAB_INCDIR ?= $(MATLAB_DIR)/matlab/extern/include
MATLAB_INCL = -I$(MATLAB_INCDIR)
SAMDIR ?= /fml/ag-raetsch/share/software/samtools
GENOME ?= ../../projects/mGene_core/genome

all: get_reads.mexa64 get_reads.mex spliced_nucleotide_db subsample

#splice_labels_from_RNA_seq:	run_splice_labels_from_RNA_seq.cpp splice_labels_from_RNA_seq.cpp get_reads_direct.cpp read.cpp
#	$(CXX) -Wall run_splice_labels_from_RNA_seq.cpp splice_labels_from_RNA_seq.cpp get_reads_direct.cpp read.cpp  ${GENOME}/genome.cpp ${GENOME}/tools.cpp -I${GENOME} -I${SAMDIR} -L${SAMDIR} -lbam -lz -lcurses -o splice_labels_from_RNA_seq

spliced_nucleotide_db:	spliced_nucleotide_db.cpp get_reads_direct.cpp read.cpp
	$(CXX) -O9 -g spliced_nucleotide_db.cpp get_reads_direct.cpp read.cpp  ${GENOME}/genome.cpp ${GENOME}/tools.cpp -I${GENOME} -I${SAMDIR} -L${SAMDIR} -lbam -lz -lcurses -o spliced_nucleotide_db

get_reads.mexa64:	get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp read.h
	${MEX_BIN} -g get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp -I${SAMDIR} -L${SAMDIR} -lbam -lz -lcurses $(MATLAB_INCL)

get_reads.mex:	get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp
	mkoctfile --mex get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp -I${SAMDIR} -L${SAMDIR} -lbam -lz -lcurses

filter_wrong_quality_string_length:	filter_wrong_quality_string_length.cpp
	$(CXX) -O9 -g filter_wrong_quality_string_length.cpp  -o filter_wrong_quality_string_length

subsample:
	$(CXX) -Wall subsample.cpp -o subsample

clean:	
	rm -f *.o *.mexa64 *.mex get_reads_direct
