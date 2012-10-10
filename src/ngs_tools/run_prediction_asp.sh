#!/bin/bash
export PYTHONPATH=/fml/ag-raetsch/home/jonas/python/usr/local/lib/python2.6/dist-packages/
export PYTHONPATH=$PYTHONPATH:/fml/ag-raetsch/home/cwidmer/svn/projects/multitask/python/base/
export PYTHONPATH=$PYTHONPATH:/fml/ag-raetsch/home/cwidmer/svn/projects/arts2/
export PYTHONPATH=$PYTHONPATH:/fml/ag-raetsch/home/jonas/svn/tools/python/
export PYTHONPATH=$PYTHONPATH:/fml/ag-raetsch/home/jonas/shogun/trunk/applications/asp/
export PYTHONPATH=$PYTHONPATH:/fml/ag-raetsch/home/jonas/svn/tools/ngs/
export LD_LIBRARY_PATH=/fml/ag-raetsch/home/jonas/shogun/trunk/src/libshogun:/fml/ag-raetsch/home/jonas/shogun/trunk/src/libshogunui


fafname=/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/genomes/japonica2/japonica2_condensed.fasta

~/shogun/trunk/applications/asp/asp --svm_type primal $fafname 
