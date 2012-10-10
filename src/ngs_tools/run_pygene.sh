#!/bin/bash
#export PYTHONPATH=/fml/ag-raetsch/home/jonas/python/usr/local/lib/python2.5/site-packages/
export PYTHONPATH=/fml/ag-raetsch/home/jonas/python/usr/local/lib/python2.6/dist-packages/
#export PYTHONPATH=/fml/ag-raetsch/home/cwidmer/lib/python2.5/site-packages/
#export PYTHONPATH=$PYTHONPATH:/fml/ag-raetsch/home/jonas/shogun/trunk/src:/fml/ag-raetsch/home/jonas/shogun/trunk/src/python_modular
export PYTHONPATH=$PYTHONPATH:/fml/ag-raetsch/home/cwidmer/svn/projects/multitask/python/base/
export PYTHONPATH=$PYTHONPATH:/fml/ag-raetsch/home/cwidmer/svn/projects/arts2/
export PYTHONPATH=$PYTHONPATH:/fml/ag-raetsch/home/jonas/svn/tools/python/
#export LD_LIBRARY_PATH=/fml/ag-raetsch/home/jonas/shogun/python/lib/
#export LD_LIBRARY_PATH=/fml/ag-raetsch/home/cwidmer/lib/ ## chris shogun fuer pygene
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/fml/ag-raetsch/home/cwidmer/lib/ ## chris shogun fuer pygene
export LD_LIBRARY_PATH=/fml/ag-raetsch/home/jonas/shogun/trunk/src/libshogun:/fml/ag-raetsch/home/jonas/shogun/trunk/src/libshogunui

python pygene.py

