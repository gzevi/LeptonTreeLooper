#!/bin/bash

# ---------------------------------------------------------------------------------- #
# Simple wrapper script to create the LeptonTree.h/cc wrappers for the looper.
# It uses the makeCMS2ClassFiles.C
#
# Usage: makeLeptonTreeClassFiles.sh lepton_tree_file_name
# ---------------------------------------------------------------------------------- #

if [ $# -ne 1 ]; then
    echo "Usage: makeLeptonTreeClassFiles.sh lepton_tree_file_name"
    exit 1
fi

# file name is the only argument
file_name=$1

# create the LeptonTree.cc/h
#pushd $TNP/tools
root_cmd="makeCMS3ClassFiles.C+(\"$file_name\", \"t\", \"LeptonTree\", \"lepton_tree\", \"lepton_tree_obj\")"
echo "[makeLeptonTreeClassFiles.sh] calling"
echo $root_cmd
root -b -l -q "$root_cmd" 

# move LeptonTree.h
#mv LeptonTree.h ..//interface/.

# move LeptonTree.cc and change include path for LeptonTree.h to conform with SCRAM convention
#sed 's/\#include "LeptonTree.h"/\#include "TagAndProbe\/Analysis\/interface\/LeptonTree.h"/g' LeptonTree.cc > ../src/LeptonTree.cc
#rm LeptonTree.cc

# remove example files
#rm doAll.C
#rm ScanChain.C
 
# done
#popd
echo "[makeLeptonTreeClassFiles.sh] finished"
