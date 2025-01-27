#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

samplename="$1"; shift # eg PD37766b, will be how html file is saved
image_file="$1"; shift # png histology image
location_tbl="$1" # table giving pixel position of each microbiopsy
clonal_prev="$2" # clonal prevelance table generated during ndp clustering
edge_tbl="$3" # edge table generated durin ndp tree generation
low_prev="$4" # whether to include low prevelence clones <0.01
outdir="$5" 
mut_tbl="$6" # optional

echo $samplename

# Run R code
  /usr/local/bin/Rscript $SCRIPT_DIR/mapscape_run.R \
  $samplename \
  $image_file \
  $location_tbl \
  $clonal_prev \
  $edge_tbl \
  $low_prev \
  $outdir \
  $mut_tbl 
