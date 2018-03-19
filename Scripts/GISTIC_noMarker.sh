#!/bin/sh
## run example GISTIC analysis

## output directory
echo --- creating output directory ---
basedir=$3
mkdir -p $basedir 

echo --- creating input gistic segment file ---
Rscript /lustre/rdi/user/hantc/tools/GISTIC/cnr4gistic.R -i $1 -o $2

echo --- running GISTIC ---
## input file definitions
segfile=$2
#markersfile=/lustre/rdi/user/hantc/Pilot-Work/171121-Zhuang/Summary/marker_noBG.txt
refgenefile=/lustre/rdi/user/hantc/app/GISTIC/refgenefiles/hg19.mat
#alf=`pwd`/examplefiles/arraylistfile.txt
#cnvfile=`pwd`/examplefiles/cnvfile.txt
## call script that sets MCR environment and calls GISTIC executable 
/lustre/rdi/user/hantc/app/GISTIC/gistic2 -b $basedir -seg $segfile -refgene $refgenefile -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.95 -armpeel 1 -savegene 1 -gcm extreme -ta 0.7 -td 0.7 -qvt 0.05
# -alf $alf -cnv $cnvfile


