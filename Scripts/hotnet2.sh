#!/usr/bin/env bash

# These commands run HotNet2 on the processed network and mutation data used for our Nature Genetics paper.

# Set location of HotNet2, number of cores, number of network permutations, and number of heat permutations.
hotnet2=/lustre/rdi/user/hantc/app/hotnet2
num_cores=5
num_network_permutations=100
num_heat_permutations=1000
samplesize=25
bed=/lustre/rdi/user/hantc/tools/bed/panel8_ctCNV.bed

# GISTIC output
gistic_out=../Results/GISTIC/$1/

# MutSigCV output
prefix=${1}_001
mafpath=../Analysis/
mutsigpath=../Analysis/$1/
outputpath=../Results/hotnet2/$1


mkdir -p ${outputpath}
bash hotnet/gistic2hotnet.sh $gistic_out/scores.gistic $bed $outputpath/gistic_${prefix}.scores
Rscript hotnet/gistic2hotnet.R -p ${prefix} -o ${outputpath} -i $outputpath/gistic_${prefix}.scores -n ${samplesize}

Rscript hotnet/maf2hotnet.R -n ${samplesize} -p ${prefix} -o ${outputpath} -m ${mafpath}/${prefix}.maf -s ${mutsigpath}/${prefix}.sig_genes.txt


# Create heat data.
python $hotnet2/makeHeatFile.py \
    scores \
    -hf ${outputpath}/${prefix}_freq.txt \
    -o  ${outputpath}/${prefix}_freq.json \
    -n  ${prefix}.freq

python $hotnet2/makeHeatFile.py \
    scores \
    -hf ${outputpath}/${prefix}_mutsig.txt \
    -o  ${outputpath}/${prefix}_mutsig.json \
    -n  ${prefix}.mutsig

python $hotnet2/makeHeatFile.py \
    scores \
    -hf ${outputpath}/${prefix}_freq_amp.txt \
    -o  ${outputpath}/${prefix}_freq_amp.json \
    -n  ${prefix}.freq_amp


python $hotnet2/makeHeatFile.py \
    scores \
    -hf ${outputpath}/${prefix}_gistic_amp.txt \
    -o  ${outputpath}/${prefix}_gistic_amp.json \
    -n  ${prefix}.gistic_amp
	
python $hotnet2/makeHeatFile.py \
    scores \
    -hf ${outputpath}/${prefix}_freq_del.txt \
    -o  ${outputpath}/${prefix}_freq_del.json \
    -n  ${prefix}.freq_del


python $hotnet2/makeHeatFile.py \
    scores \
    -hf ${outputpath}/${prefix}_gistic_del.txt \
    -o  ${outputpath}/${prefix}_gistic_del.json \
    -n  ${prefix}.gistic_del

# Run HotNet2.
python $hotnet2/HotNet2.py \
    -nf  $hotnet2/paper/data/networks/hint+hi2012/hint+hi2012_ppr_0.4.h5 \
         $hotnet2/paper/data/networks/irefindex9/irefindex9_ppr_0.45.h5 \
         $hotnet2/paper/data/networks/multinet/multinet_ppr_0.5.h5 \
    -pnp $hotnet2/paper/data/networks/hint+hi2012/permuted/hint+hi2012_ppr_0.4_##NUM##.h5 \
         $hotnet2/paper/data/networks/irefindex9/permuted/irefindex9_ppr_0.45_##NUM##.h5 \
         $hotnet2/paper/data/networks/multinet/permuted/multinet_ppr_0.5_##NUM##.h5 \
    -hf  ${outputpath}/${prefix}_freq.json \
         ${outputpath}/${prefix}_mutsig.json \
         ${outputpath}/${prefix}_freq_amp.json \
         ${outputpath}/${prefix}_gistic_amp.json \
         ${outputpath}/${prefix}_freq_del.json \
         ${outputpath}/${prefix}_gistic_del.json \
    -np  $num_network_permutations \
    -hp  $num_heat_permutations \
    -o   ${outputpath}/ \
    -c   $num_cores

