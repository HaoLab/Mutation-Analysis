#!/usr/bin/bash

MutSigCV=/work-z/user/guoh/soft/MutSigCV_1.41/run_MutSigCV.sh
MCR=/work-z/user/guoh/soft/MCR_R2016a/R2016a/v901
cov_file=/work-z/user/guoh/Data/MutSigRef/exome_full192.coverage.txt
gene_covariates=/work-z/user/guoh/Data/MutSigRef/gene.covariates.txt
mutation_type_dictionary_file=/remote/rdi/user/songx/Project/20171018_nanfang_gas/MutSigCV/mutation_type_dictionary_file.txt
chr_files_hg19=/work-z/user/guoh/Data/MutSigRef/chr_files_hg19/
AWD_NED=$1/tumor_all.maf
CD=$1/lymp_all.maf
output_AWD=$1/tumor
output_CD=$1/lymp

$MutSigCV $MCR $AWD_NED $cov_file $gene_covariates $output_AWD/tumor_all $mutation_type_dictionary_file $chr_files_hg19
$MutSigCV $MCR $CD $cov_file $gene_covariates $output_CD/lymp_all $mutation_type_dictionary_file $chr_files_hg19
