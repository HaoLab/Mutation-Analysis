#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
'''
Description: run oncotator for somatic vcf directory.
Usage: python run_oncotator.py dir
'''
import os, sys, time, re, glob, itertools, gzip
import subprocess, multiprocessing, shlex
import sh,glob
#from timeFunc import timeit
import argparse

varscan_dir = "/work-z/user/guoh/research/2017-03-23_CCES_2nd/seqdata/varscan/"
varscan_dir2 = "/work-z/user/guoh/research/2017-03-23_CCES_2nd/seqdata/varscan2/"
cfDNA_Bam = "/work-z/user/guoh/research/2017-03-23_CCES_2nd/seqdata/cfDNA_Bam/"
somatic_dir = "/remote/rdi/user/hantc/Pilot-Work/171121-Zhuang/Results/mutect/lymp_foroncotator/"
oncotator_dir = "/remote/rdi/user/hantc/Pilot-Work/171121-Zhuang/Results/mutect/lymp_oncotator/"
if not os.path.exists(oncotator_dir):
    os.mkdir(oncotator_dir)
#print "test"

#@timeit
def main():
    parser = argparse.ArgumentParser(description='run_oncotator.py, run oncotator on pilot-a.', epilog = 'Contact Hao Guo <guo.hao@genecast.com.cn> for help.')
    parser.add_argument('-i','--inputDIR', default=somatic_dir, help='Input vcf file directory.')
    parser.add_argument('-o', '--outputDIR', default=oncotator_dir, help='Output maf file directory.')
    args = parser.parse_args()
    #print args
    vcf_list = glob.glob(args.inputDIR+"/*.vcf")
    #print vcf_list
    call_concotator_mp(vcf_list, args.outputDIR)


def call_concotator_mp(vcf_list, outputDIR):
    '''Usage: call run_oncotator parallely'''
    pool = multiprocessing.Pool(20)
    status = []
    for vcf in vcf_list:
        #print vcf + "\n"
        tmpCMD = run_oncotator(vcf, outputDIR)
        command_line =  ['/bin/bash', '-c', tmpCMD]
        result = pool.apply_async(subprocess.check_call, (command_line,))
        status.append(result)
    pool.close()
    pool.join()
    rc = 0
    for st in status:
        if st.successful():
            rc += 0
        else:
            rc += 1
    return rc


def run_oncotator(inputVCF, outputDIR):
    '''Usage: return oncotator CMD string'''
    cmd1 = "oncotator --input_format=VCF --db-dir=/work-a/user/guoh/Data/oncotator/oncotator_v1_ds_Jan262014 --output_format=TCGAMAF --collapse-number-annotations "
    cmd2 = inputVCF + " "
    cmd3 = outputDIR + "/" + os.path.basename(inputVCF).split('.')[0] + ".maf hg19"
    return cmd1 + cmd2 + cmd3


if __name__ == "__main__":
    main()



