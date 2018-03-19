#/bin/bash

projectPath=$(realpath $PWD/..)
REF=/lustre/rdi/user/hantc/tools/reference/reference_panel8BC_mkdup.cnn

mkdir -p ${projectPath}/Results/cnvkit/lymp
mkdir -p ${projectPath}/Results/cnvkit/tumor

sed '1d' ../Summary/paired.info | while read BLOOD TUMOR LYMP; do
if [ ! -e "${projectPath}/Results/cnvkit/tumor/${TUMOR}.dedup.sorted.recal_reads.cns" ]
  then
    echo "cnvkit.py batch ../Data/mapped/${TUMOR}.dedup.sorted.recal_reads.bam -r ${REF} --output-dir ${projectPath}/Results/cnvkit/tumor --scatter --diagram"
fi
if [ ! -e "${projectPath}/Results/cnvkit/lymp/${LYMP}.dedup.sorted.recal_reads.cns" ]
  then
    echo "cnvkit.py batch ../Data/mapped/${LYMP}.dedup.sorted.recal_reads.bam -r ${REF} --output-dir ${projectPath}/Results/cnvkit/lymp --scatter --diagram"
fi
done > clone/cnvkit.template
perl parallel.pl clone/cnvkit.template

Rscript clone/clone_data.R
sed '1d' ../Summary/paired.info | while read BLOOD TUMOR LYMP; do
cat clone/mutect.head ${projectPath}/Results/mutect/tumor_clone/${TUMOR}.snp.oncotator.vcf | sed s/"CHROM"/"#CHROM"/ > ${projectPath}/Results/mutect/tumor_clone/${TUMOR}.vcf
cat clone/mutect.head ${projectPath}/Results/mutect/lymp_clone/${LYMP}.snp.oncotator.vcf | sed s/"CHROM"/"#CHROM"/ > ${projectPath}/Results/mutect/lymp_clone/${LYMP}.vcf
done

mkdir -p ${projectPath}/Results/cnvkit/lymp_call
mkdir -p ${projectPath}/Results/cnvkit/tumor_call

sed '1d' ../Summary/paired.info | while read BLOOD TUMOR LYMP; do
if [ ! -e "${projectPath}/Results/cnvkit/tumor_call/${TUMOR}.mkdup.call.cns" ]
  then
    echo "cnvkit.py call -o ${projectPath}/Results/cnvkit/tumor_call/${TUMOR}.mkdup.call.cns -v ${projectPath}/Results/mutect/tumor_clone/${TUMOR}.vcf -i ${TUMOR:0:5} -n NORMAL --min-variant-depth 50 ${projectPath}/Results/cnvkit/tumor/${TUMOR}.dedup.sorted.recal_reads.cns"
fi
if [ ! -e "${projectPath}/Results/cnvkit/lymp_call/${LYMP}.mkdup.call.cns" ]
  then
    echo "cnvkit.py call -o ${projectPath}/Results/cnvkit/lymp_call/${LYMP}.mkdup.call.cns -v ${projectPath}/Results/mutect/lymp_clone/${LYMP}.vcf -i ${LYMP:0:5} -n NORMAL --min-variant-depth 50 ${projectPath}/Results/cnvkit/lymp/${LYMP}.dedup.sorted.recal_reads.cns"
fi
done > clone/cnvkit.template2
perl parallel.pl clone/cnvkit.template2

