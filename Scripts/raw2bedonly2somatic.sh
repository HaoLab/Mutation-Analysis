projectPath=$(realpath $PWD/..)
BED=$1

mkdir -p ${projectPath}/Results/mutect/tumor_chr_bedonly
grep ^tumor ../Summary/mutect.list | cut -f2,3 |  while read NORMAL TUMOR; do
mkdir -p ${projectPath}/Results/mutect/tumor_chr_bedonly/${TUMOR}
bash /lustre/rdi/user/hantc/tools/vcf2bedvcf/mutect2bedvcf.sh ${projectPath}/Results/mutect/tumor_chr_raw/${TUMOR} ${BED} 135
mv ${projectPath}/Results/mutect/tumor_chr_raw/Results_bedonly/* ${projectPath}/Results/mutect/tumor_chr_bedonly/${TUMOR}
rm -rf ${projectPath}/Results/mutect/tumor_chr_raw/Results_bedonly
sed -i s/"CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1"/"CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR	NORMAL"/g ${projectPath}/Results/mutect/tumor_chr_bedonly/${TUMOR}/*
done

mkdir -p ${projectPath}/Results/mutect/lymp_chr_bedonly
grep ^lymp ../Summary/mutect.list | cut -f2,3 | while read NORMAL TUMOR; do
mkdir -p ${projectPath}/Results/mutect/lymp_chr_bedonly/${TUMOR}
bash /lustre/rdi/user/hantc/tools/vcf2bedvcf/mutect2bedvcf.sh ${projectPath}/Results/mutect/lymp_chr_raw/${TUMOR} ${BED} 135
mv ${projectPath}/Results/mutect/lymp_chr_raw/Results_bedonly/* ${projectPath}/Results/mutect/lymp_chr_bedonly/${TUMOR}
rm -rf ${projectPath}/Results/mutect/lymp_chr_raw/Results_bedonly
sed -i s/"CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1"/"CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	lymp	NORMAL"/g ${projectPath}/Results/mutect/lymp_chr_bedonly/${TUMOR}/*
done

Rscript mutect/mutect_merge.R

sed -i s/"CHROM"/"#CHROM"/g ${projectPath}/Results/mutect/tumor_bedonly/*vcf
sed -i s/"CHROM"/"#CHROM"/g ${projectPath}/Results/mutect/lymp_bedonly/*vcf

# if there are samples that has been analysed, copy *dedup.vcf to ${projectPath}/Results/mutect/tumor_bedonly (not *.snp.dedup.vcf or *.indel.dedup.vcf)
python /lustre/rdi/user/guoh/tools/annovar_VCF_subs2.py ${projectPath}/Results/mutect/tumor_bedonly
python /lustre/rdi/user/guoh/tools/annovar_VCF_subs2.py ${projectPath}/Results/mutect/lymp_bedonly

mkdir -p ../Analysis
Rscript mutect/Overlap_mutect.R
