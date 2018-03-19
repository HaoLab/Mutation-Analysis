mkdir -p ../Results/mutect

GATK=/lustre/rdi/user/songx/tools/software/GATK/GenomeAnalysisTK.jar
ref=/lustre/rdi/user/songx/tools/genome/hg19/ucsc.hg19.fasta
database_snp=/lustre/rdi/user/songx/tools/Database/dbsnp_138.hg19.vcf
database_indel=/lustre/rdi/user/songx/tools/Database/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
region=/lustre/rdi/user/hantc/Pilot-Work/180116-Zhuang/Data/bed
Cosmic=/lustre/rdi/user/songx/tools/Database/Cosmic.hg19.chr.vcf

group=$1
Normal=$2
Tumor=$3
i=$4

mkdir -p ../Results/mutect/${group}_chr_tmp
mkdir -p ../Results/mutect/${group}_chr_raw

mkdir -p ../Results/mutect/${group}_chr_raw/${Tumor}

bed=${region}/panel8_${i}.bed
java -Xmx8g -Djava.io.tmpdir=../Results/mutect/${group}_tmp -jar $GATK -T MuTect2 -R $ref -I:normal /lustre/rdi/user/hantc/Pilot-Work/180116-Zhuang/Data/mapped/${Normal}.dedup.sorted.recal_reads.bam -I:tumor /lustre/rdi/user/hantc/Pilot-Work/180116-Zhuang/Data/mapped/${Tumor}.dedup.sorted.recal_reads.bam -L ${bed} -o /lustre/rdi/user/hantc/Pilot-Work/180116-Zhuang/Results/mutect/${group}_chr_raw/${Tumor}/${Tumor}_chr${i}.dedup.vcf --dbsnp $database_snp --cosmic $Cosmic
