GATK=/lustre/rdi/user/songx/tools/software/GATK/GenomeAnalysisTK.jar
ref=/lustre/rdi/user/songx/tools/genome/hg19/ucsc.hg19.fasta
database_snp=/lustre/rdi/user/songx/tools/Database/dbsnp_138.hg19.vcf
database_indel=/lustre/rdi/user/songx/tools/Database/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
region=/lustre/rdi/user/songx/tools/genome/regions
Cosmic=/lustre/rdi/user/songx/tools/Database/Cosmic.hg19.chr.vcf
group="lymp_chr"
SAMPLE=$1


java -Xmx8g -Djava.io.tmpdir=../Data/mapped -jar $GATK -T BaseRecalibrator -R $ref -I ../Data/mapped/${SAMPLE}.mkdup.bam -o ../Data/mapped/${SAMPLE}.dedup.sorted.grp -knownSites $database_snp -knownSites $database_indel
java -Xmx8g -Djava.io.tmpdir=../Data/mapped -jar $GATK -T PrintReads -R $ref -I ../Data/mapped/${SAMPLE}.mkdup.bam -BQSR ../Data/mapped/${SAMPLE}.dedup.sorted.grp -o ../Data/mapped/${SAMPLE}.dedup.sorted.recal_reads.bam



