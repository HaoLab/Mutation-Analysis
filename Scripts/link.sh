mkdir -p ../Data/mapped

cat $1 | while read NAME FILE; do
ln -s ${FILE} /lustre/rdi/user/hantc/Pilot-Work/180116-Zhuang/Data/mapped/${NAME}.$2.bam
ln -s ${FILE}.bai /lustre/rdi/user/hantc/Pilot-Work/180116-Zhuang/Data/mapped/${NAME}.$2.bam.bai
done

