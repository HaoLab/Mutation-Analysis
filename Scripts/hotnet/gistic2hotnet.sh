input=$1
bed=$2
output=$3

bedtools intersect -a <(awk '{print "chr"$2"\t"$3"\t"$4"\t"$1" "$5" "$6" "$7" "$8}' $input | sed '1d') -b $bed -wb | awk -F"[ ]" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}'|cut -f1-3,12,4-8 > $output
