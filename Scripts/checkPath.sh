cat $1 | while read SAMPLE; do
ls -d /work/commercial_test/pipeline_V1/*/*${SAMPLE}*/qc_sample/*qc_sample.txt >> $2
done

sed -i s#"/qc_sample/".*"txt"#""#g $2