cat ../Summary/pipeline_pid.list | while read PID; do
grep ${PID} /work/commercial_test/pipeline_V1/1711*/*/qc_sample/*qc_sample.txt >> ../Summary/PID_qc_sample.txt
done
