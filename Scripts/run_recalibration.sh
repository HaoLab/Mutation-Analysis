cat $1 | while read SAMPLE FILE; do
echo "bash recalibration.sh ${SAMPLE}"
done > recalibration.template
perl parallel.pl recalibration.template
