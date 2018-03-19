bash mutect/mutect2_script.template $1 > mutect/mutect2_script.txt
perl parallel.pl mutect/mutect2_script.txt
