# get vcf for oncotator
Rscript oncotator/vcf4oncotator.R

# oncotator (must be run on pilot-a)
mkdir -p ../Results/mutect/lymp_oncotator
mkdir -p ../Results/mutect/tumor_oncotator
python oncotator/run_oncotator.py -i ../Results/mutect/lymp_foroncotator/ -o ../Results/mutect/lymp_oncotator/
python oncotator/run_oncotator.py -i ../Results/mutect/tumor_foroncotator/ -o ../Results/mutect/tumor_oncotator/

mkdir -p ../Results/mutect/lymp_oncotator_all
mkdir -p ../Results/mutect/tumor_oncotator_all
python oncotator/run_oncotator.py -i ../Results/mutect/lymp_foroncotator_all/ -o ../Results/mutect/lymp_oncotator_all/
python oncotator/run_oncotator.py -i ../Results/mutect/tumor_foroncotator_all/ -o ../Results/mutect/tumor_oncotator_all/

mkdir -p ../Results/mutect/lymp_oncotator_001
mkdir -p ../Results/mutect/tumor_oncotator_001
python oncotator/run_oncotator.py -i ../Results/mutect/lymp_foroncotator_001/ -o ../Results/mutect/lymp_oncotator_001/
python oncotator/run_oncotator.py -i ../Results/mutect/tumor_foroncotator_001/ -o ../Results/mutect/tumor_oncotator_001/


# merge maf
Rscript oncotator/Merge_maf.R

# mutsig (must be run on pilot-a)
bash oncotator/MutSigCV.sh /remote/rdi/user/hantc/Pilot-Work/180116-Zhuang/Analysis/
bash oncotator/MutSigCV_all.sh /remote/rdi/user/hantc/Pilot-Work/180116-Zhuang/Analysis/
bash oncotator/MutSigCV_001.sh /remote/rdi/user/hantc/Pilot-Work/180116-Zhuang/Analysis/


