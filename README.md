# Mutation-Analysis

## 说明
* 目前包括的模块
1. 基础模块：样本路径查找、QC
2. Calling：Mutect2、cnvkit
3. 分析模块：Mutect2结果过滤、oncotator、MutSigCV、GISTIC、HotNet2
* 需要提供的列表在Summary文件夹中均有示例
* 可以在 /lustre/rdi/user/hantc/Pilot-Work/180116-Zhuang 路径下直接查看整个流程

## 1.1 checkPath.sh ../Summary/New.list ../Summary/samplepath.list
* 该步骤是在已知sample id的情况下，在新pipeline文件夹中找到样本的路径，基本使用新id命名规则（特别长）的样本都可以通过该方法找到
* input: SAMPLE_ID
* output: SAMPLE_PATH
* run: bash checkPath.sh $input $output

## 1.2 Rscript QC_pipeline.R -f ../Summary/pipeline.list -o ../Summary
## 1.3 Rscript /lustre/rdi/user/hantc/tools/QC/QC.R ../Summary/PID.txt old
* 分别对pipeline文件夹中的样本和旧系统中的样本进行QC
* QC_pipeline.R
* f: NAME and 文件夹路径，文件夹路径需要到样本那级（即再下级为qc_sample/basic_analysis等）
* o: 输出文件夹路径
* run: Rscript QC_pipeline.R -f ../Summary/pipeline.list -o ../Summary
* /lustre/rdi/user/hantc/tools/QC/QC.R
* input: PID列表
* output: 文件前缀（输出为${output}_QC.txt）
* run: Rscript /lustre/rdi/user/hantc/tools/QC/QC.R $input $output

## 2.1 link.sh ../Summary/link_pipeline.list dedup.sorted.recal_reads
## 2.2 link.sh ../Summary/link_old.list mkdup
* 该步骤将已有的bam文件链接到新建的../Data/mapped文件夹中，后缀由传入的第二个参数指定
* 本次新加入的样本中，一部分跑过商检pipeline的流程，直接链接realign以后的bam文件；另一部分无realign bam的，链接的是sorted.filtered.mkdup.bam
* input: NAME and FILEPATH
* run: bash link.sh $input $suffix

## 3.1 run_recalibration.sh ../Summary/link_old.list
## 3.2 run_recalibration.sh ../Summary/firsttime.list
* 对无realign bam的样本先进行recalibration，此处用到了parallel.pl，需自行调整parallel.pl中设置的并行数量
* 上次分析的旧样本记录在../Summary/firsttime.list中，主要是用于后续的cnvkit分析
* 上次样本中P6372_FFG_1a有pipeline生成的realign.bam，手动链接
* input: NAME and FILEPATH
* run: bash run_recalibration.sh $input

## 4.1 run_mutect2.sh ../mutect.list
* 除input以外，调用的脚本与文件包括mutect/mutect2.sh, mutect/mutect2_script.template, /lustre/rdi/user/hantc/Pilot-Work/180116-Zhuang/Data/bed/panel8_{1..20}.bed
* 脚本中默认将panel8按大小分为20个bed，对于其他panel暂时没有现成的分割文件，需要自行分割并修改调用的mutect.sh脚本中的panel路径
* 同样用到了parallel.pl，需自行调整parallel.pl中设置的并行数量，HPC上建议为15个，不要过多
* 该步骤会生成../Results/mutect/${GROUP}_chr_raw文件夹，该文件夹中有按照${TUMOR_NAME}命名的多个文件夹，包含了chr1-chr20共计20个vcf文件（按照bed1-bed20分割）
* input: GROUP, NORMAL_NAME and TUMOR_NAME
* run: bash run_mutect2.sh $input

## 5.1 raw2bedonly2somatic.sh /lustre/rdi/user/hantc/tools/bed/panel8_ctCNV.bed
* 除input以外，调用的脚本与文件包括mutect/mutect_merge.R与mutect/Overlap_mutect.R
* 前者合并多个chr*vcf，生成文件在../Results/mutect/${GROUP}_bedonly中
* 后者进行过滤与一致性统计，生成结果在../Analysis/Overlap_snp/indel_mutect.txt，另有一些过滤后的vcf文件在../Results/mutect的各个文件夹中
* all过滤标准为：FILTER=PASS；NORMAL_VAF<0.4；ExAC_ALL<0.005；thousands<0.005；黑名单人群频率<0.1,；正常人人群频率<0.2；
* f1过滤标准为：all+TUMOR_VAF >= 0.03
* f2过滤标准为：all+TUMOR_VAF >= 0.05
* input: 完整的bed文件
* notes: 如果有旧样本，在运行脚本中的annovar注释前将合并后的mutect vcf文件拷贝至同一文件夹
* run: bash raw2bedonly2somatic.sh $input

## 6.1 oncotator2mutsig.sh
* 调用的脚本包括oncotator文件夹中的所有脚本
* 包括了oncotator注释与mutsigCV两步，均需要在pilot-a上运行
* run: bash oncotator2mutsig.sh

## 7.1 cnvkit.sh
* 调用的脚本包括clone文件夹中的全部与parallel.pl，调用文件包括panel8的正常人baseline reference与配对列表../Summary/paired.info
* 配对列表需要header: blood, tumor and lymp
* 注意本次的分组为配对分组，如果分组情况/组名不同，需要调整代码
* run: bash cnvkit.sh

## 8.1 run_GISTIC.sh ../Results/cnvkit/lymp ../Results/GISTIC/lymp
## 8.2 run_GISTIC.sh ../Results/cnvkit/tumor ../Results/GISTIC/tumor
* 调用/lustre/rdi/user/hantc/tools/GISTIC下的脚本，GISTIC安装在/lustre/rdi/user/hantc/app/GISTIC下，参数设置中amp cutoff为0.7（CN=4），del cutoff为0.7（CN=1）
* input: cnvkit的结果文件夹
* output: GISTIC的输出文件夹
* run: bash run_GISTIC.sh $input $output

## 9.1 hotnet2.sh tumor
## 9.2 hotnet2.sh lymp
* 软件安装在/lustre/rdi/user/hantc/app/hotnet2下，需要使用matlab，不确定是否需要在自己目录下安装
* 脚本中调用的其他脚本在hotnet目录下，利用maf和gistic生成的文件转化为hotnet的输入文件
* 使用的maf为过滤0.1% VAF后的MAF，如果要使用其他，可修改脚本中prefix变量
* input: group(../Summary/mutect.list中的group)
* run: bash hotnet2.sh $input
