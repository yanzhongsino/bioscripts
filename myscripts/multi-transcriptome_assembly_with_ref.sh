#!/bin/bash
# multi-transcriptome_assembly_with_ref.sh脚本用于多转录本有参组装，此脚本的目的是对测序的多对RNA-seq数据（同一个生物个体的不同器官）进行参考基因组的转录组组装。一般在为基因组测序之后做基因组基因结构注释会用到多组RNA-seq数据。最终获取这个个体的转录本和cds序列。

# 软件：提前安装好hisat2,stringtie,gffread,TransDecoder,samtools并把命令添加到环境变量。

# 数据：组装好的基因组genome，四个器官的转录组数据(rna_1)、(rna_2)、(rna_3)、(rna_4)。

# 分析流程如下:
# 1.使用hisat2将rna_seq数据比对到reference
# 2.使用stringtie进行转录本预测
# 3.使用TransDecoder进行编码区预测
# 4.使用gffread进行格式转换


# 以下需要根据您的数据更改或填入，并把行首的#删除：
# thread=4 #多个器官会使用多个thread，建议4/8线程就够了。
# sample=sampleID #物种名，或者样品ID
# genome=/path/to/genome #基因组数据，建议用绝对路径
# rna_1.1=/path/to/RNA-seq.rna_1_1 #转录组器官1的数据R1
# rna_1.2=/path/to/RNA-seq.rna_1_2 #转录组器官1的数据R2
# rna_2.1=/path/to/RNA-seq.rna_2_1 #转录组器官2的数据R1
# rna_2.2=/path/to/RNA-seq.rna_2_2 #转录组器官2的数据R2
# rna_3.1=/path/to/RNA-seq.rna_3_1 #转录组器官3的数据R1
# rna_3.2=/path/to/RNA-seq.rna_3_2 #转录组器官3的数据R2
# rna_4.1=/path/to/RNA-seq.rna_4_1 #转录组器官4的数据R1
# rna_4.2=/path/to/RNA-seq.rna_4_2 #转录组器官4的数据R2


# 0.制作结果目录
mkdir output && cd output
mkdir 1.hisat2
mkdir 2.stringtie
mkdir 3.gffread
mkdir 3.transdecoder

# 1.hisat2将rna_seq数据比对到reference
cd 1.hisat2
mkdir index && cd index 
hisat2-build ${genome} ${sample}

cd ../../
hisat2 --dta -p ${thread} -x 1.hisat2/index/${sample} -1 ${rna_1.1} -2 ${rna_1.2} |samtools sort -@ 12 >1.hisat2/rna_1.bam &
hisat2 --dta -p ${thread} -x 1.hisat2/index/${sample} -1 ${rna_2.1} -2 ${rna_2.2} |samtools sort -@ 12 >1.hisat2/rna_2.bam &
hisat2 --dta -p ${thread} -x 1.hisat2/index/${sample} -1 ${rna_3.1} -2 ${rna_3.2} |samtools sort -@ 12 >1.hisat2/rna_3.bam &
hisat2 --dta -p ${thread} -x 1.hisat2/index/${sample} -1 ${rna_4.1} -2 ${rna_4.2} |samtools sort -@ 12 >1.hisat2/rna_4.bam &

wait
echo "hisat2 done"

# 2.stringtie组装转录本
samtools merge -@ $((${thread}*4)) 1.hisat2/merged.bam 1.hisat2/rna_1.bam 1.hisat2/rna_2.bam 1.hisat2/rna_3.bam 1.hisat2/rna_4.bam ## 合并所有bam文件
stringtie -p $((${thread}*4)) -o 2.stringtie/merged.gtf 1.hisat2/merged.bam & ##组装merged的转录本
wait
echo "stringtie done"

## 3.TransDecoder预测编码区
gtf_genome_to_cdna_fasta.pl 2.stringtie/merged.gtf ${genome} >3.transdecoder/merged.transcripts.fa #把gtf文件转换成fa文件
gtf_to_alignment_gff3.pl 2.stringtie/merged.gtf >3.transdecoder/merged.transcripts.gff3 #把gtf文件转换成gff3文件
cd 3.transdecoder
TransDecoder.LongOrfs -t merged.transcripts.fa #提取转录本中长的开放阅读框，可以通过-m 参数来设定ORF的最小长度，默认是100
TransDecoder.Predict -t merged.transcripts.fa #预测编码区
cdna_alignment_orf_to_genome_orf.pl merged.transcripts.fa.transdecoder.gff3 merged.transcripts.gff3 merged.transcripts.fa >merged.transcripts.fa.transdecoder.genome.gff3 #生成基于参考基因组的gff3注释文件
cd ..
echo "transdecoder done"
## 结果merged.transcripts.fa.transdecoder.cds为最终ORF对应的cds序列。可以提供给EDTA做基因组重复序列预测的筛选步骤。
## 结果merged.transcripts.fa.transdecoder.gff3为最终ORF对应的gff3文件。
## 结果merged.transcripts.fa.transdecoder.pep为最终ORF对应的氨基酸pep序列。
## 结果merged.transcripts.fa.transdecoder.bed用于后期的IGV可视化，以BED格式存放ORF位置信息
## 结果merged.transcripts.fa.transdecoder.genome.gff3为基于参考基因组的gff3注释文件，可以提供给基因组基因预测软件做基因组注释。

# 【optional】后面两步骤选用，基于单个转录组数据做组装。
## 2plus.stringtie组装单个转录本
stringtie -p ${thread} -o 2.stringtie/rna_1.gtf 1.hisat2/rna_1.bam &
stringtie -p ${thread} -o 2.stringtie/rna_2.gtf 1.hisat2/rna_2.bam &
stringtie -p ${thread} -o 2.stringtie/rna_3.gtf 1.hisat2/rna_3.bam &
stringtie -p ${thread} -o 2.stringtie/rna_4.gtf 1.hisat2/rna_4.bam &
wait
echo "stringtie single done"

## 3plus.gffread转换单个转录本gtf格式成fa和gff格式
gffread -w 3.gffread/rna_1.fa -g ${genome} 2.stringtie/rna_1.gtf &
gffread -w 3.gffread/rna_2.fa -g ${genome} 2.stringtie/rna_2.gtf &
gffread -w 3.gffread/rna_3.fa -g ${genome} 2.stringtie/rna_3.gtf &
gffread -w 3.gffread/rna_4.fa -g ${genome} 2.stringtie/rna_4.gtf &

gffread -E 2.stringtie/rna_1.gtf -o - |sed -e "s#transcript#match#g" -e "s#exon#match_part#g" > 3.gffread/rna_1.gff &
gffread -E 2.stringtie/rna_2.gtf -o - |sed -e "s#transcript#match#g" -e "s#exon#match_part#g" > 3.gffread/rna_2.gff &
gffread -E 2.stringtie/rna_3.gtf -o - |sed -e "s#transcript#match#g" -e "s#exon#match_part#g" > 3.gffread/rna_3.gff &
gffread -E 2.stringtie/rna_4.gtf -o - |sed -e "s#transcript#match#g" -e "s#exon#match_part#g" > 3.gffread/rna_4.gff &
wait
echo "gffread single done"
