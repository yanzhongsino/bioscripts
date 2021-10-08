#/****************************************
#	filename     :extract_single_copy_genes.sh
#	date         :2021/10/07 17:02:11
#	author       :Yan Zhong
#	email        :yan.zhong.sino@gmail.com
#	description  :the script's information
#*****************************************/

#!/usr/bin/bash

# extract_single_copy_genes.sh脚本用于做完orthofinder的分析得到单拷贝直系同源基因之后提取各个物种的单拷贝直系同源基因。

## 00-1 数据准备
### 1. 单拷贝直系同源基因的信息保存在scg.ortho文件中，内容是singlecopyorthologues列表，第一列orthogroupID，后面所有列每个物种的基因序列ID，共12列，这个文件的每行数据为一个homologs单拷贝基因对应的每套数据的基因序列ID。
### 3. 所有物种的cds序列和pep序列保存在cds和pep目录下，并以物种名+cds.fa/pep.fa命名文件,序列都为一行的格式，否则应该用`seqkit seq -w 0 cds.fa > new.cds.fa`预先处理。

data_path=/path/to/data #cds和pep目录存放的目录
ortho=/path/to/scg.ortho #scg.ortho文件的绝对路径
species=(A B C D E F G) #物种名保存在species的list中，用于后面调用。物种顺序与scg.ortho文件中物种顺序一致。
desc=(AT BB C_34 d343_LL ecs F gfffe) #cds和pep序列id前缀，与species顺序一致。
thread=8 #定义提取步骤的并行运行数量

## 00-2 函数定义

#### extract_sequences用于提取orthogroup的序列
function extract_sequences(){
    a=$(echo $1)
    for i in $(seq 1 ${#species[@]});do j=`expr ${i} - 1`; grep -A 1 ${sample[${i}]} ${data_path}/$a/${species[${j}]}.$a.fa |sed "s/${sample[${i}]}.*/${species[${j}]}/" >./scg/${sample[0]}/${species[${j}]}.$a; done
    cat ./scg/${sample[0]}/*.$a > ./scg/${sample[0]}/${sample[0]}.$a.fa
    rm ./scg/${sample[0]}/*.$a
}

#### merge_sequences用于合并比对后的序列
function merge_sequences(){
    a=$(echo $1)
    b=$(echo $2)
    for i in $(echo ${species[*]}); do cat ./scg/OG*/OG*.${a}${b}fa |grep -A 1 ">${i}" |grep -v ">"|sed -E ":a;N;s/\n//g;ta" |sed "s/ //g"|sed "1i\>${i}" >./scgs/${i}.${b}${a}fa; done
    cat ./scgs/*.${b}${a}fa >./scgs/scg.${a}${b}fa
    rm ./scgs/*.${c}${b}${a}fa
}



## 00-3 为循环准备
start_time=`date +%s`              # 定义脚本运行的开始时间
echo 'start time:' ${start_time}
tmpFifo=/tmp/$$.fifo # 声明管道名称，$$表示脚本当前运行的进程PID
mkfifo ${tmpFifo} # 创建有名管道
exec 3<>${tmpFifo}                   #创建文件描述符，以可读（<）可写（>）的方式关联管道文件，这时候文件描述符3就有了有名管道文件的所有特性
rm -rf ${tmpFifo}                    #关联后的文件描述符拥有管道文件的所有特性,所以这时候管道文件可以删除，我们留下文件描述符来用就可以了
for ((i=1;i<=${thread};i++))
do
        echo "" >&3                   #&3代表引用文件描述符3，这条命令代表往管道里面放入了一个"令牌"
done

## 01 提取和比对orthogroup
echo "start extract step"
mkdir ./scg
mkdir ./scgs

while read line
do
        read -u3                           #代表从管道中读取一个令牌
        {
                echo "start ${sample[0]}"

                # 提取singlecopyorthologs.txt的每列信息保存到sample这个list中。
                for i in $(seq 0 ${#species[@]});do j=`expr ${i} + 1`; sample[${i}]=$(echo $line |awk  -v awka="$j" '{print $awka}');done

                # 提取orthogroup的cds和pep序列
	        mkdir ./scg/${sample[0]}

                extract_sequences cds
                extract_sequences pep

                # cds序列align，mafft和prank
                mafft --auto --thread 2 ./scg/${sample[0]}/${sample[0]}.cds.fa |seqkit seq -w 0 -u > ./scg/${sample[0]}/${sample[0]}.mafft.cds.fa

                prank -d=./scg/${sample[0]}/${sample[0]}.cds.fa -o=./scg/${sample[0]}/${sample[0]}.cds -codon
	        seqkit seq -w 0 -u ./scg/${sample[0]}/${sample[0]}.cds.best.fas > ./scg/${sample[0]}/${sample[0]}.prank.cds.fa
               

                # pep序列align，mafft和prank，并用pal2nal转成cds的align结果
                mafft --auto --thread 2 ./scg/${sample[0]}/${sample[0]}.pep.fa |seqkit seq -w 0 -u > ./scg/${sample[0]}/${sample[0]}.mafft.pep.fa

	        prank -d=./scg/${sample[0]}/${sample[0]}.pep.fa -o=./scg/${sample[0]}/${sample[0]}.pep -protein
	        seqkit seq -w 0 -u ./scg/${sample[0]}/${sample[0]}.pep.best.fas > ./scg/${sample[0]}/${sample[0]}.prank.pep.fa

                pal2nal.pl ./scg/${sample[0]}/${sample[0]}.mafft.pep.fa ./scg/${sample[0]}/${sample[0]}.cds.fa -output fasta |seqkit seq -w 0 -u >./scg/${sample[0]}/${sample[0]}.mafft.p2c.fa
                pal2nal.pl ./scg/${sample[0]}/${sample[0]}.prank.pep.fa ./scg/${sample[0]}/${sample[0]}.cds.fa -output fasta |seqkit seq -w 0 -u >./scg/${sample[0]}/${sample[0]}.prank.p2c.fa
                pal2nal.pl -nogap -nomismatch ./scg/${sample[0]}/${sample[0]}.mafft.pep.fa ./scg/${sample[0]}/${sample[0]}.cds.fa -output fasta |seqkit seq -w 0 -u >./scg/${sample[0]}/${sample[0]}.mafft.p2c_trim.fa
                pal2nal.pl -nogap -nomismatch ./scg/${sample[0]}/${sample[0]}.prank.pep.fa ./scg/${sample[0]}/${sample[0]}.cds.fa -output fasta |seqkit seq -w 0 -u >./scg/${sample[0]}/${sample[0]}.prank.p2c_trim.fa
                
                echo "success ${sample[0]}"       
                echo "" >&3                   #代表我这一次命令执行到最后，把令牌放回管道

        } &
done < ${ortho}

### 删除中间文件
rm ./scg/${sample[0]}/${sample[0]}.cds.best.fas
rm ./scg/${sample[0]}/${sample[0]}.pep.best.fas
rm -rf ./tmpdirprank* #删除prank运行产生的临时目录

wait

sleep 60s

exec 3<&-                       #关闭文件描述符的读
exec 3>&-                       #关闭文件描述符的写
echo "extract step end"

## 02 合并比对序列
echo "start merge step"

### 合并比对的cds和pep序列
merge_sequences mafft. cds.
merge_sequences mafft. pep.
merge_sequences prank. cds.
merge_sequences prank. pep.

#### 合并根据pep比对结果转换的CDS比对序列
merge_sequences mafft. p2c.
merge_sequences prank. p2c.
merge_sequences mafft. p2c_trim.
merge_sequences prank. p2c_trim.


echo "-------------------------** all done**--------------------"
echo "merge step end"

stop_time=`date +%s`  #定义脚本运行的结束时间
echo 'stop_time:' ${stop_time}
echo "运行用时 total time:`expr $stop_time - $start_time`"