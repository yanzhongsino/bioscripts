#/****************************************
#	filename     :extract_4dtv_from_cds.sh
#	date         :2021/10/07 16:51:00
#	author       :Yan Zhong
#	email        :yan.zhong.sino@gmail.com
#	description  :从根据codon比对的多条cds序列中提取4倍简并位点(fourfold degenerate codons, 4dtv)。
#*****************************************/

#!/usr/bin/bash

# extract_4dtv_from_cds.sh脚本是用于从根据codon比对的多条cds序列中提取4倍简并位点(fourfold degenerate codons, 4dtv)。

## 00-数据和变量定义
cds=/path/to/cds.fa
species=4dtvSpecies #以哪个物种的4dtv为标准提取

##### 提取的4dtv结果保存在4dtv.fa #####
# 以下三种方案选择一种使用即可
## 01-生成运行文件
seqkit locate -V 0 -i -d -p GCN -p CGN -p GGN -p CTN -p CCN -p TCN -p ACN -p GTN scg.m.pep2cds.fa |awk -v awka="${species}" '{if ($1 == "awka" && $6%3== 0) print $6}'|sort -k 1n |uniq |awk -v awkb="${cds}" '{print "<(seqkit subseq -r "$1":"$1" "awkb")"}' |sed -e '1i\seqkit concat' -e '$a\> 4dtv.fa'|sed -E ":a;N;s/\n/ /g;ta" >4dtv.sh #生成从cds获取4dtv的命令，很长，所以储存在文件4dtv.sh中。运行4dtv.sh就可以获取4dtv.fa文件。

sh 4dtv.sh #运行生成的文件

## 02-分布运行
### 由于一步生成运行文件占用内存，运行较慢，可以分步骤进行。
seqkit locate -V 0 -i -d -p GCN -p CGN -p GGN -p CTN -p CCN -p TCN -p ACN -p GTN scg.m.pep2cds.fa |awk -v awka="${species}" '{if ($1 == "awka" && $6%3== 0) print $6}'|sort -k 1n |uniq |awk -v awkb="${cds}" '{print "<(seqkit subseq -r "$1":"$1" "awkb")"}' >4dtv.sh 
sed -i -e '1i\seqkit concat' -e '$a\> 4dtv.fa'|sed -E ":a;N;s/\n/ /g;ta" 4dtv.sh

sh 4dtv.sh

## 03-分物种运行
for i in $(echo ${species[*]}); do cat ${species[5]}.4dtv|sed "s/${species[5]}/${i}/g" >${i}.4dtv.bed; done #为每个物种保存一份4dtv.bed文件
for i in $(echo ${species[*]}); do seqkit subseq --bed ${i}.4dtv.bed scg.m.pep2cds.fa | grep -v ">"|sed -E ":a;N;s/\n//g;ta" |sed "1i\>${i}" >${i}.4dtv.fa; done #提取4dtv位点【耗时步骤】
cat *.4dtv.fa >4dtv.fa #合并