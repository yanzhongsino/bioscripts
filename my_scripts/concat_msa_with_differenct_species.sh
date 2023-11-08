#/****************************************
#	filename     :concat_msa_with_different_species.sh
#	date         :2023/11/08
#	author       :Yan Zhong
#	email        :yan.zhong.sino@gmail.com
#	description  :the script's information
#*****************************************/

#!/usr/bin/bash

# concat_msa_with_different_species.sh脚本用于合并拥有不同数量的物种的多序列比对(mutiple sequence alignment, MSA)，使得合并后保持原有比对不变，缺乏基因序列的位置用横杠-替代。
# 待合并的单个多序列比对文件保存在input文件夹中，脚本逻辑是先统计共有多少物种，给每个多序列比对文件添加缺少的物种和用横杠-代表的空序列用来占位，输出到out/full/目录下，后缀是_full.fa。再用seqkit concat来合并序列，合并后的结果输出到out目录下all_concat.fa。


mkdir out
mkdir out/full
seqkit seq -n input/* |sort|uniq >./out/all.list # 统计完整的物种list
for i in $(ls ./input/ ) # 为每个msa序列添加缺乏的物种
do
	seqkit seq -n ./input/$i >./out/full/${i}.list # 统计msa序列的物种list
	grep -v -f ./out/full/${i}.list ./out/all.list >./out/full/${i}_v.list # 统计msa序列缺乏的物种list
	num=$(seqkit stats input/$i |awk '{print $8}'|tail -1|sed "s/,//g") # 获取msa序列的长度
	cp ./input/$i ./out/full/${i}_full.fa # 复制一份msa序列
	for j in $(cat ./out/full/${i}_v.list) # 添加每个msa序列缺乏的物种
	do
		echo ">$j" >> ./out/full/${i}_full.fa # 添加缺乏的物种名
		awk -v awka="${num}" 'BEGIN{while(n++ < awka)printf "-";print""}' >> ./out/full/${i}_full.fa # 添加msa序列长度的-(代表gap)，以保持新添加物种和原有msa序列的对齐
	done
	rm ./out/full/${i}_v.list
	rm ./out/full/${i}.list
done
rm ./out/all.list
seqkit concat -w 0 ./out/full/*_full.fa >./out/all_concat.fa # 安装序列ID合并所有含共同物种的msa序列到一个文件