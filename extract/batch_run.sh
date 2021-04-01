#!/bin/bash
# 并行运行循环指令，并控制并行数量为${thread}
# usage: 修改thread=后的值为想要设定的并行线程数，把需要并行的命令替换23行的run something。运行脚本batch_run.sh，实现并行运行循环任务。
# ref：网上找的脚本，来源找不到了，侵删。

thread=16

start_time=`date +%s`              # 定义脚本运行的开始时间
echo ${start_time}
tmpFifo=/tmp/$$.fifo # 声明管道名称，$$表示脚本当前运行的进程PID
mkfifo ${tmpFifo} # 创建有名管道
exec 3<>${tmpFifo}                   #创建文件描述符，以可读（<）可写（>）的方式关联管道文件，这时候文件描述符3就有了有名管道文件的所有特性
rm -rf ${tmpFifo}                    #关联后的文件描述符拥有管道文件的所有特性,所以这时候管道文件可以删除，我们留下文件描述符来用就可以了
for ((i=1;i<=${thread};i++))
do
        echo "" >&3                   #&3代表引用文件描述符3，这条命令代表往管道里面放入了一个"令牌"
done

for sample in $(cat sample_list.txt) #并行运行的循环
do
        read -u3                           #代表从管道中读取一个令牌
        {
                echo ${sample}
                run something # 需要并行运行的真实命令
                echo 'success' ${sample}       
                echo "" >&3                   #代表我这一次命令执行到最后，把令牌放回管道
        } &
done

wait
 
stop_time=`date +%s`  #定义脚本运行的结束时间
 
echo "TIME:`expr $stop_time - $start_time`"
exec 3<&-                       #关闭文件描述符的读
exec 3>&-                       #关闭文件描述符的写
