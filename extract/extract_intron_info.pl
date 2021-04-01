#!/usr/bin/perl
# use strict;
use warnings;

# usage: perl extract_intron_info.pl input.gff3 >intron.info
# 从gff3注释文件中提取intron的位置信息和长度，每行六列，1.chr; 2.gene_ID; 3.position_left; 4.position_right; 5.intron_length; 当基因无内含子时，3-5列写为0。
# ref: http://bioops.info/2012/11/intron-size-gff3-perl/

my $input=$ARGV[0];
my ($eachline,@exons);
my ($eachline_8,$desc,$gene,$ctg);
my $first=0;
open (IN, "<$input") or die ("no such file!");
while(defined($eachline=<IN>)){
  if($eachline=~/\tmRNA\t/){
    $first++;
    if($first != 1){
      print_intron(@exons);
      @exons=();
    }
  }elsif($eachline=~/\tCDS\t/){
    my @eachline=split(/\t/,$eachline);
    my $eachline_8=$eachline[8];
    my @desc=split(/\;/,$eachline_8);
    $gene=$desc[0];
    $ctg=$eachline[0];
    #    print "$ctg\t";
    #    print "$gene\t\n";

    push (@exons,$eachline[3],$eachline[4]);
  }
}

#print "$ctg\t";
#print "$gene\t\n";

print_intron(@exons);


sub print_intron{
  my (@exons)=@_;
  if(scalar(@exons)>2){
    my @ordered_exons=sort {$a<=>$b} @exons;
    for (my $i=1;$i<=scalar(@ordered_exons)-3;$i=$i+2){
      my $each_intron_size=$ordered_exons[$i+1]-$ordered_exons[$i]-1;
      my $intron_l=$ordered_exons[$i]+1;
      my $intron_r=$ordered_exons[$i+1]-1;
      print "$ctg\t$gene\t";
      print "$intron_l\t$intron_r\t$each_intron_size\n";
    }
  }else{print "$ctg\t$gene\t0\t0\t0\n";}
  print "\n";
}

