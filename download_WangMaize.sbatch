#!/bin/bash


PROJ=/blue/mcintyre/share/transcript_analysis

mkdir -p $PROJ/data
mkdir -p $PROJ/data/MaizeWang2020

# manually prepare file /ufgi-serrata/data/zihaoliu/WangMaize/dataInfo_WangMaize.xls

for i in `awk -F "\t" 'NR!=1 && $0!="" {print $29"|"$31}' $PROJ/data/WangMaize2020/dataInfo_WangMaize.xls`
do
    dir=$(echo $i | awk -F "|" '{print $2}')
    file=$(echo $i | awk -F "|" '{print $1}')
    ad="ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR326/${dir}/${file}"
    md5="ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR326/${dir}/${dir}.md5"
    #echo $ad
    wget -c $ad -O $PROJ/data/MaizeWang2020/${file} -t 30
    wget -c $md5 -O $PROJ/data/WangMaize2020/${dir}.md5 -t 30
done


cat $PROJ/data/WangMaize2020/*.md5 | grep -P ".*(.bam)$" > $PROJ/data/WangMaize2020/md5sum.download.txt
md5sum $PROJ/data/WangMaize2020/*.bam | awk -F "[\t/]" '{print $1"\t"$NF}' > $PROJ/data/WangMaize2020/md5sum.result.txt

awk 'BEGIN{print "fileName""\t""md5_calculated_fromFile""\t""md5_downloaded"} NR==FNR{a[$2]=$1} NR>FNR{print $2"\t"$1"\t"a[$2]}' \
    $PROJ/data/WangMaize2020/md5sum.download.txt $PROJ/data/WangMaize2020/md5sum.result.txt > \
    $PROJ/data/WangMaize2020/md5_compare.txt

rm $PROJ/data/WangMaize2020/md5sum.download.txt
rm $PROJ/data/WangMaize2020/md5sum.result.txt
