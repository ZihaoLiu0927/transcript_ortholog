#!/bin/bash

mkdir WangMaize
mkdir WangMaize/rawdata

# manually prepare file /ufgi-serrata/data/zihaoliu/WangMaize/dataInfo_WangMaize.xls

for i in `awk -F "\t" 'NR!=1 && $0!="" {print $29"|"$31}' /ufgi-serrata/data/zihaoliu//WangMaize/dataInfo_WangMaize.xls`
do
    dir=$(echo $i | awk -F "|" '{print $2}')
    file=$(echo $i | awk -F "|" '{print $1}')
    ad="ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR326/${dir}/${file}"
    md5="ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR326/${dir}/${dir}.md5"
    #echo $ad
    wget -c $ad -O /ufgi-serrata/data/zihaoliu/WangMaize/rawdata/${file} -t 30
    wget -c $md5 -O /ufgi-serrata/data/zihaoliu/WangMaize/rawdata/${dir}.md5 -t 30
done


p=$(pwd)
mv ${p}/nohup.out /ufgi-serrata/data/zihaoliu/WangMaize/nohup_WangMaize_log.out

cat /ufgi-serrata/data/zihaoliu//WangMaize/rawdata/*.md5 | grep -P ".*(.bam)$" > /ufgi-serrata/data/zihaoliu/WangMaize/md5sum.download.txt
md5sum /ufgi-serrata/data/zihaoliu/WangMaize/rawdata/md5sum *.bam > /ufgi-serrata/data/zihaoliu/WangMaize/md5sum.result.txt

awk 'BEGIN{print "fileName""\t""md5_calculated_fromFile""\t""md5_downloaded"} NR==FNR{a[$2]=$1} NR>FNR{print $2"\t"$1"\t"a[$2]}' \
    /ufgi-serrata/data/zihaoliu/WangMaize/md5sum.download.txt /ufgi-serrata/data/zihaoliu/WangMaize/md5sum.result.txt > \
    /ufgi-serrata/data/zihaoliu/WangMaize/md5_compare.txt

rm /ufgi-serrata/data/zihaoliu/WangMaize/md5sum.download.txt
rm /ufgi-serrata/data/zihaoliu/WangMaize/md5sum.result.txt
