#!bin/bash

#mkdir RoachCelegans
#mkdir RoachCelegans/rawdata

# manually prepare download info file for fastq and fast5


PROJ=/blue/mcintyre/share/transcript_distance


#form=fast5
form=fastq
touch $PROJ/data/CelegansRoach2020/md5_downloaded_$form.txt
for i in `awk -F "\t" 'NR!=1 && $0!="" {print $10"|"$11}' $PROJ/data/CelegansRoach2020/${form}_dataInfo_RoachCelegans.xls`
do
    md5=$(echo $i | awk -F "|" '{print $1}')
    ad=$(echo $i | awk -F "|" '{print $2}')
    fname=$(basename $ad)
    echo ${fname}"	"${md5} >> $PROJ/data/CelegansRoach2020/md5_downloaded_$form.txt
    wget -c $ad -O $PROJ/data/RoachCelegans2020/$form/${fname} -t 30
done

md5sum $PROJ/data/CelegansRoach2020/*.fastq.gz | awk -F "[\t/]" '{print $1"\t"$NF}'> $PROJ/data/CelegansRoach2020/md5sum.result.txt

awk 'BEGIN{print "fileName""\t""md5_calculated_fromFile""\t""md5_downloaded"} NR==FNR{a[$1]=$2} NR>FNR{print $2"\t"$1"\t"a[$2]}' \
    $PROJ/data/CelegansRoach2020/md5_downloaded_fastq.txt $PROJ/data/CelegansRoach2020/md5sum.result.txt > \
    $PROJ/data/CelegansRoach2020/md5_compare.txt

rm $PROJ/data/CelegansRoach2020/md5_downloaded_fastq.txt
rm $PROJ/data/CelegansRoach2020/md5sum.result.txt



