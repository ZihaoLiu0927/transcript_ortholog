#!/bin/bash
#SBATCH --mail-user=zihaoliu@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=download_data_MaizeWang.%A_%a
#SBATCH--output=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/download_data_MaizeWang.%A_%a.out
#SBATCH--error=/blue/mcintyre/share/transcript_distance/analysis/SLURM_LOGS/download_data_MaizeWang.%A_%a.err 
#SBATCH --cpus-per-task=12
#SBATCH --time=99:00:00
#SBATCH --ntasks=1
#SBATCH --mem=3gb
 


PROJ=/blue/mcintyre/share/transcript_analysis

mkdir -p $PROJ/data
mkdir -p $PROJ/data/MaizeWang2020

# manually prepare file $PROJ/data/MaizeWang2020/dataInfo_WangMaize.xls

for i in `awk -F "\t" 'NR!=1 && $0!="" {print $29"|"$31}' $PROJ/data/WangMaize2020/dataInfo_WangMaize.xls`
do
    dir=$(echo $i | awk -F "|" '{print $2}')
    file=$(echo $i | awk -F "|" '{print $1}')
    ad="ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR326/${dir}/${file}"
    md5="ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR326/${dir}/${dir}.md5"
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
