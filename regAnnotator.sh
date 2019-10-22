#!/bin/bash

usage() { printf "\n"; echo "Usage: $0 [-s <sample name>] [-d <input directory>]" 1>&2; printf "\n"; exit 1; }

while getopts ":s:d:" opt; do
  case $opt in
    s) sample="$OPTARG"
    ;;
    d) outdir="$OPTARG"
    ;;
    \?) usage
    exit
    ;;
  esac
done
shift $((OPTIND-1))

if [ -z "${sample}" ] || [ -z "${outdir}" ]; then
    usage
fi

mkdir -p $sample"/logs/"
samout=$(basename $sample)

getgffpath=$outdir"/bin/get_NCBI_gff.py"
parseregionspath=$outdir"/bin/parse_regions.py"
sumsamplepath=$outdir"/bin/summarize_sample.py"

annotation()
{
cat $sample"/"$samout".cov" | sed 's/"//g' | grep ">gi" | sed -e $'s/\|/\\\t/g' | awk '{print $4}' > $sample"/"$samout"_nclist.txt"
cat $sample"/"$samout".cov" | sed 's/"//g' | grep ">gi" | cut -f1 -d ' ' | sed 's/>//' | awk '{print $1}' > $sample"/"$samout"_gilist.txt"

python $getgffpath -i $sample"/"$samout"_nclist.txt" > $sample"/logs/"$samout"_get_ncbi_gff.log"


grep -v "^>" $sample"/"$samout".gt4c" | cut -d'|' -f 4,5 | perl -pe 's/\|.*?\t/\t/' | cut -f1-3 > $sample"/"$samout".txt"
cat $sample"/"$samout".txt" | awk -v var="$sample/$samout" '{print > (var"_"$1".bed"); close(var"_"$1".bed")}' RS=""


rm $sample"/"$samout*"_res.bed" || true;

for x in $sample"/"$samout*".bed";
do
NCID=$(echo $x | cut -d"_" -f2,3 | sed "s/.bed//");
$outdir"/bin/"bedmap --echo --echo-map-id --multidelim '\t' $x $outdir"/NCBI-GFF/"$NCID"_parsed.sort.bed" | tr '|' $'\t' > ${x%.bed}"_res.bed";
done

cat $sample"/"$samout*"_res.bed" > $sample"/"$samout"_res_bed.txt"
rm $sample"/"$samout*".bed" || true;

python $parseregionspath -i $sample"/"$samout"_res_bed.txt" -c $sample"/"$samout".gt4c"

cat $sample"/"$samout"_region.bed" | sort -t "|" -k 5 | uniq > $sample"/"$samout"_regions.list"
grep -v "repeat_region" $sample"/"$samout"_regions.list" > $sample"/"$samout"_genes.list"
rm $sample"/"$samout*"bed"  $sample"/"$samout*".txt" || true;

python $sumsamplepath -d $sample
rm $sample"/"$samout*".list" || true

}

if annotation ; then
    echo "Annotation Complete"
    date
else
    echo "Failed: Annotation"
    exit 1
fi

