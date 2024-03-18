#!/bin/bash
set -eou pipefail

showHelp() {
cat << EOF
$0 useage:
	$0 [args] directory
opetional args:
	-s min variant size
	-e extention
EOF
}

[ $# -eq 0 ] && showHelp && exit 1

size=10000
ext=".vcf"
write='true'

while getopts :h:s:e:w arg; do
	case $arg in
		s ) size=$OPTARG;;
		e ) ext=$OPTARG;;
		w ) write='false';;
		h ) showHelp
		    exit 1;;
	esac
done
shift $((OPTIND -1))
echo ${1}

indir=$1

echo "directory: ${indir}, size: ${size}, ext: ${ext}"


for i in ${indir}/*.vcf; 
do 
	b=$(basename ${i} ${ext})
	o=${indir}/size_filtered_${size}/${b}.vcf
	if ${write}
	then
		mkdir -p ${indir}/size_filtered_${size}
		bcftools filter -e "SVLEN<${size}" ${i} > ${o}
		echo "${b} from $(cat ${i} | wc -l) to $(cat ${o} | wc -l)"
	else
		echo "${b} from $(cat ${i} | wc -l) to" $(bcftools filter -e "SVLEN<${size}" ${i} | wc -l)
	fi
done
