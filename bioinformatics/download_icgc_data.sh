

function download_icgc {
	release=$1
	project=$2
	for e in simple_somatic_mutation.open donor sample specimen
	do
		infile="https://dcc.icgc.org/api/v1/download?fn=/${release}/Projects/${project}/${e}.${project}.tsv.gz"
		outfile=`basename $infile`
		echo $infile
		echo $outfile
		if [ ! -e $outfile ];then
	  	 wget $infile -O $outfile
		fi
	done
}

while read project
do
if [ ! -e $project ];then
	mkdir $project
fi
cd $project
download_icgc release_25 $project
cd ../

done < cancer_types_release_25.txt


## cancer_types.txt looks like:
#COCA-CN
#ESCA-CN

