#!/bin/bash

number_rep=$1
min_tr_lenght=$2
min_FPKM=$3
min_isoform_fract=$4
locus_gap=$5
transcript_merge_gap=$6

BAM_DIR=$7

SECONDS=0

###fastqc quality control - requires fastqc installed and placed in PATH

printf "Analyzing combined conditions transcriptome using:\n\n"

cd $BAM_DIR
mkdir GTFs

ls -1 *clean.bam

ls -1 *clean.bam | tr '\n' ' ' > commands.1
sed -i 's/^/bamtools merge /g' commands.1
sed -i 's/$/ -out merged.bam/g' commands.1
chmod 755 commands.1
#./commands.1

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh37_mapping/gencode.v29lift37.annotation.gtf.gz
gunzip gencode.v29lift37.annotation.gtf.gz
mv gencode.v29lift37.annotation.gtf GTFs

find . -name '*clean.bam' | {
    read firstbam
    samtools view -h "$firstbam"
    while read bam; do
        samtools view "$bam"
    done
} | samtools view -ubS - | samtools sort - -m 100000000000 -o merged.bam
samtools index merged.bam
ls -l merged.bam merged.bam.bai

echo "Stringtie analysis of combined conditions.."
echo stringtie -l MergedBAM -G GTFs/gencode.v29lift37.annotation.gtf -p 64 -m $min_tr_lenght -f $min_isoform_fract -g $locus_gap -j 10 -c 15 -o mergedBam.gtf merged.bam > commands.2.sh
chmod 755 commands.2.sh
./commands.2.sh

mv mergedBam.gtf GTFs
#find . -name mergedBam.gtf.tmp | xargs -I % sh -c 'echo %;grep chr %;' > mergedBam.gtf.tmp.gtf
#mv mergedBam.gtf.tmp.gtf GTFs/mergedBam.gtf

duration=$SECONDS
echo "Stringtie analysis of combined conditions: $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed." > mergedBam.tm


SECONDS=0

echo "Stringtie analysis of per condition transcriptome.."
number_rep=$1
for (( i=1; i<=${number_rep} ; i+=1 )) ; do
ls -1d *[$i]_clean.bam > $i.tmp.txt
done

paste -d " " *tmp.txt > per.condition.sh

echo "Analyzing per condition transcriptome using:"
cat per.condition.sh

cut -f1 -d " " per.condition.sh | sed 's/_clean.bam//g'| sed -E 's/[-0-9]//g' > merged.name

cut -f1 -d " " per.condition.sh | sed 's/_clean.bam//g'| sed -E 's/[-0-9]//g' > per.condition.sh.tmp
sed -i 's/$/.merged.bam/g' per.condition.sh.tmp
paste -d " " per.condition.sh.tmp per.condition.sh > per.condition.sh.tmp.2
sed -i 's/^/samtools merge /g' per.condition.sh.tmp.2
mv per.condition.sh.tmp.2 per.condition.sh

sed 's/$/.merged.sort.bam/g' merged.name>merged.sort.name
sed -i 's/^/samtools sort -o /g' merged.sort.name
paste -d " " merged.sort.name per.condition.sh.tmp >> per.condition.sh

sed 's/samtools sort -o /samtools index /g' merged.sort.name >> per.condition.sh
sed -i 's/samtools sort -o //g' merged.sort.name

sed 's/$/.merged.perCond.gtf/g' merged.name > merged.name.gtf
sed -i 's/^/-p 64 -m '${min_tr_lenght}' -f '${min_isoform_fract}' -g '${locus_gap}' -o GTFs\//g' merged.name.gtf 
sed 's/^/stringtie -G GTFs\/gencode.v29lift37.annotation.gtf -l /g' merged.name > begin.tmp

paste -d " " begin.tmp merged.name.gtf merged.sort.name >> per.condition.sh


rm per.condition.sh.tmp
rm merged.sort.name
rm merged.name.gtf
rm merged.name
rm *tmp.txt
rm *.tmp

chmod 755 per.condition.sh
./per.condition.sh 

duration=$SECONDS
echo "Stringtie analysis of per condition transcriptome: $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed." > percondBam.tm

SECONDS=0

printf "Downloading reference annotation..\n\n"
#wget ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/10.0/genome_coordinates/gff3/homo_sapiens.GRCh38.gff3.gz
#gunzip homo_sapiens.GRCh38.gff3.gz


printf "Merging GTFs: combined condition and condition specific..\n\n"


#mv homo_sapiens.GRCh38.gff3 GTFs
cd GTFs

stringtie --merge -l CombGTF -p 64 -G gencode.v29lift37.annotation.gtf -m $min_tr_lenght -F $min_FPKM -f $min_isoform_fract -g $transcript_merge_gap -i -o ./mergedGtf.gtf *merged.perCond.gtf
stringtie --merge -l CombGTF -p 64 -G gencode.v29lift37.annotation.gtf -m $min_tr_lenght -F $min_FPKM -f $min_isoform_fract -g $transcript_merge_gap -i -o ./CombGTF.gtf mergedBam.gtf mergedGtf.gtf 

printf "Number of discovered transcripts in combined conditions transcriptome:\n\n"
grep "	transcript	"  mergedBam.gtf | wc -l

printf "Number of discovered transcripts in per condition transcriptome:\n\n"
grep "	transcript	"  mergedGtf.gtf | wc -l

printf "Number of discovered transcripts in merged GTF - per condition and combined conditions :\n\n"
grep "	transcript	"  CombGTF.gtf | wc -l


cat mergedBam.tm
cat percondBam.tm

