#!/usr/bin/bash

##aCGH script annotation of bed files with 5 first columns of interest
##auteur: Christophe Desterke
##april 2017 version 1.0
##dependencies: need to install bedtools package, install databases files in subfolder "DATABASES_HG18/" 
##database file dependencies: hg18refseq_sorted.bed and scanDBceuHG18_sorted.bed and COSMICHG18_sorted.bed and plurinetHG18_sorted.bed and TFHG18_sorted.bed
##usage: sh acgh_subcontrol_HG18.sh example.bed control.bed

variable=${1}
if [ -z "${variable}" ]
then 
	echo "BED file not passed as parameter"
	exit 1
fi
control=${2}
if [ -z "${control}" ]
then 
	echo "BED file not passed as parameter"
	exit 1
fi

log_file="log_file.log" 

nom_fichier=$(echo $1 | sed -re 's/(.*).bed/\1/') 

cat $1 | awk -v OFS="\t" '{print $1,$2,$3,$4,$5}' > file_5c.bed

cat $2 | awk -v OFS="\t" '{print $1,$2,$3,$4,$5}' > origin.bed

sort -k 1,1 -k2,2n file_5c.bed > sample.bed

sort -k 1,1 -k2,2n origin.bed > ct_sorted.bed

bedtools subtract -a sample.bed -b ct_sorted.bed > sorted.bed

bedtools subtract -a sorted.bed -b DATABASES_HG18/scanDBceuHG18_sorted.bed > sub

bedtools intersect -a sub -b DATABASES_HG18/HG18refseq_sorted.bed -wb > annotatedWOcnv.bed

bedtools intersect -a DATABASES_HG18/scanDBceuHG18_sorted.bed -b sorted.bed > cnv.bed

bedtools intersect -a DATABASES_HG18/plurinetHG18_sorted.bed -b sub > pluripotency.bed

bedtools intersect -a DATABASES_HG18/COSMICHG18_sorted.bed -b sub > cancer.bed

bedtools intersect -a DATABASES_HG18/TFHG18_sorted.bed -b sub > TF.bed

bedtools intersect -a ct_sorted.bed -b sample.bed > locict.bed

sort -k 1,1 -k2,2n locict.bed > locict_sorted.bed

bedtools intersect -a locict_sorted.bed -b DATABASES_HG18/HG18refseq_sorted.bed -wb > controlgene.bed 

echo "---------------------------------" >> $log_file
echo "log of the aCGH annotation on HG18 human genome" >> $log_file
echo "analysis of the file = $nom_fichier" >> $log_file
echo "--------------------" >> $log_file
cal 

date >> $log_file

awk '!arr[$9]++' controlgene.bed > controlgeneuniq.bed
echo "Number of genes eliminated by the control :" >> $log_file
wc -l controlgeneuniq.bed >> $log_file

awk -v OFS="\t" '{print $9,$4,$5}' annotatedWOcnv.bed > hg18_annot_genes.tsv
echo "Number of genes annotated :" >> $log_file
wc -l hg18_annot_genes.tsv >> $log_file

awk '!arr[$1]++' hg18_annot_genes.tsv > hg18uniq_genes.tsv
echo "Number of genes uniq :" >> $log_file
wc -l hg18uniq_genes.tsv >> $log_file

echo "Number of polymorphisms eliminated :" >> $log_file
awk '!arr[$4]++' cnv.bed > cnv_uniq.bed
wc -l cnv_uniq.bed >> $log_file


echo "List of polymorphisms eliminated of the annotation : " >> $log_file
cat cnv_uniq.bed >> $log_file

echo "--------------------" >> $log_file
echo "Number of cancer genes :" >> $log_file
awk '!arr[$4]++' cancer.bed > cancer_uniq.bed
wc -l cancer_uniq.bed >> $log_file
cat cancer_uniq.bed >> $log_file

echo "--------------------" >> $log_file
echo "Number of pluripotency genes :" >> $log_file
awk '!arr[$4]++' pluripotency.bed > pluripotency_uniq.bed
wc -l pluripotency_uniq.bed >> $log_file
cat pluripotency_uniq.bed >> $log_file



echo "--------------------" >> $log_file
echo "Number of Transcrption Factors :" >> $log_file
awk '!arr[$4]++' TF.bed > TF_uniq.bed
wc -l TF_uniq.bed >> $log_file
cat TF_uniq.bed >> $log_file



rm cancer.bed
rm TF.bed
rm pluripotency.bed
rm sorted.bed
rm file_5c.bed
rm sub
rm hg18_annot_genes.tsv
rm ct_sorted.bed
rm origin.bed
rm sample.bed
rm locict.bed
rm locict_sorted.bed
rm controlgene.bed

mkdir RESULTS

mv TF_uniq.bed pluripotency_uniq.bed cancer_uniq.bed annotatedWOcnv.bed log_file.log cnv_uniq.bed cnv.bed RESULTS
mv hg18uniq_genes.tsv RESULTS
mv controlgeneuniq.bed RESULTS

cd RESULTS
mkdir CNV_removed
mkdir GENES_SCORES
mkdir ELIMINATED_CONTROL

mv controlgeneuniq.bed ELIMINATED_CONTROL

mv cnv.bed CNV_removed
mv hg18uniq_genes.tsv GENES_SCORES

mv TF_uniq.bed $(echo TF_uniq.bed | sed "s/\./".$nom_fichier"\./")
mv pluripotency_uniq.bed $(echo pluripotency_uniq.bed | sed "s/\./".$nom_fichier"\./")
mv cancer_uniq.bed $(echo cancer_uniq.bed | sed "s/\./".$nom_fichier"\./")
mv annotatedWOcnv.bed $(echo annotatedWOcnv.bed | sed "s/\./".$nom_fichier"\./")
mv cnv_uniq.bed $(echo cnv_uniq.bed | sed "s/\./".$nom_fichier"\./")
cat log_file.log

mv log_file.log $(echo log_file.log | sed "s/\./".$nom_fichier"\./")
cd .. CGHarray
mv RESULTS RESULTS_$nom_fichier
exit 0
