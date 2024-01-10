---
title: "Pasteurallaceae: Haemophilus and Aggregatibacter human oral metapangenome"
date: '2022-06-28'
output: rmdformats::readthedown
---



# 1. Set up 

##### Login

The following was used to log into a secure computer called Barhal housed at The Marine Biological Laboratory in Woods Hole, MA. All work for this analysis was executed either on Barhal or a laptop (MacBook Pro 14-inch, 2021) owned and operated by J.J. Giacomini. Please note that I omited the actual USERNAME for privacy reasons. 

```{bash, eval=FALSE}

# Choose one of the methods to log into Barhal computer
ssh USERNAME
ssh barhal-01

# Alternative option for using anvi interactive
ssh -L 8081:localhost:8081 USERNAME
ssh -L 8081:localhost:8081 barhal-01

# Alternative option for using anvi interactive
ssh -L 8080:localhost:8080 USERNAME
ssh -L 8080:localhost:8080 barhal-01

# Load anvio 
module purge
module load miniconda/3
source /bioware/miniconda3/bashrc
conda activate anvio-7.1
module load clusters/barhal
module load jbpc
```

##### Define variables

```{bash, eval=FALSE}
projectID=P_0622_Haemophilus_Aggregatibacter
mainDIR=/workspace/jmarkwelchlab
meta=$mainDIR/$projectID/DATA/Haemophilus_Aggregatibacter_Filtered.csv
ITEMS=${projectID}-add_info.items.txt
seqidHOMD=/workspace/jmarkwelchlab/HOMD_INFO/SEQID_info.txt
DIR_SITE=/storage/data/00_Public_Metagenomes/oral/HMP/
hmpID=/workspace/jmarkwelchlab/HMP_METAGENOMES/METADATA/samples_id-QC.txt
hmpMetadata=/workspace/jmarkwelchlab/HMP_METAGENOMES/METADATA/sequencing-metadata-count-QC.tsv
minContigSIZE=300
TH=15

# directories (variables)
DIR_Data=DATA
DIR_Reads=00_READS
DIR_NCBI=01_NCBI_GENOMES
DIR_MEREN=01_MEREN_GENOMES
DIR_Assemblies=02_ASSEMBLIES
DIR_Contigs=03_GENOMES_EDITED
DIR_ContigsDB=04_CONTIGS_DB
DIR_Mapping=05_MAPPING
DIR_SinglePROF=06_SINGLE_PROFILE
DIR_MergedPROF=07_MERGED_PROFILE
DIR_SummaryPROF=08_PROFILE_SUMMARY
DIR_Pangenome=09_PANGENOME
DIR_SummaryPAN=10_PANGENOME_SUMMARY
DIR_Gene_calls=11_GENE_CALLS
DIR_Annotation=12_FUNCTIONAL_ANNOTATION
DIR_DetectionGENOMES=13_DETECTED_GENOMES
DIR_Phylo=14_PHYLOGENOMICS
DIR_Derep=15_DE_REPLICATION
DIR_var=22_GENE_LEVEL
DIR_CheckM=26_CHECKM

# files (variables)
genomeMetadata=$meta
projectMetadata=$DIR_Data/02_${projectID}.csv
downloadNCBI=$DIR_Data/03_${projectID}-download_url.txt
rawNCBIid=$DIR_Data/04_${projectID}-ncbi_raw_genomes_id.txt
genusList=$DIR_Data/05_${projectID}-genus_list.txt
checkNCBIid=$DIR_Data/06_${projectID}-check_genomes_id.txt
newHOMDid=$DIR_Data/07_${projectID}-homd_id.txt
genomesNCBIid=$DIR_Data/08_${projectID}-ncbi_genomes_id.txt
duplicateStrainNames=$DIR_Data/09_${projectID}-duplicate_strain_names.txt
duplicateStrainInfo=$DIR_Data/10_${projectID}-duplicate_strain_names_info.txt
removeDuplicate=$DIR_Data/11_${projectID}-duplicate_strain_removed.txt
magsID=$DIR_Data/12_mags_id-meren.txt
nameConversions=$DIR_Data/13_${projectID}-name_conversions.txt
genomesID=$DIR_Data/id_genomes.txt
genomesALL=$DIR_Data/id_genomes-ALL.txt
genomesRefSeq=$DIR_Data/id_genomes-RefSeq.txt
RAW_ASSEMBLY=$DIR_Assemblies/${projectID}-RAW.fa
PROJECT_CONTIGS=$DIR_Contigs/${projectID}.fa
PROJECT_REPORT=$DIR_Contigs/${projectID}.report.tsv
BINNING=$DIR_Contigs/${projectID}.binning.tsv
DECOMPOSE=$DIR_Contigs/${projectID}.decompose.tsv
CONTIGS_DB=$DIR_ContigsDB/${projectID}-contigs.db
GENE_CALLS=$DIR_Gene_calls/${projectID}-gene_calls.fa
PAN_LAYERS=$DIR_Pangenome/${projectID}-add_info.layers.tsv
samplesMetadata=DATA/${projectID}-samples_metadata.txt
LAYER_ORDERS=DATA/${projectID}-layer_orders.txt

genomesALL=DATA/id_genomes-ALL.txt
nameConversions=$DIR_Data/13_${projectID}-name_conversions.txt
iDir=$DIR_Pangenome/internal_annotated_${projectID}
genomeDB=$iDir/${projectID}-GENOMES.db
panDB_detection=$iDir/${projectID}-RESULTS/${projectID}-PAN_detection.db
genomes_98ANI=DATA/id_genomes-98ANI.txt

PanProject=${projectID}-98ANI
PanDir=$DIR_Pangenome/$PanProject
GENOMES_98ANI_DEREP=$PanDir/${PanProject}.txt
GENOMES_98ANI_DEREP_DB=$PanDir/${PanProject}-GENOMES.db
PAN_DIR=$PanDir/${PanProject}-RESULTS
PAN_DB=$PanDir/${PanProject}-RESULTS/${PanProject}-PAN.db
ANI_DIR=$PanDir/${PanProject}-RESULTS/ANI_RESULTS
layersADD_98ANI=$PanDir/${PanProject}-layers_98ANI.tsv
PAN_DES=DATA/description_98ANI_derep_pangenome.txt
```


##### Make directories

```{bash, eval=FALSE}

mkdir $mainDIR/$projectID && cd $mainDIR/$projectID

mkdir $DIR_Data $DIR_Reads $DIR_NCBI $DIR_MEREN $DIR_Assemblies $DIR_Contigs $DIR_ContigsDB $DIR_Mapping $DIR_SinglePROF $DIR_MergedPROF $DIR_SummaryPROF $DIR_Pangenome $DIR_SummaryPAN $DIR_DetectionGENOMES $DIR_Gene_calls $DIR_Annotation $DIR_Phylo $DIR_Derep $DIR_CheckM
```

# 2. Download and process genomes

Family: Pasteurallacaea

Genera: Haemophilus & Aggregattibacter

Species listed in NCBI: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=416916
1. Aggregatibacter actinomycetemcomitans (serotypes a,b,c,d,e,f)
2. Aggregatibacter aphrophilus  
3. Aggregatibacter kilianii (https://pubmed.ncbi.nlm.nih.gov/29695522/)
4. Aggregatibacter pneumotropica (not validly published; identified in rats from MAGS using partial 16S rRNA genes)
5. Aggregatibacter segnis  
6. unclassified Aggregatibacter (A. sp.)

7. Haemophilus aegyptius (valid)
8. Haemophilus haemolyticus (valid)
9. Haemophilus influenzae (valid)
10. Haemophilus influenzae-murium (Unknown; from cats and dogs) 
11. Haemophilus massiliensis (valid)
12. Haemophilus paracuniculus (valid; from rabbits)
13. Haemophilus parahaemolyticus (valid)
14. Haemophilus parainfluenzae (valid)
15. Haemophilus paraphrohaemolyticus (valid)   
16. Haemophilus pittmaniae (valid; https://lpsn.dsmz.de/species/haemophilus-pittmaniae)  
17. Haemophilus quentini (not valid; https://lpsn.dsmz.de/species/haemophilus-quentini)
18. Haemophilus seminalis (valid; https://lpsn.dsmz.de/species/haemophilus-seminalis) 
19. Haemophilus simiae (not valid; unpublished study; isolated from monkey bite wound)
20. Haemophilus sputorum (valid;https://lpsn.dsmz.de/species/haemophilus-sputorum)
21. Haemophilus taxon C (not sure what this is; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC205807/?page=9)
22. Haemophilus ducreyi (valid; distantly related)
23. unclassified Haemophilus  (H. sp.)


We used the NCBI metadata file for all prokaryotes (downloaded from https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/) and filtered for our target taxa as needed. 


The HOMD metadata file can be downloaded via wget -q -P /workspace/jmarkwelchlab/SPECIES_LEVEL_PANGENOMES/DATA/HOMD_INFO http://www.homd.org/ftp/genomes/PROKKA/current/SEQID_info.txt


```{bash, eval=FALSE}
# upload NCBI metadata "prokCSV" file onto the server 
scp -r /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/DATA/Haemophilus_Aggregatibacter_Filtered.csv USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/DATA/Haemophilus_Aggregatibacter_Filtered.csv


# copy NCBI metadata for desired genomes (RefSeq)
cat $genomeMetadata | awk 'BEGIN{FS=OFS=","}NR==1{print $0}NR>1{print $0}' > $projectMetadata

# make list with gca path for download; GCA download link is in column
cat $genomeMetadata | awk -F',' -v OFS="\t" 'NR>1{print $18}' | awk 'BEGIN{FS=OFS="/"}{print $0,$NF"_genomic.fna.gz"}' | sed -e 's/"//g' > $downloadNCBI

# Download genomes
for gca_assembly in `cat $downloadNCBI`
do
gcaFile=$( echo $gca_assembly | awk -F'/' '{print $NF}' | sed -e 's/\.gz//' )
if [ ! -f $DIR_NCBI/${gcaFile} ]
then
echo "File not found - Downloading: ${gcaFile} "
wget -q -P $DIR_NCBI $gca_assembly
gunzip $DIR_NCBI/${gcaFile}.gz
echo "unzipped: ${gcaFile} "
fi
done


# Check the number of files in the directory. 
# This should be the number of genomes you expected to be downloaded, which in our case is 1226
ls $DIR_NCBI | wc -l #1226

```

##### Rename genomes

Anvio doesn't like complicated genome IDs so we have to create a list of unique IDs containing a number and the assembly ID separated by "_id_"; use GCA assembly ID column 9

```{bash, eval=FALSE}
# create a list of unique IDs containing a number and the assembly ID separated by "_id_"; use GCA assembly ID column 9
# (e.g. G_0001_id_GCA_000000000.0)
cat $projectMetadata | awk -F',' -v OFS="\t" 'NR>1{ printf "G_""%04i_id_%s\n", NR-1,$9 }' | sed -e 's/"//g'  > $rawNCBIid

# copy
cp $rawNCBIid $genomesNCBIid
cp $rawNCBIid $genomesID

# create new fasta files for each genome and 
# change names of original fasta file to number a number using the $genomesID file
# (e.g. GCA_901873365.1_Aggregatibacter_sp._BgEED05_genomic.fna -> G_0001-RAW.fa)
while IFS= read -r genomes_id
do
gca_id=$(echo $genomes_id | awk -F'_id_' '{print $2}')
old_name=$(find $DIR_NCBI -name "$gca_id*")
NEW_NAME=$(echo "$genomes_id" | awk -F'_id_' -v new_dir="$DIR_Assemblies" '{print new_dir"/"$1"-RAW.fa"}' )
cp $old_name $NEW_NAME
done < $genomesID

# remove GCA number from $genomesID file to match new fasta file names
sed -i 's/_id_.*//' $genomesID

```

##### Reformat genomes 

Using ANVIO (anvi-script-reformat-fasta) to change any non-canonical letter with N and simplify deflines. This is done in batches of 20 in a background process.

```{bash, eval=FALSE}

N=20
(
for genome in `cat $genomesID`
do
((i=i%N)); ((i++==0)) && wait
EDITED_GENOMES=$DIR_Assemblies/${genome}-RAW.fa
CONTIGS=$DIR_Contigs/${genome}.fa
REPORT=$DIR_Contigs/${genome}.report.tsv
anvi-script-reformat-fasta -l $minContigSIZE -o $CONTIGS --simplify-names --prefix ${genome} --seq-type NT -r $REPORT $EDITED_GENOMES &
done
)
```

##### Create Contig dbs

```{bash, eval=FALSE}
N=20
(
for genome in `cat $genomesID`
do
((i=i%N)); ((i++==0)) && wait
CONTIGS=$DIR_Contigs/${genome}.fa
contigsDB=$DIR_ContigsDB/${genome}-contigs.db
numThreads=5
anvi-gen-contigs-database -f $CONTIGS -n ${genome} -o $contigsDB -T $numThreads &
done
)

```

##### Annotate Contig dbs with marker genes

Marker genes: Bacterial_71, rRNAs Arch and Bac

```{bash, eval=FALSE}
N=20
(
for genome in `cat $genomesID`
do
((i=i%N)); ((i++==0)) && wait
contigsDB=$DIR_ContigsDB/${genome}-contigs.db
numThreads=5
anvi-run-hmms -c $contigsDB -T $numThreads &
done
)

# OPTIONAL
# search for tRNAs in genomes (scan-trnas)
N=20
(
for genome in `cat $genomesID`
do
((i=i%N)); ((i++==0)) && wait
contigsDB=$DIR_ContigsDB/${genome}-contigs.db
numThreads=5
anvi-scan-trnas -c $contigsDB -T $numThreads &
done
)

```


##### Contig db stats

This is an important step to make sure that each and every contig has the same hmms annotations. Otherwise there will be issues downstream. 

```{bash, eval=FALSE}

mkdir 19_Contig_db_stats
contigsDB=$DIR_ContigsDB

# run this to make sure it works
#anvi-display-contigs-stats -o 19_Contig_db_stats/G_0001-contigs.db-stats  $DIR_ContigsDB/G_0001-contigs.db

for contig in `ls $contigsDB`
  do
      anvi-display-contigs-stats  --report-as-text -o 19_Contig_db_stats/$contig-stats  $DIR_ContigsDB/$contig
done

# check each contig-stats file 
echo -e "contig_db\tBacteria_71" > Contig_bacteria_71_summary.txt

for contig in `ls 19_Contig_db_stats`
  do
    echo "##### STARTING WITH $contig #####"
    cat 19_Contig_db_stats/$contig | grep 'contigs_db\|Bacteria_71' | sed 'n; n; d' > tmp
    contig_db=$(awk 'NR==1{print $2}' tmp)
    Bacteria_71=$(awk 'NR==2{print $2}' tmp)
    echo -e "$contig_db\t$Bacteria_71" >> Contig_bacteria_71_summary.txt
    rm tmp
    echo "##### DONE WITH $contig #####"
done

```

##### Rename genomes with HOMD names

```{bash, eval=FALSE}
# get conversion names
while IFS= read -r line
do
genomeNo=$( echo $line | awk -F'_id_' '{print $1}')
assemblyID=$( echo $line | awk -F'_id_' '{print $2}')
putGenoName=$( grep "$assemblyID" $projectMetadata | awk -F',' -v OFS="\t" '{print $1,$6,$9}' | sed -e 's/"//g' | sed -e 's/ /_/g' | sed -e 's/uncultured_//' | sed -e 's/\[//' | sed -e 's/\]//' | awk -F'\t' -v OFS="|" '{print $1,$2,$3}' | awk 'BEGIN{FS=OFS="_"}{ $1= substr($1,1,1)}1' | sed -e 's/_sp\./_sp/' | sed -E 's/(.*)\./\1POINT/' | sed -e "s/[^[:alnum:]|]/_/g" | sed -E 's/(.*)POINT/\1\./' | sed -e 's/_//' | sed -e 's/|/;/' |  sed -e 's/_.*;/;/' | sed -e 's/;/_str_/' | sed -e 's/|/_id_/' | sed -e 's/^.\{1\}/&_/' | sed -e 's/_str__id_/_str_NA_id_/' | sed -e 's/\./_/')
echo -e "$putGenoName\t$genomeNo" >> $nameConversions
done < $rawNCBIid

#find HOMD ID
for hypo_hmt in `cat $rawNCBIid | awk -F'_id_' '{print $2}' `
do
cat $seqidHOMD | grep "$hypo_hmt" | awk -v gca_id="$hypo_hmt" 'BEGIN{FS=OFS="\t"}{print $3"_"$4"_str_"$5"_id_"gca_id}' | sed -e's/ //g' | sed -e 's/_sp\./_sp_/' | awk 'BEGIN{FS=OFS="_"}{ $1= substr($1,1,1)}1' | awk '{if ($1 ~/HMT/) gsub("_sp","");print}' | awk -F'_id_' -v OFS="\t" '{print $NF,$0}' | sed -e 's/-/_/g' | grep 'HMT' | sed -e 's/HMT/sp_HMT_/' >> $newHOMDid
done

# edit name conversion with HOMD sps.
while IFS= read -r line
do
gcaID=$( echo $line | awk '{print $1}' | sed -e 's/\./_/' )
oldID=$( grep "$gcaID" $nameConversions | awk '{print $1}' )
newID=$( echo $line | awk '{print $2}' | sed -e 's/\./_/' )
sed -i "s/$oldID/$newID/" $nameConversions
done < $newHOMDid

```

##### Make genome metadata 

```{bash, eval=FALSE}
# copy ITEMS info for each genome
# short version
echo -e "item\tGenome_ID\tSpecies\tRefSeq\tAssembly_ID\tAssembly_ID_norm\tStrain\tGenome_in_HOMD\tHOMD_ID\tG_ID" > $ITEMS
for gcaID in `cat $nameConversions | cut -f1 | awk -F'_id_' '{print $2}' | rev | sed -e 's/_/./' | rev `
do
gcaID_ed=$(echo $gcaID | sed -e 's/\./_/')
itemsID=$(grep "$gcaID_ed" $nameConversions | cut -f1 | awk -F'_str_' -v OFS="\t" '{print $0,$0,$1}')
strainID=$( grep "$gcaID_ed" $nameConversions | cut -f1 | awk -F'_str_' -v OFS="\t" '{print $2}' | awk -F'_id_' -v OFS="\t" '{print $1}' )
refSeq=$(grep "$gcaID" $projectMetadata | awk -F',' -v OFS="\t"  '{if ($16!="") print "Yes" ; else if ($16=="") print "No"}')
genomeHOMD=$( grep "$gcaID" $seqidHOMD)
if [ -z "$genomeHOMD" ]
then
inHOMD="no"
hmtID="NA"
else
inHOMD="yes"
hmtID=$( echo "$genomeHOMD"| awk 'BEGIN{FS=OFS="\t"}{print "HMT_"$2}')
fi
gID=$( grep "$gcaID_ed" $nameConversions | cut -f2 )
echo -e "$itemsID\t$refSeq\t$gcaID_ed\t$gcaID\t$strainID\t$inHOMD\t$hmtID\t$gID" >> $ITEMS
done


```

##### Contig db paths

```{bash, eval=FALSE}

GENOMES=DATA/$projectID-contig_paths.txt

echo -e 'name\tcontigs_db_path' > $GENOMES
for genome in `cat $genomesID`
do
cat $nameConversions | grep "$genome" | awk -v workPath="$PWD" 'BEGIN{FS=OFS="\t"}{print $1,workPath"/04_CONTIGS_DB/"$2"-contigs.db"}' >> $GENOMES
done


# check; should be 1227 rows; sometimes the code above produces duplicate rows; not sure why
cat $GENOMES | wc -l # 1229
# Looks like we have some duplicates
cat $GENOMES | awk '{print $1}' | sort | uniq -c | 2
#2 A_aphrophilus_str_CCUG_11575_id_GCA_003130375_1
#2 H_haemolyticus_str_CCUG_11096_id_GCA_006439155_1
cat $GENOMES | grep A_aphrophilus_str_CCUG_11575 #  genome ID 1036
cat $GENOMES | grep H_haemolyticus_str_CCUG_11096_id_GCA_006439155_1 # genome ID 0084
# manually removed them in nano
# check fuinal number 
cat $GENOMES | wc -l # 1227
***
```

# 3. CheckM Completness and Contamination

```{bash, eval=FALSE}
# make a directory for genome fasta files
mkdir $DIR_CheckM/BINS_DIR

# generate new nameConversions file for dereplicated genome set
nameConversions_98DEREP=DATA/13_P_0622_Haemophilus_Aggregatibacter-name_conversions-98DEREP.txt
grep -f $genomes_98ANI $nameConversions > $nameConversions_98DEREP

# copy genomes and change names $rawNCBIid
while IFS= read -r line
do
orginalName=$( echo "$line" | awk -F'\t' '{print $2}')
newName=$( echo "$line" | awk -F'\t' '{print $1}')
cp $DIR_Contigs/$orginalName.fa $DIR_CheckM/BINS_DIR/$newName.fa
done < $nameConversions_98DEREP

# taxonomy output with genomes
# check taxonomic contamination and completeness of genomes
binDir=$DIR_CheckM/BINS_DIR
outDIR=$DIR_CheckM/${projectID}-taxonomic
outFile=$outDIR/${projectID}
rank=family
taxon=Pasteurellaceae

#checkm taxon_list (for CPR using )
module load python/3.7.9-202101061157
checkm taxon_set $rank $taxon $outFile.lineage.ms.txt
#Marker set for Pasteurellaceae contains 767 marker genes arranged in 440 sets.
#Marker set inferred from 83 reference genomes.
checkm analyze -x fa -t 100 $outFile.lineage.ms.txt $binDir $outDIR
checkm qa -o 1 -f $outFile.lineage.qa.txt --tab_table -t 100 $outFile.lineage.ms.txt $outDIR

cat $DIR_CheckM/${projectID}-taxonomic/${projectID}.lineage.qa.txt | awk -F "\t" '{print $1,$12, $13}'
```

# 4. Dereplication 

##### ANIb method (98% ANI using pyANI)

To dereplicate genomes based on ANI = or > 98% identity, I used anvi-dereplicate-genomes with --program pyANI --method ANIb --similarity-threshold 0.98 --cluster-method simple_greedy --representative-method centrality.

```{bash, eval=FALSE}

clusterize -n 16 -m jgiacomini@forsyth.org -l Haem_Agg_derep_98ANI.log ./Dereplication.sh

#!/bin/bash

projectID=P_0622_Haemophilus_Aggregatibacter
GENOMES=DATA/$projectID-contig_paths.txt
DIR_Derep=15_DE_REPLICATION

module load blast+

# deduplicate Qscore (98%)
anvi-dereplicate-genomes -e $GENOMES -o $DIR_Derep/ANI_98_dereplication --skip-fasta-report --program pyANI --method ANIb --similarity-threshold 0.98 --cluster-method simple_greedy --representative-method centrality --num-threads 16 --log-file anvi-dereplicate.log


```


##### fastANI (98% ANI)

*** need to make sure fastANI verison 1.3 is loaded ***

```{bash, eval=FALSE}

module unload fastani/1.1
module load fastani/1.32
clusterize -n 8 -m jgiacomini@forsyth.org -l Haem_Agg_derep_98ANI_fastANI.log ./Dereplication-fastANI.sh

#!/bin/bash

projectID=P_0622_Haemophilus_Aggregatibacter
GENOMES=DATA/$projectID-contig_paths.txt
DIR_Derep=15_DE_REPLICATION

# deduplicate Qscore (98%)
anvi-dereplicate-genomes -e $GENOMES -o $DIR_Derep/ANI_98_dereplication_fastANI --skip-fasta-report --program fastANI --similarity-threshold 0.98 --cluster-method simple_greedy --representative-method centrality --num-threads 8 --log-file anvi-dereplicate-fastANI.log

```



Now we need to select the isolate reference genomes that we will use for mapping of the metagenomes. There are a few different methods we can use to do this. First, we can use ANI to identify clusters of genomes and unique strains. We can then manually select genomes that are type strains, or genomes that we are particularly interested in. 

H. influenzae 
ATCC_33391 CCUG_23945 CIP_102514 DSM_4690 NCTC8143 NCTC_8143
cat 15_DE_REPLICATION/ANI_98_dereplication_fastANI/CLUSTER_REPORT.txt | grep H_influenzae_str_NCTC8143_id_GCA_001457655_1

A_actinomycetemcomitans
ATCC 33384; CCUG 13227; CIP 52.106; DSM 8324; NCTC 9710
cat 15_DE_REPLICATION/ANI_98_dereplication_fastANI/CLUSTER_REPORT.txt | grep A_actinomycetemcomitans_str_NCTC_9710_id_GCA_008085305_1

H_parainfluenzae
cat 15_DE_REPLICATION/ANI_98_dereplication_fastANI/CLUSTER_REPORT.txt | grep H_parainfluenzae_str_NCTC 7857

A_aphrophilus
ATCC 33389; CCUG 3715; CIP 70.73; NCTC 5906
cat 15_DE_REPLICATION/ANI_98_dereplication_fastANI/CLUSTER_REPORT.txt | grep A_aphrophilus_str_NCTC5906

A_segnis
ATCC 33393; CCUG 10787; CCUG 12838; CIP 103292; DSM 21418; HK316; NCTC 10977
cat 15_DE_REPLICATION/ANI_98_dereplication_fastANI/CLUSTER_REPORT.txt | grep A_segnis_str_NCTC10977

H_paraphrohaemolyticus_str_NCTC10671_id_GCA_900451065_1

H_parahaemolyticus_str_NCTC_8479_id_GCA_900450675_1

A_kilianii
CCUG 70536; DSM 105094; PN_528
A_kilianii_str_PN_528_id_GCA_003130255_1

##### Send results to local

```{bash, eval=FALSE}
scp -r /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/15_DE_REPLICATION/ANI_98_dereplication_fastANI/CLUSTER_REPORT_Type_strains.txt USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/15_DE_REPLICATION/ANI_98_dereplication_fastANI/CLUSTER_REPORT_Type_strains.txt

##### pyANI all 202 dereplicated genomes
#mkdir /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/ANI

scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/09_PANGENOME/P_0622_Haemophilus_Aggregatibacter-98ANI/P_0622_Haemophilus_Aggregatibacter-98ANI-RESULTS/ANI_RESULTS/ /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/ANI
```

##### fastANI figure

fastANI all 1226 genomes. Send data to local machine.

```{bash, eval=FALSE}
scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/15_DE_REPLICATION/ANI_98_dereplication_fastANI /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/15_DE_REPLICATION/
```

Build plot using fastANI data for all 1226 genomes
```{r, eval=FALSE}
# packages
library(ggplot2)
library(tidyr)
library(reshape2)
library(phytools)
library(tibble)

library(dendextend)
library(ape)
library(dplyr)

# load ANI data
df <-read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/15_DE_REPLICATION/ANI_98_dereplication_fastANI/SIMILARITY_SCORES/fastANI_ani.txt", header = TRUE)

# load ANI tree
df_newick <-  ape::read.tree("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/15_DE_REPLICATION/ANI_98_dereplication_fastANI/SIMILARITY_SCORES/fastANI_ani.newick")


# order ANI data by ANI tree
df_newick_orders <- phytools::compute.mr(df_newick, type = "matrix")
df_newick_rownames <- rev(rownames(df_newick_orders))
data_ordered <- df[ order(match(df$key, df_newick_rownames)), ]
df_newick_rownames_df <- as.data.frame(df_newick_rownames)
cols_A<- ncol(df) -1
df_newick_rownames_df2 <- as.data.frame(matrix(df_newick_rownames_df$df_newick_rownames, ncol = cols_A, byrow = TRUE))
names(df_newick_rownames_df2) <- df_newick_rownames_df2[1,]
df_newick_rownames_df3 <- df_newick_rownames_df2[-1,]
df_newick_rownames_df3 <- df_newick_rownames_df3 %>% 
  add_column(key = NA, .before = 1)
data_ordered2<-data_ordered[names(df_newick_rownames_df3)]

cols2<- ncol(data_ordered2) 
long_df <- reshape2::melt(data_ordered2,
                          id.vars=c("key"),
                          measure.vars=colnames(data_ordered2[2:cols2]),
                          variable.name="y",
                          value.name="z")
mylevels1 <- df_newick$tip.label
long_df$key <- factor(long_df$key,levels=mylevels1)
long_df$y <- factor(long_df$y, levels=mylevels1)



# plot
plot<-ggplot(long_df, aes(key,y)) +
  geom_tile(aes(fill = z)) + 
  scale_fill_gradient2(low = "darkblue",
                       mid = "white",
                       high = "darkred",
                       midpoint = 0.95,
                       limits =c(0.90, 1))+
  #geom_text(aes(label = format(round(z, digits=3), nsmall = 3)),size=2.75) +
  ylab("") +  
  xlab("") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks.y = element_line(size = 1),
        axis.ticks.x=element_line(size = 1)) 

rows <- as.data.frame(df_newick$tip.label)
height=nrow(rows) *0.15
width=height*1.2
  
ggsave(file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/15_DE_REPLICATION/ANI_98_dereplication_fastANI/SIMILARITY_SCORES/ANI_heatmap.pdf", plot = plot, width = width, height = height, limitsize = FALSE)

```



##### ANIb figure
```{r, eval=FALSE}
# packages
library(ggplot2)
library(tidyr)
library(reshape2)
library(phytools)
library(tibble)

library(dendextend)
library(ape)
library(dplyr)

# load ANI data
df <-read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/ANI/ANI_RESULTS/ANIb_percentage_identity.txt", header = TRUE)

# load ANI tree
df_newick <-  ape::read.tree("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/ANI/ANI_RESULTS/ANIb_percentage_identity.newick")

# plot tree
plot(df_newick)

# order ANI data by ANI tree
df_newick_orders <- phytools::compute.mr(df_newick, type = "matrix")
df_newick_rownames <- rev(rownames(df_newick_orders))
data_ordered <- df[ order(match(df$key, df_newick_rownames)), ]
df_newick_rownames_df <- as.data.frame(df_newick_rownames)
cols_A<- ncol(df) -1
df_newick_rownames_df2 <- as.data.frame(matrix(df_newick_rownames_df$df_newick_rownames, ncol = cols_A, byrow = TRUE))
names(df_newick_rownames_df2) <- df_newick_rownames_df2[1,]
df_newick_rownames_df3 <- df_newick_rownames_df2[-1,]
df_newick_rownames_df3 <- df_newick_rownames_df3 %>% 
  add_column(key = NA, .before = 1)
data_ordered2<-data_ordered[names(df_newick_rownames_df3)]

cols2<- ncol(data_ordered2) 
long_df <- reshape2::melt(data_ordered2,
                          id.vars=c("key"),
                          measure.vars=colnames(data_ordered2[2:cols2]),
                          variable.name="y",
                          value.name="z")
mylevels1 <- df_newick$tip.label
long_df$key <- factor(long_df$key,levels=mylevels1)
long_df$y <- factor(long_df$y, levels=mylevels1)



# plot
plot<-ggplot(long_df, aes(key,y)) +
  geom_tile(aes(fill = z)) + 
  scale_fill_gradient2(low = "darkblue",
                       mid = "white",
                       high = "darkred",
                       midpoint = 0.95,
                       limits =c(0.90, 1))+
  #geom_text(aes(label = format(round(z, digits=3), nsmall = 3)),size=2.75) +
  ylab("") +  
  xlab("") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks.y = element_line(size = 1),
        axis.ticks.x=element_line(size = 1)) 

rows <- as.data.frame(df_newick$tip.label)
height=nrow(rows) *0.15
width=height*1.2
  
ggsave(file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/ANI/ANI_RESULTS/202_ANI_heatmap.pdf", plot = plot, width = width, height = height, limitsize = FALSE)

```


Make ANI tree for figure
```{r, eval=FALSE}
# load ANI tree
ANI_newick <-  ape::read.tree("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/ANI/ANI_RESULTS/ANIb_percentage_identity.newick")

# convert to dedrogram
ANI_dend <- ape::chronos(ANI_newick)
ANI_dend <- as.dendrogram(ANI_dend)


library(ggplot2)
library(ggdendro)

#convert cluster object to use with ggplot
dendr <- dendro_data(ANI_dend, type="rectangle") 

#your own labels (now rownames) are supplied in geom_text() and label=label
ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), size=0.25) + 
  geom_text(data=label(dendr), aes(x=x, y=y, label=label, hjust=0), size=0.65) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())

ggsave("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/FIGURES/02_ANI/ANI_dendrogram.pdf", height = 6, width = 5)
```

##### Select genomes

```{bash, eval=FALSE}

for gcaID in `cat $DIR_Derep/ANI_98_dereplication_fastANI/CLUSTER_REPORT_Type_strains.txt | awk 'NR>1{print $3}' | awk -F'_id_' '{print $2}' | rev | sed -e 's/_/\./' | rev `
do
cat $rawNCBIid | grep "$gcaID" | awk -F'_id_' '{print $1}' >> $genomesID.deduplicated
done

# order unique genomes by number
sort -u $genomesID.deduplicated > $genomesID.deduplicated.ordered

# check number of genomes (2002 expected...202 total)
wc -l $genomesID.deduplicated.ordered

#I want the full genome names for the dereplicated genomes
cat $genomesID.deduplicated.ordered | xargs -I@ -n1 grep @ $nameConversions | tee $genomesID.deduplicated_full.txt

# check list of genome IDs
cat $genomesID.deduplicated_full.txt | awk '{print $1}'


# extract genome number for deduplicated strains/species
# only include genomes for mapping
for genomeID in `cat $genomesID.deduplicated_full.txt | awk '{print $1}'`
do
cat $nameConversions | grep -w "$genomeID" | awk '{print $2}' >> $genomesID.derep98
done

```

##### Get genome count per species

```{bash, eval=FALSE}

echo -e "species\tcount" > Genomes_98ANI_derep_count.txt

for group in `cat $layersADD_98ANI | awk 'NR>1{print $3}' | sort -u`
do
      species=$(echo $group)
      count=$(cat $GENOMES_98ANI_DEREP | grep $species | wc -l)
      echo -e "$species\t$count" >> Genomes_98ANI_derep_count.txt
done
```

##### Concatenate and reformat fastas

```{bash, eval=FALSE}
# make concatenated genome file for METAPANGENOME
for genome in `cat $genomesID.derep98`
do
CONTIGS=$DIR_Contigs/${genome}.fa
cat $CONTIGS >> $RAW_ASSEMBLY
done

# Reformat concatenated-genomes using ANVIO (anvi-script-reformat-fasta)
anvi-script-reformat-fasta -l $minContigSIZE -o $PROJECT_CONTIGS --simplify-names --prefix ${projectID} -r $PROJECT_REPORT $RAW_ASSEMBLY
#Total num contigs ............................: 7,334
#Total num nucleotides ........................: 399,374,820

# Binning and decompose-file 
cat $PROJECT_REPORT | awk -F'_' -v OFS="_" 'NF{--NF};1' > $BINNING

while IFS= read -r line
do
intGenomeID=$( echo "$line" | awk '{print $2}' )
newBinID=$( echo "$line" | awk '{print $1}' )
grep "$intGenomeID" $BINNING | sed -e "s/$intGenomeID/$newBinID/" >> $DECOMPOSE
done < $nameConversions

```

Run the following commands on the concatenated fasta file you just made to make sure everything looks ok. We expect each sequence ID to have a prefix for a genome and suffix for a contig. 

```{bash, eval=FALSE}
# print number of sequences
grep -c "^>" $RAW_ASSEMBLY

# To extract ids for each sequence, just use the following:
grep -o -E "^>\w+" $RAW_ASSEMBLY | tr -d ">"

# get table of nucleotide distribution
echo -e "seq_id\tA\tU\tG\tC"; while read line; do echo $line | grep ">" | sed 's/>//g'; for i in A U G C;do echo $line | grep -v ">" | grep -o $i | wc -l | grep -v "^0"; done; done < $RAW_ASSEMBLY | paste - - - - -

```


# 5. Build initial pangenome

##### Make contig paths file

```{bash, eval=FALSE}
# Make txt file that contains genome contig file name and directory  paths for each of the genomes in the set.
genomes_98ANI=DATA/id_genomes-98ANI.txt
awk '{print $1}' $genomesID.deduplicated_full.txt > $genomes_98ANI

PanProject=${projectID}-98ANI
PanDir=$DIR_Pangenome/$PanProject
#mkdir $PanDir
GENOMES_98ANI_DEREP=$PanDir/${PanProject}.txt
layersADD_98ANI=$PanDir/${PanProject}-layers_98ANI.tsv
echo -e 'name\tcontigs_db_path' > $GENOMES_98ANI_DEREP
for genomeTypeID in `cat $genomes_98ANI`
do
cat $nameConversions | grep "$genomeTypeID" | awk -v workPath="$PWD" 'BEGIN{FS=OFS="\t"}{print $1,workPath"/04_CONTIGS_DB/"$2"-contigs.db"}' >> $GENOMES_98ANI_DEREP
done


# make meta data file to add to pangenome
head -n1 $ITEMS > $layersADD_98ANI
for genomeTypeID in `cat $genomes_98ANI`
do
cat $ITEMS | grep "$genomeTypeID" >> $layersADD_98ANI
done

# write pangenome description txt file for anvi-pan-genome --description flag
PAN_DES=DATA/description_98ANI_derep_pangenome.txt

echo "Pangenome for de-replicated set of Haemophilus and Aggregatibacter genomes (n = 202) based on 98% ANI threshold" > $PAN_DES
```


##### Run pangenome.sh script

Build pangenome for the de-replicated set of genomes using pangenome.sh script.

```{bash, eval=FALSE}

clusterize -n 20 -m jgiacomini@mbl.edu -l Haem_Agg_98ANI_derep_pangenome.log ./pangenome.sh

#!/bin/bash

export BLASTDB_LMDB_MAP_SIZE=1000000000
                             
projectID=P_0622_Haemophilus_Aggregatibacter
mainDIR=/workspace/jmarkwelchlab
DIR_Pangenome=$mainDIR/$projectID/09_PANGENOME

PanProject=${projectID}-98ANI
PanDir=$DIR_Pangenome/$PanProject
GENOMES_98ANI_DEREP=$PanDir/${PanProject}.txt
GENOMES_98ANI_DEREP_DB=$PanDir/${PanProject}-GENOMES.db
PAN_DIR=$PanDir/${PanProject}-RESULTS
PAN_DB=$PanDir/${PanProject}-RESULTS/${PanProject}-PAN.db
ANI_DIR=$PanDir/${PanProject}-RESULTS/ANI_RESULTS
layersADD_98ANI=$PanDir/${PanProject}-layers_98ANI.tsv
PAN_DES=DATA/description_98ANI_derep_pangenome.txt
TH=20

# external Genome storage
anvi-gen-genomes-storage -e $GENOMES_98ANI_DEREP -o $GENOMES_98ANI_DEREP_DB

# external pangenome
anvi-pan-genome -g $GENOMES_98ANI_DEREP_DB --use-ncbi-blast --align-with muscle  --minbit 0.5 --mcl-inflation 10 -n ${PanProject} -o $PAN_DIR --num-threads $TH --enforce-hierarchical-clustering --description $PAN_DES

# add LAYERS
anvi-import-misc-data -t layers -p $PAN_DB $layersADD_98ANI

# external genome similarity
anvi-compute-genome-similarity -e $GENOMES_98ANI_DEREP -o $ANI_DIR -p $PAN_DB --program pyANI --method ANIb --num-threads $TH --log-file pangenome_pyANIlog


```

##### Plot pangenome

```{bash, eval=FALSE}
# send pangenome to local machine

scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/09_PANGENOME/ /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/09_PANGENOME/

projectID=P_0622_Haemophilus_Aggregatibacter
DIR_Pangenome=/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/09_PANGENOME/09_PANGENOME
PanProject=${projectID}-98ANI
PanDir=$DIR_Pangenome/$PanProject
GENOMES_98ANI_DEREP_DB=$PanDir/${PanProject}-GENOMES.db
PAN_DB=$PanDir/${PanProject}-RESULTS/${PanProject}-PAN.db

anvi-display-pan -g $GENOMES_98ANI_DEREP_DB -p $PAN_DB --title 98ANI_Un_annotated_pangenome

```

# 6. Phylogeny - Bacteria_71 SCGs

For dereplicated reference genome set (N = 202).

For phylogenomic analyses we used a collection of bacteria single-copy core genes curated by Anvi'o devellpers and based on searches of ∼12000 bacterial high-quality genomes from (https://doi.org/10.1093/bioinformatics/btz188) and modified to remove Ribosomal_S20p, PseudoU_synth_1, Exonuc_VII_S, 5-FTHF_cyc-lig, YidD and Peptidase_A8 from the collection (as they were exceptionally redundant or rare among MAGs from various habitats), and added Ribosomal_S3_C, Ribosomal_L5, Ribosomal_L2 to make it more compatible with Hug et al's (https://www.nature.com/articles/nmicrobiol201648) set of ribosomal proteins.

To generate phylogenetic trees we used the Anvi'o program ‘anvi-get-sequences-for-hmm-hits’ using parameters (1) ‘--align-with muscle’ to perform alignment of protein sequences using MUSCLE (doi:10.1093/nar/gkh340 http://www.drive5.com/muscle), (2) ‘--concatenate-genes’ to concatenate separately aligned and concatenated gene sequences, (3) ‘--returnbest-hit’ to return only the most significant hit when a single HMM profile had multiple hits in one genome, (4) `--get-aa-sequences’ to output amino-acid sequence, and (5) ‘--hmmsources Bacteria_71’ to use the Bacteria 71 collection to search for genes. We used ‘--min-num-bins-gene-occurs’ to ensure that only genes that occur in at least 50% of the genomes were used for the analysis. All 71 genes from the Bacteria 71 collection were subsequrntly used for the phylogenetic tree analysis. We trimmed alignments using trimAl (Capella-Gutiérrez, Silla Martínez & Gabaldón, 2009) with the setting ‘-gt 0.5’ to remove all positions that were gaps in more than 50% of sequences, and computed maximum likelihood phylogenetic trees using IQ-TREE (Nguyen et al., 2015) with the ‘WAG’ general matrix model (Whelan & Goldman, 2001). To root phylogenetic trees we used an type strain genome for E. coli (type strain ATCC 11775; GCA_003697165.2).


##### Add outgroup

```{bash, eval=FALSE}

wget -q -P $DIR_Phylo/REDO https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/697/165/GCA_003697165.2_ASM369716v2/GCA_003697165.2_ASM369716v2_genomic.fna.gz
gunzip $DIR_Phylo/REDO/GCA_003697165.2_ASM369716v2_genomic.fna.gz

mv $DIR_Phylo/REDO/GCA_003697165.2_ASM369716v2_genomic.fna $DIR_Phylo/REDO/E_coli_str_ATCC_11775_id_GCA_003697165_2.fa

#reformat fasta file (simplify deflines)
anvi-script-reformat-fasta -l $minContigSIZE -o $DIR_Phylo/REDO/E_coli_str_ATCC_11775_id_GCA_003697165_2_reformatted.fa --simplify-names --prefix E_coli --seq-type NT -r M_micronuciformis $DIR_Phylo/REDO/E_coli_str_ATCC_11775_id_GCA_003697165_2.fa

# make contigs db
anvi-gen-contigs-database -f $DIR_Phylo/REDO/E_coli_str_ATCC_11775_id_GCA_003697165_2_reformatted.fa -o $DIR_Phylo/REDO/E_coli_str_ATCC_11775_id_GCA_003697165_2.db -T 4

#annotate marker genes
anvi-run-hmms -c $DIR_Phylo/REDO/E_coli_str_ATCC_11775_id_GCA_003697165_2.db
anvi-db-info $DIR_Phylo/REDO/E_coli_str_ATCC_11775_id_GCA_003697165_2.db

GENOMES_98ANI_DEREP_OUTGROUP=/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/09_PANGENOME/P_0622_Haemophilus_Aggregatibacter-98ANI/P_0622_Haemophilus_Aggregatibacter-98ANI_OUTGROUP.txt
cp $GENOMES_98ANI_DEREP $GENOMES_98ANI_DEREP_OUTGROUP

# Add outgroup to bac_71 fasta file
echo -e 'E_coli_str_ATCC_11775_id_GCA_003697165_2\t/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/14_PHYLOGENOMICS/REDO/E_coli_str_ATCC_11775_id_GCA_003697165_2.db' >> $GENOMES_98ANI_DEREP_OUTGROUP
```

##### Extract and concatenate aa sequences.

```{bash, eval=FALSE}
num_genomes=202
min_num_genomes=$( expr $num_genomes / 2 )

anvi-get-sequences-for-hmm-hits --e $GENOMES_98ANI_DEREP_OUTGROUP --hmm-source Bacteria_71 --min-num-bins-gene-occurs $min_num_genomes --get-aa-sequences --concatenate-genes --return-best-hit --align-with muscle -o $DIR_Phylo/REDO/Dereped_98ANI_fastani_Bac71_fasta_outgroup

```

##### Build trees 
```{bash, eval=FALSE}

# trimal removes all positions in the alignment with gaps in 50% or more of the sequences
trimal -in $DIR_Phylo/REDO/Dereped_98ANI_fastani_Bac71_fasta_outgroup -out $DIR_Phylo/REDO/Dereped_98ANI_fastani_Bac71_fasta_outgroup.clean.fa -gt 0.50 -keepheader 

clusterize -n 20 -m jgiacomini@mbl.edu -l Haem_Agg_98ANI_derep_phylogeny.log iqtree -s $DIR_Phylo/REDO/Dereped_98ANI_fastani_Bac71_fasta_outgroup.clean.fa -nt AUTO -m WAG -bb 1000 -o "E_coli_str_ATCC_11775_id_GCA_003697165_2"

```

##### Plot trees 

Send files to local
```{bash, eval=FALSE}
scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/14_PHYLOGENOMICS/REDO /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/14_PHYLOGENOMICS
```


Note: Scale bar indicates the number of substitutions per site.

```{r, eval=FALSE}
library(dendextend)
library(ape)
library(dplyr)
library(phytools)
library(phylogram)

MLtree<-ape::read.tree("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/14_PHYLOGENOMICS/REDO/Dereped_98ANI_fastani_Bac71_fasta_outgroup.clean.fa.contree")

# Opening the graphical device
pdf("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/14_PHYLOGENOMICS/Haemophilus_Aggregatibacter_202_IQTree_with_outgroup_plot.pdf",
    width = 10, height = 10)

# Creating a plot
dend_MLtree <- ape::chronos(MLtree)
plot(dend_MLtree, cex = 0.25, edge.width = 0.5)
nodelabels(text = MLtree$node.label,node=2:MLtree$Nnode+Ntip(MLtree),frame="none",adj=c(1.2,-0.5), cex = 0.2, col = "red")
add.scale.bar()

# Closing the graphical device
dev.off() 

# remove outgroup from tree file for merging into pangenome and making tanglegram
MLtree_no_outgroup <- drop.tip(MLtree, "E_coli_str_ATCC_11775_id_GCA_003697165_2")
dend_MLtree_no_outgroup <- ape::chronos(MLtree_no_outgroup)

# figure
pdf("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/14_PHYLOGENOMICS/Haemophilus_Aggregatibacter_202_IQTree_NO_outgroup_plot.pdf",
    width = 10, height = 10)
plot(dend_MLtree_no_outgroup, cex = 0.25, edge.width = 0.5)
nodelabels(text = MLtree_no_outgroup$node.label,node=2:MLtree_no_outgroup$Nnode+Ntip(MLtree_no_outgroup),frame="none",adj=c(1.2,-0.5), cex = 0.2, col = "red")
dev.off() 

write.tree(MLtree_no_outgroup, file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/14_PHYLOGENOMICS/REDO/Dereped_98ANI_fastani_Bac71_fasta_outgroup_removed.clean.fa.contree")

```

##### Add Bac71 phylogeny to pangenome

```{bash, eval=FALSE}

scp -r /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/14_PHYLOGENOMICS/REDO/Dereped_98ANI_fastani_Bac71_fasta_outgroup_removed.clean.fa.contree USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/14_PHYLOGENOMICS/REDO/Dereped_98ANI_fastani_Bac71_fasta_outgroup_removed.clean.fa.contree


# import tree to pangenome
PanProject=${projectID}-98ANI
PanDir=$DIR_Pangenome/$PanProject
PAN_DB=$PanDir/${PanProject}-RESULTS/${PanProject}-PAN.db

phyName=SCG_phylogeny
phyloTree=$DIR_Phylo/REDO/Dereped_98ANI_fastani_Bac71_fasta_outgroup_removed.clean.fa.contree
ADD_phyloTree=$DIR_Phylo/${phyName}.layer_orders.tsv

awk -v phyCount="$phyName" 'BEGIN{FS=OFS="\t"}NR==1{print "item_name","data_type","data_value"}{print phyCount,"newick",$0}' $phyloTree > $ADD_phyloTree

anvi-import-misc-data -t layer_orders -p $PAN_DB $ADD_phyloTree
```



# 7. Tanglegrams


##### Extract pangenome gene-cluster freq. tree

```{bash, eval=FALSE}
PAN_DB=/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/09_PANGENOME/09_PANGENOME/P_0622_Haemophilus_Aggregatibacter-PAN.db
PanDir=/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/09_PANGENOME/09_PANGENOME/P_0622_Haemophilus_Aggregatibacter-98ANI

anvi-export-misc-data -p $PAN_DB --target-data-table layer_orders -o $PanDir/Layers_order
cat $PanDir/Layers_order | grep 'GENE_CLUSTER_FREQ_CUSTOM'| awk -F"\t" '{print $3}' >  $PanDir/gene_cluster_frequencies_newick
```

*Run the Tanglegrams.R script to produce the figure*


# 8. Phylogeny - Pangenome SCGs

By default, this program reports amino acid sequences.

--min-geometric-homogeneity-index

This can be useful if you only want to look at gene clusters that are highly conserved in geometric configuration.


##### Extract and concatenate aa seqs.

```{bash, eval=FALSE}


PanProject=${projectID}-98ANI
PanDir=$DIR_Pangenome/$PanProject
GENOMES_98ANI_DEREP_DB=$PanDir/${PanProject}-GENOMES.db
PAN_DIR=$PanDir/${PanProject}-RESULTS
PAN_DB=$PanDir/${PanProject}-RESULTS/${PanProject}-PAN.db
N_genomes=202

anvi-get-sequences-for-gene-clusters -g $GENOMES_98ANI_DEREP_DB \
                                     -p $PAN_DB \
                                     -o $PAN_DIR/single_copy_core_genes-fasta \
                                     --max-num-genes-from-each-genome 1 \
                                     --min-num-genomes-gene-cluster-occurs $N_genomes \
                                     --min-geometric-homogeneity-index 1 \
                                     --concatenate-gene-clusters 

```

GENE CLUSTER FILTERS:

Based on --min-num-genomes-gene-cluster-occurs 202
--max-num-genomes-gene-cluster-occurs 202 
--min-num-genes-from-each-genome 0
--max-num-genes-from-each-genome 1
--min-functional-homogeneity-index -1.000
--max-functional-homogeneity-index 1.000
--min-geometric-homogeneity-index 1.000
--max-geometric-homogeneity-index 1.000
--min-combined-homogeneity-index -1.000,
--max-combined-homogeneity-index 1.000 

(some of these may be default values).

Your filters resulted in 21 gene clusters that contain a total of 4242 genes.
for downstream analyses. Just so you know.

##### Build trees

```{bash, eval=FALSE}
PanProject=${projectID}-98ANI
PanDir=$DIR_Pangenome/$PanProject
PAN_DIR=$PanDir/${PanProject}-RESULTS

iqtree -s $PAN_DIR/single_copy_core_genes-fasta -nt AUTO -m WAG -bb 1000 

```


Send data to local machine for plotting.
```{bash, eval=FALSE}
scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/09_PANGENOME/P_0622_Haemophilus_Aggregatibacter-98ANI/P_0622_Haemophilus_Aggregatibacter-98ANI-RESULTS/single_copy_core_genes-fasta.contree /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/09_PANGENOME/09_PANGENOME/P_0622_Haemophilus_Aggregatibacter-98ANI/P_0622_Haemophilus_Aggregatibacter-98ANI-RESULTS/single_copy_core_genes-fasta.contree
```

##### Plot tree 

Scale bar indicates the number of substitutions per site.
```{r, eval=FALSE}
library(dendextend)
library(ape)
library(dplyr)
library(phytools)
library(phylogram)

MLtree<-ape::read.tree("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/09_PANGENOME/09_PANGENOME/P_0622_Haemophilus_Aggregatibacter-98ANI/P_0622_Haemophilus_Aggregatibacter-98ANI-RESULTS/single_copy_core_genes-fasta.contree")

# Opening the graphical device
pdf("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/14_PHYLOGENOMICS/Haemophilus_Aggregatibacter_202_IQTree_single_copy_core_genes.pdf",
    width = 10, height = 10)

# Creating a plot
dend_MLtree <- ape::chronos(MLtree)
plot(dend_MLtree, cex = 0.25, edge.width = 0.5)
nodelabels(text = MLtree$node.label,node=2:MLtree$Nnode+Ntip(MLtree),frame="none",adj=c(1.2,-0.5), cex = 0.2, col = "red")
add.scale.bar()

# Closing the graphical device
dev.off() 

```

Tanglegram Sing-copy core gene phylogeny vs. Pangenome gen-cluster frequencies clustering
```{r, eval=FALSE}
#install.packages("dendextend")
#install.packages("phytools")
#install.packages("phylogram")

library(dendextend)
library(ape)
library(dplyr)
library(phytools)
library(phylogram)

# Step 1: load trees
PANtree<-ape::read.tree("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/09_PANGENOME/P_0622_Haemophilus_Aggregatibacter-98ANI/gene_cluster_frequencies_newick")

MLtree<-ape::read.tree("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/09_PANGENOME/09_PANGENOME/P_0622_Haemophilus_Aggregatibacter-98ANI/P_0622_Haemophilus_Aggregatibacter-98ANI-RESULTS/single_copy_core_genes-fasta.contree")

# Step 2: convert to dendrogram objects
dend_MLtree <- ape::chronos(MLtree)
dend_MLtree_2 <- as.dendrogram.phylo(dend_MLtree)
dend_PANtree <- ape::chronos(PANtree)
dend_PANtree_2 <- as.dendrogram.phylo(dend_PANtree)

# Step 3: Add colors

# Extract labels from dendrogram on the left 
labels <- dend_PANtree_2 %>% set("labels_to_char") %>% labels 
#Using a metadata table with colours create a vector of 
colourslabels <- as.data.frame(labels)
write.csv(colourslabels, "/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/DATA/tanglegram_scg_vs_pan_lables.csv", row.names = FALSE)
# edit colourslabels data frame to add colors and re-load it
metadata <- read.csv("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/DATA/tanglegram_scg_vs_pan_colors.csv", header = TRUE)

labels2 <- merge(colourslabels, metadata, by.x="labels", by.y="labels", sort=F)
cols_pan <- as.character(labels2$Colours) 
cols_phy <- metadata$Colours

# Step 4:  make dendrogram list
SCG_PAN <- dendlist(
  dend_MLtree_2 %>% 
    set("labels_col", value = cols_phy),
  dend_PANtree_2 %>% 
    set("labels_col", value = cols_pan)) 


# test untangle methods
SCG_PAN %>% dendextend::untangle(method="random") %>% entanglement() #0.07311872
SCG_PAN %>% dendextend::untangle(method="labels") %>% entanglement() #0.2584771
SCG_PAN %>% dendextend::untangle(method="ladderize") %>% entanglement() #0.08805879
SCG_PAN %>% dendextend::untangle(method="step2side") %>% entanglement() #0.02238234 ###WINNER...Usually is
SCG_PAN %>% dendextend::untangle(method="step1side") %>% entanglement() #0.03265482
SCG_PAN %>% dendextend::untangle(method="DendSer") %>% entanglement() #0.3936369

# Step 5: Call the pdf command to start the plot
rows <- as.data.frame(MLtree$tip.label)
height=nrow(rows) *0.0825

# Step 6: plot tanglegram

SCG_PAN_TANGLEGRAM <- SCG_PAN %>% dendextend::untangle(method="step2side") 

PAN_TANGLEGRAM_labels <- SCG_PAN_TANGLEGRAM[[2]] %>%  labels
PAN_TANGLEGRAM_labels <- as.data.frame(PAN_TANGLEGRAM_labels)
line_colors <- merge(PAN_TANGLEGRAM_labels, metadata, by.x="PAN_TANGLEGRAM_labels", by.y="labels", sort=F)
line_colors2 <- as.character(line_colors$Colours)

pdf(file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/Sing_copy_core_genes_vs_Pan_tanglegram_colors_2.pdf", width = 15, height = height) 

SCG_PAN_TANGLEGRAM %>% plot(common_subtrees_color_lines=FALSE,
                            highlight_distinct_edges=FALSE,
                            common_subtrees_color_branches=FALSE,
                            highlight_branches_lwd=FALSE,
                            lwd=3, 
                            lab.cex = 0.65, edge.lwd = 2, 
                            margin_inner = 20, 
                            columns_width = c(1, 0.5, 1), 
                            axes=FALSE,
                            color_lines = line_colors2,
                            main = paste("entanglement =",
                                         round(entanglement(SCG_PAN_TANGLEGRAM), 4)))

dev.off()
```


Compare entanglement for Bac71 vs SCG phylogenies against pangenome gene-cluster frequencies arrangement. Which phylogeny is better? 

RUN IN R FILE not NOTEBOOK/MARKDOWN
```{r, eval=FALSE}
pdf(file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/Bac71_vs_Pan_tanglegram_mini.pdf", width = 4, height = 4) 

bac71_SCG_PAN_TANGLEGRAM %>% plot(common_subtrees_color_lines=FALSE,
                                  highlight_distinct_edges=FALSE,
                                  common_subtrees_color_branches=FALSE,
                                  highlight_branches_lwd=FALSE,
                                  lwd=.75, 
                                  lab.cex = 0.1, edge.lwd = .75, 
                                  margin_inner = 4, 
                                  columns_width = c(1, 0.5, 1), 
                                  axes=FALSE,
                                  color_lines = bac71_line_colors2,
                                  main = paste("entanglement =",
                                               round(entanglement(bac71_SCG_PAN_TANGLEGRAM), 4)))

dev.off()

pdf(file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/SCG_vs_Pan_tanglegram_mini.pdf", width = 4, height = 4) 

SCG_PAN_TANGLEGRAM %>% plot(common_subtrees_color_lines=FALSE,
                                  highlight_distinct_edges=FALSE,
                                  common_subtrees_color_branches=FALSE,
                                  highlight_branches_lwd=FALSE,
                                  lwd=.75, 
                                  lab.cex = 0.1, edge.lwd = .75, 
                                  margin_inner = 4, 
                                  columns_width = c(1, 0.5, 1), 
                                  axes=FALSE,
                                  color_lines = line_colors2,
                                  main = paste("entanglement =",
                                               round(entanglement(SCG_PAN_TANGLEGRAM), 4)))

dev.off()

```


 
# 9. Mapping (TEST SET)  

Metagenomic read recruitment (TEST SET) to dereplicated reference genome set.

Send metagenome sample info files to remote server
```{bash, eval=FALSE}
scp -r /Users/home/Veillonella/phasei.txt USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/DATA/phasei.txt

scp -r /Users/home/Veillonella/phaseii.txt USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/DATA/phaseii.txt

scp -r /Users/home/Veillonella/sequencing-metadata-count-QC_subset10.txt USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/DATA/sequencing-metadata-count-QC_subset10.txt

```


##### (TEST SET): Set up 10 largest metagenomes per oral site.

Make a symlink for the directories on Barhal that house the zipped metagenomes.
```{bash, eval=FALSE}

# Let's check to see how many metagenomes there are to work with
# There are 268 metagenomes in this diretcory
find /storage/data/00_Public_Metagenomes/oral/HMP/phaseii/01_QC/*R1.fastq.gz -type f | wc -l
find /storage/data/00_Public_Metagenomes/oral/HMP/phaseii/01_QC/*R2.fastq.gz -type f | wc -l

# There are 418 metagenomes in this diretcory
find /storage/data/00_Public_Metagenomes/oral/HMP/phasei/01_QC/*R1.fastq.gz -type f | wc -l
find /storage/data/00_Public_Metagenomes/oral/HMP/phasei/01_QC/*R2.fastq.gz -type f | wc -l



# Now we should create symbolic links for HMP metagenomes. This may be tricky since there are two different directories and one list. So I went ahead and filtered the list in R by phaseii and phasei  to create two txt files that list the sample IDs. 

# 1) create symlinks for phasei smaples
for metagenome in `cat DATA/phasei.txt`
do
ln -s $DIR_SITE/phasei/01_QC/${metagenome}-QUALITY_PASSED_R1.fastq.gz 00_READS/${metagenome}_R1.fastq.gz
ln -s $DIR_SITE/phasei/01_QC/${metagenome}-QUALITY_PASSED_R2.fastq.gz 00_READS/${metagenome}_R2.fastq.gz
done

# 2) create symlinks for phaseii smaples
for metagenome in `cat DATA/phaseii.txt`
do
ln -s $DIR_SITE/phaseii/01_QC/${metagenome}-QUALITY_PASSED_R1.fastq.gz 00_READS/${metagenome}_R1.fastq.gz
ln -s $DIR_SITE/phaseii/01_QC/${metagenome}-QUALITY_PASSED_R2.fastq.gz 00_READS/${metagenome}_R2.fastq.gz
done
 
cd 00_READS/
ls -lt TD_HC_HMP_S059_02_R2.fastq.gz
cd ..
```

Build a  list of metagenome samples for preliminary mapping with the 10 largest metagenomes for each of the sites.

```{bash, eval=FALSE}

## this is the master list of metagenome IDs from Julian (only fo ruse on Forsyth machine)
#head DATA/samples_id-QC.txt

# This is the subset of metagenome IDs that I (JJG) made on my local machine and uploaded to the remote server above
head -n1 DATA/sequencing-metadata-count-QC_subset10.txt

# make samples ID list 
cat DATA/sequencing-metadata-count-QC_subset10.txt | awk 'NR>1{print $1}' > DATA/samples_id-QC_IDs_subset10.txt

#Also make sure samples_id-QC_IDs_subset10.txt is in main directory!!!
# note that we will need to re-copy the master list into the main directory later on when we we map the full set of metagenomes

cp DATA/samples_id-QC_IDs_subset10.txt samples_id-QC_IDs.txt
```

##### (TEST SET): Run mapping, contigDB, profiling scripts

Make sure samples_id-QC_IDs.txt is in main directory!

```{bash, eval=FALSE}
# script for MAPPING in parallel (notice there are two scripts)
clusterize -n 20 -m jgiacomini@mbl.edu -log mapping_top_10.log /workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-02-indexing_and_mapping_parallel.sh P_0622_Haemophilus_Aggregatibacter 


# script for contigsdb, gene calling and hmms annotation (within anvio); need the single contigs db for the profiling steps; note that we have yet to annotate the contig db with gene functions. since this is just the first mapping run. 
clusterize -n 20 -m jgiacomini@mbl.edu -log contigdb.log /workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-01-contigsDB_annotation_hmms.sh P_0622_Haemophilus_Aggregatibacter 


# single profiling (notice there are two scripts)
clusterize -n 15 -m jgiacomini@mbl.edu -log profiling_top_10.log /workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-03-single_profiling_parallel.sh P_0622_Haemophilus_Aggregatibacter

# merged profiles (make sure "site" is in the correct position. Should match the sampleList=samples_id-QC_IDs.txt. I think if using Julian's metagenomes then site is a suffix, but with Meren's metagenomes it is a prefix.) 
# Julian's Single profiles per site syntax
#SINGLE_PROFILES_DB=$DIR_SinglePROF/*_${site}/PROFILE.db
# Meren's Single profiles per site syntax
#SINGLE_PROFILES_DB=$DIR_SinglePROF/${site}_*/PROFILE.db

# need to amend the merge profile for single site HP. I just manually copied the HP single profile directory into the merged profile directory.  
mkdir 07_MERGED_PROFILE/P_0622_Haemophilus_Aggregatibacter_HP 
cp 06_SINGLE_PROFILE/HP_HC_HMP_S082_01/* 07_MERGED_PROFILE/P_0622_Haemophilus_Aggregatibacter_HP/

clusterize -n 15 -m jgiacomini@mbl.edu -log merge_top_10.log /workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-04-merged_profile_parallel.sh P_0622_Haemophilus_Aggregatibacter 


# decompose merged profile, estimate SCG taxonomy and summarize
clusterize -n 15 -m jgiacomini@mbl.edu -log decompose_top_10.log /workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-05-decompose_summarize_profile_parallel.sh P_0622_Haemophilus_Aggregatibacter 

```

##### (TEST SET): Get detection, mean_coverage, and mean_coverage_Q2Q3 data

Relative abundance based on mean coverage of Q2-Q3.

First, we need to concatenate mean coverage (depth) and detection (breadth) for all oral sites.

```{bash, eval=FALSE}

for data in detection mean_coverage mean_coverage_Q2Q3
do
#variables
ref_site=SA
detTableTrans=$DIR_DetectionGENOMES/${projectID}-${data}.transposed.txt
detTableTransHeadless=$DIR_DetectionGENOMES/${projectID}-${data}.transposed.noheader.txt
detTSV=$DIR_DetectionGENOMES/${projectID}-${data}.txt
# transpose detection-tables to combine values from 9 oral sites
for site in PP PB BM KG TD PT TH HP SA 
do
projSite=${projectID}_${site}
siteDetTSV=08_PROFILE_SUMMARY/${projSite}-profile/bins_across_samples/${data}.txt
siteDetTSVtrans=$DIR_DetectionGENOMES/${projSite}-${data}.transposed.txt
anvi-script-transpose-matrix -o $siteDetTSVtrans $siteDetTSV
awk 'BEGIN{FS=OFS="\t"}NR>1{print $0}' $siteDetTSVtrans >> $detTableTransHeadless
done
# add header to detection-table (I choose ${ref_site}), transpose table, remove intermediate files
awk 'BEGIN{FS=OFS="\t"}NR==1{print $0}' $DIR_DetectionGENOMES/${projectID}_${ref_site}-${data}.transposed.txt > $detTableTrans
cat $detTableTransHeadless >> $detTableTrans
rm $DIR_DetectionGENOMES/${projectID}_*-${data}*txt $DIR_DetectionGENOMES/${projectID}-${data}*noheader*
anvi-script-transpose-matrix -o $detTSV $detTableTrans
rm $detTableTrans
done

```

Export total reads mapped data

```{bash, eval=FALSE}
MAPPING_STATS=$DIR_DetectionGENOMES/${projectID}-mapping_stats.txt
# add header
echo -e "layers\tnum_SNVs_reported\ttotal_reads_mapped\tnum_SCVs_reported\tnum_INDELs_reported\ttotal_reads_kept" >> $MAPPING_STATS

for site in BM TD PP PT PB TH KG SA HP 
do
MAPPING_STATS_SITE=08_PROFILE_SUMMARY/${projectID}_${site}-profile/misc_data_layers/default.txt 
# add data
awk 'BEGIN{FS=OFS="\t"}NR>1{print $0}' $MAPPING_STATS_SITE >> $MAPPING_STATS
done
```

Download data from remote server.

```{bash, eval=FALSE}
scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/13_DETECTED_GENOMES /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/

```


##### (TEST SET): Relative abundance, mean coverage and detection plots

See 'Haemophilus_Aggregatibacter_realtive_abundance_plots.Rmd' 


# 10. Annotate contig DB with gene functions.

```{bash, eval=FALSE}
# script for annotation with pfams and cogs and kofams (within anvio) 
clusterize -n 20 -m jgiacomini@mbl.edu -log annotation.log /workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-01-contigsDB_annotation_pfams_cogs.sh P_0622_Haemophilus_Aggregatibacter
```

Gene functions ...............................: 699,452 function calls from 1
                                                source (Pfam) for 347,867 unique
                                                gene calls have been added to
                                                the contigs database.

Gene functions ...............................: 793,414 function calls from 3   
                                                sources (COG20_PATHWAY,
                                                COG20_CATEGORY, COG20_FUNCTION)
                                                for 345,636 unique gene calls
                                                have been added to the contigs
                                                database.
Gene functions ...............................: 306,665 function calls from 1
                                                source (KOfam) for 291,888
                                                unique gene calls have been
                                                added to the contigs database.
Gene functions ...............................: 56,717 function calls from 1
                                                source (KEGG_Module) for 55,749
                                                unique gene calls have been
                                                added to the contigs database.
Gene functions ...............................: 56,717 function calls from 1
                                                source (KEGG_Class) for 55,749
                                                unique gene calls have been
                                                added to the contigs database.

# 11. Mapping (FULL SET)

Metagenome read recruitment to dereplicated reference genome set (FULL SET).

Set up samples.txt file with full set of metagenomes
```{bash, eval=FALSE}

# copy full list of metagenome IDs from Veillonella metapangenome directory
cp /workspace/jmarkwelchlab/P_1221_Veillonella/samples_id-QC_IDs_686.txt /workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/DATA/samples_id-QC_IDs_686.txt


# copy samples_id-QC_IDs_686.txt to the samples_id-QC_IDs.txt file in the main directory for mapping
cp DATA/samples_id-QC_IDs_686.txt samples_id-QC_IDs.txt

# check; should be 686 rows representing 686 metagenome sample IDs
cat samples_id-QC_IDs.txt | wc -l

```

Make sure samples_id-QC_IDs.txt is in the main directory. This will be the metagenome sample IDs that will be mapped

Directories that are needed to be cleared (i.e., move all stuff that is already in there from a previous into new directory named TEST_RUN05: 
05_MAPPING
06_SINGLE_PROFILE
07_MERGED_PROFILE
08_PROFILE_SUMMARY
13_DETECTED_GENOMES

##### (FULL SET) Run scripts

```{bash, eval=FALSE}
clusterize -n 15 -m jgiacomini@mbl.edu -log mapping_686_full_set.log /workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-02-indexing_and_mapping_parallel.sh P_0622_Haemophilus_Aggregatibacter 


# single profiling (notice there are two scripts)
clusterize -n 10 -m jgiacomini@mbl.edu -log profiling_686_full_set.log /workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-03-single_profiling_parallel.sh P_0622_Haemophilus_Aggregatibacter


# need to amend the merge profile for single site HP. I just manually copied the HP single profile directory into the merged profile directory.  
mkdir 07_MERGED_PROFILE/P_0622_Haemophilus_Aggregatibacter_HP 
cp 06_SINGLE_PROFILE/HP_HC_HMP_S082_01/* 07_MERGED_PROFILE/P_0622_Haemophilus_Aggregatibacter_HP/

clusterize -n 15 -m jgiacomini@mbl.edu -log merge_686_full_set.log /workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-04-merged_profile_parallel.sh P_0622_Haemophilus_Aggregatibacter 


# decompose merged profile, estimate SCG taxonomy and summarize
clusterize -n 15 -m jgiacomini@mbl.edu -log decompose_686_full_set.log /workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-05-decompose_summarize_profile_parallel.sh P_0622_Haemophilus_Aggregatibacter 


```


##### (FULL SET) Get detection, mean_coverage and mean_coverage_Q2Q3

Concatenate mean coverage (depth) and detection (breadth) for all oral sites.

```{bash, eval=FALSE}

for data in detection mean_coverage mean_coverage_Q2Q3
do
#variables
ref_site=SA
detTableTrans=$DIR_DetectionGENOMES/${projectID}-${data}.transposed.txt
detTableTransHeadless=$DIR_DetectionGENOMES/${projectID}-${data}.transposed.noheader.txt
detTSV=$DIR_DetectionGENOMES/${projectID}-${data}.txt
# transpose detection-tables to combine values from 9 oral sites
for site in PP PB BM KG TD PT TH HP SA 
do
projSite=${projectID}_${site}
siteDetTSV=08_PROFILE_SUMMARY/${projSite}-profile/bins_across_samples/${data}.txt
siteDetTSVtrans=$DIR_DetectionGENOMES/${projSite}-${data}.transposed.txt
anvi-script-transpose-matrix -o $siteDetTSVtrans $siteDetTSV
awk 'BEGIN{FS=OFS="\t"}NR>1{print $0}' $siteDetTSVtrans >> $detTableTransHeadless
done
# add header to detection-table (I choose ${ref_site}), transpose table, remove intermediate files
awk 'BEGIN{FS=OFS="\t"}NR==1{print $0}' $DIR_DetectionGENOMES/${projectID}_${ref_site}-${data}.transposed.txt > $detTableTrans
cat $detTableTransHeadless >> $detTableTrans
rm $DIR_DetectionGENOMES/${projectID}_*-${data}*txt $DIR_DetectionGENOMES/${projectID}-${data}*noheader*
anvi-script-transpose-matrix -o $detTSV $detTableTrans
rm $detTableTrans
done

```



##### (FULL SET) Export total reads mapped data

* Need to fix this; Each site has a different header (order changed). *

BM = layers  total_reads_kept        num_SCVs_reported       num_SNVs_reported       total_reads_mapped      num_INDELs_reported

TD = layers  total_reads_kept        num_SNVs_reported       total_reads_mapped      num_SCVs_reported       num_INDELs_reported

```{r, eval=FALSE}

# First need to arrange each data frame into the same order of columns
library(dplyr)

# Define the desired column order
new_cols <- c("layers", "total_reads_mapped", "total_reads_kept", "num_SNVs_reported", "num_SCVs_reported", "num_INDELs_reported")


# Loop through the input files and reorder columns BM TD PP PT PB TH KG SA HP 
for (file_name in c("/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/08_PROFILE_SUMMARY/P_0622_Haemophilus_Aggregatibacter_TD-profile/misc_data_layers/default.txt",
                    "/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/08_PROFILE_SUMMARY/P_0622_Haemophilus_Aggregatibacter_BM-profile/misc_data_layers/default.txt",
                    "/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/08_PROFILE_SUMMARY/P_0622_Haemophilus_Aggregatibacter_PP-profile/misc_data_layers/default.txt",
                    "/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/08_PROFILE_SUMMARY/P_0622_Haemophilus_Aggregatibacter_PT-profile/misc_data_layers/default.txt",
                    "/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/08_PROFILE_SUMMARY/P_0622_Haemophilus_Aggregatibacter_PB-profile/misc_data_layers/default.txt",
                    "/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/08_PROFILE_SUMMARY/P_0622_Haemophilus_Aggregatibacter_TH-profile/misc_data_layers/default.txt",
                    "/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/08_PROFILE_SUMMARY/P_0622_Haemophilus_Aggregatibacter_KG-profile/misc_data_layers/default.txt",
                    "/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/08_PROFILE_SUMMARY/P_0622_Haemophilus_Aggregatibacter_SA-profile/misc_data_layers/default.txt",
                    "/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/08_PROFILE_SUMMARY/P_0622_Haemophilus_Aggregatibacter_HP-profile/misc_data_layers/default.txt"
                    )) {
  # Read in the data
  df <- read.delim(file_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  # Reorder the columns based on the desired order
  df <- df[, new_cols]
  
  # Write out the updated data to a new file
  write.table(df, file = paste0(file_name, "_new.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
}


# Check

df <-read.table("/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/08_PROFILE_SUMMARY/P_0622_Haemophilus_Aggregatibacter_TD-profile/misc_data_layers/default.txt_new.txt", header = TRUE)


```


Now we can combine files
```{bash, eval=FALSE}
DIR_DetectionGENOMES=/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/13_DETECTED_GENOMES
projectID=P_0622_Haemophilus_Aggregatibacter
MAPPING_STATS=$DIR_DetectionGENOMES/${projectID}-mapping_stats_new.txt

# add header
echo -e "layers\ttotal_reads_mapped\ttotal_reads_kept\tnum_SNVs_reported\tnum_SCVs_reported\tnum_INDELs_reported" >> $MAPPING_STATS

for site in BM TD PP PT PB TH KG SA HP 
do
MAPPING_STATS_SITE=/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/08_PROFILE_SUMMARY/${projectID}_${site}-profile/misc_data_layers/default.txt_new.txt 
# add data
awk 'BEGIN{FS=OFS="\t"}NR>1{print $0}' $MAPPING_STATS_SITE >> $MAPPING_STATS
done
```

Download data from remote server
```{bash, eval=FALSE}
scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/13_DETECTED_GENOMES/P_0622_Haemophilus_Aggregatibacter-mapping_stats_new.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/DATA/P_0622_Haemophilus_Aggregatibacter-mapping_stats_new.txt

```

##### (FULL SET) Detection and Relative abundance figures

*See Haemophilus_Aggregatibacter_realtive_abundance_plots.Rmd for detection, relative abundance and mean coverage plots*


# 12. Final pangenome (annotated with gene functions)

```{bash, eval=FALSE}
export BLASTDB_LMDB_MAP_SIZE=10000000000


site="SA"
WORK_DIR=$PWD
TH=10
COLLECTION=Genomes
CONTIGS_DB=$WORK_DIR/$DIR_ContigsDB/${projectID}-contigs.db
PROFILE_DB=$WORK_DIR/$DIR_MergedPROF/${projectID}_${site}/PROFILE.db
PAN_SITE=$WORK_DIR/$DIR_Pangenome/internal_annotated_${projectID}
INTERNAL=$PAN_SITE/internal_${projectID}.txt
GENOMES_STORAGE=$PAN_SITE/${projectID}-GENOMES.db
PAN_RESULTS=$PAN_SITE/${projectID}-RESULTS
PAN_DB=$PAN_RESULTS/${projectID}-PAN.db
ANIb_RESULTS=$PAN_SITE/ANIb-RESULTS
ANIt_RESULTS=$PAN_SITE/ANIt-RESULTS

mkdir $PAN_SITE


genomes_98ANI=DATA/id_genomes-98ANI.txt
awk '{print $1}' $genomesID.deduplicated_full.txt > $genomes_98ANI
echo -e 'name\tbin_id\tcollection_id\tprofile_db_path\tcontigs_db_path' > $INTERNAL
cat $genomes_98ANI | awk -v collection="$COLLECTION" -v profile_db="$PROFILE_DB" -v contigs_db="$CONTIGS_DB" 'BEGIN{FS=OFS="\t"}{print $1,$1,collection,profile_db,contigs_db}' >> $INTERNAL

# build annotated genome data base
clusterize -n 20 -m jgiacomini@mbl.edu -l genomes_storage_686.log anvi-gen-genomes-storage -i $INTERNAL --gene-caller prodigal -o $GENOMES_STORAGE

# build annoated pagenome
clusterize -n 20 -m jgiacomini@mbl.edu -l pangenome_686.log anvi-pan-genome -g $GENOMES_STORAGE --align-with muscle --use-ncbi-blast --minbit 0.5 --mcl-inflation 10 -n $projectID -o $PAN_RESULTS --num-threads 20 --enforce-hierarchical-clustering 

# create reference genome metadata $PAN_LAYERS from dereplicated set of genomes (n = 202)
PAN_LAYERS=$DIR_Pangenome/${projectID}-add_info.layers.tsv
#To extract the lines from data.txt with the genes listed in genelist.txt:
grep -w -F -f $genomes_98ANI -e item P_0622_Haemophilus_Aggregatibacter-add_info.items.txt > $PAN_LAYERS

# import annotated reference genome metadata into pangenome
anvi-import-misc-data -p $PAN_DB -t layers $PAN_LAYERS

# add ANI (ANIb method) to annotated pangenome
module load blast+/2.10.1

clusterize -n 20 -m jgiacomini@mbl.edu -l ANIb_annotated_pangenome.log anvi-compute-genome-similarity -i $INTERNAL -o $ANIb_RESULTS --pan-db $PAN_DB --program pyANI --method ANIb --num-threads 20


# add ANI (tetra method) to annotated pangenome
#clusterize -n 20 -m jgiacomini@mbl.edu -l tetra.log anvi-compute-genome-similarity -i $INTERNAL -o $ANIt_RESULTS --pan-db $PAN_DB --program pyANI --method TETRA --num-threads 20

# hmms matrix
#for hmms_source in Transfer_RNAs Ribosomal_RNA_28S Archaea_76 Ribosomal_RNA_16S Bacteria_71 Protista_83 Ribosomal_RNA_18S Ribosomal_RNA_23S Ribosomal_RNA_5S Ribosomal_RNA_12S
#do
#HMMS_MATRIX=$DIR_Phylo/${projectID}-${hmms_source}-hits.tsv
#anvi-script-gen-hmm-hits-matrix-across-genomes -i $INTERNAL --hmm-source ${hmms_source} -o $HMMS_MATRIX
#done
```

##### Add Pangenome SCG phylogeny to pangenome

```{bash, eval=FALSE}
# Note that we built the SCG tree earlier
# import tree to pangenome
phyName=SCG_phylogeny_iqtree
phyloTree=$DIR_Pangenome/P_0622_Haemophilus_Aggregatibacter-98ANI/P_0622_Haemophilus_Aggregatibacter-98ANI-RESULTS/single_copy_core_genes-fasta.contree
ADD_phyloTree=$PAN_SITE/${phyName}.layer_orders.tsv

awk -v phyCount="$phyName" 'BEGIN{FS=OFS="\t"}NR==1{print "item_name","data_type","data_value"}{print phyCount,"newick",$0}' $phyloTree > $ADD_phyloTree

anvi-import-misc-data -t layer_orders -p $PAN_DB $ADD_phyloTree
```

##### Add Bacteria 71 SCG phylogeny to pangenome

```{bash, eval=FALSE}
# Note that we built the SCG tree earlier
# import tree to pangenome
phyName_bac71=Bac71_phylogeny_iqtree
phyloTree_bac71=$DIR_Phylo/REDO/Dereped_98ANI_fastani_Bac71_fasta_outgroup_removed.clean.fa.contree
ADD_phyloTree_bac71=$PAN_SITE/${phyName_bac71}.layer_orders.tsv

awk -v phyCount="$phyName_bac71" 'BEGIN{FS=OFS="\t"}NR==1{print "item_name","data_type","data_value"}{print phyCount,"newick",$0}' $phyloTree_bac71 > $ADD_phyloTree_bac71

anvi-import-misc-data -t layer_orders -p $PAN_DB $ADD_phyloTree_bac71
```


##### Gene cluster collections: Core, ACC and Singletons

```{bash, eval=FALSE}
iDir=$DIR_Pangenome/internal_annotated_${projectID}
genomeDB=$iDir/${projectID}-GENOMES.db
PAN_DB=$iDir/${projectID}-RESULTS/${projectID}-PAN.db
numGenomes="202"

# get all gene-clusters
anvi-get-sequences-for-gene-clusters -g $genomeDB -p $PAN_DB --min-num-genomes-gene-cluster-occurs 1 -o $iDir/${projectID}-accessoryGC.temp.txt
# get genusCore gene-clusters
anvi-get-sequences-for-gene-clusters -g $genomeDB -p $PAN_DB --min-num-genomes-gene-cluster-occurs $numGenomes -o $iDir/${projectID}-genusCoreGC.temp.txt
# get singleton gene-clusters
anvi-get-sequences-for-gene-clusters -g $genomeDB -p $PAN_DB --max-num-genomes-gene-cluster-occurs 1 -o $iDir/${projectID}-singletonsGC.temp.txt

# make list with unique gene-cluster IDs for genusCore & singletons
# add bin name in second column
# append to cas3 collection file
cat $iDir/${projectID}-genusCoreGC.temp.txt | grep ">" | awk -F'|' '{print $2}' | sed -e 's/gene_cluster\://' | sort | uniq | awk 'BEGIN{FS=OFS="\t"}{print $1, "Genus_core"}' > $iDir/${projectID}-collection.cas3.txt
cat $iDir/${projectID}-singletonsGC.temp.txt | grep ">" | awk -F'|' '{print $2}' | sed -e 's/gene_cluster\://' | sort | uniq | awk 'BEGIN{FS=OFS="\t"}{print $1, "Singletons"}' >> $iDir/${projectID}-collection.cas3.txt

# remove genus core & singletons gene-cluster IDs from accessory file
# print only missing gene clusters from accessoryID
# add bin name in second column
# append to cas3 collection file
comm -23 <( cat $iDir/${projectID}-accessoryGC.temp.txt | grep ">" | awk -F'|' '{print $2}' | sed -e 's/gene_cluster\://' | sort | uniq ) <( cat $iDir/${projectID}-collection.cas3.txt | cut -f1 | sort | uniq) | awk 'BEGIN{FS=OFS="\t"}{print $1, "Accessory"}' >> $iDir/${projectID}-collection.cas3.txt

# make collection info file
echo -e "Genus_core\tUNKNOWN_SOURCE\t#870000" > $iDir/${projectID}-collection.cas3-info.txt
echo -e "Accessory\tUNKNOWN_SOURCE\t#f2f2f2" >> $iDir/${projectID}-collection.cas3-info.txt
echo -e "Singletons\tUNKNOWN_SOURCE\t#db9e04" >> $iDir/${projectID}-collection.cas3-info.txt

# import cas3 collection & info to pangenome
anvi-import-collection -p $PAN_DB -C cas3 --bins-info $iDir/${projectID}-collection.cas3-info.txt $iDir/${projectID}-collection.cas3.txt

rm $iDir/*temp.txt


cat 09_PANGENOME/internal_annotated_P_0622_Haemophilus_Aggregatibacter/P_0622_Haemophilus_Aggregatibacter-collection.cas3.txt | grep -e "Genus_core" | wc -l #727
cat 09_PANGENOME/internal_annotated_P_0622_Haemophilus_Aggregatibacter/P_0622_Haemophilus_Aggregatibacter-collection.cas3.txt | grep -e "Accessory" | wc -l #8342
cat 09_PANGENOME/internal_annotated_P_0622_Haemophilus_Aggregatibacter/P_0622_Haemophilus_Aggregatibacter-collection.cas3.txt | grep -e "Singletons" | wc -l #4173


```

# 13. Metapangenome figure

Generate oral site layers for metapangenome; Run after profile summary complete.
```{bash, eval=FALSE}
# LAYERS from metadata and merge with mapping reads and relative mapping
mapIntermediate=mapping-reads.intermediate.tsv
samplesMetadata=DATA/P_0622_Haemophilus_Aggregatibacter-686-samples_metadata.txt
LAYER_ORDERS=DATA/P_0622_Haemophilus_Aggregatibacter-686-layer_orders.txt

#sequencing-metadata-count-QC.txt needs to be copied into the DATA directory
cp /workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/sequencing-metadata-count-QC.txt DATA/sequencing-metadata-count-QC.txt 



### Extract total reads mapped
for site in PP PB BM KG TD PT HP TH SA 
do
siteLayers=$DIR_SummaryPROF/${projectID}_${site}-profile/misc_data_layers/default.txt
awk -F'\t' -v OFS="\t" 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} { print $(f["layers"]), $(f["total_reads_mapped"])}' $siteLayers | awk 'BEGIN{FS=OFS="\t"}NR>1{print $1,$2}' >> $mapIntermediate
done

### Build meta data for all samples 
echo -e "Sample_ID\tOral_site\tSubject_gender\tTotal_reads\tMapping_reads\tRel_map_reads" > $samplesMetadata
for metagenome in `cat samples_id-QC_IDs.txt`
do
META=$(cat DATA/sequencing-metadata-count-QC.txt | grep -w "$metagenome" | awk 'BEGIN{FS=OFS="\t"}{print $1,$3,$4,$10}')
MAP=$(cat $mapIntermediate | grep -w "$metagenome" | awk 'BEGIN{FS=OFS="\t"}{print $2}')
echo -e "$META\t$MAP" | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6=((100 * $5 )/ $4)}' >> $samplesMetadata
done

rm $mapIntermediate

# generate LAYER_ORDERS file
# order samples using multiple conditions (intermediate files)
for site in PP PB BM KG TD PT HP TH SA 
do
# Abundance
awk 'BEGIN{FS=OFS="\t"}NR>1{ print $1,$2,$6}' $samplesMetadata | sort -k2,2 -k3,3n | cut -f1 > Abundance.intermediate.layer_orders.tsv
# Site_abundance
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$6}' $samplesMetadata | sort -k2,2 -k3,3n | cut -f1 >> Site_abundance.intermediate.layer_orders.tsv
# Site_abundance_reverse
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$6}' $samplesMetadata | sort -k2,2 -k3,3nr | cut -f1 >> Site_abundance_reverse.intermediate.layer_orders.tsv
# Site_reverse_abundance
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$6}' $samplesMetadata | sort -k2,2r -k3,3n | cut -f1 >> Site_reverse_abundance.intermediate.layer_orders.tsv
# Site_reverse_abundance_reverse
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$6}' $samplesMetadata | sort -k2,2r -k3,3nr | cut -f1 >> Site_reverse_abundance_reverse.intermediate.layer_orders.tsv
# Site_mapped_reads
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$5}' $samplesMetadata | sort -k2,2 -k3,3n | cut -f1 >> Site_mapped_reads.intermediate.layer_orders.tsv
# Site_mapped_reads_reverse
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$5}' $samplesMetadata | sort -k2,2 -k3,3nr | cut -f1 >> Site_mapped_reads_reverse.intermediate.layer_orders.tsv
# Site_reverse_mapped_reads
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$5}' $samplesMetadata | sort -k2,2r -k3,3n | cut -f1 >> Site_reverse_mapped_reads.intermediate.layer_orders.tsv
# Site_reverse_mapped_reads_reverse
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$5}' $samplesMetadata | sort -k2,2r -k3,3nr | cut -f1 >> Site_reverse_mapped_reads_reverse.intermediate.layer_orders.tsv
# Site_total_reads
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$4}' $samplesMetadata | sort -k2,2 -k3,3n | cut -f1 >> Site_total_reads.intermediate.layer_orders.tsv
# Site_total_reads_reverse
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$4}' $samplesMetadata | sort -k2,2 -k3,3nr | cut -f1 >> Site_total_reads_reverse.intermediate.layer_orders.tsv
# Site_reverse_total_reads
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$4}' $samplesMetadata | sort -k2,2r -k3,3n | cut -f1 >> Site_reverse_total_reads.intermediate.layer_orders.tsv
# Site_reverse_total_reads_reverse
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$4}' $samplesMetadata | sort -k2,2r -k3,3nr | cut -f1 >> Site_reverse_total_reads_reverse.intermediate.layer_orders.tsv
done

# create LAYER ORDERS file
echo -e 'item_name\tdata_type\tdata_value' > $LAYER_ORDERS

# linearize order and append
for file in *intermediate.layer_orders.tsv
do
FILE_NAME=$(echo $file | awk -F'.' '{print $1}')
paste -sd"," $file | awk -v file_name="$FILE_NAME" 'BEGIN{FS=OFS="\t"}{print file_name,"basic",$0}' >> $LAYER_ORDERS
done

# transpose, linearize order and append
for file in *intermediate.layer_orders.tsv
do
FILE_NAME=$(echo $file | awk -F'.' '{print "Inverted_"$1}')
cat $file | tac | paste -sd","| awk -v file_name="$FILE_NAME" 'BEGIN{FS=OFS="\t"}{print file_name,"basic",$0}' >> $LAYER_ORDERS
done

# remove intermediate files
rm *intermediate.layer_orders.tsv

```



We first did an initial run with just 10 samples for each site, so we need to move some files from the following directories for a final run with 686 samples:
 - need to clear 10_PANGENOME_SUMMARY
 - need to clear 09_PANGENOME of HMP_site pangenomes
 
 MAKE SURE TO USE script-06.2 for the annotated pangenomes
 

##### Script to build metpangenome

```{bash, eval=FALSE}

clusterize -n 20 -log metapangenome.log /workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-06.2-metapangenome_summarize_pan_parallel.sh cas3 P_0622_Haemophilus_Aggregatibacter

```

##### Ordering of genomes

```{bash, eval=FALSE}
# export order from pangenome
iDir=$DIR_Pangenome/internal_annotated_${projectID}
genomeDB=$iDir/${projectID}-GENOMES.db
panDB=$iDir/${projectID}-RESULTS/${projectID}-PAN.db
panLayersOrder=$iDir/${projectID}-layer_orders.txt

# enter anvi-interactive and edit genome arrangement as desired
anvi-display-pan -p $panDB -g $genomeDB

# export layer orders from pangenome
anvi-export-misc-data -p $panDB -o $panLayersOrder -t layer_orders

# show list of layer orders
cut -f1 $panLayersOrder


# make layer orders tree files for those we are interested in...
for pantree in GENE_CLUSTER_FREQ_CUSTOM SCG_phylogeny_iqtree Bac71_phylogeny_iqtree ANIb_percentage_identity
do
panLayersOrder=$iDir/${projectID}-layer_orders.txt
cat $panLayersOrder | grep "$pantree" | awk 'BEGIN{FS=OFS="\t"}{print $3}' > $DIR_DetectionGENOMES/${projectID}-${pantree}.tree
done

```

##### ECG and EAG

Collect and concatenate environmental core/accessory data.

```{bash, eval=FALSE}

iDir=$DIR_Pangenome/internal_annotated_${projectID}
genomeDB=$iDir/${projectID}-GENOMES.db
panDB=$iDir/${projectID}-RESULTS/${projectID}-PAN.db
panLayersOrder=$iDir/${projectID}-layer_orders.txt

# combining metapangenome
# get environmental core/accessory data
environmentalGenes=$iDir/environmental_genes.tsv
for site in BM TD PP PT PB TH KG SA HP 
do
siteDir=$DIR_Pangenome/HMP_${site}
panDB=$siteDir/${projectID}-RESULTS/${projectID}-PAN.db
itemsTable=$siteDir/${projectID}-${site}-items.txt
anvi-export-misc-data -p $panDB -t items -o $itemsTable
cat $itemsTable | awk -F'\t' -v OFS="\t" 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} { print $(f["ECGs_and_EAGs!EAG"]),$(f["ECGs_and_EAGs!ECG"]),$(f["ECGs_and_EAGs!NA"]),$(f["EAG_ECG_ratio"])}' | awk -v SITE="${site}_" -v OFS="\t" 'NR==1{print SITE$1,SITE$2,SITE$3,SITE$4}NR>1{print $1,$2,$3,$4}' > ${itemsTable}.intermediate
done

cat $DIR_Pangenome/HMP_BM/${projectID}-BM-items.txt | cut -f1 > $environmentalGenes.intermediate
paste $environmentalGenes.intermediate $DIR_Pangenome/HMP_*/${projectID}-*-items.txt.intermediate > $environmentalGenes
rm $environmentalGenes.intermediate

# import EAG/ECG info to pangenomeDB
panDB=$iDir/${projectID}-RESULTS/${projectID}-PAN.db
anvi-import-misc-data -p $panDB -t items $environmentalGenes
```

##### Import Mean coverage Q2Q3

*May want to sort by total number of reads and then choose a subset for import*

Mean coverage for BM, TD and SUPP top 30 metagenomes by Total reads (Quality checked). Total of 90 metagenomes 
```{r, eval=FALSE}

mean_cov_sorted <- read.csv("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/13_DETECTED_GENOMES/P_0622_Haemophilus_Aggregatibacter-EDITED_mean_coverage_strain_level_Q2Q3.csv", header = TRUE)


mean_cov_sorted$Internal_ID <- as.factor(mean_cov_sorted$Internal_ID)
# filter to only select BM, TD and SUPP ortal sites
mean_cov_sorted_BM_SUPP_TD <- mean_cov_sorted %>% 
  dplyr::filter(site == "BM" | site == "TD" | site == "SUPP")

# create data frame for top 30 samples per oral site base don total reads
BM_top30 <- mean_cov_sorted_BM_SUPP_TD %>% 
  dplyr::filter(site == "BM") %>% 
  dplyr::select(Internal_ID, QC_total_reads) %>% 
  dplyr::distinct() %>% 
  slice_max(QC_total_reads, n = 30)


BM_top30 <- mean_cov_sorted_BM_SUPP_TD %>% 
  dplyr::filter(site == "BM") %>% 
  dplyr::select(Internal_ID, QC_total_reads) %>% 
  dplyr::distinct() %>% 
  slice_max(QC_total_reads, n = 30)

TD_top30 <- mean_cov_sorted_BM_SUPP_TD %>% 
  dplyr::filter(site == "TD") %>% 
  dplyr::select(Internal_ID, QC_total_reads) %>% 
  dplyr::distinct() %>% 
  slice_max(QC_total_reads, n = 30)

SUPP_top30 <- mean_cov_sorted_BM_SUPP_TD %>% 
  dplyr::filter(site == "SUPP") %>% 
  dplyr::select(Internal_ID, QC_total_reads) %>% 
  dplyr::distinct() %>% 
  slice_max(QC_total_reads, n = 30)

# Combine
top30<-rbind(BM_top30,TD_top30,SUPP_top30)
View(top30)

#
mean_cov_sorted_BM_SUPP_TD_top30 <- mean_cov_sorted_BM_SUPP_TD %>% 
  filter(Internal_ID %in% top30$Internal_ID)

# pivot wider (rownames = bins, col names = Internal_ID, values = value)
mean_cov_sorted_BM_SUPP_TD_top30_wide <- pivot_wider(mean_cov_sorted_BM_SUPP_TD_top30, id_cols = bins, names_from = Internal_ID, values_from = value)

View(mean_cov_sorted_BM_SUPP_TD_top30_wide)

#Save data as text file 
write.table(mean_cov_sorted_BM_SUPP_TD_top30_wide, "/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/13_DETECTED_GENOMES/P_0622_Haemophilus_Aggregatibacter-Top30_BM_TD_SUPP_mean_coverage_strain_level_Q2Q3.txt", row.names = FALSE, quote = FALSE, sep = "\t")


```

send sorted top 30 (BM,TD,SUPP) mean coverage data to remote server.

```{bash, eval=FALSE}

scp -r /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/13_DETECTED_GENOMES/P_0622_Haemophilus_Aggregatibacter-Top30_BM_TD_SUPP_mean_coverage_strain_level_Q2Q3.txt USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/13_DETECTED_GENOMES/P_0622_Haemophilus_Aggregatibacter-Top30_BM_TD_SUPP_mean_coverage_strain_level_Q2Q3.txt

```

```{bash, eval=FALSE}
panDB=$iDir/${projectID}-RESULTS/${projectID}-PAN.db

anvi-import-misc-data -t layers -D Mean_coverage_Q2Q3 -p $panDB $DIR_DetectionGENOMES/${projectID}-mean_coverage_Q2Q3.txt
anvi-import-misc-data -t layers -D Top30_Mean_coverage_Q2Q3 -p $panDB $DIR_DetectionGENOMES/${projectID}-Top30_BM_TD_SUPP_mean_coverage_strain_level_Q2Q3.txt



```

##### Visualize metapangenome

```{bash, eval=FALSE}
iDir=$DIR_Pangenome/internal_annotated_${projectID}
genomeDB=$iDir/${projectID}-GENOMES.db
panDB=$iDir/${projectID}-RESULTS/${projectID}-PAN.db

anvi-display-pan -p $panDB -g $genomeDB
```


```{bash, eval=FALSE}
# Not working for some strange reason; just stays stuck in loading screen on google chrome (11/19/2022). Lets move it to local machine and run

scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/09_PANGENOME/internal_annotated_P_0622_Haemophilus_Aggregatibacter/P_0622_Haemophilus_Aggregatibacter-RESULTS/P_0622_Haemophilus_Aggregatibacter-PAN.db  /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/09_PANGENOME/09_PANGENOME/P_0622_Haemophilus_Aggregatibacter-PAN.db

scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/09_PANGENOME/internal_annotated_P_0622_Haemophilus_Aggregatibacter/P_0622_Haemophilus_Aggregatibacter-GENOMES.db  /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/09_PANGENOME/09_PANGENOME/P_0622_Haemophilus_Aggregatibacter-GENOMES.db


Local_GENOMES=/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/09_PANGENOME/09_PANGENOME/P_0622_Haemophilus_Aggregatibacter-GENOMES.db
Local_PAN=/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/09_PANGENOME/09_PANGENOME/P_0622_Haemophilus_Aggregatibacter-PAN.db
#anvi-migrate -p $Local_PAN --migrate-safely

anvi-display-pan -p $Local_PAN -g $Local_GENOMES
```


# 14. Gene-level coverage analysis

The goal is to extrcat gene-level coverage for each ref genome in each site.

##### Get genome coverage by site

...for downstream gene level coverage inspection and functional enrichment analyses using anvi-display-functions.

```{bash, eval=FALSE}
N=9
(
for site in BM HP KG PB PP PT SA TD TH
do
((i=i%N)); ((i++==0)) && wait
projectID=P_0622_Haemophilus_Aggregatibacter
MERGED_PROFILE=/workspace/jmarkwelchlab/${projectID}/07_MERGED_PROFILE
CONTIGS=/workspace/jmarkwelchlab/${projectID}/04_CONTIGS_DB/${projectID}-contigs.db
anvi-split -p $MERGED_PROFILE/${projectID}_${site}/PROFILE.db -c $CONTIGS -C Genomes -o 25_SPLITS/$site &
done
)

```

##### Method 1: anvi-script-gen-distribution-of-genes-in-a-bin

```{bash, eval=FALSE}

clusterize -n 20  -m jgiacomini@mbl.edu -l gene_level_coverage.log ./SCRIPTS/gene_level.sh

#!/bin/bash

SUPP_profile=07_MERGED_PROFILE/P_0622_Haemophilus_Aggregatibacter_PP/PROFILE.db
TD_profile=07_MERGED_PROFILE/P_0622_Haemophilus_Aggregatibacter_TD/PROFILE.db
BM_profile=07_MERGED_PROFILE/P_0622_Haemophilus_Aggregatibacter_BM/PROFILE.db
contigs=04_CONTIGS_DB/P_0622_Haemophilus_Aggregatibacter-contigs.db
DIR_var=22_GENE_LEVEL
genomes_98ANI=DATA/id_genomes-98ANI.txt

#mkdir $DIR_var/SUPP $DIR_var/BM $DIR_var/TD

for bin in `cat $genomes_98ANI`
do 
anvi-script-gen-distribution-of-genes-in-a-bin -c $contigs -p $SUPP_profile -C Genomes -b $bin --fraction-of-median-coverage 0.25 
mv ${bin}-ENV-DETECTION.txt $DIR_var/SUPP/${bin}-ENV-DETECTION.txt
mv ${bin}-GENE-COVs.txt $DIR_var/SUPP/${bin}-GENE-COVs.txt
done


for bin in `cat $genomes_98ANI`
do 
anvi-script-gen-distribution-of-genes-in-a-bin -c $contigs -p $TD_profile -C Genomes -b $bin --fraction-of-median-coverage 0.25 
mv ${bin}-ENV-DETECTION.txt $DIR_var/TD/${bin}-ENV-DETECTION.txt
mv ${bin}-GENE-COVs.txt $DIR_var/TD/${bin}-GENE-COVs.txt
done

for bin in `cat $genomes_98ANI`
do 
anvi-script-gen-distribution-of-genes-in-a-bin -c $contigs -p $BM_profile -C Genomes -b $bin --fraction-of-median-coverage 0.25 
mv ${bin}-ENV-DETECTION.txt $DIR_var/BM/${bin}-ENV-DETECTION.txt
mv ${bin}-GENE-COVs.txt $DIR_var/BM/${bin}-GENE-COVs.txt
done


```

##### Visualize Method 1

```{bash, eval=FALSE}

#### H. parainfluenzae TD vs SUPP specialists ####

# H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1 is a SUPP specialist strain 
#COG0156

# send to local
scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/22_GENE_LEVEL/SUPP/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1-GENE-COVs.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/SUPP/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1-GENE-COVs.txt

# View on MBL server
WD=/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter
DIR_var=$WD/22_GENE_LEVEL
SUPP_profile=$WD/07_MERGED_PROFILE/P_0622_Haemophilus_Aggregatibacter_PP/PROFILE.db
TD_profile=$WD/07_MERGED_PROFILE/P_0622_Haemophilus_Aggregatibacter_TD/PROFILE.db
BM_profile=$WD/07_MERGED_PROFILE/P_0622_Haemophilus_Aggregatibacter_BM/PROFILE.db

bin=H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1

anvi-interactive -d $DIR_var/SUPP/${bin}-GENE-COVs.txt \
                 -A $DIR_var/SUPP/${bin}-ENV-DETECTION.txt \
                 --manual \
                 -p 07_MERGED_PROFILE/P_0622_Haemophilus_Aggregatibacter_PP/tmp.db

# find gene caller IDs to search for (Run on local)
cat /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/24_FUNCTIONAL_ENRICHMENT/Gene_caller_functional_summary.txt | grep "H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1"  | grep "GC_00001314"

# or try for a way to search functional annotations and manually inspect gene level coverage 
anvi-interactive -p $WD/07_MERGED_PROFILE/P_0622_Haemophilus_Aggregatibacter_PP/PROFILE.db -c $WD/04_CONTIGS_DB/P_0622_Haemophilus_Aggregatibacter-contigs.db -C Genomes -b H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1 --gene-mode
```


##### Method 2: anvi-export-gene-coverage-and-detection.

```{bash, eval=FALSE}

num_processes=100

(
for genome in `cat $genomes_98ANI`
do
((i=i%$num_processes)); ((i++==0)) && wait
bin=$genome
for site in PP TD BM
do
tempDIR=25_SPLITS/$site/$bin
CONTIGS=$tempDIR/CONTIGS.db
PROFILE=$tempDIR/PROFILE.db
anvi-export-gene-coverage-and-detection -c $CONTIGS -p $PROFILE -O $tempDIR/${bin}_${site}_GENE_LEVEL
done &
done
)

```



Each geneome has a profile summary for each oral site that contains a gene-detection.txt file with the proportion of nucleotides with 1x coverage for each gene. We should be able to filter and combine the gene-detection.txt files from the SP, TD and KG oral sites to produce a gene-level figure of detection. 

H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1, H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1 and A_sp_HMT_458_str_W10330_id_GCA_000466335_1


##### Method 3: (prefered) Gene-level detection for genomes in their preferred oral sites

Run on MBL Barhal-01 server...

###### Concatenate, sort and filter each file...

```{bash, eval=FALSE}

WD=/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter
DIR_SummaryPROF=$WD/08_PROFILE_SUMMARY
GENE_LEVEL_DETECTION=$WD/23_GENE_LEVEL_DETECTION
projectID=P_0622_Haemophilus_Aggregatibacter

# copy required files to a new directory named 23_GENE_LEVEL_DETECTION
for bin in `cat $GENE_LEVEL_DETECTION/Gene_level_bin_list.txt`
  do
  for site in TD PP BM
    do
      mkdir $GENE_LEVEL_DETECTION/${site}
      cp $DIR_SummaryPROF/${projectID}_${site}-profile/bin_by_bin/$bin/${bin}-gene_detection.txt $GENE_LEVEL_DETECTION/${site}/${bin}-gene_detection.txt
  done
done

# concatenate, sort and filter each file 
clusterize -n 2 -l indivMetaCombiner_DET_A.log ./SCRIPTS/indivMetaCombiner_DETECTION.py  H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1 30
clusterize -n 2 -l indivMetaCombiner_DET_B.log ./SCRIPTS/indivMetaCombiner_DETECTION.py  H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1 30
clusterize -n 2 -l indivMetaCombiner_DET_C.log ./SCRIPTS/indivMetaCombiner_DETECTION.py  A_sp_HMT_458_str_W10330_id_GCA_000466335_1 30

```

###### Make combo-30 files for all genomes in pangenome.

Run on MBL Barhal-01 server...

```{bash, eval=FALSE}
clusterize -n 20 -l LOGS/COMBO30.log /workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-COMBO30.sh
```

Move files to local machine for plotting...Run on local...
```{bash, eval=FALSE}
scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/23_GENE_LEVEL_DETECTION/COMBO_30_DATA /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/07_GENE_LEVEL/
```


###### Sort files

Sort by rank order based on the proportion of samples in which the gene is detected (threshold >= 0.90); Run on local...

```{r, eval=FALSE}

library(dplyr)


# Define a function to calculate the proportion of detected samples for a row
prop_detected <- function(row) {
  detected <- sum(row[-1] >= 0.90) # Count the number of samples where the value is >= 0.90
  prop <- detected / (length(row) - 1) # Calculate the proportion
  return(prop)
}


# Apply the function to each row of the data frame and add the results to a new column

df <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/23_GENE_LEVEL_DETECTION/A_sp_HMT_458_str_W10330_id_GCA_000466335_1-gene_detection.txt-COMBO-30", header = TRUE, sep = "\t")
df$prop_detected <- apply(df, 1, prop_detected)
df_sorted <- df %>% 
  select(gene_callers_id, prop_detected) %>% 
  arrange(desc(prop_detected)) %>% 
  select(gene_callers_id) 
write.table(df_sorted, "/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/23_GENE_LEVEL_DETECTION/A_sp_HMT_458_str_W10330_id_GCA_000466335_1-gene_detection_items_order.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)


df <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/23_GENE_LEVEL_DETECTION/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1-gene_detection.txt-COMBO-30", header = TRUE, sep = "\t")
df$prop_detected <- apply(df, 1, prop_detected)
df_sorted <- df %>% 
  select(gene_callers_id, prop_detected) %>% 
  arrange(desc(prop_detected)) %>% 
  select(gene_callers_id) 
write.table(df_sorted, "/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/23_GENE_LEVEL_DETECTION/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1-gene_detection_items_order.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)

df <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/23_GENE_LEVEL_DETECTION/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1-gene_detection.txt-COMBO-30", header = TRUE, sep = "\t")
df$prop_detected <- apply(df, 1, prop_detected)
df_sorted <- df %>% 
  select(gene_callers_id, prop_detected) %>% 
  arrange(desc(prop_detected)) %>% 
  select(gene_callers_id) 
write.table(df_sorted, "/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/23_GENE_LEVEL_DETECTION/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1-gene_detection_items_order.txt", row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)

```

###### FIGURES - gene-level analysis 

*See Haemophilus and Aggregatibacter gene-level analysis and figures.Rmd for remainder of code used to build gene-level figures*


# 15. Metabolic pathway enrichment

##### Make metadata file

Run on MBL Barhal-01 server...

```{bash, eval=FALSE}

# we can create a meta data file if needed using anvi-estimate-genome-completeness
extg202=09_PANGENOME/P_0622_Haemophilus_Aggregatibacter-98ANI/P_0622_Haemophilus_Aggregatibacter-98ANI.txt
clusterize -n 10 -m jgiacomini@mbl.edu -l anvio_geneome_completeness.log anvi-estimate-genome-completeness -e $extg202 -o P_0622_Haemophilus_Aggregatibacter_G202_completness.txt

# check that gene functions sources for KOfams
CONTIGS_DB=04_CONTIGS_DB/${projectID}-contigs.db
anvi-db-info $CONTIGS_DB

# make metabolism directory
#mkdir 23_METABOLISM

# note that we need to run both of these ; the second because we need the modules.txt file for estimating module enrichment
intg202=09_PANGENOME/internal_annotated_P_0622_Haemophilus_Aggregatibacter/internal_P_0622_Haemophilus_Aggregatibacter.txt

clusterize -n 10 -m jgiacomini@mbl.edu -l G202_Estimate_metabolism.log anvi-estimate-metabolism -i $intg202 -O 23_METABOLISM/${projectID}-MetabolicPaths --matrix-format --kegg-output-modes modules

clusterize -n 10 -m jgiacomini@mbl.edu -l G202_Estimate_metabolism.log anvi-estimate-metabolism -i $intg202 -O 23_METABOLISM/${projectID}-MetabolicPaths --kegg-output-modes modules
```

We need to create txt files that list the genome IDs for groups of interest. Column names need to be "name" and "group". I think that we can use the strain-level relative abundance data frame to choose genomes and groups based on a threshold: 

1: The reference genome shows up in at least half of the samples for the oral site.
2: The reference genome shows up in at least half of the samples for the oral site above X% relative abundance.
3: The reference genome shows up in at 2/3 of the samples for the oral site.

##### Build list of groups

```{r, eval=FALSE}
# build txt list with groups of interest
#### group 1: TD H. parainfluenzae vs SUPP H. parainfluenzae


# Load mean coverage Q2Q3 data into R
strain_level_mean_covQ2Q3 <- read.csv("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/13_DETECTED_GENOMES/P_0622_Haemophilus_Aggregatibacter-EDITED_mean_coverage_strain_level_Q2Q3.csv", header = TRUE)

strain_level_mean_covQ2Q3$site <- as.factor(strain_level_mean_covQ2Q3$site)
strain_level_mean_covQ2Q3$bins <- as.factor(strain_level_mean_covQ2Q3$bins)
strain_level_mean_covQ2Q3$layers <- as.factor(strain_level_mean_covQ2Q3$layers)

################# TD ################# 
# all genomes with at least 1X coverage TD (i.e., coverage > 1) in at least half of the samples
# how many samples?
strain_level_mean_covQ2Q3_TD <- strain_level_mean_covQ2Q3 %>% 
  filter(site == "TD") %>% 
  select(layers, bins, value) %>% 
  droplevels
n_distinct(strain_level_mean_covQ2Q3_TD$layers) # 220

# create variable for if cov > 1x; call it Cov1X
strain_level_mean_covQ2Q3_TD <- strain_level_mean_covQ2Q3_TD %>% 
  mutate(Cov1X = case_when(value >= 1 ~ 1, 
                            value < 1 ~ 0))

# now we have a col that tells us if a sample passes a coverage threshold or not. 
#select bins that occur at least X times (i.e., 197/2 = 98.5 times)
Cov1X_df <- strain_level_mean_covQ2Q3_TD %>% 
  group_by(bins) %>% 
  filter(sum(Cov1X) >= (n_distinct(strain_level_mean_covQ2Q3_TD$layers)/2)) %>%
  droplevels
TD_genomes <- as.data.frame(levels(as.factor(Cov1X_df$bins)))
TD_genomes_1x_in_half<- TD_genomes %>% rename(name = 'levels(as.factor(Cov1X_df$bins))')

# manually inspect list and remove non-H. parainfluenzae genomes H_haemolyticus_str_M26166_id_GCA_003492745_1
TD_genomes_1x_in_half <- TD_genomes_1x_in_half %>% 
  filter(name != "H_influenzae_str_841_HINF_id_GCA_001058575_1")

################# SUPP ################# 
# all SUPP genomes with at least 1X coverage SUPP (i.e., coverage > 1) in at least half of the samples
# how many samples?
strain_level_mean_covQ2Q3_SUPP <- strain_level_mean_covQ2Q3 %>% 
  filter(site == "SUPP") %>% 
  select(layers, bins, value) %>% 
  droplevels
n_distinct(strain_level_mean_covQ2Q3_SUPP$layers) # 210

# create variable for if cov > 1x; call it Cov1X
strain_level_mean_covQ2Q3_SUPP <- strain_level_mean_covQ2Q3_SUPP %>% 
  mutate(Cov1X = case_when(value >= 1 ~ 1, 
                            value < 1 ~ 0))

# now we have a col that tells us if a sample passes a coverage threshold or not. 
#select bins that occur at least X times (i.e., 197/2 = 98.5 times)
SUPP_Cov1X_df <- strain_level_mean_covQ2Q3_SUPP %>% 
  group_by(bins) %>% 
  filter(sum(Cov1X) >= (n_distinct(strain_level_mean_covQ2Q3_SUPP$layers)/2)) %>%
  droplevels
SUPP_genomes <- as.data.frame(levels(as.factor(SUPP_Cov1X_df$bins)))
SUPP_genomes_1x_in_half<- SUPP_genomes %>% rename(name = 'levels(as.factor(SUPP_Cov1X_df$bins))')

# manually inspect list and remove non-H. parainfluenzae genomes H_haemolyticus_str_M26166_id_GCA_003492745_1
SUPP_genomes_1x_in_half <- SUPP_genomes_1x_in_half %>% 
  filter(name != "H_haemolyticus_str_M26166_id_GCA_003492745_1")

################# Final lists
SUPP_genomes_1x_in_half
TD_genomes_1x_in_half

# add group ID for each
SUPP_genomes_1x_in_half <- SUPP_genomes_1x_in_half %>% 
  mutate(group = "SUPP")

TD_genomes_1x_in_half <- TD_genomes_1x_in_half %>% 
  mutate(group = "TD")

# combine 
Hpara_SUPP_TD_list <- rbind(SUPP_genomes_1x_in_half, TD_genomes_1x_in_half)

write.table(Hpara_SUPP_TD_list, "/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/DATA/Hpara_SUPP_TD_Cov1x_in_half.txt", row.names = FALSE, quote = FALSE, sep = "\t")

```


##### Group 1: TD H. parainfluenzae vs SUPP H. parainfluenzae

```{bash, eval=FALSE}

# send list to server
scp -r /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/DATA/Hpara_SUPP_TD_Cov1x_in_half.txt USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/23_METABOLISM/Hpara_SUPP_TD_Cov1x_in_half.txt

#### group 1: TD H. parainfluenzae vs SUPP H. parainfluenzae (genomes with 1x or greater coverage in at least half of the samples)
23_METABOLISM/Hpara_SUPP_TD_Cov1x_in_half.txt

anvi-compute-metabolic-enrichment -M 23_METABOLISM/${projectID}-MetabolicPaths_modules.txt -G 23_METABOLISM/Hpara_SUPP_TD_Cov1x_in_half.txt -o 23_METABOLISM/HparaSupp_vs_HparaTD-metabolic_enrichment.txt

# send results to local
scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/23_METABOLISM /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/
```


##### HMT_036 versus H. haemolyticus and H. parainlfuenzae

metabolic enrichment analyses...

```{bash, eval=FALSE}

cd $mainDIR/$projectID
anvi-compute-metabolic-enrichment -M 23_METABOLISM/${projectID}-MetabolicPaths_modules.txt -G 23_METABOLISM/HMT_036_vs_haemolyticus.txt -o 23_METABOLISM/HMT_036_vs_haemolyticus-metabolic_enrichment.txt


anvi-compute-metabolic-enrichment -M 23_METABOLISM/${projectID}-MetabolicPaths_modules.txt -G 23_METABOLISM/HMT_036_vs_H_para.txt -o 23_METABOLISM/HMT_036_vs_H_para-metabolic_enrichment.txt


cat 23_METABOLISM/HMT_036_vs_haemolyticus.txt <(awk 'NR>1' 23_METABOLISM/Hpara_SUPP_TD_Cov1x_in_half.txt) > 23_METABOLISM/HMT_036_vs_H_para.txt

anvi-compute-metabolic-enrichment -M 23_METABOLISM/${projectID}-MetabolicPaths_modules.txt -G 23_METABOLISM/HMT_036_vs_H_para.txt -o 23_METABOLISM/HMT_036_vs_H_para-metabolic_enrichment.txt

```

Send reuslts to local.
```{bash, eval=FALSE}

scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/23_METABOLISM/HMT_036_vs_haemolyticus-metabolic_enrichment.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/23_METABOLISM/HMT_036_vs_haemolyticus-metabolic_enrichment.txt

scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/23_METABOLISM/HMT_036_vs_H_para-metabolic_enrichment.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/23_METABOLISM/HMT_036_vs_H_para-metabolic_enrichment.txt

```


# 16. Functional enrichment: COG20, Pfam and KOfam 

COG20_CATEGORY, KOfam, Pfam, COG20_FUNCTION, COG20_PATHWAY, KEGG_Module, KEGG_Class

We made use of anvi’o script anvi-compute-functional-enrichment to identify functional annotations that are differentially enriched or depleted in one set of genomes compared to another. First each genome was assigned to a group (oral site). Then, the script associated each gene cluster with the most frequently annotated function and generated a frequency table of functions across genomes. Finally, the enrichment test was done using a generalized linear model with a logit linkage function to obtain the enrichment score and a p-value. This analysis was performed for COG20, Pfams and KEGG annotations independently and the results were combined based on the gene cluster id.

##### Combine pangenome summaries

Combine pangenome summaries into one file for functional enrichment analysis.

```{bash, eval=FALSE}

FUN=$DIR_SummaryPAN/Functional_Enrichment.txt

head -n1 $DIR_SummaryPAN/HMP_BM/P_0622_Haemophilus_Aggregatibacter_gene_clusters_summary.txt | sed -e "s/$/\tsite/"> $FUN

for site in BM HP KG PB PP PT SA TD TH
do
cat $DIR_SummaryPAN/HMP_${site}/P_0622_Haemophilus_Aggregatibacter_gene_clusters_summary.txt | awk "NR>1" | sed -e "s/$/\t$site/" >> $FUN
done

# directories
iDir=$DIR_Pangenome/internal_annotated_${projectID}
genomeDB=$iDir/${projectID}-GENOMES.db
panDB=$iDir/${projectID}-RESULTS/${projectID}-PAN.db

# double check list of groups
cat 23_METABOLISM/Hpara_SUPP_TD_Cov1x_in_half.txt

# loaded groups into pangenome data base 
anvi-import-misc-data -p $panDB -t layers 23_METABOLISM/Hpara_SUPP_TD_Cov1x_in_half.txt

# COG20 Functions 
anvi-compute-functional-enrichment-in-pan -g $genomeDB -p $panDB --category-variable group --annotation-source COG20_FUNCTION -o 24_FUNCTIONAL_ENRICHMENT/H_para_SUPP_vs_TD-enriched-functions-COG.txt 

# Pfams
# Sept 6, 2022: Pfams warning because a function name was repeated. Weird error to get
anvi-compute-functional-enrichment-in-pan -g $genomeDB -p $panDB --category-variable group --annotation-source Pfam -o 24_FUNCTIONAL_ENRICHMENT/H_para_SUPP_vs_TD-enriched-functions-Pfam.txt 

# KOfams
anvi-compute-functional-enrichment-in-pan -g $genomeDB -p $panDB --category-variable group --annotation-source KOfam -o 24_FUNCTIONAL_ENRICHMENT/H_para_SUPP_vs_TD-enriched-functions-KOfam.txt 

# enrichment for individual gene clusters (i.e., some functional annotations are repetative becasue some gene clusters share them)
anvi-compute-functional-enrichment-in-pan -g $genomeDB -p $panDB --category-variable group --annotation-source IDENTITY -o 24_FUNCTIONAL_ENRICHMENT/H_para_SUPP_vs_TD-enriched-functions-IDENTITY.txt  --include-gc-identity-as-function


# COG20 Pathways
anvi-compute-functional-enrichment-in-pan -g $genomeDB -p $panDB --category-variable group --annotation-source COG20_PATHWAY -o 24_FUNCTIONAL_ENRICHMENT/H_para_SUPP_vs_TD-enriched-functions-COG_PATHWAY.txt 


scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/24_FUNCTIONAL_ENRICHMENT /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/

scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/24_FUNCTIONAL_ENRICHMENT/H_para_SUPP_vs_TD-enriched-functions-COG_PATHWAY.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/24_FUNCTIONAL_ENRICHMENT/H_para_SUPP_vs_TD-enriched-functions-COG_PATHWAY.txt 



```

##### Generate external genomes file

Generate external genomes txt file for the SUPP H para and TD H para specialists.

```{bash, eval=FALSE}
extg202=09_PANGENOME/P_0622_Haemophilus_Aggregatibacter-98ANI/P_0622_Haemophilus_Aggregatibacter-98ANI.txt

extg202_split_contigs=25_SPLITS/external_genomes_g202_split_contigs.txt
echo -e "name\tcontigs_db_path" > $extg202_split_contigs

for genome in `cat $extg202 | awk 'NR>1{print $1}'`
do
path=25_SPLITS/SA/${genome}/CONTIGS.db
echo -e "${genome}\t${path}" >> $extg202_split_contigs
done


Hpara_SUPP_TD_split_contigs=25_SPLITS/Hpara_SUPP_TD_Cov1x_in_half_split_contigs.txt
echo -e "name\tcontigs_db_path" > $Hpara_SUPP_TD_split_contigs

for genome in `cat 23_METABOLISM/Hpara_SUPP_TD_Cov1x_in_half.txt | awk 'NR>1{print $1}'`
do
path=SA/${genome}/CONTIGS.db
echo -e "${genome}\t${path}" >> $Hpara_SUPP_TD_split_contigs
done

```

##### Anvi-dsiplay-functions 

Anvi-dsiplay-functions for each functional aannotation source (COG20, Pfams, KOfams, KEGG modules).

```{bash, eval=FALSE}
extg202_split_contigs=25_SPLITS/external_genomes_g202_split_contigs.txt
projectID=P_0622_Haemophilus_Aggregatibacter
DIR_Pangenome=09_PANGENOME
iDir=$DIR_Pangenome/internal_annotated_${projectID}
genomeDB=$iDir/${projectID}-GENOMES.db
```


##### COG20 FUNCTIONS: H.parainfluenzae subgroups

```{bash, eval=FALSE}

projectID=P_0622_Haemophilus_Aggregatibacter
DIR_Pangenome=09_PANGENOME
iDir=$DIR_Pangenome/internal_annotated_${projectID}
genomeDB=$iDir/${projectID}-GENOMES.db

anvi-display-functions -g $genomeDB \
         --groups-txt 23_METABOLISM/Hpara_SUPP_TD_Cov1x_in_half.txt \
         --annotation-source COG20_FUNCTION \
         --profile-db 24_FUNCTIONAL_ENRICHMENT/GENOMES_db_Hpara_SUPP_TD_Cov1x_in_half-COG20_FUNCTION-PROFILE.db \
         | tee 24_FUNCTIONAL_ENRICHMENT/GENOMES_db_Hpara_SUPP_TD_Cov1x_in_half-COG20_FUNCTION.log
```



Functional occurrence stats input file path:  : /usr/local/tmp/tmpb7wjbw0e/FUNC_OCCURENCE_STATS.txt
Functional enrichment output file path:  .....: /usr/local/tmp/tmpb7wjbw0e/FUNC_ENRICHMENT_OUTPUT.txt
Temporary log file (use `--debug` to keep):  .: /usr/local/tmp/tmp0tz8h433

Num genomes ..................................: 202
Function annotation source ...................: COG20_FUNCTION
Num unique keys ..............................: 2,174
Keys correspond to ...........................: 'functions' rather than 'accession IDs'
Only the best hits are considered ............: True

```{bash, eval=FALSE}

mv /usr/local/tmp/tmpb7wjbw0e/FUNC_OCCURENCE_STATS.txt 24_FUNCTIONAL_ENRICHMENT/COG20/GENOMES_db_FUNC_OCCURENCE_STATS.txt
mv /usr/local/tmp/tmpb7wjbw0e/FUNC_ENRICHMENT_OUTPUT.txt 24_FUNCTIONAL_ENRICHMENT/COG20/GENOMES_db_FUNC_ENRICHMENT_OUTPUT.txt
```

##### COG20 FUNCTIONS: HMT-036 vs H. haemolyticus

```{bash, eval=FALSE}

projectID=P_0622_Haemophilus_Aggregatibacter
WD_DIR=/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter
genomeDB=$WD_DIR/09_PANGENOME/09_PANGENOME/${projectID}-GENOMES.db
groups_txt=$WD_DIR/23_METABOLISM/HMT_036_vs_haemolyticus.txt
out=$WD_DIR/23_METABOLISM/HMT_036_vs_haemolyticus-COG20_FUNCTION-enrichment


anvi-script-gen-function-matrix-across-genomes -g $genomeDB \
         --groups-txt $groups_txt \
         --annotation-source COG20_FUNCTION \
         --output-file-prefix $out

```


##### COG20 FUNCTIONS: HMT_036 vs H_para

```{bash, eval=FALSE}

projectID=P_0622_Haemophilus_Aggregatibacter
DIR_Pangenome=09_PANGENOME
iDir=$DIR_Pangenome/internal_annotated_${projectID}
genomeDB=$iDir/${projectID}-GENOMES.db
  
anvi-display-functions -g $genomeDB \
         --groups-txt 23_METABOLISM/HMT_036_vs_H_para.txt \
         --annotation-source COG20_FUNCTION \
         --profile-db 24_FUNCTIONAL_ENRICHMENT/COG20/HMT_036_vs_H_para_COG20_FUNCTION-PROFILE.db 
         
         
mv /usr/local/tmp/tmpnw0wkpso/FUNC_OCCURENCE_STATS.txt 24_FUNCTIONAL_ENRICHMENT/COG20/HMT_036_vs_H_para_COG20_FUNCTIONQ_OCCURENCE_STATS.txt
mv /usr/local/tmp/tmpnw0wkpso/FUNC_ENRICHMENT_OUTPUT.txt 24_FUNCTIONAL_ENRICHMENT/COG20/HMT_036_vs_H_para_COG20_FUNCTION_ENRICHMENT_OUTPUT.txt

```
Functional occurrence stats input file path:  : /usr/local/tmp/tmpnw0wkpso/FUNC_OCCURENCE_STATS.txt
Functional enrichment output file path:  .....: /usr/local/tmp/tmpnw0wkpso/FUNC_ENRICHMENT_OUTPUT.txt
Temporary log file (use `--debug` to keep):  .: /usr/local/tmp/tmp0mairo4k


Send to local
```{bash, eval=FALSE}
scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/24_FUNCTIONAL_ENRICHMENT/COG20/HMT_036_vs_H_para_COG20_FUNCTION_ENRICHMENT_OUTPUT.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/24_FUNCTIONAL_ENRICHMENT/

scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/24_FUNCTIONAL_ENRICHMENT/COG20/HMT_036_vs_H_para_COG20_FUNCTIONQ_OCCURENCE_STATS.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/24_FUNCTIONAL_ENRICHMENT/
```

Combine files
```{r, eval=FALSE}

df_1 <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/24_FUNCTIONAL_ENRICHMENT/HMT_036_vs_H_para_COG20_FUNCTION_ENRICHMENT_OUTPUT.txt", header = TRUE, sep = "\t", fill = TRUE)

df_2 <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/24_FUNCTIONAL_ENRICHMENT/HMT_036_vs_H_para_COG20_FUNCTIONQ_OCCURENCE_STATS.txt", header = TRUE, sep = "\t", fill = TRUE)

df_3 <- merge(df_1, df_2, by = "key")

write.table(df_3, "/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/24_FUNCTIONAL_ENRICHMENT/HMT_036_vs_H_para_COG20_FUNCTION-enrichment.txt", row.names = FALSE, quote = FALSE, sep = "\t")


```


##### Pfams : H parainfluenzae subgroups

```{bash, eval=FALSE}

projectID=P_0622_Haemophilus_Aggregatibacter
DIR_Pangenome=09_PANGENOME
iDir=$DIR_Pangenome/internal_annotated_${projectID}
genomeDB=$iDir/${projectID}-GENOMES.db

anvi-display-functions -g $genomeDB \
         --groups-txt 23_METABOLISM/Hpara_SUPP_TD_Cov1x_in_half.txt \
         --annotation-source Pfam \
         --profile-db 24_FUNCTIONAL_ENRICHMENT/GENOMES_db_Hpara_SUPP_TD_Cov1x_in_half-Pfam-PROFILE.db \
         | tee 24_FUNCTIONAL_ENRICHMENT/GENOMES_db_Hpara_SUPP_TD_Cov1x_in_half-Pfam.log
         
mv /usr/local/tmp/tmp3am_izgy/FUNC_OCCURENCE_STATS.txt 24_FUNCTIONAL_ENRICHMENT/Pfams/GENOMES_db_FUNC_OCCURENCE_STATS.txt
mv /usr/local/tmp/tmp3am_izgy/FUNC_ENRICHMENT_OUTPUT.txt 24_FUNCTIONAL_ENRICHMENT/Pfams/GENOMES_db_FUNC_ENRICHMENT_OUTPUT.txt
```


Functional occurrence stats input file path:  : /usr/local/tmp/tmp3am_izgy/FUNC_OCCURENCE_STATS.txt
Functional enrichment output file path:  .....: /usr/local/tmp/tmp3am_izgy/FUNC_ENRICHMENT_OUTPUT.txt
Temporary log file (use `--debug` to keep):  .: /usr/local/tmp/tmpwnry5u6q

Num genomes ..................................: 202
Function annotation source ...................: Pfam
Num unique keys ..............................: 2,570
Keys correspond to ...........................: 'functions' rather than 'accession IDs'
Only the best hits are considered ............: True


##### KOfams : H parainfluenzae subgroups

```{bash, eval=FALSE}

projectID=P_0622_Haemophilus_Aggregatibacter
DIR_Pangenome=09_PANGENOME
iDir=$DIR_Pangenome/internal_annotated_${projectID}
genomeDB=$iDir/${projectID}-GENOMES.db

anvi-display-functions -g $genomeDB \
         --groups-txt 23_METABOLISM/Hpara_SUPP_TD_Cov1x_in_half.txt \
         --annotation-source KOfam \
         --profile-db 24_FUNCTIONAL_ENRICHMENT/GENOMES_db_Hpara_SUPP_TD_Cov1x_in_half-KOfams-PROFILE.db \
         | tee 24_FUNCTIONAL_ENRICHMENT/GENOMES_db_Hpara_SUPP_TD_Cov1x_in_half-KOfams.log
         
mv /usr/local/tmp/tmp54mvedji/FUNC_OCCURENCE_STATS.txt 24_FUNCTIONAL_ENRICHMENT/KOfams/GENOMES_db_FUNC_OCCURENCE_STATS.txt
mv /usr/local/tmp/tmp54mvedji/FUNC_ENRICHMENT_OUTPUT.txt 24_FUNCTIONAL_ENRICHMENT/KOfams/GENOMES_db_FUNC_ENRICHMENT_OUTPUT.txt
```


Functional occurrence stats input file path:  : /usr/local/tmp/tmp54mvedji/FUNC_OCCURENCE_STATS.txt
Functional enrichment output file path:  .....: /usr/local/tmp/tmp54mvedji/FUNC_ENRICHMENT_OUTPUT.txt
Temporary log file (use `--debug` to keep):  .: /usr/local/tmp/tmp1xsk5v7k

Num genomes ..................................: 202
Function annotation source ...................: KOfam
Num unique keys ..............................: 2,052
Keys correspond to ...........................: 'functions' rather than 'accession IDs'
Only the best hits are considered ............: True


##### KEGG Modules : H parainfluenzae subgroups

```{bash, eval=FALSE}

projectID=P_0622_Haemophilus_Aggregatibacter
DIR_Pangenome=09_PANGENOME
iDir=$DIR_Pangenome/internal_annotated_${projectID}
genomeDB=$iDir/${projectID}-GENOMES.db

anvi-display-functions -g $genomeDB \
         --groups-txt 23_METABOLISM/Hpara_SUPP_TD_Cov1x_in_half.txt \
         --annotation-source KEGG_Module \
         --profile-db 24_FUNCTIONAL_ENRICHMENT/GENOMES_db_Hpara_SUPP_TD_Cov1x_in_half-KEGG_Module-PROFILE.db \
         | tee 24_FUNCTIONAL_ENRICHMENT/GENOMES_db_Hpara_SUPP_TD_Cov1x_in_half-KEGG_Module.log
         
mv /usr/local/tmp/tmpuqwb865z/FUNC_OCCURENCE_STATS.txt 24_FUNCTIONAL_ENRICHMENT/KEGG_Module/GENOMES_db_FUNC_OCCURENCE_STATS.txt
mv /usr/local/tmp/tmpuqwb865z/FUNC_ENRICHMENT_OUTPUT.txt 24_FUNCTIONAL_ENRICHMENT/KEGG_Module/GENOMES_db_FUNC_ENRICHMENT_OUTPUT.txt
```

Functional occurrence stats input file path:  : /usr/local/tmp/tmpuqwb865z/FUNC_OCCURENCE_STATS.txt
Functional enrichment output file path:  .....: /usr/local/tmp/tmpuqwb865z/FUNC_ENRICHMENT_OUTPUT.txt
Temporary log file (use `--debug` to keep):  .: /usr/local/tmp/tmpnt9kz3t5

Num genomes ..................................: 202
Function annotation source ...................: KEGG_Module
Num unique keys ..............................: 125
Keys correspond to ...........................: 'functions' rather than 'accession IDs'
Only the best hits are considered ............: True

##### KEGG Modules : HMT-036 vs H. haemolyticus
```{bash, eval=FALSE}

projectID=P_0622_Haemophilus_Aggregatibacter
WD_DIR=/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter
genomeDB=$WD_DIR/09_PANGENOME/09_PANGENOME/${projectID}-GENOMES.db
groups_txt=$WD_DIR/23_METABOLISM/HMT_036_vs_haemolyticus.txt
out=$WD_DIR/23_METABOLISM/HMT_036_vs_haemolyticus-KEGG_Module-enrichment

anvi-script-gen-function-matrix-across-genomes -g $genomeDB \
         --groups-txt $groups_txt \
         --annotation-source KEGG_Module \
         --output-file-prefix $out

```


##### COG20_PATHWAYS: H parainfluenzae subgroups

```{bash, eval=FALSE}

projectID=P_0622_Haemophilus_Aggregatibacter
DIR_Pangenome=09_PANGENOME
iDir=$DIR_Pangenome/internal_annotated_${projectID}
genomeDB=$iDir/${projectID}-GENOMES.db

anvi-display-functions -g $genomeDB \
         --groups-txt 23_METABOLISM/Hpara_SUPP_TD_Cov1x_in_half.txt \
         --annotation-source COG20_PATHWAY \
         --profile-db 24_FUNCTIONAL_ENRICHMENT/GENOMES_db_Hpara_SUPP_TD_Cov1x_in_half-COG20_PATHWAYS-PROFILE.db \
         | tee 24_FUNCTIONAL_ENRICHMENT/GENOMES_db_Hpara_SUPP_TD_Cov1x_in_half-COG20_PATHWAYS.log

mkdir 24_FUNCTIONAL_ENRICHMENT/COG20_PATHWAYS        
mv /usr/local/tmp/tmpi3cj775j/FUNC_OCCURENCE_STATS.txt 24_FUNCTIONAL_ENRICHMENT/COG20_PATHWAYS/GENOMES_db_FUNC_OCCURENCE_STATS.txt
mv /usr/local/tmp/tmpi3cj775j/FUNC_ENRICHMENT_OUTPUT.txt 24_FUNCTIONAL_ENRICHMENT/COG20_PATHWAYS/GENOMES_db_FUNC_ENRICHMENT_OUTPUT.txt
```


Functional occurrence stats input file path:  : /usr/local/tmp/tmpi3cj775j/FUNC_OCCURENCE_STATS.txt
Functional enrichment output file path:  .....: /usr/local/tmp/tmpi3cj775j/FUNC_ENRICHMENT_OUTPUT.txt
Temporary log file (use `--debug` to keep):  .: /usr/local/tmp/tmp_wa0yl33


Num genomes ..................................: 202
Function annotation source ...................: COG20_PATHWAY
Num unique keys ..............................: 64
Keys correspond to ...........................: 'functions' rather than 'accession IDs'
Only the best hits are considered ............: True


##### COG20_PATHWAYS: HMT-036 vs H. haemolyticus

```{bash, eval=FALSE}

projectID=P_0622_Haemophilus_Aggregatibacter
WD_DIR=/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter
genomeDB=$WD_DIR/09_PANGENOME/09_PANGENOME/${projectID}-GENOMES.db
groups_txt=$WD_DIR/23_METABOLISM/HMT_036_vs_haemolyticus.txt
out=$WD_DIR/23_METABOLISM/HMT_036_vs_haemolyticus-COG20_PATHWAY-enrichment


  
anvi-script-gen-function-matrix-across-genomes -g $genomeDB \
         --groups-txt $groups_txt \
         --annotation-source COG20_PATHWAY \
         --output-file-prefix $out

```

##### Clean up files

```{bash, eval=FALSE}
# Move data from tmp directory into 24_FUNCTIONAL_ENRICHMENT diectory
scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/24_FUNCTIONAL_ENRICHMENT /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/24_FUNCTIONAL_ENRICHMENT/Genomes_Storage_AnviDisplayFuctions/

# easier to process data manually with excel. R is having trouble reading columns properly. not sure why.
# simply sort both files by "key" and then copy and paste one into the other. manually inspect for errors.
# ALSO combined all sources into one excel file called FUNCTIONAL_ENRICHMENT_TABLE_COMBO.xlsx

# send back to server
scp -r /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/24_FUNCTIONAL_ENRICHMENT/MBL_SERVER USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/24_FUNCTIONAL_ENRICHMENT/


# move files around; only need to keep files named *FUNCTIONAL_ENRICHMENT_TABLE.txt in each source sub-directory and the COMBO file named FUNCTIONAL_ENRICHMENT_TABLE_COMBO.xlsx

tmp_dir=/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/24_FUNCTIONAL_ENRICHMENT/MBL_SERVER/24_FUNCTIONAL_ENRICHMENT
final_dir=/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/24_FUNCTIONAL_ENRICHMENT

for source in COG20 Pfams KOfams KEGG_Module
do
mv $tmp_dir/$source/*FUNCTIONAL_ENRICHMENT_TABLE.txt $final_dir/$source/FUNCTIONAL_ENRICHMENT_TABLE.txt
done

mv $tmp_dir/FUNCTIONAL_ENRICHMENT_TABLE_COMBO.xlsx $final_dir/FUNCTIONAL_ENRICHMENT_TABLE_COMBO.xlsx

rm -r /workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/24_FUNCTIONAL_ENRICHMENT/MBL_SERVER
```


Send to local machine gene caller functional enrichment data frame
```{bash, eval=FALSE}
scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/10_PANGENOME_SUMMARY/Functional_Enrichment.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/24_FUNCTIONAL_ENRICHMENT/Gene_caller_functional_summary.txt

```

# 17. GTDB-tk - classification

##### Set up environment

```{bash, eval=FALSE}
screen -S gtdb

conda deactivate # ignore it if it complains
module purge
module load miniconda/3
source /bioware/miniconda3/bashrc 
conda activate /bioware/gtdbtk-2.3.0
gtdbtk check_install

module load clusters/barhal
module load jbpc
```


##### Create directories
```{bash, eval=FALSE}
mainDIR=/workspace/jmarkwelchlab
projectID=P_0622_Haemophilus_Aggregatibacter

mkdir 29_GTDBTK
mkdir 29_GTDBTK/Genomes
mkdir 29_GTDBTK/classify_out 
mkdir 29_GTDBTK/tmp
```

##### Run scripts
```{bash, eval=FALSE}

# make custom_taxonomy_file for gtdbtk decorate
cat P_0622_Haemophilus_Aggregatibacter-add_info.items.txt | awk -F"\t" 'NR>1 {print $10, $1}' | grep -Fwf bin_list.txt > 29_GTDBTK/custom_taxonomy_file.txt

# move genomes to 29_GTDBTK/Genomes directory
for G_ID in `cat 29_GTDBTK/custom_taxonomy_file.txt | awk '{print $1}'`; do
cp 03_GENOMES_EDITED/${G_ID}.fa 29_GTDBTK/Genomes/${G_ID}.fa
done

# run script

# Added --mash_db gtdb-tk-r214.msh
clusterize -n 10 -m jgiacomini@forsyth.org -l LOGS/GTDBTK.log /workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-gtdbtk-classify_wf.sh

# send resutls to local (RUN ON LOCAL)
scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/29_GTDBTK/classify_out/gtdbtk.bac120.summary.tsv /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/29_GTDBTK/classify_out/gtdbtk.bac120.summary.tsv

scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/29_GTDBTK/custom_taxonomy_file.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/29_GTDBTK/custom_taxonomy_file.txt

```


##### Organize results

```{r, eval=FALSE}

gtdb_results <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/29_GTDBTK/classify_out/gtdbtk.bac120.summary.tsv", header = TRUE, sep = "\t")

meta <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/29_GTDBTK/custom_taxonomy_file.txt", header = FALSE)
colnames(meta) <- c("G_ID", "pangenome_ID")

gtdb_results_merged <- merge(gtdb_results, meta, by.x = "user_genome", by.y = "G_ID")

gtdb_results_merged <- gtdb_results_merged %>% 
  mutate(gtdb_species_classification = classification)

# Remove all before and up to ";s_"
gtdb_results_merged$gtdb_species_classification <- gsub(".*;s__","",gtdb_results_merged$gtdb_species_classification)

# reorder columns
gtdb_results_merged <- gtdb_results_merged %>%
  select(user_genome,pangenome_ID,classification,gtdb_species_classification,classification_method, fastani_reference,fastani_reference_radius,fastani_taxonomy, fastani_ani,fastani_af, note, other_related_references.genome_id.species_name.radius.ANI.AF.)

# write results
write.table(gtdb_results_merged, "/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/29_GTDBTK/classify_out/gtdbtk.bac120.merged_results.txt", quote = FALSE, row.names = FALSE, sep = "\t")
```

##### GTDB-tk without ANI screen
```{bash, eval=FALSE}

clusterize -n 10 -m jgiacomini@forsyth.org -l LOGS/No_ANI_GTDBTK.log /workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-gtdbtk-classify_wf_no_ANI.sh


mainDIR=/workspace/jmarkwelchlab/$projectID
OUT=$mainDIR/29_GTDBTK/classify_out_no_ANI/classify
CUSTOM=$mainDIR/29_GTDBTK/custom_taxonomy_file.txt

sed 's/ /\t/g' $CUSTOM > $mainDIR/29_GTDBTK/custom_taxonomy_file_tab.txt

gtdbtk decorate --input_tree $OUT/gtdbtk.bac120.classify.tree.1.tree \
                --output_tree $OUT/Custom_gtdbtk.bac120.classify.tree \
                --custom_taxonomy_file $mainDIR/29_GTDBTK/custom_taxonomy_file_tab.txt


# send resutls to local
scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/29_GTDBTK/classify_out_no_ANI /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/29_GTDBTK/

```



# 18. Gene and Nucleotide coverage 

The purpose of the scripts below is to visualize gene-level coverage of specific ref genomes of interest for a subset of metagenomes and targeted genes of interest.

We need a profile database and a contigs database that was generated from the anvi-split program earlier. The program anvi-interatcive will generate gene-level distrubution plots for the genome and all of the smaples in the profile db, but will also randomly select a contig split. If we want to target specific genes we may need to figure out which contig split the gene is in and then simply specify the split ID when running anvi-interactive.

1) Use anvi-interactive to identify contig splits that contain genes of interest. Search functions and copy contig name and nucleotide positions
2) Use anvi-gen-variability-profile to get data for SNVs. Be sure to use --include-split-names and --include-contig-names
3) Use anvi-get-split-coverages to get coverages of the split for a specific ref genome and a specific oral site
4) Use anvi-script-visualize-split-coverages to generate a pdf of the split coverages for a subset of samples. 
        Use --sample-data sample_data.txt with header sample_name	sample_group sample_color (colors are #000080 and #8b0000 for blue and red)
        Use --free-y-scale TRUE
        Use --snv-data if we want to include SNV data
        Use --snv-marker-transparency 0.9
        Use --snv-marker-width 0.1

Optional: Use anvi-inspect to preview split coverages before running anvi-script-visualize-split-coverages 

anvi-inspect -p PROFILE.db \
             -c CONTIGS.db \
             --split-name XXX_split_00003

Build sample_data.txt file for anvi-script-visualize-split-coverages.

```{bash, eval=FALSE}

# column 1, 2 and 10
cat DATA/sequencing-metadata-count-QC.txt | awk -v OFS="\t" '{print $1,$3,$10}' > DATA/QC_Total_reads_metagenomes.txt

```

```{r, eval=FALSE}
library(dplyr)

df <- read.table("/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/DATA/QC_Total_reads_metagenomes.txt", header = TRUE, sep = "\t")

df_TD <- df %>% filter(Sub_site == "TD") %>% droplevels()

df_TD_sorted <- df_TD %>% arrange(desc(QC_total_reads))

df_TD_top10 <- df_TD_sorted %>%  slice(1:10)

df_SUPP <- df %>% filter(Sub_site == "SUPP") %>% droplevels()

df_SUPP_sorted <- df_SUPP %>% arrange(desc(QC_total_reads))

df_SUPP_top10 <- df_SUPP_sorted %>%  slice(1:10)

dfSUPP_TD_Top10 <- rbind(df_SUPP_top10, df_TD_top10) 

dfSUPP_TD_Top10 <- dfSUPP_TD_Top10 %>% mutate(sample_color = case_when(grepl('TD', Sub_site) ~ '#8b0000',
                                                                       grepl('SUPP', Sub_site) ~ '#000080'))

write.table(dfSUPP_TD_Top10, "/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/DATA/SUPP_TD_Top10_QC_Total_reads_metagenomes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

```

```{bash, eval=FALSE}

echo -e "sample_name\tsample_group\tsample_color" > 22_GENE_LEVEL/DATA/SUPP_TD_sample_data.txt

cat DATA/SUPP_TD_Top10_QC_Total_reads_metagenomes.txt | awk -v OFS="\t" 'NR>1{print $1,$2,$4}' >> 22_GENE_LEVEL/DATA/SUPP_TD_sample_data.txt
```

```{bash, eval=FALSE}

#### H. parainfluenzae TD vs SUPP specialists ####

# H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1 is a SUPP specialist strain 
#COG0156


#1)
anvi-interactive -p 25_SPLITS/PP/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1/PROFILE.db \
-c 25_SPLITS/PP/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1/CONTIGS.db \
-P 8081

#1b)
#split = P_0622_Haemophilus_Aggregatibacter_000000005393_split_00093

#nucleotide positions: start = 1851203	stop = 1852346

# 10,200 ~ 11,400

echo -e "P_0622_Haemophilus_Aggregatibacter_000000005393_split_00093" > 25_SPLITS/PP/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1/COG0156.txt 

#2)
anvi-gen-variability-profile -p 25_SPLITS/PP/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1/PROFILE.db \
                             -c 25_SPLITS/PP/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1/CONTIGS.db \
                             --splits-of-interest 25_SPLITS/PP/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1/COG0156.txt  \
                             --include-split-names \
                             --include-contig-names \
                             -o 25_SPLITS/PP/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1/P_0622_Haemophilus_Aggregatibacter_000000005393_split_00093_SNVs.txt

#3)
anvi-get-split-coverages -p 25_SPLITS/PP/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1/PROFILE.db \
                         -c 25_SPLITS/PP/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1/CONTIGS.db \
                         --split-name P_0622_Haemophilus_Aggregatibacter_000000005393_split_00093 \
                         -o 25_SPLITS/PP/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1/P_0622_Haemophilus_Aggregatibacter_000000005393_split_00093_coverage.txt

#4)
anvi-script-visualize-split-coverages -i 25_SPLITS/PP/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1/P_0622_Haemophilus_Aggregatibacter_000000005393_split_00093_coverage.txt \
                                      -o 25_SPLITS/PP/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1/P_0622_Haemophilus_Aggregatibacter_000000005393_split_00093_inspect.pdf \
                                      --sample-data 22_GENE_LEVEL/DATA/SUPP_TD_sample_data.txt \
                                      --free-y-scale TRUE \
                                      --snv-data 25_SPLITS/PP/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1/P_0622_Haemophilus_Aggregatibacter_000000005393_split_00093_SNVs.txt \
                                      --snv-marker-transparency 0.9 \
                                      --snv-marker-width 0.1
#5) send to local
scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/25_SPLITS/PP/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1/P_0622_Haemophilus_Aggregatibacter_000000005393_split_00093_inspect.pdf /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/

```


```{bash, eval=FALSE}

echo -e "sample_name\tsample_group\tsample_color" > 22_GENE_LEVEL/DATA/SUPP_TD_sample_data.txt

cat DATA/SUPP_TD_Top10_QC_Total_reads_metagenomes.txt | awk -v OFS="\t" 'NR>1{print $1,$2,$4}' >> 22_GENE_LEVEL/DATA/SUPP_TD_sample_data.txt
```

```{bash, eval=FALSE}

#### H. parainfluenzae TD vs SUPP specialists ####

#"H. parainfluenzae str: M1C142-1 (GCA_014931375.1)" is a TD specialist
#"Pyruvate/oxaloacetate carboxyltransferase (OadA1) (PDB:2NX9)")))

#1) ID the split that contains the gene
anvi-interactive -p 25_SPLITS/TD/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1/PROFILE.db \
-c 25_SPLITS/TD/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1/CONTIGS.db \
-P 8082

#1b)
#split = P_0622_Haemophilus_Aggregatibacter_000000005386_split_00090
# COG5016


echo -e "P_0622_Haemophilus_Aggregatibacter_000000005386_split_00090" > 25_SPLITS/TD/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1/COG5016.txt 

#2)
anvi-gen-variability-profile -p 25_SPLITS/TD/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1/PROFILE.db \
                             -c 25_SPLITS/TD/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1/CONTIGS.db \
                             --splits-of-interest 25_SPLITS/TD/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1/COG5016.txt  \
                             --include-split-names \
                             --include-contig-names \
                             -o 25_SPLITS/TD/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1/P_0622_Haemophilus_Aggregatibacter_000000005386_split_00090_SNVs.txt

#3)
anvi-get-split-coverages --p 25_SPLITS/TD/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1/PROFILE.db \
                         -c 25_SPLITS/TD/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1/CONTIGS.db \
                         --split-name P_0622_Haemophilus_Aggregatibacter_000000005386_split_00090 \
                         -o 25_SPLITS/TD/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1/P_0622_Haemophilus_Aggregatibacter_000000005386_split_00090_coverage.txt

#4)
anvi-script-visualize-split-coverages -i 25_SPLITS/TD/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1/P_0622_Haemophilus_Aggregatibacter_000000005386_split_00090_coverage.txt \
                                      -o 25_SPLITS/TD/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1/P_0622_Haemophilus_Aggregatibacter_000000005386_split_00090_inspect.pdf \
                                      --sample-data 22_GENE_LEVEL/DATA/SUPP_TD_sample_data.txt \
                                      --free-y-scale TRUE \
                                      --snv-data 25_SPLITS/TD/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1/P_0622_Haemophilus_Aggregatibacter_000000005386_split_00090_SNVs.txt \
                                      --snv-marker-transparency 0.9 \
                                      --snv-marker-width 0.1
#5) send to local
scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/25_SPLITS/TD/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1/P_0622_Haemophilus_Aggregatibacter_000000005386_split_00090_inspect.pdf /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/

```



Send data to local machine

```{bash, eval=FALSE}

#H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1 is a SUPP specialist
scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/22_GENE_LEVEL/TD/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1-GENE-COVs.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/TD/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1-GENE-COVs.txt



#H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1 is a TD specialist

scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/25_SPLITS/PP/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1_PP_GENE_LEVEL-GENE-COVERAGES.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/SUPP/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1_PP_GENE_LEVEL-GENE-COVERAGES.txt

scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/25_SPLITS/TD/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1_TD_GENE_LEVEL-GENE-COVERAGES.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/TD/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1_TD_GENE_LEVEL-GENE-COVERAGES.txt

scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/25_SPLITS/PP/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1_PP_GENE_LEVEL-GENE-DETECTION.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/SUPP/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1_PP_GENE_LEVEL-GENE-DETECTION.txt

scp -r USERNAME:/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/25_SPLITS/TD/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1_TD_GENE_LEVEL-GENE-DETECTION.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/TD/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1_TD_GENE_LEVEL-GENE-DETECTION.txt
```

Run on local!

Mean nucleotide coverage of genes of interest for ref genome in oral site.

```{r, eval=FALSE}
library(tidyr)
library(dplyr)
library(egg) 

# Load mapping stats into R
mapping_stats_df <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/DATA/QC_Total_reads_metagenomes.txt", header=TRUE)


### SUPP
# subset mapping stats df by oral site of interest
SUPP_mapping_stats_df <- mapping_stats_df %>% filter(Sub_site == "SUPP") %>% droplevels()

# load coverage data
SUPP_cov_df <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/SUPP/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1-GENE-COVs.txt", header = TRUE)

# filter out gene of interest 208827
SUPP_cov_df_208827 <- SUPP_cov_df %>% filter(key == 208827) %>% droplevels()

# convert wide to long form
SUPP_melted_cov_df_208827 <- SUPP_cov_df_208827 %>% 
  gather(key, value)

# merge 
SUPP_melted_cov_df_208827 = merge.data.frame(SUPP_melted_cov_df_208827, SUPP_mapping_stats_df, by.x = "key", by.y = "Internal_ID") 


dodge <- position_dodge(width=0.9)
SUPP_plot <- ggplot(data = SUPP_melted_cov_df_208827, aes(x=reorder(key, -QC_total_reads), y=value)) + 
  geom_bar(stat = 'identity',position = dodge, color = "blue") +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 25), breaks=seq(0,25,by=5), oob=squish, labels = scales::comma) +
  xlab("Supragingival plaque sample")+
  ylab("Mean nucleotide coverage") +
  theme_classic() + 
  theme(text = element_text(size = 12),
        axis.ticks.y = element_line(size = 1),
        axis.text.x=element_blank(),
        axis.text.y =  element_text(face = "italic"),
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"))
SUPP_plot

#
SUPP_Total_reads_plot <- ggplot(data = SUPP_melted_cov_df_208827, aes(x=reorder(key, -QC_total_reads), y=QC_total_reads)) + 
  geom_bar(stat = 'identity',position = dodge, color = "blue") +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 200000000), breaks=seq(0,200000000,by=50000000), oob=squish, labels = scales::comma) +
  xlab(NULL)+
  ylab("Total reads") +
  theme_classic() + 
  theme(text = element_text(size = 8),
        axis.ticks.y = element_line(size = 1),
        axis.text.x=element_blank(),
        axis.text.y =  element_text(face = "italic"),
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"))
SUPP_Total_reads_plot

## TD
# subset mapping stats df by oral site of interest
TD_mapping_stats_df <- mapping_stats_df %>% filter(Sub_site == "TD") %>% droplevels()

# load coverage data
TD_cov_df <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/TD/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1-GENE-COVs.txt", header = TRUE)

# filter out gene of interest 208827
TD_cov_df_208827 <- TD_cov_df %>% filter(key == 208827) %>% droplevels()

# convert wide to long form
TD_melted_cov_df_208827 <- TD_cov_df_208827 %>% 
  gather(key, value)

# merge 
TD_melted_cov_df_208827 = merge.data.frame(TD_melted_cov_df_208827, TD_mapping_stats_df, by.x = "key", by.y = "Internal_ID") 


dodge <- position_dodge(width=0.9)
TD_plot <- ggplot(data = TD_melted_cov_df_208827, aes(x=reorder(key, -QC_total_reads), y=value)) + 
  geom_bar(stat = 'identity',position = dodge, color = "red") +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 25), breaks=seq(0,25,by=5), oob=squish, labels = scales::comma) +
  xlab("Tongue dorsum sample")+
  ylab("Mean nucleotide coverage") +
  theme_classic() + 
  theme(text = element_text(size = 12),
        axis.ticks.y = element_line(size = 1),
        axis.text.x=element_blank(),
        axis.text.y =  element_text(face = "italic"),
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"))
TD_plot

#
TD_Total_reads_plot <- ggplot(data = TD_melted_cov_df_208827, aes(x=reorder(key, -QC_total_reads), y=QC_total_reads)) + 
  geom_bar(stat = 'identity',position = dodge, color = "red") +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 200000000), breaks=seq(0,200000000,by=50000000), oob=squish, labels = scales::comma) +
  xlab(NULL)+
  ylab("Total reads") +
  theme_classic() + 
  theme(text = element_text(size = 8),
        axis.ticks.y = element_line(size = 1),
        axis.text.x=element_blank(),
        axis.text.y =  element_text(face = "italic"),
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white")) 
TD_Total_reads_plot


library(ggpubr)  
SUPP_and_TD_208827 <- ggpubr::ggarrange(SUPP_Total_reads_plot,SUPP_plot, TD_Total_reads_plot, TD_plot, ncol = 1, heights = c(0.3,1,0.3,1), align = "v")



title <- expression(atop(bold("H. parainfluenzae str: NCTC_7857 (GCA_900450845.1)"), scriptstyle("7-keto-8-aminopelargonate synthetase or related enzyme (BioF) (PDB:1BS0)")))
Annotated_SUPP_and_TD_208827 <- annotate_figure(SUPP_and_TD_208827,
                top=text_grob(title),
                fig.lab.size = 18)

ggsave("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/FIGURES/H_parainfluenzae_str_NCTC_7857_id_GCA_900450845_1_SUPP_vs_TD_BioF.pdf", Annotated_SUPP_and_TD_208827, width = 8, height = 10)

```


Obtain gene caller IDs manually

```{bash, eval=FALSE}

#for the H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1 genome
cat SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/24_FUNCTIONAL_ENRICHMENT/Gene_caller_functional_summary.txt | grep -e "H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1" | grep -e "COG5016"

# gene caller ID is 20021


#for the H_parainfluenzae_str_M1C160_1_id_GCA_014931275_1 genome - OadB gene
cat 24_FUNCTIONAL_ENRICHMENT/Gene_caller_functional_summary.txt | grep -e "H_parainfluenzae_str_M1C160_1_id_GCA_014931275_1" | grep -e "OadB" | awk -F'\t' '{print $5,$28}'
```


```{r, eval=FALSE}
# Load mapping stats into R
mapping_stats_df <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/DATA/QC_Total_reads_metagenomes.txt", header=TRUE)


### SUPP
# subset mapping stats df by oral site of interest
SUPP_mapping_stats_df <- mapping_stats_df %>% filter(Sub_site == "SUPP") %>% droplevels()

# load coverage data
M1C142_SUPP_cov_df <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/TD/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1_PP_GENE_LEVEL-GENE-COVERAGES.txt", header = TRUE)

# filter out gene of interest 20021
M1C142_SSUPP_cov_df_20021 <- M1C142_SUPP_cov_df %>% filter(key == 20021) %>% droplevels()

# convert wide to long form
M1C142_SSUPP_melted_cov_df_20021 <- M1C142_SSUPP_cov_df_20021 %>% 
  gather(key, value)

# merge 
M1C142_SSUPP_melted_cov_df_20021 = merge.data.frame(M1C142_SSUPP_melted_cov_df_20021, SUPP_mapping_stats_df, by.x = "key", by.y = "Internal_ID") 


dodge <- position_dodge(width=0.9)
SUPP_plot <- ggplot(data = M1C142_SSUPP_melted_cov_df_20021, aes(x=reorder(key, -QC_total_reads), y=value)) + 
  geom_bar(stat = 'identity',position = dodge, color = "blue") +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 25), breaks=seq(0,25,by=5), oob=squish, labels = scales::comma) +
  xlab("Supragingival plaque sample")+
  ylab("Mean nucleotide coverage") +
  theme_classic() + 
  theme(text = element_text(size = 12),
        axis.ticks.y = element_line(size = 1),
        axis.text.x=element_blank(),
        axis.text.y =  element_text(face = "italic"),
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"))
SUPP_plot

#
SUPP_Total_reads_plot <- ggplot(data = M1C142_SSUPP_melted_cov_df_20021, aes(x=reorder(key, -QC_total_reads), y=QC_total_reads)) + 
  geom_bar(stat = 'identity',position = dodge, color = "blue") +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 200000000), breaks=seq(0,200000000,by=50000000), oob=squish, labels = scales::comma) +
  xlab(NULL)+
  ylab("Total reads") +
  theme_classic() + 
  theme(text = element_text(size = 8),
        axis.ticks.y = element_line(size = 1),
        axis.text.x=element_blank(),
        axis.text.y =  element_text(face = "italic"),
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"))
SUPP_Total_reads_plot

## TD
# subset mapping stats df by oral site of interest
TD_mapping_stats_df <- mapping_stats_df %>% filter(Sub_site == "TD") %>% droplevels()

# load coverage data
M1C142_TD_cov_df <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/TD/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1_TD_GENE_LEVEL-GENE-COVERAGES.txt", header = TRUE)

# filter out gene of interest 20021
M1C142_TD_cov_df_20021 <- M1C142_TD_cov_df %>% filter(key == 20021) %>% droplevels()

# convert wide to long form
M1C142_TD_melted_cov_df_20021 <- M1C142_TD_cov_df_20021 %>% 
  gather(key, value)

# merge 
M1C142_TD_melted_cov_df_20021 = merge.data.frame(M1C142_TD_melted_cov_df_20021, TD_mapping_stats_df, by.x = "key", by.y = "Internal_ID") 


dodge <- position_dodge(width=0.9)
TD_plot <- ggplot(data = M1C142_TD_melted_cov_df_20021, aes(x=reorder(key, -QC_total_reads), y=value)) + 
  geom_bar(stat = 'identity',position = dodge, color = "red") +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 25), breaks=seq(0,25,by=5), oob=squish, labels = scales::comma) +
  xlab("Tongue dorsum sample")+
  ylab("Mean nucleotide coverage") +
  theme_classic() + 
  theme(text = element_text(size = 12),
        axis.ticks.y = element_line(size = 1),
        axis.text.x=element_blank(),
        axis.text.y =  element_text(face = "italic"),
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"))
TD_plot

#
TD_Total_reads_plot <- ggplot(data = M1C142_TD_melted_cov_df_20021, aes(x=reorder(key, -QC_total_reads), y=QC_total_reads)) + 
  geom_bar(stat = 'identity',position = dodge, color = "red") +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 200000000), breaks=seq(0,200000000,by=50000000), oob=squish, labels = scales::comma) +
  xlab(NULL)+
  ylab("Total reads") +
  theme_classic() + 
  theme(text = element_text(size = 8),
        axis.ticks.y = element_line(size = 1),
        axis.text.x=element_blank(),
        axis.text.y =  element_text(face = "italic"),
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white")) 
TD_Total_reads_plot


library(ggpubr)  
SUPP_and_TD_20021 <- ggpubr::ggarrange(SUPP_Total_reads_plot,SUPP_plot, TD_Total_reads_plot, TD_plot, ncol = 1, heights = c(0.3,1,0.3,1), align = "v")



title <- expression(atop(bold("H. parainfluenzae str: M1C142-1 (GCA_014931375.1)"), scriptstyle("Pyruvate/oxaloacetate carboxyltransferase (OadA1) (PDB:2NX9)")))
Annotated_SUPP_and_TD_20021 <- annotate_figure(SUPP_and_TD_20021,
                top=text_grob(title),
                fig.lab.size = 18)

ggsave("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/FIGURES/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1_SUPP_vs_TD_OadA1.pdf", Annotated_SUPP_and_TD_20021, width = 8, height = 10)

```



```{r, eval=FALSE}
#### TD
# subset mapping stats df by oral site of interest
TD_mapping_stats_df <- mapping_stats_df %>% filter(Sub_site == "TD") %>% droplevels()

# load coverage data
M1C142_TD_cov_df <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/TD/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1_TD_GENE_LEVEL-GENE-COVERAGES.txt", header = TRUE)

# filter out gene of interest 20021
M1C142_TD_cov_df_20021 <- M1C142_TD_cov_df %>% filter(key == 20021) %>% droplevels()

# convert wide to long form
M1C142_TD_melted_cov_df_20021 <- M1C142_TD_cov_df_20021 %>% 
  gather(key, value)

# merge 
M1C142_TD_melted_cov_df_20021 = merge.data.frame(M1C142_TD_melted_cov_df_20021, TD_mapping_stats_df, by.x = "key", by.y = "Internal_ID") 


# load detection coverage data
M1C142_TD_detection_df <- read.table("//Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/TD/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1_TD_GENE_LEVEL-GENE-DETECTION.txt", header = TRUE)

# filter out gene of interest 20021
M1C142_TD_detection_df_20021 <- M1C142_TD_detection_df %>% filter(key == 20021) %>% droplevels()

# convert wide to long form
M1C142_TD_melted_detection_df_20021 <- M1C142_TD_detection_df_20021 %>% 
  gather(key, value)

# rename detection value 
M1C142_TD_melted_detection_df_20021 <- M1C142_TD_melted_detection_df_20021 %>% 
  rename(breadth = value)

# merge with coverage data
M1C142_TD_melted_breadth_and_cov_df_20021 = merge.data.frame(M1C142_TD_melted_cov_df_20021, M1C142_TD_melted_detection_df_20021, by = "key") 

M1C142_TD_melted_breadth_and_cov_df_20021 <- M1C142_TD_melted_breadth_and_cov_df_20021 %>% 
  rename(depth = value,
         breadth = breadth.y, 
         site = Sub_site) %>% 
  select(depth, breadth, key, site, QC_total_reads)



M1C142_TD_melted_breadth_and_cov_df_20021_plot <- ggplot(M1C142_TD_melted_breadth_and_cov_df_20021, aes(x = depth, y = breadth)) +
  geom_point() + 
  theme_classic() +
  xlab("Mean nucleotide depth of coverage")+
  ylab("Gene breadth of coverage")


title <- expression(atop(bold("H. parainfluenzae str: M1C142-1 (GCA_014931375.1)"), scriptstyle("Pyruvate/oxaloacetate carboxyltransferase (OadA1) (PDB:2NX9)")))


Annotated_M1C142_TD_melted_breadth_and_cov_df_20021_plot <- annotate_figure(M1C142_TD_melted_breadth_and_cov_df_20021_plot,
                                               top=text_grob(title),
                                               fig.lab.size = 20)


ggsave("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/FIGURES/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1_SUPP_vs_TD_OadA1_depth_vs_breadth.pdf", 
       Annotated_M1C142_TD_melted_breadth_and_cov_df_20021_plot, 
       width = 5, height = 4)



#################
##### SUPP ######
# subset mapping stats df by oral site of interest
SUPP_mapping_stats_df <- mapping_stats_df %>% filter(Sub_site == "SUPP") %>% droplevels()

# load coverage data
M1C142_SUPP_cov_df <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/SUPP/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1_PP_GENE_LEVEL-GENE-COVERAGES.txt", header = TRUE)

# filter out gene of interest 20021
M1C142_SUPP_cov_df_20021 <- M1C142_SUPP_cov_df %>% filter(key == 20021) %>% droplevels()

# convert wide to long form
M1C142_SUPP_melted_cov_df_20021 <- M1C142_SUPP_cov_df_20021 %>% 
  gather(key, value)

# merge 
M1C142_SUPP_melted_cov_df_20021 = merge.data.frame(M1C142_SUPP_melted_cov_df_20021, SUPP_mapping_stats_df, by.x = "key", by.y = "Internal_ID") 


# load detectionerage data
M1C142_SUPP_detection_df <- read.table("//Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/SUPP/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1_PP_GENE_LEVEL-GENE-DETECTION.txt", header = TRUE)

# filter out gene of interest 20021
M1C142_SUPP_detection_df_20021 <- M1C142_SUPP_detection_df %>% filter(key == 20021) %>% droplevels()

# convert wide to long form
M1C142_SUPP_melted_detection_df_20021 <- M1C142_SUPP_detection_df_20021 %>% 
  gather(key, value)

# rename detection value 
M1C142_SUPP_melted_detection_df_20021 <- M1C142_SUPP_melted_detection_df_20021 %>% 
  rename(breadth = value)

# merge with coverage data
M1C142_SUPP_melted_breadth_and_cov_df_20021 = merge.data.frame(M1C142_SUPP_melted_cov_df_20021, M1C142_SUPP_melted_detection_df_20021, by = "key") 

M1C142_SUPP_melted_breadth_and_cov_df_20021 <- M1C142_SUPP_melted_breadth_and_cov_df_20021 %>% 
  rename(depth = value,
         site = Sub_site) %>% 
  select(depth, breadth, key, site, QC_total_reads)



M1C142_SUPP_melted_breadth_and_cov_df_20021_plot <- ggplot(M1C142_SUPP_melted_breadth_and_cov_df_20021, aes(x = depth, y = breadth)) +
  geom_point() + 
  theme_classic() +
  xlab("Mean nucleotide depth of coverage")+
  ylab("Gene breadth of coverage")


title <- expression(atop(bold("H. parainfluenzae str: M1C142-1 (GCA_014931375.1)"), scriptstyle("Pyruvate/oxaloacetate carboxyltransferase (OadA1) (PDB:2NX9)")))


Annotated_M1C142_SUPP_melted_breadth_and_cov_df_20021_plot <- annotate_figure(M1C142_SUPP_melted_breadth_and_cov_df_20021_plot,
                                                                            top=text_grob(title),
                                                                            fig.lab.size = 20)



library(ggpubr)  
M1C142_SUPP_and_TD_20021_breath_vs_depth <- ggpubr::ggarrange(Annotated_M1C142_TD_melted_breadth_and_cov_df_20021_plot,
                                        Annotated_M1C142_SUPP_melted_breadth_and_cov_df_20021_plot,
                                        ncol = 1, heights = c(1,1), align = "v", labels = c("TD", "SUPP"))



#####combine the data and make a single plot with colored points for oral site

M1C142_SUPP_melted_breadth_and_cov_df_20021
M1C142_TD_melted_breadth_and_cov_df_20021

M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021 <- rbind(M1C142_SUPP_melted_breadth_and_cov_df_20021, M1C142_TD_melted_breadth_and_cov_df_20021)


M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021_plot <- ggplot(M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021, aes(x = depth, y = breadth, color = site, shape = site)) +
  geom_point(size = 3, alpha = 0.5) + 
  theme_classic() +
  xlab("Mean nucleotide depth of coverage")+
  ylab("Gene breadth of coverage") +
  scale_color_manual(values=c("#56B4E9","#E69F00")) +
  scale_shape_manual(values = c(17, 19)) #+
  #geom_hline(yintercept=0.9, linetype="dashed", 
             #color = "red", size=0.5, alpha = 0.25) +
  #geom_vline(xintercept=1, linetype="dashed", 
             #color = "blue", size=0.5, alpha = 0.25)

M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021_plot



title <- expression(atop(bold("H. parainfluenzae str: M1C142-1 (GCA_014931375.1)"), scriptstyle("Pyruvate/oxaloacetate carboxyltransferase (OadA1) (PDB:2NX9)")))

Annotated_M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021_plot <- annotate_figure(M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021_plot,
                                                                              top=text_grob(title),
                                                                              fig.lab.size = 20)

ggsave("/Users/home/SPECIES_LEVEL_PANGENOMES/Haemophilus_and_Aggregatibacter/22_GENE_LEVEL/FIGURES/H_parainfluenzae_str_M1C142_1_id_GCA_014931375_1_SUPP_vs_TD_OadA1_depth_vs_breadth.pdf", 
       Annotated_M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021_plot, 
       width = 6, height = 5)

```


```{r, eval=FALSE}
# How many samples with breadth above 0.9?
TD_samples_above_90percent <- M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021 %>% filter(site == "TD", breadth > 0.9) %>% summarise(n = n())
total_TD_smaples <- M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021 %>% filter(site == "TD") %>% summarise(n = n())
Prop_above_90_TD <- (TD_samples_above_90percent/total_TD_smaples)*100

TD_samples_above_25percent <- M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021 %>% filter(site == "TD", breadth > 0.25) %>% summarise(n = n())
total_TD_smaples <- M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021 %>% filter(site == "TD") %>% summarise(n = n())
Prop_above_25_TD <- (TD_samples_above_25percent/total_TD_smaples)*100

TD_samples_above_50percent <- M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021 %>% filter(site == "TD", breadth > 0.5) %>% summarise(n = n())
total_TD_smaples <- M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021 %>% filter(site == "TD") %>% summarise(n = n())
Prop_above_50_TD <- (TD_samples_above_50percent/total_TD_smaples)*100

TD_samples_above_1X <- M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021 %>% filter(site == "TD", depth > 1) %>% summarise(n = n())
total_TD_smaples <- M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021 %>% filter(site == "TD") %>% summarise(n = n())
Prop_TD_samples_above_1X <- (TD_samples_above_1X/total_TD_smaples)*100


SUPP_samples_above_90percent <- M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021 %>% filter(site == "SUPP", breadth > 0.9) %>% summarise(n = n())
total_SUPP_smaples <- M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021 %>% filter(site == "SUPP") %>% summarise(n = n())
Prop_above_90_SUPP <- (SUPP_samples_above_90percent/total_SUPP_smaples)*100


SUPP_samples_above_50percent <- M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021 %>% filter(site == "SUPP", breadth > 0.5) %>% summarise(n = n())
total_SUPP_smaples <- M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021 %>% filter(site == "SUPP") %>% summarise(n = n())
Prop_above_50_SUPP <- (SUPP_samples_above_50percent/total_SUPP_smaples)*100


SUPP_samples_above_10percent <- M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021 %>% filter(site == "SUPP", breadth > 0.1) %>% summarise(n = n())
total_SUPP_smaples <- M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021 %>% filter(site == "SUPP") %>% summarise(n = n())
Prop_above_10_SUPP <- (SUPP_samples_above_10percent/total_SUPP_smaples)*100


SUPP_samples_above_1X <- M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021 %>% filter(site == "SUPP", depth > 1) %>% summarise(n = n())
total_SUPP_smaples <- M1C142_SUPP_and_TD_melted_breadth_and_cov_df_20021 %>% filter(site == "SUPP") %>% summarise(n = n())
Prop_SUPP_samples_above_1X <- (SUPP_samples_above_1X/total_SUPP_smaples)*100


```


