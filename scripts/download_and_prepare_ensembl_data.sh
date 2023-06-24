#!/usr/bin/bash
## Download the content from ensembl

MUS_CHR='Mus_musculus.GRCm38.dna.chromosome.'
HS_CHR='Homo_sapiens.GRCh38.dna.chromosome.'
EXTENSION='.fa.gz'
DIRECTORIES=(cds ncrna cdna pep dna)
RELEASE='release-90'
BASE_URL='ftp://anonymous@ftp.ensembl.org/pub/'$RELEASE'/fasta/'
DIR_PATH=${BASE_URL:16}
MOUSE_URL=$BASE_URL'mus_musculus/'
HS_URL=$BASE_URL'homo_sapiens/'
OTHERS=(MT X Y)
DOWNLOAD_DIR='/ontology/databaseImport/ensembl/data/testing/'
CURDIR=$(pwd)

echo $CURDIR
echo "Fetching mouse sequence data from ensembl ..."
rm -rf $DOWNLOAD_DIR'data'
mkdir $DOWNLOAD_DIR'data'
for dir in "${DIRECTORIES[@]}";
do
    echo "Current directory: "$dir
    if [[ "dna" == "$dir" ]]
    then
        DNA_DIR_PATH=$DOWNLOAD_DIR"/data/"$DIR_PATH"mus_musculus/dna/"
        mkdir -p $DNA_DIR_PATH
        echo $DNA_DIR_PATH
        for i in {1..19};
        do
            wget -P $DNA_DIR_PATH $MOUSE_URL"dna/"$MUS_CHR$i$EXTENSION
        done
        for i in "${OTHERS[@]}";
        do
            wget -P $DNA_DIR_PATH $MOUSE_URL"dna/"$MUS_CHR$i$EXTENSION
        done
    else
        wget -P $DOWNLOAD_DIR'/data/' -r $MOUSE_URL$dir'/'
    fi
done

echo "Fetching Human chromosomal data..."
for dir in "${DIRECTORIES[@]}";
do
    echo "Current Directory"$dir
    if [[ "dna" == "$dir" ]]
    then
        DNA_DIR_PATH=$DOWNLOAD_DIR"/data/"$DIR_PATH'homo_sapiens/dna/'
        mkdir -p $DNA_DIR_PATH
        for i in {1..22};
        do
           wget -P $DNA_DIR_PATH $HS_URL"dna/"$HS_CHR$i$EXTENSION
        done
        for i in "${OTHERS[@]}";
        do
            wget -P $DNA_DIR_PATH $HS_URL"dna/"$HS_CHR$i$EXTENSION
        done
    else
        # For all the other directories
        wget -P $DOWNLOAD_DIR'/data/' -r $HS_URL$dir'/'
    fi
done

## Moving the files to the right location
mkdir -p $DOWNLOAD_DIR'/data/ensembl/hs'
mkdir -p $DOWNLOAD_DIR'/data/ensembl/mus'
mv $DOWNLOAD_DIR'/data/'ftp.ensembl.org/pub/$RELEASE/fasta/homo_sapiens/* $DOWNLOAD_DIR'/data/ensembl/hs/'
mv $DOWNLOAD_DIR'/data/'ftp.ensembl.org/pub/$RELEASE/fasta/mus_musculus/* $DOWNLOAD_DIR'/data/ensembl/mus/'

## Possible enhancement. Delete only if the directories are empty.
rm -rf $DOWNLOAD_DIR'/data/ftp.ensembl.org'

echo "Unzipping the downloaded files"
for i in "${DIRECTORIES[@]}";
do
  gunzip $DOWNLOAD_DIR'/data/ensembl/mus/'$i'/'*.gz
  gunzip $DOWNLOAD_DIR'/data/ensembl/hs/'$i'/'*.gz
done