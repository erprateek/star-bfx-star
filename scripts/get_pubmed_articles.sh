#!/bin/bash

echo "Getting files"
#URL='ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/'
URL='ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/'
BASE_NAME='articles.txt.'
EXTENSION='.tar.gz'
#NAME_LIST=("O-Z" "I-N" "C-H" "0-9A-B")
NAME_LIST=("comm_use.0-9A-B.txt.tar.gz" "comm_use.A-B.xml.tar.gz" "comm_use.C-H.txt.tar.gz" "comm_use.C-H.xml.tar.gz" "comm_use.I-N.txt.tar.gz" "comm_use.I-N.xml.tar.gz" "comm_use.O-Z.txt.tar.gz" "comm_use.O-Z.xml.tar.gz" "non_comm_use.0-9A-B.txt.tar.gz" "non_comm_use.A-B.xml.tar.gz" "non_comm_use.C-H.txt.tar.gz" "non_comm_use.C-H.xml.tar.gz" "non_comm_use.I-N.txt.tar.gz" "non_comm_use.I-N.xml.tar.gz" "non_comm_use.O-Z.txt.tar.gz" "non_comm_use.O-Z.xml.tar.gz" )
for name in "${NAME_LIST[@]}"
do
    FILE_NAME=$URL$name
    #LOCAL_FILE=$BASE_NAME$name$EXTENSION
    # If not exists file or unzipped file
    echo "Getting ... "$FILE_NAME
    wget $FILE_NAME 
done
echo "Done getting files"

# Extract tar file contents
for name in "${NAME_LIST[@]}"
do
    LOCAL_FILE=$name
    gunzip $LOCAL_FILE
done
