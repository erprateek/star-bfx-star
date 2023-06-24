#!/bin/bash

url_file=${1}
download_location=${2}

echo "Making directory"
mkdir -p $download_location

cd $download_location

echo "getting"
wget --reject-regex ".*\.pdf\..*" -i ../$url_file
