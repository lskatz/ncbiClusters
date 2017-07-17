#!/bin/bash

set -e

TAXON=$1
OUT=$2

if [ "$OUT" == "" ]; then
  echo "Usage: $(basename $0) taxon outdir/"
  echo 
  echo "  Where taxon is something listed on the NCBI pathogens page"
  exit 0;
fi

cd $OUT
wget --continue -r \
  -X/pathogen/Results/$TAXON/latest_snps/SNP_trees \
  -X/pathogen/Results/$TAXON/latest_snps/Trees \
  ftp://ftp.ncbi.nlm.nih.gov/pathogen/Results/$TAXON/latest_kmer/

echo "Results can be found in $OUT"
