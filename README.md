# ncbiClusters

Get a report of new clusters via the NCBI Pathogen Detection Pipeline

## Usage
    
    perl scripts/downloadClusters.pl --resultsSet Listeria latest
    R --no-save < scripts/snpLmNCBI.R

## Requirements

* Perl v5.12 or greater
* BioPerl
* ps2pdf (required for pdf instead of ps reports)

## Getting help
    
    perl scripts/downloadClusters.pl -h

