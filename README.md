# ncbiClusters

Get a report of new clusters via the NCBI Pathogen Detection Pipeline

## Usage
    
    perl downloadClusters.pl resultsSet remoteDir
    
Where `remoteDir` is the Pathogen Detection Pipeline directory from which to retrieve results and `resultsSet`
is the taxon to download, e.g., Listeria


## Requirements

* Perl v5.12 or greater
* BioPerl
* ps2pdf (required for pdf instead of ps reports)

## Getting help
    
    perl scripts/downloadClusters.pl -h

