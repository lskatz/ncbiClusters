#!/usr/bin/env perl 
#
# Author: Error Strain
# Modified: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Getopt::Long;
use Net::FTP;
use File::Basename qw/basename/;
use File::Temp qw/tempdir/;

# Import local modules
#use FindBin;
#use lib "$FindBin::RealBin/../lib";
#use Docopt;

local $0=basename($0);
sub logmsg { print STDERR "$0: @_\n"; }

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s set|results|resultsSet=s)) or die $!;
  $$settings{domain}||="ftp.ncbi.nlm.nih.gov";
  $$settings{tempdir}||=tempdir("$0.XXXXXX",TMPDIR=>1,CLEANUP=>1);
  $$settings{set}||="Listeria";

  my $remoteDir=$ARGV[0];

  if(!$remoteDir || $$settings{help}){
    die usage();
  }

  logmsg "Temporary directory is $$settings{tempdir}";
  logmsg "Downloading for $$settings{set}";

  downloadAll($remoteDir,$settings);

  return 0;
}

sub downloadAll{
  my($remoteDir,$settings)=@_;

  my $ftp = Net::FTP->new($$settings{domain}, Debug => 0)
    or die "Cannot connect to $$settings{domain}: $@";
  $ftp->login("anonymous",'-anonymous@')
      or die "Cannot login ", $ftp->message;


  #Retrieve the metadata file
  logmsg "Retrieving metadata";
  $ftp->cwd("//pathogen/Results/$$settings{set}/$remoteDir/Metadata/")
    or die "Cannot change working directory ", $ftp->message;
  my @metafiles = $ftp->ls("*.tsv");
  foreach(@metafiles) {
    $ftp->get($_,"$$settings{tempdir}/$_")
      or die "get failed", $ftp->message;
  }

  #Retrieve SNP distances
  logmsg "Retrieving SNP distances";
  $ftp->cwd("//pathogen/Results/$$settings{set}/$remoteDir/Clusters/")
    or die "Cannot change working directory ", $ftp->message;
  my @distfiles = $ftp->ls("*.SNP_distances.tsv");
  foreach(@distfiles) {
    $ftp->get($_,"$$settings{tempdir}/$_")
      or die "get failed", $ftp->message;
  }

  #Retrieve SNP trees
  logmsg "Retrieving trees";
  $ftp->cwd("//pathogen/Results/$$settings{set}/$remoteDir/SNP_trees/")
    or die "Cannot change working directory ", $ftp->message;
  my @snpfiles = $ftp->ls();
  logmsg scalar(@snpfiles)." trees to download";
  for(my $i=0;$i<@snpfiles;$i++){
    $ftp->cwd($snpfiles[$i])
      or die "Cannot change working directory ", $ftp->message;
    my @newfiles = $ftp->ls("*.newick");
    $ftp->ascii();
    $ftp->get($newfiles[0],"$$settings{tempdir}/$newfiles[0]")
      or die "get failed", $ftp->message;
    $ftp->cwd("..")
      or die "Cannot change working directory ", $ftp->message;

    if($i % 10 == 0 && $i>0){
      logmsg "Finished downloading $i trees";
    }
  }

  $ftp->quit;

  return scalar(@snpfiles);
}

sub usage{
  "$0: downloads NCBI Pathogen Detection Pipeline results
  Usage: $0 latest
    where 'latest' is the Pathogen Detection Pipeline 
    directory from which to retrieve results

  --tempdir     /tmp      Where temporary files go including
                          trees, metadata, etc
  --resultsSet  Listeria  Which NCBI results set to download
  "
}
