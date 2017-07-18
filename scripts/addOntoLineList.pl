#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper qw/Dumper/;
use File::Temp qw/tempdir/;
use File::Basename qw/basename/;
use Getopt::Long qw/GetOptions/;
use List::Util qw/uniq/;
use Bio::DB::EUtilities;

local $0 = basename $0;
sub logmsg{print STDERR "$0: @_\n";}

exit main();

sub main{
  my $settings={};
  GetOptions($settings, qw(help taxon=s out|output|outfile=s linelist=s indir|input|in=s)) or die $!;
  die usage() if($$settings{help});
  for (qw(linelist out indir taxon)){
    $$settings{$_} //= die "ERROR: need --$_";
  }

  # TODO sort by version number to get the latest file
  my $metadataFile=(glob("$$settings{indir}/pathogen/Results/$$settings{taxon}/latest_snps/Metadata/*.metadata.tsv"))[0];
  die "ERROR: could not locate metadata file in $$settings{indir}/pathogen/Results/$$settings{taxon}/latest_snps/Metadata/*.metadata.tsv" if(!$metadataFile);
  die "ERROR: metadata file not found ($metadataFile)" if(!-e $metadataFile);
  my $clusterFile =(glob("$$settings{indir}/pathogen/Results/$$settings{taxon}/latest_snps/Clusters/*.SNP_distances.tsv"))[0];
  die "ERROR: could not locate cluster file in $$settings{indir}/pathogen/Results/$$settings{taxon}/latest_snps/Clusters/*.SNP_distances.tsv" if(!$clusterFile);
  die "ERROR: cluster file not found ($clusterFile)" if(!-e $clusterFile);


  my $pnIsolates=readLineList($$settings{linelist}, $settings);
  my $metadata=readMetadata($metadataFile,$settings);

  my $relatedBiosamples=findCloselyRelatedIsolates($pnIsolates, $clusterFile, $settings);

  printNewLineList($relatedBiosamples, $metadata, $settings);

  return 0;
}

sub readMetadata{
  my($infile,$settings)=@_;
  
  my %metadata;
  logmsg "Reading $infile";

  open(my $metadataFh,"<",$infile) or die "ERROR: could not read $infile: $!";
  # Grab the header and turn it into an array
  my $header=<$metadataFh>; chomp($header);
  $header=~s/^\s+|^#|\s+$//g;  # trim whitespace and leading hash
  my @header=split(/\t/,$header);
  # Read the values in the file to associate them with the header(column) names
  while(my $line=<$metadataFh>){
    chomp($line);
    my @field=split(/\t/,$line);
    my %F;
    @F{@header}=@field; # get an index of header_name => field_value
    
    # Add onto the metadata hash.
    for(@header){
      # The ID of this hash will be the PDT identifier
      # because it is the same one used in the newick trees.
      # Because there could potentially be multiple metadata
      # spreadsheets, use the "equals if not blank"
      # operator to set each field. NOTE: this means that
      # in case of conflict, the first spreadsheet
      # that is read will take priority in setting the value.
      $metadata{$F{biosample_acc}}{$_}||=$F{$_};
    }
  }
  close $metadataFh;
  
  return \%metadata;
}

sub readLineList{
  my($linelist,$settings)=@_;

  my %isolate;

  open(my $lineListFh, "<", $linelist) or die "ERROR: could not read $linelist: $!";
  # Grab the header and turn it into an array
  my $header=<$lineListFh>; chomp($header);
  $header=~s/^\s+|^#|\s+$//g;  # trim whitespace and leading hash
  my @header=split(/\t/,$header);
  # Some spreadsheet quality control
  if(@header != uniq(@header)){
    my %uniq;
    for(@header){
      if($uniq{$_}++){
        logmsg "Found duplicate header $_";
      }
    }
    logmsg "WARNING: some headers are duplicated in $linelist";
  }
  while(my $line=<$lineListFh>){
    $line=~s/^\s+|\s+$//g;
    my @field=split(/\t/,$line);

    # Set %F carefully in an if-defined method
    my %F;
    for(my $i=0;$i<@header;$i++){
      $F{$header[$i]}||=$field[$i]||"";
    }

    # Figure out what I want this key to be. 
    # Preferably, the BioSample accession
    my $index;
    for my $possibleIndexName(qw(NCBI_ACCESSION biosample_acc SAMN)){
        # Need to save this variable to avoid uninitialized value in hash within hash warning
        my $possibleIndexFromLine=$F{$possibleIndexName};
        next if(!defined($possibleIndexFromLine));

        $index=$F{$possibleIndexName};
        last;
    }
    # TODO try to use eutils to find a biosample accession if there isn't one
    # in the spreadsheet.
    if(!$index){
      $F{biosample_acc} = wgsIdToBiosample($F{WGS_id}, $settings);
      $index = $F{biosample_acc};
    }

    # If there isn't a biosample for the hash to be
    # index, then put out a warning and set undefined
    # fields to empty string.
    if(!$index){
      $_ ||= "" for(@field);
      logmsg "WARNING: Could not find a biosample accession in the line list for @field[0..2]";
      next;
    }

    # Copy over values carefully by checking if they are defined
    for my $key(keys(%F)){
      $isolate{$index}{$key} ||= $F{$key};
    }
  }
  close $lineListFh;
  return \%isolate;
}

sub findCloselyRelatedIsolates{
  my($pnIsolates, $clusterFile, $settings)=@_;

  logmsg "Finding closely related isolates from $clusterFile";

  my @filteredIsolate;

  open(my $clusterFh, $clusterFile) or die "ERROR: could not read $clusterFile: $!";
  my $header=<$clusterFh>;
  chomp($header);
  my @header=split(/\t/,$header);
  while(<$clusterFh>){
    my %F;
    @F{@header}=split(/\t/,$_);

    # At least one member must be in the line list
    if(!$$pnIsolates{$F{biosample_acc_1}} && !$$pnIsolates{$F{biosample_acc_2}}){
      next;
    }

    # These two isolates should be kept if they are close enough
    # to each other.
    if($F{compatible_distance} > 25){
      next;
    } 
    
    # Passed all filters: keep these two targets
    push(@filteredIsolate, $F{biosample_acc_1}, $F{biosample_acc_2});
  }
  close $clusterFh;
  
  # There are a ton of duplicates; remove them.
  @filteredIsolate=uniq(@filteredIsolate);
  return \@filteredIsolate;
}

sub printNewLineList{
  my ($relatedBiosamples, $metadata, $settings) = @_;
  #print Dumper $metadata; # SAMN02990449 SAMN03754731

  open(my $outFh, ">", $$settings{out}) or die "ERROR: could not write to $$settings{out}: $!";

  # Get the correct keys in order from metadata
  my($firstKey, $firstEntry)=each(%$metadata);
  my @metadataKey=keys(%$firstEntry);

  # Reorder some headers by plucking them out and then
  # placing them at the beginning
  for my $specificHeader(reverse qw(label Run attribute_package isolation_source bioproject_acc bioproject_center strain collection_date serovar)){
    @metadataKey = grep{!/^$specificHeader$/} @metadataKey;
    unshift(@metadataKey, $specificHeader);
  }

  print $outFh join("\t", @metadataKey)."\n";
  for my $biosample(@$relatedBiosamples){
    if(!$$metadata{$biosample}){
      logmsg "WARNING: metadata not found for $biosample";
      next;
    }
    for my $key(@metadataKey){
      print $outFh "$$metadata{$biosample}{$key}\t";
    }
    print $outFh "\n";
  }

  close $outFh;
}

# figure out the biosample accession from the WGS_id.
sub wgsIdToBiosample{
  my ($WGS_id, $settings) = @_;

  logmsg "Using edirect to find the biosample for $WGS_id";
  
  my $biosample = `esearch -query '$WGS_id' -db biosample | esummary | xtract -pattern DocumentSummary -element Accession`;
  if($?){
    die "ERROR: could not run esearch | esummary | xtract on $WGS_id";
  }
  $biosample//="";
  $biosample=~s/^\s+|\s+$//g;

  if(!$biosample){
    die "ERROR: could not extract biosample accession from WGS_id $WGS_id";
  }

  return $biosample;
}



sub usage{
  "$0: create a new line list of related isolates
  Usage: $0 --linelist in.tsv --indir ./ftp.ncbi.nlm.nih.gov/ --out newlinelist.tsv

  --indir       The input directory from NCBI. Is a mirror 
                of the original directory structure.
  --linelist    From PulseNet.  Must have biosample_acc as a header
  --out         The new line list.
  --taxon       The taxon as listed in the NCBI FTP site
  "
}
