#!/usr/bin/env perl 
#
# AUTHORS: Lee Katz <lkatz@cdc.gov> and Errol Strain <Errol.Strain@fda.hhs.gov>
# Run this script with --help for usage information.

# Need Perl version 5.12 or greater
require 5.12.0;

use strict;
use warnings;
use Getopt::Long;
use Net::FTP;
use File::Basename qw/basename/;
use File::Temp qw/tempdir/;
use POSIX qw/strftime/;
use File::Copy qw/mv cp/;
use List::Util qw/min max sum shuffle/;
use List::MoreUtils qw/uniq/;

# Import local modules
use FindBin;
use lib "$FindBin::RealBin/../lib/perl5";
use Bio::Tree::Draw::Cladogram; # requires PostScript/TextBlock.pm in the lib dir and so this is 'local'
use Bio::TreeIO;
use PostScript::Simple;
use Time::Piece; # for parsing dates from Metadata.tsv
use Config::Simple;
use Algorithm::Combinatorics qw/variations_with_repetition/;
use NcbiClusters qw/logmsg listSets downloadAll parseDate Dumper colorScheme randColor/;

local $0=basename($0);

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help new-isolates from=s to=s tempdir=s outdir=s list maxTrees=i line-list=s@)) or die $!;
  $$settings{domain}||="ftp.ncbi.nlm.nih.gov";
  $$settings{tempdir}||=tempdir("$0.XXXXXX",TMPDIR=>1,CLEANUP=>1);
  $$settings{outdir}||="out";
  $$settings{maxTrees}||=0;
  $$settings{'line-list'}||=[];

  # Parse the date parameters
  if($$settings{from}){
    $$settings{from}=parseDate($$settings{from});
  } else {
    $$settings{from}=parseDate(strftime("%m/%d/%Y",localtime(0)));
  }
  if($$settings{to}){
    $$settings{to}=parseDate($$settings{to});
  } else {
    $$settings{to}||=parseDate(strftime("%m/%d/%Y",localtime()));
  }

  die usage($settings) if($$settings{help});

  ($$settings{set},$$settings{remoteDir})=@ARGV;

  if($$settings{list}){
    my $list=listSets($settings);
    for(@$list){
      print $_."\n";
    }
    return 0;
  }

  # Can only die on error only after the program has
  # had a chance to run listSets().
  die "ERROR: need taxon such as Listeria\n".usage($settings) if(!$$settings{set});
  die "ERROR: need remote directory\n".usage($settings) if(!$$settings{remoteDir});

  mkdir($$settings{tempdir}) if(!-e $$settings{tempdir});
  logmsg "Temporary directory is $$settings{tempdir}";
  logmsg "Downloading for $$settings{set}";
  mkdir($$settings{outdir}) if(!-e $$settings{outdir});

  downloadAll($$settings{remoteDir},$settings);

  # Sort the metadata by order of latest result first.
  # This should only be one file unless the tmp folder
  # has been used more than once.
  my @metadataFiles=sort { 
      my $aInt=basename($a);
      my $bInt=basename($b);
      $aInt=~s/PDG\d+\.(\d+)\.metadata.tsv/$1/;
      $bInt=~s/PDG\d+\.(\d+)\.metadata.tsv/$1/;
      return $bInt <=> $aInt;
    } glob("$$settings{tempdir}/Metadata/*.metadata.tsv");
  my $PDG=basename($metadataFiles[0],".metadata.tsv");
  my $metadata=readMetadata($metadataFiles[0],$settings);

  # Copies anything that passes the filters to the 
  # output directory.
  my $filteredIsolates = filteredIsolates($metadata,$settings);
  
  printReport($filteredIsolates,$settings);
  downloadSra($filteredIsolates, $settings);
  
  return 0;
}

sub readMetadata{
  my($infile,$settings)=@_;
  
  # Be able to map back and forth from pdt to biosample
  my %pdt_biosample;
  my %biosample_pdt;

  my %metadata;
  logmsg "Reading $infile";
  # Make this file available to the user but read from the temp
  # file because the user would probably rather mess with
  # the output one (however unlikely).
  cp($infile,"$$settings{outdir}/Metadata/");

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
      $metadata{$F{target_acc}}{$_}||=$F{$_};
    }

    if($F{target_acc} && $F{biosample_acc}){
      $pdt_biosample{$F{target_acc}}=$F{biosample_acc};
      $biosample_pdt{$F{biosample_acc}}=$F{target_acc};
    }
  }
  close $metadataFh;

  # Read the config file for spreadsheet headers
  my $headerMappings=new Config::Simple("config/headerMappings.ini")->get_block("header_mappings");

  # Added information from PulseNet line lists
  for my $infile(@{$$settings{'line-list'}}){
    open(my $lineList, "<", $infile) or die "ERROR: could not read line list $infile: $!";
    # Grab the header and turn it into an array
    my $header=<$lineList>; chomp($header);
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
      die "ERROR: some headers are duplicated in $infile";
    }
    while(my $line=<$lineList>){
      chomp $line;
      my @field=split(/\t/,$line);
      my %F;
      @F{@header}=@field;

      # Figure out what I want this key to be. 
      # Preferably, the PDT identifier (target_acc)
      my $index;
      for my $possibleIndexName(qw(target_acc NCBI_ACCESSION biosample_acc SAMN strain sample_name)){
        # Need to save this variable to avoid uninitialized value in hash within hash warning
        my $possibleIndexFromLine=$F{$possibleIndexName};
        next if(!defined($possibleIndexFromLine));

        $index=$F{$possibleIndexName};
        # Don't break the loop here because we'd prefer the biosample_pdt index
        
        # Although the index can be found in the line list,
        # it'd be better to find the PDT value in the NCBI metadata.
        if(defined($biosample_pdt{$possibleIndexFromLine})){
          $index=$biosample_pdt{$possibleIndexFromLine};
          last;
        }
      }
      # TODO try to use eutils to find a biosample accession if there isn't one
      # in the spreadsheet.
      if(!$index){
        logmsg "WARNING: Could not find a biosample accession in the line list for \n".Dumper \%F;
        next;
      }

      # Since there is a match on metadata from the line list,
      # mark this isolate as having come from the line list.
      $metadata{$index}{'is_in_lineList'}=1;

      # Add onto the large metadata hash.
      for my $pnHeader(@header){
        # If this header can be translated to an NCBI header,
        # do it.
        if(my $ncbiHeader=$$headerMappings{$pnHeader}){
          $metadata{$index}{$ncbiHeader}=$F{$pnHeader};
        }
        # Add the pulsenet value under the pulsenet header
        $metadata{$index}{$pnHeader}=$F{$pnHeader};
      }     
    }
    close $lineList;
  }

  # Get the rules for changing the label for each sample
  my $customLabel=new Config::Simple("config/headerMappings.ini")->get_block("custom_label");
  my @labelFields=split(/\s*\|\s*/,$$customLabel{label});
  # For each isolate, modify its label
  for my $target_acc(keys(%metadata)){
    my $newLabel="$target_acc | ";
    for my $field(@labelFields){
      if(
           !defined($metadata{$target_acc}{$field})
        || $metadata{$target_acc}{$field}=~/^NULL$|^\s*$|^missing$/i
      ){
        $metadata{$target_acc}{$field}="";
      }
      $newLabel.=$metadata{$target_acc}{$field}.' | ';
    }
    $newLabel=~s/ \| $//; # remove trailing pipe as a result of the above loop
    $metadata{$target_acc}{label}=$newLabel;
  }

  return \%metadata;
}

sub filteredIsolates{
  my($metadata,$settings)=@_;

  logmsg "Filtering by time and/or by line list. If none were given, then no trees will be filtered by this process. If --new-isolates was given, we are also filtering by that!";
  
  my %filteredIsolate;
  
  # Index the list of new isolates if requested to 
  # filter by new isolates
  if($$settings{'new-isolates'}){
  # TODO sort by the version number
    my $newisolatesFile=(glob("$$settings{tempdir}/Clusters/*.new_isolates.tsv"))[0];
    open(my $fh,"<",$newisolatesFile) or die "ERROR: cannot read $newisolatesFile: $!";
    my $header=<$fh>; chomp($header);
    my @header=split(/\t/,$header);
    while(<$fh>){
      chomp;
      my @F=split(/\t/,$_);
      my %F;
      @F{@header}=@F;
      $filteredIsolate{$F{target_acc}}=$$metadata{$F{target_acc}};
    }
    close $fh;
  }
  
  # open the cluster information but only keep recent isolates so that
  # the memory footprint is low
  # TODO sort by version number and get the latest.
  my $clusterFile = (glob("$$settings{tempdir}/Clusters/*.SNP_distances.tsv"))[0];
  open(my $clusterFh, $clusterFile) or die "ERROR: could not read $clusterFile: $!";
  my $header=<$clusterFh>;
  chomp($header);
  my @header=split(/\t/,$header);
  while(<$clusterFh>){
    my %F;
    @F{@header}=split(/\t/,$_);

    # At least one member must be in the line list
    if(!$$metadata{$F{target_acc_1}}{'is_in_lineList'} && 
       !$$metadata{$F{target_acc_2}}{'is_in_lineList'}
      ){
      next;
    }

    die Dumper \%F;
    
    # Discard this pairwise comparison if either isolate is older than
    # X days.
    my $date_1 = parseDate($$metadata{$F{target_acc_1}}{target_creation_date});
    my $date_2 = parseDate($$metadata{$F{target_acc_2}}{target_creation_date});
    if(
      $date_1 < $$settings{from} || $date_1 > $$settings{to} ||
      $date_2 < $$settings{from} || $date_2 > $$settings{to}
      ){
      next;
    }
    
    # These two isolates should be kept if they are close enough
    # to each other.
      if($F{compatible_distance} > 50){
      next;
    }  
    
    # Passed all filters: keep these two targets
    $filteredIsolate{$F{target_acc_1}}=$$metadata{$F{target_acc_1}};
    $filteredIsolate{$F{target_acc_2}}=$$metadata{$F{target_acc_2}};
  }
  close $clusterFh;
  return \%filteredIsolate;
}

sub printReport{
  my($filteredIsolates,$settings)=@_;
  
  my $report="$$settings{outdir}/ncbiLineList.tsv";
  open(my $fh, ">", $report) or die "ERROR: could not write to $report: $!";
  # Make the new line list
  my @header=qw(label  HHS_region      LibraryLayout   PFGE_PrimaryEnzyme_pattern      PFGE_SecondaryEnzyme_pattern    Platform        Run     asm_acc asm_display_name        asm_level       asm_stats_contig_n50    asm_stats_length_bp     asm_stats_n_contig   assembly_method attribute_package       bioproject_acc  bioproject_center       bioproject_title        biosample_acc   collected_by    collection_date complete_fl     fullasm_id      geo_loc_name    host    host_disease    isolation_source     lat_lon outbreak        sample_name     scientific_name serovar species_taxid   sra_center      sra_release_date        strain  sub_species     target_acc      target_creation_date    tax-id  wgs_acc_prefix  wgs_master_acc);
  print $fh join("\t", @header)."\n";
  for my $isolate(keys(%$filteredIsolates)){
    for my $h(@header){
      print $fh $$filteredIsolates{$isolate}{$h}."\t";
    }
    print $fh "\n";
  }
  close $fh;
}

sub downloadSra{
  my($filteredIsolates, $settings)=@_;
  
  my $outdir="$$settings{outdir}/sra";
  my $tempdir="$$settings{tempdir}/sra"; 
  mkdir $tempdir;
  mkdir $outdir;

  # option(s) for fastq-dump
  my $deflineSeq='@$ac.$si/$ri';
  
  # Make a SneakerNet-style folder by adding in a SampleSheet.csv
  my $sh = "$$settings{outdir}/sra/SampleSheet.csv";
  open(my $shFh, ">", $sh) or die "ERROR: could not write to $sh: $!";
  print $shFh "[Data]\n";
  print $shFh join("\t", qw(Sample_ID Sample_Name Sample_Plate Sample_Well I7_Index_ID index I5_Index_ID index2 Sample_Project Description))."\n";
  
  # Download isolates and also add in samplesheet entries
  my $isolateCounter=0;
  while(my($name, $isolate) = each(%$filteredIsolates)){
    if(!$$isolate{Run}){
      logmsg "Warning: could not find SRA entry for $$isolate{label}";
      my $emptyFile="$outdir/$$isolate{strain}.null.fastq";
      open(my $fh, ">", $emptyFile) or die "ERROR: could not write to $emptyFile: $!";
      close $fh;
      system("gzip -f '$emptyFile'"); # make it compatible
    }
    
    # Download into the tmp dir and then move it over.
    # That way, we can check if the file was downloaded completely.
    # Check if the file already exists
    elsif(-e "$outdir/$$isolate{Run}_1.fastq.gz"){
      logmsg "$$isolate{label} ($$isolate{Run}) was already downloaded; will not download again";
    } else {
      # Download into the temp dir
      logmsg "Downloading $$isolate{label} into $$isolate{Run}";
      
      system("fastq-dump --gzip --accession $$isolate{Run} --outdir $tempdir --defline-seq $deflineSeq --defline-qual '+' --split-files --skip-technical --dumpbase --clip");
      if($?){
  die "ERROR: could not download $$isolate{label} using SRA ID $$isolate{Run}";
      }
      for(glob("$tempdir/*")){
  mv($_, $outdir);
      }
    }
    
    my $scientificName=$$isolate{scientific_name};
    $scientificName=~s/\s+/_/g;
    print $shFh join("\t", $$isolate{Run}, ($$isolate{isolate_name} || $$isolate{strain} || 'strain'),
                           "DownloadNcbiClustersScript", 
         ++$isolateCounter, 'fake_index_id', 'fake_index',
         'fake_index_id', 'fake_index', 'fake_project', 
         "Species=$scientificName;Route=CalcEngine");
    print $shFh "\n";
  }
  
  close $shFh;
}

sub usage{
  my($settings)=@_;
  my $to=$$settings{to}->strftime("%m/%d/%Y");
  my $from=$$settings{from}->strftime("%m/%d/%Y");
  "$0: downloads NCBI Pathogen Detection Pipeline results
  Usage: $0 taxon remoteDir
    where 'remoteDir' is the Pathogen Detection Pipeline 
    directory from which to retrieve results

  --tempdir     /tmp        Where temporary files go including
                            trees, metadata, etc. Useful for
                            running this pipeline multiple times
                            and downloading only once.
  --outdir      ./out       Where output files go
  --list                    List the options for results sets, 
                            i.e., taxa, and then exit.
                            If a taxon parameter is already 
                            given, then it will list all 
                            possible remoteDirs.

  FILTERING
  --line-list   ''          A tab-delimited spreadsheet of a
                            PulseNet line list. Only trees
                            containing these isolates will
                            be reported. The headers may
                            not contain special characters
                            such as '#'.
                            Multiple --line-list flags are
                            allowed; one per line list.
  --from        $from   Include trees with any isolates
                            as early as this date. Dates
                            can be in the format of either
                            YYYY-MM-DD or MM/DD/YYYY
  --to          $to   Include trees with isolates
                            as late as this date.
  --new-isolates            Only include trees that contain
                            newly uploaded isolates.
  --maxTrees    0           Download at most this many trees
  "
}
