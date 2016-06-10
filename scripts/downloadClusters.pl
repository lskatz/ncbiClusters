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
use Data::Dumper;
use POSIX qw/strftime/;
use Time::Piece; # for parsing dates from Metadata.tsv
use File::Copy qw/mv/;
use List::Util qw/min max/;

# Import local modules
use FindBin;
use lib "$FindBin::RealBin/../lib";
use Bio::Tree::Draw::Cladogram; # requires PostScript/TextBlock.pm in the lib dir
use Bio::TreeIO;
use PostScript::Simple;
#use Docopt;

local $0=basename($0);
sub logmsg { print STDERR "$0: @_\n"; }

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help report from=s to=s tempdir=s outdir=s set|results|resultsSet=s list maxTrees=i)) or die $!;
  $$settings{domain}||="ftp.ncbi.nlm.nih.gov";
  $$settings{tempdir}||=tempdir("$0.XXXXXX",TMPDIR=>1,CLEANUP=>1);
  $$settings{outdir}||="out";
  $$settings{set}||="Listeria";
  $$settings{maxTrees}||=0;

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

  if($$settings{list}){
    listSets($settings);
    return 0;
  }

  my $remoteDir=$ARGV[0];

  if(!$remoteDir || $$settings{help}){
    die usage($settings);
  }

  mkdir($$settings{tempdir}) if(!-e $$settings{tempdir});
  logmsg "Temporary directory is $$settings{tempdir}";
  logmsg "Downloading for $$settings{set}";

  downloadAll($remoteDir,$settings);
  my $metadata=readMetadata($settings);

  # Remove older trees and anything else that doesn't
  # pass the user-defined filters.
  filterTrees($metadata,$settings);

  makeReport($metadata,$settings) if($$settings{report});

  # Move anything over that passed the filter
  mkdir $$settings{outdir};
  for(glob("$$settings{tempdir}/*.{newick,tsv,pdf}")){
    mv($_,$$settings{outdir});
  }
  mkdir "$$settings{outdir}/images";
  for(glob("$$settings{tempdir}/*.{eps,gif}")){
    mv($_,"$$settings{outdir}/images");
  }

  return 0;
}

sub listSets{
  my($settings)=@_;

  my $ftp = Net::FTP->new($$settings{domain}, Debug => 0)
    or die "Cannot connect to $$settings{domain}: $@";
  $ftp->login("anonymous",'-anonymous@')
      or die "Cannot login ", $ftp->message;
  
  $ftp->cwd("//pathogen/Results")
    or die "Cannot change working directory ", $ftp->message;

  my @resultSets=$ftp->ls("");

  for(@resultSets){
    print $_."\n";
  }
  
  return \@resultSets;
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
  @snpfiles=reverse(@snpfiles);
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
    if($$settings{maxTrees} && $i >= $$settings{maxTrees}){
      last;
    }
  }

  $ftp->quit;

  return scalar(@snpfiles);
}

sub readMetadata{
  my($settings)=@_;

  # Assume the metadata file is the only one in the temp directory
  my @metadataFiles=glob("$$settings{tempdir}/*.metadata.tsv");
  my $infile=$metadataFiles[0];
  if(@metadataFiles > 1){
    logmsg "WARNING: there is more than one metadata file in $$settings{tempdir}! Assuming $infile.";
  }

  my %metadata;
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
    # The ID of this hash will be the PDT identifier
    # because it is the same one used in the newick trees
    $metadata{$F{target_acc}}=\%F;
  }
  close $metadataFh;

  return \%metadata;
}

sub filterTrees{
  my($metadata,$settings)=@_;

  for my $tree(glob("$$settings{tempdir}/*.newick")){
    # Read in the tree and its leaves
    my $treeObj=Bio::TreeIO->new(-file=>$tree)->next_tree;
    my @sample=$treeObj->get_leaf_nodes(); # returns node objects

    # Determine whether this tree has any isolates in the 
    # time window that the user supplied
    my $is_in_timewindow=0; # assume false until proven true by >0 isolates
    for my $s(@sample){
      my $target_acc=$s->id;
      $target_acc=~s/^\'|\'$//g; # remove single quotes that might appear around taxa in the tree
      
      # Read in the date to an object, assuming the date is in Y-m-d format.
      # The best date to read is the collection date, but if not that, then
      # the target_creation_date is a good assumption (although sometimes not
      # accurate for DoC).
      my $dateToRead=$$metadata{$target_acc}{collection_date};
      $dateToRead=~s/^\s+|\s+//g; # whitespace trim so I don't have to worry about it below
      if($dateToRead=~/^$|^missing$|^null$|^0$/i){ # regex for blank or 'missing' or 0 or 'null'
        $dateToRead=$$metadata{$target_acc}{target_creation_date};
      }
      my $dateObj=parseDate($dateToRead);
      if($dateObj >= $$settings{from} && $dateObj <= $$settings{to}){
        $is_in_timewindow=1;
        last;
      }
    }

    # Delete the tree if it doesn't occur within the time window.
    if($is_in_timewindow){
      logmsg "Tree passed: $tree";
    } else {
      logmsg "Tree is not in the time window: $tree";
      unlink($tree);
    }

    # TODO any future filters?
  }
}

sub makeReport{
  my($metadata,$settings)=@_;
  
  my $reportFile="$$settings{tempdir}/report.pdf";

  for my $tree(glob("$$settings{tempdir}/*.newick")){
    my $eps="$$settings{tempdir}/".basename($tree,'.newick').".eps";
    my $gif="$$settings{tempdir}/".basename($tree,'.newick').".gif";
    my $treeObj=Bio::TreeIO->new(-file=>$tree)->next_tree; # assume only one tree in the file
    #
    # Avoid a random divide by zero error in the cladogram module
    if($treeObj->get_root_node->height==0){
      #$treeObj->get_root_node->branch_length(1e-8);
      next;  # Skip; not sure what to do about this right now
    }

    addMetadataToTree($treeObj,$metadata,$settings);
    my $cladogram=Bio::Tree::Draw::Cladogram->new(
      -tree       => $treeObj,
      #-font       => "sans-serif",
      -size       => 12,
      -colors     => 1,
      -bootstrap  => 0,  # draw bootstraps?
      
      # no margins
      -top        => 0,  
      -right      => 0,
      -bottom     => 0,
      -left       => 0,
    );
    $cladogram->print(-file=>$eps);
    
    # Convert the eps to a gif for a larger report
    #system("convert -density 300 -resize '1024x1024' -colorspace RGB -flatten $eps $gif");
    #die "ERROR with imagemagick 'convert'" if $?;
  }

  # Make the full report
  my $p=new PostScript::Simple(
    papersize   => "Letter",
    colour      => 1,
    eps         => 0,
    units       => "pt",        # Must be in pt units instead of in because the tree eps files are in pt
    #coordorigin => "LeftTop",   # coordinate origin top-left
    #direction   => "RightDown", # direction for x-y coordinates
  );

  my $pageNumber=1;
  $p->newpage($pageNumber);
  $p->setcolour("black");
  $p->setfont("Times-Roman",32);
  $p->text({align=>"centre"},72*4,72*10,"Report from NCBI Pathogen Pipeline");
  $p->text({align=>"centre"},72*4,72*9.5,"generated today"); # TODO fill in the date
  #$p->setfont("Times-Roman",24);
  #$p->text(1,7,"Generated by EDLB/CFSAN software");

  # Make a new page per tree
  #$p->setcolour(30,30,30);
  for my $eps(glob("$$settings{tempdir}/*.eps")){
    $p->newpage(++$pageNumber); # the PS module increments page numbers automatically starting with 1
    #$p->{direction}="RightDown";
    #$p->{coordorigin} = "LeftTop";   # coordinate origin bottom-left
    $p->setfont("Times-Roman",16);
    $p->setcolour("black");
    logmsg "Writing $eps to file";

    $p->text({align=>"centre"},4*72,10*72,basename($eps,".newick_tree.eps"));
    if($p->err()){
      die $p->err();
    }

    # Bio::Tree::Draw::Cladogram puts decimals into the width
    # but PostScript::Simple cannot understand decimals.
    # Therefore I have to find the bbox dimensions and write
    # them correctly to a new temporary file, then replace
    # the original with the temporary file.
    my @treeCoordinates;  #x1,y1,x2,y2
    open(my $epsFh,"<",$eps) or die "ERROR: could not read $eps: $!";
    open(my $epsNewFh,">","$eps.tmp") or die "ERROR: could not write to $eps.tmp: $!";
    while(my $line=<$epsFh>){
      if ($line=~/^\%\%BoundingBox:\s+(.+)\s*$/){
        my $dimensions=$1;
        @treeCoordinates=split(/\s+/,$dimensions);
        $_=int($_) for(@treeCoordinates);

        print $epsNewFh '%%BoundingBox: '.join(" ",@treeCoordinates)."\n";
      } else {
        print $epsNewFh $line;
      }
    }
    close $epsFh;

    die "ERROR: could not find bbox dimensions of $eps" if(!defined($treeCoordinates[3]));

    # Cement the new bbox dimensions into the file by
    # making a move (rename) command.
    system("mv $eps.tmp $eps");

    # Set up the direction of flow to bottom to top so that the
    # tree image doesn't get reversed and
    # set up pt because the image is in pt units.
    #$p->{direction}="RightUp";
    #$p->{coordorigin} = "LeftBottom";   # coordinate origin bottom-left

    # 72 points per inch.
    #$_=sprintf("%0.4f",($_/72)) for(@treeCoordinates);
    $_=($_+=72*0.5) for(@treeCoordinates[0,2]); # shift X
    $_=($_+=72*0.5) for(@treeCoordinates[1,3]); # shift Y
    #@treeCoordinates=(200,200,400,400);

    $p->importepsfile($eps, @treeCoordinates);
    if($p->err()){
      die "ERROR with $eps:\n". $p->err();
    }
  }
  $p->output("report.ps");
  logmsg "Wrote $pageNumber pages to report.ps";
  
}
# A hacky way to parse a variety of dates
sub parseDate{
  my($date,$settings)=@_;
  if($date=~m|\d{1,2}/\d{1,2}/\d{2,4}|){
    return Time::Piece->strptime($date,"%m/%d/%Y");
  } elsif($date=~m|\d{2,4}\-\d{1,2}\-\d{1,2}|){
    return Time::Piece->strptime($date,"%Y-%m-%d");
  } else {

    # IF MMYYDD formats don't work, try just the year
    if($date=~/^\d{4}$/){
      return Time::Piece->strptime($date,"%Y");
    } elsif($date=~/^(\d{4})/){ # just try first four digits of a number as a year
      return Time::Piece->strptime($1,"%Y");
    }
      
  }

  logmsg "WARNING: I could not parse $date for dates. Using 12/31/1969.";
  return Time::Piece->strptime("%Y-%m-%d","12/31/1969");
}

sub addMetadataToTree{
  my($treeObj,$metadata,$settings)=@_;

  for my $node($treeObj->get_leaf_nodes){
    # remove leading and lagging quotes
    my $target_acc=$node->id;
    $target_acc=~s/^\'|\'$//g;
    logmsg $node->id ." -> $target_acc";

    # rename the node
    $node->id($$metadata{$target_acc}{label});

    # add RGB color
    my($red,$blue,$green)=(0.1,0.1,0.1); # almost black by default
    if($$metadata{$target_acc}{attribute_package}=~/environmental|food/i){
      $blue=0.9;
      $red=0.3;
      $green=0.3;
    } elsif($$metadata{$target_acc}{attribute_package}=~/clinical|host/i){
      $blue=0.3;
      $red=0.9;
      $green=0.3;
    }
    logmsg "added $red $green $blue to $target_acc";
    $node->add_tag_value("Rcolor",$red);
    $node->add_tag_value("Gcolor",$green);
    $node->add_tag_value("Bcolor",$blue);
  }
}
    
sub usage{
  my($settings)=@_;
  my $to=$$settings{to}->strftime("%m/%d/%Y");
  my $from=$$settings{from}->strftime("%m/%d/%Y");
  my $epoch=strftime("%m/%d/%Y",localtime(0));
  "$0: downloads NCBI Pathogen Detection Pipeline results
  Usage: $0 latest
    where 'latest' is the Pathogen Detection Pipeline 
    directory from which to retrieve results

  --tempdir     /tmp        Where temporary files go including
                            trees, metadata, etc
  --outdir      ./out       Where output files go
  --resultsSet  Listeria    Which NCBI results set to download
  --list                    List the options for results sets
  --report                  Create a report in outdir/report.pdf
                            and then exit.
  FILTERING
  --from        $from   Include trees with any isolates
                            as early as this date. Dates
                            can be in the format of either
                            YYYY-MM-DD or MM/DD/YYYY
  --to          $to   Include trees with isolates
                            as late as this date.
  --maxTrees    0           Download at most this many trees
  "
}
