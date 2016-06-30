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
use Data::Dumper;
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

local $0=basename($0);
sub logmsg { print STDERR "$0: @_\n"; }

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
  die "ERROR: need results set (ie, taxon) such as Listeria\n".usage() if(!$$settings{set});
  die "ERROR: need remote directory\n".usage() if(!$$settings{remoteDir});

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
  filterTrees($metadata,$settings);

  makeReport($metadata,$PDG,$settings);

  return 0;
}

sub listSets{
  my($settings)=@_;

  if($$settings{set}){
    return listRemoteDirs($settings);
  }

  my $ftp = Net::FTP->new($$settings{domain}, Debug => 0)
    or die "Cannot connect to $$settings{domain}: $@";
  $ftp->login("anonymous",'-anonymous@')
      or die "Cannot login ", $ftp->message;
  
  $ftp->cwd("//pathogen/Results")
    or die "Cannot change working directory ", $ftp->message;

  my @resultSets=$ftp->ls("");

  $ftp->quit;
  
  return \@resultSets;
}

sub listRemoteDirs{
  my($settings)=@_;
  my $ftp = Net::FTP->new($$settings{domain}, Debug => 0)
    or die "Cannot connect to $$settings{domain}: $@";
  $ftp->login("anonymous",'-anonymous@')
      or die "Cannot login ", $ftp->message;
  
  $ftp->cwd("//pathogen/Results/$$settings{set}")
    or die "Cannot change working directory ", $ftp->message;

  my @remoteDir=sort {$a cmp $b} $ftp->ls("");

  $ftp->quit;
  
  return \@remoteDir;
}

sub downloadAll{
  my($remoteDir,$settings)=@_;

  my $ftp = Net::FTP->new($$settings{domain}, Debug => 0)
    or die "Cannot connect to $$settings{domain}: $@";
  $ftp->login("anonymous",'-anonymous@')
      or die "Cannot login ", $ftp->message;
  $ftp->ascii(); # download ascii encoding


  #Retrieve the metadata file
  logmsg "Retrieving metadata";
  mkdir "$$settings{tempdir}/Metadata";
  mkdir "$$settings{outdir}/Metadata";
  $ftp->cwd("//pathogen/Results/$$settings{set}/$remoteDir/Metadata/")
    or die "Cannot change working directory. It is possible that '$$settings{set}' ($remoteDir) does not exist. ", $ftp->message;
  my @metafiles = $ftp->ls("*.tsv");
  foreach(@metafiles) {
    $ftp->get($_,"$$settings{tempdir}/Metadata/$_")
      or die "get failed", $ftp->message;
  }

  #Clusters directory files
  mkdir "$$settings{tempdir}/Clusters";
  mkdir "$$settings{outdir}/Clusters";
  #logmsg "Retrieving SNP distances";
  #$ftp->cwd("//pathogen/Results/$$settings{set}/$remoteDir/Clusters/")
  #  or die "Cannot change working directory ", $ftp->message;
  #my @distfiles = $ftp->ls("*.SNP_distances.tsv");
  #foreach(@distfiles) {
  #  $ftp->get($_,"$$settings{tempdir}/Clusters/$_")
  #    or die "get failed", $ftp->message;
  #}
  #
  # Retrieve the list of newest isolates, also in this directory
  logmsg "Retrieving a list of the latest isolates";
  $ftp->cwd("//pathogen/Results/$$settings{set}/$remoteDir/Clusters/")
    or die "Cannot change working directory to Clusters. It is possible that '$$settings{set}' ($remoteDir) does not exist. ", $ftp->message;
  my @newIsolateFiles=$ftp->ls("*.new_isolates.tsv")
    or die "ERROR with command while in Clusters directory: ls *.new_isolates.tsv: ", $ftp->message;
  foreach(@newIsolateFiles){
    $ftp->get($_,"$$settings{tempdir}/Clusters/$_")
      or die "get failed", $ftp->message;
  }
  

  #Retrieve SNP trees
  logmsg "Reading the remote SNP_trees directory";
  mkdir "$$settings{tempdir}/SNP_trees";
  mkdir "$$settings{outdir}/SNP_trees";
  $ftp->cwd("//pathogen/Results/$$settings{set}/$remoteDir/SNP_trees/")
    or die "Cannot change working directory ", $ftp->message;
  # Have to ls on the local folder because it would be too many results
  # for the server to return.  Cannot ls on */*.newick.
  my @SNP_dir = $ftp->ls("")
    or die "ERROR: cannot get directory contents: ", $ftp->message;
  die "INTERNAL ERROR: I only saw ".scalar(@SNP_dir) if(@SNP_dir < 1);
  @SNP_dir=reverse(@SNP_dir);        # Reverse order so that we get mostly new trees first
  @SNP_dir=grep(/^PDS\d+/,@SNP_dir); # Keep only directories that begin with PDS
  # TODO filter trees to download, to save tons of time. Do the filtering
  # step here instead of after all trees are downloaded.
  logmsg scalar(@SNP_dir)." trees to download";
  for(my $i=0;$i<@SNP_dir;$i++){
    if($i % 100 == 0 && $i>0){
      logmsg "Finished downloading $i trees out of ".scalar(@SNP_dir);
    }

    # If it is already here and downloaded, then don't
    # download again.
    my $localtree="$$settings{tempdir}/SNP_trees/$SNP_dir[$i].newick_tree.newick"; 
    next if(-e $localtree);

    # Look for the tree file
    my @cluster_files=$ftp->ls($SNP_dir[$i]) or die "ERROR: could not read $SNP_dir[$i]: ",$ftp->message;
    my $remotetree=(grep(/\.newick$/,@cluster_files))[0];
    # there should be a newick tree in each dir, but don't bother downloading
    # a null tree in case it happens.
    next if(!$remotetree);

    # Download the tree file.
    $ftp->get($remotetree,$localtree)
      or die "get failed", $ftp->message;

    if($$settings{maxTrees} && $i >= $$settings{maxTrees}){
      last;
    }
  }

  $ftp->quit;

  return scalar(@SNP_dir);
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
      for my $possibleIndexName(qw(NCBI_ACCESSION biosample_acc SAMN)){
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
        logmsg "WARNING: Could not find a biosample accession for \n".Dumper \%F;
        next;
      }

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
  for my $target_acc(keys(%metadata)){
    my $newLabel="$target_acc | ";
    for my $field(@labelFields){
      next if(!defined($metadata{$target_acc}{$field}));
      next if($metadata{$target_acc}{$field}=~/^NULL$|^\s*$|^missing$/);
      $newLabel.=$metadata{$target_acc}{$field}.' | ';
    }
    $newLabel=~s/ \| $//; # remove trailing pipe as a result of the above loop
    $metadata{$target_acc}{label}=$newLabel;
  }

  return \%metadata;
}

sub filterTrees{
  my($metadata,$settings)=@_;

  # Index the list of new isolates if requested to 
  # filter by new isolates
  my %treeWithNewIsolate;
  if($$settings{'new-isolates'}){
    logmsg "Filtering trees by new isolates, identified by $$settings{tempdir}/Clusters/*.new_isolates.tsv";
    my $newisolatesFile=(glob("$$settings{tempdir}/Clusters/*.new_isolates.tsv"))[0];
    open(my $fh,"<",$newisolatesFile) or die "ERROR: cannot read $newisolatesFile: $!";
    my $header=<$fh>; chomp($header);
    my @header=split(/\t/,$header);
    while(<$fh>){
      chomp;
      my @F=split(/\t/,$_);
      my %F;
      @F{@header}=@F;
      $treeWithNewIsolate{$F{PDS_acc}}=\%F;
    }
    close $fh;
  }

  # Start the filtering process
  logmsg "Filtering trees by time and/or by line list. If none were given, then no trees will be filtered by this process.";
  for my $tree(glob("$$settings{tempdir}/SNP_trees/*.newick")){
    # If the tree passes filters, then it gets moved here
    my $outtree="$$settings{outdir}/SNP_trees/".basename($tree);
    my $tree_passed_filters=1;  # innocent until guilty
    my @whyFail; #Keep track of why a tree didn't pass filters for debugging

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

    # If it doesn't occur within the time window.
    if($is_in_timewindow){
      #logmsg "Tree passed: $tree";
    } else {
      #logmsg "Tree is not in the time window: $tree";
      $tree_passed_filters=0;
      push(@whyFail,"Not in time window");
    }
    
    # Filter for new isolates only if requested
    if($$settings{'new-isolates'}){
      my $PDS_acc=basename($tree,".newick_tree.newick");
      if($treeWithNewIsolate{$PDS_acc}){
        
      } else {
        $tree_passed_filters=0;
        push(@whyFail,"Tree does not have isolates in the new isolates list");
      }
    }  

    # Only keep trees in the line lists, if line lists
    # were supplied.
    if(@{$$settings{'line-list'}} > 0){
      
      my $tree_has_linelist=0;
      for my $s(@sample){
        my $target_acc=$s->id;
        $target_acc=~s/^\'|\'$//g; # remove single quotes that might appear around taxa in the tree
        
        # Check for anything that pulsenet would put into
        # the metadata that NCBI wouldn't have. That will
        # be the indication that it is in the line list.
        # TODO: I should probably just read the line list
        # directly.
        if(
             defined($$metadata{$target_acc}{Key}) # state ID
          || defined($$metadata{$target_acc}{SourceCity}) 
          || defined($$metadata{$target_acc}{SourceCounty})
          || defined($$metadata{$target_acc}{PatientSex})
          || defined($$metadata{$target_acc}{SeroType})
          || defined($$metadata{$target_acc}{PatientAge})
        ){
          $tree_has_linelist=1;
          last;
        }
      }
      if($tree_has_linelist==0){
        $tree_passed_filters=0;
        push(@whyFail,"Tree does not have isolates in the line list");
      }
    }

    # TODO any future filters?

    # If the tree passed all filters, then copy it over
    if($tree_passed_filters){
      cp($tree,$outtree) or die "ERROR copying $tree to $tree:\n  $!";
    } else {
      #logmsg "Tree failed the filters because ".join("\n",@whyFail);
    }
  }

  logmsg "Done filtering";

}

sub makeReport{
  my($metadata,$PDG,$settings)=@_;

  logmsg "Making the PDF report";
  
  my $postscript="$$settings{outdir}/report.ps";
  my $PDF="$$settings{outdir}/report.pdf";

  # Figure out a coloring scheme
  my @colorPercentage=(0.1,0.3,0.5,0.7,0.9);
  my @availableColor = variations_with_repetition(\@colorPercentage,3);
  # Sort colors brightest to darkest instead of by RGB. Randomize colors
  # that have an equal score in the sort.
  @availableColor=sort { sum(map{exp($_)} @$b) <=> sum(map{exp($_)} @$a) } shuffle(@availableColor);
  # Remove colors that are too bright
  @availableColor=grep{!($$_[0]>0.7 && $$_[1] >0.7 && $$_[2] > 0.7)} @availableColor;
  #die Dumper [map{join(", ",@$_)} @availableColor];

  # Figure out what categories the user wants to color by
  my $colorBy=new Config::Simple("config/colorBy.ini")->get_block("global");
  my @colorBy=keys(%$colorBy);
  my %colorCoding=(''=>[0,0,0]); # holds color coding combinations defined by the config file
  $colorCoding{''} = pop(@availableColor); # when there is no color, choose black

  logmsg "Creating phylogeny images";
  mkdir "$$settings{outdir}/images";
  for my $tree(glob("$$settings{outdir}/SNP_trees/*.newick")){
    my $PDS=basename($tree,".newick_tree.newick");
    my $eps="$$settings{outdir}/images/".basename($tree,'.newick').".eps";
    my $treeObj=Bio::TreeIO->new(-file=>$tree)->next_tree; # assume only one tree in the file
    
    # Avoid a random divide by zero error in the cladogram module
    if($treeObj->get_root_node->height==0){
      #$treeObj->get_root_node->branch_length(1e-8);
      next;  # Skip; not sure what to do about this right now
    }

    # Adds colors into the tree object (and in the future, anything else).
    addMetadataToTree($treeObj,$metadata,$colorBy,\%colorCoding,\@availableColor,$settings);

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
  $p->setfont("Times-Roman",24);
  $p->text({align=>"left"},72*0.5,72*9.5,"NCBI dataset: $$settings{set} ($$settings{remoteDir})");
  $p->text({align=>"left"},72*0.5,72*9.0,"Filters:");
  $p->setfont("Times-Roman",18);
  $p->text({align=>"left"},72*0.5,72*8.5,"  from: $$settings{from}");
  $p->text({align=>"left"},72*0.5,72*8.0,"  to: $$settings{to}");

  # Add a color legend in a black box 7"x2"
  $p->setcolour("black");
  $p->setlinewidth(1);
  my $wholeLegendHeight=4*72;
  my $wholeLegendY1=4*72;
  $p->box(72*0.5,$wholeLegendY1,72*7.5,72*4+$wholeLegendHeight);

  # The height/width of each color legend box is 2". For
  # keeping a margin, we will double the effective size
  # and add one.
  my $ptPerColor=($wholeLegendHeight/(1+2*scalar(keys(%colorCoding))));
  my $legendYMarker=$wholeLegendY1; # the marker for where legend boxes are starts with the outer box
  # Alphabetize the legend
  for my $category(sort {$b cmp $a} keys(%colorCoding)){
    my $color=$colorCoding{$category} || [0,0,0];
    $legendYMarker+=$ptPerColor; # adding margin
    # Convert the decimal color used in BioPerl into the 0-255 range for PostScript::Simple.
    $p->setcolour($$color[0]*255, $$color[1]*255, $$color[2]*255);
    $p->box({filled=>1},
      72*1,                       # x1
      $legendYMarker,             # y1
      72*1.0+$ptPerColor,         # x2
      $legendYMarker+$ptPerColor, # y2
    );

    $p->setfont("Times-Roman",int($ptPerColor));
    $p->setcolour("black");
    $category||="No info";
    #$p->text(72*1.0+$ptPerColor*2, $legendYMarker, $category."  ". join(", ",@$color)); # To the right of the box
    $p->text(72*1.0+$ptPerColor*2, $legendYMarker, $category); # To the right of the box

    $legendYMarker+=$ptPerColor; # advance marker past the legend box
  }

  # Footer for time generaged
  $p->setcolour("black");
  $p->setfont("Times-Roman",16);
  $p->text({align=>"centre"},72*4,72*1,"generated ".localtime()); 

  # Make a new page per tree
  #$p->setcolour(30,30,30);
  logmsg "Adding phylogeny images to larger report";
  for my $eps(glob("$$settings{outdir}/images/*.eps")){
    $p->newpage(++$pageNumber); # the PS module increments page numbers automatically starting with 1
    #$p->{direction}="RightDown";
    #$p->{coordorigin} = "LeftTop";   # coordinate origin bottom-left
    $p->setfont("Times-Roman",16);
    $p->setcolour("black");
    #logmsg "Writing to $eps";

    my $tree_acc=basename($eps,".newick_tree.eps");
    $p->text({align=>"centre"},4*72,10.5*72,$tree_acc);
    if($p->err()){
      die $p->err();
    }

    # Make URL to NCBI results directory
    $p->setfont("Times-Roman",12);
    $p->text({align=>"left"},72*0.3,72*10,"http://www.ncbi.nlm.nih.gov/pathogens/$$settings{set}/$PDG/$tree_acc");
    $p->text({align=>"left"},72*0.3,72*9.7,"ftp://ftp.ncbi.nlm.nih.gov/pathogen/Results/$$settings{set}/$$settings{remoteDir}/SNP_trees/$tree_acc");

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

    # Make a 0.5 inch margin from the bottom left
    # 72 points per inch.
    $_=($_+=72*0.5) for(@treeCoordinates[0,2]); # shift X
    $_=($_+=72*0.5) for(@treeCoordinates[1,3]); # shift Y

    # Rescale if wider than the page or taller than the page.
    # Dimensions inside of the margin: 7.5" x 10"
    if($treeCoordinates[2] > 7.5*72 || $treeCoordinates[3] > 10*72){
      my $scale=min(7.5*72/$treeCoordinates[2], 10*72/$treeCoordinates[3]);
      $_ *= $scale for(@treeCoordinates);
    }

    $p->importepsfile($eps, @treeCoordinates);
    if($p->err()){
      die "ERROR with $eps:\n". $p->err();
    }
  }
  $p->output($postscript);
  logmsg "Wrote $pageNumber pages to $postscript";
  
  system("ps2pdf $postscript $PDF");
  if($?){
    logmsg "WARNING: could not use `ps2pdf` to convert $postscript to $PDF";
  } else {
    logmsg "Converted $postscript to $PDF";
  }
}
# A hacky way to parse a variety of dates
sub parseDate{
  my($date,$settings)=@_;

  return Time::Piece->strptime("12/31/1969","%m/%d/%Y") if($date=~/^\s*$/);

  # Don't die just because Time::Piece craps out
  my $timePiece;
  eval{
    if($date=~m|\d{1,2}/\d{1,2}/\d{2,4}|){
      $timePiece=Time::Piece->strptime($date,"%m/%d/%Y");
    } elsif($date=~m|\d{2,4}\-\d{1,2}\-\d{1,2}|){
      $timePiece=Time::Piece->strptime($date,"%Y-%m-%d");
    } else {

      # IF MMYYDD formats don't work, try just the year
      if($date=~/^\d{4}$/){
        $timePiece=Time::Piece->strptime($date,"%Y");
      } elsif($date=~/^(\d{4})/){ # just try first four digits of a number as a year
        $timePiece=Time::Piece->strptime($1,"%Y");
      }
        
    }
    return $timePiece;
  };
  return $timePiece if(defined $timePiece);

  logmsg "WARNING: I could not parse $date for dates. Using 12/31/1969.";
  return Time::Piece->strptime("12/31/1969","%m/%d/%Y");
}

sub addMetadataToTree{
  my($treeObj,$metadata,$colorBy,$colorCoding,$availableColor,$settings)=@_;

  for my $node($treeObj->get_leaf_nodes){
    # remove leading and lagging quotes
    my $target_acc=$node->id;
    $target_acc=~s/^\'|\'$//g;
    #logmsg $node->id ." -> $target_acc";

    # rename the node
    $node->id($$metadata{$target_acc}{label});

    # Color by user-defined categories
    my $color;
    # Currently this loop only works if there is only one category in the ini file.
    # TODO allow for more categories.
    for my $colorKey (keys(%$colorBy)){
      $$metadata{$target_acc}{$colorKey}//="";
      if($$metadata{$target_acc}{$colorKey} =~/^NULL$|^\s*$|^missing$/){
        next;
      }

      $color=$$colorCoding{ $$metadata{$target_acc}{$colorKey} };
      if(!$color){
        $color=shift(@$availableColor);
        $$colorCoding{ $$metadata{$target_acc}{$colorKey} } = $color;
      }
    }

    if(!$color){
      $color=$$colorCoding{''};
    }

    $node->add_tag_value("Rcolor",$$color[0]);
    $node->add_tag_value("Gcolor",$$color[1]);
    $node->add_tag_value("Bcolor",$$color[2]);
  }
}
    
sub usage{
  my($settings)=@_;
  my $to=$$settings{to}->strftime("%m/%d/%Y");
  my $from=$$settings{from}->strftime("%m/%d/%Y");
  "$0: downloads NCBI Pathogen Detection Pipeline results
  Usage: $0 resultsSet remoteDir
    where 'remoteDir' is the Pathogen Detection Pipeline 
    directory from which to retrieve results
    and 'resultsSet' is the taxon to download, e.g., Listeria

  --tempdir     /tmp        Where temporary files go including
                            trees, metadata, etc. Useful for
                            running this pipeline multiple times
                            and downloading only once.
  --outdir      ./out       Where output files go
  --list                    List the options for results sets, 
                            i.e., taxa.
                            If a taxon/resultsSet parameter
                            is already given, then it will list
                            all possible remoteDirs.

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
