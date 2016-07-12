#!/usr/bin/env perl
package NcbiClusters;
use strict;
use warnings;
use Data::Dumper qw/Dumper/;
# Alphabetize data dump keys
$Data::Dumper::Sortkeys=sub{my($hash)=@_; return [sort {lc $a cmp lc $b} keys(%$hash)]};

use Net::FTP;
use Time::Piece;

use Exporter qw/import/;

our @EXPORT_OK = qw(logmsg listSets downloadAll parseDate findBiosample Dumper colorScheme randColor);

# Keep track of the colors that are used
our @_colorScheme=();
our %_colorScheme=();

sub logmsg { print STDERR "$0: @_\n"; }

=head1 NcbiClusters

A module for reporting on NCBI SNP clusters

=head1 AUTHORS

Lee Katz L<lkatz@cdc.gov>

Errol Strain L<Errol.Strain@fda.hhs.gov>

=head2 Methods

=over 12

=item C<Dumper>

  Runs a Data Dump on a variable using L<Data::Dumper>

=back 

=cut

=over 12

=item C<listSets>

  List the possible values for sets of data that can be
  retrieved. If C<$$settings{set}> is given, then a list 
  of possible remote directories will be listed instead.

  Args:   $$settings{set}:    The taxon to look at. If set,
                              the remote directories of this
                              taxon will be reported instead.
          $$settings{domain}: The ncbi ftp domain

=back

=cut
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

=over 12

=item C<listRemoteDirs>

  List the remote directories of the taxon given.
  Args:   $$settings{set}:    The taxon to look at. If set,
                              the remote directories of this
                              taxon will be reported instead.
          $$settings{domain}: The ncbi ftp domain

=back

=cut
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

=over 12

=item C<downloadAll>

  Downloads everything from the NCBI FTP site

  Args:   $$settings{set}:     The taxon to look at. If set,
                               the remote directories of this
                               taxon will be reported instead.
          $$settings{domain}:  The ncbi ftp domain
          $$settings{tempdir}: Where things are being
                               downloaded to
          $$settings{outdir}:  Where things will eventually be
                               moved to (but not by this sub)
          $$settings{maxTrees} Maximum number of trees to
                               download. Usually newer trees
                               are downloaded first.


=back

=cut 

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

=over 12

=item C<parseDate>
  
  Parses the date with C<Time::Piece>

  Arguments:   $date      string in either MM/DD/YYYY format or YYYY-MM-DD format
               $settings  (not used right now)

=back

=cut

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

=over 12

=item C<findBiosample>

  Finds the biosample information for a given query

  Args:  $biosample  string
         $settings   hash     (not used right now)

=back

=cut

sub findBiosample{
  my($query,$settings)=@_;
  
  
}

=over 12

=item C<colorScheme>

  Defines a set of colors
  Returns: array of RGB values, zero to one.
           Each RGB value is a three-element list.
  Palette was defined at http://paletton.com/#uid=62X0Z0k7xFM3Ur75UujeqG4j0KZ

=back

=cut

sub colorScheme{
  my($type,$settings)=@_;

  $type||="rgb0";

  # TODO allow for other types of color coding:
  #   hex, rgb(0-255), rgba(0-255 + transparency), etc
  # This currently only gives rgb0.
  
  my @color=(
    [0,0,0],              # This will be the default color
    [0.627,0.82,0.698],
    [0.541,0.616,0.569],
    [0.561,0.686,0.608],
    [0.451,0.82,0.588],
    [0.349,0.863,0.541],
    [0.631,0.651,0.792],
    [0.522,0.529,0.584],
    [0.549,0.561,0.651],
    [0.486,0.522,0.796],
    [0.408,0.459,0.843],
    [0.973,0.992,0.757],
    [0.827,0.835,0.733],
    [0.922,0.933,0.761],
    [0.957,0.992,0.545],
    [0.949,0.992,0.404],
    [1,0.808,0.765],
    [0.847,0.761,0.745],
    [0.945,0.804,0.773],
    [1,0.627,0.549],
    [1,0.51,0.408],
  );

  # Interleave the different colors, but leave the zeroth
  # color alone since it is for 'missing'.
  # Also, leave the first color alone since we don't really
  # need to shuffle it.
  @color = sort { $$a[1] <=> $$b[1] } @color;

  @_colorScheme=@color;
  for(@color){
    my $key=join(" ",@$_);
    $_colorScheme{$key}=1;
  }

  return @color if wantarray;
  return \@color;
}

sub randColor{
  my($settings)=@_;

  die "ERROR: need to call colorScheme() before this subroutine" if(!@_colorScheme);

  # Make up a new color
  my $color=[rand(1),rand(1),rand(1)];
  $_=sprintf("%0.3f",$_) for(@$color);

  while($_colorScheme{"@$color"}){
    logmsg "WARNING: I have already seen color @$color";
    $color=randColor($settings);
  }

  $_colorScheme{join(" ",@$color)}=1;
  
  return $color;
}


1;
