#!/usr/bin/env perl
package NcbiClusters;
use strict;
use warnings;

use Net::FTP;

use Exporter qw/import/;

our @EXPORT_OK = qw(logmsg listSets);

sub logmsg { print STDERR "$0: @_\n"; }

=encoding UTF-8
=head1 NcbiClusters

A module for reporting on NCBI SNP clusters

=head2 Methods

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


1;
