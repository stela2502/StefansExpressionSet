#! /usr/bin/perl -w

#  Copyright (C) 2015-11-18 Stefan Lang

#  This program is free software; you can redistribute it
#  and/or modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation;
#  either version 3 of the License, or (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <http://www.gnu.org/licenses/>.

=head1 CreateRandomForestScripts.pl

This script creates all required R scripts to be run on the server.

To get further help use 'CreateRandomForestScripts.pl -help' at the comman line.

=cut

use Getopt::Long;
use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';

use stefans_libs::root;

my (
	$help,  $debug,   $database, $infile,
	$trees, $forests, $splits,   $outfile,
);

Getopt::Long::GetOptions(
	"-infile=s"     => \$infile,
	"-trees=s"      => \$trees,
	"-forests=s"    => \$forests,
	"-splits=s"     => \$splits,
	"-outfile=s"    => \$outfile,

	"-help"  => \$help,
	"-debug" => \$debug
);

my $warn  = '';
my $error = '';

unless ( defined $infile ) {
	$error .= "the cmd line switch -infile is undefined!\n";
}

unless ( defined $trees ) {
	$warn .= "trees set to 433\n";
	$trees = 433;
}
unless ( defined $forests ) {
	$warn .= "forests set to 7\n";
	$forests = 7;
}
unless ( defined $splits ) {
	$warn .= "splits set to 32\n";
	$splits = 32;
}
unless ( defined $outfile ) {
	$error .= "the cmd line switch -outfile is undefined!\n";
}

if ($help) {
	print helpString();
	exit;
}

if ( $error =~ m/\w/ ) {
	print helpString($error);
	exit;
}

sub helpString {
	my $errorMessage = shift;
	$errorMessage = ' ' unless ( defined $errorMessage );
	return "
 $errorMessage
 command line switches for CreateRandomForestScripts.pl
 
 This script creates a set of R scripts, that can be run using the Sun Grid Engine to
 create the random forest based distance matrix needed from the random forest 
 unsupervised clustering.

   -infile       :the RData object containing the datRF data frame object to be analyzed

   -trees        :how many trees to use
   -forests      :how many forests to plant
   -splits       :how many scripts to run in paralele (one processor each)

   -outfile      :the outfile to be transfered to the server (tar.gz)

   -help         :print this help
   -debug        :verbose output
   

";
}

my ($task_description);

$task_description .= 'perl '
  . root->perl_include() . ' '
  . $plugin_path
  . '/CreateRandomForestScripts.pl';
$task_description .= " -infile $infile"         if ( defined $infile );
$task_description .= " -trees $trees"           if ( defined $trees );
$task_description .= " -forests $forests"       if ( defined $forests );
$task_description .= " -splits $splits"         if ( defined $splits );
$task_description .= " -outfile $outfile"       if ( defined $outfile );

open( LOG, ">$outfile.log" ) or die $!;
print LOG $task_description . "\n";
close(LOG);

my $fm = root->filemap($outfile);
mkdir( $fm->{'path'} ) unless ( -d $fm->{'path'} );
system("rm -Rf $fm->{'path'}tmp/") if ( -d "$fm->{'path'}tmp/" );
mkdir("$fm->{'path'}tmp/");

my ( @files, @gene_files, $ifm );

$ifm =  root->filemap($infile);

for ( my $i = 0 ; $i < $splits ; $i++ ) {
	open( RSCRIPT, ">$fm->{'path'}tmp/randomForest_worker_$i.R" )
	  or Carp::confess(
"I could not create the R script '$fm->{'path'}tmp/randomForest_worker_$i.R'\n$!\n"
	  );
	system("touch $fm->{'path'}tmp/randomForest_worker_$i.Rdata.lock");
	system("touch $fm->{'path'}tmp/randomForest_worker_genes_$i.Rdata.lock");
	print RSCRIPT "## this is using work published in\n"
	  . "## Tao Shi and Steve Horvath (2006) Unsupervised Learning with Random Forest Predictors. \n"
	  . "## Journal of Computational and Graphical Statistics. Volume 15, Number 1, March 2006, pp. 118-138(21)\n"
	  . "initial.options <- commandArgs(trailingOnly = FALSE)\n"
	  . "script.dir <- dirname ( initial.options[ grep( 'randomForest', initial.options ) ] )\n"
	  . "print ( paste( 'working directory = ',script.dir))\n"
	  . "setwd(script.dir)\n"
	  . "lock.name <- 'randomForest_worker_$i.Rdata'\n"
	  . "lock.name2 <- 'randomForest_worker_genes_$i.Rdata'\n"
	  . "source ('libs/Tool_RandomForest.R')\n"
	  . "load ('$ifm->{'filename'}')\n"
	  . "no.forests=$forests\n"
	  . "no.trees=$trees\n"
	  . "system.time (RF <- calculate_RF(  datRF, $forests , $trees ,imp=T, oob.prox1=T, mtry1=3 ))\n"
	  . "save_RF(RF, 'randomForest_worker_$i.Rdata' )\n"
	  ## The gene based randomForests
	  . "datRF <- data.frame(t( datRF ))\n"
	  . "attach(datRF)\n"
	  . "system.time (RF <- calculate_RF( datRF, $forests , $trees ,imp=T, oob.prox1=T, mtry1=3 ))\n"
	  . "save_RF(RF, 'randomForest_worker_genes_$i.Rdata' )\n"
	  . "release_lock (lock.name)\n"
	  . "release_lock (lock.name2)\n";
	close(RSCRIPT);
	$files[$i]      = "randomForest_worker_" . $i . ".Rdata";
	$gene_files[$i] = "randomForest_worker_genes_" . $i . ".Rdata";
}
## calculate the distance matrix
open( RSCRIPT, ">$fm->{'path'}tmp/randomForest_finisher.R" )
	  or Carp::confess(
"I could not create the R script '$fm->{'path'}tmp/randomForest_worker_finisher.R'\n$!\n"
	  );

print RSCRIPT "## this is using work published in\n"
	  . "## Tao Shi and Steve Horvath (2006) Unsupervised Learning with Random Forest Predictors. \n"
	  . "## Journal of Computational and Graphical Statistics. Volume 15, Number 1, March 2006, pp. 118-138(21)\n"
	  . "initial.options <- commandArgs(trailingOnly = FALSE)\n"
	  . "script.dir <- dirname ( initial.options[ grep( 'randomForest', initial.options ) ] )\n"
	  . "print ( paste( 'working directory = ',script.dir))\n"
	  . "setwd(script.dir)\n"
	  . "source ('libs/Tool_RandomForest.R')\n"
	  . "load ('$ifm->{'filename'}')\n"
  . "Rf.data <- read_RF ( c('"
  . join( "', '", @gene_files ) . "'),"
  . scalar(@gene_files) . "  )\n"
  . "datRF <- data.frame(t( datRF ))\n"
  . "distRF = RFdist(Rf.data, datRF, no.tree= $trees, imp=F)\n"
  . "save( distRF, file=\"RandomForestdistRFobject_genes.RData\" )\n"
  . "load ('$infile')\n"
  . "attach(datRF)\n"
  . "Rf.data <- read_RF ( c('"
  . join( "',\n '", @files ) . "'),"
  . scalar( @files )
  . " )\n"    ## 25*20 sec wait for finish of the worker scripts
  . "distRF = RFdist(Rf.data, datRF, no.tree= $trees, imp=F)\n"
  . "save( distRF, file=\"RandomForestdistRFobject.RData\" )\n";

close(RSCRIPT);

open ( SCRIPT ,">$fm->{'path'}tmp/single_proc.sh" ) or die $!;
print SCRIPT map { '/bin/bash -c \"DISPLAY=:7 R CMD BATCH --no-save --no-restore --no-readline -- '.$_."\n" } map{ $_ =~ s/Data$//; $_ } @files, "randomForest_finisher.RData";
close ( SCRIPT);


## now copy all important things to the folder
system( "cp $plugin_path/submit_RF_to_psub.pl $fm->{'path'}tmp/" );
system( "cp $infile $fm->{'path'}/tmp/$ifm->{'filename'}");
mkdir( "$fm->{'path'}tmp/libs");
system( "cp $plugin_path/../R/Tool_RandomForest.R fm->{'path'}tmp/libs");
system( "cp $plugin_path/submit_RF_to_psum.pl $fm->{'path'}tmp/" );
chdir( "$fm->{'path'}tmp/" );
system( "tar -cf $fm->{'filename_core'}.tar *" );
system( "gzip2 -9 $fm->{'filename_core'}.tar");
system( "cp $fm->{'filename_core'}.tar.gz ../" );
unless ( $debug ){
	chdir( "$fm->{'path'}" );
	system( "rm -Rf $fm->{'path'}/tmp" );
}
print "please copy the file $fm->{'path'}$fm->{'filename_core'}.tar.gz to the calculation server extract all files and run 'perl submit_RF_to_psub.pl -email <your email address> -files randomForest_*.R' in the folder\n";


