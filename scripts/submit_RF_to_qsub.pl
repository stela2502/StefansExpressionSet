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

=head1 submit_RF_to_qsub.pl

This script submits all scripts to the qsub, waits for them to finish and sums up the data in the end.

To get further help use 'submit_RF_to_qsub.pl -help' at the comman line.

=cut

use Getopt::Long;
use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

use Cwd;
my $wd = getcwd;

my $VERSION = 'v1.0';


my ( $help, $debug, $database, @files, $email);

Getopt::Long::GetOptions(
	 "-files=s{,}"    => \@files,
	 "-email=s" => \$email,
	 

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $files[0]) {
	$error .= "the cmd line switch -files is undefined!\n";
}
unless ( defined $email ) {
	$error .= "the cmd line switch -email is undefined!\n";
}

if ( $help ){
	print helpString( ) ;
	exit;
}

if ( $error =~ m/\w/ ){
	print helpString($error ) ;
	exit;
}

sub helpString {
	my $errorMessage = shift;
	$errorMessage = ' ' unless ( defined $errorMessage); 
 	return "
 $errorMessage
 command line switches for submit_RF_to_qsub.pl

   -files   :all R files you want to getrun on SGE
   -email   :your mail for the SGE

   -help           :print this help
   -debug          :verbose output
   

"; 
}


my ( $task_description);

$task_description .= 'perl '.$plugin_path .'/submit_RF_to_qsub.pl';
$task_description .= ' -files '.join( ' ', @files ) if ( defined $files[0]);

for ( my $i = 0; $i < @files; $i ++ ){ ## the last script is the sum up script and that has to be run after all others have finished!
	open ( QSCRIPT, ">$files[$i].sh" ) or die $!;
	print QSCRIPT 
	'#!/bin/bash'."\n"
	.'#$ -S /bin/bash'."\n"
	.'#$ -M '.$email."\n"
	.'#$ -m eas'."\n"
	.'#$ -pe orte 1'."\n"
	.'#$ -l mem_free=20G'."\n"
	.'cd '.$wd."\n"
	.'R CMD BATCH --no-save --no-restore --no-readline -- '.$files[$i]."\n";
	close ( QSCRIPT);
}
for ( my $i = 0; $i < @files-1; $i ++ ){
	system( "qsub $files[$i].sh");
}
while ( 1 ) {
	open ( CHECK, "ls *.lock|" );
	if ( split( ' ', <CHECK>) == 0 ){
		close ( CHECK);
		last;
	}
	close(CHECK);
}
system( "qsub $files[$#files].sh");

