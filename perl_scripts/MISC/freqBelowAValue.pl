#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use ConfigPDB;
use MyGeom;

use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($reverse,$infile,$outfile,$which_tech,$listfile);
my (@expressions);
my $idx = 0;
my $cutoff = 0;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "infile=s"=>\$infile ,
            "reverse"=>\$reverse ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "idx=i"=>\$idx ,
            "cutoff=i"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -cutoff ") if(!defined $cutoff);
my $ifh = util_read($infile);
my $ofh = util_write($outfile);
my $CNT = 0 ; 
my ($RESULTDIR,$PDBDIR,$FASTADIR) = util_SetEnvVars();


my $lessoreq = 0 ;
my $greater = 0 ;
while(<$ifh>){
     next if(/^\s*$/);
	 my (@l) = split ;
	 my $v = $l[$idx];
	 if($v <= $cutoff){
	 	$lessoreq++;
	 }
	 else{
	 	$greater++;
	 }
}

my $tot = $lessoreq + $greater;
my $perL = 100* $lessoreq/$tot;
my $perG = 100* $greater/$tot;
print "% lesser = $perL, % greater $perG\n";


close($ifh);


print STDERR "Output written in $outfile\n";

chmod 0777, $outfile ;
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
