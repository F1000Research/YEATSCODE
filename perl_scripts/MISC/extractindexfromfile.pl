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
my ($forR,$reverse,$infile,$outfile,$which_tech,$seperator,$listfile);
my (@idx);
$seperator = ",";
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "infile=s"=>\$infile ,
            "reverse"=>\$reverse ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "forR=s"=>\$forR ,
            "seperator=s"=>\$seperator ,
            "idx=i"=>\@idx,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
die if(!@idx);

my $idx1 = $idx[0];
$outfile = "$infile.$idx1" if (! defined $outfile);
my $ofh = util_write($outfile);
my $CNT = 0 ; 
my ($RESULTDIR,$PDBDIR,$FASTADIR) = util_SetEnvVars();


print $ofh "$forR = ( " if(defined $forR);
my $first = 1 ;
while(<$ifh>){
     next if(/^\s*$/);
	 my (@l) = split " " ,$_ ;

	 if($first){
	 	$seperator = "";
	 }
	 else{
	 }

	 foreach my $i (@idx){
	     print $ofh "$l[$i]$seperator";
	 }
	 print $ofh "\n" if(!defined $forR);
}

print $ofh ")\n" if(defined $forR);
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
