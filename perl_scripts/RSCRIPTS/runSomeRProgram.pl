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
my ($varname,$absolute,$reverse,$infile,$cutoff,$outfile,$which_tech,$listfile);
my (@expressions);
my $idx = 0;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "infile=s"=>\$infile ,
            "varname=s"=>\$varname ,
            "absolute"=>\$absolute ,
            "reverse"=>\$reverse ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "idx=i"=>\$idx ,
            "cutoff=i"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -idx ") if(!defined $idx);
usage( "Need to give a input file name => option -varname ") if(!defined $varname);
my $CNT = 0 ; 
my ($RESULTDIR,$PDBDIR,$FASTADIR) = util_SetEnvVars();

my $ofh = util_open_or_append($outfile);


my ($list,$mean,$sd,$assign) = util_ExtractOneIdxFromFile($infile,$idx,$varname,100);



print $ofh " #$infile\n";
print $ofh " $assign\n";
print $ofh "wilcox.test(A,B,paired=TRUE)\n";
#print $ofh " t.test($varname,mu=$mean) \n";
#print $ofh " t.test($varname) \n";
#print $ofh " shapiro.test($varname)\n";

#print "wilcox.test(C, Y, paired = T)\n";

chmod 0777, $outfile ;
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
