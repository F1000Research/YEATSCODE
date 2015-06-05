#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use ConfigPDB;
use MyGeom;
use PDB;

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($chainID,$infile,$outfile,$onecolor,$color,$listfile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "color=s"=>\$color ,
            "chainID=s"=>\$chainID ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a chainID pdb id -option -chainID  ") if(!defined $chainID);
usage( "Need to give a color pdb id -option -color  ") if(!defined $color);
my $CNT = 0 ; 
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $PWD = cwd;

my @colors = qw ( magenta cyan yellow blue black red green gray orange  ) ;


my $N = @colors ; 

## lsit contains list of 2 numbers - start and end
my @list= util_read_list_words($listfile);
my $list = {};
map { s/\s*//g ; $list->{$_} = 1 ; } @list ;

my $cnt = 0 ; 
while(@list){

   my $a = shift @list ;
   print $ofh "select block_query$cnt, /PDBB//$chainID/$a\n";
   print $ofh "color $color, block_query$cnt\n";
   print $ofh "show sticks, block_query$cnt \n";

}


chmod 0777, $outfile ;

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
