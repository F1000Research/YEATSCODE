#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use PDB;
use ConfigPDB;
use Math::Geometry ;
use Math::Geometry::Planar;


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($duplicate,$blastout,$infile,$p1,$p2,$outfile,$trs,$cutoff,$which_tech,$listfile,$protein);
my ($chooseone,$ignorefile,$tag,$donone,@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "trs=s"=>\$trs ,
            "tag=s"=>\$tag ,
            "blastout=s"=>\$blastout ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "donone"=>\$donone ,
            "duplicate"=>\$duplicate ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a infile -option -infile  ") if(!defined $infile);
usage( "Need to give a tag -option -tag  ") if(!defined $tag);
$outfile = "$tag.$infile";
my $ofh = util_write($outfile);
my $ifh = util_read($infile);

my $table = {};
{
my @list= util_read_list_sentences($listfile);
map {  my @ll = split ;  $table->{$ll[0]} = 1 ; } @list ;
}

my $cnt = 0 ;
while(<$ifh>){

   my $orog = $_;

   # just hardcoded, essentially we take the first
   s/ln -s//;
   s/.ORF.*//;
   my ($trs) = split ;
   if(exists $table->{$trs}){
   	   $cnt++;
   	   print $ofh "$orog";
   }
}

print "Wrote $cnt to $outfile\n";
   
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
