#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use PDB;
use ConfigPDB;

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($infile,$outfile,$which_tech,$listfile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
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

my @list= util_read_list_sentences($listfile);
my $list = {};
map { s/\s*//g ; $list->{$_} = 1 ; } @list ;

my $DB = {};
my $DBDONE = {};
my $done = {};
my $indices = {};
foreach my $i (@list){
   my $ifh = util_read($i);
   print "Merging $i\n";
   while(<$ifh>){
     next if(/^\s*$/);
     next if(/^\s*#/);
	 my @l  = split ; 
	 my $nm = shift @l ;
	 foreach my $index (@l){
		    if(!exists $DB->{$nm}){
		        $DB->{$nm} = [];	
		    }
		    if(! exists $DBDONE->{$nm}->{$index}){
			    $DBDONE->{$nm}->{$index} = 1 ;
		        push @{$DB->{$nm}}, $index; 
		    }
	      }
	   }
   close($ifh);
}

my @L = keys %{$indices};
my $N = @L ;
die "There were $N indices\n";


foreach my $k (keys %{$DB}){
	 my @list = @{$DB->{$k}} ; 
	 print $ofh "$k @list \n";
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
