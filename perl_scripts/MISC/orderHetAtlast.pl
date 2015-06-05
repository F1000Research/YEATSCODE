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
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);



my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();

my @list= util_read_list_sentences($listfile);
my $list = {};
map { s/\s*//g ; $list->{$_} = 1 ; } @list ;

foreach my $i (@list){
    my $pdb = "$PDBDIR/$i.pdb";
	my $ifh = util_read($pdb);
	my $ofh = util_write("temp.$i");

	my @het ;
    while(<$ifh>){
	    if(/^HET/){
	       push @het, $_ ;	
	    }
		else{
			print $ofh "$_";
		}
    }
    close($ifh);

	foreach my $h (@het){
		print $ofh $h ;
	}
	close($ofh);
	system("mv -f temp.$i $pdb");
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
