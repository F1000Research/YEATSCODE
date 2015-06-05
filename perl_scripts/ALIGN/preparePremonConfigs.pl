#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use PDB;
use ConfigPDB;
use MyGeom;
# just test

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($infile,$outfile,$anndir,$polar,$listfile,$protein);
my (@expressions,$config);
my $size ;
my $verbose = 1 ;
my $LARGESTDISTALLOWED = 30 ;
GetOptions(
            "anndir=s"=>\$anndir ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "config=s"=>\$config,
            "size=i"=>\$size ,
            "polar=s"=>\$polar ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
#usage( "Need to give a howmany -option -howmany  ") if(!defined $howmany);
usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
usage( "Need to give a size pdb id -option -size  ") if(!defined $size);
usage( "Need to give a polar pdb id -option -polar  ") if(!defined $polar);
usage( "Need to give a anndir pdb id -option -anndir  ") if(!defined $anndir);
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $PWD = cwd;

my  ($seconds, $microseconds) = gettimeofday;

my $pdb = "$PDBDIR/$protein.pdb";
my $pdb1 = new PDB();
$pdb1->ReadPDB($pdb);

my @list= util_read_list_sentences($listfile);
my $list = {};
my @reallist ; 
#map {my @j = split ; push @reallist, @j;  } @list ;
map { s/,/ /g ; my @j = split ; push @reallist, @j;  } @list ;

my $sorting = {};
foreach my $i (@reallist){
	my ($nm,$type) = split "/", $i ;
	my ($name,$number) = ($nm =~ /([a-zA-Z]+)([0-9]+)/);
	my $x = $pdb1->GetSingleLetter($name);
	#print "$x $number\n";
	$x = $x.$number;
	$sorting->{$x} = $i ; 
}
my @tmplist ;

my @backboneatoms = qw( O N CA);
my $backboneatomstable = util_make_table(\@backboneatoms);

my $IDX = 0 ;
foreach my $k (sort keys %{$sorting}){
	my $i = $sorting->{$k};
	my ($nm,$type) = split "/", $i ;
	$IDX++;
	my ($name,$number) = ($nm =~ /([a-zA-Z]+)([0-9]+)/);
	my $x = $pdb1->GetSingleLetter($name);
	#print "$x $number\n";
	my $newnm = $name . $number . $type ; 
	if(exists $backboneatomstable->{$type}){
	     system( "cp $SRC/PREMONCONFIGS/aalist.allinone $anndir/$protein.$size.$polar.premon.in.aalist$IDX\n");
	     system( "cp $SRC/PREMONCONFIGS/config.grp.all$type $anndir/$protein.$size.$polar.premon.in.config$IDX\n");
	}
	else{
	     system( "cp $SRC/PREMONCONFIGS/aalist.normal $anndir/$protein.$size.$polar.premon.in.aalist$IDX\n");
	     system( "cp $SRC/PREMONCONFIGS/config.grp $anndir/$protein.$size.$polar.premon.in.config$IDX\n");
	}
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
