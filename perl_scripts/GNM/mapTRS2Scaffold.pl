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
my ($infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();

my $info = {};
my $mapTRS = {};
while(<$ifh>){
	$CNT++;
	 my ($trs,$Subjectname,$querylength,$qS,$qE,$sS,$sE,$diffQ,$diffS,$diffMatchedQ,$diffMatchedQ2S,$ignoredfar) = split ;
	 if(!defined $info->{$Subjectname}){
	 	$info->{$Subjectname} = [];
	 	$mapTRS->{$Subjectname} = {};
	 }
	 die "$sS $sE" if ($sS > $sE);

	 my $done = {};
	 foreach my $i ($sS..$sE){
	 	if(exists $mapTRS->{$Subjectname}->{$i}){
			my $oldtrs = $mapTRS->{$Subjectname}->{$i};
			#print "found two $trs and $oldtrs mapping to $Subjectname at $i \n" if(! exists $done->{$oldtrs});
			print $ofh "$trs $oldtrs\n";
			$done->{$oldtrs} = 1 ;
		}
		$mapTRS->{$Subjectname}->{$i} = $trs ;
	 }
	 push @{$info->{$Subjectname}} , $sS;
	 push @{$info->{$Subjectname}} , $sE;

}

#foreach my $k (keys %{$

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
