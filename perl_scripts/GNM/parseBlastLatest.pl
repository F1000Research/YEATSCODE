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
my ($checkforsame,$infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions,$trs);
my $howmany = 100000 ;
my $verbose = 0 ;

my $percentlength = 10;
my $percentmatched = 70;
my $percentidentity = 30;
my $expectlimit = 0.00000000001;

GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "trs=s"=>\$trs ,
            "checkforsame"=>\$checkforsame ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "percentlength=i"=>\$percentlength ,
            "percentmatched=i"=>\$percentmatched ,
            "percentidentity=i"=>\$percentidentity ,
            "expectlimit=f"=>\$expectlimit ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_open_or_append($outfile);
my $ofhonlyone = util_open_or_append("onlyone");

print "-percentidentity $percentidentity -percentmatched $percentmatched -percentlength $percentlength\n" if($verbose);

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -trs ") if(!defined $trs);

## this is for the findgene iterative prog
my $ofhTRS = util_write("$trs.trs.out");

my ($info,$querylength) = util_PARSEBLAST($infile);

die "$trs - querylength not defined in $infile" if(!defined $querylength);


my $done = {};
my @others ;
my $found100percent = 0;
foreach my $k (@{$info}){
	my $org = $k ;
    my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;
	print "$name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect\n" if($verbose);
	$name =~ s/;//g;
	if($name eq $trs && $iden_percent eq 100){
		$found100percent  = 1 ;
		$querylength = $subjectlength;
		$done->{$name} = 1 ;
	}
	else{
		push @others, $k ;
	}
}
die if (!$found100percent && defined $checkforsame);
die if (!defined $querylength && defined $checkforsame);

my $percentdiffofL = ($querylength * $percentlength) /100 ;
my $maxL = $querylength + $percentdiffofL ;
my $minL = $querylength - $percentdiffofL ;

my $percentfMatched = ($querylength * $percentmatched) /100 ;

my $grouping = {};
my $foundone = 0 ;
my $cntfound = 0 ;
foreach my $k (@{$info}){
	my $org = $k ;
    my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;
	next if($expect > $expectlimit);
	#print  "EE $expect\n";

	$name =~ s/;//g;
	next if($trs eq $name);

	print "This one has $subjectlength $subjectmatched $iden_percent. MaxL = $maxL and minL = $minL, and querylength = $querylength\n" if($verbose);
	if($subjectlength > $minL && $subjectlength < $maxL){
		if($subjectmatched > $percentfMatched){
		    if($iden_percent > $percentidentity){
				 print "\t chosen \n" if($verbose);
		    	 $foundone = 1 ;
			     print $ofh "$trs $name 0\n" if(! exists $done->{$name});
			     $done->{$name} = 1;
				 $cntfound++;

				 $name =~ s/.ORF.*//;
				 print $ofhTRS "$name\n";
			}
		}
	}
}

print "Info: found $cntfound\n";
if(!$foundone){
	print $ofhonlyone "$trs $trs\n";
    print $ofh "$trs $trs 0 # only one \n";
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
