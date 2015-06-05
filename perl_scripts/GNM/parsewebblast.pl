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
my ($pairwise,$infile,$outfile,$which_tech,$listfile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
my $cutoff = 0.00005 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "pairwise"=>\$pairwise ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_append($outfile);

system ("touch logclose") if(! -e "logclose");
my $ofhlog = util_append("logclose");

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;

my $info = {};
my ($junk,$percent);
my @l ; 
my $push = 0 ; 
my $length  ;
my $Identities  ;
while(<$ifh>){
     next if(/^\s*$/);
	 chomp ;
     if(/Length/ && ! defined $length){
	    ($length) = (/Length=(.*)/)	;
	 }
	 elsif(/Identities/ && ! defined $Identities){
	 	($Identities,$junk,$percent) = (/Identities\s*=\s*(\d+)\/(\d+)\s*\((\d+)\%\)/);
		 last ;
	 }
	 elsif(/Sequences producing significant alignments/){
		$push = 1 ; 
		next ;
	 }
	 if($push && !/(chromosome|hypothe|unnamed|uncharacterized)/i){
	     push @l, $_ ;
	 }
}

my $percentdiff = 100 - $percent ;

($infile =~ s/.*BLASTOUT..//);

my $v = "";
if(!defined $pairwise){
    $v = $l[0] ;
    if(!defined $v || $v =~ /Length/){
	    print "Warning:Did not find any characterized\n";
	    print $ofh "$infile Warning:Did not find any characterized\n";
	    exit;
    }
}

my $diff = $length - $Identities ;
my @lll = split " ", $v ;
my $N = @lll ;
my $evalue = defined $pairwise? 0 : $lll[$N-1];
if($evalue > $cutoff){
	print "Warning; $infile has evalue $evalue more than cutoff $cutoff\n";
}

$v =~ s/            //g;
print $ofh "$infile $length $Identities/$junk $diff $percentdiff  $v\n";
print  "$v\n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
