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
my ($pairwise,$infile,$outfile,$trs,$which_tech,$listfile,$protein);
my (@expressions);
my $checklength ;
my $verbose = 1 ;
my $cutoff = 30 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "pairwise"=>\$pairwise ,
            "outfile=s"=>\$outfile ,
            "trs=s"=>\$trs ,
            "expr=s"=>\@expressions,
            "checklength=i"=>\$checklength ,
            "cutoff=i"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_append($outfile);
my $ofhTRS = util_write("$trs.trs.out");

system ("touch logclose") if(! -e "logclose");
my $ofhlog = util_append("logclose");

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;


my ($fulllength,$TRSnm,$STRS,$LLLs,$SCORES,$identities) = util_ParseWebBlast($infile);

my @TRSnm = @{$TRSnm};
my @STRS = @{$STRS};
my @LLLs = @{$LLLs};
my @identities = @{$identities};
   
my $NN = @STRS -1 ;
#warn "Not same for $infile : $N and $NN..." if($N ne $NN);

print $ofh "$infile\n";
my $wrote = 0 ;
#:wqprint "$NN = nmnmm\n";

foreach my $idx (0..$NN){
	#my $a = $l[$idx];
	my $str = $STRS[$idx];
	my $TRS = $TRSnm[$idx];
	$TRS =~ s/.ORF.*//;
	my $length = $LLLs[$idx];

	my @l = split " ", $str ;
	#my $b = @l[@l -1] ;
	my $b = $identities[$idx];
	#print "$b = cutoff \n";
	if($b > $cutoff){

	    my $diff = abs($fulllength - $length) ;
	    #print "$TRS $b $diff = $fulllength -$length lll\n";
		next if(defined $checklength && $diff >   $checklength);
		#next if(defined $checklength && ($fulllength <  460 || $fulllength > 560));
	    $wrote++;
	    print $ofh "$str $diff  \n";
		print $ofhTRS "$TRS\n";
	}
}
print $ofh "--\n";


print "cutoff = $cutoff, wrote $wrote \n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
