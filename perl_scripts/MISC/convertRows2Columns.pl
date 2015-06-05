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
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;

my $cnt = 0 ; 
my $max = 0 ;
my @L ; 
while(<$ifh>){
     next if(/^\s*$/);
	 my (@l) = split ; 
	 push @L,\@l ;
	 my $N = @l ;
	 $max = $N if($N > $max);
}
print "Max = $max \n";

my @LL ; 
my $CNT = 0 ; 
foreach my $i (@L){
   my @l = @{$i};
   my $N = @l ;
   my $diff = $max - $N ; 
   print "Adding $diff \n";
   foreach my $x (1..$diff){
   	  push @l, "";
   }
   $CNT++;
   my @TMPL = @l ;
   my $WHAT = shift @TMPL ;

   my $idx = 1 ;
   my $nm = "$WHAT.data.$CNT";
   #my $OFH = util_write($nm);
   #foreach my $i (@TMPL){
        #print $OFH "$idx $i \n";
		#$idx++;
   #}
   push @LL,\@l;
   
}

my @strs ; 
foreach my $i (1..$max){
	my $s = "";
	push @strs, $s ; 
}

while(@LL){
	my $l = shift @LL;
	my @l = @{$l};
	my $N = @l -1 ; 
	foreach my $idx (0..$N){
		$strs[$idx] =  $strs[$idx]  . $l[$idx] . "\t";
	}
}


print "Writing to $outfile\n";
foreach my $s (@strs){
	print $ofh "$s \n";
}

close($ifh);

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
