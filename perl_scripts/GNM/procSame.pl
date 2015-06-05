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
my ($infile,$p1,$p2,$outfile,$which_tech,$listfile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
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

my $info = {};
my $CNT = 0 ;
my $cnt = 0 ;
while(<$ifh>){
	next if(/Warning/i);
	next if(/ORF/i);
	$cnt++;
	 my ($nm,$two,$three,$diff,$percentdiff,$junk) = split ; 
	 #next if($percentdiff > 2);
	 $CNT++;
	 $info->{$junk} = {} if(!defined $info->{$junk});
	 $info->{$junk}->{$nm} = 1 ;
}

print "Process $CNT out of $cnt\n";

my $sortn = {};
foreach my $k (keys %{$info}){
	 my $val = $info->{$k} ; 
	 my $N = (keys %{$val});
	 $sortn->{$k} = $N ;
}

my @s = sort { $sortn->{$b} <=> $sortn->{$a}} (keys %{$sortn});



my $revtrs = {};
my $postrs = {};
my $ifhREVTRS = util_read("posrev.TRS");
while(<$ifhREVTRS>){
	 my ($nm,$what,$val) = split ; 
	 if($what){
	     $postrs->{$nm} = $val ;
	 }
	 else{
	     $revtrs->{$nm} = $val ;
	 }
}

my $DIFFDIR = 0 ;
my $ofhdiffdir = util_write("diffdir");
foreach my $s (@s){
	my $v = $sortn->{$s};
	print $ofh  "$s $v \n";

	my $val = $info->{$s};

	my $dir ;
	my $prev ;

	my @POS ; 
	my @REV ; 
	foreach my $i (keys %{$val}){
		#if ($s eq "super521"){
			#print "$i lll\n";
		#}
		#next ;
	     my ($subject,$isrev) = util_Blastout_isrev("BLASTOUT_2WGS/$i.blast.nt");
		 die "$subject $s BLASTOUT_2WGS/$i.blast.nt" if($subject ne $s);


         #if($isrev && exists $revtrs->{$i}){
		 	if(exists $postrs->{$i} && $postrs->{$i} eq 1){
		 		$isrev = 1 ;		
			}
		 	elsif(exists $revtrs->{$i} && $revtrs->{$i} eq 1){
		 		$isrev = 0 ;		
			}

		 if(1){
		 $dir = $isrev if (!defined $dir);
		 if($dir ne $isrev){
		    print $ofhdiffdir "$s has different dirs.. for $prev ($dir) and $i ($isrev)\n" ;
			$DIFFDIR++;
		    #last ;
		 }
		 $prev = $i ;
		 }

	}

	my $NP = @POS ;
	my $NR = @REV ;
	my $M = $NP * $NR ;
	if($M){
	   print $ofhdiffdir "$s has pos $NP and neg $NR\n";
	}

}
print "There were $DIFFDIR with diff dirs, output in diffdir\n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
