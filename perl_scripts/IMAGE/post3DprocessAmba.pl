#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use ConfigPDB;
use MyGeom;
use PDB;
use MyGeom;
use MyMagick;
use MyUtils;
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use ConfigPDB;
use PDB;
use POSIX ;
use Algorithm::Combinatorics qw(combinations) ;
use Math::Geometry ;
use Math::Geometry::Planar;

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($infile,$outfile,$which_tech,$listfile,$resolution,$protein,$randomsampling);
my (@expressions);
my $threshold ;
my $verbose = 1 ;
my $ignoredist ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "threshold=i"=>\$threshold ,
            "randomsampling=i"=>\$randomsampling ,
            "ignoredist=i"=>\$ignoredist ,
            "resolution=f"=>\$resolution ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofhmin = util_write("mindist");
my $ofhminzero = util_write("mindistzero");
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a resolution -option -resolution  ") if(!defined $resolution);
usage( "Need to give a threshold -option -threshold  ") if(!defined $threshold);
usage( "Need to give a ignoredist -option -ignoredist  ") if(!defined $ignoredist);
usage( "Need to give a randomsampling -option -randomsampling  ") if(!defined $randomsampling);
my $CNT = 0 ; 
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();

my $PWD = cwd;

my @list= util_read_list_sentences($listfile);
my $list = {};
map { s/\s*//g ; $list->{$_} = 1 ; } @list ;
my $TOTAL = @list ;

my $info = {};
my $cnt = 0 ; 
my $CCC = 0 ;
my $CNT = 0 ; 
foreach my $i (@list){

	if($randomsampling){
	    if(util_round(rand())){
			next ;
		}

	}

    my $ifh = util_read($i);

	my $Ngroups = 0 ;
	my $mindist = 1000 ;
	my $mindistNozero = 1000 ;
	my @groups ; 
    my $groupinfo = {};
    while(<$ifh>){


		if(/GROUP/){
			my @l = split ; 
			shift @l ;
			my $N = @l ; 


			$Ngroups++;
			$groupinfo->{$Ngroups} = \@l ; 
			push @groups, \@l ;

		}
		else{
			my ($p,$dist) = split ; 
			$info->{$p} = $dist ; 
			$mindist = $dist if($mindist > $dist);
			$mindistNozero = $dist if($mindistNozero > $dist && $dist > 0 );
		}


	
    }

	my $validgroups = 0 ; 
	foreach my $k (sort {$a <=> $b}  keys %{$groupinfo}){
		my $l = $groupinfo->{$k} ; 
		my @l = @{$l};
		my $N = @l;
		next if($N < $threshold);
		foreach my $p (@l){
			  my $dist = $info->{$p} ;
			  ################# ignoring 0 dist#####
			  next if($ignoredist && !$dist);
			  $dist = $dist * $resolution ;
			  print $ofhmin "$CNT $dist\n";
			  $CNT++ ; 
		}
		$validgroups++;
	}
	print "validgroup for $i = $validgroups\n";
	if($validgroups){
		print "$i is good !! \n";
	}

	if($Ngroups > 0){
		$CCC++;
	}

	$mindistNozero = $mindistNozero * $resolution;
	$mindist= $mindist* $resolution;
	if($Ngroups eq 1 || $Ngroups eq 2){
		$cnt++ ;

		#print $ofhmin "$cnt $mindist\n";
		#print $ofhminzero "$cnt  $mindistNozero\n";

		#foreach my $g (@groups){
			#my @l = @{$g};
			#my ($x, $y, $z) ;
			#$x = $y = $z = 0 ; 
			#foreach my $p (@l){
				#print "P = $p \n";
	            #my ($X,$Y,$Z) = geom_MakeCoordFromKey($p);
		        #$Z = $RATIO * $Z ; 
				#$x = $x + $X ; X
				#$y = $y + $Y ; 
				#$z = $z + $Z ; 
			#}
		#}


	}
    close($ifh);
}
#print "Processed $cnt images out of $TOTAL\n";
print "found something in $CCC of $TOTAL \n";


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
