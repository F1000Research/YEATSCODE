#!/usr/bin/perl -w 
use strict ;
use PDB;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use ConfigPDB;
use MyGeom;
use MyPymol;
use Math::Geometry ;
use Math::VectorReal qw(:all);  # Include O X Y Z axis constant vectors

use AAConfig;

my $aaconfig = new AAConfig("/home/sandeepc/aa.config");



my $padding = 3 ; 
my $DISTSURFACE = 1 ; 
my $DOSURFACE = 0 ; 
my $DISTATOMS   = 1 ; 
my $BOUNDARYDELTA   = 1000 ; 
my $RECALIBRATECOUNT= 100 ; 
my $CNT1000 = 3 ;
my $ADDATOMSCLOSE2 = 0 ;

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($debugresidue,$ann,$config,$p1,$p2,$infile,$threshPD,$threshsign,$threshDist,$outfile,$readpotential,$which_tech,$listfile,$protein);
my $DISTANCEWITHOUTSEQMATCH = 1 ;
my $verbose = 1 ;

my ($onlypolar,$radii,$before1,$before2);
$readpotential = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "debugresidue=s"=>\$debugresidue ,
            "onlypolar=i"=>\$onlypolar ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "ann=s"=>\$ann ,
            "config=s"=>\$config,
            "radii=i"=>\$radii ,
            "padding=i"=>\$padding ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a config file name => option -config ") if(!defined $config);
usage( "Need to give a p1 file name => option -p1 ") if(!defined $p1);



my $ofh = util_write($outfile);
my $fhatoms = util_write("coords.atoms");
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $PWD = cwd;

#ConfigPDB_Init($config,$ofh);

my $pdb = "$PDBDIR/$p1.pdb";
my $pdb1 = new PDB();
$pdb1->ReadPDB($pdb);

my @atoms1 = $pdb1->GetAtoms();
my ($minx,$miny,$minz);
my ($maxx,$maxy,$maxz);
$minx = $miny = $minz = 100000 ;
$maxx = $maxy = $maxz = -10000 ;
foreach my $a1 (@atoms1){
	my ($x,$y,$z) = $a1->Coords();
	next if(!defined $x);
	$maxx = $x if($x > $maxx); $maxy = $y if($y > $maxy); $maxz = $z if($z > $maxz);
	$minx = $x if($x < $minx); $miny = $y if($y < $miny); $minz = $z if($z < $minz);
}

### round off 
$maxx = int($maxx); $maxy = int($maxy); $maxz = int($maxz);
$minx = int($minx); $miny = int($miny); $minz = int($minz);

### pad 
$maxx = $maxx + $padding ; $maxy = $maxy + $padding ; $maxz = $maxz + $padding ;
$minx = $minx - $padding ; $miny = $miny - $padding ; $minz = $minz - $padding ;

print "padding = $padding \n";
print "maxx maxy maxz\n";
print "$maxx $maxy $maxz\n";
print "minx miny minz\n";
print "$minx $miny $minz\n";
my $diffx = $maxx - $minx ; my $diffy = $maxy - $miny  ; my $diffz = $maxz - $minz  ;
my $CNT = 0 ; 
my $pseudoatom = new Atom();
$pseudoatom->SetIdx(10000);

my $natoms = @atoms1 ;
my $CNT = 0 ;
my $DIST = 5 ;

my $ignore = {};
foreach my $a1 (@atoms1){
	my ($x,$y,$z) = $a1->Coords();

	my $X = int($x);
	my $Y = int($y);
	my $Z = int($z);
	my $str = geom_MakeKeyFromCoord($X,$Y,$Z);
	if(! exists $ignore->{$str}){
	    $ignore->{$str} = 1;
	    $CNT++;
	}

	my @l = MakeCube($str,$DIST -2);
	foreach my $k (@l){
		if(! exists $ignore->{$k}){
		   $ignore->{$k} = 1;
		   $CNT++;
		}
	}
}

my $cnt = 0 ;
my @res ; 
my $ALL = {};
my @ALL = ();
foreach my $x ($minx..$maxx){
    foreach my $y ($miny..$maxy){
        foreach my $z ($minz..$maxz){
	        my $str = geom_MakeKeyFromCoord($x,$y,$z);
			$ALL->{$str} = 1 ;
			push @ALL, $str ;
        }
    }
}

my $NALL = @ALL ;
print "There are $NALL ignore= $CNT \n";
my $ccc = 0 ;
foreach my $str (@ALL){
			next if(exists $ignore->{$str});
			my @l = MakeCube($str,$DIST) ;
            my $isempty =  CheckIfEmpty(\@l,$ignore);
			if($isempty){
			    print "$str is empty\n";
			}
			$ccc++ ;
			#print "$ccc\n";
}


foreach my $r (@res){
	print $ofh $r, "\n"; 
}
my $gridresilution = 1 ; 




sub usage{
my ($msg) = @_ ;
print $msg , "\n" ;
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
die ;
}


sub MakeCube{
  my ($k,$D) = @_ ;
  my ($x,$y,$z) = geom_MakeCoordFromKey($k);

  my $minx = $x - $D ;
  my $maxx = $x + $D ;
  my $miny = $y - $D ;
  my $maxy = $y + $D ;
  my $minz = $z - $D ;
  my $maxz = $z + $D ;

  my @ret ; 
  foreach my $x ($minx..$maxx){
    foreach my $y ($miny..$maxy){
        foreach my $z ($minz..$maxz){
	        my $str = geom_MakeKeyFromCoord($x,$y,$z);
			push @ret, $str ;
		}
	}
  }
  return @ret ;
}


sub CheckIfEmpty{
	my ($l,$filled) = @_;

	my @l = @{$l}; 
	foreach my $l (@l){
		return 0 if(!exists $ALL->{$l});
		return 0 if(exists $filled->{$l});
	}
	return 1 ;
	
}
