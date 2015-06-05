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




use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($helix,$ann,$config,$p1,$p2,$infile,$outfile,$which_tech,$listfile,$protein);
my $maxdist ;
my $DISTANCEWITHOUTSEQMATCH = 1 ;
my $verbose = 1 ;
my $isnew = 0 ;

my ($verify,$radii,$before1,$before2);
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "verify"=>\$verify ,
            "p1=s"=>\$p1 ,
            "helix=s"=>\$helix ,
            "p2=s"=>\$p2 ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "ann=s"=>\$ann ,
            "maxdist=f"=>\$maxdist ,
            "config=s"=>\$config,
            "radii=i"=>\$radii ,
            "isnew=i"=>\$isnew ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a config file name => option -config ") if(!defined $config);
usage( "Need to give a radii file name => option -radii ") if(!defined $radii);
usage( "Need to give a radii file name => option -protein ") if(!defined $protein);


my $ofh = util_write($outfile);
my $ofhstats = util_append("statsCLASP.AH");
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $PWD = cwd;

ConfigPDB_Init($config);


my $info = {};
my @resultlines ;
my $ifh = util_read($infile);
my @proteins ;
push @proteins, $protein; 


my @info = util_ReadPdbs($PDBDIR,$APBSDIR,1,@proteins) ; 
my $info = shift @info ;
my $pdb1 = $info->{PDBOBJ};
my $pqr1 = $info->{PQR};
my $pots1 = $info->{POTS};


$, = " ";

while(<$ifh>){
     next if(/^\s*$/);
     next if(/RESULT/);
     chop ;

	 my @l = split ;
	 while(@l>1){

		 my $a = shift @l ;
		 my $b = $l[0];

		 my $STR = " $a $b ";
	     my ($dist,$pots,$name) = util_ProcessSingleLine($pdb1,$pqr1,$pots1,$STR,$isnew);
	     my @dist = @{$dist}; 
	     my @pots = @{$pots}; 
    
	     print $ofh "$name, D " , @dist, "\n";
	     print $ofh " PD ", @pots, "\n";

		 my @l = split ",", $name ;

		 my $FINALNM = "";
		 my $REALNM = "";
		 foreach my $l (@l){
		     my @l1 = split "/", $l ;
			 $REALNM = $REALNM . $l1[0];
			 my ($NAME) = ($l1[0] =~ /(...)/);
			 print "$l1[0] $NAME ,,,,,,,,,,,,,,,,,,,\n";;
			 $FINALNM = $FINALNM . $NAME ; 
		 }

		 print $ofhstats "$protein $FINALNM @dist @pots $REALNM $helix \n";
	 }

}


print "Wriitng to $outfile\n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
