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
my ($holo2apo,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "holo2apo=s"=>\$holo2apo ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);

system ("touch logclose") if(! -e "logclose");
my $ofhlog = util_append("logclose");

usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
#usage( "Need to give a holo2apo -option -holo2apo  ") if(!defined $holo2apo);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;


my @NONPOLAR = qw(PRO ALA GLY ILE LEU MET VAL);
my $NONPOLAR = {};
foreach my $i (@NONPOLAR){
	$NONPOLAR->{$i}= 1;
}

my @backboneatoms = qw( O N CA);
my $backboneatomstable = util_make_table(\@backboneatoms);

my $Allowed = {};
$Allowed->{"ON"} = 1 ; $Allowed->{"NO"} = 1 ;
$Allowed->{"OH"} = 1 ; $Allowed->{"HO"} = 1 ;
$Allowed->{"OO"} = 1 ;
$Allowed->{"FN"} = 1 ;
$Allowed->{"PO"} = 1 ; $Allowed->{"OP"} = 1 ;
$Allowed->{"NN"} = 1 ;
$Allowed->{"NH"} = 1 ; $Allowed->{"HN"} = 1 ;
$Allowed->{"SH"} = 1 ; $Allowed->{"HS"} = 1 ;

my @list= util_read_list_sentences($listfile);
my $list = {};
map { s/\s*//g ; $list->{$_} = 1 ; } @list ;

my ($holo2apoinfo) ;
if(defined $holo2apo){
    ($holo2apoinfo) = util_parseHolo2Apo($holo2apo);
}

foreach my $protein (@list){

    my $orig = defined $holo2apo? $holo2apoinfo->{$protein} : $protein ;

    next if(! -e "$protein.allinfo.sorted");
    my $ifh = util_read("$protein.allinfo.sorted");
	my @lines ; 
    while(<$ifh>){
		chomp ;
		push @lines, $_;
    }
    close($ifh);

	if(@lines < 4){
		print "There are no contacts for $protein\n";
		next ;
	}

	## first find how many below 2.8
	my $cntreallyclose = 0 ;

	my $required = 4 ;

	my $doneclose = {};
	my $FINAL2PRINT = "";
	my $FINAL2PRINT1 = "";
	foreach my $line (@lines){
		my ($fullnm,$atom,$type,$dist) = split " ",$line ;

		
	    my ($name,$number,$atomnm) = ($fullnm =~ /([a-zA-Z]+)([0-9]+)([a-zA-Z]+)/);
		## sometimes numbers are negative
		next if(!defined $name);

		my $NAMENUM = $name . $number ;
		next if(exists $doneclose->{$NAMENUM});

		if($dist < 2.8 ){
		    $doneclose->{$NAMENUM} = 1 ;
			$cntreallyclose++;
			$required--;
		    $FINAL2PRINT = $FINAL2PRINT .  " $NAMENUM/$atomnm/$atom/$dist ";
		    $FINAL2PRINT1 = $FINAL2PRINT1 .  " $NAMENUM ";
		    #print "$protein $name/$number/$atomnm $atom,$type,$dist\n";
		}
	}
	my $str  = "$protein ";
	$str = $str . " $cntreallyclose ";

	my $donePolar = {};
	my $cntPolar = 0 ;
	foreach my $line (@lines){
	    last if(!$required);

		my ($fullnm,$atom,$type,$dist) = split " ",$line ;
		
	    my ($name,$number,$atomnm) = ($fullnm =~ /([a-zA-Z]+)([0-9]+)([a-zA-Z]+)/);

		## sometimes numbers are negative
		next if(!defined $name);

		my $NAMENUM = $name . $number ;

		next if(exists $doneclose->{$NAMENUM});
		next if(exists $donePolar->{$NAMENUM});
		

		next if(exists $NONPOLAR->{$name});
		next if(! exists $Allowed->{$type});
		next if(exists $backboneatomstable->{$atomnm});

		$donePolar->{$NAMENUM} = 1 ;
		if($dist < 3.8 ){
		    $donePolar->{$NAMENUM} = 1 ;
			$cntPolar++;
			$required--;
		    $FINAL2PRINT = $FINAL2PRINT .  " $NAMENUM/$atomnm/$atom/$dist ";
		    $FINAL2PRINT1 = $FINAL2PRINT1 .  " $NAMENUM ";
		    #print "$protein $name/$number/$atomnm $atom,$type,$dist\n";
		}
	}

	$str = $str . " $cntPolar ";

	my $sum = $cntreallyclose + $cntPolar ; 
	if($sum > 3){
		my $OFH = util_write("$orig.4.1.table.in");
		my $OFH1 = util_write("$orig.4.1.clasp.in");
		print $OFH "$FINAL2PRINT \n";
		print $OFH1 "$FINAL2PRINT1 \n";
		close($OFH);
	     print $ofh  " $str  \n";
	}

}

chmod 0777, $outfile ;

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
