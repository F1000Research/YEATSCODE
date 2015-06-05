#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use PDB;
use ConfigPDB;
use MyGeom;
# just test

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($infile,$outfile,$which_tech,$listfile,$protein);
my (@expressions,$config);
my $howmany ;
my $verbose = 1 ;
my $LARGESTDISTALLOWED = 30 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "config=s"=>\$config,
            "howmany=i"=>\$howmany ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a config file name => option -config ") if(!defined $config);
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
#usage( "Need to give a howmany -option -howmany  ") if(!defined $howmany);
usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $PWD = cwd;

my  ($seconds, $microseconds) = gettimeofday;


    my @proteins ; 
    push @proteins, $protein ; 
	my $i = $protein;
    die "No apbs" if( ! -e "$APBSDIR/$i/$i.pqr");
    die "No apbs"  if( ! -e "$APBSDIR/$i/pot1.dx.atompot" && ! -e "$APBSDIR/$i/pot1.dx.atompot");
    my @info = util_ReadPdbs($PDBDIR,$APBSDIR,1,@proteins) ; 
    my $info = shift @info ;
    my $pdb1 = $info->{PDBOBJ};
    my $pqr1 = $info->{PQR};
    my $pots1 = $info->{POTS};

my @list= util_read_list_sentences($listfile);
my $list = {};
my @reallist ; 
#map {my @j = split ; push @reallist, @j;  } @list ;
map { s/,/ /g ; my @j = split ; push @reallist, @j;  } @list ;

############ Sort based on single letter only 
my @tmplist ;

my @LIST2CHECK1 ;
foreach my $i (@reallist){
	my ($name,$number) = ($i =~ /([a-zA-Z]+)([0-9]+)/);
	$name = uc($name);

	my ($res) = $pdb1->GetResidueIdx($number);
	my $RESNAME = $res->GetName();
	if($name ne $RESNAME){
		unlink $outfile ;
		die "Residue names at number $number do not match - $name and $RESNAME. Removing $outfile";
	}

	my $len = length($name); 
	die "Wrong length $len" if($len != 1 && $len != 3);
	$name = $pdb1->GetSingleLetter($name);
	my $newnm = $name . $number ; 
	push @tmplist, $newnm ; 
}
my @sortlist = sort @tmplist ;


my @chain;
my @atoms ; 
my $done ;
my $original ;
my $CNT = 0 ; 
my $finalstringlist = "";
my $realstring = "";
my $resultstyle = "";
my @LIST2CHECK2 ;
foreach my $i (@sortlist){
	print "$i\n";
	$realstring = $realstring .  $i . " " ;
	$CNT++ ;
	$i =~ s/,//g;
	my ($name,$number) = ($i =~ /([a-zA-Z]+)([0-9]+)/);
	$name = uc($name);
	#print "$name,$number \n";
	my $len = length($name); 
	die "Wrong length $len" if($len != 1);
	$finalstringlist = $finalstringlist . $name ;
	if($len == 1){
	      $name = $pdb1->GetThreeLetter($name);
	}


	my $configtmp = "$outfile.config$CNT";
	if(-e $configtmp){
		print "using nonstandard \n";
         ConfigPDB_Init($configtmp);
	}
	else{
		print "using standard \n";
        ConfigPDB_Init($config);
	}

	my ($res) = $pdb1->GetResidueIdx($number);

	my $RESNAME = $res->GetName();
	push @LIST2CHECK2, $RESNAME ;
	my $type = ConfigPDB_GetAtom($res->GetName()) or die;
	my ($atom1) = $pdb1->GetAtomFromResidueAndType($number,$type);
	die "$number,$type not found" if(!defined $atom1);
	$resultstyle = $resultstyle .   $name . "/" . $number . "/" . $type  .  " " ;
	$atom1->Print();
	push @atoms, $atom1 ;
}


## keep different scope - sorted lists
if(1){
    my @distlist = sort @{$pdb1->DistanceInGivenSetOfAtoms(\@atoms)};
    my @pdlist =sort  @{$pdb1->PDInGivenSetOfAtoms(\@atoms,$pqr1,$pots1)};
    my $N = @distlist - 1 ;
    my $largestdist = $distlist[$N];
    if($largestdist > $LARGESTDISTALLOWED){
	    $, = " , " ;
	    print "Distances are @distlist \n";
        die "Error: $protein: Largest distance = $largestdist  more than $LARGESTDISTALLOWED"  ;
    }
}

$, =  " " ;


my @distlist = @{$pdb1->DistanceInGivenSetOfAtoms(\@atoms)};
my @pdlist = @{$pdb1->PDInGivenSetOfAtoms(\@atoms,$pqr1,$pots1)};

my $ofh = util_write($outfile);
my $ofhresultstyle = util_write("$protein.resultstyle");
print $ofh  "STR $finalstringlist \n";
print $ofh  "REALSTRING $realstring \n";
print $ofh  "RESULTSTYLE $resultstyle \n";
#print $ofh  "RESULTSTYLE \n";
print $ofh  "D @distlist \n";;
print $ofh  "PD @pdlist \n"; ;


print $ofhresultstyle  "$resultstyle \n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
