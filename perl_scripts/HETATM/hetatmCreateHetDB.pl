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
my ($holo2apo,$pdbseqres,$hetatm,$het2pdb,$pdbseqres,$infile,$outfile,$which_tech,$listfile,$protein);
my (@expressions);
my $size ;
my $verbose = 0 ;
GetOptions(
            "het2pdb=s"=>\$het2pdb ,
            "pdbseqres=s"=>\$pdbseqres ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "hetatm=s"=>\$hetatm ,
            "holo2apo=s"=>\$holo2apo ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "size=i"=>\$size ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
print "Writing to $outfile - the commands to convert premon.in's\n";
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -pdbseqres ") if(!defined $pdbseqres);

usage( "Need to give a input file name => option -het2pdb ") if(!defined $het2pdb);
#usage( "Need to give a input file name => option -hetatm ") if(!defined $hetatm);
#usage( "Need to give a input file name => option -holo2apo ") if(!defined $holo2apo);

#usage( "Need to give a input file name => option -pdbseqres ") if(!defined $pdbseqres);
#print "READING $pdbseqres\n";
#my ($info,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRES($pdbseqres,0);

my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;


print "READING $pdbseqres\n";
my ($info,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRES($pdbseqres,0);

print "READING $het2pdb\n";
my ($NOHET,$YESHET,$HET2PDB,$HET2PDBSIZE) = util_parseHETPDB($het2pdb);
my ($holo2apoinfo) ;
if(defined $holo2apo){
    ($holo2apoinfo) = util_parseHolo2Apo($holo2apo);
}

my $ifh = util_read($het2pdb);
my $found = 0 ;

my $samepdb = {};
while(<$ifh>){
	next if(/YESHET/);
	next if(/NOHET/);
	my ($hetatm,$SIZE,@list) = split ;
	next if(defined $size && $size > $SIZE);


	my $sequencesdone = {};
	foreach my $protein (@list){
		my $numhets = $YESHET->{$protein};
		next if($numhets > 1);


		my $seq = $mapChainedName2Name->{$protein};
		next if(!defined $seq || exists $sequencesdone->{$seq});
		$sequencesdone->{$seq} = 1 ;

		if(! defined $holo2apo || exists $holo2apoinfo->{$protein}){

            print $ofh "hetatmProc.pl -outf ooo  -ispolar 1 -ini 3.9 -max 3.9  -het $hetatm -p $protein\n";
            print $ofh "sort.pl -idx 3 -in $protein.allinfo -out $protein.allinfo.sorted \n";


			my $orig = defined $holo2apo ?  $holo2apoinfo->{$protein} : $protein ;
			if(exists $samepdb->{$orig}){
				$samepdb->{$orig} = $samepdb->{$orig} + 1;
			}
			else{
				$samepdb->{$orig} = 1 ;
			}
			$found++;
		}
	}
}

print "Found = $found\n";
foreach my $k (keys %{$samepdb}){
	my $v = $samepdb->{$k};
	if($v ne 1){
		print "SAME $k => $v\n";
   }
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
