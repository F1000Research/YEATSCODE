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
my ($het2pdb,$pdbseqres,$infile,$outfile,$which_tech,$listfile,$protein);
my (@expressions);
my $size ;
my $verbose = 0 ;
GetOptions(
            "het2pdb=s"=>\$het2pdb ,
            "pdbseqres=s"=>\$pdbseqres ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "size=i"=>\$size ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);

usage( "Need to give a input file name => option -het2pdb ") if(!defined $het2pdb);
usage( "Need to give a input file name => option -pdbseqres ") if(!defined $pdbseqres);

usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
my $ofhlistNOHET = util_write("$listfile.nohet");
my $ofhlistHET = util_write("$listfile.het");

my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;

my @list= util_read_list_sentences($listfile);
my $list = {};
map { s/\s*//g ; $list->{$_} = 1 ; } @list ;


print "READING $het2pdb\n";
my ($NOHET,$YESHET) = util_parseHETPDB($het2pdb);

print "READING $pdbseqres\n";
my ($info,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRES($pdbseqres,0);


foreach my $protein (@list){
   if(!exists $mapChainedName2Name->{$protein}){
	   die "protein $protein does not exist - have you specified the chain?\n";
   }
   my $seq = $mapChainedName2Name->{$protein};

   my @pdbswithseq  = @{$infoSeq2PDB->{$seq}} ;
   my $N = @pdbswithseq ;
   my $pdbstr = join ",", @pdbswithseq ;
   print "There are $N proteins with the same seq as $protein $pdbstr \n" ;
   foreach my $p (@pdbswithseq){
   	   print $ofh "$p\n";
   }

}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
