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
my ($het2pdb,$pdbseqres,$infile,$size,$outfile,$which_tech,$listfile,$protein);
my (@expressions);
my $verbose = 0 ;
GetOptions(
            "het2pdb=s"=>\$het2pdb ,
            "pdbseqres=s"=>\$pdbseqres ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "verbose=i"=>\$verbose ,
            "size=i"=>\$size ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
print "Writing to $outfile - the commands to convert premon.in's\n";
my $ofh = util_write($outfile);

usage( "Need to give a input file name => option -het2pdb ") if(!defined $het2pdb);
usage( "Need to give a input file name => option -size ") if(!defined $size);
usage( "Need to give a input file name => option -pdbseqres ") if(!defined $pdbseqres);

my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;


print "READING $het2pdb\n";
my ($NOHET,$YESHET) = util_parseHETPDB($het2pdb);

print "READING $pdbseqres\n";
my ($info,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRES($pdbseqres,0);


my $Nnotfound = 0 ;
my $Nfound = 0 ;
my $done = {};

foreach my $protein (keys %{$YESHET}){

   print "Processing $protein \n" if($verbose);

   my $NUMHETS = $YESHET->{$protein};
   next if($NUMHETS > $size);

   if(!exists $mapChainedName2Name->{$protein}){
	   warn "protein $protein does not exist - have you specified the chain?\n";
	   next ;
   }
   my $seq = $mapChainedName2Name->{$protein};

   my @pdbswithseq  = @{$infoSeq2PDB->{$seq}} ;
   my $N = @pdbswithseq ;
   my $pdbstr = join ",", @pdbswithseq ;
   print "There are $N proteins with the same seq as $protein $pdbstr \n" if($verbose);



   my $found = 0 ; 


   ### Try Chain A first
   my $TOPRINT ; 
   foreach my $pdb (@pdbswithseq){
   	    next if(!($pdb =~ /....A/));
        if(exists $NOHET->{$pdb}){
			print "Found $pdb in NOHET\n" if($verbose);
		   ## just choose one
		   $TOPRINT = $pdb ;
		   $found = 1 ;
		   
		   last ;
	    }
   }

   if(!$found){
      foreach my $pdb (@pdbswithseq){
        if(exists $NOHET->{$pdb}){
		   ## just choose one
		   $TOPRINT = $pdb ;
		   $found = 1 ;
		   last ;
	    }
      }
   }

   if(!$found){
   	   $Nnotfound++; 
   }
   else{
		if(! exists $done->{$TOPRINT}){


	 	    print $ofh "$protein $TOPRINT $NUMHETS \n";
   	        $Nfound++;
		    $done->{$TOPRINT} = 1;
		}
   }
}
print "Found = $Nfound, not found = $Nnotfound\n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
