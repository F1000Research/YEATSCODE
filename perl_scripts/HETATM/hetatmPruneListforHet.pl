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
print "Writing to $outfile - the commands to convert premon.in's\n";
my $ofh = util_write($outfile);

usage( "Need to give a input file name => option -het2pdb ") if(!defined $het2pdb);
usage( "Need to give a input file name => option -pdbseqres ") if(!defined $pdbseqres);
usage( "Need to give a input file name => option -size ") if(!defined $size);

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


my $Nnotfound = 0 ;
my $Nfound = 0 ;
my $done = {};

my $IGNOREPDBS = {};
if(0){
$IGNOREPDBS->{"4J4RA"} = 1;
$IGNOREPDBS->{"1G40A"} = 1;
$IGNOREPDBS->{"1G44A"} = 1;
$IGNOREPDBS->{"1G40B"} = 1;
$IGNOREPDBS->{"1G44B"} = 1;
$IGNOREPDBS->{"1G44C"} = 1;
}
foreach my $protein (@list){
   if(!exists $mapChainedName2Name->{$protein}){
	   die "protein $protein does not exist - have you specified the chain?\n";
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
   	    next if( exists $IGNOREPDBS->{$pdb});
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
   	    next if( exists $IGNOREPDBS->{$pdb});
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

	    my $ANNDIR = "ANNOTATE.4";

	 	print $ofh "preparePremonConfigs.pl -outf lll -lis $protein.$size.1.table.in -pr $TOPRINT -pol \$1 -size $size -anndir $ANNDIR\n";
	 	print $ofh "createPremoninput.pl -out $ANNDIR/$TOPRINT.4.\$1.premon.in -con \$CONFIGGRP -li $ANNDIR/$protein.in -pr $TOPRINT\n";

	 	print $ofh "printPairwise.pl -out 1 -c \$CONFIGGRP -ra 222 -pr $protein -in $TOPRINT.resultstyle\n";
	 	print $ofh "printPairwise.pl -out 2 -c \$CONFIGGRP -ra 222 -pr $TOPRINT -in $TOPRINT.resultstyle\n";
	 	print $ofh "cat 1 2 > \! $protein.$TOPRINT.table\n\n\n" ;

	 	print $ofhlistNOHET "$TOPRINT\n";
	 	print $ofhlistHET "$protein\n";
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
