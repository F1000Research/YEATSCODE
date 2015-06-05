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
my ($infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein,$appendfile);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "appendfile=s"=>\$appendfile ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a output file name => option -listfile ") if(!defined $listfile);
usage( "Need to give a output file name => option -appendfile ") if(!defined $appendfile);
my $ofh = util_write($outfile);
my $ofhLIST = util_append($appendfile);

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();


my $mapping = {};
my $revmapping = {};
while(<$ifh>){
	 my ($nm,$junk) = split ; 
	 $nm =~ s/.ORF.*//;
	 $junk =~ s/.ORF.*//;
	 $mapping->{$nm} = $junk ;
	 $revmapping->{$junk} = $nm ;
}

$ifh = util_read($listfile);

my $PRINTED = 0 ;
my $DIDNOTCONCAT = 0 ;
while(<$ifh>){
	if(/2 =N/){
	   #s/.* =N//;
	   #print $ofh "$_";
	   #next ;
		
	}
	s/.* =N//;
	my @l = split ;

	foreach my $i (@l){
	    $i =~ s/.ORF.*//;
	}

    my $first ;
	my $last ;
	foreach my $i (@l){
		if(exists $mapping->{$i} && ! exists $revmapping->{$i}){
			$first = $i ;
		}
		if(exists $revmapping->{$i} && ! exists $mapping->{$i}){
			$last = $i ;
		}
	}


	my $x = $first ;
	$x =~ s/.ORF.*//;
	my ($fe,$JJJJ) =  util_readfasta("$FASTADIR/$x.ALL.1.fasta");
	
    my $newnm = "$x.MERGED" ; 
	print $ofhLIST "$newnm\n";
	my $PRINT = 1;
	while(1){
		print $ofh "$x   ";
		last if(!exists $mapping->{$x});
		print "$x ....\n";
		$x = $mapping->{$x};
		$x =~ s/.ORF.*//;
		my ($fb,$QQQQ) =  util_readfasta("$FASTADIR/$x.ALL.1.fasta");

		$fe = util_ConcatTwoStringsWithCommonBeginAndEnd($fe,$fb,100);
        #$newnm = $newnm . ".$x";
		if($fe eq -1){
			print "$newnm does not concat \n";
			$PRINT = 0 ;
			$DIDNOTCONCAT++;
			last  ;
		}


	}
	if($PRINT){
	#print $ofh "\n";
	$PRINTED++;
    my $FH = util_write("$newnm.ALL.1.fasta");
    print $FH ">$newnm\n";
    print $FH "$fe\n";
    close ($FH);
	}

}

print "DIDNOTCONCAT = $DIDNOTCONCAT, PRINTED = $PRINTED\n";


chmod 0777, $outfile ;

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
