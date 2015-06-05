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
my ($infile,$p1,$p2,$outfile,$trs,$cutoff,$which_tech,$listfile,$protein);
my ($chooseone,$ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "trs=s"=>\$trs ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "chooseone"=>\$chooseone ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);

usage( "Need to give a trs -option -trs  ") if(!defined $trs);
usage( "Need to give a cutoff -option -cutoff  ") if(!defined $cutoff);

my $ofhannotate = util_open_or_append("genome.annotated.csv");
my $ofhmorethanone = util_open_or_append("genome.annotated.morethanone.csv");
my $ofhnone = util_open_or_append("genome.annotated.none.csv");
my $ofhcommands = util_open_or_append("genome.commands.csh");

my @l = <YY/$trs*>;
#print "@l \n";

if(@l ne 3){
	print "Not all found for $trs\n";
	exit ;
}
foreach my $infile (@l){
	if(-z $infile){
	    print "Not all found\n";
	    exit ;
	}
}


my $sort = {};
my $sortEVALS2SCORES = {};
my $sortEVALS2INFILE = {};
my $MAXSCORE = 0 ;
my $ALLSTRS = "";
foreach my $infile (@l){
   my ($fulllength,$TRSnm,$STRS,$LLLs,$SCORES,$IDENTITIES,$EVALUES) = util_ParseWebBlast($infile);


   if(!defined $fulllength){
   	  print "---File not complete, rerun $trs\n";
   	  exit ;
   }
   if(!defined $TRSnm){
   	  $infile =~ s/^YY.//;
      $ALLSTRS = $ALLSTRS . "\n" . "\t$infile ". "\t".  "***** No hits found ****";
   	  next ;
   }

   my @TRSnm = @{$TRSnm};
   my @STRS = @{$STRS};
   my @LLLs = @{$LLLs};
   my @SCORES = @{$SCORES};
   my @EVALUES = @{$EVALUES};
   my @IDENTITIES = @{$IDENTITIES};

   $infile =~ s/^YY.//;
   $ALLSTRS = $ALLSTRS . "\n" . "\t$infile ". $STRS[0];

   #$sort->{$SCORES[0]} = $STRS ;
   $sort->{$EVALUES[0]} = $STRS ;
   $sortEVALS2SCORES->{$EVALUES[0]} = $SCORES[0] ;
   $sortEVALS2INFILE->{$EVALUES[0]} = $infile ;
   if($MAXSCORE < $SCORES[0]){
   	  $MAXSCORE = $SCORES[0];
   }
   #print "$SCORES[0] kkkkkkkkk\n"
   
}

my $max = 0   ;
foreach my $k (sort { $b <=> $a} keys %{$sort}){
	if($k > $max){
		$max = $k ;
	}
}

my $cnt = 0 ;
my $onlyone  ;
my $trsrealname  ;
my $sortTRSLen = {};
foreach my $k (sort { $b <=> $a} keys %{$sort}){
	#if($k > $cutoff && (defined $chooseone && $k > $max/2)){
	my $SCORE = $sortEVALS2SCORES->{$k};
	my $infile = $sortEVALS2INFILE->{$k};
	$infile =~ s/\.blast/.ALL.1.fasta/;
	$infile =~ s/\._/.ORF_/;
	

	my ($str,$firstline) = util_readfasta("/home/sandeepc/DATA/rafael/walnut/2/FASTADIR_ORFNEW/$infile");
	my $len = length($str);
	$sortTRSLen->{$infile} = $len;

	#print "$infile $len\n";

	if(($k < $cutoff && $SCORE > $MAXSCORE/2) && defined $chooseone ){
		$cnt++;
	    my @l = @{$sort->{$k}};
		foreach my $i (@l){
			if(!($i =~ /(chromosome|hypothe|unnamed|uncharacterized)/i)){
				my @l = split " ", $i ;
				my $EEE = $l[@l -1 ];
				if($EEE > $cutoff){
					print "Warning: $trs Annotated, but with high E value $EEE $cutoff\n";
				    $cnt--;
				}
				else{
			        $onlyone = $i ;
			        $trsrealname = $infile ;
				}
				last ;
			}

		}

	}
	#print "$v $k \n";
}

if($cnt eq 1){
	print $ofhannotate "$trs $onlyone \n";
	print $ofhcommands "ln -s $trsrealname $trs.ALL.1.fasta # unique\n";
}
else{
	if($cnt eq 0){
	    print $ofhnone "$trs $ALLSTRS\n";
		my @sorted = sort {$sortTRSLen->{$b} <=> $sortTRSLen->{$a}} (keys %{$sortTRSLen});
		my $trsrealname = $sorted[0];
		my $len = $sortTRSLen->{$trsrealname};
	    print $ofhcommands "ln -s $trsrealname $trs.ALL.1.fasta # none\n";
	}
	elsif($cnt > 1){
	    print $ofhmorethanone "$trs $ALLSTRS\n";
	}
}

   
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
