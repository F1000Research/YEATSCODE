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
my ($duplicate,$blastout,$infile,$p1,$p2,$outfile,$trs,$cutoff,$which_tech,$listfile,$protein);
my ($chooseone,$ignorefile,$donone,@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "trs=s"=>\$trs ,
            "blastout=s"=>\$blastout ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "donone"=>\$donone ,
            "duplicate"=>\$duplicate ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);

usage( "Need to give a trs -option -trs  ") if(!defined $trs);
usage( "Need to give a blastout -option -blastout  ") if(!defined $blastout);
usage( "Need to give a cutoff -option -cutoff  ") if(!defined $cutoff);

my $tableuniquelyannotated = {};
{
system("touch genome.list");
my @list= util_read_list_sentences("genome.list");
map {  my @ll = split ;  $tableuniquelyannotated->{$ll[0]} = 1 ; } @list ;
}
if(exists $tableuniquelyannotated->{$trs}){
	exit ;
}

my $ofhannotate = util_open_or_append("genome.annotated.csv");
my $ofhmorethanone = util_open_or_append("genome.annotated.morethanone.csv");
my $ofhnone = util_open_or_append("genome.annotated.none.csv");
my $ofhcommands = util_open_or_append("genome.commands.csh");
my $ofhuniquelyannotated = util_open_or_append("genome.list");
my $ofhwarn = util_open_or_append("genome.warn");


print "Writing bunch of genome* files\n" if($verbose);

my $trsfasta = "$trs.ALL.1.fasta";

my @l = <$blastout/$trs.*>;
#print "@l \n";

if(@l eq 4){
	print $ofhwarn "$trs repeats?\n";
	shift @l ;
}

if(@l ne 3){
	print "Not all found for $trs\n";
	exit ;
}
foreach my $infile (@l){
	if(-z $infile){
	    print "Not all found $trs\n";
	    exit ;
	}
}


my $sort = {};
my $sortEVALS2SCORES = {};
my $sortEVALS2ALLSTRS = {};
my $sortEVALS2INFILE = {};
my $MAXSCORE = 0 ;
my $ALLSTRS = "";

my $JUNKEVAL = 1 ;
foreach my $infile (@l){
   my ($fulllength,$TRSnm,$STRS,$LLLs,$SCORES,$IDENTITIES,$EVALUES) = util_ParseWebBlast($infile);


   if(!defined $fulllength){
   	  print "---File not complete, rerun $trs\n";
   	  exit ;
   }
   if(!defined $TRSnm){
   	  $infile =~ s/^$blastout.//;
          $ALLSTRS = $ALLSTRS . "\n" . "\t$infile ". "\t".  "***** No hits found ****";
   $sort->{$JUNKEVAL} = {};
   $sortEVALS2INFILE->{$JUNKEVAL} = $infile ;
   $JUNKEVAL++;
   	  next ;
   }

   my @TRSnm = @{$TRSnm};
   my @STRS = @{$STRS};
   my @LLLs = @{$LLLs};
   my @SCORES = @{$SCORES};
   my @EVALUES = @{$EVALUES};
   my @IDENTITIES = @{$IDENTITIES};

   $infile =~ s/^$blastout.//;
   $ALLSTRS = $ALLSTRS . "\n" . "\t$infile ". $STRS[0];

   #$sort->{$SCORES[0]} = $STRS ;
   $sort->{$EVALUES[0]} = $STRS ;
   $sortEVALS2SCORES->{$EVALUES[0]} = $SCORES[0] ;
   $sortEVALS2ALLSTRS->{$EVALUES[0]} = $STRS ;
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
	

	#my ($str,$firstline) = util_readfasta("/home/tsbutter/FASTADIR_ORFNEW/$infile");
	my ($str,$firstline) = util_readfasta("FASTADIR_ORFNEW/$infile");
	my $len = length($str);
	$sortTRSLen->{$infile} = $len;

	#print "$infile $len\n";

	if(($k < $cutoff && $SCORE > $MAXSCORE/2)){
		$cnt++;
	    my @l = @{$sort->{$k}};
		my $firstone ;
		foreach my $i (@l){
			$firstone = $i if(!defined $firstone);
			if(!($i =~ /(chromosome|hypothe|unnamed|uncharacterized)/i)){
				my @l = split " ", $i ;
				my $EEE = $l[@l -1 ];
				if($EEE > $cutoff){
					print $ofhwarn "$trs Annotated, but with high E value $EEE $cutoff\n";
					$onlyone = $firstone . "\tNo function found";
			                 $trsrealname = $infile ;
				}
				else{
			        $onlyone = $i ;
			        $trsrealname = $infile ;
				}
				last ;
			}

		}
		if(!defined $onlyone){
					print $ofhwarn "$trs not annotated at all\n";
					$onlyone = $firstone;
			                 $trsrealname = $infile ;
			
		}

	}
	#print "$v $k \n";
}

print "$cnt lll \n" if(defined $donone);

if($cnt eq 1){
	print $ofhannotate "$trs $onlyone \n";
	print $ofhcommands "ln -s $trsrealname $trsfasta # unique\n";
	print $ofhuniquelyannotated "$trs unique-$cutoff\n";
}
else{
	#if($cnt eq 0 && defined $donone)
	if(defined $donone){
	        #print $ofhnone "$trs $ALLSTRS\n";
		my @sorted = sort {$sortTRSLen->{$b} <=> $sortTRSLen->{$a}} (keys %{$sortTRSLen});

		my $trsrealname = $sorted[0];
		my $len = $sortTRSLen->{$trsrealname};
	        print $ofhcommands "ln -s $trsrealname $trsfasta # none\n";
	        print $ofhuniquelyannotated "$trs none-$cutoff\n";
	}
	elsif($cnt > 1 && defined $duplicate){
		my @sorted = (sort {$a <=> $b} keys %{$sortEVALS2INFILE});

		my $X = shift @sorted;
		my $Y = shift @sorted;
		my @SCORESA = @{$sortEVALS2ALLSTRS->{$X}};
		my @SCORESB = @{$sortEVALS2ALLSTRS->{$Y}};
		my $Atable = {};

		## make table of lower
		foreach my $str (@SCORESB){
			$str =~ s/\|/ /g;
			my @l = split " ", $str ;
			$Atable->{$l[1]} = 1 ;
		}

		my $foundone = 0 ;
		foreach my $str (@SCORESA){
			my $orig = $str ;
			$str =~ s/\|/ /g;
			my @l = split " ", $str ;
			my $N = @l - 1;
			my $EEE = $l[$N];
			if(exists $Atable->{$l[1]} && $EEE < 0.0000000001){
				print "Found $l[1] $str \n";
				$foundone = 1 ;
	            print $ofhannotate "$trs $onlyone \n";
	            print $ofhcommands "ln -s $trsrealname $trsfasta # unique fix $l[1]\n";
	            print $ofhuniquelyannotated "$trs unique-assemblyerror$cutoff\n";
				last ;
			}
		}



		if(!$foundone){
	        print $ofhmorethanone "$trs $ALLSTRS\n";
		    my $a = $sortEVALS2INFILE->{$X};
		    my $b = $sortEVALS2INFILE->{$Y};
	        my $trsa = $trs. "_a." . "ALL.1.fasta";
	        my $trsb = $trs. "_b." . "ALL.1.fasta";
                $a =~ s/\./.ORF/;
                $b =~ s/\./.ORF/;

                $a =~ s/\.blast/\.ALL.1.fasta/;
                $b =~ s/\.blast/\.ALL.1.fasta/;
		    print $ofhcommands "ln -s $a $trsa # duplicate 2 \n";
		    print $ofhcommands "ln -s $b $trsb # duplicate 2 \n";
		    if($cnt eq 2){
	            print $ofhuniquelyannotated "$trs duplicate-2$cutoff\n";
    
    
		    }
		    else{
			    die if($cnt ne 3);
		        my $c = $sortEVALS2INFILE->{shift @sorted};
                        $c =~ s/\./.ORF/;
                        $c =~ s/\.blast/\.ALL.1.fasta/;
	                my $trsc = $trs . "_c." . "ALL.1.fasta";
		        print $ofhcommands "ln -s $c $trsc # duplicate 3\n";
	                print $ofhuniquelyannotated "$trs duplicate-3$cutoff\n";
		    }
		}
	}
}

   
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
