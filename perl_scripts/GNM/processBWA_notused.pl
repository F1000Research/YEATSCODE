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
my ($infile,$p1,$ignoreband,$cutoff,$p2,$outfile,$which_tech,$listfile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "cutoff=i"=>\$cutoff ,
            "ignoreband=i"=>\$ignoreband ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -listfile ") if(!defined $listfile);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a output file name => option -cutoff ") if(!defined $cutoff);
usage( "Need to give a output file name => option -ignoreband ") if(!defined $ignoreband);


my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);

my @list= util_read_list_sentences($listfile);
my $list = {};
map { s/\s*//g ; $list->{lc($_)} = 1 ; } @list ;
#my @NAMES = qw(CE CI CK EM FL HC HL HP HU IF LE LM LY PK PL PT RT SE TZ VB);
my @NAMES = qw (CE  CI  CK  EM  FL  HC  HL  HP  HU  IF  LE  LM  LY  PK  PL  PT  RT  SE  TZ  VB);


my @FFF ; 
my @FFFNEG ; 
foreach my $n (@NAMES){
	my $ooo = "SPLIT.$n";
	my $ofhmm = util_write($ooo);
	push @FFF, $ofhmm;

	$ooo = "NEGSPLIT.$n";
	$ofhmm = util_write($ooo);
	push @FFFNEG, $ofhmm;
}

my $HIGHCUTOFF = $cutoff + $ignoreband;
my $POS = {};
my $NEG = {};
my $MERGED = {};
my $IGN = 0 ;
my $NAME2TOTALCNT = {};
while(<$ifh>){
	 my (@l) = split ; 
	 my ($name) = shift @l ;
	 if(exists $list->{$name}){
	 	#print "Ignoring not existent $name\n";
		$IGN++;
		next ;
	 }
	 my $N = @NAMES -1 ;
	 foreach my $idx (0..$N){
	 	my $code = $NAMES[$idx];
		my $a = shift @l ;
		#my $b = shift @l ;

	    $MERGED->{$name} = {} if(! defined $MERGED->{$name});
	    $MERGED->{$name}->{$code} = 1 ;

		$NAME2TOTALCNT->{$name} = 0 if(!defined $NAME2TOTALCNT->{$name});
		$NAME2TOTALCNT->{$name} = $NAME2TOTALCNT->{$name} + $a ;

		if($a > $HIGHCUTOFF){
	 	    my $OFH = $FFF[$idx];
		    print $OFH "$name $a \n";
			$POS->{$name} = {} if(! defined $POS->{$name});
			$POS->{$name}->{$code} = 1 ;
		}
		elsif($a < $cutoff){
	 	    my $OFH = $FFFNEG[$idx];
		    print $OFH "$name $a \n";
			$NEG->{$name} = {} if(! defined $NEG->{$name});
			$NEG->{$name}->{$code} = 1 ;
		}
	 }
	 die if(@l);
}


my $PW = {};
foreach my $name (keys %{$POS}){
	my $tab = $POS->{$name};
	my @codes = (sort keys %{$tab});

	while(@codes){
		my $a = shift @codes ;
	    foreach my $b (@codes){
			$PW->{$a} = {} if(!exists $PW->{$a});
			if(! exists $PW->{$a}->{$b}){
				$PW->{$a}->{$b} = 1 ;
			}
			else{
				$PW->{$a}->{$b} = $PW->{$a}->{$b} + 1 ;
			}
			
	    }
	}
}




if(0){
foreach my $name (keys %{$NEG}){
	my $tab = $NEG->{$name};
	my @codes = (sort keys %{$tab});

	while(@codes){
		my $a = shift @codes ;
	    foreach my $b (@codes){
			$PW->{$a} = {} if(!exists $PW->{$a});
			if(! exists $PW->{$a}->{$b}){
				$PW->{$a}->{$b} = 1 ;
			}
			else{
				$PW->{$a}->{$b} = $PW->{$a}->{$b} + 1 ;
			}
			
	    }
	}
}
}

foreach my $k1 (keys %{$PW}){
	my $tab = $PW->{$k1};
    foreach my $k2 (keys %{$tab}){
		my $v = $tab->{$k2};
		print $ofh "$k1 $k2 $v\n";
	}
	
}

my $FHHALL = util_write("transcripts.bycode.ALL.$HIGHCUTOFF");
my $FHHFEW = util_write("transcripts.bycode.FEW.$HIGHCUTOFF");
foreach my $name (keys %{$POS}){
	my $tab = $POS->{$name};
	my @codes = (sort keys %{$tab});
	my $NUM = @codes ;
	my $totalcnt = $NAME2TOTALCNT->{$name};
	if($NUM < 4){
	     print $FHHFEW "$name $NUM $totalcnt @codes\n";
	}
	elsif($NUM > 17){
	     print $FHHALL "$name $NUM $totalcnt @codes\n";
	}
}

print "Ignored $IGN\n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
