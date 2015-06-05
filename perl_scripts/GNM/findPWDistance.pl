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
my ($annotatedfile,$infile,$p1,$mapfile,$cutoff,$p2,$outfile,$which_tech,$ignorefile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "annotatedfile"=>\$annotatedfile ,
            "cutoff=i"=>\$cutoff ,
            "mapfile=s"=>\$mapfile ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -ignorefile ") if(!defined $ignorefile);


$outfile = "$infile.PW";
print "Writing to $outfile\n";
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);

my $meaninfo = {};
my $sdindfo = {};

my @ignorelist= util_read_list_sentences($ignorefile);
my $ignorelist = {};
map { s/\s*//g ; $ignorelist->{lc($_)} = 1 ; } @ignorelist ;
#my @NAMES = qw(CE CI CK EM FL HC HL HP HU IF LE LM LY PK PL PT RT SE TZ VB);
my @NAMES = qw (CE  CI  CK  EM  FL  HC  HL  HP  HU  IF  LE  LM  LY  PK  PL  PT  RT  SE  TZ  VB);
my @NAMESSAVED = qw (CE  CI  CK  EM  FL  HC  HL  HP  HU  IF  LE  LM  LY  PK  PL  PT  RT  SE  TZ  VB);
		

my $ofhpw = util_write("pairwisenames");

my @PWnames ;
my $NAMESTABLE = {};
while(@NAMES){
	my $a = shift @NAMES ;
	$NAMESTABLE->{$a} = 0;
	foreach my $i (@NAMES){
		print $ofhpw "$a-$i\n";

		my @l;
		push @l, $a;
		push @l, $i;
		my @sorted = sort @l;

		push @PWnames, "$sorted[0]-$sorted[1]";
	}
}

my @NUMS ;
my $firstone = 1 ;
while(<$ifh>){
	 my (@l) = split ; 
	 my ($name) = shift @l ;

	 if(defined $annotatedfile){
	 	shift @l ;
	 	shift @l ;
	 	shift @l ;
	 	shift @l ;
	 	shift @l ;
	 }
	 if(exists $ignorelist->{$name}){
	 	print "Ignoring not existent $name\n";
		next ;
	 }
	 my $N = @l - 1;

	 foreach my $idx (0..$N){

		my $CNT = 0 ;
		while(@l){
		    my $a = shift @l ;
			foreach my $i (@l){
				my $diff = abs ($a - $i);
				print "$diff=diff $CNT\n" if($verbose);
				if($firstone){
					push @NUMS, $diff ;
				}
				else{
					$NUMS[$CNT] = $NUMS[$CNT] + $diff ;
				}
			    $CNT++;
			}
		}
	 }
	 $firstone = 0 ;
}


my $NN = @NUMS ;
die "$NN " if($NN ne @PWnames);
print "There were $NN items\n";
$, = "\n";

my $VALUES = {};
my $max = 0 ;
foreach my $idx (1..$NN){
	my $realidx  = $idx -1 ;
    #print $ofh "$PWnames[$realidx] $NUMS[$realidx] \n"; 
	$VALUES->{$PWnames[$realidx]} = $NUMS[$realidx]  ;
	if($NUMS[$realidx] > $max){
		$max = $NUMS[$realidx] ;
	}
}

my @sorted = sort { $VALUES->{$a} <=> $VALUES->{$b}} (keys %{$VALUES});

my $GROUPCNT = 0 ;
my $GROUPS = {};


my $NEWVALUES = {};
my $cnt = 0 ;
foreach my $k (@sorted){
	$cnt++;
	my $v = int(100*($VALUES->{$k}/$max));
	print $ofh "$k $v\n";
	$NEWVALUES->{$k} = $v ;

	if($cnt < $cutoff){
		#print "$k $v lll \n";
	    InsertOrCreateGroup($k)  ;
	}
}

print "There were $GROUPCNT groups\n";

foreach my $key (keys %{$GROUPS}){
	my $group = $GROUPS->{$key};
	print "$group: ";

	my @l ;
	foreach my $nm (keys %{$group}){
		print " $nm ";
		push @l, $nm ;
		$NAMESTABLE->{$nm} = 1 ;
	}
	print " \n";


	while(@l){
		my $a = shift @l ;
		foreach my $i (@l){
		   my @ll;
		   push @ll, $a;
		   push @ll, $i;
		   my @sorted = sort @ll;
		   my $newnm = "$sorted[0]-$sorted[1]";
		   my $newval = $NEWVALUES->{$newnm};
		   print " $newnm $newval ...";
		}
		print "\n";
	}
}

foreach my $k (keys %{$NAMESTABLE}){
	if(! $NAMESTABLE->{$k}){
		print "$k not there \n";
	}
}

sub InsertOrCreateGroup{
	my ($k) = @_ ;
	my @l = split "-", $k;

	my $found = 0 ;
	foreach my $key (keys %{$GROUPS}){
			my $group = $GROUPS->{$key};
			if(exists $GROUPS->{$key}->{$l[0]} || exists $GROUPS->{$key}->{$l[1]}){
				if($found){
					if(exists $GROUPS->{$key}->{$l[0]}){
					    print "$l[0] is a link\n";
					 }
					 if(exists $GROUPS->{$key}->{$l[1]}){
					    print "$l[1] is a link\n";
					 }
				}
				$GROUPS->{$key}->{$l[0]} = 1 ;
				$GROUPS->{$key}->{$l[1]} = 1 ;
				$found = 1 ;
			}
			#last if($found);
	}
	if(!$found){
        $GROUPCNT++;
	    $GROUPS->{$GROUPCNT} = {};
	    $GROUPS->{$GROUPCNT}->{$l[0]} = 1 ;
	    $GROUPS->{$GROUPCNT}->{$l[1]} = 1 ;
	}
}



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
