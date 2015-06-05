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
my ($mergenames,$infile,$p1,$ignoreband,$cutoff,$p2,$mappingfile,$outfile,$which_tech,$ignorefile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "mergenames"=>\$mergenames ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "mappingfile=s"=>\$mappingfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "cutoff=i"=>\$cutoff ,
            "ignoreband=i"=>\$ignoreband ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -mappingfile ") if(!defined $mappingfile);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);


my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);

my $ignoretable = {};
if(defined $ignorefile){
   my @lll= util_read_list_sentences($ignorefile);
   map { s/\s*//g ; $ignoretable->{$_} = 1 ; } @lll ;
}

my $ifhmap = util_read($mappingfile);
my $mappingtable = {};

my $REALNAMEMAP = {};
while(<$ifhmap>){
	next if(/^\s*#/);
	 my (@l) = split ; 
	 my ($first) = $l[0];
	 my $realname = $first ;

	 $first =~ s/.ORF.*//;
	 $REALNAMEMAP->{$first} = $realname;

	 foreach my $i (@l){
	    my $realname = $i ;
	    $i =~ s/.ORF.*//;
		$REALNAMEMAP->{$i} = $realname;
	 	$mappingtable->{$i} = $first ;
	 }
}


my $NUM ; 
my $MAPPED = {};
my $CNT = 0 ;
while(<$ifh>){
	$CNT++;
	 my (@l) = split ; 
	 my ($name) = $l[0];
	 $name = uc($name);
	 $name =~ s/.ORF//;
	 next if(exists $ignoretable->{$name});

	 

	 if(!defined $NUM){
	 	$NUM = @l ;
	 }
	 else{
	 	die "$name not same" if(@l ne $NUM);
	 }

	 if(! exists $mappingtable->{$name}){
	 	die if(exists $MAPPED->{$name});
	 	$MAPPED->{$name} = [];
	 }
	 else{
	 	$name = $mappingtable->{$name} ;
	 	$MAPPED->{$name} = [] if(!exists $MAPPED->{$name});
	 }

	 push @{$MAPPED->{$name}}, \@l;
}


my $NN = (keys %{$MAPPED});
print "Started with $CNT, but ended with $NN\n";

my $warned = 0 ;
foreach my $k (keys %{$MAPPED}){

	my @list = @{$MAPPED->{$k}};
	
	my $first = shift @list ;
	my @first = @{$first};
	my $ZZZ = @first - 1 ;

	if(defined $mergenames && exists $REALNAMEMAP->{$k}){
		$k = $REALNAMEMAP->{$k};
		my $NUMM =  @list ;
		if(@list >  0){
			if(!$warned){
		       warn  "Merging more than one $k $NUMM, just printing one "   ;
		   }
		   $warned++;
		   @list  = ();
		}
		#my $second = $list[0];
		#$second = $REALNAMEMAP->{$second};
		#print "$k $second\n";

		
	}

	foreach my $l (@list){
		my @ll = @{$l};
		my $N = @ll -1 ;
		die if($N ne $ZZZ);
		foreach my $idx (1..$N){
			$first[$idx] = $first[$idx] + $ll[$idx];
		}
	}

	print $ofh "@first \n";
	
}

print "Merged $warned \n";


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
