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
my ($infile,$p1,$p2,$outfile,$tag,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "tag=s"=>\$tag ,
            "p1=s"=>\$p1 ,
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
usage( "Need to give a tag file name => option -tag ") if(!defined $tag);
my $ofh = util_write($outfile);

system ("touch logclose") if(! -e "logclose");
my $ofhlog = util_append("logclose");

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();

my $PWD = cwd;

my $ignoretable = {};
if(defined $ignorefile){
   my @lll= util_read_list_sentences($ignorefile);
   map { s/\s*//g ; $ignoretable->{$_} = 1 ; } @lll ;
}

my $name = "";
my $CNT = 0 ; 

my $ofhonlyone = util_write("onlyone");
my $hasmorethanoneofh = util_write("hasmorethanone");
my $ORIGLENGTH ;

my $total = 0 ;
my $hasmorethanone = 0 ;

my $info = {};

my $seenforeach; 

my $ALL = {};
my $hasmorethanonemap = {};
my $onlyonemap = {};
while(<$ifh>){
	if(/--/){
		die "$name $CNT $ORIGLENGTH" if(!$CNT);
		if($CNT eq 1){
		   print $ofhonlyone "$name\n";
		   $onlyonemap->{$name} = 1;
		}
		if(!exists $seenforeach->{$name}){
			die "There is no $name found \n";
		}
		next ;
	}
	chomp;
    if(/$tag/){
		s/$tag.//;
		s/.blast.out.*//;
		$name = $_ ;
		$ALL->{$name} = 1 ;
		$CNT = 0 ;
		$ORIGLENGTH = 0 ;
		$total++;
		$seenforeach = {};

    }
	else{
		my @l = split ;
		my $N = @l - 1 ;
		my $t = $l[0];
		#print "$t llll\n";
		#$t =~ s/;//;
		$t =~ s/.ORF.*;//;
		if(($t =~ /ORF.*ORF/)){
		   die ;
		}
		#print "$t kkkk\n";


	    #C53436_G2_I1.ORF_35;  unnamed protein product                           882   0.0   425 / 425 100
		my $percent = $l[$N];
		my $NUM  = $l[$N-3];

		#print "$t $name \n";
		if(!$CNT || $name eq $t){
			$ORIGLENGTH = $NUM;
		    $CNT++;
			$seenforeach->{$t} = 1 ;
		}
		elsif($percent > 95 && abs($ORIGLENGTH - $NUM) < 200 && $NUM > 100){
		    #print "$name $t $percent $NUM $SCORE \n";
			$info->{$name} = {} if(!defined $info->{$name});
			$info->{$name}->{$t} = 1 ;
		    $CNT++;
			$seenforeach->{$t} = 1 ;
			$hasmorethanone++;

		}
	}

}
close($ifh);


foreach my $k (sort keys %{$info}){
	foreach my $j (sort keys %{$info->{$k}}){
		if(exists $info->{$j}->{$k}){
			print $ofh "$k $j 0 \n";
			print $hasmorethanoneofh "$k\n";
			print $hasmorethanoneofh "$j\n";
			$hasmorethanonemap->{$k} = 1;
			$hasmorethanonemap->{$j} = 1;
		}
    }
}

foreach my $i (keys %{$ALL}){
	if(! exists $hasmorethanonemap->{$i} && ! exists $onlyonemap->{$i}){
	    print $ofhonlyone "$i\n";
	}
}

print "processed $total hasmorethanone = $hasmorethanone\n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
