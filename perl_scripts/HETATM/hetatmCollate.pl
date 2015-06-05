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
my ($infile,$outfile,$which_tech,$listfile,$hetatm);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "hetatm=s"=>\$hetatm ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
usage( "Need to give a hetatm pdb id -option -hetatm  ") if(!defined $hetatm);
my $CNT = 0 ; 
my $info = {};

my $stats = {};
my $fhs = {};
while(<$ifh>){
     next if(/^\s*$/);
	 my ($atom,$pdb,$res,$num,$type,$d) = split ; 
	 my $what = "$res,$num,$type,$d";
	 $info->{$atom} = {} if(!defined $info->{$atom});
	 $info->{$atom}->{PDB}->{$pdb} = [] if(!defined $info->{$atom}->{PDB}->{$pdb});

	 push @{$info->{$atom}->{PDB}->{$pdb}}, $what ;

	 ## save just the unique residue - and not each atoms 
	 my $str = "$pdb.$num";

	 $fhs->{$pdb} = util_write("$pdb.clasp.in") if(!defined $fhs->{$pdb});
	 my $fhclasp = $fhs->{$pdb};

	 if(!defined $stats->{$str}){
	     $stats->{$str} = $res ;
		 print $fhclasp "$res$num \n";
	 }

}

### does not make sense - since atoms numbers are not conserved
foreach my $atom (keys %{$info}){
	my $tab = $info->{$atom}->{PDB};
	print $ofh "$atom  \n";
	foreach my $pdb (keys %{$tab}){
		my $l = $tab->{$pdb};
		foreach my $what (@{$l}){
		     print $ofh "\t$what $pdb \n";
		}
	}
}


my $cntstats = {};
my $whatstr = {};
foreach my $k (keys %{$stats}){
	my $what = $stats->{$k};
	if(!defined $cntstats->{$what}){
		$cntstats->{$what} = 1;
		$whatstr->{$what} = $k;
	}
	else{
		$cntstats->{$what} = $cntstats->{$what} + 1; 
		$whatstr->{$what} = $whatstr->{$what} . "-" . $k;
	}
}


$, = "";
foreach my $what (keys %{$cntstats}){
	my $v = $cntstats->{$what};
	my $str = $whatstr->{$what};
	my @l = split "-", $str ;
	#print @l ;
	#die ;
	print "$what $v $str\n";
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
