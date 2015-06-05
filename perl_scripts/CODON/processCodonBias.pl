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
my ($getAAfreq,$infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "getAAfreq"=>\$getAAfreq ,
            "infile=s"=>\$infile ,
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
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);

my $sort = {};
while(<$ifh>){
	 next if(/^\s*$/);
	 my ($nm,@l) = split ; 
	$sort->{$nm} = $_ ;
}


my @ignoreforPYPU = qw (C D E F H K N Q Y );
my $ignoreforPYPU = {};
foreach my $q (@ignoreforPYPU){
	$ignoreforPYPU->{$q} = 1 ;
}
my $info = {};
my $AT = 0 ;
my $GC = 0 ;
my $PY = 0 ;
my $PU = 0 ;

my $fhAT = util_write("vAT");
my $fhGC = util_write("vGC");
my $fhPY = util_write("vPY");
my $fhPU = util_write("vPU");
my $fhdiffATGC = util_write("vDIFFATGC");

my $mapnm2sum = {};
my $total = 0 ;
foreach my $k (sort keys %{$sort}){
	$_ = $sort->{$k};
    my $at = 0 ;
    my $gc = 0 ;
    my $py = 0 ;
    my $pu = 0 ;
	 next if(/^\s*$/);
	 my ($nm,@l) = split ; 
	 #next if(@l eq 1);
	 my $sum = 0 ;
	 foreach my $i (@l){
	 	$i =~ s/=/ /;
		my @ll = split " ",$i ;
		$sum = $sum + $ll[1];
	 }
	 #my $QQ = int($sum/10000) ;
	 #print $fhnmsum "($nm,$QQ)\n";
	 $mapnm2sum->{$nm} = $sum ;
	 $total = $total + $sum ;


	 print $ofh "$nm ";
	 foreach my $i (@l){
	 	$i =~ s/=//;
		my @ll = split " ",$i ;
		my $percent = util_format_float((($ll[1]/$sum)*100),1);
		$percent =~ s/ //g;
		print $ofh " $ll[0]=$percent ";

	    if($ll[0] =~ /..T/ || $ll[0] =~ /..A/){
	    #if($ll[0] =~ /T../ || $ll[0] =~ /A../){
	    #if($ll[0] =~ /.T./ || $ll[0] =~ /.A./){
	 	   $AT = $AT + $ll[1] ;
	 	   $at = $at + $ll[1] ;
	    }
		else{
	 	   $GC = $GC + $ll[1] ;
	 	   $gc = $gc + $ll[1] ;
		}


		if(! exists $ignoreforPYPU->{$nm}){
	     if($ll[0] =~ /..T/ || $ll[0] =~ /..C/){
	 	   $PY = $PY + $ll[1] ;
	 	   $py = $py + $ll[1] ;
	     }
		else{
	 	   $PU = $PU + $ll[1] ;
	 	   $pu = $pu + $ll[1] ;
		 }
		}
		
	 }
	 print $ofh "\n";


	 # for TIKZ 
	 if(1){
	 	my $diffatgc = util_format_float($at - $gc,1) ;
	    print $fhAT "($nm, $at) \n";
	    print $fhGC "($nm, $gc) \n";
	    print $fhPY "($nm, $py) \n";
	    print $fhPU "($nm, $pu) \n";
	    print $fhdiffATGC "($nm, $diffatgc) \n";
	 }
	 

}

if(defined $getAAfreq){
   my $fhnmsum = util_write("aafreq");
   foreach my $k (sort keys %{$mapnm2sum}){
	   my $v = $mapnm2sum->{$k} ;
	   my $percent = util_format_float(100 * ($v/$total),1);
	print $fhnmsum "($k, $percent)\n";
   }
}

my $PYPU = util_format_float($PY/$PU,1);
my $ATGC = util_format_float($AT/$GC,1);

my $ATDIFFGC = $AT - $GC ;

my $ofhatdiffgc = util_write("atdiffgc");
print $ofhatdiffgc "$ATDIFFGC \n";


print "Total PY/PU (ignoring M and W): $PY $PU $PYPU \n";
print "Total AT/GC (ignoring M and W): $AT $GC $ATGC \n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
