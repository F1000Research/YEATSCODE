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
my ($infile,$outfile,$which_tech,$listfile,$protein);
my ($only2list,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
my $cutoff = 60 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "only2list"=>\$only2list ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
my $ifh = util_read($infile);


usage( "Need to give a input file name => option -infile ") if(!defined $infile);

my $infosimi = {};
my $infoiden = {};
while(<$ifh>){
	my ($a,$b,$simi,$iden) = split ;
	$infosimi->{$a} = {} if(!defined $infosimi->{$a});
	$infosimi->{$a}->{$b} = $simi;

	if(! defined $only2list){
	    $infosimi->{$b} = {} if(!defined $infosimi->{$b});
	    $infosimi->{$b}->{$a} = $simi;
	}
}


my $SETS = {};
my $FASTA2SET = {};
my $setnumber = 0 ; 
foreach my $k (keys %{$infosimi}){

	my $table = $infosimi->{$k}; 
	my @sorted = sort { $table->{$b} <=> $table->{$a} } (keys %{$table});
	foreach my $s (@sorted){
		my $v = $table->{$s};
		if($v < $cutoff){

	       if(! exists $FASTA2SET->{$k} && ! exists $FASTA2SET->{$s}){
	           $setnumber++;
	           $SETS->{$setnumber} = {};
	           $SETS->{$setnumber}->{$k} = 1 ;
	           $FASTA2SET->{$k} = $setnumber ;
            }
            else{
	            $setnumber = $FASTA2SET->{$k} if(defined $FASTA2SET->{$k});
	            $setnumber = $FASTA2SET->{$s} if(defined $FASTA2SET->{$s});
             }
			 $SETS->{$setnumber}->{$s} = 1 ;
	         $FASTA2SET->{$s} = $setnumber ;
		 }
		
		#print "$s/$v ..";
	}
	#print "\n";
}


foreach my $i (1..$setnumber){
	my $set = $SETS->{$i} ; 

	my @l = (keys %{$set});
	print $ofh "$i ";
	foreach my $x (keys %{$set}){
		print  $ofh "$x ";
	}
	print $ofh "\n";


	print "Checking...\n";
	foreach my $j (@l){
		my $atleastone = 0 ;
	    foreach my $k (@l){
		   next if($j eq $k);
		   my $v = $infosimi->{$j}->{$k};
		   if($v < $cutoff){
		   	   $atleastone = 1 ;
		   }
		   #print "$j $k $v\n";
	    }
		die if(!$atleastone);
	}
}

print "There are $setnumber with cutoff of $cutoff\n";
system("cat $outfile");

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
