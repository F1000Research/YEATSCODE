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
my ($fastadir,$infile,$p1,$p2,$orfdir,$outfile,$cutoff,$listfile,$orf,$trs);
my ($fastafile,@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
GetOptions(
            "listfile=s"=>\$listfile ,
            "trs=s"=>\$trs ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "orf=s"=>\$orf ,
            "fastadir=s"=>\$fastadir ,
            "outfile=s"=>\$outfile ,
            "orfdir=s"=>\$orfdir ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a output file name => option -orfdir ") if(!defined $orfdir);
my $ofh = util_write($outfile);
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a fastadir -option -fastadir  ") if(!defined $fastadir);

my $info = {};

my $ignoredduetounknown = 0 ;
my $ignoredduetosize = 0 ;
my $actuallprocessed = 0 ;
my $IFH = util_read($listfile);
while(<$IFH>){
	next if(/^\s*$/);
	 my ($trs,$orf) = split ; 
	 my $infile = "ORF/$trs.orf" ;
	 my $fastafile = "$fastadir/$trs.ALL.1.fasta";
	 ProcessOne($trs,$infile,$fastafile,$orf);
}


sub ProcessOne{
  my ($trs,$infile,$fastafile,$orf) = @_ ;
  my $ifh = util_read($infile);
  while(<$ifh>){
	  if(/$orf/){
	     s/-//g;
	     s/\[//g;
	     s/\]//g;
		 my @l = split ;
		 my $start = $l[1];
		 my $end = $l[2];

		 my $isreverse = 0 ;
		 if(/REVERSE/){
		 	$isreverse = 1 ;
		 }

		 my $SSS = "";
         while(<$ifh>){
		 	chomp;
		 	if(/^\s*>/){
				last ;
			}
			$SSS = $SSS . $_ ;
			
		 }
		 my $LLL = length($SSS);
		 if(defined $cutoff && $LLL < $cutoff){
		 	#print "Ignoring $trs\n";
			$ignoredduetosize++;
			next ;
		 }
		 
		 my ($str,$firstline) = util_readfasta($fastafile);
		 if($str =~ /(R|Y|K|M|S|W|B|D|H|V|N)/){
		 	$ignoredduetounknown++;
			next ;
		 }
		 $actuallprocessed++;


         if($isreverse){
		 	my $tmp = $start ;
			$start = $end ;
			$end = $tmp ;
		 }


		 my $retstr = util_extractSliceFromFasta($str,$start,$end);
		 my $len = length($retstr);



         if($isreverse){
			$retstr = util_getComplimentaryString($retstr) ;
		 }

		 my $rlen = length($retstr);
		 die "Something wrong" if($LLL*3 ne $len);
		 print "Extracting $start $end , got $LLL for aa and $len for nt. isreverse = $isreverse and $rlen = rlen\n" if($verbose);


		 my $codontable = util_getCodonTable();

		 my @NT = ($retstr =~ /(...)/g);
		 my @AA = ($SSS =~ /(.)/g);
		 my $N = @AA ;
		 my $NTN = @NT ;
		 die "$trs $N =n an $NTN = NTN" if($N ne $NTN);
		 while(@NT){
		 	my $a = shift @NT ;
		 	my $b = shift @AA ;
			if(!defined $codontable->{$a}){
				die "$a";
			}

			die "$a $b while this is in $codontable->{$a} " if($codontable->{$a} ne $b);
			#print "$a $b\n";
			if(!exists $info->{$b}){
				$info->{$b} = {};
			}
			if(!exists $info->{$b}->{$a}){
				$info->{$b}->{$a} = 0 ;
			}
			$info->{$b}->{$a} = $info->{$b}->{$a} + 1 ;
		 }

		 last ;

	  }
  }
  close($ifh);
}

foreach my $aa (sort keys %{$info}){
	my $tab = $info->{$aa} ;
    print $ofh "$aa ";
	foreach my $k (keys %{$tab}){
		my $v = $tab->{$k} ;
		print $ofh " $k=$v ";
	}
	print $ofh "\n";
}

print "actuallprocessed ignoredduetounknown ignoredduetosize \n";
print "$actuallprocessed $ignoredduetounknown $ignoredduetosize \n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
