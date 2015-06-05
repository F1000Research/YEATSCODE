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
use POSIX qw(floor);
use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);

my $bufferintron = 9 ;
my $buffer = 0 ;
my $dointron = 0 ;

my $commandline = util_get_cmdline("",\@ARGV) ;
my ($query,$nameofsubject,$findifreverse,$subject,$exact,$infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($trs,$fasta,$ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "nameofsubject=s"=>\$nameofsubject ,
            "exact"=>\$exact ,
            "findifreverse"=>\$findifreverse ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "trs=s"=>\$trs ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "query=s"=>\$query ,
            "subject=s"=>\$subject ,
            "expr=s"=>\@expressions,
            "dointron=i"=>\$dointron ,
            "bufferintron=i"=>\$bufferintron ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);


usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -query ") if(!defined $query);
usage( "Need to give a input file name => option -subject ") if(!defined $subject);
usage( "Need to give a input file name => option -trs ") if(!defined $trs);
my $ifh = util_read($infile);

my $ofhlog = util_append("infointrons.csv");


my $ofhscaffoldlog = util_open_or_append("LOGS/$nameofsubject.log");
print "Appending to LOGS/$nameofsubject.log \n";
print $ofhscaffoldlog "\n---------- $trs ------------\n";



print "bufferintron = $bufferintron, and buffer = $buffer: you can change it if you want or give -exact which will make buffer 0 \n\n";

my $info = {};
my $CNT = 0 ; 
my $qstart ;
my $qend ;
my $sstart ;
my $send ;

my $DONOW = 0 ;
my $iden ; 

my @snumbers ;
my @qnumbers ;

my $SUB = {};
my $QUERY = {};

my $origlength ;
my $hasintrons=1;
my $TOTALDIFF = 0 ;
my $identityline ;

my $IDXS = {};
my $STRINGSSORTED = {};
while(<$ifh>){
	chomp ;
	if(/Length=/ && ! defined $origlength){
		($origlength) = /Length=(\d+)/ ;
		print "origlength = $origlength\n" if($verbose);
	}


     if(/^\s*Identities/){
	    print "$_\n" if($verbose);
	 	$DONOW = 1 ;
	 	$CNT++;
		$identityline = $_ if(!defined $identityline);
		
		if(defined $qstart){
			my @llll = split " ", $iden ;
			my $IDEN = $llll[3];
			$IDEN =~ s/\(//;
			$IDEN =~ s/\)//;
			$IDEN =~ s/\%//;
			$IDEN =~ s/\,//;
		    my $SSSS = "Query ($qstart-> $qend) Subject ($sstart-> $send) $iden ";
			$STRINGSSORTED->{$qstart} = $SSSS ;
			if($IDEN ne 100){
			   my @IDXS ;
			   push @IDXS,$sstart;
			   push @IDXS,$send;
			   push @IDXS,$qstart;
			   push @IDXS,$qend;
			   push @IDXS,$IDEN;
			   $IDXS->{$sstart} = \@IDXS ;
			}
			if($sstart > $send){
				die "$sstart > $send - should not be here" if(!defined $findifreverse);
				system("touch $trs.isrev");
				die "Please comp the subject - it is reverse\n";
			}
			else{
				if(defined $findifreverse){
					die "Trial run - this is in the fwd direction\n";
				}
			}
			push @snumbers, $sstart ;
			push @snumbers, $send ;

			foreach my $i ($sstart..$send){
				$SUB->{$i} = 1 ;
			}
			foreach my $i ($qstart..$qend){
				$QUERY->{$i} = 1 ;
			}
			push @qnumbers, $qstart ;
			push @qnumbers, $qend ;

			# just doing the first ....may need to change
            my $diff = abs($qstart - $qend)+1;
			$TOTALDIFF = $diff + $TOTALDIFF ;
			#print $ofhscaffoldlog "$diff \n";
			if($diff eq $origlength){
				print "There are no introns - exact match\n";
				$dointron = 0 ;
				$hasintrons = 0 ;
				last ;
			}

		}


		undef $qstart ;
		undef $qend ;
		undef $sstart ;
		undef $send ;

		$iden = $_;
	}
	if($DONOW && /Query/){
		my ($one,$two,$three,$four) = split ;
		$qstart = $two if(!defined $qstart);
		$qend = $four ;
		print "qstart $qstart $_ \n" if($verbose);
	}
	if($DONOW && /Sbjct/){
		my ($one,$two,$three,$four) = split ;
		$sstart = $two if(!defined $sstart);
		$send = $four ;
		print "sstart $sstart $_ \n" if($verbose);
	}
}
close($ifh);



#if(!@snumbers){
if(1){
		if(defined $qstart){
			my @llll = split " ", $iden ;
			my $IDEN = $llll[3];
			$IDEN =~ s/\(//;
			$IDEN =~ s/\)//;
			$IDEN =~ s/\%//;
			$IDEN =~ s/\,//;
		    my $SSSS = "Query ($qstart-> $qend) Subject ($sstart-> $send) $iden ";
			$STRINGSSORTED->{$qstart} = $SSSS ;
			if($IDEN ne 100){
			my @IDXS ;
			push @IDXS,$sstart;
			push @IDXS,$send;
			push @IDXS,$qstart;
			push @IDXS,$qend;
			push @IDXS,$IDEN;
			$IDXS->{$sstart} = \@IDXS ;
			}
			if($sstart > $send){
				system("touch $trs.isrev");
				die "Please comp the subject - it is reverse\n";
			}

			push @snumbers, $sstart ;
			push @snumbers, $send ;

			foreach my $i ($sstart..$send){
				$SUB->{$i} = 1 ;
			}
			foreach my $i ($qstart..$qend){
				$QUERY->{$i} = 1 ;
			}
			push @qnumbers, $qstart ;
			push @qnumbers, $qend ;

			# just doing the first ....may need to change
            my $diff = abs($qstart - $qend)+1;
			$TOTALDIFF = $diff + $TOTALDIFF ;
			#print $ofhscaffoldlog "$diff \n";
			if($diff eq $origlength){
				print "There are no introns - exact match\n";
				$hasintrons = 0 ;
				$dointron = 0 ;
			}

		}
}

foreach my $i (sort {$a <=> $b} keys %{$STRINGSSORTED}){
	my $v = $STRINGSSORTED->{$i} ;
	print $ofhscaffoldlog "$v \n";
}


my $MARKED = {};
my $ERRORS = 0 ;
if(1){

	my $strS = "";
	my $strQ = "";
	if(1){
	   ($strS) = util_readfasta($subject);
	   my $lenS = length($strS);
	   ($strQ) = util_readfasta("../../FASTADIR/$trs.ALL.1.fasta");
	   my $lenQ = length($strQ);
	}

    my @sl = sort {$a <=> $b}  (keys %{$IDXS});
	
	my @l ;
	foreach my $i (@sl){
		my @x = @{$IDXS->{$i}};
		push @l, @x ;
	}


	if(!@l){
		print "Hundred percent match $nameofsubject $trs \n";
	}
	else{

	   my $totalSlen = length($strS);

	   my @SLICES ; 
       my $sliceFirst;
	   my $BEGIN ; 
	   while(@l){
	       my $sS = shift @l	;
	       my $sE = shift @l	;
	       my $qS = shift @l	;
	       my $qE = shift @l	;
	       my $IDEN = shift @l	;
   
		   foreach my $xxx ($sS...$sE){
			   if(exists $MARKED->{$xxx}){
			      $ERRORS = 1;
			   }
			   else{
				   $MARKED->{$xxx} = 1 ;
			   }
		   }
   
		   if(!defined $sliceFirst){
			   $BEGIN = 1 ;
		   }
           $sliceFirst = util_extractSliceFromFastaString($strS,$BEGIN,$sS-1);
		   push @SLICES,$sliceFirst;
           my $sliceQ = util_extractSliceFromFastaString($strQ,$qS,$qE);
	       push @SLICES,$sliceQ;
		   $BEGIN = $sE+1;
   
   
	   }
       my $sliceLast = util_extractSliceFromFastaString($strS,$BEGIN,$totalSlen);
	   push @SLICES,$sliceLast;
   
       if(!$ERRORS){
		   my $final = "";
		   foreach my $str (@SLICES){
		       $str =~ s/\\n//;
			   $final = $final . $str ;
		   }
	       my $fname = "NEWFASTA/$nameofsubject.$trs.ALL.1.fasta";
		   print "Writing to $fname \n";
	       my $ofhnewfasta = util_write($fname);
	       print $ofhnewfasta ">$nameofsubject\n";
	       print $ofhnewfasta "$final\n";
	   }
   }

}

if(defined $exact){
	$buffer = 0 ;
}

my @sl = sort {$a <=> $b}  @snumbers ;
my @ql = sort {$a <=> $b}  @qnumbers ;

my $N = @sl - 1 ;
my $sS = $sl[0] - $buffer ;
my $sE = $sl[$N] + $buffer ;

my $qS = $ql[0] ;
my $qE = $ql[$N];

my $qDIFF = abs($qS - $qE)+1;
my $sDIFF = abs($sS - $sE)+1;
$identityline = $hasintrons ? "": $identityline ;
print $ofhlog "$trs $subject $hasintrons $qS $qE $qDIFF $sS $sE $sDIFF $identityline \n";


my $first = 0 ; 
my $begin ;
my $introncnt = 0 ;
if($dointron){
    print "Getting intron slices\n";
    foreach my $x ($sS..$sE){
	    if(! exists $SUB->{$x} ){
		    $begin = $x -1 if(!defined $begin);
		    $first = 1 ;
	    }
	    else{
		    if($first){
		       my $end = $x - 1;
		       $introncnt++;
		       my $beginplus = $begin + $bufferintron ;
		       my $endminus = $end - $bufferintron ;
               system(" extractslicefromfasta.pl -out INTRONS/$trs.intron.$introncnt.start -inf $subject -sta $begin -end $beginplus");
               system(" extractslicefromfasta.pl -out INTRONS/$trs.intron.$introncnt.end -inf $subject -sta $endminus -end $end");
		       $first = 0 ;
		       undef $begin;
		    }
		    
	    }
    }
}

if(1){
   print "Getting real slices\n";
   system ("extractslicefromfasta.pl -out $outfile.subject -inf $subject -sta $sS -end $sE");
   system ("extractslicefromfasta.pl -out $outfile.query -inf $query -sta $qS -end $qE");


   my $DIFF = (($origlength - $TOTALDIFF));
   my $PERCENT = (($DIFF)/$origlength)*100 ;
   print "origlength = $origlength, TOTALDIFF = $TOTALDIFF, DIFF = $DIFF, PERCENT = $PERCENT \n";
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
