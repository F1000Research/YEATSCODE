#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($all,$infile,$outfile,$or,$silent,$groupinfo);
my ($DIR,$listfile);
my $howmany = 600000 ; 
my $cutofflength = 0 ; 
my @types = (); 
my $protein ;
my @motifs = (); 
GetOptions(
            "all"=>\$all ,
            "groupinfo"=>\$groupinfo ,
            "silent"=>\$silent ,
            "infile=s"=>\$infile ,
            "dir=s"=>\$DIR ,
            "listfile=s"=>\$listfile ,
            "howmany=i"=>\$howmany ,
            "or=i"=>\$or ,
            "cutofflength=i"=>\$cutofflength ,
            "type=s"=>\@types,
            "protein=s"=>\$protein,
            "motif=s"=>\@motifs,
            "outfile=s"=>\$outfile 
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);

print STDERR "Info: parsing file $infile - might take some time\n";
my ($info,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEWTMP($infile,1,1);


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
sub util_parsePDBSEQRESNEWTMP{
	my ($infile,$all,$writedata) = @_ ; 
    my $ifh = util_read($infile);
    my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();

	my $seenNm ; 
	my $info  = {}; 
	my $infoSeq2PDB = {} ; 
	my $mapChainedName2Name = {} ; 
	my $uniqueseqs = {} ; 
	
	my $nm ;
	my $CNT = 0 ;
    while(<$ifh>){
         next if(/^\s*$/);
	     if(/>/){
	        ($nm) = />(.*)\s*/;
			$nm = uc($nm);
			die if(!defined $nm); 
			die "Repeat $nm " if(exists $info->{$nm});
			$CNT++;
			$info->{$nm} = "" ;
			next ;
		}


		## add to seq
		die if(/>/);
		if(!defined $_){
		    die "$nm $CNT " ;
		}
		chomp;
		$info->{$nm} = $info->{$nm} . $_ ;

	 }

	 foreach my $k (keys %{$info}){
	 	my $seq = $info->{$k};
		my $length = length($seq);
		#print "$length $k \n";
		$infoSeq2PDB->{$seq} = [] if(!defined $infoSeq2PDB->{$seq}) ;
		push @{$infoSeq2PDB->{$seq}}, $k ;
	 }

	 my $N = (keys %{$infoSeq2PDB});
	 print "There were $CNT seqeunces to start with, and there are $N unique \n";

	 if(defined $writedata && $writedata){
         my $ofh = util_write("list.unique");

		 my $SORTING = {};
         foreach my $k (keys %{$infoSeq2PDB}){
		 	my $len = length($k);
			$SORTING->{$k} = $len ;
		 }


		 my @sl = sort { $SORTING->{$b} <=> $SORTING->{$a}} (keys %{$SORTING});
		 my $CNTwrite = 0 ;
         foreach my $k (@sl){
	        my @pdbswithseq  = sort @{$infoSeq2PDB->{$k}} ;
		 	my $len = length($k);
			if($CNTwrite < 3 ){
			#if(1){
			   #print "$len \n";
		       foreach my $pdb (@pdbswithseq){
				   ($pdb) = split " ", $pdb ;
			       print "KKKKKKKK $pdb\n";
				   exit ;
			       my $fastanm = "$FASTADIR/$pdb.ALL.1.fasta";
				   $CNTwrite++;
			       my $FH = util_write($fastanm);
                  print $FH "\>$pdb;\n";
                  print $FH "$k\n";
		          close($FH);
       
		       }
		       my $v = $pdbswithseq[0];
	 	       print $ofh "$v\n";
			}
	     }
		 print "Wrote $CNTwrite finally\n";
             close($ofh);


         $ofh = util_write("list.allpdbschained");
         foreach my $k (keys %{$mapChainedName2Name}){
	 	    print $ofh "$k\n";
	     }
         close($ofh);

         $ofh = util_write("list.allpdbs");
         foreach my $k (keys %{$info}){
	 	    $k = uc($k);
	 	    print $ofh "$k\n";
	     }
         close($ofh);
	 }


	 return ($info,$infoSeq2PDB,$mapChainedName2Name) ;
}
