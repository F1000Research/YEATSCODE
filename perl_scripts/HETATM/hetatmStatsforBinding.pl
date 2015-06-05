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
my ($infile,$outfile,$special,$listfile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "special"=>\$special ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);

system ("touch logclose") if(! -e "logclose");
my $ofhlog = util_append("logclose");

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
#usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;

my  ($seconds, $microseconds) = gettimeofday;

#my $pdb = "$PDBDIR/$protein.pdb";
#my $pdb1 = new PDB();
#$pdb1->ReadPDB($pdb);



my $info = {};
my $done = {};
my $CNT = 0 ;
my $GROUPS = {} ; 
while(<$ifh>){
     next if(/^\s*$/);
	 my (@l) = split ; 
	 my $PROT = shift @l ;
	 my $HETATM = shift @l ;
	 my $SIZE = shift @l ;
     my $sorted = {};

	 my $ATOMSTR =  "";
	 foreach my $r (@l){
	 	my (@l) = split "/", $r ;
		my $N = @l - 1; 
		my $x = $l[0];

		my $RESATOM = $l[1];
		my $HETATOM = $l[2];

		my ($firstA) = ($RESATOM =~ /(.)/);
		my ($firstB) = ($HETATOM =~ /(.)/);

		my $PAIR = $firstA. $firstB;
		$ATOMSTR = $ATOMSTR . $PAIR ;

		my $XXX = $x . "/". $l[1];
        my $num = $x ;
	    $num =~ s/...//;

	 	$sorted->{$num} = $XXX ;

		my $dist = $l[$N];
		next if($dist ge 3.8);
		#print "$num $dist lllll\n";
		if(! exists $info->{$x}){
			$info->{$x} = 1 ;
		}
		else{
			$info->{$x} = $info->{$x} + 1 ;
		}
	 	
	 }
	 #print "$ATOMSTR \n";

	 my $final = "";
	 foreach my $i (sort {$a <=> $b} keys %{$sorted}){
	 	my $val = $sorted->{$i};
	 	$final = $final . "XXX" . $val ;
	 }
	 $final =~ s/GLU203/GLU205/;
	 $final =~ s/GLU204/GLU206/;
	 $final =~ s/SER631/SER630/;
	 $final =~ s/TYR548/TYR547/;
	 $final =~ s/TYR663/TYR662/;
	 $final =~ s/TYR632/TYR631/;
	 if(! exists $done->{$final}){
	 	 $CNT++;
		 my $MOTIFNAME = "2OQVA$CNT";

	     $GROUPS->{$ATOMSTR} = [] if(!defined $GROUPS->{$ATOMSTR});
	     push @{$GROUPS->{$ATOMSTR}}, $PROT ;

		 my $fname = "${protein}MOTIF.$CNT";
	     #print "$PROT $MOTIFNAME $final\n";
		 if(0){
		 my $OFH = util_write($fname);
	     print $OFH "$final\n";
		 close($OFH);
		 }
		 print $ofh "preparePremonConfigs.pl -outf lll -lis $PROT.4.1.table.in -pr 2OQVA -pol 1 -size 4 -anndir ANNOTATE.4\n";
		 print $ofh "cp -f $PROT.4.clasp.in ANNOTATE.4/\n";
		 print $ofh "createPremoninput.pl -out ANNOTATE.4/2OQVA$CNT.4.1.premon.in -con \$CONFIGGRP -li ANNOTATE.4/$PROT.4.clasp.in -pr 2OQVA\n";
		 print $ofh "ln -s ANNOTATE.4/2OQVA$CNT.4.1.premon.in ANNOTATE.4/$PROT.4.1.premon.in \n";

		 #print $ofh "mv -f ANNOTATE.4/2OQVA.4.1.premon.in ANNOTATE.4/2OQVA$CNT.4.1.premon.in \n";
		 foreach my $IDX (1..4){
		      print $ofh "mv -f ANNOTATE.4/2OQVA.4.1.premon.in.config$IDX ANNOTATE.4/2OQVA$CNT.4.1.premon.in.config$IDX \n";
		      print $ofh "mv -f ANNOTATE.4/2OQVA.4.1.premon.in.aalist$IDX ANNOTATE.4/2OQVA$CNT.4.1.premon.in.aalist$IDX \n";
		 }

		 $done->{$final} = 1 ;
	 }

}

my $GROUPFH = util_write("groups");

foreach my $k (sort keys %{$GROUPS}){
	my @l = @{$GROUPS->{$k}};
	print $GROUPFH "$k @l \n";
}

foreach my $k (sort keys %{$info}){
     $k = uc($k);

     my $num = $k ;
	 $num =~ s/...//;
	 my $val = $info->{$k} ; 
	 print $ofh "$num $k $val\n";
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
