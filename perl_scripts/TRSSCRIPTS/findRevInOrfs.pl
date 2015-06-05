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
my ($infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
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

system ("touch logclose") if(! -e "logclose");
my $ofhlog = util_append("logclose");

usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();

my $ignoretable = {};
if(defined $ignorefile){
   my @lll= util_read_list_sentences($ignorefile);
   map { s/\s*//g ; $ignoretable->{$_} = 1 ; } @lll ;
}


my @list= util_read_list_sentences($listfile);
my $list = {};
map { s/\s*//g ; $list->{$_} = 1 ; } @list ;

foreach my $i (@list){

    my $ifh = util_read("FASTADIR_ORFNEW/$i.ALL.1.fasta");
    my $info = {};
    while(<$ifh>){
	    if(/>/){
		    my ($ORF) = (/ORF_(\d+)/);
		    $ORF  = "ORF_$ORF";
		    #print "ORF $ORF\n";
            my $ifh1 = util_read("ORF/$i.orf");
			while(<$ifh1>){
				if(/$ORF/){
					#print "$ORF $_\n";
					if(/REVERSE/){
	                    my ($subject,$isrev) = util_Blastout_isrev("BLASTOUT_2WGS/$i.blast.nt");
						if(!$isrev){
						    print $ofh "$i 0 1 \n";
						}
						else{
						   print $ofh "$i 0 0\n";
						}
					}
					else{
	                    my ($subject,$isrev) = util_Blastout_isrev("BLASTOUT_2WGS/$i.blast.nt");
						if($isrev){
						    print $ofh "$i 1 1 \n";
						}
						else{
						    print $ofh "$i 1 0 \n";
						}

					}
					last ;
				}
			}
			close($ifh1);
	    }
    }
	close($ifh);
}


chmod 0777, $outfile ;

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
