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
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
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
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;


my $info = {};
my $first = 1 ;
my $Q ; 
my $X ; 
my $ONE ; 
my $TWO ; 
my $CCC = 0 ; 
while(<$ifh>){
     next if(/^\s*$/);
	 chop; 
     if(/---/ && defined $ONE){
	 	PrintNewQ($CCC++);
	 }
     if(/^\s*Q/){
	    s/\s*..//;
	 	$Q = $_ ;
	 }
     elsif(/^\s*X/){
	    s/\s*..//;
	 	$X = $_ ;
	 }
     elsif(/^\s*(A|B|C|D|E)/){
	    s/\s*..//;
	 	if(!defined $ONE){
	 	    $ONE = $_ ;
		}
		else{
			$TWO = $_ ;
		}
	 }
}

if(defined $ONE){
	PrintNewQ($CCC++);
}
close($ifh);

foreach my $k (keys %{$info}){
     $k = uc($k);
	 my $val = $info->{$k} ; 
	 print $ofh "$k\n";
}

sub PrintNewQ{
	my ($counter) = @_ ; 
    print "sub Q$counter { \n";
	print "my \$level = \" 9 \"; \n";
	print "my \$subject = \"MATHS\"; \n";
    print "\tmy \@WA ; \n";
    print "\tmy \$question = \"$Q \";\n";
    print "\tmy \$answer = \"$X\" ; \n";
    print "\tmy \$hint = \"\"; \n";
    print "\tpush \@WA , \" $ONE \" ;  \n";
    print "\tpush \@WA , \" $TWO\" ;  \n";
    print "\tmy (\$ANSLIST,\$correctindex) = QP_InsertCorrectAnswer(\$answer,\@WA); \n";
    print "\t(\$question,\$answer,\$hint,\$ANSLIST,\$correctindex); \n";
    print "}\n";
	undef $ONE;
	undef $TWO;
}

print STDERR "Output written in $outfile\n";


chmod 0777, $outfile ;

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
