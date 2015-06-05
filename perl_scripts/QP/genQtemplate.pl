#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use PDB;
use ConfigPDB;

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

my $info = {};
while(<$ifh>){
     next if(/^\s*$/);
	 chomp;
	 my ($nm,$q,$a,$hint,$WA1,$WA2,$subject,$level) = split "#", $_ ;
	 Template($nm,$q,$a,$hint,$WA1,$WA2,$subject,$level) ;
}
close($ifh);

system ("cat $outfile");
	
sub Template{
	my ($nm,$q,$a,$hint,$WA1,$WA2,$subject,$level) = @_ ; 
	 print $ofh "sub $nm {\n";
	 print $ofh "\tmy \$level = \"$level\"; \n"; 
     print $ofh "\tmy (\$n1,\$n2,\$n3,\$n4) =  QP_GetNRandomNumbersBelowValue(4,20); \n";
     print $ofh "\tmy (\$n1,\$n2,\$n3,\$n4) =  QP_GetNRandomNumbersBelowValueBothSign(4,20); \n";
	 print $ofh "\tmy \$subject = \"$subject\"; \n"; 
	 print $ofh "\tmy \$a =   \n";
	 print $ofh "\tmy \@WA ; \n";
	 print $ofh "\tmy \$WA1 = $a + 5 ; \n";
	 print $ofh "\tmy \$WA2 =  $a - 5; \n";
     print $ofh "\tmy \$question =\"$q \" ; \n";
     print $ofh "\tmy \$answer = \"$a \" ; \n";
	 print $ofh "\tmy \$hint = \"$hint \";\n";
	 print $ofh "\tpush \@WA , \$WA1  ;\n";
	 print $ofh "\tpush \@WA , \$WA2  ;\n";
	 print $ofh "\tmy (\$ANSLIST,\$correctindex) = QP_InsertCorrectAnswer(\$answer,\@WA);\n";
     print $ofh "\t(\$question,\$answer,\$hint,\$ANSLIST,\$correctindex);\n";
     print $ofh "}\n";
}
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
