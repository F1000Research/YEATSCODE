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
my ($trs,$infile,$p1,$ignoreband,$cutoff,$p2,$outfile,$which_tech,$ignorefile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "trs=s"=>\$trs ,
            "p2=s"=>\$p2 ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "cutoff=i"=>\$cutoff ,
            "ignoreband=i"=>\$ignoreband ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -trs ") if(!defined $trs);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);


my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);

#my @NAMES = qw(CE CI CK EM FL HC HL HP HU IF LE LM LY PK PL PT RT SE TZ VB);
my @NAMES = qw (CE  CI  CK  EM  FL  HC  HL  HP  HU  IF  LE  LM  LY  PK  PL  PT  RT  SE  TZ  VB);



my $STR2PRINT= "";

my $max = 0 ;
while(<$ifh>){
	 my (@l) = split ; 
	 my ($name) = shift @l ;

	 next if($name ne $trs);
	 my $N = @NAMES -1 ;

	 foreach my $idx (0..$N){
	 	my $code = $NAMES[$idx];
		my $a = shift @l ;
		#print $ofh "$idx $a\n";
		print $ofh "($code,$a)\n";
		$STR2PRINT = $STR2PRINT . "($code,$a)\n";
		$max = $a if($a > $max);
	

	 }
}

$max = $max + 100 ;

$trs =~ s/_/\\_/g;
PrintTexFile();

sub PrintTexFile{
print  << "ENDOFUSAGE" ; 
\\documentclass[border=10pt]{standalone}
\\usepackage{pgfplots}
\\pgfplotsset{width=7cm,compat=1.8}
\\usepackage{pgfplotstable}
\\renewcommand*{\\familydefault}{\\sfdefault}
\\usepackage{sfmath}
\\begin{document}
\\begin{tikzpicture}
  \\centering
  \\begin{axis}[
        ybar, axis on top,
        title={$trs},
        height=8cm, width=15.5cm,
        bar width=0.4cm,
        %ymajorgrids, tick align=inside,
        major grid style={draw=white},
        %enlarge y limits={value=.1,upper},
        ymin=0, ymax=$max,
        axis x line*=bottom,
        axis y line*=right,
        y axis line style={opacity=0},
        tickwidth=0pt,
        enlarge x limits=true,
        legend style={
            at={(0.5,-0.2)},
            anchor=north,
            legend columns=-1,
            /tikz/every even column/.append style={column sep=0.5cm}
        },
        ylabel={Count},
      symbolic x coords={CE,CI,CK,EM,FL,HC,HL,HP,HU,IF,LE,LM,LY,PK,PL,PT,RT,SE,TZ,VB},
       xtick=data,
       nodes near coords={
        \\pgfmathprintnumber[precision=0]{\\pgfplotspointmeta}
       }
    ]
	\\addplot [draw=none, fill=blue!30] coordinates {
ENDOFUSAGE
print  $STR2PRINT ;
print  << "ENDOFUSAGE1" ; 
	  };

    \\legend{}
  \\end{axis}
  \\end{tikzpicture}
\\end{document}

ENDOFUSAGE1
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
