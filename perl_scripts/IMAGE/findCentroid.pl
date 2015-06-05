#!/usr/bin/perl -w 
use MyMagick;
use MyUtils;
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use ConfigPDB;
use MyGeom;
use PDB;

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($infile,$outfile,$display,$which_tech,$debug,$listfile,$protein);
my ($color,@expressions,$from,$to);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "color=s"=>\$color ,
            "to=s"=>\$to ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "display"=>\$display ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "debug"=>\$debug ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

usage( "Need to give a color -option -color  ") if(!defined $color);
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);

my $ofh = util_write("centroid.out");


my @list= util_read_list_sentences($listfile);


my @MX ;
my @MY ;


my $fx ;
my $fy ;

my $prevX;
my $prevY;

my $currdist = -1 ;
my $cumudist = -1 ;
foreach my $infile (@list){
   my $ifh = util_read($infile);
   $outfile  ="$infile.centroid.png";
   my $image = new MyMagick($infile);
   
   $image->DebugInfo() if(defined $debug);
   
   my $w = $image->GetWidth();
   my $h = $image->GetHeight();
   
   print "W = $w H = $h\n";
   
   my $obj = $image->GetObj();
   
   my (@colorpoints) = $image->FindCoordsForColor($color);
   

   my @X ;
   my @Y ;
   while(@colorpoints){
	   my $x = shift @colorpoints ;
	   my $y = shift @colorpoints ;

	   push @X, $x ;
	   push @Y, $y ;
   }


   my ($meanX,$sdX) = util_GetMeanSD(\@X);
   my ($meanY,$sdY) = util_GetMeanSD(\@Y);

	   push @MX, $meanX ;
	   push @MY, $meanY ;


   if(!defined $fx){
   	$fx = $meanX ;
   	$fy = $meanY ;

   }
   else{

       $currdist = geom_Distance_2D ($meanX,$meanY,$prevX,$prevY);
       $cumudist = geom_Distance_2D ($meanX,$meanY,$fx,$fy);
   }
   $prevX = $meanX ;
   $prevY = $meanY ;

   $image->SetColorToPixelByColorName($meanX,$meanY,"white");
   print $ofh "meanX=$meanX meanY=$meanY currdist=$currdist cumudist=$cumudist\n";

   $image->Write($outfile,$display);
   close($ifh);
}





sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ;
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}

