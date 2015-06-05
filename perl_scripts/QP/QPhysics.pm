
package QPhysics;
use Carp ;
use MyUtils ;
use MyQPaper ;
use POSIX ;
require Exporter;
use Algorithm::Combinatorics qw(combinations) ;
use Math::NumberCruncher;
use Math::MatrixReal;  # Not required for pure vector math as above
use Math::Geometry ; 
use Math::VectorReal qw(:all);  # Include O X Y Z axis constant vectors
#use Math::Trig;
#use Math::Trig ':radial';
no warnings 'redefine';
my $EPSILON = 0.01;

  local $SIG{__WARN__} = sub {};

@ISA = qw(Exporter);
@EXPORT = qw( 
MOMENTUMDEF
NEWTON2LAW

MOMENTUM1_CHANGEINVEL
MOMENTUM2_CHANGEINMOM
MOMENTUM3_RATEOFCHANGE
MOMENTUM4_RIFLERECOIL
MOMENTUM5_BULLETWOOD

MOTION1_AIRPLANEDIST
MOTION2_FORCEHOWLONG
);

use strict ;
use FileHandle ;
use Getopt::Long;

sub MOMENTUMDEF  {
	my @WA ; 
	my $question ="Momentum is defined as -  " ; 
	my $answer = "velocity * mass  " ; 
	my $hint = "";
	push @WA , "mass * acceleration  " ;
	push @WA , "velocity * acceleration " ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub NEWTON2LAW  {
	my @WA ; 
	my $question ="According to Newton's second law, force is proportional to  " ; 
	my $answer = "rate of change of momentum  " ; 
	my $hint = "";
	push @WA , "rate of change of speed  " ;
	push @WA , "rate of change of acceleration  " ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}


sub MOMENTUM1_CHANGEINVEL  {
	my $level = " 9 "; 
    my ($x,$y,$z) =  QP_GetNRandomNumbersBelowValue(3,15);
	my $subject = " Physics "; 
	my @WA ; 
	my $question ="The initial velocity of a ball of mass $x kg is $y metre/s and the final velocity is $z metre/sec. What is the change in velocity  " ; 
	my $ANS = $z -$y;
	my $answer = $ANS . " metre/sec ";
	my $WA1 = $ANS + 2 ;
	my $WA2 = $ANS - 2 ;
	my $hint = "  ";
	push @WA , $WA1 . " metre/sec ";
	push @WA , $WA2 . " metre/sec ";
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub MOMENTUM2_CHANGEINMOM  {
	my $level = " 9 "; 
    my ($x,$y,$z) =  QP_GetNRandomNumbersBelowValue(3,15);
	my $subject = " Physics "; 
	my @WA ; 
	my $question ="The initial velocity of a ball of mass $x kg is $y metre/s and the final velocity is $z metre/sec. What is the change in momentum?  " ; 
	my $ANS = ($z -$y)*$x;
	my $answer = $ANS . " (kg*metre)/sec ";
	my $WA1 = ($z -$y)*$x + 2 ;
	my $WA2 = $ANS*$x - 2 ;
	my $hint = "  ";
	push @WA , $WA1 . " (kg*metre)/sec ";
	push @WA , $WA2 . " (kg*metre)/sec ";
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub MOMENTUM3_RATEOFCHANGE{
	my $level = " 9 "; 
    my ($x,$y,$z,$t) =  QP_GetNRandomNumbersBelowValue(4,15);
	my $subject = " Physics "; 
	my @WA ; 
	my $question ="The initial velocity of a ball of mass $x kg is $y metre/s, and after $t secs the final velocity is $z metre/sec. What is the rate of change in momentum?" ; 
	my $ANS = util_format_float((($z -$y)*$x)/$t,1);
	my $answer = $ANS . " metre/(sec*sec) ";
	my $WA1 = ($z -$y)*$x + 2 ;
	my $WA2 = $ANS*$x - 2 ;
	my $hint = "  ";
	push @WA , $WA1 . " (kg*metre)/sec ";
	push @WA , $WA2 . " (kg*metre)/sec ";
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}


sub MOMENTUM4_RIFLERECOIL  {
	my $level = " 9"; 
	my $subject = " Physics "; 
    my ($M1,$M2,$V) =  QP_GetNRandomNumbersBelowValue(3,15);
	my $a = util_format_float(($M2*$V)/$M1,1);
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  - 5; 
	my $question ="A $M1 kg rifle fires a $M2 g bullet at a velocity of $V m/s find the recoil velocity of the rifle? " ; 
	my $answer = " $a  " ; 
	my $hint = " Use the conservation of momentum concept  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);

}


sub MOMENTUM5_BULLETWOOD  {
	my $level = " 9"; 
	my $subject = " Physics "; 
    my ($M1,$M2,$V) =  QP_GetNRandomNumbersBelowValue(3,15);
	my $a = util_format_float(($M2*$V)/$M1,1);
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  - 5; 
	my $question ="A $M2 g bullet at a travelling at a velocity of $V m/s hits a stationary plank of wood of mass $M1 g and gets embedded in the plank. What is the velocity with which the combined mass of the wood and the bullet moves ? " ; 
	my $answer = " $a  " ; 
	my $hint = " Use the conservation of momentum concept  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub MOTION1_AIRPLANEDIST  {
	my $level = " 9"; 
	my $subject = " Physics "; 
    my ($A,$T) =  QP_GetNRandomNumbersBelowValue(2,15);
    my $a = 1/2*$A*$T*$T ; 
	my @WA ; 
	my $WA1 = $a - 5 ;
	my $WA2 = $a + 5 ;
	my $question ="An airplane accelerates down a runway at $A m/s*s for $T s until it finally lifts off the ground. Determine the distance traveled before takeoff " ; 
	my $answer = " $a  " ; 
	my $hint = " Use equation s = ut + 1/2*a*t*t. Remember to check unit.  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub MOTION2_FORCEHOWLONG  {
	my $level = " 9"; 
	my $subject = " Physics "; 
    my ($V,$F,$M) =  QP_GetNRandomNumbersBelowValue(3,15);
my $A = $F/$M ; 
my $T = $V/$A; 
my $answer = util_format_float($T,1) ; 
	my @WA ; 
	my $WA1 = $answer + 3 ; 
	my $WA2 =  $answer - 3; 
	my $question ="In order to gain a velocity of $V m/s how long should a force of $F N be exerted on a body of mass $M kg that is initially at rest? " ; 
	my $hint = " Find the acceleration first from the formuala: F = m*a  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}
