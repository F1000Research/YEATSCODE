#!/usr/bin/perl -w 

use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;

use MyUtils;
use MyQPaper;
use QPhysics;
use POSIX qw(floor);

my $commandline = util_get_cmdline("",\@ARGV) ;

my ($easy, $nodecimals, $infile,$class,$subset,$qpaper,$howmany,$design,$multiple,$answer,$clocked);
$howmany = 0 ;
$multiple = 1 ;
GetOptions(
            "howmany=s"=>\$howmany ,
            "infile=s"=>\$infile ,
            "clocked:s"=>\$clocked ,
            "design=s"=>\$design ,
            "answer=s"=>\$answer ,
            "multiple=i"=>\$multiple ,
            "class=i"=>\$class ,
            "subset=i"=>\$subset ,
            "nodecimals"=>\$nodecimals ,
            "easy"=>\$easy ,
            "qpaper=s"=>\$qpaper );

die "Dont recognize command line arg @ARGV " if(@ARGV);

usage( "Need to give a infile, => option -infile ") if(!defined $infile);

my $id = floor(100000*rand());
$qpaper = "qpaper.$id";
$answer = "answer.$id";
while(-e $qpaper){
   $id = floor(100000*rand());
   $qpaper = "qpaper.$id";
   $answer = "answer.$id";
}

my @PLACES = qw( house bridge road garden );
my @PPL = qw(boys girls men women);
my @THINGS = qw( oranges apples mangoes books );
my @alpha = qw( a b c d e f g h i j k l m n o p q r s t u v w x y z );
my @ALPHA = qw( A B C D E F G H I J K L M N O P Q R S T U V W X Y Z );
my @names = qw(Ravi Rishi Sanjay Amit Pooja Ritu Amrita Amba Shalini);
my @numbersabove80 = qw(80 90 100 110 120 130 140 150 160 170 );
my @SUBJECTS = qw(Maths Science History Geography Life-science);
my @SUMFORTRIANGLE = qw(18 20 10 9 36 90 45 );

my $ofh = util_write($qpaper);
my $ofhDB = util_write("DB");
my $answerfh = util_write($answer);
print STDERR "Writing qpaper file $qpaper and answers in $answer\n";
unlink ("QP");
unlink ("AP");
system("ln -s $qpaper QP");
system("ln -s $answer AP");

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}

my @qtypes ;
my $leveltable = {};
my $subjecttable = {};
my $ifh = util_read($infile);
while(<$ifh>){
     next if(/^\s*$/);
     next if(/^\s*#/);
     chop ;
	 my ($nm,$level,$subject) = split ; 
	 push @qtypes, $nm ;
	 $leveltable->{$nm} = $level ;
	 $subjecttable->{$nm} = $subject ;
}


my $cnt = 1 ;
while($multiple){
    #print $ofh      " QUESTIONS id=$id  $multiple \\\\\n";
    #print $ofh      " ================================ \\\\\n";

    print $ofh      "\\section{Question set}\n";
    #print $ofh      " ================================ \\\\\n";

    print $answerfh " ANSWERS id=$id $multiple \\\\\n";
    print $answerfh      " ================================ \\\\\n";
    
    $multiple--;


	my @qtypestmp;
    if(defined $subset){
       my $SMALL = util_pick_n_random_from_list(\@qtypes,$subset);
       @qtypestmp = @{$SMALL};
    }
	else{
       @qtypestmp = @qtypes ;
	}

        foreach my $qtype (@qtypestmp){

            my ($q,$a,$hint,$MQA,$correctindex) ;
		    my $sub = \&{$qtype}; 
			my $level = $leveltable->{$qtype} or die ; 
			my $subject = $subjecttable->{$qtype} or die ; 
			next if(defined $class && $level != $class);
			print STDERR "Processing $qtype \n";
		    ($q,$a,$hint,$MQA,$correctindex) = $sub->();
			die "$qtype not proper for answer" if (!defined $a);
			die "$qtype not proper for question" if (!defined $q);
			$hint = "No hints, sorry..." if($hint eq "");

            print $ofh      " $cnt ) {\\bf $q}\\\\ \n";
            print $answerfh " $cnt ) {\\bf $a}\\\\ \n";
			
			print $ofhDB "myArray[$cnt]=\" ";
			print $ofhDB "$q#";
			print $ofhDB "$a#";
			my $NMQA = @{$MQA} ; 
			die "Need 3 and only 3" if($NMQA ne 3);
			foreach my $mqa (@{$MQA}){
				print $ofhDB "$mqa #";
			}
			print $ofhDB "$correctindex#";
			print $ofhDB "$hint#";
			print $ofhDB "$qtype#";
			print $ofhDB "$level#";
			print $ofhDB "$subject#";

			print $ofhDB "\" ; \n ";

            $cnt++;
       }
}


sub ADDMULTIPLE{
    my ($div) = @_ ;
	my $sum = 0 ; 
	my $str = "";
	while($sum >= 0){
		$sum = 0 ; 
		$str = "";
        my $n1 = 5 ; 
        while($n1){
           my $celprob = util_round(rand());
           my $p = floor(20*rand()) + 1 ; 
	       $str = $celprob ? "$str + $p " : "$str - $p"; # do this first
	       $p = $celprob ? $p : $p * (-1) ; 
	       $sum = $sum + $p ;
	       $n1--;
	    }
	}


	### The common parts 
	my @WA ; 
    my $answer = $sum;
    my $question = "Find the sum of: $str ";
	my $hint = "Group positive and negative numbers together";
	push @WA , $answer + floor(10*rand()) + 1 ;
	push @WA , $answer - floor(10*rand()) - 1 ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($sum,@WA);
    ($question,$answer,$hint,$ANSLIST,$correctindex);
}


sub DIV{
	return MULT(1);
}

sub MULT{
    my ($div) = @_ ;
    my $n1 = 20;
    my $n2 = 10 ;
    if(defined $div){
        $nodecimals = 1 ;
    }

    $n2 = 1 if(defined $nodecimals);
    my ($p,$q);
   
    while(!($p = floor($n1*rand()+1)/$n2)){}
    while(!($q = floor($n1*rand()+1)/$n2)){}
    
    my $a = $p * $q ;
    my $oper = " X " ;
    if(defined $div){
         $oper = " / " ;
         my $t = $a ; 
         $a = $p ; 
         $p = $t ; 
    }
    undef $nodecimals ;

	### The common parts 
	my @WA ; 
    my $answer = $a;
    my $question = "$p $oper $q\ = ?";
	my $hint = defined $div? "See if there are common factors": "";
	push @WA , $answer + floor(10*rand()) + 1 ;
	push @WA , $answer - floor(10*rand()) - 1 ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($a,@WA);
    ($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub TEMPERATURE{
    my $n1 = 100 ;
    my $cel = floor($n1*rand());
    my $fahren = (9*$cel)/5 + 32 ; 
    
    my $celprob = util_round(rand());
	my ($q, $a);
    if($celprob){
        $q = "Convert $cel degree celsius into fahrenheit" ;
		$a = $fahren ;
    }
    else {
        $q = "Convert $fahren degree fahrenheit into celsius" ;
		$a = $cel ;
    }

	### The common parts 
	my @WA ; 
    my $answer = $a;
    my $question = $q ;
	my $hint = "Equation for conversion of celsius to fahrenheit is: C/5 = (F -32)/9";
	push @WA , $answer + floor(10*rand()) + 1 ;
	push @WA , $answer - floor(10*rand()) - 1 ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($a,@WA);
    ($question,$answer,$hint,$ANSLIST,$correctindex);
}


sub FRACTIONS{
    my $n1 = 15 ;
	my ($a,$b,$c,$d) = QP_GetNRandomNumbersBelowValue(4,$n1);
    my $answer = ($a/$b + $c/$d) ; 
    $answer = QP_round2place($answer);
    my $q = $a . "/" . $b . " + " . $c . "/" . $d ; 

	my @WA ; 
    my $question = $q ;
	my $hint = "Find the LCM of $b and $d";
	push @WA , $answer + floor(5*rand()) + 1 ;
	push @WA , $answer - floor(5*rand()) - 1 ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
    ($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub ASCDESC{
    my $n1 = 15 ;
	my ($a,$b,$c,$d,$e,$f) = QP_GetNRandomNumbersBelowValue(6,$n1);

    my $p = $a . "/" . $b  ;
    my $q = $c . "/" . $d  ;
    my $r = $e . "/" . $f  ;
 
    my $x = QP_round2place($a/$b * 1.00);
    my $y = QP_round2place($c/$d * 1.00);
    my $z = QP_round2place($e/$f * 1.00);
    if($x > $y ){
        ($x,$y) = ($y,$x) ;
        ($p,$q) = ($q,$p) ;
    }
    if($y > $z ){
        ($y,$z) = ($z,$y) ;
        ($q,$r) = ($r,$q) ;
    }
   
    my $celprob = util_round(rand());
    my $answer = $celprob ?  " $p, $q, $r " : " $r, $q, $p";
    my $WA1 = $celprob ?  " $r, $p, $q " : " $q, $p, $r";
    my $WA2 = $celprob ?  " $r, $q, $p " : " $q, $r, $p";
    my $what = $celprob ? " ascending " : " descending ";
    my $expr = $a . "/" . $b . " , " . $c . "/" . $d . " , " . $e . "/" . $f  ; 

	my @WA ; 
    my $question = " Arrange in $what order $expr ";
	my $hint = "Like fractions, make the denominator the same by finding LCM";
	push @WA , $WA1 ;
	push @WA , $WA2;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
    ($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub CIRCLEAreaCircum{
    my $n1 = 9 ;
    my $r = (floor($n1*rand()) + 1)*7;
    my $a = (22/7) * $r * $r ;
    my $c = (22/7) * 2 * $r ;
    my $celprob = util_round(rand());
	my ($q);
    if($celprob){
        $q = " Radius of a circle is $r cm. Find its area";
	    $a = $a ; 
	}
	else{
        $q = " Radius of a circle is $r cm. Find its circumference" ;
	    $a = $c ; 
	}
	my @WA ; 
    my $question = $q ; 
    my $answer = $a ; 
	my $hint = "Area = pi *r *r ; Circumference = 2 * pi * r ; pi = 22/7 ; ";
	push @WA , $answer + floor(5*rand()) + 10 ;
	push @WA , $answer - floor(5*rand()) - 10 ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
    ($question,$answer,$hint,$ANSLIST,$correctindex);
}


sub UNITARYWORK{
    my $n1 = 15 ;
    my $n2 = 15 ;
    my ($days,$boys) =  QP_GetNRandomNumbersBelowValue(2,15);
    my ($X1,$X2) =  QP_GetNRandomNumbersBelowValue(2,10);
    my $days1 = $days * ($X1+1);
    my $boys1 =  $boys * ($X2+1);

	my $WHOM = util_pick_random_from_list(\@PPL);
    my $what = util_pick_random_from_list(\@PLACES);
    
    
    my $celprob = util_round(rand());
	my ($q,$a) ; 
    if($celprob){
        $q = "$boys $WHOM can build a $what in $days days. How many days will $boys1 $WHOM build it?"; 
        $a = QP_round2place(($boys*$days)/$boys1);
    }
    else {
        $q = "$boys $WHOM can build a $what in $days1 days. How many $WHOM will make it in $days days ?"; 
        $a = QP_round2place(($boys*$days1)/$days);
    }

	my @WA ; 
    my $question = $q ; 
    my $answer = $a ; 
	my $hint = $celprob ?  "Find out how many days will 1 $WHOM take to make a $what. Think whether 1 $WHOM will take more days or less days - and divide or multiply accordingly" : "Find out in 1 day how may $WHOM will be required to make the $what. Will it be more or less $WHOM required. multiply or divide accordingly";
	push @WA , $answer + floor(5*rand()) + 10 ;
	push @WA , $answer - floor(5*rand()) - 10 ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
    ($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub UNITARYCOST{
    my ($cost,$num) =  QP_GetNRandomNumbersBelowValue(2,15);
    my ($X1,$X2) =  QP_GetNRandomNumbersBelowValue(2,10);
    my $cost1 = $cost* ($X1+1);
    my $num1 =  $num * ($X2+1);

    my $what = util_pick_random_from_list(\@THINGS);
    
    
    my $celprob = util_round(rand());
	my ($q,$a) ; 
    if($celprob){
        $q = "$num $what costs Rs $cost . What is the cost of $num1 $what?";
        $a = QP_round2place(($cost/$num)*$num1);
    }
    else {
        $q = "$num $what costs Rs $cost . How many $what can you buy for Rs $cost1?";
        $a = QP_round2place(($num/$cost)*$cost1);
    }
	my @WA ; 
    my $question = $q ; 
    my $answer = $a ; 
	my $hint = $celprob ?  "Find out how what is the cost of 1 $what" : "Find out how many $what can you buy for Rs 1.";
	push @WA , $answer + floor(5*rand()) + 10 ;
	push @WA , $answer - floor(5*rand()) - 10 ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
    ($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub SQUAREROOTIRRATIONAL{


    my ($num,$num1) =  QP_GetNRandomNumbersBelowValue(2,10);
	my $sum = $num + $num1 ; 
	my $prod  = 4* $num * $num1 ; 
    
    my $celprob = util_round(rand());
	my ($q,$a);

	my $WA1;
	my $WA2;
    my ($r1,$r2,$r3,$r4) =  QP_GetNRandomNumbersBelowValue(4,10);
    if($celprob){
        $q = "Find the square root of $sum + sqrt($prod)";
        $a = "sqrt($num) + sqrt($num1)";
        $WA1 = "sqrt($r1) + sqrt($r3)";
        $WA2 = "sqrt($r2) + sqrt($r4)";
    }
    else {
        $q = "Find the square root of $sum - sqrt($prod)";
        $a = "sqrt($num) - sqrt($num1), or sqrt($num1) - sqrt($num)";
        $WA1 = "sqrt($r1) - sqrt($r4), or sqrt($r2) - sqrt($r4)";
        $WA2 = "sqrt($r2) - sqrt($r3), or sqrt($r3) - sqrt($r1)";
    }
	my @WA ; 
    my $question = $q ; 
    my $answer = $a ; 
	my $hint = "Remember (a+b)^2 = a^2 + 2ab + b^2. Try to factor the expresion under the root into a 4. Thus, a 2 will come out of the square root.";
	push @WA , $WA1 ;
	push @WA , $WA2 ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
    ($question,$answer,$hint,$ANSLIST,$correctindex);
	
}


### Not really used now 
sub SUMOFSQUARES{
    my ($n1,$n2) =  QP_GetNRandomNumbersBelowValue(2,10);
	my $square1 = $n1 *$n1 ; 
	my $square2 = $n2 *$n2 ; 
	my $sum = $square1 + $square2 ;


    my $q = "Express $sum as the sum of two squares";
    my $a = "$n1 x $n1 + $n2 x $n2 " ;
    return $q , "$q = $a";
	
}

sub tmpMULTEXPR{
    my @l1 ; 
    my @l2 ; 
    my @l ; 
    
    push @l1 , genTerm("a");
    push @l1 , genTerm("b");
    
    push @l2 , genTerm("a");
    push @l2 , genTerm("b");
    
    push @l , \@l1 ;
    push @l , \@l2;
    
    my $q = printExpr(@l);
    my $a = multExpr(@l);
	return ($q,$a);
}


sub MULTEXPR{
    my ($q,$a) = tmpMULTEXPR();
    my ($q1,$WA1) = tmpMULTEXPR();
    my ($q2,$WA2) = tmpMULTEXPR();

	my @WA ; 
    my $question = "Multiply $q";
    my $answer = $a ; 
	my $hint = "";
	push @WA , $WA1 ;
	push @WA , $WA2 ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
    ($question,$answer,$hint,$ANSLIST,$correctindex);

}

sub FACTOREXPR{
    my ($a,$q) = tmpMULTEXPR();
    my ($WA1,$q1) = tmpMULTEXPR();
    my ($WA2,$q2) = tmpMULTEXPR();

	my @WA ; 
    my $question = "Factorise $q";
    my $answer = $a ; 
	my $hint = "";
	push @WA , $WA1 ;
	push @WA , $WA2 ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
    ($question,$answer,$hint,$ANSLIST,$correctindex);

}
sub FACTOREXPROLDFORcheck{
    my @l1 ; 
    my @l2 ; 
    my @l ; 
    
    push @l1 , genTerm("a");
    push @l1 , genTerm("b");
    
    push @l2 , genTerm("a");
    push @l2 , genTerm("b");
    
    push @l , \@l1 ;
    push @l , \@l2;
    
    my $a = printExpr(@l);
    my $q = multExpr(@l);
    return "Factorise $q" , "$q = $a";
}

sub genTerm{
	my ($a) = @_ ; 
    my $celprob = util_round(rand());
    my $n1 =  floor(5*rand())+1;
	my @term ; 
	if($celprob){
		 push @term , "+";
		 push @term , $n1 ; 
		 push @term , $a ; 
	}
	else{
		 push @term , "-";
		 push @term , $n1 ; 
		 push @term , $a ; 
	}
	return \@term ; 
}

sub printExpr{
	my (@l) = @_ ; 
	my $finalexpr = "";
	foreach my $l (@l){
		my @exprs = @{$l} ; 
        my $expr = "( ";
		foreach my $e (@exprs){
		    my $term = join "",  @{$e};
			$expr = $expr . $term ; 
		}
		$expr = $expr . ") " ; 
		$finalexpr = $finalexpr . $expr ; 
	}
	return $finalexpr ;
}

sub multExpr{
	my (@l) = @_ ; 
	my $l1 = shift (@l);
	my @exprs1 = @{$l1} ; 
	my $l2 = shift (@l);
	my @exprs2 = @{$l2} ; 

	my $done = {};
	foreach my $e1 (@exprs1){
	    my ($sign1, $coeff1, $nm1) = @{$e1};
		
	    foreach my $e2 (@exprs2){
			 my ($sign2, $coeff2, $nm2) = @{$e2};
			 my $SIGN = $sign1 eq $sign2 ? "+" : "-";
			 my $COEFF = $coeff1 * $coeff2 ;
			 my $NM =  $nm1 lt $nm2 ? $nm1 . $nm2 :$nm2 . $nm1  ; 
			 $done->{$NM} = 0 if(!defined $done->{$NM});

			 $done->{$NM} = $sign1 eq $sign2 ?  $COEFF + $done->{$NM} : - $COEFF + $done->{$NM}  ; 
		}
	}
	my @XXX ; 
	foreach my $k (keys %{$done}){
		my $v = $done->{$k} ;
		my $expr = "$v$k" ;
		$expr = "+" . $expr if(!($expr =~ /^\s*-/));
		##$finalexpr = $finalexpr . " " .  $expr ; 
		push @XXX, $expr ;
	}

	## just to set the order correct
	my $AA = shift @XXX ;
	my $BB = shift @XXX ;
	push @XXX, $AA ; 
	push @XXX, $BB ; 

    my $finalexpr = "";
	foreach my $expr (@XXX){
		$finalexpr = $finalexpr . " " .  $expr ; 
	}
	
	return $finalexpr ;
}





sub EQNWhatnumberaddedgives  {
	my $level = " 6 "; 
	my ($x,$y,$n3,$n4) =  QP_GetNRandomNumbersBelowValueBothSign(4,50); 
	my $subject = "MATHS"; 
	my $a = $y - $x ; 
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  - 5; 
	my $question =" What number added to $x gives $y? " ; 
	my $answer = " $a  " ; 
	my $hint = " HINT  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}
sub EQNwhatnumbermultiplied  {
	my $level = " 6 "; 
	my ($x,$y,$n3,$n4) =  QP_GetNRandomNumbersBelowValueBothSign(4,50); 
	my $subject = "MATHS"; 
	my $a =  $y/$x ;
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  - 5; 
	my $question =" What number multiplied by $x gives $y? " ; 
	my $answer = " $a  " ; 
	my $hint = " HINT  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}
sub EQNwhatnumberdivided  {
	my $level = " 6 "; 
	my $subject = "MATHS"; 
	my ($x,$y,$n3,$n4) =  QP_GetNRandomNumbersBelowValueBothSign(4,50); 
	my $a =  $x * $y ; 
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  - 5; 
	my $question =" What number divided by $x gives $y? " ; 
	my $answer = " $a  " ; 
	my $hint = " HINT  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub EQN_MONEY{
	my ($nm1,$nm2) = @{util_pick_n_random_from_list(\@names,2)} ; 

    my ($M1,$M2)         =  QP_GetNRandomNumbersBelowValue(2,15);
    my ($n1,$n2,$n3,$n4) =  QP_GetNRandomNumbersBelowValue(4,15);

	my $V1 = $n1 * $M1 ; 
	my $V2 = $n2 * $M2 ; 
	my $S1 = $V1 + $V2 ; 

	my $V3 = $n3 * $M1 ; 
	my $V4 = $n4 * $M2 ; 
	my $S2 = $V3 + $V4 ; 

    my $q = "$nm1 and $nm2 have some money. If you multiply $nm1 's money by $n1 and $nm2 's money by $n2 and add them up, you get Rs $S1. If you multiply $nm1 's money by $n3 and $nm2 's money by $n4 and add them up, you get Rs $S2. How much money does each person have? " ;
 
    my ($r1,$r2,$r3,$r4,$r5) =  QP_GetNRandomNumbersBelowValue(5,15);
	my $a = "$nm1 = $M1, $nm2 = $M2";
	my $WA1 = "$nm1 = $r1, $nm2 = $r2";
	my $WA2 = "$nm1 = $r3, $nm2 = $r4";

	my @WA ; 
    my $question = $q ; 
    my $answer = $a ; 
	my $hint = "Set up the simulataneous equation using x and y. Eliminate one variable by making the coefficient same, and then subtracting";
	push @WA , $WA1 ;
	push @WA , $WA2 ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
    ($question,$answer,$hint,$ANSLIST,$correctindex);

}

sub FACTOREXPR_CUBE{
	my ($nm1,$nm2) = @{util_pick_n_random_from_list(\@alpha,2)} ; 

    my ($M1,$M2)         =  QP_GetNRandomNumbersBelowValue(2,6);
	
	my $cube1 = $M1 * $M1 * $M1 ;
	my $cube2 = $M2 * $M2 * $M2 ;


	my $sq1 = $M1 * $M1 ;
	my $sq2 = $M2 * $M2 ;


	#print "$nm1*$nm2 $M1 $M2\n";

    my $random = util_round(rand()) ;
    my $celprob = ! $random ? "-" : "+";
    my $celprobNot = ! $random ? "+" : "-" ;
    my $ql = "Factorize \$ ${cube1}$nm1\^\{3\} $celprob ${cube2}$nm2\^\{3\} \$\n";
    my $q = "Factorize ${cube1}*$nm1*$nm1*$nm1 $celprob ${cube2}*$nm2*$nm2*$nm2";
 
 	my $MULT = $M1 * $M2 ; 
	my $al = "\$ (${M1}$nm1 $celprob ${M2}$nm2) (${sq1}$nm1\^\{2\} $celprobNot ${MULT}${nm1}$nm2 + ${sq2}$nm2\^\{2\}) \$  ";
	my $a = "(${M1}$nm1 $celprob ${M2}$nm2) (${sq1}$nm1*$nm1 $celprobNot ${MULT}${nm1}$nm2 + ${sq2}$nm2*$nm2)  ";

    my ($r1,$r2,$r3,$r4,$r5) =  QP_GetNRandomNumbersBelowValue(5,15);
	my $WA1 = "(${r1}$nm1*$nm1 $celprobNot ${r2}$nm2*$nm2) (${r3}$nm1*$nm1 $celprob ${r4}${nm1}$nm2 + ${r5}$nm2*$nm2)  ";
    ($r1,$r2,$r3,$r4,$r5) =  QP_GetNRandomNumbersBelowValue(5,15);
	my $WA2 = "(${r1}$nm1*$nm1 $celprobNot ${r2}$nm2*$nm2) (${r3}$nm1*$nm1 $celprob ${r4}${nm1}$nm2 + ${r5}$nm2*$nm2)  ";

	my @WA ; 
    my $question = $q ; 
    my $answer = $a ; 
	my $hint = "a*a*a + b*b*b = (a - b)(a*a + ab + b*b) and a*a*a - b*b*b = (a + b)(a*a - ab + b*b)";
	push @WA , $WA1 ;
	push @WA , $WA2 ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
    ($question,$answer,$hint,$ANSLIST,$correctindex);

}


sub PERCENT1  {
	my $level = " 7"; 
	my $subject = " MATHS "; 
	my $y = 100 ;
    while($y eq 100) {
        $y = util_pick_random_from_list(\@numbersabove80);
    }
    my $x = floor($y*rand()) ;
	my $subject = util_pick_random_from_list(\@SUBJECTS);
	my $a =   
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  - 5; 
	my $question ="If you have scored $x marks in $subject out of $y marks, what is the percentage of marks scored by you " ; 
	my $answer = " $a  " ; 
	my $hint = " Use unitary method - find out how much you scored in if the exam was of 1 marks, and then if the exam was of 100 marks. That is the percentage  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub TRIANGLEANGLE1  {
	my $level = " 8"; 
	my $subject = "MATHS"; 
    my $finalsum = util_pick_random_from_list(\@SUMFORTRIANGLE);
	my ($x, $y, $z) = QP_Find_N_numbers_that_addupto($finalsum,3);
	my $common = 180/$finalsum ; 
	my $A = $common * $x ; 
	my $B = $common * $y ; 
	my $C = $common * $z ; 
	my $BWA = $common * $y -1 ; 
	my $CWA = $common * $z + 1 ; 
	my $a = " $A $B $C ";
	my @WA ; 
	my $WA1 =  " $A $BWA $C ";
	my $WA2 =  " $A $B $CWA ";
	my $question ="The angles of a triangle are in the ratio of $x, $y and $z. Find their absolute values " ; 
	my $answer = " $a  " ; 
	my $hint = "   ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}



sub GEOMDistancePoints  {
	my $level = " 9 "; 
	my $subject = "MATHS"; 
    my ($x1,$x2,$y1,$y2)         =  QP_GetNRandomNumbersBelowValueBothSign(4,6);

	my $X = ($x1 - $x2) ;
	my $Y = ($y1 - $y2); 
    my $a = $X*$X + $Y *$Y;
	$a = util_format_float($a);
	my @WA ; 
	my $A1 =  $a  + 5 ; 
	my $A2 =   $a  - 5; 
	my $WA1 =  "sqrt($A1)";
	my $WA2 =  "sqrt($A2)";
	my $question =" What is the distance between the points ($x1,$y1) and ($x2, $y2)? " ; 
	my $answer = " sqrt ($a)  " ; 
	my $hint = " HINT  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}


sub GEOMDividingInternally  {
	my $level = " 9 "; 
	my ($x1,$x2,$y1,$y2) =  QP_GetNRandomNumbersBelowValueBothSign(4,20); 
	my ($n1,$n2,$n3,$n4) =  QP_GetNRandomNumbersBelowValue(4,20); 
	my $subject = "MATHS"; 
	my $X =  util_format_float(($n1*$x2 + $n2*$x1)/($n1+$n2),1);
	my $Y =  util_format_float(($n1*$y2 + $n2*$y1)/($n1+$n2),1);
    my $XX = $X + 5 ;
    my $YY = $Y + 5 ;
	my $a = "($X,$Y)";
	my @WA ; 
	my $WA1 =  "($X,$YY)";
	my $WA2 =  "($XX,$Y)";
	my $question =" Find the coordinate of the point that divides the points ($x1,$y1) and ($x2, $y2) internally in the ratio of $n1:$n2? " ; 
	my $answer = " $a  " ; 
	my $hint = " HINT  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub RATIOProportionalitConst  {
	my $level = " 8"; 
	my ($n1,$n2,$n3,$n4) =  QP_GetNRandomNumbersBelowValueBothSign(4,20); 
	my $subject = "MATHS"; 
    my $k = $n2/$n1;
	my $a = $k * $n3 ; 
	$a = util_format_float($a);
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  - 5; 
	my $question ="x is directly proportional to y. When x is $n1, y is $n2.  Find x when y is $n3.  " ; 
	my $answer = " $a  " ; 
	my $hint = " HINT  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}


sub RATIOSimpleEquation  {
	my $level = " 8"; 
	my ($n1,$n2,$n3,$n4) =  QP_GetNRandomNumbersBelowValueBothSign(4,20); 
	my $subject = "MATHS"; 
	my $a =  ($n1 * $n3 )/ $n2 ; 
	$a = util_format_float($a);
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  - 5; 
	my $question =" Find the unknown value in the proportion:  $n1 : x = $n2 : $n3.  " ; 
	my $answer = " $a  " ; 
	my $hint = " HINT  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub RATIOCutpipeinRatio  {
	my $level = " 8"; 
    my ($n1,$n2,$n3,$n4) =  QP_GetNRandomNumbersBelowValueSorted(4,50);
	my ($x1,$x2) =  QP_GetNRandomNumbersBelowValue(4,20); 
	my $subject = "MATHS"; 
    my $k = ($n4)/($x1 + $x2);
	my $P = $x1 * $k ;
	$P = util_format_float($P);
	my $Q = $x2 * $k ;
	$Q = util_format_float($P);
	my $PW = $x1 * $k + 5 ;
	my $QW = $x2 * $k -5 ;
	my $a =   "$P and $Q";
	my @WA ; 
	my $WA1 =  "$PW and $Q";
	my $WA2 =  "$P  and $QW";
	my $question ="One piece of pipe $n4 meters long is to be cut into two pieces, with the lengths of the pieces being in a $x1 : $x2 ratio. What are the lengths of the pieces? " ; 
	my $answer = " $a  " ; 
	my $hint = " HINT  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}




sub GEOMAreaofRightAngleTriangle  {
	my $level = " 8"; 
	my ($n1,$n2,$n3,$n4) =  QP_GetNRandomNumbersBelowValue(4,20); 
	my $subject = "MATHS"; 
	my $a = 1/2 * $n1 * $n2 ; 
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  - 5; 
	my $question ="What is the area of a right angled triangle with base = $n1 cm and height = $n2 cm ?  " ; 
	my $answer = " $a  " ; 
	my $hint = " HINT  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}
sub GEOMAreaofIsocelesTriangle  {
	my $level = " 8"; 
	my ($n1,$n2,$n3,$n4) =  QP_GetNRandomNumbersBelowValue(4,20); 
	my $subject = "MATHS"; 
	my $a =  $n1 * $n2 ; 
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  - 5; 
	my $question ="What is the area of a isosceles triangle with base = $n1 cm and height = $n2 cm ?  " ; 
	my $answer = " $a  " ; 
	my $hint = " HINT  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}
sub GEOMCircleFindRadiiFromChord  {
	my $level = " 8"; 
	my ($n1,$n2,$n3,$n4) =  QP_GetNRandomNumbersBelowValue(4,20); 
	my $subject = "MATHS"; 
	my $a = util_format_float(sqrt ( $n1 * $n1 + $n2 * $n2),1);
	my @WA ; 
	my $WA1 =  $a  + 1 ; 
	my $WA2 =   $a  - 2; 
	my $question ="What is the radius of a circle with a chord = $n1 cm and distance of the chord from centre = $n2? cm  " ; 
	my $answer = " $a  " ; 
	my $hint = " HINT  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}



sub EQN_FATHERAGE  {
	my $level = " 9"; 
	my ($n1) =  QP_GetNRandomNumbersBetweenValues(1,3,6); 
	my ($n2) =  QP_GetNRandomNumbersBetweenValues(1,8,15); 
	my ($Sage) =  QP_GetNRandomNumbersBetweenValues(1,10,20); 
	my ($Fage) =  $n1 * $Sage ;
	my $Ratio = util_format_float(($Fage + $n2)/($Sage + $n2),1);

	my $subject = "MATHS"; 
	my $a =  $Fage ;
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  - 5; 
	my $question ="A father is $n1 times as old as his son, after $n2 years, the father will be $Ratio as old as his son at that time, therefore, the fathers current age is: " ; 
	my $answer = " $a  " ; 
	my $hint = " HINT  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub EQN_TWONUMBERRATIO  {
	my $level = " 9"; 
	my ($n1) =  QP_GetNRandomNumbersBetweenValues(1,3,6); 
	my ($n2) =  QP_GetNRandomNumbersBetweenValues(1,8,15); 
	my ($Sage) =  QP_GetNRandomNumbersBetweenValues(1,10,20); 
	my ($Fage) =  $n1 * $Sage ;
	my $Ratio = util_format_float(($Fage + $n2)/($Sage + $n2),1);

	my $subject = "MATHS"; 
	my $a =  $Fage ;
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  - 5; 
	#my $question ="A father is $n1 times as old as his son, after $n2 years, the father will be $Ratio as old as his son at that time, therefore, the fathers current age is: " ; 
	my $question ="Two numbers are in the ratio $n1. If you add $n2 to both, their ratio becomes $Ratio. Find the numbers : " ; 

	my $answer = " $a  " ; 
	my $hint = " HINT  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub DENSITY{
    my ($n1,$n2,$n3,$n4) =  QP_GetNRandomNumbersBelowValue(4,20);
    my $question = "The volume of $n1 Kg of material M1 is $n2 litres, and the volume of $n3 KG of material M2 is $n4 litres. Which material is heavier?";
	my $heavier = $n1/$n2 >= $n3/$n4 ? "M1" : "M2";
	my $notheavier = $n1/$n2 < $n3/$n4 ? "M1" : "M2";
    my $a = $heavier ; 

	my $level = " 8"; 
	my $subject = "SCIENCE"; 
	my @WA ; 
	my $WA1 =  " $notheavier";
	my $WA2 =   "neither";
	my $answer = " $a  " ; 
	my $hint = " HINT  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub MOLE{
    my ($n1,$n2) =  QP_GetNRandomNumbersBelowValue(2,30);
    my $question = "The molecular mass of a compound C is $n1. How many atoms are there in $n2 grams of the compound?";
	my $a = "$n2/$n1 * 6.023 * \$10^{23}\$" ;

	my $level = " 9"; 
	my $subject = "SCIENCE"; 
	my @WA ; 
	my $WA1 = "$n1/$n2 * 6.023 * \$10^{23}\$" ;
	my $WA2 = "$n1 * 6.023 * \$10^{23}\$" ;
	my $answer = " $a  " ; 
	my $hint = " HINT  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub ATOMICNUMBER{
    my ($n1,$n2) =  QP_GetNRandomNumbersBelowValueSorted(2,30);
    my $question = "The atomic number of an element E is $n1, and its mass number is $n2. How many neutrons does the element have?";
	my $a = $n2 - $n1 ;

	my $level = " 9"; 
	my $subject = "SCIENCE"; 
	my @WA ; 
	my $WA1 = "$n2";
	my $WA2 = "$n1";
	my $answer = " $a  " ; 
	my $hint = " HINT  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub SET_INTERSECTION{
	my $U = util_pick_n_random_from_list(\@alpha,10);
	my $A = util_pick_n_random_from_list($U,5);
	my $B = util_pick_n_random_from_list($U,6);

	my $strU = join ", ", @{$U};
	my $strA = join ", ", @{$A};
	my $strB = join ", ",  @{$B};

    my $celprob = util_round(rand());
	my $what ; 
	my $strC ; 
    if($celprob){
		$what = "intersection";
	    my $tabA = util_make_table($A);
	    my @intersection = ();
	    foreach my $p (@{$B}){
	      push @intersection,$p if(exists $tabA->{$p});	
	    }
	    $strC = join",",  @intersection ;
	}
	else{
		$what = "union";
	    my $tabA = util_make_table($A);
	    foreach my $p (@{$B}){
	      $tabA->{$p} = 1 ;
	    }
	    my @intersection = (keys %{$tabA});
	    $strC = join",",  @intersection ;
	}

    my $question = "Universal set = ( $strU ) set A = ( $strA ), set B = ( $strB ) . Find A $what B by drawing venn diagrams";
	my $a = " ( $strC ) " ;
	my $WA1 = " ($strA )";
	my $WA2 = " ($strB )";

	my $level = " 9"; 
	my $subject = "MATHS"; 
	my @WA ; 
	my $answer = " $a  " ; 
	my $hint = " HINT  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);

}

sub APGenerate{
	my ($n) =  QP_GetNRandomNumbersBetweenValues(1,20,30); 
	my ($a) =  QP_GetNRandomNumbersBetweenValues(1,1,5); 
	my ($d) =  QP_GetNRandomNumbersBelowValueBothSign(1,5); 

	my @sums ; 
	my @terms ; 
	my $sum = 0 ;
	foreach my $i (1..$n){
		my $term = $a + ($i -1)* $d ;
		$sum = $sum + $term ; 
		push @terms, $term ;
		push @sums, $sum ;
	}

	return ($n,$d,$a,\@sums, \@terms);
	
}


sub AP_FINDNTHTERM  {
	my $level = " 10"; 
	my ($n1) =  QP_GetNRandomNumbersBetweenValues(1,10,20); 
	my ($n,$d,$first,$sums, $terms) = APGenerate();
	my @sums = @{$sums};
	my @terms = @{$terms};
	my $subject = "MATHS"; 
	my $a =  $first + ($n1 -1)*$d;
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  + 12; 
	my $seq  = "$terms[0], $terms[1], $terms[2]....";
	my $question ="Find the $n1 th term of the AP given by $seq." ; 
	my $answer = " $a  " ; 
	my $hint = "  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub AP_SOLVEEQN  {
	my $level = " 10"; 
	my ($n,$d,$first,$sums, $terms) = APGenerate();
	my ($n1,$n2) =  QP_GetNRandomNumbersBetweenValues(2,10,19); 
	my @sums = @{$sums};
	my @terms = @{$terms};

	my $N1 = $terms[$n1];
	my $N2 = $terms[$n2];
	my $subject = "MATHS"; 
	my $a =  "first = $first, difference = $d ";
	my @WA ; 
	my $AAA =  $first  + 5 ; 
	my $DDD =   $d  + 12; 
	my $WA1 =  "first = $AAA, difference = $d ";
	my $WA2 =  "first = $first, difference = $DDD ";
	my $seq  = "$terms[0], $terms[1], $terms[2]....";
	my $question ="The $n1 th term of an AP is $N1, and the $n2 th term of the AP is $N2. Find the first number and the difference. " ; 
	my $answer = " $a  " ; 
	my $hint = "  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub AP_LOAN  {
	my $level = " 10"; 
	my ($n,$d,$first,$sums, $terms) = APGenerate();
	my $FINALSUM = -1 ;
	while($FINALSUM < 0){
	   ($n,$d,$first,$sums, $terms) = APGenerate();
	   my @sums = @{$sums};
	   $FINALSUM = $sums[$n-1];
	}
	my @sums = @{$sums};
	my @terms = @{$terms};
	my $subject = "MATHS"; 
	my $a =  $first ;
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  + 12; 
	my $seq  = "$terms[0], $terms[1], $terms[2]....\n";
	my $what = $d > 0 ? "increasing": "decreasing";
	my $question ="A man takes a loan of $FINALSUM, and repays it by $what the payment by $d every month, and finishing the payment in $n months. How much should he pay in the first month?";
	my $answer = " $a  " ; 
	my $hint = "  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub AP_SUMN  {
	my $level = " 10"; 
	my ($n,$d,$first,$sums, $terms) = APGenerate();
	my @sums = @{$sums};
	my $FINALSUM = $sums[$n-1];
	my @terms = @{$terms};
	my $subject = "MATHS"; 
	my $a =  $FINALSUM ;
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  + 12; 
	my $seq  = "$terms[0], $terms[1], $terms[2]....\n";
	my $what = $d > 0 ? "increaing": "decreasing";
	my $question ="Find Sum(n) for n=$n of an AP with first term = $first and difference = $d.";
	my $answer = " $a  " ; 
	my $hint = "  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub GPGenerate{
	my ($n) =  QP_GetNRandomNumbersBetweenValues(1,10,15); 
	my ($a) =  QP_GetNRandomNumbersBetweenValues(1,1,5); 
	my ($d) =  QP_GetNRandomNumbersBelowValueBothSign(1,5); 

	my @sums ; 
	my @terms ; 
	my $sum = 0 ;
	my $prevtem = $a ;
	foreach my $i (1..$n){
		my $term = $i == 1 ? $prevtem : $prevtem * $d ; 
		$prevtem = $term ;
		$sum = $sum + $term ; 
		push @terms, $term ;
		push @sums, $sum ;
	}

	return ($n,$d,$a,\@sums, \@terms);
	
}

sub GP_FINDNTHTERM  {
	my $level = " 10"; 
	my ($n1) =  QP_GetNRandomNumbersBetweenValues(1,4,6); 
	my ($n,$d,$first,$sums, $terms) = GPGenerate();
	my @sums = @{$sums};
	my @terms = @{$terms};
	my $subject = "MATHS"; 
	my $a =  $terms[$n1 -1];
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  + 12; 
	my $seq  = "$terms[0], $terms[1], $terms[2]....";
	my $question ="Find the $n1 th term of the GP given by $seq." ; 
	my $answer = " $a  " ; 
	my $hint = "  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub GP_SUMN  {
	my $level = " 10"; 
	my ($n1) =  QP_GetNRandomNumbersBetweenValues(1,4,6); 
	my ($n,$d,$first,$sums, $terms) = GPGenerate();
	my @sums = @{$sums};
	my @terms = @{$terms};
	my $subject = "MATHS"; 
	my $a =  $sums[$n1 -1];
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  + 12; 
	my $seq  = "$terms[0], $terms[1], $terms[2]....";
	my $question ="Find the Sum($n1) for n=$n1 of the GP given by $seq." ; 
	my $answer = " $a  " ; 
	my $hint = "  ";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}


sub _ElecGenerate{
	my ($q) =  QP_GetNRandomNumbersBetweenValues(1,5,10); 
	my ($t) =  QP_GetNRandomNumbersBetweenValues(1,5,10); 
	my ($i) =  util_format_float($q/$t,1) ;
	my ($r) =  QP_GetNRandomNumbersBetweenValues(1,5,10); 
	my $v = $i * $r ; 
	my $w = $v * $q ; 
	my $p = util_format_float($w/$t,1) ; 
	return ($q,$t,$i,$r,$v,$w,$p); 
}

sub ELEC_I_from_QT  {
	my $level = " 10"; 
	my ($q,$t,$i,$r,$v,$w,$p) = _ElecGenerate();
	my $question ="Find the current when $q Coulomb charge passes through a cross section in $t secs\n";
	my $a =  $i ;
	$a = util_format_float($a,1);
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  + 12; 
	my $answer = " $a  " ; 
	my $hint = "I = q/t";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub ELEC_Q_from_IT  {
	my $level = " 10"; 
	my ($q,$t,$i,$r,$v,$w,$p) = _ElecGenerate();
	my $question ="Find the charge flowing when $i Amp currnet passes through a cross section for $t secs\n";
	my $a =  $q ;
	$a = util_format_float($a,1);
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  + 12; 
	my $answer = " $a  " ; 
	my $hint = "I = q/t";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}


sub ELEC_V_from_IR  {
	my $level = " 10"; 
	my ($q,$t,$i,$r,$v,$w,$p) = _ElecGenerate();
	my $question ="Find the potential difference when $i Amp current passes through conductor with $r ohm resistance. \n";
	my $a =  $v ;
	$a = util_format_float($a,1);
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  + 12; 
	my $answer = " $a  " ; 
	my $hint = "V = IR";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub ELEC_V_from_WQ  {
	my $level = " 10"; 
	my ($q,$t,$i,$r,$v,$w,$p) = _ElecGenerate();
	my $question ="Find the potential difference when $w joules of work is done in moving $q Coulombs\n";
	my $a =  $v ;
	$a = util_format_float($a,1);
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  + 12; 
	my $answer = " $a  " ; 
	my $hint = "V = W/Q";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub ELEC_P_from_WT{
	my $level = " 10"; 
	my ($q,$t,$i,$r,$v,$w,$p) = _ElecGenerate();
	my $question ="Find the power when $w joules of work is done in $t secs";
	my $a =  $p ;
	$a = util_format_float($a,1);
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  + 12; 
	my $answer = " $a  " ; 
	my $hint = "P = W/t";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub ELEC_P_from_IR{
	my $level = " 10"; 
	my ($q,$t,$i,$r,$v,$w,$p) = _ElecGenerate();
	my $question ="Find the power when $i current flows through a conductor with resistance $r";
	my $a =  $p ;
	$a = util_format_float($a,1);
	my @WA ; 
	my $WA1 =  $a  + 5 ; 
	my $WA2 =   $a  + 12; 
	my $answer = " $a  " ; 
	my $hint = "P = I*R*R";
	push @WA , $WA1  ;
	push @WA , $WA2  ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
	($question,$answer,$hint,$ANSLIST,$correctindex);
}

sub QUADROOTS_FACTOR{
	my $level = " 10"; 
    my ($a,$q) = tmpMULTEXPR();
	$q =~ s/aa/xx/;
	$q =~ s/ab/x/;
	$q =~ s/bb//;

	$a =~ s/b//g;
	$a =~ s/a/x/g;
    my ($WA1,$q1) = tmpMULTEXPR();
    my ($WA2,$q2) = tmpMULTEXPR();

	my @WA ; 
	my @methods = qw ("factorization", "square", "roots");
	my $WHATMETHOD = util_pick_n_random_from_list(\@methods,1);
    my $question = "Find the roots of $q using factorization method";
    my $answer = $a ; 
	my $hint = "";
	push @WA , $WA1 ;
	push @WA , $WA2 ;
	my ($ANSLIST,$correctindex) = QP_InsertCorrectAnswer($answer,@WA);
    ($question,$answer,$hint,$ANSLIST,$correctindex);

}
