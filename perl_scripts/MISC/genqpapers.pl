#!/usr/bin/perl -w 

use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;

use MyUtils;
use MyQPaper;
use POSIX qw(floor);

my $commandline = util_get_cmdline("",\@ARGV) ;

my ($easy, $nodecimals, $infile,$qpaper,$howmany,$design,$multiple,$answer,$clocked);
$howmany = 0 ;
$multiple = 1 ;
GetOptions(
            "howmany=s"=>\$howmany ,
            "infile=s"=>\$infile ,
            "clocked:s"=>\$clocked ,
            "design=s"=>\$design ,
            "answer=s"=>\$answer ,
            "multiple=i"=>\$multiple ,
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

my $ofh = util_write($qpaper);
my $answerfh = util_write($answer);
print STDERR "Writing qpaper file $qpaper and answers in $answer\n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}

my @qtypes = util_read_list_words($infile);

while($multiple){
    my $cnt = 0 ;
    #print $ofh      " ================================\n";
    print $ofh      " QUESTIONS id=$id  $multiple \\\\\n";
    print $ofh      " ================================ \\\\\n";
    #print $answerfh      " ================================\n";
    print $answerfh " ANSWERS id=$id $multiple \\\\\n";
    print $answerfh      " ================================ \\\\\n";
    
    $multiple--;

    foreach my $iii (0..$howmany){
        foreach my $qtype (@qtypes){
            $cnt++;
            my ($q,$a) ;
		    my $sub = \&{$qtype}; 
		    ($q,$a) = $sub->();
            print $ofh      " $cnt ) {\\bf $q}\\\\ \n";
            print $answerfh " $cnt ) {\\bf $a}\\\\ \n";
        }
    }
}

sub ADD{
    my ($div) = @_ ;
    my $n1 = defined $easy ? 100 : 10000 ;
    my $n2 = defined $easy ? 10 : 100.0  ;

    $n2 = 1 if(defined $nodecimals);
    my ($p,$q);
   
    while(!($p = floor($n1*rand())/$n2)){}
    while(!($q = floor($n1*rand())/$n2)){}
 
    my $addorsub = 1 * rand() > 0.5 ? 0 : 1 ; 
    
    my $a = $addorsub ? $p + $q : $p - $q  ;
    my $oper = $addorsub ? " + " :  " - "  ;
    my $answer = "$p $oper $q\ = $a";
    my $question = "$p $oper $q\ = ? ";
    ($question,$answer);
}

sub ADDMULTIPLE{
    my ($div) = @_ ;
    my $n1 = 8 ; 
	my $sum = 0 ; 
	my $str = "";
    while($n1){
       my $celprob = util_round(rand());
       my $p = floor(20*rand()) + 1 ; 
	   $str = $celprob ? "$str + $p " : "$str - $p"; # do this first
	   $p = $celprob ? $p : $p * (-1) ; 
	   $sum = $sum + $p ;
	   $n1--;
	}
    my $answer = "$str = $sum";
    my $question = "Find the sum of: $str  ";
    ($question,$answer);
}

sub DIV{
	return MULT(1);
}

sub MULT{
    my ($div) = @_ ;
    my $n1 = defined $easy ? 100 : 10000 ;
    my $n2 = defined $easy ? 10 : 100.0  ;
    if(defined $div){
        $nodecimals = 1 ;
         $n1 = defined $easy ? 20 : 100 ; 
         $n2 = defined $easy ? 10.0 : 10.0;
    }

    $n2 = 1 if(defined $nodecimals);
    my ($p,$q);
   
    while(!($p = floor($n1*rand())/$n2)){}
    while(!($q = floor($n1*rand())/$n2)){}
    
    my $a = $p * $q ;
    my $oper = " X " ;
    if(defined $div){
         $oper = " / " ;
         my $t = $a ; 
         $a = $p ; 
         $p = $t ; 
    }
    my $answer = "$p $oper $q\ = $a";
    my $question = "$p $oper $q\ = ? ";
    undef $nodecimals ;
    ($question,$answer);
}

sub TEMPERATURE{
    my $n1 = 1000 ;
    my $cel = floor($n1*rand());
    my $fahren = (9*$cel)/5 + 32 ; 
    
    my $celprob = util_round(rand());
    if($celprob){
        return "$cel degree C = ? degree F" , " $cel degree C = $fahren F";
    }
    else {
        return "$fahren degree F = ? degree C" , " $fahren degree F = $cel C";
    }
}


sub FRACTIONS{
    my $n1 = 20 ;
	my ($a,$b,$c,$d) = QP_GetNRandomNumbersBelowValue(4,$n1);
    my $answer = ($a/$b + $c/$d) ; 
    $answer = QP_round2place($answer);
    my $q = $a . "/" . $b . "+" . $c . "/" . $d ; 
    return " $q ", " $q = $answer";
}

sub ASCDESC{
    my $n1 = 20 ;
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
    my $answer = $celprob ?  " $p $q $r " : " $r $q $p";
    my $what = $celprob ? " ascending " : " descending ";
    my $question = $a . "/" . $b . " , " . $c . "/" . $d . " , " . $e . "/" . $f  ; 
    return " Arrange in $what order $question ", " $answer";
}

sub CIRCLE{
    my $n1 = 50 ;
    my $r = floor($n1*rand());
    my $a = 3.14 * $r * $r ;
    my $c = 3.14 * 2 * $r ;
    return " Radius of a circle is $r cm. Find its area and circumference" , 
	       "Area = $a sq cm , Circumference = $c cm";
}


sub UNITARYWORK{
    my $n1 = 30 ;
    my $n2 = $n1 ;

	my $WHOM = util_pick_random_from_list(\@PPL);
    my $what = util_pick_random_from_list(\@PLACES);
    
    my $days = floor($n1*rand())+1;
    my $days1 = $days * (floor($n1*rand())+1);

    my $boys =  floor($n2*rand())+2;
    my $boys1 =  $boys * (floor($n2*rand())+1);
    
    my $celprob = util_round(rand());
    if($celprob){
        my $q = "$boys $WHOM can build a $what in $days days. How many days will $boys1 $WHOM build it?"; 
        my $a = QP_round2place(($boys*$days)/$boys1);
        return $q , "$q = $a";
    }
    else {
        my $q = "$boys $WHOM can build a $what in $days days. How many $WHOM will make it in $days1 days ?"; 
        my $a = QP_round2place(($boys*$days)/$days1);
        return $q , "$q = $a";
    }
}

sub UNITARYCOST{
    my $n1 = 50 ;
    my $n2 = 50 ;

    my $what = util_pick_random_from_list(\@THINGS);
    
    my $cost = floor($n1*rand())+1;
    my $cost1 = $cost* (floor(20*rand())+1);

    my $num =  floor($n2*rand())+1;
    my $num1 =  $num * (floor(20*rand())+1);
    
    my $celprob = util_round(rand());
    if($celprob){
        my $q = "$num $what costs Rs $cost . What is the cost of $num1 $what?";
        my $a = QP_round2place(($cost/$num)*$num1);
        return $q , "$q = $a";
    }
    else {
        my $q = "$num $what costs Rs $cost . How many $what can you buy for Rs $cost1?";
        my $a = QP_round2place(($num/$cost)*$cost1);
        return $q , "$q = $a";
    }
}

sub SQUAREROOTIRRATIONAL{


    my ($num,$num1) =  QP_GetNRandomNumbersBelowValue(2,10);
	my $sum = $num + $num1 ; 
	my $prod  = 4* $num * $num1 ; 
    
    my $celprob = util_round(rand());
    if($celprob){
        my $q = "Find the square root of $sum + sqrt($prod)";
        my $a = "sqrt($num) + sqrt($num1)";
        return $q , "$q = $a";
    }
    else {
        my $q = "Find the square root of $sum - sqrt($prod)";
        my $a = "sqrt($num) - sqrt($num1), or sqrt($num1) - sqrt($num)";
        return $q , "$q = $a";
    }
	
}


sub SUMOFSQUARES{
    my ($n1,$n2) =  QP_GetNRandomNumbersBelowValue(2,10);
	my $square1 = $n1 *$n1 ; 
	my $square2 = $n2 *$n2 ; 
	my $sum = $square1 + $square2 ;


    my $q = "Express $sum as the sum of two squares";
    my $a = "$n1 x $n1 + $n2 x $n2 " ;
    return $q , "$q = $a";
	
}



sub MULTEXPR{
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
    return "multiply $q" , "$q = $a";
}
sub FACTOREXPR{
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
    my $n1 =  floor(10*rand())+1;
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

    my $finalexpr = "";
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
	foreach my $k (keys %{$done}){
		my $v = $done->{$k} ;
		my $expr = "$v$k" ;
		$expr = "+" . $expr if(!($expr =~ /^\s*-/));
		$finalexpr = $finalexpr . " " .  $expr ; 
	}
	return $finalexpr ;
}
