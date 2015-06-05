#!/usr/bin/perl -w 
use strict ;
use PDB;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use ConfigPDB;
use MyGeom;
use MyPymol;
use Math::Geometry ;
use Math::VectorReal qw(:all);  # Include O X Y Z axis constant vectors




use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($ann,$resnum,$config,$p1,$p2,$infile,$outfile,$which_tech,$listfile,$protein);
my $maxdist ;
my $DISTANCEWITHOUTSEQMATCH = 1 ;
my $verbose = 0 ;

my ($all,$radii,$before1,$before2);
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "all"=>\$all ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "ann=s"=>\$ann ,
            "maxdist=f"=>\$maxdist ,
            "config=s"=>\$config,
            "radii=i"=>\$radii ,
            "resnum=i"=>\$resnum ,
            "verbose=i"=>\$verbose ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a config file name => option -config ") if(!defined $config);
usage( "Need to give a file name => option -protein ") if(!defined $protein);
usage( "Need to give a file name => option -radii ") if(!defined $radii);


my $ofh = util_write($outfile);
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $PWD = cwd;

ConfigPDB_Init($config);

my $pdb1 = "$PDBDIR/$protein.pdb";

my $pdb = new PDB();
$pdb->ReadPDB($pdb1);



my @fullatomlist = $pdb->GetAtoms();
my @atomlist; 
foreach my $a (@fullatomlist){
		#$a->Print();
		my $atomstr = $a->GetAtomStr();
		next if($atomstr eq "HETATM");

		my $type = $a->GetType();
		next if($type ne "CA");

		push @atomlist,$a;
		my $resnum = $a->GetResNum();
		print "$type $resnum \n";
}

my @tmplist = @atomlist ;

my $disttable = {};
while(@tmplist){
	my $a = shift @tmplist;
	foreach my $b (@tmplist){
	        my $d = util_format_float($pdb->DistanceAtoms($a,$b),1) ;
			$disttable->{$a} = {} if(!exists $disttable->{$a});
			$disttable->{$b} = {} if(!exists $disttable->{$b});
			if($d <= $radii){
			     $disttable->{$a}->{$b} = $d; 
			     $disttable->{$b}->{$a} = $d; 
			}

	}
}

my @FullMatrix = ();

foreach my $a (@atomlist){
	my ($X,$Y,$Z) = $a->Coords();

	foreach my $IDX (1..3){
         my @row ; 		

	    my $sumX = 0 ;
	    my $sumY = 0 ;
	    my $sumZ = 0 ;
        foreach my $b (@atomlist){
	        my ($x,$y,$z) = $b->Coords();
		    my $val;
		    if($a eq $b){
			    $val = -10000 ;
		        push @row, $val ;
		        push @row, $val ;
		        push @row, $val ;
		    }
		    else{
				my $FIRST = 0 ;
			    if(exists $disttable->{$a}->{$b}){
				    $val = 1 ;
				    my $dist = $disttable->{$a}->{$b} ;
				    my $distsquared  = $dist* $dist;
                    if($IDX eq 1){
						$FIRST = ($X-$x);
	                }
	                elsif ($IDX eq 2){
						$FIRST = ($Y-$y);
	                }
	                else{
						$FIRST = ($Z-$z);
	                }
		            my $val1 = util_format_float(($FIRST*($X -$x))/$distsquared,1);
		            my $val2 = util_format_float(($FIRST*($Y -$y))/$distsquared,1);
		            my $val3 = util_format_float(($FIRST*($Z -$z))/$distsquared,1);

					if($verbose){
		            print "my $val1 = util_format_float(($FIRST*($X -$x))/$distsquared,1);\n";
		            print "my $val2 = util_format_float(($FIRST*($Y -$y))/$distsquared,1);\n";
		            print "my $val3 = util_format_float(($FIRST*($Z -$z))/$distsquared,1);\n";
					exit ;
					}

		            push @row, $val1 ;
		            push @row, $val2 ;
		            push @row, $val3 ;
					$sumX = $sumX + $val1 ;
					$sumY = $sumY + $val2 ;
					$sumZ = $sumZ + $val3 ;
    
			    }
			    else{
				    $val = 0 ;
		             push @row, $val ;
		             push @row, $val ;
		             push @row, $val ;
			    }
			    
		    }
    
        }

	    my $sum = 0 ;
	    my @tmprow = @row;
	    while(@tmprow){
		    my $x = shift @tmprow ;
		    my $y = shift @tmprow ;
		    my $z = shift @tmprow ;
		    next if($x eq -10000) ;
		    $sum = $sum + $x + $y+ $z ; 
	    }

		my @SUMARR ;
		push @SUMARR, $sumX; 
		push @SUMARR, $sumY; 
		push @SUMARR, $sumZ; 

	    my $N = @row ;
	    foreach my $idx (0..$N-1){
		    if($row[$idx] eq -10000){
		        die  if(!@SUMARR);
			    $row[$idx] = - (shift @SUMARR) ;
		    }
	    }
		die  if(@SUMARR);
    
	    push @FullMatrix,\@row ;
	} ### foreach x,y,z column

}

print $ofh "use strict;\n";
print $ofh "use PDL; \n";
print $ofh "use PDL::Matrix; \n";
print $ofh "use PDL::MatrixOps; \n";
print $ofh "my \$mat = PDL::Matrix->pdl([" ; 
foreach my $row (@FullMatrix){
	my $str = join ",", @{$row};
    print $ofh "[$str], ";
}
print $ofh "]); \n";
print $ofh "my (\$ev,\$e) = eigens \$mat; \n";
print $ofh "print \$ev ; \n";
print $ofh "print \$e ; \n";


print "@$_\n" for @FullMatrix;

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
