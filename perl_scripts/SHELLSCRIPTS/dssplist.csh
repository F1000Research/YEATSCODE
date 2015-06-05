#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set list = $1 

touch HTH
foreach i (`cat $list`)
	if(! -e "$DSSP/$i.dssp") then
	    echo mkdssp -i $PDBDIR/$i.pdb -o $DSSP/$i.dssp
	    mkdssp -i $PDBDIR/$i.pdb -o $DSSP/$i.dssp
	endif 
	if(! -e $HELIXDIR/$i.HELIXLIST) then 
	    helixparseDSSPoutput.pl -outf $HELIXDIR/$i.HELIXLIST -p $i -inf $DSSP/$i.dssp
    endif 
end 

 

