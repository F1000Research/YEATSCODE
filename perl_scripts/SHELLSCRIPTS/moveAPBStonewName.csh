#!/bin/csh -f

if($#argv != 2  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set nm1 = $1 
set nm2 = $2


if(! -e $APBSDIR/$nm2) then
   mv $APBSDIR/$nm1 $APBSDIR/$nm2 
   cd $APBSDIR/$nm2
   mv $nm1.csv $nm2.csv
   mv $nm1.pqr $nm2.pqr
endif 

