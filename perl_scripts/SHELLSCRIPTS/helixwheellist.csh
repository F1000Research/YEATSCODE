#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set list = $1 

touch allvalues
touch commands
#setenv PDBDIR $HELIXDIR
foreach i (`cat $list`)
      #helixwheel.pl -out ooooooo -con $CONFIGGRP -aa ~/aalist -prote $i > & ! /dev/null
      helixwheel.pl -out ooooooo -con $CONFIGGRP -aa ~/aalist -prote $i 
end 

 

