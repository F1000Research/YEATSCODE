#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set list = $1

source $SRC/setup.csh 

newfile.csh allnames

foreach i (`cat $list`)
echo -n "$i  "
$SRC/ALIGN/helixidentify.pl -outf ooooooo -aal ~/aalist -con $CONFIGGRP -pr $i > & ! /dev/null
#$SRC/ALIGN/helixidentify.pl -outf ooooooo -aal ~/aalist -con $CONFIGGRP -pr $i 
end 

echo 
