#!/bin/csh 

if($#argv != 3  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set list = `getTmpFileName.pl`

newfile.csh $list
 
echo $1 >> $list
echo $2 >> $list

set simi=$3

$SRC/MISC/needleFastalist.pl  -out qwqwqq -simi $simi -needleout oooooooooo -arg ~/needle.arg -list $list

unlink qwqwqq
unlink $list
cat oooooooooo

