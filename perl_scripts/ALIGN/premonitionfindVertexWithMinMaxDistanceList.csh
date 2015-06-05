#!/bin/csh -f

if($#argv != 2  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set PWD = ` pwd`
set list = $1
set anndir = $2 
set out = premon.config
newfile.csh $out
foreach i (`cat $list`)
    findVertexWithMinMaxDistance.pl -out $out -in $anndir/$i.outconf.annotated -pr $i 
end 
