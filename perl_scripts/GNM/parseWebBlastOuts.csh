#!/bin/csh 

if($#argv != 2  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set list = $1 
newfile.csh genome.annotated.csv
newfile.csh genome.annotated.morethanone.csv
newfile.csh genome.annotated.none.csv
newfile.csh genome.commands.csh
newfile.csh genome.list
newfile.csh genome.filenotfound
newfile.csh genome.warn


set cutoff=0.000000000001
# do the duplicate one only for a low Evalue
foreach i (`cat $list`)
   parseWebBlastOuts.pl -outf jjj -tr $i -cutoff $cutoff -blast $2 -duplicate
end 


## this just gets the uniquely annotated for a larger Evalue
set cutoff=0.00000001
foreach i (`cat $list`)
   parseWebBlastOuts.pl -outf jjj -tr $i -cutoff $cutoff -blast $2
end 

set cutoff=0.00000001
# this just gets the none - see it has the previous E value
foreach i (`cat $list`)
   parseWebBlastOuts.pl -outf jjj -tr $i -cutoff $cutoff -blast $2 -donone
end 
