#!/bin/csh -f

if($#argv != 3  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 

set PWD = ` pwd`
#set WORKDIR = $PWD/PREMONITION
set list = $PWD/$1
set mlen = $2
set config = $PWD/$3

#mkdir -p $WORKDIR/$mlen
#cd $WORKDIR/$mlen
foreach i (`cat $list`)
	if(! -e "$i.premon.out.zip") then 
       $SRC/ALIGN/premonition.pl -p1 $i  -ml $mlen -list $config
	   $SRC/ALIGN/premonitionMerge.pl -out $i.premon.out -p1 $i
	   zip -r $i.premon.out.zip $i.premon.out
	   \rm $i*singleout
	   \rm $i.premon.out
   endif 
end 


