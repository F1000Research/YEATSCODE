#!/bin/csh 

if($#argv != 4  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 

cat *trs.out > ! list 
UNIQ list

onlyuniquePrefix.pl -in list.uniq -out list.uniq.prefix

wc -l list.uniq.prefix
#foreach i (`cat list.uniq.prefix`)
foreach i (`cat list.uniq`)
    if(! -e "$i.trs.out") then 
		echo $i
    endif
end 

foreach i (`cat list.uniq`)
    if(! -e "$i.trs.out") then 
         setenv FASTADIR ../FASTADIR_ORFNEW
		 $SRC/BLAST/findgene.csh $i $1 $2 $3 $4 
    endif
end 

\rm *reflection*
