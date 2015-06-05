#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set list = $1 

if(-e "F")then 
	echo "please remove file F"
   exit 
endif 
touch F 
foreach i (`cat $list`)
    
	writeFasta.pl -p $i
	cat $i.ALL*fasta >> F 
end 

$SRC/MISC/processStringFromFasta.pl -outf ooo -in F -cutof 0 
cat ooo | grep "^>" > ! ppp 
cat ppp | perl -ne ' s/.//;  my @l = split ; print $l[0], "\n" ;' > ! list.pdb 

 

