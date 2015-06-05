#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 

set PWD = ` pwd`
set WORKDIR = $PWD/PREMONITION
mkdir -p $WORKDIR
set list = $PWD/$1
foreach dist ( 15 20 )
	mkdir -p $WORKDIR/$dist
	cd $WORKDIR/$dist 
    foreach mlen ( 3 4 )
	    mkdir -p $WORKDIR/$dist/$mlen
	    cd $WORKDIR/$dist/$mlen
        foreach i (`cat $list`)
		$SRC/ALIGN/premonition.pl -outf $i.out -p1 $i -dis $dist -ml $mlen
        end 
		cd -
    end
	cd - 
end


