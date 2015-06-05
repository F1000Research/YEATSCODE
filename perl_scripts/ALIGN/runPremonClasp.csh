#!/bin/csh -f

if($#argv != 8  ) then 
  echo "Usage : "
  exit 
endif 

set PWD = ` pwd`
set listref = $1
set listquery = $2
set readpotential = $3
set close2active = $4
set size = $5
set aalist = $6
set annotate = $7 
set tag = $8 

set premondir = $PREMONITION/$size/
set results = Results.$close2active.premon.$listref.$listquery
mkdir -p $results
\cp -f $aalist $results

set annotatefixed = $results/ANNOTATE.$size
mkdir -p $results/ANNOTATE.$size

\cp -f mapfileHolo2Apo $annotatefixed


set errlog=$PWD/premon.log
newfile.csh $errlog 
 
set fixedreflist=$PWD/$results/fixedreflist
newfile.csh $fixedreflist 

echo "CONVERTING premon.in"
convertPremoninBasedonAAlist.pl -anndir1 $annotate -anndir2 $annotatefixed -aalist $aalist -size $size -outf $fixedreflist -lis $listref  -tag $tag
#foreach i (`cat $listref`)
	#if(! -e "$annotate/$i.$size.premon.in") then 
		#echo "$annotate/$i.$size.premon.in does not exist" >> $errlog
	#else
		#echo $i >> $fixedreflist
	#endif
#end 


foreach i (`cat $listquery`)
	if(! -e "$premondir/$i.premon.out.zip") then 
		echo "$premondir/$i.premon.out.zip does not exist" >> $errlog
	else
		mkdir -p $results/$i
		cd $results/$i 

		if(! -e "$i.cumu.scores") then 
		     \rm -f $i.premon.out
		     \rm -f $premondir/$i.premon.out ## just remove this too...
		     unzip $premondir/$i.premon.out.zip 
     
		     echo "Running premon on query $i"
		     $SRC/ALIGN/runPremonClasp.pl -protein $i -premon $i.premon.out -list $fixedreflist -size $size -rundir $PWD -errlog $errlog -conf $CONFIGGRP -readpo $readpotential -close2active $close2active -anndir $annotatefixed -tag $tag
		     $SRC/ALIGN/getPremonCumu.csh $i finallist 


		     \rm -f $i.premon.out
		else
		   echo "$i already done" >> $errlog
		endif 

		cd $PWD
	endif
end 


echo "-=================  error log in $errlog ========================"
