#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : ./apbs.csh  <list>  "
  exit 
endif 

set PWD = ` pwd`
set list = $1

foreach i (`cat $list`)
	if(-e "$APBSDIR/$i.zip") then
		cd $APBSDIR
	    unzip $APBSDIR/$i.zip
		unlink $APBSDIR/$i.zip
		cd -
	else if(! -e $APBSDIR/$i) then 
       $SRC/APBS/apbssingle.csh $i 
       cd $i 
       moveapbsfiles.csh
       cd -
	   mv $i $APBSDIR
   endif
end 

