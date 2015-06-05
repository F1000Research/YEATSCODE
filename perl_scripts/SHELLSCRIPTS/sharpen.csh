#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : "
  exit 
endif 

set listref = $1


foreach i ( ` cat $listref` )
  convert $i -sharpen 0x3 $i.out
  mv -f $i.out $i 
end

