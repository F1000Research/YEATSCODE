#!/bin/csh 

if($#argv != 1  ) then 
  echo "exec> DB "
  exit 
endif 

## all files are created within the MERGED dir
set db=$1 

mkdir -p MERGED

newfile.csh alreadydone
foreach i ( 5 6 7 8 9 10 11 12 13 14 15 )
  if(! -e MERGED/data.MERGED.OVERLAP$i) then 
      joinTRS.pl -inf $db -outf data.MERGED.OVERLAP -size $i
 endif 
end


foreach i ( 5 6 7 8 9 10 11 12 13 14 15 )
      #joinTRSPostProcess.pl -inf data.MERGED.OVERLAP -size $i -anno annotations.txt
      joinTRSPostProcess.pl -inf data.MERGED.OVERLAP -size $i 
end
