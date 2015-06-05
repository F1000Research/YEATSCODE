#!/bin/csh -f

if($#argv != 5  ) then 
  echo "Usage : "
  exit 
endif 

set PWD = ` pwd`
set listref = $PWD/$1
set listquery = $PWD/$2
set dir = $PWD/$3
set anndir = $PWD/$4
set close2activesite = $5

mkdir -p $dir

foreach ref ( ` cat $listref` )
  mkdir -p $dir/$ref 	
  cd $dir/$ref 	
  cp -rf  $CONFIGGRP . 

  if( ! -e $ref.pdb.out   ) then 
      cp -r  $anndir/$ref.outconf.annotated . 
      newfile.csh LLL
      echo $ref $ref $anndir/$ref.outconf.annotated $dir/$ref/$ref.pdb.out >> LLL
      $SRC/3DMatch -annd $anndir -outf llllll -pdb $ref -list LLL -newfindresidues -close2activesite $close2activesite 
  endif 

  cd $PWD 
end

echo "Starting queries "
set FFF = $dir/LLL
foreach query ( ` cat $listquery` )


  newfile.csh $FFF >> /dev/null
  echo -n "$query "
  foreach ref ( ` cat $listref` )
       cd $dir/$ref 	
	   if(! -e $query.pdb.out) then 
               echo $ref $query $anndir/$ref.outconf.annotated $dir/$ref/$query.pdb.out >> $FFF
	   endif 
	   cd $PWD
  end
  
   if (! -z $FFF) then 
     echo -n "$query "
     $SRC/3DMatch -annd $anndir -outf llllll -pdb $query -listfile $FFF -newfindresidues -close2activesite $close2activesite 
  endif



  cd $PWD 
end
echo ""
