#!/bin/csh -f

if($#argv != 4  ) then 
  echo "Usage wrong for extractRefQuertList.cshctRefQuertList.csh as $#argv is not 4  "
  exit 
endif 

set PWD = ` pwd`
set listref = $PWD/$1
set listquery = $PWD/$2
set dir = $PWD/$3
set resultsdir = $PWD/$4

set runfile = $PWD/run.csh

mkdir -p $dir


echo processing lists $listquery and $listref. Please wait...

foreach query ( ` cat $listquery` )

  mkdir -p $dir/$query 	 
  cd $dir/$query/
  newfile.csh list2run  > & ! /dev/null
  echo -n "."
  foreach ref ( ` cat $listref` )
      if( ! -e $query.$ref.pdb.out || -z $query.$ref.pdb.out ) then  
	      if( (! -e $dir/$query/$query.$ref.pdb.out)  && (  -e $resultsdir/$ref/$query.pdb.out ) && \
			          ( ! -z $resultsdir/$ref/$query.pdb.out) && (  -e $resultsdir/$ref/$ref.pdb.out ) && ( ! -z  $resultsdir/$ref/$ref.pdb.out) ) then
		  	echo $ref >> list2run 
		  endif 
	  endif 

  end

  cd -
end

echo "."

echo "Running::: $SRC/CLASP -q $query -dontrunapbs  -relax  -reflist list2run -link $resultsdir -querylist $listquery -maindir $dir -outfile error.log. Logged into dev null  "
$SRC/CLASP -q $query -dontrunapbs  -relax  -reflist list2run -link $resultsdir -querylist $listquery -maindir $dir -outfile error.log  
