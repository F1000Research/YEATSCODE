
sub PrintABSequence{
	my ($self) = @_ ;
	my $retlist = {};

    my @res = $self->GetResidues();
    my $N = @res;
	
	my $prevnum ; 
	my $numtable = {};

	my $state = 0 ;
	
	my $fname = $self->{REALNAME} . ".AB.fasta";
	my $ofh = util_write($fname);
	print $ofh ">$self->{REALNAME}AB\n";

	my $donepocketH = {};
	my $donepocketB = {};
    while(@res){
		my $res = shift @res ;
	    next if($res->GetAtomStr() eq "HETATM");
	    next if($res->GetName() eq "HOH");
		my $num = $res->GetResNum();
		my $isonhelx  = $self->IsResidueNumOnHelix($num);
		my $isonbeta  = $self->IsResidueNumOnBETA($num);
		die if ($isonhelx > -1  && $isonbeta > -1);
		if(!$state){
		   if($isonhelx > -1 || $isonbeta > -1  ){
		   	  my $what2print = $isonhelx > -1 ? "H" : "B";
			  print $ofh "$what2print";
		   	  $state = $isonhelx > -1 ? 1 : 2 ;
		   }
		}
		else{
		   if($isonhelx eq -1 && $isonbeta eq -1){
		   	  $state = 0 ;
		   }
		}
		if($state){
			my $IsDisulphide = $self->IsDisulphide($num);
			if($IsDisulphide){
				print $ofh "C";
			}
		}
	}
	print $ofh "\n";
}
