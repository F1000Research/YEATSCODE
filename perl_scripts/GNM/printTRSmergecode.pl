
foreach my $i (19..100){
    print "elsif(\$size eq $i){\n";
       print "(\$begin) = (\$seq =~/^(" ; 
foreach my $j (1..$i){
print ".";
}
       print ")/); \n";
       print "(\$end) = (\$seq =~/(" ; 
foreach my $j (1..$i){
print ".";
}
       print ")\$/); \n";
    print "}\n";
}
