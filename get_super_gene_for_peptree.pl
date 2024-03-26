#!/usr/bin/perl -w
use FindBin qw($Bin $Script);

BEGIN{
     unshift (@INC,"$Bin/Soft/");
}

my %hash;
use Bio::SeqIO;
my $inputfile1=$ARGV[0];
my $inputfile2=$ARGV[1];
my $species=$ARGV[2];
my @name;
open(NAME,$species);
while(<NAME>){
   chomp;
   my @tmp=split/\s+/;
   push @name,$tmp[0];
}
close NAME;

open (OUT, ">$inputfile2")|| die "cannot open $inputfile1:$!";
my @file=glob("$inputfile1/OG*.fa.mafft.fa-gb");
for(my $i=0;$i<=$#file;$i++){
my $ina = Bio::SeqIO->new(-file => $file[$i], -format => 'fasta');
        while(my $obj = $ina->next_seq()){
              my $id = $obj->id;
		 my $seq = $obj->seq;
         push @{$hash{$file[$i]}},[($id,$seq)];
}
}
my @array2;
my $num;
foreach my $key(keys %hash){
  my @array=@{$hash{$key}};
  for(my $i=0;$i<=$#array;$i++){
    $array2[$i].=$array[$i][1];
  }
  $num=$#array+1;
}
my $len=length($array2[0]);
print OUT "$num  $len\n";
for(my $i=0;$i<=$#array2;$i++){
$len=length($array2[$i]);
print OUT "$name[$i]  $array2[$i]\n";
}
close OUT;
