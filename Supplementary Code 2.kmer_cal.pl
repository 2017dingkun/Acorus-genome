#!/usr/bin/perl

open ALL,"$ARGV[0]";
while(<ALL>){
   chomp;
   if(/>(\d+)/){
      $kmer=<ALL>;
      chomp($kmer);
      $t_hash{$kmer}=$1;
   }
}
close ALL;

open LIST,"$ARGV[1]";
while($chr=<LIST>){
   chomp($chr);
   push @chr,$chr;
   open FA,"$chr.jf.fa";
   while(<FA>){
      chomp;
      if(/>(\d+)/){
         $kmer=<FA>;
         chomp($kmer);
         if(defined($t_hash{$kmer})){
            $hash{$kmer}{$chr}=$1;
         }
      }
   }
   close FA;
}
close LIST;

print "13mer\tALL";
for($i=0;$i<=$#chr;$i++){
   print "\t$chr[$i]";
}
print "\n";

foreach $kmer (keys %t_hash){
   print "$kmer\t$t_hash{$kmer}";
   for($i=0;$i<=$#chr;$i++){
      if(!defined($hash{$kmer}{$chr[$i]})){
         $hash{$kmer}{$chr[$i]}=0;
      }
      print "\t$hash{$kmer}{$chr[$i]}";
   }
   print "\n";
}
