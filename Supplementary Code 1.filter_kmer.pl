#!/usr/bin/perl

open FREQ,"$ARGV[0]";
$_=<FREQ>;
print "$_";
while(<FREQ>){
   chomp;
   $line=$_;
   @item=split;
   $otherA=0;$otherB=0;
   @array_A=();
   @array_B=();
   for($i=2;$i<=23;$i++){
      if($item[$i]==0){$item[$i]=0.01;}
   }
   for($i=4;$i<=11;$i++){
      push @array_A,$item[$i];
      $otherA+=$item[$i];
   }
   for($i=12;$i<=21;$i++){
      push @array_B,$item[$i];
      $otherB+=$item[$i];
   }

   if($item[2]>$item[23]){$max_1=$item[2];$min_1=$item[23];}
      else{$min_1=$item[2];$max_1=$item[23];}
   if($item[3]>$item[22]){$max_2=$item[3];$min_2=$item[22];}
      else{$min_2=$item[3];$max_2=$item[22];}

   @sort_A=sort {$a <=> $b} @array_A;
   @sort_B=sort {$a <=> $b} @array_B;
   if($otherA>$otherB){
      $max_3=$otherA;$min_3=$otherB;
      if($sort_A[2]<=$sort_B[7]){next;}
      if($sort_A[7]/$otherA>0.6){next;}
      if($sort_B[9]>$sort_A[2]){next;}
      if($sort_A[0]<$sort_B[7]){next;}
#      $max_33=$sort_A[2];$min_33=$sort_B[7];
   }
   else{
      $min_3=$otherA;$max_3=$otherB;
      if($sort_B[2]<=$sort_A[5]){next;}
      if($sort_B[9]/$otherB>0.6){next;}
      if($sort_A[7]>$sort_B[2]){next;}
      if($sort_B[0]<$sort_A[5]){next;}
#      $max_33=$sort_B[2];$min_33=$sort_A[5];
   }

   if(($max_1/$min_1>2)&&($max_2/$min_2>2)&&($max_3/$min_3>2)){
      if($max_1>$max_2){$max=$max_2;}
         else{$max=$max_1;}
      if($min_1<$min_2){$min=$min_2;}
         else{$min=$min_1;}
#      if(($max>$min)&&($max>$max_33)){#&&($min<$min_33)){
      if($max>$min){
         print "$line\n";
         print STDERR "$item[0]\t$min_1\t$max_1\t$min_2\t$max_2\t$sort_A[1]\t$sort_A[6]\t$sort_B[1]\t$sort_B[8]\n";
      }
   }
}
close FREQ;
