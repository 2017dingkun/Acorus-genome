#!/usr/bin/perl
$name="start";
open FA, "$ARGV[0]";
while(<FA>){
   chomp;
   if(/>Chr(\d+)/i){
      $name="Chr".$1;
      unless($name eq "start"){
         close CHR;
      }
      open CHR,">$name.fa";
      print CHR ">$name\n";
   }
   else{
      print CHR "$_\n";
   }
}
close FA;
close CHR;
