while read line
do
   jellyfish count -m 13 -s 200M -C -t 20 ${line}.fa -o ${line}.jf
   rm ${line}.fa
   jellyfish dump ${line}.jf > ${line}.jf.fa
   rm ${line}.jf
done<Chr.lst
