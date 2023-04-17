jellyfish count -m 13 -s 200M -C -t 20 Acorus.chromosome_level.genome.fa -o Acorus.chromosome_level.genome_counts.jf 
jellyfish dump -L 100 Acorus.chromosome_level.genome_counts.jf > Acorus.chromosome_level.genome_counts.fa
perl split2chr.pl Acorus.chromosome_level.genome.fa
sh run_jellyfish.sh > run_jellyfish.log 2> run_jellyfish.err
perl kmer_cal.pl Acorus.chromosome_level.genome_counts.fa Chr.lst >Kmer_freq.txt
perl filter_kmer.pl Kmer_freq.txt >filter_kmer.txt 2>err
