database="?/uniprot_sprot.fasta"
program="?/swsharpdbh"
out="?"
params="--nocache --max-aligns 500 --outfmt bm0"

{ time $program -j $database -i ?/humdiv-neutral.fa $params -s 3 > divneu.3.out ; } 2>> divneu.3.err

{ time $program -j $database -i ?/humdiv-neutral.fa $params -s 3 -p > divneu.3p.out ; } 2>> divneu.3p.err

#mv *.err *.out $out

