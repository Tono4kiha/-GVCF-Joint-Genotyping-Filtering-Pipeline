#!/bin/bash
set +o posix
set -eo pipefail

genome=
Kmer=25
Number_of_errors=2
cpu=40
genmap=
bedtools=
pigz=

${genmap} index -F ${genome} -I index_genmap
[ -d out_genmap ] || mkdir out_genmap
${genmap} map -K ${Kmer} -E ${Number_of_errors} -I index_genmap -O out_genmap -bg --threads ${cpu}
awk -F "\t" '($4+0.0)>0.499' out_genmap/*.genmap.bedgraph > out_genmap/high_mappability.bed
${bedtools} merge -i out_genmap/high_mappability.bed > out_genmap/high_mappability.merge.bed
${pigz} --best -p ${cpu} out_genmap/high_mappability.bed
${pigz} --best -p ${cpu} out_genmap/*.genmap.bedgraph
rm -rf index_genmap 
