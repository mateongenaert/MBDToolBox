# bedtools to create necessary BED files from BAM to work with BALM
# some processing, sorting etc.

BEDTools-Version-2.15.0/bin/bamToBed -i 11_meth_nodups.bam > 11_meth.bed
BEDTools-Version-2.15.0/bin/bamToBed -i 11_input_nodups.bam > 11_input.bed
sort 11_meth.bed -k1,1 -k2,2n -k3,3n -k6,6 -u > 11_meth_sorted.bed
sort 11_input.bed -k1,1 -k2,2n -k3,3n -k6,6 -u > 11_input_sorted.bed
cut -f1,2,3,6- 11_meth_sorted.bed > 11_meth.bed
cut -f1,2,3,6- 11_input_sorted.bed > 11_input.bed

# Run BALM

BALM1.0.1_linux64_multithread -g hg19 -p 0.99 -b -m 0.1 -c 11_input_sorted.bed -i -f 11_meth_sorted.bed
