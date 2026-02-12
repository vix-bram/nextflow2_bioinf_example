# nextFlow bioinformatics example use case

This workflow performs the trimming and assembly of raw compressed fastq reads (R1.fastq.gz and R2.fastq.gz). Then, quality assessment and genotyping is performed on the assembly. The tools needed to run this workflow are:
*fastp
skesa
quast
mlst*



You can run it on the terminal with the following command:

nextflow run bioinf\_example\_nextflow2.nf --quast\_path /path/to/quast.py

**Please hardcode the paths of the two (forward and reverse strands) input raw fastq files into the nextflow script.**

