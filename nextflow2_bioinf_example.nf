#!/usr/bin/env nf
nextflow.enable.dsl=2

params.quast_path = '' // User must provide this path

if (params.quast_path == '') {
    println("""\nERROR: Missing required parameters.
    Usage:
        nextflow run nextflow2_bioinf_example.nf --quast_path /path/to/quast.py
    Please provide all required file paths.\n""")
    System.exit(1) // Exit if parameters are not provided
}
// Trimming the reads
process trim_fastp {
  // maxForks 1
  // Input files
  input:
  file reads1
  file reads2
  
  // Output files
  output:
  file "trimmed_R1.fq.gz" 
  file "trimmed_R2.fq.gz"
  
  // Script to run fastp
  script:
  """
  fastp -i $reads1 -I $reads2 -o trimmed_R1.fq.gz -O trimmed_R2.fq.gz
  """
}

// File assebly
process assemble_skesa {
  // maxForks 1
  // Input files
  input:
  file "trimmed_R1.fq.gz" 
  file "trimmed_R2.fq.gz"
  
  // Output files
  output:
  path 'assembly.fasta'
  path 'skesa.stdout.txt'
  path 'skesa.stderr.txt'
  
  // Script to run skesa
  script:
  """
  skesa --fastq trimmed_R1.fq.gz trimmed_R2.fq.gz --contigs_out assembly.fasta 1> skesa.stdout.txt 2> skesa.stderr.txt
  """
}


// Quality Assessment
process qa_quast {
  // Input files
  input:
  path "assembly.fasta" 
  path "skesa.stdout.txt"
  path "skesa.stderr.txt"
    
  // Output files
  output:
  path "quast_output"

  // Performing quality assessment on the assemblies
  """
  python ${params.quast_path} assembly.fasta -o quast_output
  """
}

// Using mlst to characterize the bacterial strain
process proc_mlst {
  // Input files
  input:
  file "assembly.fasta" 
  path "skesa.stdout.txt"
  path "skesa.stderr.txt"

  // Output files
  output:
  path "mlst.tsv"
    
  // Running mlst
  script:
  """
  mlst assembly.fasta > mlst.tsv
  """
}


workflow {

  // Print instructions about input files
  println("Please ensure that the input reads R1.fastq.gz and R2.fastq.gz are in the current directory where this Nextflow script is being run.")

  // Hardcoding the input files (taking in inputs didnt work)
  reads_ch_R1 = file("E0475741_R1.fastq.gz")
  reads_ch_R2 = file("E0475741_R2.fastq.gz")

  trim_fastp(reads_ch_R1, reads_ch_R2)
  | assemble_skesa
  | qa_quast & proc_mlst
  | mix
  | view
}
