params {
  input = "s3://aws-batch-genomics-shared/secondary-analysis/example-files"
  reference = "s3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
  sample_id = "NIST7035"
  output = "s3://pwyming-tmp-us-east-1"
}

ref_name = file(params.reference).name
ref_indices = Channel.fromPath("${params.reference}*").toList()

fastq = file("${params.input}/fastq")
Channel
  .fromFilePairs("${fastq}/${params.sample_id}_*{1,2}*{fastq.gz}")
  .set { fastqs }

process bwa_mem {
    container "bwa:latest"
    cpu 8
    memory "64 GB"

  input:
    file '*' from ref_indices
    set key, file(reads) from fastqs
  
  output:
    file "${params.sample_id}.sam" into sam
  
  script:
  """
  bwa mem -t 16 -p \
        ${ref_name} \
        ${params.sample_id}_*1*.fastq.gz \
        > ${params.sample_id}.sam
  """
}

/*
process samtools_sort {
    container "samtools:latest"
    cpu 8
    memory "32 GB"

  input:
  
  output:
  
  script:
  """
  samtools sort \
        -@ 16 \
        -o $OUTPUT_PATH/${SAMPLE_ID}.bam \
        $INPUT_PATH/${SAMPLE_ID}.sam
  """
}

process samtools_index {
    container "samtools:latest"
    cpu 8
    memory "32 GB"

  input:
  
  output:
  
  script:
  """
  samtools index \
        $INPUT_PATH/${SAMPLE_ID}.bam
  """
}

process bcftools_mpileup {
    container "bcftools:latest"
    cpu 8
    memory "32 GB"

  input:
  
  output:
  
  script:
  """
  bcftools mpileup \
        --threads 16 \
        -Oz \
        -r chr21 \
        -f $REFERENCE_PATH/${REFERENCE_NAME}.fasta \
        $INPUT_PATH/${SAMPLE_ID}.bam \
        > $OUTPUT_PATH/${SAMPLE_ID}.mpileup.vcf.gz
  """
}

process bcftools_call {
    container "bcftools:latest"
    cpu 8
    memory "32 GB"

  input:
  
  output:
  
  script:
  """
  bcftools call \
        -m \
        --threads 16 \
        -t chr21 \
        -o $OUTPUT_PATH/${SAMPLE_ID}.vcf \
        $INPUT_PATH/${SAMPLE_ID}.mpileup.vcf.gz
  """
}
*/