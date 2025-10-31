// Move the UMI from I2 to the R1/R2 read ID
process extract_umi {
  tag "Extract UMI for ${sample_id}"
  
  output:
  tuple val(sample_id), path(r1_out_umi), path(r2_out_umi)
  
  input:
  tuple val(sample_id), path(r1), path(r2), path(i2)
  
  script:
  """
  if [ ${params.umi_side} = "5" ]; then
        cutadapt -j ${task.cpus} -u ${params.umi_length} --rename='{id}_{r1.cut_prefix} {comment}' -o i2_out_umi -p r1_out_umi $i2 $r1 > extract_umi_r1.log
        cutadapt -j ${task.cpus} -u ${params.umi_length} --rename='{id}_{r1.cut_prefix} {comment}' -o i2_out_umi -p r2_out_umi $i2 $r2 > extract_umi_r2.log
    else
        cutadapt -j ${task.cpus} -u -${params.umi_length} --rename='{id}_{r1.cut_suffix} {comment}' -o i2_out_umi -p r1_out_umi $i2 $r1 > extract_umi_r1.log
        cutadapt -j ${task.cpus} -u -${params.umi_length} --rename='{id}_{r1.cut_suffix} {comment}' -o i2_out_umi -p r2_out_umi $i2 $r2 > extract_umi_r2.log
  fi
  """
}

// Trim the dsODN tag from the 5' region of the R2 read
process trim_odn {
  tag "Trim ODN for ${sample_id}"
  
  input:
  tuple val(sample_id), path(r1), path(r2)
  
  output:
  tuple val(sample_id), path(r1_out_odn), path(r2_out_odn)
  
  script:
  """
  cutadapt -j ${task.cpus} \
  -G "negative=${params.primer.negative.R2_leading};max_error_rate=0;rightmost" \
  -G "positive=${params.primer.positive.R2_leading};max_error_rate=0;rightmost" \
  --discard-untrimmed \
  --rename='{id}_{r2.adapter_name} {comment}' \
  -o r1_out_odn -p r2_out_odn $r1 $r2 > trim_odn.log
  """
}

// Trim adapters from the 3' regions of the R1 and R2 reads
// for read R2, trim the R1 primer binding region (if present)
// for read R1, trim the dsODN tag (if present)
process trim_adapters {
  tag "Trim adapters for ${sample_id}"
  
  input:
  tuple val(sample_id), path(r1), path(r2)
  
  output:
  tuple val(sample_id), path(r1_out_adapter_trim), path(r2_out_adapter_trim)
  
  script:
  """
  cutadapt -j ${task.cpus} \
  -A "${params.primer.positive.R2_trailing};min_overlap=6;max_error_rate=0.1" \
  -A "${params.primer.negative.R2_trailing};min_overlap=6;max_error_rate=0.1" \
  -a "${params.primer.positive.R1_trailing};min_overlap=6;max_error_rate=0.1" \
  -a "${params.primer.negative.R1_trailing};min_overlap=6;max_error_rate=0.1" \
  -o r1_out_adapter_trim -p r2_out_adapter_trim $r1 $r2 > adapter_trim.log
  """
}

// filter for reads that are sufficiently long (i.e., both r1 and r2 are greater than or equal to 25 bp in length)
process filter_reads {
  tag "Filter reads on length for ${sample_id}"
  
  input:
  tuple val(sample_id), path(r1), path(r2)
  
  output:
  tuple val(sample_id), path(r1_length_filtered), path(r2_length_filtered), emit: filter_good
  tuple val(sample_id), path(r2_too_short), emit: filter_bad
  
  script:
  """
  cutadapt -j ${task.cpus} \
  --pair-filter=any \
  --minimum-length ${params.min_length} \
  --too-short-output r1_too_short \
  --too-short-paired-output r2_too_short \
  -o r1_length_filtered \
  -p r2_length_filtered \
  $r1 $r2 > filter_length.log
  """
}


// Align paired end reads to genome
process align_to_genome {
    tag "Align ${sample_id} to genome"
  
    input:
    tuple val(sample_id), path(r1), path(r2)
    path index
    
    output:
    tuple val(sample_id), path("alignment.sam"), emit: alignment
    tuple val(sample_id), path("unaligned_reads.R1.fastq.gz"), path("unaligned_reads.R2.fastq.gz"), emit: unmapped
    
    script:
    """
    bowtie2 -p ${task.cpus} --no-unal \
    -I ${params.min_frag_length} \
    -X ${params.max_frag_length} \
    --dovetail \
    --no-mixed \
    --no-discordant \
    --un-conc-gz unaligned_reads.R%.fastq.gz \
    -x $index \
    alignment.sam -1 ${r1} -2 ${r2} 2> alignment.log
    """
}

// WORKFLOW
workflow {
// Set up fastq file channels
  Channel.fromPath(params.samplesheet)
      .splitCsv(header: true, sep: '\t') // Read the TSV file
      .map { row ->
          def sample_id = row.sample_id
          def r1 = file(row.R1, checkIfExists: true)
          def r2 = file(row.R2, checkIfExists: true)
          def i2 = file(row.I2, checkIfExists: true)
          return [ sample_id, r1, r2, i2 ]
      }
      .set { ch_input_reads }
  ch_extract_umi = extract_umi(ch_input_reads)
  ch_trim_odn = trim_odn(ch_extract_umi)
  ch_trim_adapter = trim_adapters(ch_trim_odn)
  ch_filter_reads = filter_reads(ch_trim_adapter)
  // ch_filter_reads.filter_bad
  align_to_genome(ch_filter_reads.filter_good, params.index)
}
