/*
* PART I: BASIC READ PROCESSING
*/

process run_initial_read_processing {
  tag "Running initial read processing and alignment for ${sample_id}"
  
  cpus 4
  memory 6.GB
  time 20.m
  
  input:
  tuple val(sample_id), path("r1"), path("r2"), path("i2")
  
  output:
  tuple val(sample_id), path("paired_alignment.sam"), emit: aligned
  tuple val(sample_id), path("r2_too_short"), path("unaligned_reads.R2.fastq"), emit: bad_r2
  
  script:
  """
  # step A
  
  if [ ${params.umi_side} = "5" ]; then
      cutadapt -j ${task.cpus} -u ${params.umi_length} --rename='{id}_{r1.cut_prefix} {comment}' -o i2_out_umi -p r1_out_umi $i2 $r1 > extract_umi_r1.log
      cutadapt -j ${task.cpus} -u ${params.umi_length} --rename='{id}_{r1.cut_prefix} {comment}' -o i2_out_umi -p r2_out_umi $i2 $r2 > extract_umi_r2.log
  else
      cutadapt -j ${task.cpus} -u -${params.umi_length} --rename='{id}_{r1.cut_suffix} {comment}' -o i2_out_umi -p r1_out_umi $i2 $r1 > extract_umi_r1.log
      cutadapt -j ${task.cpus} -u -${params.umi_length} --rename='{id}_{r1.cut_suffix} {comment}' -o i2_out_umi -p r2_out_umi $i2 $r2 > extract_umi_r2.log
  fi
  
  # step B
  
  cutadapt -j ${task.cpus} \
  -G "negative=${params.primer.negative.R2_leading};max_error_rate=0;rightmost" \
  -G "positive=${params.primer.positive.R2_leading};max_error_rate=0;rightmost" \
  --discard-untrimmed \
  --rename='{id}_{r2.adapter_name} {comment}' \
  -o r1_out_odn -p r2_out_odn r1_out_umi r2_out_umi > trim_odn.log
  
  # step C
  
  cutadapt -j ${task.cpus} \
  -A "${params.primer.positive.R2_trailing};min_overlap=6;max_error_rate=0.1" \
  -A "${params.primer.negative.R2_trailing};min_overlap=6;max_error_rate=0.1" \
  -a "${params.primer.positive.R1_trailing};min_overlap=6;max_error_rate=0.1" \
  -a "${params.primer.negative.R1_trailing};min_overlap=6;max_error_rate=0.1" \
  -o r1_out_adapter_trim -p r2_out_adapter_trim r1_out_odn r2_out_odn > adapter_trim.log
  
  # step D
  
  cutadapt -j ${task.cpus} \
  --pair-filter=any \
  --minimum-length ${params.min_length} \
  --too-short-output r1_too_short \
  --too-short-paired-output r2_too_short \
  -o r1_length_filtered \
  -p r2_length_filtered \
  r1_out_adapter_trim r2_out_adapter_trim > filter_length.log
  
  # step E
  
  bowtie2 -p ${task.cpus} --no-unal \
  -I ${params.min_frag_length} \
  -X ${params.max_frag_length} \
  --dovetail \
  --no-mixed \
  --no-discordant \
  --un-conc unaligned_reads.R%.fastq \
  -x ${params.index} \
  -S paired_alignment.sam -1 r1_length_filtered -2 r2_length_filtered 2> alignment.log
  """
}


process process_paired_end_alignments {
  tag "Processing paired end alignments for sample ${sample_id}"
  
  cpus 1
  memory 6.GB
  time 5.m
  
  input:
  tuple val(sample_id), path("alignment")
  
  output:
  tuple val(sample_id), path("output.umi")

  """
  # step A
  
  samtools view ${alignment} -e '((mapq == 1 && [AS] == [XS]) || (mapq >=${params.min_mapq}))' | cut -f1 | sort | uniq  > filtered_reads.txt

  # step B
  
  samtools sort -@ ${task.cpus} ${alignment} > output.bamPos
  samtools view -hb --qname-file filtered_reads.txt output.bamPos > output.bam
  samtools index output.bam
  
  # step C
  
  bedtools bamtobed -bedpe -mate1 -i output.bam > output.tmp
  awk 'BEGIN{OFS="\\t";FS="\\t"} (\$1 == \$4) {split(\$7,a,"_"); if(\$10=="+") print \$4,\$5,\$5,a[1],a[2],a[3],\$8,\$10,\$3-\$5; else print \$4,\$6-1,\$6-1,a[1],a[2],a[3],\$8,\$10,\$6-\$2}' output.tmp | sort -k1,1 -k2,3n -k6,6 -k5,5 -k8,8 | bedtools groupby -g 1,2,3,6,5,8 -c 4,7 -o count_distinct,mean  > output.umi
  """
}


process rescue_r2_reads {
  tag "Rescuing R2 reads for sample ${sample_id}"
    
  cpus 4
  memory 6.GB
  time 20.m
  
  input:
  tuple val(sample_id), path("r2_short"), path("r2_unmapped")
  
  output:
  tuple val(sample_id), path("output_r2.umi")
  
  
  script:
  """
  # step A
  
  bowtie2 -p ${task.cpus} --no-unal \
  -x ${params.index} \
  -U ${r2_short},${r2_unmapped} -S r2_alignment.sam 2> alignment.log
  
  # step B
  
  samtools view -b -h r2_alignment.sam -e '((mapq == 1 && [AS] == [XS]) || (mapq >=${params.min_mapq}))' | samtools sort - > r2_alignment_sorted
  samtools index r2_alignment_sorted
  
  # step C
  
  bedtools bamtobed -i r2_alignment_sorted > output.tmp
  awk 'BEGIN{OFS="\\t";FS="\\t"} {split(\$4,a,"_"); if(\$6=="+") print \$1,\$2,\$2,a[1],a[2],a[3],\$5,\$6,\$3-\$2; else print \$1,\$3-1,\$3-1,a[1],a[2],a[3],\$5,\$6,\$3-\$2}' output.tmp | sort -k1,1 -k2,3n -k6,6 -k5,5 -k8,8 | bedtools groupby -g 1,2,3,6,5,8 -c 4,7 -o count_distinct,mean > output_r2.umi
  """
}


// produce combined data frame
// combine paired-end and r2 read data frames; collapse by UMI
process produce_combined_collapsed_data_frame {
  cpus 1
  memory 6.GB
  time 5.m
  
  tag "Producing combined/collapsed data frame from sample ${sample_id}"
  
  input:
  tuple val(sample_id), path("paired_df"), path("r2_df")

  output:
  tuple val(sample_id), path("${sample_id}_count_table.rds")
  
  script:
  """
  collapse_umis.R ${paired_df} ${r2_df} ${sample_id}_count_table.rds
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
  ch_run_initial_read_processing = run_initial_read_processing(ch_input_reads)
  ch_process_paired_end_alignments = process_paired_end_alignments(ch_run_initial_read_processing.aligned)
  ch_rescue_r2_reads = rescue_r2_reads(ch_run_initial_read_processing.bad_r2)
  joined_ch = ch_process_paired_end_alignments.join(ch_rescue_r2_reads)
  produce_combined_collapsed_data_frame(joined_ch).view()
}
