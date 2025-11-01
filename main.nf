/*
* PART I: BASIC READ PROCESSING
*/

process initial_read_processing {
  tag "Run initial read processing for ${sample_id}"
  
  input:
  tuple val(sample_id), path("r1"), path("r2"), path("i2")
  
  output:
  tuple val(sample_id), path("r1_length_filtered"), path("r2_length_filtered"), emit: filter_good
  tuple val(sample_id), path("r2_too_short"), emit: filter_short
  
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
  
  
  """
}



// Move the UMI from I2 to the R1/R2 read ID
process extract_umi {
  tag "Extract UMI for ${sample_id}"
  
  output:
  tuple val(sample_id), path("r1_out_umi"), path("r2_out_umi")
  
  input:
  tuple val(sample_id), path("r1"), path("r2"), path("i2")
  
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
  tuple val(sample_id), path("r1"), path("r2")
  
  output:
  tuple val(sample_id), path("r1_out_odn"), path("r2_out_odn")
  
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
  tuple val(sample_id), path("r1"), path("r2")
  
  output:
  tuple val(sample_id), path("r1_out_adapter_trim"), path("r2_out_adapter_trim")
  
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
  tuple val(sample_id), path("r1"), path("r2")
  
  output:
  tuple val(sample_id), path("r1_length_filtered"), path("r2_length_filtered"), emit: filter_good
  tuple val(sample_id), path("r2_too_short"), emit: filter_short
  
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
    cpus 2
  
    input:
    tuple val(sample_id), path("r1"), path("r2")
    
    output:
    tuple val(sample_id), path("alignment.sam"), emit: alignment
    tuple val(sample_id), path("unaligned_reads.R2.fastq"), emit: unmapped
    
    script:
    """
    bowtie2 -p ${task.cpus} --no-unal \
    -I ${params.min_frag_length} \
    -X ${params.max_frag_length} \
    --dovetail \
    --no-mixed \
    --no-discordant \
    --un-conc unaligned_reads.R%.fastq \
    -x ${params.index} \
    -S alignment.sam -1 ${r1} -2 ${r2} 2> alignment.log
    """
}

// Filter reads. We retain reads for which (1) the mapq is sufficiently high, or (2) the read mapped to 2+ locations with 0 or 1 mismatches. (In the latter case, choose randomly.)

// (From https://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html#bt2expt)
// In bowtie2, the best a true multiread (AS=XS) can get is MAPQ=1 regardless of how low or high its multiplicity.
// This occurs when there are 0 or 1 mismatches over perfect base calls in the read, or when AS=XS goes down to -6.
// When there are 2-5 mismatches over perfect base calls (or the AS=XS <= -12 ---- i.e. -12 to -30.6), the MAPQ becomes 0.
process filter_alignments {
  tag "Filter alignment for sample ${sample_id}"
  
  input:
  tuple val(sample_id), path("alignment")
  
  output:
  tuple val(sample_id), path("alignment"), path("filtered_reads.txt")
  
  script:
  """
  samtools view ${alignment} -e '((mapq == 1 && [AS] == [XS]) || (mapq >=${params.min_mapq}))' | cut -f1 | sort | uniq  > filtered_reads.txt
  """
}

// Sort and index alignments for downstream processing
process sort_alignments {
  tag "Sort alignments for sample ${sample_id}"
  
  input:
  tuple val(sample_id), path("alignment"), path("filtered_reads")
  
  output:
  tuple val(sample_id), path("output.bamName")
  
  script:
  """
  samtools sort -@ ${task.cpus} ${alignment} > output.bamPos
  samtools view -hb --qname-file ${filtered_reads} output.bamPos > output.bam
  samtools index output.bam
  samtools sort -n -@ ${task.cpus} output.bam > output.bamName
  """
}

// call the integration site (the first base of the R2 read)
process call_integration_sites {
  tag "Call integration sites for sample ${sample_id}"
  
  input:
  tuple val(sample_id), path("input_bam")
  
  output:
  tuple val(sample_id), path("output.umi")
  
  """
  bedtools bamtobed -bedpe -mate1 -i ${input_bam} > output.tmp
  awk 'BEGIN{OFS="\\t";FS="\\t"} (\$1 == \$4) {split(\$7,a,"_"); if(\$10=="+") print \$4,\$5,\$5,a[1],a[2],a[3],\$8,\$10,\$3-\$5; else print \$4,\$6-1,\$6-1,a[1],a[2],a[3],\$8,\$10,\$6-\$2}' output.tmp | sort -k1,1 -k2,3n -k6,6 -k5,5 -k8,8 | bedtools groupby -g 1,2,3,6,5,8 -c 4,7 -o count_distinct,mean  > output.umi
  """
}

// align R2 reads that (i) were filtered out because they were too short or (ii) failed to align as a pair
process rescue_r2_reads {
  tag "Rescue R2 reads for sample ${sample_id}"
  cpus 2
  
  input:
  tuple val(sample_id), path(r2_short)
  tuple val(sample_id), path(r2_unmapped)
  
  output:
  tuple val(sample_id), path("r2_alignment.sam")
  
  script:
  """
  bowtie2 -p ${task.cpus} --no-unal \
  -x ${params.index} \
  -U ${r2_short},${r2_unmapped} -S r2_alignment.sam 2> alignment.log
  """
}

// sort the aligned R2 reads
process sort_rescued_r2_reads {
  tag "Sort rescued R2 reads for sample ${sample_id}"
  
  input:
  tuple val(sample_id), path("r2_alignment")
  
  output:
  tuple val(sample_id), path("r2_alignment_sorted")
  
  shell:
  """
  samtools view -b -h  ${r2_alignment} -e '((mapq == 1 && [AS] == [XS]) || (mapq >=${params.min_mapq}))' | samtools sort - > r2_alignment_sorted
  samtools index r2_alignment_sorted
  """
}

// call integration sites for the rescued R2 reads
process call_integration_sites_r2_reads {
  tag "Call integration sites for rescued R2 reads for sample ${sample_id}"
  
  input:
  tuple val(sample_id), path("r2_alignment")
  
  output:
  tuple val(sample_id), path("output_r2.umi")
  
  """
  bedtools bamtobed -i ${r2_alignment} > output.tmp
  awk 'BEGIN{OFS="\\t";FS="\\t"} {split(\$4,a,"_"); if(\$6=="+") print \$1,\$2,\$2,a[1],a[2],a[3],\$5,\$6,\$3-\$2; else print \$1,\$3-1,\$3-1,a[1],a[2],a[3],\$5,\$6,\$3-\$2}' output.tmp | sort -k1,1 -k2,3n -k6,6 -k5,5 -k8,8 | bedtools groupby -g 1,2,3,6,5,8 -c 4,7 -o count_distinct,mean > output_r2.umi
  """
}

// produce combined data frame
// combine paired-end and r2 read data frames; collapse by UMI
process produce_combined_collapsed_data_frame {
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
  ch_extract_umi = extract_umi(ch_input_reads)
  ch_trim_odn = trim_odn(ch_extract_umi)
  ch_trim_adapter = trim_adapters(ch_trim_odn)
  ch_filter_reads = filter_reads(ch_trim_adapter)
  // MAIN LOGIC: align good paired end reads
  ch_align_genome = align_to_genome(ch_filter_reads.filter_good)
  ch_filter_alignments = filter_alignments(ch_align_genome.alignment)
  ch_sort_alignments = sort_alignments(ch_filter_alignments)
  ch_call_integration_sites = call_integration_sites(ch_sort_alignments)
  // PARALLEL RESCUE LOGIC: try to rescue R2 reads
  ch_rescue_r2_reads = rescue_r2_reads(ch_filter_reads.filter_short, ch_align_genome.unmapped)
  ch_sort_rescued_r2_reads = sort_rescued_r2_reads(ch_rescue_r2_reads)
  ch_call_integration_sites_r2_reads = call_integration_sites_r2_reads(ch_sort_rescued_r2_reads)
  // PRODUCE FINAL UMI-COLLAPSED DATA FRAME
  joined_ch = ch_call_integration_sites.join(ch_call_integration_sites_r2_reads)
  produce_combined_collapsed_data_frame(joined_ch)
}
