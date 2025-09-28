# 1. Fastq to BAM conversion
bwa mem -K 100000000 -v 3 -Y -R "${RG}" -t ${nThreads} ${ref_fasta} ${sample_fq_1} ${sample_fq_2} 2>${sample_tag}.bwa.log | samtools view -1 -b - -o ${sample_tag}.bwa.bam

# 2. Sort and mark duplicates
gatk --java-options -Xmx15000m SortSam --INPUT ${bam_raw} --OUTPUT ${bam_sorted} --SORT_ORDER "coordinate"

gatk --java-options -Xmx15000m MarkDuplicates --INPUT ${bam_sorted} --OUTPUT ${bam_markdup} --METRICS_FILE ${txt_metrics} \
			--CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --ASSUME_SORT_ORDER "coordinate"


# 3. Base quality score recalibration
gatk --java-options -Xmx12000m BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --use-original-qualities \
      -O ${recalibration_report_filename} \
      --known-sites ${dbSNP_vcf} \
      --known-sites ${known_indels_sites_VCF}

# 4. Apply BQSR

gatk --java-options -Xmx15000m ApplyBQSR \
	      	-R ${ref_fasta} \
	      	-I ${input_bam} \
	      	-O ${sample_tag}.bam \
	      	-bqsr ${recalibration_report} \
	      	--add-output-sam-program-record \
	      	--create-output-bam-md5 \
	      	--use-original-qualities
