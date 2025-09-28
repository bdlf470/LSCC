# 1. CRAM to BAM conversion if needed
samtools view -h -T ${ref_fasta} ${cram} | samtools view -b -o ${name}.bam -
samtools index -b ${name}.bam

# 2. Split intervals
gatk --java-options -Xmx3000m SplitIntervals \
            -R ${ref_fasta} \
            ${"-L " + intervals} \
            -scatter ${scatter_count} \
            -O interval-files \
            ${split_intervals_extra_args}

# 3. Run Mutect2 on each interval
        touch bamout.bam
        touch f1r2.tar.gz
        echo "" > normal_name.txt

        gatk --java-options -Xmx8000m GetSampleName -R ${ref_fasta} -I ${tumor_bam} -O tumor_name.txt -encode
        tumor_command_line="-I ${tumor_bam} -tumor `cat tumor_name.txt`"

        if [[ ! -z "${normal_bam}" ]]; then
            gatk --java-options -Xmx8000m GetSampleName -R ${ref_fasta} -I ${normal_bam} -O normal_name.txt -encode
            normal_command_line="-I ${normal_bam} -normal `cat normal_name.txt`"
        fi

        gatk --java-options -Xmx8000m Mutect2 \
            -R ${ref_fasta} \
            $tumor_command_line \
            $normal_command_line \
            ${"--germline-resource " + gnomad} \
            ${"-pon " + pon} \
            ${"-L " + intervals} \
            ${"--alleles " + gga_vcf} \
            -O "${output_vcf}" \
            ${true='--bam-output bamout.bam' false='' make_bamout} \
            ${true='--f1r2-tar-gz f1r2.tar.gz' false='' run_ob_filter} \
            ${m2_extra_args}

        ### GetPileupSummaries
        # These must be created, even if they remain empty, as cromwell doesn't support optional output
        touch tumor-pileups.table
        touch normal-pileups.table

        if [[ ! -z "${variants_for_contamination}" ]]; then
            gatk --java-options -Xmx8000m GetPileupSummaries -R ${ref_fasta} -I ${tumor_bam} ${"--interval-set-rule INTERSECTION -L " + intervals} \
                -V ${variants_for_contamination} -L ${variants_for_contamination} -O tumor-pileups.table

            if [[ ! -z "${normal_bam}" ]]; then
                gatk --java-options -Xmx8000m GetPileupSummaries -R ${ref_fasta} -I ${normal_bam} ${"--interval-set-rule INTERSECTION -L " + intervals} \
                    -V ${variants_for_contamination} -L ${variants_for_contamination} -O normal-pileups.table
            fi
        fi

# 4. Merge VCFs
        gatk --java-options -Xmx4000m MergeVcfs -I ${sep=' -I ' input_vcfs} -O ${output_vcf}

# 5. Merge stats
gatk --java-options -Xmx7000m GatherBamFiles \
            -I ${sep=" -I " bam_outs} -O unsorted.out.bam -R ${ref_fasta}

        # We must sort because adjacent scatters may have overlapping (padded) assembly regions, hence
        # overlapping bamouts

        gatk --java-options -Xmx7000m SortSam -I unsorted.out.bam \
            -O ${output_vcf_name}.out.bam \
            --SORT_ORDER coordinate -VALIDATION_STRINGENCY LENIENT
        gatk --java-options -Xmx7000m BuildBamIndex -I ${output_vcf_name}.out.bam -VALIDATION_STRINGENCY
        
        gatk --java-options -Xmx3000m MergeMutectStats -stats ${sep=" -stats " stats} -O merged.stats

# 6. Calculate contamination if variants provided
        gatk --java-options -Xmx4000m GatherPileupSummaries \
        --sequence-dictionary ${ref_dict} \
        -I ${sep=' -I ' input_tables} \
        -O ${output_name}.tsv

# 7. Learn read orientation model if requested
        gatk --java-options -Xmx8000m LearnReadOrientationModel -I ${sep=" -I " f1r2_tar_gz} -O "artifact-priors.tar.gz"
        
                gatk --java-options -Xmx3000m CalculateContamination -I ${tumor_pileups} \
        -O contamination.table --tumor-segmentation segments.table ${"-matched " + normal_pileups}

# 8. Filter variants
gatk --java-options -Xmx7000m FilterMutectCalls -V ${unfiltered_vcf} \
            -R ${ref_fasta} \
      	    -O ${output_vcf} \
      	    ${"--contamination-table " + contamination_table} \
      	    ${"--tumor-segmentation " + maf_segments} \
      	    ${"--ob-priors " + artifact_priors_tar_gz} \
      	    ${"-stats " + mutect_stats} \
      	    --filtering-stats filtering.stats \
      	    ${m2_extra_filtering_args}

# 9. Filter alignment artifacts if requested
gatk --java-options -Xmx9000m FilterAlignmentArtifacts \
            -V ${input_vcf} \
            -I ${bam} \
            --bwa-mem-index-image ${realignment_index_bundle} \
            ${realignment_extra_args} \
            -O ${output_vcf}

