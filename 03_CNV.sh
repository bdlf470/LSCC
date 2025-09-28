

# 1. Collect read counts
        gatk --java-options -Xmx4g CollectReadCounts \
            -L ${intervals} \
            --input ${bam} \
            --reference ${ref_fasta} \
            --format HDF5 \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ${counts_filename}

# 2. Collect allelic counts
        gatk --java-options -Xmx4g CollectAllelicCounts \
            -L ${common_sites} \
            --input ${bam} \
            --reference ${ref_fasta} \
            --minimum-base-quality ${default="20" minimum_base_quality} \
            --output ${allelic_counts_filename}

# 3. Denoise read counts

        gatk --java-options -Xmx4g DenoiseReadCounts \
            --input ${read_counts} \
            --count-panel-of-normals ${read_count_pon} \
            ${"--number-of-eigensamples " + number_of_eigensamples} \
            --standardized-copy-ratios ${entity_id}.standardizedCR.tsv \
            --denoised-copy-ratios ${entity_id}.denoisedCR.tsv

# 4. Model segments

gatk --java-options -Xmx8g ModelSegments \
            --denoised-copy-ratios ${denoised_copy_ratios} \
            --allelic-counts ${allelic_counts} \
            ${"--normal-allelic-counts " + normal_allelic_counts} \
            --minimum-total-allele-count-case ${min_total_allele_count_} \
            --minimum-total-allele-count-normal ${default="30" min_total_allele_count_normal} \
            --genotyping-homozygous-log-ratio-threshold ${default="-10.0" genotyping_homozygous_log_ratio_threshold} \
            --genotyping-base-error-rate ${default="0.05" genotyping_base_error_rate} \
            --maximum-number-of-segments-per-chromosome ${default="1000" max_num_segments_per_chromosome} \
            --kernel-variance-copy-ratio ${default="0.0" kernel_variance_copy_ratio} \
            --kernel-variance-allele-fraction ${default="0.025" kernel_variance_allele_fraction} \
            --kernel-scaling-allele-fraction ${default="1.0" kernel_scaling_allele_fraction} \
            --kernel-approximation-dimension ${default="100" kernel_approximation_dimension} \
            --window-size ${sep=" --window-size " window_sizes} \
            --number-of-changepoints-penalty-factor ${default="1.0" num_changepoints_penalty_factor} \
            --minor-allele-fraction-prior-alpha ${default="25.0" minor_allele_fraction_prior_alpha} \
            --number-of-samples-copy-ratio ${default="100" num_samples_copy_ratio} \
            --number-of-burn-in-samples-copy-ratio ${default="50" num_burn_in_copy_ratio} \
            --number-of-samples-allele-fraction ${default="100" num_samples_allele_fraction} \
            --number-of-burn-in-samples-allele-fraction ${default="50" num_burn_in_allele_fraction} \
            --smoothing-credible-interval-threshold-copy-ratio ${default="2.0" smoothing_threshold_copy_ratio} \
            --smoothing-credible-interval-threshold-allele-fraction ${default="2.0" smoothing_threshold_allele_fraction} \
            --maximum-number-of-smoothing-iterations ${default="10" max_num_smoothing_iterations} \
            --number-of-smoothing-iterations-per-fit ${default="0" num_smoothing_iterations_per_fit} \
            --output ${output_dir_} \
            --output-prefix ${entity_id}


# 5. Call copy ratio segments
    gatk --java-options -Xmx8g CallCopyRatioSegments \
            --input ${copy_ratio_segments} \
            --neutral-segment-copy-ratio-lower-bound ${default="0.9" neutral_segment_copy_ratio_lower_bound} \
            --neutral-segment-copy-ratio-upper-bound ${default="1.1" neutral_segment_copy_ratio_upper_bound} \
            --outlier-neutral-segment-copy-ratio-z-score-threshold ${default="2.0" outlier_neutral_segment_copy_ratio_z_score_threshold} \
            --calling-copy-ratio-z-score-threshold ${default="2.0" calling_copy_ratio_z_score_threshold} \
            --output ${entity_id}.called.seg

