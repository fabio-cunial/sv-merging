version 1.0


#
workflow VisualizeSVs {
    input {
        File remote_vcf_dir
        String chr
        Int pos_from
        Int pos_to
        File show_svs_in_intervals
    }
    parameter_meta {
        show_svs_in_intervals: "File containing a possibly empty list of intervals (as in bcftools view). The program prints to STDOUT all the SVs that start inside each interval, from every file."
    }
    call VisualizeSVsImpl {
        input:
            remote_vcf_dir = remote_vcf_dir,
            chr = chr,
            pos_from = pos_from,
            pos_to = pos_to,
            show_svs_in_intervals = show_svs_in_intervals
    }
    output {
        File log = VisualizeSVsImpl.log
    }
}


task VisualizeSVsImpl {
    input {
        File remote_vcf_dir
        String chr
        Int pos_from
        Int pos_to
        File show_svs_in_intervals
    }
    
    String docker_dir = "/sv-merging"
    String work_dir = "/cromwell_root/sv-merging"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        
        # Visualizing SVs
        while : ; do
            TEST=$(gsutil -m cp ~{remote_vcf_dir}/merged.1.vcf.gz ~{remote_vcf_dir}/merged.2.vcf ~{remote_vcf_dir}/merged.3.vcf ~{remote_vcf_dir}/merged.4.vcf ~{remote_vcf_dir}/merged.5.vcf ~{remote_vcf_dir}/hg38.sorted.fa.out.cleaned.orderedByChromosome  ~{remote_vcf_dir}/grch38_noalt-human_GRCh38_no_alt_analysis_set.trf.bed.orderedByChromosome . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading files. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        gunzip merged.1.vcf.gz
        OUTPUT_FILE="mergedSVs_~{chr}_~{pos_from}_~{pos_to}.png"
        REPEAT_MASKER_FILE="hg38.sorted.fa.out.cleaned.orderedByChromosome"
        TRF_FILE="grch38_noalt-human_GRCh38_no_alt_analysis_set.trf.bed.orderedByChromosome"
        ${TIME_COMMAND} java -cp ~{docker_dir} PrintPopulationSVs merged.1.vcf merged.2.vcf merged.3.vcf merged.4.vcf merged.5.vcf ${OUTPUT_FILE} ~{chr} ~{pos_from} ~{pos_to} ${REPEAT_MASKER_FILE} $(wc -l < ${REPEAT_MASKER_FILE}) ${TRF_FILE} $(wc -l ${TRF_FILE})
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ${OUTPUT_FILE} ~{remote_vcf_dir} && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading image. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        
        # Printing SVs that start in specified intervals (if any).
        if [ $(wc -l < ~{show_svs_in_intervals}) -gt 0 ]; then
            for i in $(seq 1 5); do
                bgzip merged.${i}.vcf
                tabix merged.${i}.vcf.gz
            done
            while read REGION; do
                for i in $(seq 1 5); do
                    echo ""
                    echo "==== SVs that start inside ${REGION} in file merged.${i}.vcf.gz:"
                    bcftools view --regions ${REGION} --output-type v merged.${i}.vcf.gz
                    echo ""
                done
            done < ~{show_svs_in_intervals}
        fi
    >>>
    output {
        File log = stdout()
    }
    runtime {
        docker: "fcunial/sv-merging"
        cpu: 1  # bcftools and truvari are sequential
        memory: "64GB"  # Arbitrary
        disks: "local-disk 250 HDD"  # Arbitrary
        preemptible: 0
    }
}
