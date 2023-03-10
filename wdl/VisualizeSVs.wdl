version 1.0


#
workflow VisualizeSVs {
    input {
        String remote_vcf_dir
        String chr
        Int pos_from
        Int pos_to
        File show_svs_in_intervals
    }
    parameter_meta {
        show_svs_in_intervals: "File containing a possibly empty list of intervals (as in bcftools view). The program prints to STDOUT all the SVs that intersect each interval, from every file."
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
    }
}


task VisualizeSVsImpl {
    input {
        String remote_vcf_dir
        String chr
        Int pos_from
        Int pos_to
        File show_svs_in_intervals
    }
    
    Int ram_gb = 64  # Arbitrary
    Int ram_gb_effective = ram_gb - 4
    String docker_dir = "/sv-merging"
    String work_dir = "/cromwell_root/sv-merging"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        # Downloading and sorting VCFs
        while : ; do
            TEST=$(gsutil -m cp ~{remote_vcf_dir}/merged.1.vcf.gz ~{remote_vcf_dir}/merged.2.vcf ~{remote_vcf_dir}/merged.3.vcf ~{remote_vcf_dir}/merged.4.vcf ~{remote_vcf_dir}/merged.5.vcf \
                ~{remote_vcf_dir}/survivor.vcf \
                ~{remote_vcf_dir}/jasmine.vcf \
                ~{remote_vcf_dir}/hg38.sorted.fa.out.cleaned.orderedByChromosome ~{remote_vcf_dir}/grch38_noalt-human_GRCh38_no_alt_analysis_set.trf.bed.orderedByChromosome \
                . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading files. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        gunzip merged.1.vcf.gz
        for i in $(seq 1 5); do
            bcftools sort --output-type v --output merged.${i}.sorted.vcf merged.${i}.vcf
        done
        bcftools sort --output-type v --output survivor.sorted.vcf survivor.vcf
        bcftools sort --output-type v --output jasmine.sorted.vcf jasmine.vcf
        
        # Visualizing SVs
        OUTPUT_FILE="mergedSVs_~{chr}_~{pos_from}_~{pos_to}.png"
        REPEAT_MASKER_FILE="hg38.sorted.fa.out.cleaned.orderedByChromosome"
        TRF_FILE="grch38_noalt-human_GRCh38_no_alt_analysis_set.trf.bed.orderedByChromosome"
        ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx~{ram_gb_effective}G PrintPopulationSVs \
            ~{chr} \
            ~{pos_from} \
            ~{pos_to} \
            ${REPEAT_MASKER_FILE} \
            $(wc -l < ${REPEAT_MASKER_FILE}) \
            ${TRF_FILE} \
            $(wc -l < ${TRF_FILE}) \
            merged.1.sorted.vcf merged.2.sorted.vcf merged.3.sorted.vcf merged.4.sorted.vcf merged.5.sorted.vcf \
            survivor.sorted.vcf \
            null \
            jasmine.sorted.vcf \
            ${OUTPUT_FILE}
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ${OUTPUT_FILE} ~{remote_vcf_dir} && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading image. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        
        # Printing SVs that start inside the given intervals (if any).
        if [ $(wc -l < ~{show_svs_in_intervals}) -gt 0 ]; then
            for i in $(seq 1 5); do
                bgzip --threads ${N_THREADS} merged.${i}.sorted.vcf
                tabix merged.${i}.sorted.vcf.gz
            done
            while read REGION; do
                for i in $(seq 1 5); do
                    bcftools view --no-header --regions ${REGION} --output-type v merged.${i}.sorted.vcf.gz > merged.${i}.in.${REGION}.txt
                done
            done < ~{show_svs_in_intervals}
            
            bgzip --threads ${N_THREADS} survivor.sorted.vcf
            tabix survivor.sorted.vcf.gz
            while read REGION; do
                bcftools view --no-header --regions ${REGION} --output-type v survivor.sorted.vcf.gz > survivor.in.${REGION}.txt
            done < ~{show_svs_in_intervals}
            
            bgzip --threads ${N_THREADS} jasmine.sorted.vcf
            tabix jasmine.sorted.vcf.gz
            while read REGION; do
                bcftools view --no-header --regions ${REGION} --output-type v jasmine.sorted.vcf.gz > jasmine.in.${REGION}.txt
            done < ~{show_svs_in_intervals}
            
            while : ; do
                TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp "*.in.*.txt" ~{remote_vcf_dir} && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading query intervals. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        fi
    >>>
    output {
    }
    runtime {
        docker: "fcunial/sv-merging"
        cpu: 1  # bcftools and truvari are sequential
        memory: ram_gb + "GB"  # Arbitrary
        disks: "local-disk 250 HDD"  # Arbitrary
        preemptible: 0
    }
}
