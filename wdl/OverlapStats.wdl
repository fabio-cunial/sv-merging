version 1.0


#
workflow OverlapStats {
    input {
        String bucket_dir
        String merger
        String genotyper
        Int n_samples
        Int max_overlap
    }
    parameter_meta {
    }
    call OverlapStatsImpl {
        input:
            bucket_dir = bucket_dir,
            merger = merger,
            genotyper = genotyper,
            n_samples = n_samples,
            max_overlap = max_overlap
    }
    output {
        File report = OverlapStatsImpl.report
    }
}


task OverlapStatsImpl {
    input {
        String bucket_dir
        String merger
        String genotyper
        Int n_samples
        Int max_overlap
    }
    parameter_meta {
        genotyper: "null: no regenotyping."
    }
    
    String docker_dir = "/sv-merging"
    String work_dir = "/cromwell_root/sv-merging"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        seq 0 $(( ~{n_samples}-1 )) > list.txt
        N_ROWS=$(( ~{n_samples}/${N_THREADS} ))
        split -d -l ${N_ROWS} list.txt chunk-
        if [ ~{genotyper} != "null" ]; then
            ID="regenotyped_~{genotyper}_standardized"
        else
            ID="standardized"
        fi    
        VCF_GZ_FILE="${ID}.vcf.gz"
        gsutil cp ~{bucket_dir}/~{merger}/${VCF_GZ_FILE} .
        gunzip ${VCF_GZ_FILE}
        for SAMPLES_FILE in $(ls chunk-*); do
            java -Xmx2g -cp ~{docker_dir} OverlapStats ${ID}.vcf ${SAMPLES_FILE} 1 . ~{max_overlap} &
        done
        wait
        if [ ~{genotyper} != "null" ]; then
            INFIX="_~{genotyper}"
        else
            INFIX=""
        fi
        cat *_chr21_del_histogram.txt > ~{merger}${INFIX}_chr21_del_histogram.matrix
        cat *_chr21_dup_histogram.txt > ~{merger}${INFIX}_chr21_dup_histogram.matrix
        cat *_chr21_inv_histogram.txt > ~{merger}${INFIX}_chr21_inv_histogram.matrix
        cat *_chr22_del_histogram.txt > ~{merger}${INFIX}_chr22_del_histogram.matrix
        cat *_chr22_dup_histogram.txt > ~{merger}${INFIX}_chr22_dup_histogram.matrix
        cat *_chr22_inv_histogram.txt > ~{merger}${INFIX}_chr22_inv_histogram.matrix
        cat *_gtCounts.txt > ~{merger}${INFIX}_gtCounts.matrix
        rm -f *.txt ${ID}.vcf
        tar -czvf report.tar.gz *.matrix
    >>>
    output {
        File report = work_dir + "/report.tar.gz"
    }
    runtime {
        docker: "fcunial/sv-merging"
        cpu: 16  # Arbitrary
        memory: "32GB"  # Arbitrary
        disks: "local-disk 250 HDD"  # Arbitrary
        preemptible: 0
    }
}
