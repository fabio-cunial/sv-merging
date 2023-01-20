version 1.0


#
workflow GetVCFs {
    input {
        File vcf_addresses
        String region
        String output_dir
        File reference_fa
        File reference_fai
        Int n_cpus
    }
    parameter_meta {
        vcf_addresses: "File containing a list of .vcf.gz bucket addresses."
        region: "Region of the genome to be extracted from each VCF (as in bcftools view)."
        output_dir: "The output is stored in this destination directory in a bucket."
    }
    call GetVCFsImpl {
        input:
            vcf_addresses = vcf_addresses,
            region = region,
            output_dir = output_dir,
            reference_fa = reference_fa,
            reference_fai = reference_fai
    }
    output {
    }
}


task GetVCFsImpl {
    input {
        File vcf_addresses
        String region
        String output_dir
        File reference_fa
        File reference_fai
        Int n_cpus
    }
    parameter_meta {
        vcf_addresses: "File containing a list of .vcf.gz bucket addresses."
        region: "Region of the genome to be extracted from each VCF (as in bcftools view)."
        output_dir: "The output is stored in this destination directory in a bucket."
    }
    
    String docker_dir = "/sv-merging"
    String work_dir = "/cromwell_root/sv-merging"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        
        # Downloading VCF regions
        touch list.txt
        while read VCF_FILE; do
            while : ; do
                TEST=$(gsutil cp ${VCF_FILE} . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <${VCF_FILE}>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            LOCAL_FILE=$(basename -s .vcf.gz ${VCF_FILE})
            bcftools view --threads ${N_THREADS} --regions "~{region}" --output-type z ${LOCAL_FILE}.vcf.gz | bcftools sort --output-type z --output ${LOCAL_FILE}.region.vcf.gz
            rm -f ${LOCAL_FILE}.vcf.gz
            tabix ${LOCAL_FILE}.region.vcf.gz
            echo ${LOCAL_FILE}.region.vcf.gz >> list.txt
        done < ~{vcf_addresses}
        
        # Merging SVs
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --file-list list.txt --output-type v --output merged.1.vcf
        ${TIME_COMMAND} truvari collapse --input merged.1.vcf --output merged.2.vcf --collapsed-output collapsed.2.vcf --reference ~{reference_fa} --keep first --passonly
        ${TIME_COMMAND} truvari collapse --input merged.1.vcf --output merged.3.vcf --collapsed-output collapsed.3.vcf --reference ~{reference_fa} --keep maxqual --passonly
        ${TIME_COMMAND} truvari collapse --input merged.1.vcf --output merged.4.vcf --collapsed-output collapsed.4.vcf --reference ~{reference_fa} --keep common --passonly
        ${TIME_COMMAND} truvari collapse --input merged.1.vcf --output merged.5.vcf --collapsed-output collapsed.5.vcf --reference ~{reference_fa} --keep common --chain --passonly
        
        # Outputting
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp "merged.*.vcf" "collapsed.*.vcf" ~{output_dir} && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading file <merged.*.vcf>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>
    output {
    }
    runtime {
        docker: "fcunial/sv-merging"
        cpu: n_cpus
        memory: "64GB"  # Arbitrary
        disks: "local-disk 100 HDD"  # Arbitrary
        preemptible: 0
    }
}
