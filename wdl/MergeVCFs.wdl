version 1.0


#
workflow MergeVCFs {
    input {
        File vcf_addresses
        Array[String] regions
        File reference_fa
        File reference_fai
        String output_dir
    }
    parameter_meta {
        vcf_addresses: "File containing a list of .vcf.gz bucket addresses."
        region: "Region of the genome to be extracted from each VCF (as in bcftools view)."
        output_dir: "The output is stored in this destination directory in a bucket."
    }
    
    scatter (region in regions) {
        call MergeVCFsImpl {
            input:
                vcf_addresses = vcf_addresses,
                region = region,
                output_dir = output_dir,
                reference_fa = reference_fa,
                reference_fai = reference_fai
        }
    }
    
    output {
    }
}


task MergeVCFsImpl {
    input {
        File vcf_addresses
        String region
        String output_dir
        File reference_fa
        File reference_fai
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
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
        TIME_COMMAND="/usr/bin/time --verbose"
        
        function uploadVCF() {
            local FILE1=$1
            local FILE2=$2
            
            while : ; do
                TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ${FILE1} ${FILE2} ~{output_dir} && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading file <${FILE1}>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        }
        
        # Downloading VCF region
        touch list.txt
        while read REMOTE_VCF_GZ_FILE; do
            LOCAL_FILE=$(basename -s .vcf.gz ${REMOTE_VCF_GZ_FILE})
            TEST=$(bcftools view --threads ${N_THREADS} --include "FILTER=\"PASS\" || FILTER=\".\"" --regions "~{region}" --output-type v ${REMOTE_VCF_GZ_FILE} > tmp1.vcf  && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
                bcftools view --threads ${N_THREADS} --include "FILTER=\"PASS\" || FILTER=\".\"" --regions "~{region}" --output-type v ${REMOTE_VCF_GZ_FILE} > tmp1.vcf
            fi
            python3 /sv-merging/preprocess_vcf.py tmp1.vcf ~{reference_fa} > tmp2.vcf
            rm -f tmp1.vcf
            bcftools sort --output-type v tmp2.vcf > ${LOCAL_FILE}.vcf
            rm -f tmp2.vcf
            echo ${LOCAL_FILE}.vcf >> list.txt
        done < ~{vcf_addresses}
        
        # JASMINE, default parameters.
        MERGED_VCF="jasmine_~{region}.vcf"
        TEST=$(gsutil -q stat ~{output_dir}/${MERGED_VCF}.gz && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            source activate jasmine
            ${TIME_COMMAND} jasmine --output_genotypes threads=${N_THREADS} file_list=list.txt out_file=${MERGED_VCF}
            conda deactivate
            bcftools sort --output-type z ${MERGED_VCF} > ${MERGED_VCF}.gz
            tabix ${MERGED_VCF}.gz
            uploadVCF ${MERGED_VCF}.gz ${MERGED_VCF}.gz.tbi
        fi
    >>>
    output {
    }
    runtime {
        docker: "fcunial/sv-merging"
        cpu: 8  # Arbitrary
        memory: "64GB"  # Arbitrary
        disks: "local-disk 128 HDD"  # Arbitrary
        preemptible: 0
    }
}
