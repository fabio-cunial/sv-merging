version 1.0


# regions = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"] 
#
workflow PhasedMerge {
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
        call MergePAV {
            input:
                vcf_addresses = vcf_addresses,
                region = region,
                output_dir = output_dir,
                reference_fa = reference_fa,
                reference_fai = reference_fai
        }
        call PangenieMerge {
            input: 
                vcf_gz = MergePAV.merged_vcf,
                reference_fa = reference_fa,
                reference_fai = reference_fai
        }
    }
    
    output {
    }
}


# Merges the phased PAV VCFs (which contain both SVs and SNPs) by merging SVs
# with $truvari collapse$ and SNPs with $bcftools merge$.
#
# Remark: the output VCF does not contain mulatiallelics, since we can use
# PanGenie's tools to create a multiallelic VCF later.
#
task MergePAV {
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
        MIN_SV_LENGTH="50"  # Arbitrary
        
        # Downloading and sorting all VCFs
        touch list_svs.txt list_snps.txt
        while read REMOTE_VCF_GZ_FILE; do
            LOCAL_FILE=$(basename -s .vcf.gz ${REMOTE_VCF_GZ_FILE})
            TEST=$(bcftools view --threads ${N_THREADS} --include "FILTER=\"PASS\" || FILTER=\".\"" --regions "~{region}" --output-type v ${REMOTE_VCF_GZ_FILE} > tmp.vcf  && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
                bcftools view --threads ${N_THREADS} --include "FILTER=\"PASS\" || FILTER=\".\"" --regions "~{region}" --output-type v ${REMOTE_VCF_GZ_FILE} > tmp.vcf
            fi
            java -cp ~{docker_dir} PAV2SVs tmp.vcf ${MIN_SV_LENGTH} svs.vcf snps.vcf
            rm -f tmp.vcf
            python3 /sv-merging/preprocess_vcf.py svs.vcf ~{reference_fa} > svs_cleaned.vcf
            rm -f svs.vcf
            bcftools sort --output-type z svs_cleaned.vcf > ${LOCAL_FILE}_svs.vcf.gz
            tabix -f ${LOCAL_FILE}_svs.vcf.gz
            rm -f svs_cleaned.vcf
            echo ${LOCAL_FILE}_svs.vcf.gz >> list_svs.txt
            bcftools sort --output-type z snps.vcf > ${LOCAL_FILE}_snps.vcf.gz
            tabix -f ${LOCAL_FILE}_snps.vcf.gz
            rm -f snps.vcf
            echo ${LOCAL_FILE}_snps.vcf.gz >> list_snps.txt
        done < ~{vcf_addresses}
        
        # Remark: $bcftools merge$ makes the following changes to GT, which we
        # have to reverse:
        # .|1 -> ./1
        # 1|. -> 1/.
        #
        # Remark: we convert every '.' to a 0, since PanGenie's multiallelic VCF
        # tool discards a haplotype if it contains even a single '.'.
        REPLACEMENT_COMMAND='s@\./1@0|1@g;s@1/\.@1|0@g;s@\./\.@0|0@g;s@\.|1@0|1@g;s@1|\.@1|0@g;s@\.|\.@0|0@g'
        
        # Merging SVs
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --file-list list_svs.txt | sed ${REPLACEMENT_COMMAND} | bgzip > bcftools_svs.vcf.gz
        tabix -f bcftools_svs.vcf.gz
        bcftools view --no-header bcftools_svs.vcf.gz | head -n 20 && echo 0 || echo 1
        ${TIME_COMMAND} truvari collapse --input bcftools_svs.vcf.gz --output truvari_collapse.vcf --reference ~{reference_fa}
        bcftools sort --output-type z truvari_collapse.vcf > truvari_collapse_sorted.vcf.gz
        tabix -f truvari_collapse_sorted.vcf.gz
        rm -f truvari_collapse.vcf
        bcftools view --no-header truvari_collapse_sorted.vcf.gz | head -n 20 && echo 0 || echo 1
        
        # Naive merging of SNPs
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --file-list list_snps.txt | sed ${REPLACEMENT_COMMAND} | bgzip > bcftools_snps.vcf.gz
        tabix -f bcftools_snps.vcf.gz
        bcftools view --no-header bcftools_snps.vcf.gz | head -n 20 && echo 0 || echo 1
        
        # Combining SVs and SNPs
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --output-type z --output final.vcf.gz truvari_collapse_sorted.vcf.gz bcftools_snps.vcf.gz
        tabix -f final.vcf.gz
        bcftools view --no-header final.vcf.gz | head -n 20 && echo 0 || echo 1
    >>>
    
    output {
        File merged_vcf = work_dir + "/final.vcf.gz"
        File merged_tbi = work_dir + "/final.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/sv-merging"
        cpu: 16  # Arbitrary
        memory: "64GB"  # Arbitrary
        disks: "local-disk 128 HDD"  # Arbitrary
        preemptible: 0
    }
}


# Creates a multiallelic VCF representing PanGenie's bubbles.
#
task PangenieMerge {
    input {
        File vcf_gz
        File reference_fa
        File reference_fai
    }
    parameter_meta {
    }
    
    Int disk_size_gb = ceil(size(vcf_gz, "GB"))*20 + ceil(size(reference_fa, "GB")) + 50
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

        gunzip -c ~{vcf_gz} > input.vcf
        cp -r ~{docker_dir}/vcf-merging/pangenome-graph-from-callset/ ./
        cd pangenome-graph-from-callset/
        rm -f config.yaml
        echo "vcf: ~{work_dir}/input.vcf" >> config.yaml
        echo "reference: ~{reference_fa}" >> config.yaml
        echo "outdir: ~{work_dir}/results" >> config.yaml
        source activate svpop
        ${TIME_COMMAND} snakemake -j ${N_THREADS}
        conda deactivate
        cd ~{work_dir}
        bgzip -@ ${N_THREADS} ./results/pangenome/pangenome.vcf
        tabix -f ./results/pangenome/pangenome.vcf.gz
    >>>
    
    output {
        File pangenome_vcf = work_dir + "/results/pangenome/pangenome.vcf.gz"
        File pangenome_tbi = work_dir + "/results/pangenome/pangenome.vcf.gz.tbi"
        File log_prepare_vcf = work_dir + "/results/input-vcf/prepare-vcf.log"
        File log_callset = work_dir + "/results/input-vcf/callset.log"
        File log_callset_biallelic = work_dir + "/results/input-vcf/callset-biallelic.log"
        File log_pangenome = work_dir + "/results/pangenome/pangenome.log"
    }
    runtime {
        docker: "fcunial/sv-merging"
        cpu: 8  # Arbitrary
        memory: "64GB"  # Arbitrary
        disks: "local-disk " + disk_size_gb + " HDD"  # Arbitrary
        preemptible: 0
    }
}