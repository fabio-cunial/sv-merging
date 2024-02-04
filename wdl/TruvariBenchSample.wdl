# Truvari intra-merge benchmarking for AoU SV
version 1.0


# Workflow for intra-sample evaluation for AoU.
#
workflow TruvariBenchSample {
    input {
        String sample_id
        File baseline_variants
        File baseline_variants_tbi
        File comparison_variants
        File comparison_variants_tbi
        File includebed
        Boolean inter
    }
    parameter_meta {
    }
    
    call BenchImpl {
        input:
            sample_id = sample_id,
            baseline_variants = baseline_variants,
            baseline_variants_tbi = baseline_variants_tbi,
            comparison_variants = comparison_variants,
            comparison_variants_tbi = comparison_variants_tbi,
            includebed = includebed,
            inter = inter
    }
    
    output {
        File truvari_bench_summary = BenchImpl.bench_summary
        File bench_results = BenchImpl.bench_results
    }
}


task BenchImpl {
    input {
        String sample_id
        File baseline_variants
        File baseline_variants_tbi
        File comparison_variants
        File comparison_variants_tbi
        File includebed
        Boolean inter
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*( ceil(size(baseline_variants,"GB")) + ceil(size(comparison_variants,"GB")) ) + 50
    String docker_dir = "/truvari_intrasample"
    String work_dir = "/cromwell_root/truvari_intrasample"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        extra=""
        if [ "~{inter}" = true ]; then
           extra="--no-ref a"
        fi

        # Updating date of TBI files
        touch ~{baseline_variants_tbi} ~{comparison_variants_tbi}

        # Making sure that the baseline has the right sample name (e.g. this is
        # not true for dipcall).
        echo ~{sample_id} > samples.txt
        bcftools reheader --samples samples.txt ~{baseline_variants} > reheaded.vcf.gz
        tabix reheaded.vcf.gz
        
        # Making sure that the baseline has no multiallelic records (e.g. this
        # is not true for dipcall).
        bcftools norm --multiallelics - --output-type z reheaded.vcf.gz > baseline_fixed.vcf.gz
        tabix baseline_fixed.vcf.gz
        rm -f reheaded.vcf.gz*

        truvari bench -b baseline_fixed.vcf.gz \
           -c ~{comparison_variants} \
           --includebed ~{includebed} \
           --bSample ~{sample_id} \
           --cSample ~{sample_id} \
           --pick multi ${extra} \
           -o ./bench_results/
        tar czf bench_results.tar.gz ./bench_results/
    >>>
    
    output {
        File bench_summary = work_dir + "/bench_results/summary.json"
        File bench_results = work_dir + "/bench_results.tar.gz"
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
        cpu: 1
        memory: "16GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
