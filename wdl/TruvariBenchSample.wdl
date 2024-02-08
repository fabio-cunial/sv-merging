# Truvari intra-merge benchmarking and region summary for AoU SV
version 1.0


# Workflow for intra-sample evaluation for AoU.
#
workflow TruvariBenchSample {
    input {
        String sample_id
        File truth_vcf
        File truth_vcf_idx
        File truth_bed
        File reference
        File bcftools_merged
        File bcftools_merged_idx
        File jasmine_default
        File jasmine_default_idx
        File svimmer
        File svimmer_idx
        File svmerger
        File svmerger_idx
        File truvari_collapsed
        File truvari_collapsed_idx
    }
    parameter_meta {
    }
    
    call BenchImpl {
        input:
            sample_id = sample_id,
            truth_vcf = truth_vcf,
            truth_vcf_idx = truth_vcf_idx,
            truth_bed = truth_bed,
            reference = reference,
            bcftools_merged = bcftools_merged,
            bcftools_merged_idx = bcftools_merged_idx,
            jasmine_default = jasmine_default,
            jasmine_default_idx = jasmine_default_idx,
            svimmer = svimmer,
            svimmer_idx = svimmer_idx,
            svmerger = svmerger,
            svmerger_idx = svmerger_idx,
            truvari_collapsed = truvari_collapsed,
            truvari_collapsed_idx = truvari_collapsed_idx
    }
    
    output {
        File bench_results = BenchImpl.bench_results
        File region_table = BenchImpl.region_table
        File region_summary = BenchImpl.region_summary
    }
}


task BenchImpl {
    input {
        String sample_id
        File truth_vcf
        File truth_vcf_idx
        File truth_bed
        File reference
        File bcftools_merged
        File bcftools_merged_idx
        File jasmine_default
        File jasmine_default_idx
        File svimmer
        File svimmer_idx
        File svmerger
        File svmerger_idx
        File truvari_collapsed
        File truvari_collapsed_idx
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*( ceil( size(truth_vcf,"GB")+size(bcftools_merged,"GB")+size(jasmine_default,"GB")+size(svimmer,"GB")+size(svmerger,"GB")+size(truvari_collapsed,"GB") ) ) + 50
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

        declare -A all_vcfs
        all_vcfs["baseline"]="~{truth_vcf}"
        all_vcfs["bcftools"]="~{bcftools_merged}"
        all_vcfs["jasmine"]="~{jasmine_default}"
        all_vcfs["svimmer"]="~{svimmer}"
        all_vcfs["svmerger"]="~{svmerger}"
        all_vcfs["truvari"]="~{truvari_collapsed}"

        # Making sure that the baseline has the right sample name (e.g. this is
        # not true for dipcall).
        echo ~{sample_id} > samples.txt
        bcftools reheader --samples samples.txt ${all_vcfs["baseline"]} > reheaded.vcf.gz
        tabix reheaded.vcf.gz

        # Making sure that the baseline has no multiallelic records (e.g. this
        # is not true for dipcall).
        bcftools norm --multiallelics - --output-type z reheaded.vcf.gz > baseline_fixed.vcf.gz
        tabix baseline_fixed.vcf.gz
        rm -f reheaded.vcf.gz*

        echo "#### Fix GT for svimmer"
        python3 ~{docker_dir}/add_gt.py ${all_vcfs["svimmer"]} ~{sample_id} | bcftools sort -O z -o fixed_svimmer.vcf.gz
        tabix -f fixed_svimmer.vcf.gz
        all_vcfs["svimmer"]="fixed_svimmer.vcf.gz"

        echo "#### Resolving variants"
        comparison_vcfs=("jasmine" "svimmer" "svmerger")
        for key in ${comparison_vcfs[@]};
        do
            val=${all_vcfs[$key]}
            echo "Key: $key, Value: $val"
            outname="resolved_$(basename $val)"
            python ~{docker_dir}/resolve.py $val ~{reference} \
                | bcftools view -i "SVTYPE != 'BND'" \
                | bcftools sort -O z -o $outname
            tabix -f $outname
            all_vcfs[$key]="$outname"
        done

        echo "#### Run Truvari bench"
        mkdir -p bench_results/
        comparison_vcfs+=("truvari" "bcftools")
        for key in "${comparison_vcfs[@]}"; do
            val=${all_vcfs[$key]}
            echo "Key: $key, Value: $val"

            truvari bench -b baseline_fixed.vcf.gz \
               -c $val \
               --includebed ~{truth_bed} \
               --pick multi \
               -o ./bench_results/${key}
        done
        tar czf ~{sample_id}_bench_results.tar.gz ./bench_results/

        echo "#### Run kdp on each of the vcfs"
        for key in "${comparison_vcfs[@]}"; do
            val=${all_vcfs[$key]}
            echo "Key: $key, Value: $val"
            outname="kdp_$(basename ${val})"
            python ~{docker_dir}/kfdphase_cosine.py --sizemin 50 -b baseline_fixed.vcf.gz --kmer 4 --pg -c ${val} -o ${outname}
            all_vcfs[$key]="$outname"
        done

        echo "#### Run merge region analysis"
        python ~{docker_dir}/merge_report_table.py \
            --truth baseline_fixed.vcf.gz \
            --bcftools ${all_vcfs["bcftools"]} \
            --regions  ~{truth_bed} \
            --collapsed ${all_vcfs["truvari"]} \
            --jasmine ${all_vcfs["jasmine"]} \
            --svimmer ${all_vcfs["svimmer"]} \
            --svmerger ${all_vcfs["svmerger"]} \
            --output region.table.tsv

        echo "### Making region summary"
        python ~{docker_dir}/summarize_merge_report_table.py region.table.tsv > region.summary.txt
    >>>
    
    output {
        File bench_results = work_dir + "/" + sample_id + "_bench_results.tar.gz"
        File region_table = work_dir + "/region.table.tsv"
        File region_summary = work_dir + "region.summary.txt"
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
        cpu: 1
        memory: "16GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
