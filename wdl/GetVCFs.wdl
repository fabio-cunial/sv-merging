version 1.0


#
workflow GetVCFs {
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
        
        function uploadVCF() {
            local FILE1=$1
            local FILE2=$2
            
            TEST=$(gsutil -q stat ~{output_dir}/${FILE1} && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                while : ; do
                    TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ${FILE1} ${FILE2} ~{output_dir} && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error uploading file <${FILE1}>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
            fi
        }
        
        # Downloading VCF regions
        TEST=$(cat ~{vcf_addresses} | gsutil -m cp -I . && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            echo "Error downloading VCF files. Trying again..."
            sleep ${GSUTIL_DELAY_S}
        fi
        touch list.txt
        while read VCF_FILE; do
            LOCAL_FILE=$(basename -s .vcf.gz ${VCF_FILE})
            tabix ${LOCAL_FILE}.vcf.gz
            bcftools filter --threads ${N_THREADS} --regions "~{region}" --include "FILTER=\"PASS\"" --output-type z ${LOCAL_FILE}.vcf.gz > tmp.vcf.gz
            rm -f ${LOCAL_FILE}.vcf.gz ${LOCAL_FILE}.vcf.gz.tbi
            tabix tmp.vcf.gz
            bcftools sort --output-type z tmp.vcf.gz > ${LOCAL_FILE}.region.vcf.gz
            rm -f tmp.vcf.gz tmp.vcf.gz.tbi
            tabix ${LOCAL_FILE}.region.vcf.gz
            echo ${LOCAL_FILE}.region.vcf.gz >> list.txt
        done < ~{vcf_addresses}
        
        # BCFTOOLS MERGE
        #${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --apply-filters PASS --merge none --file-list list.txt --output-type z --output merged.1.vcf.gz
        #uploadVCF "merged.1.vcf.gz" " "
        #tabix merged.1.vcf.gz
        
        # TRUVARI COLLAPSE
        TEST=$(gsutil -q stat ~{output_dir}/merged.2.vcf && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            ${TIME_COMMAND} truvari collapse --input merged.1.vcf.gz --output merged.2.vcf --collapsed-output collapsed.2.vcf --reference ~{reference_fa} --keep first --passonly
            uploadVCF merged.2.vcf collapsed.2.vcf
            rm -f merged.2.vcf collapsed.2.vcf
        fi
        
        TEST=$(gsutil -q stat ~{output_dir}/merged.3.vcf && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            ${TIME_COMMAND} truvari collapse --input merged.1.vcf.gz --output merged.3.vcf --collapsed-output collapsed.3.vcf --reference ~{reference_fa} --keep maxqual --passonly
            uploadVCF merged.3.vcf collapsed.3.vcf
            rm -f merged.3.vcf collapsed.3.vcf
        fi
        
        TEST=$(gsutil -q stat ~{output_dir}/merged.4.vcf && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            ${TIME_COMMAND} truvari collapse --input merged.1.vcf.gz --output merged.4.vcf --collapsed-output collapsed.4.vcf --reference ~{reference_fa} --keep common --passonly
            uploadVCF merged.4.vcf collapsed.4.vcf
            rm -f merged.4.vcf collapsed.4.vcf
        fi
        
        TEST=$(gsutil -q stat ~{output_dir}/merged.5.vcf && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            ${TIME_COMMAND} truvari collapse --input merged.1.vcf.gz --output merged.5.vcf --collapsed-output collapsed.5.vcf --reference ~{reference_fa} --keep common --chain --passonly
            uploadVCF merged.5.vcf collapsed.5.vcf
            rm -f merged.5.vcf collapsed.5.vcf
        fi
        
        # SVPOP. Parameters are set to approximate truvari's defaults.
        #TEST=$(gsutil -q stat ~{output_dir}/svpop.vcf.gz && echo 0 || echo 1)
        #TEST=0
        #if [ ${TEST} -eq 1 ]; then
        #    echo -e "chr22\t50818468\t2864531079\t70\t71" > chr22only.fai
        #    rm -rf config/; mkdir config/
        #    touch config/samples.tsv
        #    echo -e "NAME\tSAMPLE\tTYPE\tDATA\tVERSION\tPARAMS\tCOMMENT" >> config/samples.tsv
        #    echo -e "sniffles2\tDEFAULT\tsniffles\t~{work_dir}/{sample}.vcf.gz\t2\t\t"  >> config/samples.tsv
        #    cat config/samples.tsv
        #    touch config/config.json
        #    echo "{" >> config/config.json
        #    echo "\"reference\": \"~{reference_fa}\"," >> config/config.json
        #    echo "\"reference_fai\": \"chr22only.fai\"," >> config/config.json
        #    echo "\"ucsc_ref_name\": \"hg38\"," >> config/config.json
        #    # samplelist section
        #    echo "\"samplelist\": {" >> config/config.json
        #    echo "\"mySamples\": [" >> config/config.json
        #    i=0
        #    while read VCF_FILE; do
        #        if [ $i -eq 0 ]; then 
        #            echo -n "\"${VCF_FILE%.vcf.gz}\"" >> config/config.json
        #            i=1
        #        else
        #            echo -en ",\n\"${VCF_FILE%.vcf.gz}\"" >> config/config.json
        #        fi
        #    done < list.txt
        #    echo -e "\n]" >> config/config.json
        #    echo "}," >> config/config.json
        #    # sampleset section
        #    echo "\"sampleset\": {" >> config/config.json
        #    echo "\"myMerge\": {" >> config/config.json
        #    echo "\"sourcetype\": \"caller\"," >> config/config.json
        #    echo "\"sourcename\": \"sniffles2\"," >> config/config.json
        #    echo "\"merge\": \"nr::szro(szro=0.7,dist=500,match(score=0.7,limit=4000,ksize=9))\"," >> config/config.json
        #    echo "\"name\": \"myMerge\"," >> config/config.json
        #    echo "\"description\": \"myMerge\"" >> config/config.json
        #    echo "}," >> config/config.json
        #    echo "}" >> config/config.json
        #    # end of config file
        #    echo "}" >> config/config.json
        #    cat config/config.json
        #    source activate svpop
        #    # INS return error at line: <https://github.com/EichlerLab/svpop/blob/1d62b72187172a7898443c5b1a27f3838621a199/svpoplib/svmerge.py#L822>
        #    ${TIME_COMMAND} snakemake -s ~{docker_dir}/svpop/Snakefile --cores ${N_THREADS} results/variant/sampleset/myMerge/mySamples/all/all/bed/sv_ins.bed.gz
        #    # DEL returns error at line: <https://github.com/EichlerLab/svpop/blob/1d62b72187172a7898443c5b1a27f3838621a199/svpoplib/svmerge.py#L822>
        #    ${TIME_COMMAND} snakemake -s ~{docker_dir}/svpop/Snakefile --cores ${N_THREADS} results/variant/sampleset/myMerge/mySamples/all/all/bed/sv_del.bed.gz
        #    ${TIME_COMMAND} snakemake -s ~{docker_dir}/svpop/Snakefile --cores ${N_THREADS} results/variant/sampleset/myMerge/mySamples/all/all/bed/sv_inv.bed.gz
        #    ${TIME_COMMAND} snakemake -s ~{docker_dir}/svpop/Snakefile --cores ${N_THREADS} results/variant/sampleset/myMerge/mySamples/all/all/bed/sv_dup.bed.gz
        #    conda deactivate
        #    zcat results/variant/sampleset/myMerge/mySamples/all/all/bed/sv_ins.bed.gz \
        #         results/variant/sampleset/myMerge/mySamples/all/all/bed/sv_del.bed.gz \
        #         results/variant/sampleset/myMerge/mySamples/all/all/bed/sv_inv.bed.gz \
        #         results/variant/sampleset/myMerge/mySamples/all/all/bed/sv_dup.bed.gz | gzip > svpop.bed.gz
        #    uploadVCF svpop.bed.gz " "
        #fi
        
        
        
        
        
        
        
        
        
        # ----------- Decompressing VCFs for the following tools ---------------
        touch list_decompressed.txt
        while read VCF_GZ_FILE; do
            gunzip ${VCF_GZ_FILE}
            echo ${VCF_GZ_FILE%.gz} >> list_decompressed.txt
        done < list.txt
        
        # SURVIVOR. Parameters are set to truvari's defaults. survivor does not
        # work on .vcf.gz files.
        TEST=$(gsutil -q stat ~{output_dir}/survivor.vcf && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            ${TIME_COMMAND} SURVIVOR merge list_decompressed.txt 500 1 1 1  0  50 survivor.vcf
            uploadVCF survivor.vcf " "
            rm -f survivor.vcf
        fi
        
        # JASMINE. Parameters are set to truvari's defaults, except:
        # --clique_merging: removed because too slow.
        # min_seq_id=0.7 --use_edit_dist: removed because too slow.
        TEST=$(gsutil -q stat ~{output_dir}/jasmine.vcf && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            source activate jasmine
            ${TIME_COMMAND} jasmine threads=${N_THREADS} max_dist=500 file_list=list_decompressed.txt out_file=jasmine.vcf
            conda deactivate
            uploadVCF jasmine.vcf " "
            rm -f jasmine.vcf
        fi
        
        # SVIMMER. Parameters are set to approximate truvari's defaults, except
        # that SVIMMER uses an absolute size difference instead of a relative
        # one.
        # Remark: according to its paper, "svimmer ignores the samples' genotype
        # information to reduce compute time and memory, as only SV site
        # information is needed for GraphTyper’s graph construction." Otherwise
        # they say it's similar to SURVIVOR.
        TEST=$(gsutil -q stat ~{output_dir}/svimmer.vcf && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            ${TIME_COMMAND} svimmer --threads ${N_THREADS} --ids --max_distance 500 --output svimmer.vcf list_decompressed.txt chr22
            uploadVCF svimmer.vcf " "
            rm -f svimmer.vcf
        fi
        
        # SVTOOLS
        TEST=$(gsutil -q stat ~{output_dir}/svtools.vcf && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            ${TIME_COMMAND} svtools lsort list_decompressed.txt svtools.1.vcf
            ${TIME_COMMAND} svtools lmerge svtools.1.vcf > svtools.vcf
            rm -f svtools.1.vcf
            uploadVCF svtools.vcf " "
            rm -f svtools.vcf
        fi
        
    >>>
    output {
    }
    runtime {
        docker: "fcunial/sv-merging"
        cpu: 1  # bcftools and truvari are sequential
        memory: "64GB"  # Arbitrary
        disks: "local-disk 250 HDD"  # Arbitrary
        preemptible: 0
    }
}
