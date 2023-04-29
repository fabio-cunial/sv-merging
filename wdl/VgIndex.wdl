version 1.0


#
workflow VgIndex {
    input {
        String caller
        Array[Int] coverages
        Int max_sv_length
        String remote_dir
        File reference_fa
        File reference_fai
    }
    parameter_meta {
    }
    scatter(c in coverages) {
        call VgIndexImpl {
            input:
                caller = caller,
                coverage = c,
                max_sv_length = max_sv_length,
                remote_dir = remote_dir,
                reference_fa = reference_fa,
                reference_fai = reference_fai
        }
    }
    output {
    }
}


task VgIndexImpl {
    input {
        String caller
        Int coverage
        Int max_sv_length
        String remote_dir
        File reference_fa
        File reference_fai
    }
    parameter_meta {
        remote_dir: "Containing the VCF to be indexed."
        max_sv_length: "Used for cleaning."
    }
    
    String docker_dir = "/infogain"
    String work_dir = "/cromwell_root/infogain"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        VG_COMMAND="~{docker_dir}/vg"
        
        cat ~{reference_fa}| awk '{ if (substr($0,1,1)==">") { filename=(substr($1,2) ".fa") } print $0 >> filename; close(filename) }'
        gsutil -m cp ~{remote_dir}/~{caller}-~{coverage}.vcf.gz .
        gunzip ~{caller}-~{coverage}.vcf.gz
        java -cp ~{docker_dir} CleanVCF ~{caller}-~{coverage}.vcf . ~{max_sv_length} 1 ~{caller}-~{coverage}-vg.vcf
        ${TIME_COMMAND} ${VG_COMMAND} construct --threads ${N_THREADS} --handle-sv --alt-paths --vcf ~{caller}-~{coverage}-vg.vcf --reference ~{reference_fa} > ~{caller}-~{coverage}.vg
        mkdir ./vgtmp
        ${TIME_COMMAND} ${VG_COMMAND} index --threads ${N_THREADS} --temp-dir ./vgtmp --progress --xg-alts --xg-name ~{caller}-~{coverage}.xg --gcsa-out ~{caller}-~{coverage}.gcsa ~{caller}-~{coverage}.vg
        rm -rf ./vgtmp
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp ~{caller}-~{coverage}-vg.vcf ~{caller}-~{coverage}.vg ~{caller}-~{coverage}.xg ~{caller}-~{coverage}.gcsa ~{remote_dir} && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading indexes. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/infogain"
        cpu: 32
        memory: "128GB"
        disks: "local-disk 512 HDD"
        preemptible: 0
    }
}
