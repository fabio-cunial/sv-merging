version 1.0


workflow JointCalling {
    input {
        File snfs
    }

    call JointCallingImpl { 
        input: 
            snfs = snfs
    }
    
    output {
        File joint_vcf = JointCallingImpl.joint_vcf
        File joint_tbi = JointCallingImpl.joint_tbi
        File counts = JointCallingImpl.counts
    }
}


task JointCallingImpl {
    input {
        File snfs
    }

    command <<<
        set -euxo pipefail
        
        GSUTIL_DELAY_S="600"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        
        rm -f list.txt; touch list.txt
        rm -f counts.txt; touch counts.txt
        while read SNF_FILE; do
            LOCAL_FILE=$(basename ${SNF_FILE})
            while : ; do
                TEST=$(gsutil -m cp ${SNF_FILE} ${LOCAL_FILE} && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <${SNF_FILE}>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            echo "${LOCAL_FILE}\t${LOCAL_FILE%.snf}" >> list.tsv
        done < ~{snfs}
        head list.tsv
        ${TIME_COMMAND} sniffles --threads ${N_THREADS} --combine-separate-intra True --input list.tsv --vcf joint.vcf
        N_INS=$(grep "SVTYPE=INS" joint.vcf | awk '{ if ($7=="PASS") print $0; }' | wc -l)
        N_DEL=$(grep "SVTYPE=DEL" joint.vcf | awk '{ if ($7=="PASS") print $0; }' | wc -l)
        echo "${N_INS},${N_DEL}" >> counts.txt
        bgzip --threads ${N_THREADS} joint.vcf
        tabix joint.vcf.gz
    >>>

    output {
        File joint_vcf = "joint.vcf.gz"
        File joint_tbi = "joint.vcf.gz.tbi"
        File counts = "counts.txt"
    }
    runtime {
        cpu: 4
        memory: "16 GiB"
        disks: "local-disk 100 HDD"
        preemptible: 0
        docker: "fcunial/simulation"
    }
}