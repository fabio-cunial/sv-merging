#!/bin/bash
#
DOCKER_DIR=$1
WORK_DIR=$2
VCF_GZ=$3
REFERENCE_FA=$4

GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
GSUTIL_DELAY_S="600"
N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
TIME_COMMAND="/usr/bin/time --verbose"

gunzip -c ${VCF_GZ} > input.vcf
cp -r ${DOCKER_DIR}/vcf-merging/pangenome-graph-from-callset/ ./
cd pangenome-graph-from-callset/
rm -f config.yaml
echo "vcf: ${WORK_DIR}/input.vcf" >> config.yaml
echo "reference: ${REFERENCE_FA}" >> config.yaml
echo "outdir: ${WORK_DIR}/results" >> config.yaml
source activate svpop
${TIME_COMMAND} snakemake -j ${N_THREADS}
conda deactivate
cd ${WORK_DIR}
bgzip -@ ${N_THREADS} ./results/pangenome/pangenome.vcf
tabix -f ./results/pangenome/pangenome.vcf.gz
