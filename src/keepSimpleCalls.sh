#!/bin/bash
#
CALLER="jasmine"
PASS_ONLY="1"
TRF_BED="human_GRCh38_no_alt_analysis_set.trf.bed"
IDENTITY_THRESHOLD="100"
SV_LENGTH_THRESHOLD="10000"


for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
    bcftools view ${CALLER}_chr${CHR}.vcf.gz > ${CALLER}_chr${CHR}.vcf
    java -Xmx12G KeepSimpleCalls ${CALLER}_chr${CHR}.vcf ${PASS_ONLY} ${TRF_BED} ${IDENTITY_THRESHOLD} ${SV_LENGTH_THRESHOLD} ${CALLER}_chr${CHR}_short.vcf ${CALLER}_chr${CHR}_long.vcf
    rm -f ${CALLER}_chr${CHR}.vcf
done
