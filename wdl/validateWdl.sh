#!/bin/bash
#
set -x
WOMTOOL_PATH="/Users/fcunial/apps/cromwell/womtool-84.jar"

java -jar ${WOMTOOL_PATH} validate -l GetVCFs.wdl -i inputs/private_GetVCFs.json
java -jar ${WOMTOOL_PATH} validate -l VisualizeSVs.wdl -i inputs/private_VisualizeSVs.json
java -jar ${WOMTOOL_PATH} validate -l JointCalling.wdl -i inputs/private_JointCalling.json
java -jar ${WOMTOOL_PATH} validate -l OverlapStats.wdl -i inputs/private_OverlapStats.json