#!/bin/bash

# tell me if I am analysing data or MC
MC=0

JSONCONFIG=dpl-config_bmes_mc.json
if [ $MC -eq 0 ]; then
    JSONCONFIG=dpl-config_bmes.json
fi
COMMONCONFIGS="-b --aod-memory-rate-limit 1000000000 --configuration json://$JSONCONFIG --shm-segment-size 7500000000"

# a complex way to assign a suffix to our outputs
JSONDATA=$(cat $JSONCONFIG)
FILENAMEWITHAT=$(echo "$JSONDATA" | jq -r '."internal-dpl-aod-reader"."aod-file-private"')
FILENAME=$(echo "$FILENAMEWITHAT" | sed 's/@//g')
AODNAME=$(head -n 1 $FILENAME)
IFS='/' read -ra AODNAMENOPATH <<< $AODNAME
LASTIDX=$(( ${#AODNAMENOPATH[@]} - 1 ))
AODNAMENOPATH=${AODNAMENOPATH[$LASTIDX]}
IFS='_' read -ra SUFFIX <<< $AODNAMENOPATH
SUFFIX=${SUFFIX[0]}
SUFFIX=$(echo "$AODNAMENOPATH" | sed "s/$SUFFIX//g")

o2-analysis-hf-task-b0-reduced $COMMONCONFIGS |
o2-analysis-hf-candidate-selector-b0-to-d-pi-reduced $COMMONCONFIGS |
o2-analysis-hf-candidate-creator-b0-reduced $COMMONCONFIGS --aod-writer-json OutputDirector.json --resources-monitoring 2 --fairmq-ipc-prefix .

mv AnalysisResults.root AnalysisResultsTask$SUFFIX
touch ao2ds_to_merge.txt
echo AO2D.root > ao2ds_to_merge.txt
o2-aod-merger --input ao2ds_to_merge.txt --max-size 10000000000000 --output Tree$SUFFIX
rm AO2D.root
rm ao2ds_to_merge.txt
rm localhost*