#!/bin/bash

MC=0                    # Default value for MC
PARTICLE="b0"           # Default value for PARTICLE
CHARMDAUGHTER="d"

while [[ $# -gt 0 ]]; do
    case $1 in
        --mc)
            if [[ -z "$2" ]]; then
                echo "Error: Missing argument for --mc"
                exit 1
            fi
            MC="$2"
            shift # Move past argument value
            shift # Move past argument name
            ;;
        --particle)
            if [[ -z "$2" ]]; then
                echo "Error: Missing argument for --particle"
                exit 1
            fi
            PARTICLE="$2"
            if [[ "$PARTICLE" != "bplus" && "$PARTICLE" != "b0" ]]; then
                echo "Wrong particle: only 'bplus' and 'b0' are accepted"
                exit 1
            fi
            shift # Move past argument value
            shift # Move past argument name
            ;;
        *)    # Unknown option
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Output the final values
echo "MC: $MC"
echo "Particle: $PARTICLE"

if [ "$PARTICLE" == "bplus" ]; then
    CHARMDAUGHTER="d0"
fi

JSONCONFIG=dpl-config_${PARTICLE}_mc.json
if [ $MC -eq 0 ]; then
    JSONCONFIG=dpl-config_${PARTICLE}.json
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

o2-analysis-hf-task-${PARTICLE}-reduced $COMMONCONFIGS |
o2-analysis-hf-candidate-selector-${PARTICLE}-to-${CHARMDAUGHTER}-pi-reduced $COMMONCONFIGS |
o2-analysis-hf-candidate-creator-${PARTICLE}-reduced $COMMONCONFIGS --aod-writer-json OutputDirector_${PARTICLE}.json --resources-monitoring 2 --fairmq-ipc-prefix .

mv AnalysisResults.root AnalysisResultsTask$SUFFIX
touch ao2ds_to_merge.txt
echo AO2D.root > ao2ds_to_merge.txt
o2-aod-merger --input ao2ds_to_merge.txt --max-size 10000000000000 --output Tree$SUFFIX
rm AO2D.root
rm ao2ds_to_merge.txt
rm localhost*