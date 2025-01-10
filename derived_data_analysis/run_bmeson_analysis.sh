#!/bin/bash

MC=0                    # Default value for MC
PARTICLE="b0"           # Default value for PARTICLE
CHARMDAUGHTER="d"
JOBS=1
DELETE=0

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
        -j)
            if [[ -z "$2" ]]; then
                echo "Error: Missing argument for -j"
                exit 1
            fi
            JOBS="$2"
            shift # Move past argument value
            shift # Move past argument name
            ;;
        -d)
            DELETE=1
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
FILENAME="input_mc.txt"
if [ $MC -eq 0 ]; then
    JSONCONFIG=dpl-config_${PARTICLE}.json
    FILENAME="input_data.txt"
fi
COMMONCONFIGS="-b --aod-memory-rate-limit 1000000000 --configuration json://../$JSONCONFIG --shm-segment-size 7500000000"

# a complex way to assign a suffix to our outputs
AODNAME=$(head -n 1 $FILENAME)
IFS='/' read -ra AODNAMENOPATH <<< $AODNAME
LASTIDX=$(( ${#AODNAMENOPATH[@]} - 1 ))
AODNAMENOPATH=${AODNAMENOPATH[$LASTIDX]}
IFS='_' read -ra SUFFIX <<< $AODNAMENOPATH
SUFFIX=${SUFFIX[0]}
SUFFIX=$(echo "$AODNAMENOPATH" | sed "s/$SUFFIX//g")

# Add new line to the end of the file, if there isn't one
sed -i -e '$a\' "$FILENAME"

# Count the number of AODs in the file
NAO2D=$(grep -cve '^$' "$FILENAME")
NFILESPERJOB=$(($NAO2D / $JOBS))
if (( NAO2D % JOBS != 0 )); then
    (( NFILESPERJOB++ )) # Increment num by 1
fi

AODCOUNTER=0
NFILES=0
declare -a INPUTNAMES=()
while IFS= read -r line; do

    # Add the file name to the list of inputs
    if [[ $AODCOUNTER -eq 0 ]]; then
        INPUTNAMES+=("inputs_file_temp_$NFILES.txt")
    fi
    # Save the file name into a grouped file
    echo "$line" >> "inputs_file_temp_$NFILES.txt"
    AODCOUNTER=$((AODCOUNTER + 1))
    if [[ $AODCOUNTER -eq $NFILESPERJOB ]]; then
        AODCOUNTER=0
        NFILES=$((NFILES + 1))
    fi
done < "$FILENAME"

parallel -j "$JOBS" 'mkdir -p output_{#} && cd output_{#} && o2-analysis-hf-task-b0-reduced '"$COMMONCONFIGS"' | o2-analysis-hf-candidate-selector-'"$PARTICLE"'-to-d-pi-reduced '"$COMMONCONFIGS"' | o2-analysis-hf-candidate-creator-'"$PARTICLE"'-reduced '"$COMMONCONFIGS"' --aod-writer-json ../OutputDirector_'"$PARTICLE"'.json --resources-monitoring 2 --fairmq-ipc-prefix . --aod-file @../{1} > log_{#}.txt 2>&1 && cd ..' ::: ${INPUTNAMES[@]}

# loop on all the output files and merge them
declare -a ANALYSISRESULTS=()
for i in $(seq 1 $((NFILES))); do
    if [ -f "output_$i/AO2D.root" ]; then
        echo "output_$i/AO2D.root" >> ao2ds_to_merge.txt
        ANALYSISRESULTS+=("output_$i/AnalysisResults.root")
    fi
done

o2-aod-merger --input ao2ds_to_merge.txt --max-size 10000000000000 --output Tree$SUFFIX
hadd -f AnalysisResultsTask$SUFFIX ${ANALYSISRESULTS[@]}
rm inputs_file_temp_*.txt
rm ao2ds_to_merge.txt
rm ./**/localhost*

if [ $DELETE -eq 1 ]; then
    rm -rf output_*
fi