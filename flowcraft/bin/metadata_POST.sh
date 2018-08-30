#!/usr/bin/env sh

set -ex

projectid=$1
pipelineid=$2
processid=$3
sample=$4
url=$5
username=$6
userid=$7
task=$8
species=$9

metadata_str="{}"

# If a .report.json file was populated, set the json_str variable
if [ -s .metadata.json ];
then
    metadata_str=$(cat $(pwd)/.metadata.json | sed 's/ /%20/g' | sed s/\"/\'/g)
fi

# If a .versions OR .report.json file was populated send the request
if [ ! "$metadata_str" = "{}" ];
then
    workdir=$(pwd)
    json="{'projectid':'$projectid','pipelineId':'$pipelineid','processId':'nfMetadata','sample_name':'$sample','nfMetadata':$metadata_str,'username':'$username','userId':'$userid','workdir':'$workdir','task':'nfMetadata','processName':'nfMetadata','species':'$species','overwrite':'false'}"
    echo \"${json}\" > .final.json
    {
        cat .final.json | curl -H  "Content-Type: application/json" -k -L -X POST -d @- $url > /dev/null
    } || {
        echo Curl request failed
    }

fi
