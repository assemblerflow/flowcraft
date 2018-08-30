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
overwrite=${10}

json_str="{}"
version_str="[]"
trace_str=""

# If a .report.json file was populated, set the json_str variable
if [ -s .report.json ];
then

    # Modification of the JSON string should be different for chewbbaca
    # output
    if [ $task = "chewbbaca" ];
    then
        json_str=$(cat $(pwd)/.report.json | sed 's/ //g' | sed s/\"/\'/g)
    else
        json_str=$(cat $(pwd)/.report.json | sed 's/ /%20/g' | sed s/\"/\'/g)
    fi
fi

# If a .versions file was populated, set the version_str variable
if [ -s .versions ];
then
    version_str=$(< $(pwd)/.versions sed 's/ /%20/g' | sed s/\"/\'/g)
fi

if [ -s .command.trace ];
then
    trace_str="$(< $(pwd)/.command.trace tr "\n" ";")"
fi

# If a .versions OR .report.json file was populated send the request
if [ ! "$json_str" = "{}" ] || [ ! "$version_str" = "[]" ] || [ ! "$trace_str" = "" ];
then
    workdir=$(pwd)
    json="{'projectid':'$projectid','pipelineId':'$pipelineid','processId':'$processid','sample_name':'$sample','reportJson':$json_str,'username':'$username','userId':'$userid','workdir':'$workdir','task':'$task','processName':'$task','species':'$species','versions':$version_str,'trace':'$trace_str', 'overwrite': '$overwrite'}"
    echo \"${json}\" > .final.json
    {
        cat .final.json | curl -H  "Content-Type: application/json" -k -L -X POST -d @- $url > /dev/null
    } || {
        echo Curl request failed
    }

fi
