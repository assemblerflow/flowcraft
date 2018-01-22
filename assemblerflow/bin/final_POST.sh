#!/usr/bin/env sh

st=$(cat $(pwd)/.status)

json="{'project_id':'$1','pipeline_id':'$2','process_id':'$3','run_info':'None','run_output':'None','warnings':'$(pwd)/.warning','log_file':'$(pwd)/.command.log','status':'$st','type':'output'}"

{
    curl -H  "Content-Type: application/json" -L -X POST -d \"$json\" $4 > /dev/null
} || {
    echo Curl request failed
}