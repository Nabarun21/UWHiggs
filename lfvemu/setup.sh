#!/bin/bash
export OVERRIDE_META_TREE_data_EM='em/metaInfo'

export IGNORE_LUMI_ERRORS=1

source jobid.sh
export jobid=$jobid13

echo $jobid
#export datasrc=/hdfs/store/user/$USER/  #$(ls -d /scratch/*/data/$jobid | awk -F$jobid '{print $1}')
#export MEGAPATH=/hdfs/store/user/$USER
export datasrc=/hdfs/store/user/taroni  #$(ls -d /scratch/*/data/$jobid | awk -F$jobid '{print $1}')
export MEGAPATH=/hdfs/store/user/taroni

./make_proxies.sh
rake "meta:getinputs[$jobid, $datasrc,em/metaInfo]"
rake "meta:getmeta[inputs/$jobid, em/metaInfo, 13]"
#rake "meta:getinputs[$jobid, $datasrc,ee/metaInfo]"
#rake "meta:getmeta[inputs/$jobid, ee/metaInfo, 13]"

unset OVERRIDE_META_TREE_data_EM
unset IGNORE_LUMI_ERRORS
