#!/bin/bash
# Run all of the analysis

set -o nounset
set -o errexit

export MEGAPATH=/hdfs/store/user/fmeng/
source jobid.sh
#jobid13=LFV_8013v1
export jobid=$jobid13

#rake genkin
#rake recoplots
#rake recoplotsMVA
#rake controlplots
#rake controlplotsMVA
#rake fakeeet
#rake fakeeetMVA
#rake  efits
#rake drawTauFakeRate
#export jobid=$jobidmt
rake recoplotsMVA_C
#rake stitched
#rake reco2
#rake emuOptim
##rake drawplots
#rake genkinEMu
#rake genkinMuTau
#rake fakemmmMVA
#rake fakeeemMVA
#rake fakeeeeMVA
#rake fakemmeMVA
#rake  fits
#rake efits
#rake zmm
#rake test
