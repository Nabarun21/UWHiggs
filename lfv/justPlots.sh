#!/bin/bash                                                                                                                                   

usage='Usage: -a <analyzer name>  -lumi <luminosity in pb> -ns <num_samples> -nf <sampling frequency> -ph <phase position> (-cpt <old pulse type> -cns\
 <old no. of samples> -cnf<old sampling freq>) '

args=`getopt rdlp: -- "$@"`
if test $? != 0
     then
         echo $usage
         exit 1
fi

eval set -- "$args"


for i
 do
    case "$i" in
      -analyzer) shift; analyzer=$2;shift;;
      -lumi) shift; luminosity=$2;shift;;
      -jobid) shift;jobid=$2;shift;;
      -analtype) shift;analtype=$2;shift;;
      -kinplots) shift;kinplots=$2;shift;;
    esac
done

echo $analyzer
echo $luminosity

echo $jobid
echo $analtype


source jobid.sh
#plot opposite sign region Using fakes from MC
python plotterMCBDT.py $analyzer $luminosity $jobid $analtype $kinplots "True"

#plot same sign Using fakes from MC
python plotterMCssBDT.py $analyzer $luminosity $jobid $analtype $kinplots

#plot same sign Using fakes from data
python plotterDataFakesBDT.py $analyzer $luminosity $jobid $analtype $kinplots

#plot same sign Using fakes from data
python plotterDataFakesMethod2BDT.py $analyzer $luminosity $jobid $analtype $kinplots

#plot anti-isolated OS region
python plotterAntiIsolatedOSBDT.py $analyzer $luminosity $jobid $analtype $kinplots

#plot anti-isolated SS region
python plotterAntiIsolatedSSBDT.py $analyzer $luminosity $jobid $analtype $kinplots

#plot weighted anti-isolated SS region all fakes -- muon + electron -electron muon
python plotterAntiIsolatedWeightedSSBDT.py  $analyzer $luminosity $jobid $analtype $kinplots 

#plot weighted anti-isolated OS region
python plotterAntiIsolatedWeightedOSBDT.py  $analyzer $luminosity $jobid $analtype $kinplots

#plot weighted anti-isolated OS region electron fakes
python plotterAntiIsolatedWeightedElectronOSBDT.py $analyzer $luminosity $jobid $analtype $kinplots


#plot weighted anti-isolated SS region electron fakes
python plotterAntiIsolatedWeightedElectronSSBDT.py $analyzer $luminosity $jobid $analtype $kinplots


#plot weighted anti-isolated OS region muon fakes                                                                                                                               
python plotterAntiIsolatedWeightedMuonOSBDT.py $analyzer $luminosity $jobid $analtype $kinplots


#plot weighted anti-isolated SS region muon fakes                                                                                                                               
python plotterAntiIsolatedWeightedMuonSSBDT.py $analyzer $luminosity $jobid $analtype $kinplots


#plot weighted anti-isolated OS region muon electron fakes                                                                                                                               
python plotterAntiIsolatedWeightedMuonElectronOS.py $analyzer $luminosity $jobid $analtype $kinplots


#plot weighted anti-isolated SS region muon electron fakes                                                                                                                               
python plotterAntiIsolatedWeightedMuonElectronSS.py $analyzer $luminosity $jobid $analtype $kinplots
