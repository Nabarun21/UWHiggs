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
      -sys) shift;sys=$2;shift;;
    esac
done

echo $analyzer
echo $luminosity

echo $jobid
echo $analtype


rm -r LFVHEMuAnalyzerMVA$analyzer$luminosity*


#copy results into current directory in the form wanted
cp -r results/$jobid/LFVHEMuAnalyzerMVA$analyzer LFVHEMuAnalyzerMVA$analyzer$luminosity
cp move.sh LFVHEMuAnalyzerMVA$analyzer$luminosity 
cd LFVHEMuAnalyzerMVA$analyzer$luminosity
source move.sh
rm move.sh
cd -


#get QCD (data-MC) in ss *2.30(SF)
python getQCD.py $analyzer $luminosity $jobid $analtype $sys 
python getQCDforcombine.py $analyzer $luminosity $jobid $analtype $sys   #get the same without sysuncertainty in histos, combine script already adds them


#construct  Fakes histos from anti-isolated region
#python getFakesBDT.py $analyzer $luminosity $jobid $analtype $sys
#python getFakesforcombineBDT.py $analyzer $luminosity   $jobid $analtype $sys  #get the same without sysuncertainty in histos, combine script already adds them

#construct  Fakes histos from anti-isolated region:method 2 
#python getFakesmethod2BDT.py $analyzer $luminosity $jobid $analtype $sys
#python getFakesmethod2forcombineBDT.py $analyzer $luminosity $jobid $analtype $sys   #get the same without sysuncertainty in histos, combine script already adds them


#folder for plotting
cp -r LFVHEMuAnalyzerMVA$analyzer$luminosity LFVHEMuAnalyzerMVA$analyzer$luminosity'plot'
#cp FAKES$analyzer.root LFVHEMuAnalyzerMVA$analyzer$luminosity'plot'/FAKES.root
#cp FAKESmethod2$analyzer.root LFVHEMuAnalyzerMVA$analyzer$luminosity'plot'/FAKESmethod2.root
cp QCD$analyzer.root LFVHEMuAnalyzerMVA$analyzer$luminosity'plot'/QCD.root

#folders for constructing limits
#cp -r LFVHEMuAnalyzerMVA$analyzer$luminosity LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromdata'
#cp -r LFVHEMuAnalyzerMVA$analyzer$luminosity LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromdatamethod2'
cp -r LFVHEMuAnalyzerMVA$analyzer$luminosity LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromMC'
#cp FAKESforcombine$analyzer.root LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromdata'/FAKES.root
#cp FAKESmethod2forcombine$analyzer.root LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromdatamethod2'/FAKES.root
#rm LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromdatamethod2'/WJETSMC.root
#rm LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromdata'/QCD.root
cp QCDforcombine$analyzer.root LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromMC'/QCD.root
#rm LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromMC'/*ETau*
#rm LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromdata'/*ETau*

echo 'hoiyagese'
echo 'hoiyagese'
echo 'hoiyagese'
echo 'hoiyagese'
echo 'hoiyagese'
echo 'hoiyagese'
echo 'hoiyagese'
echo 'hoiyagese'
echo 'hoiyagese'
echo 'hoiyagese'
echo 'hoiyagese'
echo 'hoiyagese'
echo 'hoiyagese'
echo 'hoiyagese'
echo 'hoiyagese'
echo 'hoiyagese'
echo 'hoiyagese'
echo 'hoiyagese'

#in plotting folder remove neg bins + add systematic uncertainities
python prepHistos.py $analyzer $luminosity $jobid $analtype

#plotting replace weight files by current lumi (I am a superfluos remnant feature from ice ages, get rid of me)
python replaceFileWeightBy1.py $jobid $luminosity


source jobid.sh
#plot opposite sign region Using fakes from MC
python plotterMCBDT.py $analyzer $luminosity $jobid $analtype $kinplots "True"

##plot same sign Using fakes from MC
#python plotterMCssBDT.py $analyzer $luminosity $jobid $analtype $kinplots 
#
##plot same sign Using fakes from data
#python plotterDataFakesBDT.py $analyzer $luminosity $jobid $analtype $kinplots
#
##plot same sign Using fakes from data
#python plotterDataFakesMethod2BDT.py $analyzer $luminosity $jobid $analtype $kinplots
#
##plot anti-isolated OS region
#python plotterAntiIsolatedOSBDT.py $analyzer $luminosity $jobid $analtype $kinplots
#
##plot anti-isolated SS region
#python plotterAntiIsolatedSSBDT.py $analyzer $luminosity $jobid $analtype $kinplots
#
##plot weighted anti-isolated SS region all fakes -- muon + electron -electron muon
#python plotterAntiIsolatedWeightedSSBDT.py  $analyzer $luminosity $jobid $analtype $kinplots 
#
##plot weighted anti-isolated OS region
#python plotterAntiIsolatedWeightedOSBDT.py  $analyzer $luminosity $jobid $analtype $kinplots
#
##plot weighted anti-isolated OS region electron fakes
#python plotterAntiIsolatedWeightedElectronOSBDT.py $analyzer $luminosity $jobid $analtype $kinplots
#
#
##plot weighted anti-isolated SS region electron fakes
#python plotterAntiIsolatedWeightedElectronSSBDT.py $analyzer $luminosity $jobid $analtype $kinplots
#
#
##plot weighted anti-isolated OS region muon fakes                                                                                                                               
#python plotterAntiIsolatedWeightedMuonOSBDT.py $analyzer $luminosity $jobid $analtype $kinplots
#
#
##plot weighted anti-isolated SS region muon fakes                                                                                                                               
#python plotterAntiIsolatedWeightedMuonSSBDT.py $analyzer $luminosity $jobid $analtype $kinplots
#
#
##plot weighted anti-isolated OS region muon electron fakes                                                                                                                               
#python plotterAntiIsolatedWeightedMuonElectronOS.py $analyzer $luminosity $jobid $analtype $kinplots
#
#
##plot weighted anti-isolated SS region muon electron fakes                                                                                                                               
#python plotterAntiIsolatedWeightedMuonElectronSS.py $analyzer $luminosity $jobid $analtype $kinplots
