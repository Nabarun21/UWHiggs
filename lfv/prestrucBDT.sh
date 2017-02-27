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
    esac
done

echo $analyzer
echo $luminosity

echo $jobid

#get QCD (data-MC) in ss *2.30(SF)
python getQCDBDT.py $analyzer $luminosity $jobid
python getQCDforcombineBDT.py $analyzer $luminosity $jobid   #get the same without sysuncertainty in histos, combine script already adds them

#copy results into current directory in the form wanted
cp -r results/$jobid/LFVHEMuAnalyzerMVA$analyzer LFVHEMuAnalyzerMVA$analyzer$luminosity
cp move.sh LFVHEMuAnalyzerMVA$analyzer$luminosity 
cd LFVHEMuAnalyzerMVA$analyzer$luminosity
source move.sh
rm move.sh
cd -


#construct  Fakes histos from anti-isolated region
#python getFakes.py $analyzer $luminosity
#python getFakesforcombine.py $analyzer $luminosity   #get the same without sysuncertainty in histos, combine script already adds them


#folder for plotting
cp -r LFVHEMuAnalyzerMVA$analyzer$luminosity LFVHEMuAnalyzerMVA$analyzer$luminosity'plot'
#cp FAKES.root LFVHEMuAnalyzerMVA$analyzer$luminosity'plot'
cp QCDBDT.root LFVHEMuAnalyzerMVA$analyzer$luminosity'plot'/QCD.root

#folders for constructing limits
#cp -r LFVHEMuAnalyzerMVA$analyzer$luminosity LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromdata'
cp -r LFVHEMuAnalyzerMVA$analyzer$luminosity LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromMC'
#cp FAKESforcombine.root LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromdata'/FAKES.root
#rm LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromdata'/WJETSMC.root
#rm LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromdata'/QCD.root
cp QCDforcombineBDT.root LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromMC'/QCD.root



#in plotting folder remove neg bins + add systematic uncertainities
python prepHistosBDT.py $analyzer $luminosity

#plotting replace weight files by current lumi (I am a superfluos remnant feature from ice ages, get rid of me)
python replaceFileWeightBy1.py $luminosity


source jobid.sh
#plot opposite sign region Using fakes from MC
#python plotterMC.py $analyzer $luminosity $jobid

#plot same sign Using fakes from MC
#python plotterMCss.py $analyzer $luminosity $jobid

#plot same sign Using fakes from data
#python plotterDataFakes.py $analyzer $luminosity $jobid

#plot anti-isolated OS region
#python plotterAntiIsolatedOS.py $analyzer $luminosity $jobid

#plot anti-isolated SS region
#python plotterAntiIsolatedSS.py $analyzer $luminosity $jobid

#plot weighted anti-isolated SS region
#python plotterAntiIsolatedWeightedSS.py  $analyzer $luminosity $jobid

#plot weighted anti-isolated OS region
#python plotterAntiIsolatedWeightedOS.py  $analyzer $luminosity $jobid


python plotterBDT.py $analyzer $luminosity $jobid