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


rm -r LFVHEMuAnalyzerMVA$analyzer$luminosity*

#get QCD (data-MC) in ss *2.30(SF)
python getQCDBDT.py $analyzer $luminosity $jobid $analtype 
python getQCDforcombineBDT.py $analyzer $luminosity $jobid $analtype   #get the same without sysuncertainty in histos, combine script already adds them

#copy results into current directory in the form wanted
cp -r results/$jobid/LFVHEMuAnalyzerMVA$analyzer LFVHEMuAnalyzerMVA$analyzer$luminosity
cp move.sh LFVHEMuAnalyzerMVA$analyzer$luminosity 
cd LFVHEMuAnalyzerMVA$analyzer$luminosity
source move.sh
rm move.sh
cd -


#construct  Fakes histos from anti-isolated region
python getFakesBDT.py $analyzer $luminosity $jobid $analtype
python getFakesforcombineBDT.py $analyzer $luminosity   $jobid $analtype  #get the same without sysuncertainty in histos, combine script already adds them

#construct  Fakes histos from anti-isolated region:method 2 
python getFakesmethod2BDT.py $analyzer $luminosity $jobid $analtype  
python getFakesmethod2forcombineBDT.py $analyzer $luminosity $jobid $analtype   #get the same without sysuncertainty in histos, combine script already adds them


#folder for plotting
cp -r LFVHEMuAnalyzerMVA$analyzer$luminosity LFVHEMuAnalyzerMVA$analyzer$luminosity'plot'
cp FAKES$analtype.root LFVHEMuAnalyzerMVA$analyzer$luminosity'plot'/FAKES.root
cp FAKESmethod2$analtype.root LFVHEMuAnalyzerMVA$analyzer$luminosity'plot'/FAKESmethod2.root
cp QCD$analtype.root LFVHEMuAnalyzerMVA$analyzer$luminosity'plot'/QCD.root

#folders for constructing limits
cp -r LFVHEMuAnalyzerMVA$analyzer$luminosity LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromdata'
cp -r LFVHEMuAnalyzerMVA$analyzer$luminosity LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromdatamethod2'
cp -r LFVHEMuAnalyzerMVA$analyzer$luminosity LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromMC'
cp FAKESforcombine$analtype.root LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromdata'/FAKES.root
cp FAKESmethod2forcombine$analtype.root LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromdatamethod2'/FAKES.root
rm LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromdata'/WJETSMC.root
#rm LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromdata'/QCD.root
cp QCDforcombine$analtype.root LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromMC'/QCD.root
rm LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromMC'/*ETau*
rm LFVHEMuAnalyzerMVA$analyzer$luminosity'fakesfromdata'/*ETau*

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
