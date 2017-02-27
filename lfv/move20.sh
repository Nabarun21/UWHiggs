
mv GluGluHToTauTau_M125_13TeV_powheg_pythia8.root  ggHTauTau.root             
mv VBFHToTauTau_M125_13TeV_powheg_pythia8.root vbfHTauTau.root
mv TT_TuneCUETP8M1_13TeV-powheg-pythia8*.root TT.root

mv GluGlu_LFV_HToMuTau_M125_13TeV_powheg_pythia8.root LFVGG125.root
mv VBF_LFV_HToMuTau_M125_13TeV_powheg_pythia8.root LFVVBF125.root

hadd -f Dibosons.root WW_TuneCUETP8M1_13TeV-pythia8.root WZ_TuneCUETP8M1_13TeV-pythia8.root ZZ_TuneCUETP8M1_13TeV-pythia8.root 

hadd -f DY.root DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root

hadd -f WJETSMC.root WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root

hadd -f  data_obs.root   data_SingleMuon_Run2016B*.root  data_SingleMuon_Run2016C*.root  data_SingleMuon_Run2016D*.root  data_SingleMuon_Run2016E*.root data_SingleMuon_Run2016F*.root

hadd -f WG.root WGstarToLNuEE_012Jets_13TeV-madgraph.root WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root  WGstarToLNuMuMu_012Jets_13TeV-madgraph.root

hadd -f T.root ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root


rm WW_TuneCUETP8M1_13TeV-pythia8.root
rm WZ_TuneCUETP8M1_13TeV-pythia8.root
rm ZZ_TuneCUETP8M1_13TeV-pythia8.root
rm WG*mad*
rm WG*amcat*
rm ST*
rm DY*Jets*
rm W*Jets*
rm data_SingleMuon_*