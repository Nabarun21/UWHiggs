hadd -f SignalTree.root results/LFV_feb18_mc/LFVHEMuAnalyzerMVAmakeBDTtrees/GluGlu_LFV_HToMuTau_M125*.root   results/LFV_feb18_mc/LFVHEMuAnalyzerMVAmakeBDTtrees/VBF_LFV_HToMuTau_M125*.root

hadd -f BackgroundTreeTT.root results/LFV_feb18_mc/LFVHEMuAnalyzerMVAmakeBDTtrees/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_v6-v1.root
hadd -f trainingTreeTT.root BackgroundTreeTT.root SignalTree.root

hadd -f BackgroundTreeDY.root results/LFV_feb18_mc/LFVHEMuAnalyzerMVAmakeBDTtrees/DY*JetsToLL*TuneCUETP8M1_13TeV-*.root 
hadd -f trainingTreeDY.root BackgroundTreeDY.root SignalTree.root


hadd -f BackgroundTree.root BackgroundTreeDY.root BackgroundTreeTT.root
hadd -f trainingTree.root BackgroundTree.root SignalTree.root
