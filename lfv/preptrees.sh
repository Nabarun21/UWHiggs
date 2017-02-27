hadd -f SignalTree.root results/Oct30/LFVHEMuAnalyzerMVAmakeBDTtrees/GluGlu_LFV_HToMuTau_M125_13TeV_powheg_pythia8.root   results/Oct30/LFVHEMuAnalyzerMVAmakeBDTtrees/VBF_LFV_HToMuTau_M125_13TeV_powheg_pythia8.root

hadd -f BackgroundTreeTT.root results/Oct30/LFVHEMuAnalyzerMVAmakeBDTtrees/TT_TuneCUETP8M1_13TeV-powheg-pythia8.root
hadd -f trainingTreeTT.root BackgroundTreeTT.root SignalTree.root

hadd -f BackgroundTreeDY.root results/Oct30/LFVHEMuAnalyzerMVAmakeBDTtrees/DY*JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root
hadd -f trainingTreeDY.root BackgroundTreeDY.root SignalTree.root


hadd -f BackgroundTree.root BackgroundTreeDY.root BackgroundTreeTT.root
hadd -f trainingTree.root BackgroundTree.root SignalTree.root
