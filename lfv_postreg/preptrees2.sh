hadd -f BackgroundTree2.root results/SMHTT_aug16_v2/LFVHEMuAnalyzerMVA_makeBDTtreesvbfunchanged/*HToTau*   results/SMHTT_aug16_v2/LFVHEMuAnalyzerMVA_makeBDTtreesvbfunchanged/*TuneCUETP8M1*
hadd -f SignalTree2.root results/SMHTT_aug16_v2/LFVHEMuAnalyzerMVA_makeBDTtreesvbfunchanged/*LFV*
hadd -f trainingTree2.root BackgroundTree2.root SignalTree2.root
