hadd -f BackgroundTreeD.root results/SMHTT_aug16_v2/LFVHEMuAnalyzerMVA_makeBDTtrees/*DY*
hadd -f SignalTreeD.root results/SMHTT_aug16_v2/LFVHEMuAnalyzerMVA_makeBDTtrees/*LFV*
hadd -f trainingTreeD.root BackgroundTreeD.root SignalTreeD.root
