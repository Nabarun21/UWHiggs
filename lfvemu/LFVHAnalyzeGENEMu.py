from ETauTree import ETauTree
import os
import ROOT
import math
import glob
import array
from FinalStateAnalysis.PlotTools.decorators import  memo_last
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
from math import sqrt, pi

def deltaPhi(phi1, phi2):
    PHI = abs(phi1-phi2)
    if PHI<=pi:
        return PHI
    else:
        return 2*pi-PHI
    
class LFVHAnalyzeGENEMu(MegaBase):
    tree = 'em/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel='ET'
        super(LFVHAnalyzeGENEMu, self).__init__(tree, outfile, **kwargs)
        self.out=outfile
        self.histograms = {}


    def begin(self):

        #GenQuantities 
        self.book('gen',"eGenMotherPdgId", "gen Mother pdgId", 200, 0, 200)
        self.book('gen',"eGenEta", "gen Eta, electron", 50, -2.5, 2.5)
        self.book('gen',"eGenPhi", "gen Phi, electron", 100, -3.2, 3.2)
        self.book('gen',"eGenEnergy", "gen Energy, electron", 40, 0, 200)
        self.book('gen',"eGenPt", "gen p_{T}, electron", 40, 0, 200)
        
        self.book('gen',"mGenEta", "gen Eta, muon", 50, -2.5, 2.5)
        self.book('gen',"mGenPhi", "gen Phi, muon",  100, -3.2, 3.2)
        self.book('gen',"mGenEnergy", "gen Energy, muon", 40, 0, 200)
        self.book('gen',"mGenPt", "gen p_{T}, muon", 40 , 0, 200)
 
        self.book('gen',"emGenDeltaPhi", "gen emu, delta phi",  50, 0, 3.2)

        
        self.book('gen',"eGenMotherPdgId_all", "gen Mother pdfId, all electrons", 200, 0, 200)
        self.book('gen',"eGenEta_all", "gen Eta, all electrons", 50, -2.5, 2.5)
        self.book('gen',"eGenPhi_all", "gen Phi, all electrons",  100, -3.2, 3.2)
        self.book('gen',"eGenPt_all", "gen p_{T}, all electrons", 40, 0, 200)
        self.book('gen',"eGenEnergy_all", "gen Energy, all electrons", 40, 0, 200)
    
        self.book('gen',"mGenEta_all", "gen Eta, all muons", 50, -2.5, 2.5)
        self.book('gen',"mGenPhi_all", "gen Phi, all muons",  100, -3.2, 3.2)
        self.book('gen',"mGenEnergy_all", "gen Energy, all muons", 40, 0, 200)
        self.book('gen',"mGenPt_all", "gen p_{T}, all muons",  40, 0, 200)
#        self.book('gen',"mGenDecayMode_all", "all gen Taus decay mode", 20, 0, 20)
        self.book('gen',"emGenDeltaPhi_all", "all gen e mu delta phi",  50, 0, 3.2)

        self.book('gen', 'higgsPt', 'higgs p_{T}', 40, 0, 200); 

            
    def fill_histos(self, row, folder='gen', fakeRate = False):
        histos = self.histograms

        histos[folder+'/eGenMotherPdgId_all'].Fill(row.eGenMotherPdgId) 
        histos[folder+'/eGenEta_all'].Fill(row.eGenEta) 
        histos[folder+'/eGenPhi_all'].Fill(row.eGenPhi)
        histos[folder+'/eGenPt_all'].Fill(row.eGenPt)
        histos[folder+'/eGenEnergy_all'].Fill(row.eGenEnergy)
        
        histos[folder+'/mGenEta_all'].Fill(row.mGenEta)
        histos[folder+'/mGenPhi_all'].Fill(row.mGenPhi)
        histos[folder+'/mGenPt_all'].Fill(row.mGenPt)
        histos[folder+'/mGenEnergy_all'].Fill(row.mGenEnergy)
#        histos[folder+'/mGenDecayMode_all'].Fill(row.mGenDecayMode)
 
        histos[folder+'/emGenDeltaPhi_all'].Fill(deltaPhi(row.eGenPhi, row.mGenPhi))

 
#        if row.eGenMotherPdgId == 25: 
        if row.eComesFromHiggs == True:
            histos[folder+'/eGenMotherPdgId'].Fill(row.eGenMotherPdgId) 
            histos[folder+'/eGenEta'].Fill(row.eGenEta) 
            histos[folder+'/eGenPhi'].Fill(row.eGenPhi)
            histos[folder+'/eGenPt'].Fill(row.eGenPt)
            histos[folder+'/eGenEnergy'].Fill(row.eGenEnergy)
        if row.mComesFromHiggs  == True:            
            histos[folder+'/mGenEta'].Fill(row.mGenEta)
            histos[folder+'/mGenPhi'].Fill(row.mGenPhi)
            histos[folder+'/mGenPt'].Fill(row.mGenPt)
            histos[folder+'/mGenEnergy'].Fill(row.mGenEnergy)
#            histos[folder+'/mGenDecayMode'].Fill(row.mGenDecayMode)
            
        if row.eComesFromHiggs == True and row.mComesFromHiggs == True: 
            #print row.eGenMotherPdgId, " ", row.mGenMotherPdgId
            histos[folder+'/emGenDeltaPhi'].Fill(deltaPhi(row.eGenPhi, row.mGenPhi))
            
            histos[folder+'/higgsPt'].Fill(sqrt((row.eGenPx+row.mGenPx)**2 + (row.eGenPy+row.mGenPy)**2)) # correct  for LFVHiggsToETau only, add the neutrinos
            
    
 
    def process(self):
        for row in self.tree:
        #os.environ['megatarget'] ## to get the sample name
        #for i, row in enumerate(self.tree):
         #   if  i >= 1000:
          #      return
            self.fill_histos(row, 'gen')
    def finish(self):
        self.write_histos()
