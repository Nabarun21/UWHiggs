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
    
class LFVHAnalyzeGENMuTau(MegaBase):
    tree = 'mt/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel='ET'
        super(LFVHAnalyzeGENMuTau, self).__init__(tree, outfile, **kwargs)
        self.out=outfile
        self.histograms = {}


    def begin(self):
        folder = ['gen', 'gen_presel', 'gen_sel']
        #GenQuantities 
        for d in folder:
            
            self.book(d,"mGenMotherPdgId", "gen Mother pdgId", 200, 0, 200)
            self.book(d,"mGenEta", "gen Eta, muon", 50, -2.5, 2.5)
            self.book(d,"mGenPhi", "gen Phi, muon", 100, -3.2, 3.2)
            self.book(d,"mGenEnergy", "gen Energy, muon", 40, 0, 200)
            self.book(d,"mGenPt", "gen p_{T}, muon", 40, 0, 200)
        
            self.book(d,"tGenEta", "gen Eta, tau", 50, -2.5, 2.5)
            self.book(d,"tGenPhi", "gen Phi, tau",  100, -3.2, 3.2)
            self.book(d,"tGenEnergy", "gen Energy, tau", 40, 0, 200)
            self.book(d,"tGenPt", "gen p_{T}, tau", 40 , 0, 200)
            self.book(d,"tGenDecayMode", "gen Tau, decay mode", 20, 0, 20)
            self.book(d,"mtGenDeltaPhi", "gen e tau, delta phi",  50, 0, 3.2)

            
            self.book(d,"mGenMotherPdgId_all", "gen Mother pdfId, all muons", 200, 0, 200)
            self.book(d,"mGenEta_all", "gen Eta, all muons", 50, -2.5, 2.5)
            self.book(d,"mGenPhi_all", "gen Phi, all muons",  100, -3.2, 3.2)
            self.book(d,"mGenPt_all", "gen p_{T}, all muons", 40, 0, 200)
            self.book(d,"mGenEnergy_all", "gen Energy, all muons", 40, 0, 200)
    
            self.book(d,"tGenEta_all", "gen Eta, all taus", 50, -2.5, 2.5)
            self.book(d,"tGenPhi_all", "gen Phi, all taus",  100, -3.2, 3.2)
            self.book(d,"tGenEnergy_all", "gen Energy, all taus", 40, 0, 200)
            self.book(d,"tGenPt_all", "gen p_{T}, all taus",  40, 0, 200)
            self.book(d,"tGenDecayMode_all", "all gen Taus decay mode", 20, 0, 20)
            self.book(d,"mtGenDeltaPhi_all", "all gen e tau delta phi",  50, 0, 3.2)
            
            self.book(d, 'higgsPt', 'higgs p_{T}', 40, 0, 200); 

            
    def fill_histos(self, row, folder='gen', fakeRate = False):
        histos = self.histograms

        histos[folder+'/mGenMotherPdgId_all'].Fill(row.mGenMotherPdgId) 
        histos[folder+'/mGenEta_all'].Fill(row.mGenEta) 
        histos[folder+'/mGenPhi_all'].Fill(row.mGenPhi)
        histos[folder+'/mGenPt_all'].Fill(row.mGenPt)
        histos[folder+'/mGenEnergy_all'].Fill(row.mGenEnergy)
        
        histos[folder+'/tGenEta_all'].Fill(row.tGenEta)
        histos[folder+'/tGenPhi_all'].Fill(row.tGenPhi)
        histos[folder+'/tGenPt_all'].Fill(row.tGenPt)
        histos[folder+'/tGenEnergy_all'].Fill(row.tGenEnergy)
        histos[folder+'/tGenDecayMode_all'].Fill(row.tGenDecayMode)
 
        histos[folder+'/mtGenDeltaPhi_all'].Fill(deltaPhi(row.mGenPhi, row.tGenPhi))

 
#        if row.mGenMotherPdgId == 25: 
        if row.mComesFromHiggs == True:
            histos[folder+'/mGenMotherPdgId'].Fill(row.mGenMotherPdgId) 
            histos[folder+'/mGenEta'].Fill(row.mGenEta) 
            histos[folder+'/mGenPhi'].Fill(row.mGenPhi)
            histos[folder+'/mGenPt'].Fill(row.mGenPt)
            histos[folder+'/mGenEnergy'].Fill(row.mGenEnergy)
        if row.tComesFromHiggs  == True:            
            histos[folder+'/tGenEta'].Fill(row.tGenEta)
            histos[folder+'/tGenPhi'].Fill(row.tGenPhi)
            histos[folder+'/tGenPt'].Fill(row.tGenPt)
            histos[folder+'/tGenEnergy'].Fill(row.tGenEnergy)
            histos[folder+'/tGenDecayMode'].Fill(row.tGenDecayMode)
            
        if row.mComesFromHiggs == True and row.tComesFromHiggs == True: 
            #print row.mGenMotherPdgId, " ", row.tGenMotherPdgId
            histos[folder+'/mtGenDeltaPhi'].Fill(deltaPhi(row.mGenPhi, row.tGenPhi))
            
            histos[folder+'/higgsPt'].Fill(sqrt((row.mGenPx+row.tGenPx)**2 + (row.mGenPy+ row.tGenPy)**2)) # correct  for LFVHiggsToETau only, add the neutrinos
            
    
 
    def process(self):
        for row in self.tree:
        #os.environ['megatarget'] ## to get the sample name
        #for i, row in enumerate(self.tree):
         #   if  i >= 1000:
          #      return
            self.fill_histos(row, 'gen')
            if row.tGenPt < 20 : continue
            if row.tGenEta < -2.5 or row.tGenEta > 2.5 : continue 
            if row.mGenPt < 10 : continue
            if row.mGenEta  < -2.5 or row.mGenEta > 2.5 : continue 
            self.fill_histos(row, 'gen_presel')
            if row.tGenPt < 25 : continue
            if row.tGenEta < -2.5 or row.tGenEta > 2.5 : continue 
            if row.mGenPt < 17 : continue
            if row.mGenEta  < -2.5 or row.mGenEta > 2.5 : continue 
            if row.mMtToMVAMET > 70 : continue
            self.fill_histos(row, 'gen_sel')


    def finish(self):
        self.write_histos()
