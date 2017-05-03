from EMuTree import EMuTree
import os
import ROOT
from math import *
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
    
class LFVHAnalyzeGEN(MegaBase):
    tree = 'em/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel='EMu'
        super(LFVHAnalyzeGEN, self).__init__(tree, outfile, **kwargs)
        self.out=outfile
        self.histograms = {}


    def begin(self):

        #GenQuantities 
        self.book('gen',"eGenMotherPdgId", "gen Mother pdgId", 200, 0, 200)
        self.book('gen',"eGenEta", "gen Eta, electron", 50, -2.5, 2.5)
        self.book('gen',"eGenPhi", "gen Phi, electron", 100, -3.2, 3.2)
        self.book('gen',"eGenEnergy", "gen Energy, electron", 40, 0, 400)
        self.book('gen',"eGenPt", "gen p_{T}, electron", 40, 0, 400)
        
        self.book('gen',"mGenEta", "gen Eta, mu", 50, -2.5, 2.5)
        self.book('gen',"mGenPhi", "gen Phi, mu",  100, -3.2, 3.2)
        self.book('gen',"mGenEnergy", "gen Energy, mu", 40, 0, 400)
        self.book('gen',"mGenPt", "gen p_{T}, mu", 40 , 0, 400)
#        self.book('gen',"mGenDecayMode", "gen Mu, decay mode", 20, 0, 20)
        self.book('gen',"emGenDeltaPhi", " gen e mu delta phi",  50, 0, 3.2)

        self.book('gen',"eGenMotherdeltaPhi", "gen e tau(e_mother), delta phi",  50, 0, 3.2)
        self.book('gen',"eGenMotherdeltaR", "gen e tau(e_mother), delta R",  100, 0, 6.5)

        
        self.book('gen',"eGenMotherPdgId_all", "gen Mother pdfId, all electrons", 200, 0, 200)
        self.book('gen',"eGenEta_all", "gen Eta, all electrons", 50, -2.5, 2.5)
        self.book('gen',"eGenPhi_all", "gen Phi, all electrons",  100, -3.2, 3.2)
        self.book('gen',"eGenPt_all", "gen p_{T}, all electrons", 40, 0, 400)
        self.book('gen',"eGenEnergy_all", "gen Energy, all electrons", 40, 0, 400)
    
        self.book('gen',"mGenEta_all", "gen Eta, all mus", 50, -2.5, 2.5)
        self.book('gen',"mGenPhi_all", "gen Phi, all mus",  100, -3.2, 3.2)
        self.book('gen',"mGenEnergy_all", "gen Energy, all mus", 40, 0, 400)
        self.book('gen',"mGenPt_all", "gen p_{T}, all mus",  40, 0, 400)
 #       self.book('gen',"mGenDecayMode_all", "all gen Mus decay mode", 20, 0, 20)
        self.book('gen',"emGenDeltaPhi_all", "all gen e mu delta phi",  50, 0, 3.2)

        self.book('gen', 'higgsPt', 'higgs p_{T}', 40, 0, 400); 

            
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
  #      histos[folder+'/mGenDecayMode_all'].Fill(row.mGenDecayMode)
 
        histos[folder+'/emGenDeltaPhi_all'].Fill(deltaPhi(row.eGenPhi, row.mGenPhi))

 
#        if row.eGenMotherPdgId == 25: 
#        print "electron mother's particle ID=  ",row.eGenMotherPdgId

        if row.eComesFromHiggs == True:
            if abs(row.eGenMotherPdgId)==15:
                print "found tau mother"
                histos[folder+'/eGenMotherdeltaPhi'].Fill(deltaPhi(row.eGenPhi, row.eGenMotherPhi))
                histos[folder+'/eGenMotherdeltaR'].Fill(sqrt(deltaPhi(row.eGenPhi, row.eGenMotherPhi)**2+(eGenEta-eGenMotherEta)**2))

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
   #         histos[folder+'/mGenDecayMode'].Fill(row.mGenDecayMode)
            
        if row.eComesFromHiggs == True and row.mComesFromHiggs == True: 
            #print row.eGenMotherPdgId, " ", row.mGenMotherPdgId
            histos[folder+'/emGenDeltaPhi'].Fill(deltaPhi(row.eGenPhi, row.mGenPhi))
            
            histos[folder+'/higgsPt'].Fill(row.eGenMotherPt) 
            
    
 
    def process(self):
        for row in self.tree:
        #os.environ['megatarget'] ## to get the sample name
        #for i, row in enumerate(self.tree):
         #   if  i >= 1000:
          #      return
            #print row.eGenMotherPdgId, row.eGenPdgId, row.eComesFromHiggs, row.tComesFromHiggs
            #if row.eGenMotherPdgId==25:
            #print row.eComesFromHiggs, row.eGenMotherPdgId, row.eGenPdgId
            self.fill_histos(row, 'gen')
    def finish(self):
        self.write_histos()
