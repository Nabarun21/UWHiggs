from MuTauTree import MuTauTree
import os
import ROOT
import math
import glob
import array
import baseSelections as selections
import FinalStateAnalysis.PlotTools.pytree as pytree
from FinalStateAnalysis.PlotTools.decorators import  memo_last
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
from math import sqrt, pi, cos

def collmass(row, met, metPhi):
    ptnu = met*cos(deltaPhi(metPhi, row.tPhi))
    visfrac = row.tPt/(row.tPt+abs(ptnu))
    #print met, cos(deltaPhi(metPhi, row.tPhi)), ptnu, visfrac
    return (row.m_t_Mass / sqrt(visfrac))

def deltaPhi(phi1, phi2):
    PHI = abs(phi1-phi2)
    if PHI<=pi:
        return PHI
    else:
        return 2*pi-PHI
def deltaR(phi1, ph2, eta1, eta2):
    deta = eta1 - eta2
    dphi = abs(phi1-phi2)
    if (dphi>pi) : dphi = 2*pi-dphi
    return sqrt(deta*deta + dphi*dphi);

class LFVHMuTauAnalyzer(MegaBase):
    tree = 'mt/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel='MT'
        super(LFVHMuTauAnalyzer, self).__init__(tree, outfile, **kwargs)
        self.tree = MuTauTree(tree)
        self.out=outfile
        self.histograms = {}

    @staticmethod 
    def tau_veto(row):
        if not row.tAntiMuonTight2 or not row.tAntiElectronMVA3Loose or not row.tDecayFinding :
            return False

    @staticmethod
    def obj1_matches_gen(row):
        return row.eGenPdgId == -1*row.eCharge*11
    @staticmethod 
    def obj3_matches_gen(row):
        return t.genDecayMode != -2 

    
 
## 
    def begin(self):

        processtype=['gg', 'vbf']
        threshold=['mpt0', 'mpt40']
        sign=['os', 'ss']
        jetN = [0, 1, 2, 3]
        for i in sign:
            for j in processtype:
                for k in threshold:
                    for jn in jetN: 

                        folder =i+'/'+j+'/'+k +'/'+str(jn)
                   
                        self.book(folder,"tPt", "tau p_{T}", 200, 0, 200)
                        self.book(folder,"tPhi", "tau phi", 100, -3.2, 3.2)
                        self.book(folder,"tEta", "tau eta",  50, -2.5, 2.5)
                        
                        self.book(folder,"mPt", "m p_{T}", 200, 0, 200)
                        self.book(folder,"mPhi", "m phi",  100, -3.2, 3.2)
                        self.book(folder,"mEta", "m eta", 50, -2.5, 2.5)
                        
                        self.book(folder, "mt_DeltaPhi", "m-tau DeltaPhi" , 50, 0, 3.2)
                        self.book(folder, "mt_DeltaR", "m-tau DeltaR" , 50, 0, 3.2)
                        
                        self.book(folder, "h_collmass_pfmet",  "h_collmass_pfmet",  32, 0, 320)
                        self.book(folder, "h_collmass_mvamet",  "h_collmass_mvamet",  32, 0, 320)
                        
                        self.book(folder, "h_vismass",  "h_vismass",  32, 0, 320)
                        
                        
                        self.book(folder, "tPFMET_DeltaPhi", "tau-PFMET DeltaPhi" , 50, 0, 3.2)
                        self.book(folder, "tPFMET_Mt", "tau-PFMET M_{T}" , 200, 0, 200)
                        self.book(folder, "tMVAMET_DeltaPhi", "tau-MVAMET DeltaPhi" , 50, 0, 3.2)
                        self.book(folder, "tMVAMET_Mt", "tau-MVAMET M_{T}" , 200, 0, 200)
                        
                        self.book(folder, "mPFMET_DeltaPhi", "m-PFMET DeltaPhi" , 50, 0, 3.2)
                        self.book(folder, "mPFMET_Mt", "m-PFMET M_{T}" , 200, 0, 200)
                        self.book(folder, "mMVAMET_DeltaPhi", "m-MVAMET DeltaPhi" , 50, 0, 3.2)
                        self.book(folder, "mMVAMET_Mt", "m-MVAMET M_{T}" , 200, 0, 200)
                        
                        self.book(folder, "jetN_20", "Number of jets, p_{T}>20", 10, -0.5, 9.5) 
                        self.book(folder, "jetN_30", "Number of jets, p_{T}>30", 10, -0.5, 9.5) 

    def fill_histos(self, row, folder='os/gg/mpt0/0', fakeRate = False):

        histos = self.histograms

        histos[folder+'/tPt'].Fill(row.tPt)
        histos[folder+'/tEta'].Fill(row.tEta)
        histos[folder+'/tPhi'].Fill(row.tPhi) 

        histos[folder+'/mPt'].Fill(row.mPt)
        histos[folder+'/mEta'].Fill(row.mEta)
        histos[folder+'/mPhi'].Fill(row.mPhi)

        histos[folder+'/mt_DeltaPhi'].Fill(deltaPhi(row.mPhi, row.tPhi))
        histos[folder+'/mt_DeltaR'].Fill(row.m_t_DR)
        
        
        histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.pfMetEt, row.pfMetPhi))
        histos[folder+'/h_collmass_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi))

        histos[folder+'/h_vismass'].Fill(row.m_t_Mass)
               

        histos[folder+'/mPFMET_DeltaPhi'].Fill(deltaPhi(row.mPhi, row.pfMetPhi))
        histos[folder+'/mMVAMET_DeltaPhi'].Fill(deltaPhi(row.mPhi, row.mva_metPhi))
        histos[folder+'/mPFMET_Mt'].Fill(row.mMtToPFMET)
        histos[folder+'/mMVAMET_Mt'].Fill(row.mMtToMVAMET)

        histos[folder+'/tPFMET_DeltaPhi'].Fill(deltaPhi(row.tPhi, row.pfMetPhi))
        histos[folder+'/tMVAMET_DeltaPhi'].Fill(deltaPhi(row.tPhi, row.mva_metPhi))
        histos[folder+'/tPFMET_Mt'].Fill(row.tMtToPFMET)
        histos[folder+'/tMVAMET_Mt'].Fill(row.tMtToMVAMET)


        histos[folder+'/jetN_20'].Fill(row.jetVeto20) 
        histos[folder+'/jetN_30'].Fill(row.jetVeto30) 
        


        

    def process(self):
        for row in self.tree:
            
        #for i, row in enumerate(self.tree):
            #if  i >= 100:
                #return
            if not bool(row.isoMu24eta2p1Pass) : continue
            
 
            sign = 'ss' if row.m_t_SS else 'os'
            processtype = '' ## use a line as for sign when the vbf when selections are defined            
            ptthreshold = [0,40]
            if row.vbfJetVeto30:
                processtype = 'vbf'
            else:
                processtype ='gg'##changed from 20
            
            jn = row.jetVeto30
            if jn > 3 : jn = 3
            if jn != 0 and row.bjetCSVVeto30!=0 : continue 

            if not selections.muSelection(row, 'm'): continue
            if not selections.lepton_id_iso(row, 'm', 'mutauiso012'): continue#'mid12Tight_mtauiso012'): continue
            if not selections.tauSelection(row, 't'): continue
            
            #if not selections.vetos(row) : continue
            
            if not row.tTightIso3Hits : continue
            if not row.tAntiElectronMVA3Loose : continue
            if not row.tAntiMuonTight2 : continue
            #print row.mPt     
            if row.tauVetoPt20EleTight3MuLoose : continue #to be change to Pt20EleLooseMutight
            if row.muVetoPt5IsoIdVtx : continue
            if row.eVetoCicTightIso : continue                     
            
           
            for j in ptthreshold:
                folder = sign+'/'+processtype+'/mpt'+str(j)+'/'+str(int(jn))
                    
                if row.mPt < j : continue
                 
               # if row.jetVeto30_DR05 > 0 : continue 
                    
                self.fill_histos(row, folder)
                    

           
            
            
    def finish(self):
        self.write_histos()


