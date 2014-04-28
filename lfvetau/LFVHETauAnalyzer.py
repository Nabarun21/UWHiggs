from ETauTree import ETauTree
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
    return (row.e_t_Mass / sqrt(visfrac))

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

class LFVHETauAnalyzer(MegaBase):
    tree = 'et/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel='ET'
        super(LFVHETauAnalyzer, self).__init__(tree, outfile, **kwargs)
        self.tree = ETauTree(tree)
        self.out=outfile
        self.histograms = {}

    @staticmethod 
    def tau_veto(row):
        if not row.tAntiMuonLoose2 or not row.tAntiElectronMVA3Tight or not row.tDecayFinding :
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
        threshold=['ept0', 'ept40']
        sign=['os', 'ss']
        jetN = [0, 1, 2, 3]
        for i in sign:
            for j in processtype:
                for k in threshold:
                    for jn in jetN: 
                        folder=[]

                        folder.append(i+'/'+j+'/'+k +'/'+str(jn))
                        folder.append(i+'/'+j+'/'+k +'/'+str(jn)+'/selected')
                        
                        for f in folder: 
                            self.book(f,"tPt", "tau p_{T}", 200, 0, 200)
                            self.book(f,"tPhi", "tau phi", 100, -3.2, 3.2)
                            self.book(f,"tEta", "tau eta",  50, -2.5, 2.5)
                    
                            self.book(f,"ePt", "e p_{T}", 200, 0, 200)
                            self.book(f,"ePhi", "e phi",  100, -3.2, 3.2)
                            self.book(f,"eEta", "e eta", 50, -2.5, 2.5)
                            
                            self.book(f, "et_DeltaPhi", "e-tau DeltaPhi" , 50, 0, 3.2)
                            self.book(f, "et_DeltaR", "e-tau DeltaR" , 50, 0, 3.2)
                            
                            self.book(f, "h_collmass_pfmet",  "h_collmass_pfmet",  32, 0, 320)
                            self.book(f, "h_collmass_mvamet",  "h_collmass_mvamet",  32, 0, 320)
                            
                            self.book(f, "h_vismass",  "h_vismass",  32, 0, 320)
                        
                            
                            self.book(f, "tPFMET_DeltaPhi", "tau-PFMET DeltaPhi" , 50, 0, 3.2)
                            self.book(f, "tPFMET_Mt", "tau-PFMET M_{T}" , 200, 0, 200)
                            self.book(f, "tMVAMET_DeltaPhi", "tau-MVAMET DeltaPhi" , 50, 0, 3.2)
                            self.book(f, "tMVAMET_Mt", "tau-MVAMET M_{T}" , 200, 0, 200)
                            
                            self.book(f, "ePFMET_DeltaPhi", "e-PFMET DeltaPhi" , 50, 0, 3.2)
                            self.book(f, "ePFMET_Mt", "e-PFMET M_{T}" , 200, 0, 200)
                            self.book(f, "eMVAMET_DeltaPhi", "e-MVAMET DeltaPhi" , 50, 0, 3.2)
                            self.book(f, "eMVAMET_Mt", "e-MVAMET M_{T}" , 200, 0, 200)
                            
                            self.book(f, "jetN_20", "Number of jets, p_{T}>20", 10, -0.5, 9.5) 
                            self.book(f, "jetN_30", "Number of jets, p_{T}>30", 10, -0.5, 9.5) 

                    
    def fill_histos(self, row, folder='os/gg/ept0/0', fakeRate = False):

        histos = self.histograms

        histos[folder+'/tPt'].Fill(row.tPt)
        histos[folder+'/tEta'].Fill(row.tEta)
        histos[folder+'/tPhi'].Fill(row.tPhi) 

        histos[folder+'/ePt'].Fill(row.ePt)
        histos[folder+'/eEta'].Fill(row.eEta)
        histos[folder+'/ePhi'].Fill(row.ePhi)

        histos[folder+'/et_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.tPhi))
        histos[folder+'/et_DeltaR'].Fill(row.e_t_DR)
        
        
        histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.pfMetEt, row.pfMetPhi))
        histos[folder+'/h_collmass_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi))

        histos[folder+'/h_vismass'].Fill(row.e_t_Mass)
               

        histos[folder+'/ePFMET_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.pfMetPhi))
        histos[folder+'/eMVAMET_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.mva_metPhi))
        histos[folder+'/ePFMET_Mt'].Fill(row.eMtToPFMET)
        histos[folder+'/eMVAMET_Mt'].Fill(row.eMtToMVAMET)

        histos[folder+'/tPFMET_DeltaPhi'].Fill(deltaPhi(row.tPhi, row.pfMetPhi))
        histos[folder+'/tMVAMET_DeltaPhi'].Fill(deltaPhi(row.tPhi, row.mva_metPhi))
        histos[folder+'/tPFMET_Mt'].Fill(row.tMtToPFMET)
        histos[folder+'/tMVAMET_Mt'].Fill(row.tMtToMVAMET)


        histos[folder+'/jetN_20'].Fill(row.jetVeto20) 
        histos[folder+'/jetN_30'].Fill(row.jetVeto30) 



        

    def process(self):
        for row in self.tree:
#        for i, row in enumerate(self.tree):
#            if  i >= 100:
#                return
#           if not bool(row.isoMu24eta2p1Pass) : continue
 
            sign = 'ss' if row.e_t_SS else 'os'
            processtype = '' ## use a line as for sign when the vbf when selections are defined            
            ptthreshold = [0,40]
            if row.vbfJetVeto30:
                processtype = 'vbf'
            else:
                processtype ='gg'##changed from 20
            jn = row.jetVeto30
            if jn > 3 : jn = 3
            if not bool(row.singleEPass) : continue
            if jn != 0 and row.bjetCSVVeto30!=0 : continue 
            
            if not selections.eSelection(row, 'e'): continue
            if not selections.lepton_id_iso(row, 'e', 'eid12Loose_etauiso012'): continue
            if not selections.tauSelection(row, 't'): continue
            #if not selections.vetos(row) : continue
            if not row.tTightIso3Hits : continue
            if not row.tAntiElectronMVA3Tight : continue
            if not row.tAntiMuon2Loose : continue
            #if row.tJetOverlap  : continue 
            #if  row.eJetOverlap : continue
                
            if row.tauVetoPt20EleTight3MuLoose : continue 
            if row.muVetoPt5IsoIdVtx : continue
            if row.eVetoCicLooseIso : continue # change it with Loose
  
            for j in ptthreshold:
                folder = sign+'/'+processtype+'/ept'+str(j)+'/'+str(int(jn))
                
                if row.ePt < j : continue
                  
                self.fill_histos(row, folder)
                #selections
                if jn == 0 :
                    if row.tPt < 35: continue 
                    if row.ePt < 40 : continue
                    if deltaPhi(row.ePhi, row.tPhi) < 2.7 : continue
                    if row.tMtToPFMET > 50 : continue
                if jn == 1 :
                    if row.tPt < 40: continue 
                    if row.ePt < 35 : continue
                    if row.tMtToPFMET > 35 : continue
                if jn == 2 :
                    if row.tPt < 40: continue 
                    if row.ePt < 30 : continue # no cut as only electrons with pt>30 are in the ntuples
                    if row.tMtToPFMET > 35 : continue
                    if row.vbfMass < 550 : continue
                    if row.vbfDeta < 3.5 : continue
                folder = sign+'/'+processtype+'/ept'+str(j)+'/'+str(int(jn))+'/selected'
                self.fill_histos(row, folder)
                    
             
            
    def finish(self):
        self.write_histos()


