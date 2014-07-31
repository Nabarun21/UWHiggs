##Correction Factor still to add
from EETree import EETree
import os
import ROOT
import math
import glob
import array
import baseSelections as selections
import mcCorrections
import FinalStateAnalysis.PlotTools.pytree as pytree
from FinalStateAnalysis.PlotTools.decorators import  memo_last
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
from math import sqrt, pi, sin, cos, acos, sinh

def Z(row):
    e1p=ROOT.TVector3(row.e1Pt*cos(row.e1Phi),row.e1Pt*sin(row.e1Phi),row.e1Pt*sinh(row.e1Eta))
    e2p=ROOT.TVector3(row.e2Pt*cos(row.e2Phi),row.e2Pt*sin(row.e2Phi),row.e2Pt*sinh(row.e2Eta))
    e1FourVector= ROOT.TLorentzVector(e1p, sqrt(e1p.Mag2()+row.e1Mass*row.e1Mass))
    e2FourVector= ROOT.TLorentzVector(e2p, sqrt(e2p.Mag2()+row.e2Mass*row.e2Mass))
    zFourVector = e1FourVector+e2FourVector
    return zFourVector

def zPhi(pt1, eta1, phi1, pt2, eta2, phi2):
    px1 = pt1*cos(phi1)
    py1 = pt1*sin(phi1)
    pz1 = pt1*sinh(eta1)
    px2 = pt2*cos(phi2)
    py2 = pt2*sin(phi2)
    pz2 = pt2*sinh(eta2)
    
    px = px1+px2
    py = py1+py2
    pt = sqrt(px*px+py*py)
    phi = acos(px/pt)
    return phi

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

class EEAnalyzerMVA(MegaBase):
    tree = 'ee/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel='EE'
        super(EEAnalyzerMVA, self).__init__(tree, outfile, **kwargs)
        self.tree = EETree(tree)
        self.out=outfile
        self.histograms = {}
        self.pucorrector = mcCorrections.make_puCorrector('singlee')

        
    def event_weight(self, row):
        if row.run > 2: #FIXME! add tight ID correction
            return 1.
        
        etrig = 'e1'
        if row.e2Pt > row.e1Pt : etrig = 'e2'
        if bool(row.e1MatchesSingleE27WP80) and  not bool(row.e2MatchesSingleE27WP80) : etrig = 'e1'
        if not bool(row.e1MatchesSingleE27WP80) and  bool(row.e2MatchesSingleE27WP80) :  etrig = 'e2'

        return  self.pucorrector(row.nTruePU) * \
            mcCorrections.get_electronId_corrections13_MVA(row, 'e1') * \
            mcCorrections.get_electronIso_corrections13_MVA(row, 'e1') * \
            mcCorrections.get_electronId_corrections13_MVA(row, 'e2') * \
            mcCorrections.get_electronIso_corrections13_MVA(row, 'e2')


 
## 
    def begin(self):
        sign=['os', 'ss']
        jets = [0,1,2,3]
        folder=[]
        for i in sign:
            folder.append(i)
            for j in jets:
                folder.append(i+'/'+str(j))
        

        for f in folder: 
            self.book(f, "evtInfo", "evtInfo", "run/l:lumi/l:evt/l", type=pytree.PyTree)
            
            self.book(f,"e1Pt", "e1 p_{T}", 200, 0, 200)
            self.book(f,"e1Phi", "e1 phi",  100, -3.2, 3.2)
            self.book(f,"e1Eta", "e1 eta", 50, -2.5, 2.5)
            self.book(f,"e2Pt", "e2 p_{T}", 200, 0, 200)
            self.book(f,"e2Phi", "e2 phi",  100, -3.2, 3.2)
            self.book(f,"e2Eta", "e2 eta", 50, -2.5, 2.5)
            
            self.book(f, "e1e2_DeltaPhi", "e1-e2 DeltaPhi" , 50, 0, 3.2)
            self.book(f, "e1e2_DeltaR", "e1-e2 DeltaR" , 50, 0, 3.2)
                
            self.book(f, "e1e2Mass",  "e1e2 Inv Mass",  320, 0, 320)
                
            self.book(f, "pfMetEt", "pfMetEt",  50, 0, 100)
            self.book(f, "type1_pfMetEt", "type1_pfMetEt",  50, 0, 100)
            #self.book(f, "mvaMetEt", "mvaMetEt", 50, 0, 100)
            self.book(f, "pfMetPhi", "pfMetPhi", 100, -3.2, 3.2)
            self.book(f, "type1_pfMetPhi", "type1_pfMetPhi", 100, -3.2, 3.2)
            #self.book(f, "mvaMetPhi", "mvaMetPhi", 100, -3.2, 3.2)
            self.book(f, "type1_pfMetEt_par", "type1_pfMetEt_par", 100, -100, 100)
            self.book(f, "type1_pfMetEt_perp", "type1_pfMetEt_perp", 50, 0, 100)
            self.book(f, "pfMetEt_par", "pfMetEt_par", 100, -100, 100)
            self.book(f, "pfMetEt_perp", "pfMetEt_perp", 50, 0, 100)
            #self.book(f, "mvaMetEt_par", "mvaMetEt_par", 100, -100, 100)
            #self.book(f, "mvaMetEt_perp", "mvaMetEt_perp", 50, 0, 100)
                
                
            self.book(f, "e1PFMET_DeltaPhi", "e1-type1PFMET DeltaPhi" , 50, 0, 3.2)
            self.book(f, "e1PFMET_Mt", "e1-type1PFMET M_{T}" , 200, 0, 200)
            #self.book(f, "e1MVAMET_DeltaPhi", "e1-MVAMET DeltaPhi" , 50, 0, 3.2)
            #self.book(f, "e1MVAMET_Mt", "e1-MVAMET M_{T}" , 200, 0, 200)
            
            self.book(f, "e2PFMET_DeltaPhi", "e2-type1PFMET DeltaPhi" , 50, 0, 3.2)
            self.book(f, "e2PFMET_Mt", "e2-type1PFMET M_{T}" , 200, 0, 200)
            #self.book(f, "e2MVAMET_DeltaPhi", "e2-MVAMET DeltaPhi" , 50, 0, 3.2)
            #self.book(f, "e2MVAMET_Mt", "e2-MVAMET M_{T}" , 200, 0, 200)
                
                    
    def fill_histos(self, row, folder='os', fakeRate = False):
        
        weight = self.event_weight(row)

        histos = self.histograms
        ##histos[folder+'/evtInfo'].Fill(row)

        histos[folder+'/e1Pt'].Fill(row.e1Pt, weight)
        histos[folder+'/e1Eta'].Fill(row.e1Eta, weight)
        histos[folder+'/e1Phi'].Fill(row.e1Phi, weight) 

        histos[folder+'/e2Pt'].Fill(row.e2Pt, weight)
        histos[folder+'/e2Eta'].Fill(row.e2Eta, weight)
        histos[folder+'/e2Phi'].Fill(row.e2Phi, weight)

        histos[folder+'/e1e2_DeltaPhi'].Fill(deltaPhi(row.e1Phi, row.e2Phi), weight)
        histos[folder+'/e1e2_DeltaR'].Fill(row.e1_e2_DR, weight)

        histos[folder+'/e1e2Mass'].Fill(row.e1_e2_Mass, weight)

        histos[folder+'/e1PFMET_DeltaPhi'].Fill(deltaPhi(row.e1Phi, row.type1_pfMetPhi), weight)
        #histos[folder+'/e1MVAMET_DeltaPhi'].Fill(deltaPhi(row.e1Phi, row.mva_metPhi), weight)
        histos[folder+'/e1PFMET_Mt'].Fill(row.e1MtToPFMET, weight)
        #histos[folder+'/e1MVAMET_Mt'].Fill(row.e1MtToMVAMET, weight)

        histos[folder+'/type1_pfMetEt'].Fill(row.type1_pfMetEt, weight)
        histos[folder+'/pfMetEt'].Fill(row.type1_pfMetEt, weight)
        #histos[folder+'/mvaMetEt'].Fill(row.mva_metEt, weight)
        histos[folder+'/type1_pfMetPhi'].Fill(row.type1_pfMetPhi, weight)
        histos[folder+'/pfMetPhi'].Fill(row.pfMetPhi, weight)
        #histos[folder+'/mvaMetPhi'].Fill(row.mva_metPhi, weight)

        zphi = Z(row).Phi()
        histos[folder+'/type1_pfMetEt_par'].Fill(row.pfMetEt*cos(deltaPhi(zphi, row.type1_pfMetPhi)), weight)
        histos[folder+'/type1_pfMetEt_perp'].Fill(row.pfMetEt*sin(deltaPhi(zphi, row.type1_pfMetPhi)), weight)
        histos[folder+'/pfMetEt_par'].Fill(row.pfMetEt*cos(deltaPhi(zphi, row.type1_pfMetPhi)), weight)
        histos[folder+'/pfMetEt_perp'].Fill(row.pfMetEt*sin(deltaPhi(zphi, row.type1_pfMetPhi)), weight)
        #histos[folder+'/mvaMetEt_par'].Fill(row.mva_metEt*cos(deltaPhi(zphi, row.mva_metPhi)), weight)
        #histos[folder+'/mvaMetEt_perp'].Fill(row.mva_metEt*sin(deltaPhi(zphi, row.mva_metPhi)), weight)
 
        histos[folder+'/e2PFMET_DeltaPhi'].Fill(deltaPhi(row.e2Phi, row.type1_pfMetPhi), weight)
        #histos[folder+'/e2MVAMET_DeltaPhi'].Fill(deltaPhi(row.e2Phi, row.mva_metPhi), weight)
        histos[folder+'/e2PFMET_Mt'].Fill(row.e2MtToPFMET, weight)
        #histos[folder+'/e2MVAMET_Mt'].Fill(row.e2MtToMVAMET, weight)


    def process(self):
        event = ()
        for row in self.tree:
#        for i, row in enumerate(self.tree):
#            if  i >= 100:
#                return

            
            jn = row.jetVeto30
            if jn > 3 : jn = 3
            sign = 'ss' if row.e1_e2_SS else 'os'
            if row.run > 2 :
                if not bool(row.singleE27WP80Pass) : continue
                if  not  bool(row.e1MatchesSingleE27WP80) and   not bool(row.e1MatchesSingleE27WP80) : continue
           
            if jn != 0 and row.bjetCSVVeto30!=0 : continue 
            
            if row.e1Pt < 30 : continue
            if row.e2Pt < 30 : continue
            #if not abs(row.e1_e2_Mass-91.2) < 20: continue
            
            if not selections.eSelection(row, 'e1'): continue
            if not selections.lepton_id_iso(row, 'e1', 'eid13Tight_etauiso01'): continue
            if abs(row.e1Eta) > 1.4442 and abs(row.e1Eta) < 1.566 : continue
            
            if not selections.eSelection(row, 'e2'): continue
            if not selections.lepton_id_iso(row, 'e2', 'eid13Tight_etauiso01'): continue
            if abs(row.e2Eta) > 1.4442 and abs(row.e2Eta) < 1.566 : continue
           

            #if not selections.vetos(row) : continue
            if row.muVetoPt5IsoIdVtx : continue
            if row.eVetoCicLooseIso : continue # change it with Loose
#            if row.tauVetoPt20EleTight3MuLoose : continue           
            if row.tauHpsVetoPt20 : continue
        

            folder = sign
            self.histograms[folder+'/evtInfo'].Fill(row)
            new_event = (row.run, row.lumi, row.evt) 
            if event != new_event :
                event = new_event 
                self.fill_histos(row, folder)
                folder = sign+'/'+str(int(jn)) 
                self.fill_histos(row, folder)
             
            
    def finish(self):
        self.write_histos()


