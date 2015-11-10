##Correction Factor still to add
from MMTTree import MMTTree
import os
import ROOT
import math
import optimizer
import glob
import array
import mcCorrections
import baseSelections as selections
import FinalStateAnalysis.PlotTools.pytree as pytree
from FinalStateAnalysis.PlotTools.decorators import  memo_last
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
from math import sqrt, pi, sin, cos, acos, sinh
from cutflowtracker import cut_flow_tracker
#Makes the cut flow histogram
from inspect import currentframe

def get_linenumber():
    cf = currentframe()
    return cf.f_back.f_lineno
def deltaPhi(phi1, phi2):
    PHI = abs(phi1-phi2)
    if PHI<=pi:
        return PHI
    else:
        return 2*pi-PHI
def deltaR(phi1, phi2, eta1, eta2):
    deta = eta1 - eta2
    dphi = abs(phi1-phi2)
    if (dphi>pi) : dphi = 2*pi-dphi
    return sqrt(deta*deta + dphi*dphi);
    
def mtDR(row):
    return row.m1_t_DR if row.m1_t_DR < row.m2_t_DR else row.m2_t_DR

def mtDPhi(row):
    m1tDPhi=deltaPhi(row.m1Phi, row.tPhi)
    m2tDPhi=deltaPhi(row.m2Phi, row.tPhi)
    return m1tDPhi if m1tDPhi < m2tDPhi else m2tDPhi

def Z(row):
    m1p=ROOT.TVector3(row.m1Pt*cos(row.m1Phi),row.m1Pt*sin(row.m1Phi),row.m1Pt*sinh(row.m1Eta))
    m2p=ROOT.TVector3(row.m2Pt*cos(row.m2Phi),row.m2Pt*sin(row.m2Phi),row.m2Pt*sinh(row.m2Eta))
    m1FourVector= ROOT.TLorentzVector(m1p, sqrt(m1p.Mag2()+row.m1Mass*row.m1Mass))
    m2FourVector= ROOT.TLorentzVector(m2p, sqrt(m2p.Mag2()+row.m2Mass*row.m2Mass))
    zFourVector = m1FourVector+m2FourVector
    return zFourVector


class TauFakeRateAnalyzerMVA(MegaBase):
    tree = 'mmt/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel='MMT'
        super(TauFakeRateAnalyzerMVA, self).__init__(tree, outfile, **kwargs)
        self.tree = MMTTree(tree)
        self.out=outfile
        self.histograms = {}
        #self.pucorrector = mcCorrections.make_puCorrector('singlem')


        if  ('RUN_OPTIMIZATION' in os.environ) and eval(os.environ['RUN_OPTIMIZATION']):
            optimizer_keys   = [ i for i in optimizer.grid_search.keys() if i.startswith(self.channel) ]
            self.grid_search = {}
            if len(optimizer_keys) > 1:
                for key in optimizer_keys:
                    self.grid_search[key] = optimizer.grid_search[key]
            else:
                self.grid_search[''] = optimizer.grid_search[optimizer_keys[0]]


    def event_weight(self, row):
        if row.run > 2: #FIXME! add tight ID correction
            return 1.

        etrig = 'm1'
        if row.m2Pt > row.m1Pt : etrig = 'm2'
        #if bool(row.m1MatchesSingleE27WP80) and  not bool(row.e2MatchesSingleE27WP80) : etrig = 'e1'
        #if not bool(row.m1MatchesSingleE27WP80) and  bool(row.e2MatchesSingleE27WP80) :  etrig = 'e2'

            
        return 1.    
#        return self.pucorrector(row.nTruePU) * \
#            mcCorrections.eid_correction( row, 'e1', 'e2') * \
#            mcCorrections.eiso_correction(row, 'e1', 'e2') * \
#            mcCorrections.trig_correction(row, etrig     )

           
    def begin(self):
        
        tauiso = ['tNoCuts', 'tSuperSuperLoose', 'tSuperLoose', 'tLoose', 'tTigh']
        folder = []
        sign = ['ss','os']
        for iso in tauiso:
            for s in sign:
                folder.append(s+'/'+iso)
                j=0
                while j < 4 :
                    folder.append(s+'/'+iso+'/'+str(j))
                    j+=1
                    
        for f in folder: 
            
            self.book(f,"m1Pt", "m1 p_{T}", 200, 0, 200)
            self.book(f,"m1Phi", "m1 phi",  100, -3.2, 3.2)
            self.book(f,"m1Eta", "m1 eta", 50, -2.5, 2.5)

            self.book(f,"m2Pt", "m2 p_{T}", 200, 0, 200)
            self.book(f,"m2Phi", "m2 phi",  100, -3.2, 3.2)
            self.book(f,"m2Eta", "m2 eta", 50, -2.5, 2.5)

            self.book(f, "m1m2Mass",  "m1m2 Inv Mass",  32, 0, 320)

            self.book(f, "tMtToPFMET", "#tau Met MT", 100, 0, 100)
            
            self.book(f,"tByPileupWeightedIsolationRaw3Hits", "tByPileupWeightedIsolationRaw3Hits", 500, 0, 100) 
            self.book(f,"tPt", "t p_{T}", 200, 0, 200)
            self.book(f,"tPtbarrel", "t p_{T} barrel", 200, 0, 200)
            self.book(f,"tPtendcap", "t p_{T} endcap", 200, 0, 200)
            self.book(f,"tPhi", "t phi",  100, -3.2, 3.2)
            self.book(f,"tEta", "t eta", 50, -2.5, 2.5)
            self.book(f,"tAbsEta", "t abs eta", 50, -2.5, 2.5)
 
            self.book(f,"mtDR", "m t DR", 50, 0, 10)
            self.book(f,"mtDPhi", "m t DPhi", 32, 0, 3.2)

            self.book(f,"ztDR", "Z #tau DR", 50, 0, 10)
            self.book(f,"ztDPhi", "Z #tau DPhi", 32, 0, 3.2)
            self.book(f,"Zpt", "Z p_{T}", 200, 0, 200)

            

            self.book(f, "type1_pfMetEt", "type1_pfMetEt",200, 0, 200) 
            self.book(f, "jetN_30", "Number of jets, p_{T}>30", 10, -0.5, 9.5)
            self.book(f, "bjetCSVVeto30", "number of bjets", 10, -0.5, 9.5)

                      
    def fill_histos(self, row, folder='os/tSuperLoose', fakeRate = False):
        weight = self.event_weight(row)
        histos = self.histograms
 
        histos[folder+'/m1Pt'].Fill(row.m1Pt, weight)
        histos[folder+'/m1Eta'].Fill(row.m1Eta, weight)
        histos[folder+'/m1Phi'].Fill(row.m1Phi, weight) 

        histos[folder+'/m2Pt'].Fill(row.m2Pt, weight)
        histos[folder+'/m2Eta'].Fill(row.m2Eta, weight)
        histos[folder+'/m2Phi'].Fill(row.m2Phi, weight)

        histos[folder+'/m1m2Mass'].Fill(row.m1_m2_Mass, weight)
        histos[folder+'/tMtToPFMET'].Fill(row.tMtToPfMet_type1,weight)
    
        histos[folder+'/tByPileupWeightedIsolationRaw3Hits'].Fill(row.tByPileupWeightedIsolationRaw3Hits, weight)
        histos[folder+'/tPt'].Fill(row.tPt, weight)
        if abs(row.tEta) < 1.5 :  histos[folder+'/tPtbarrel'].Fill(row.tPt, weight)
        if abs(row.tEta) > 1.5 :  histos[folder+'/tPtendcap'].Fill(row.tPt, weight)
        histos[folder+'/tEta'].Fill(row.tEta, weight)
        histos[folder+'/tAbsEta'].Fill(abs(row.tEta), weight)
        histos[folder+'/tPhi'].Fill(row.tPhi, weight) 
 
        histos[folder+'/type1_pfMetEt'].Fill(row.type1_pfMetEt)
        histos[folder+'/mtDR'].Fill(mtDR(row)) 
        histos[folder+'/mtDPhi'].Fill(mtDPhi(row)) 
        histos[folder+'/jetN_30'].Fill(row.jetVeto30, weight) 
        histos[folder+'/bjetCSVVeto30'].Fill(row.bjetCISVVeto30Medium, weight) 
       
        histos[folder+'/ztDR'].Fill(deltaR(Z(row).Phi(), row.tPhi, Z(row).Eta(), row.tEta))
        histos[folder+'/ztDPhi'].Fill(deltaPhi(Z(row).Phi(), row.tPhi))
        histos[folder+'/Zpt'].Fill(Z(row).Pt())
            

    def process(self):


        myevent =()
        for row in self.tree:
            jn = row.jetVeto30
            if jn > 3 : jn = 3
            
            if not bool(row.singleIsoMu20Pass) : continue

            if row.bjetCISVVeto30Medium!=0 : continue 
            if row.m1Pt < 30 : continue
            if row.m2Pt < 30 : continue
            
#        for i, row in enumerate(self.tree):
#            if  i >= 100:
#                return
 
            if not selections.muSelection(row, 'm1'): continue
            #print row.m1Pt, row.m2Pt, row.tPt, row.m1RelPFIsoDBDefault
            if not selections.lepton_id_iso(row, 'm1', 'muId_etauiso01'): continue
            if not selections.muSelection(row, 'm2'): continue
            if not selections.lepton_id_iso(row, 'm2', 'muId_etauiso01'): continue
            
            if not abs(row.m1_m2_Mass-91.2) < 20: continue
            if not selections.tauSelection(row, 't'): continue
            if row.tPt < 30 : continue

            if not row.tAgainstMuonTight3: continue
            if not row.tAgainstElectronTightMVA5: continue #was 3

            if row.tauVetoPt20Loose3HitsNewDMVtx : continue 
            if row.muVetoPt5IsoIdVtx : continue
            if row.eVetoMVAIso: continue # change it with Loose

            
            #            if  etDR(row) < 1. : continue 
            if (row.run, row.lumi, row.evt, row.m1Pt, row.m2Pt)==myevent: continue
            myevent=(row.run, row.lumi, row.evt, row.m1Pt, row.m2Pt)

            tauiso = 'tNoCuts'
            sign = 'ss' if row.m1_m2_SS else 'os'
            folder = sign+'/'+tauiso
          
            self.fill_histos(row, folder)
            folder=folder+'/'+str(int(jn))
            self.fill_histos(row, folder)
            
            if not row.tByPileupWeightedIsolationRaw3Hits< 10 : continue
            tauiso = 'tSuperSuperLoose'
            folder = sign+'/'+tauiso
            self.fill_histos(row, folder)
            folder=folder+'/'+str(int(jn))
            self.fill_histos(row, folder)                
            if not row.tByPileupWeightedIsolationRaw3Hits < 5 : continue
            tauiso = 'tSuperLoose'
            folder = sign+'/'+tauiso
            self.fill_histos(row, folder)
            folder=folder+'/'+str(int(jn))
            self.fill_histos(row, folder)

            if  row.tByLooseIsolationMVA3newDMwLT : 
                tauiso = 'tLoose'
                folder = sign+'/'+tauiso
                self.fill_histos(row,  folder)
                
                folder=folder+'/'+str(int(jn))
                self.fill_histos(row, folder)
               
            if row.tByTightIsolationMVA3newDMwLT :
                tauiso = 'tTigh' 
                folder = sign+'/'+tauiso
                self.fill_histos(row,  folder)

                folder=folder+'/'+str(int(jn))
                self.fill_histos(row, folder)
              
    
    def finish(self):
        self.write_histos()


