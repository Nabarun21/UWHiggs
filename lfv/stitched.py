# check in https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2015#MET when the mva met receipe is available.
from EMTree import EMTree
import os
import ROOT
import math
import glob
import array
#import mcCorrections
import baseSelections as selections
import FinalStateAnalysis.PlotTools.pytree as pytree
from FinalStateAnalysis.PlotTools.decorators import  memo_last
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
from cutflowtracker import cut_flow_tracker
from math import sqrt, pi, cos
#from fakerate_functions import fakerate_central_histogram, fakerate_p1s_histogram, fakerate_m1s_histogram
import FinalStateAnalysis.TagAndProbe.PileupWeight as PileupWeight
import FinalStateAnalysis.TagAndProbe.MuonPOGCorrections as MuonPOGCorrections
from inspect import currentframe




cut_flow_step=['allEvents','HLTIsoPasstrg','esel','eiso','musel','muiso','bjetveto','DR_e_mu','surplus_mu_veto','jet0sel','jet1sel','jet2sel']

def collmass(row, met, metPhi):
    ptnu =abs(met*cos(deltaPhi(metPhi, row.ePhi)))
    visfrac = row.ePt/(row.ePt+ptnu)
    #print met, cos(deltaPhi(metPhi, row.tPhi)), ptnu, visfrac
    return (row.e_m_Mass / sqrt(visfrac))

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

pu_distributions = glob.glob(os.path.join('inputs', os.environ['jobid'], 'data_SingleMu*pu.root'))
pu_corrector = PileupWeight.PileupWeight('MC_Spring16', *pu_distributions)
id_corrector  = MuonPOGCorrections.make_muon_pog_PFTight_2016B()
iso_corrector = MuonPOGCorrections.make_muon_pog_TightIso_2016B()
tr_corrector  = MuonPOGCorrections.make_muon_pog_IsoMu20oIsoTkMu20_2015()


class stitched(MegaBase):
    tree = 'em/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel='EMu'
        super(stitched, self).__init__(tree, outfile, **kwargs)
        target = os.path.basename(os.environ['megatarget'])
        self.target=target
        self.is_WJet=('WJetsToLNu' in target or 'W1JetsToLNu' in target or 'W2JetsToLNu' in target or 'W3JetsToLNu' in target or 'W4JetsToLNu' in target)
        self.is_DYJet= ('DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target or  'DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target or 'DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target or 'DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target or 'DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target) 
        self.tree = EMTree(tree)
        self.out=outfile
        self.histograms = {}
        self.mym1='m'
        self.mye1='e'
        if self.is_WJet:
            self.binned_weight=[0.003079413,0.00035568,0.000181734,0.000123275,0.000293479]
        elif self.is_DYJet:
            self.binned_weight=[0.000280512026,0.000027968676,0.000026799986,0.000033257606,0.000181650478]
        else:
            self.binned_weight=[1,1,1,1,1]

#        print self.binned_weight
        """need to think about this"""
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

    def mc_corrector_2015(self, row):
        
        pu = pu_corrector(row.nTruePU)
 #       pu=1
        muidcorr = id_corrector(getattr(row, self.mym1+'Pt'), abs(getattr(row, self.mym1+'Eta')))
        muisocorr = iso_corrector('Tight', getattr(row, self.mym1+'Pt'), abs(getattr(row, self.mym1+'Eta')))
#        print "id corr", muidcorr
 #       print "iso corr", muisocorr
 #       mutrcorr = tr_corrector(getattr(row, self.mym1+'Pt'), abs(getattr(row, self.mym1+'Eta'))) 
        mutrcorr=1
     # if pu*muidcorr1*muisocorr1*muidcorr2*muisocorr2*mutrcorr==0: print pu, muidcorr1, muisocorr1, muidcorr2, muisocorr2, mutrcorr
#        print "pileup--------   =",pu
        return pu*muidcorr*muisocorr*mutrcorr
      

    def correction(self,row):
	return self.mc_corrector_2015(row)
        
    def event_weight(self, row):
 
        if row.run > 2: #FIXME! add tight ID correction
            return 1.
       # if row.GenWeight*self.correction(row) == 0 : print 'weight==0', row.GenWeight*self.correction(row), row.GenWeight, self.correction(row), row.m1Pt, row.m2Pt, row.m1Eta, row.m2Eta
       # print row.GenWeight, "lkdfh"
        return row.GenWeight*self.correction(row) 




    def begin(self):
        sign=[ 'ss','os']

        for s in sign:
            self.book(s, "jetN_30", "Number of jets, p_{T}>30", 10, -0.5, 9.5) 
            self.book(s, "NUP", "Number of Partons", 12, -0.5, 11.5) 
            self.book(s, "numGenJets", "Number of Gen Level Jets", 12, -0.5, 11.5) 
            self.book(s+'/tNoCuts', "CUT_FLOW", "Cut Flow", len(cut_flow_step), 0, len(cut_flow_step))
            self.book(s+'all', "jetN_30", "Number of jets, p_{T}>30", 10, -0.5, 9.5) 
            self.book(s+'all', "NUP", "Number of Partons", 12, -0.5, 11.5) 
            self.book(s+'all', "numGenJets", "Number of Gen Level Jets", 12, -0.5, 11.5) 

            
            xaxis = self.histograms[s+'/tNoCuts/CUT_FLOW'].GetXaxis()
            self.cut_flow_histo = self.histograms[s+'/tNoCuts/CUT_FLOW']
            self.cut_flow_map   = {}
            for i, name in enumerate(cut_flow_step):
                xaxis.SetBinLabel(i+1, name)
                self.cut_flow_map[name] = i+0.5


    def fill_jet_histos(self,row,sign):

        if self.is_WJet or self.is_DYJet:
            weight = self.event_weight(row)*self.binned_weight[int(row.numGenJets)]
        else:
            weight = self.event_weight(row)

        self.histograms[sign+'/jetN_30'].Fill(row.jetVeto30,weight)
        self.histograms[sign+'/NUP'].Fill(row.NUP,weight)
        self.histograms[sign+'/numGenJets'].Fill(row.numGenJets,weight)

    def fill_jet_histos_all(self,row,sign):

        if self.is_WJet or self.is_DYJet:
            weight = self.event_weight(row)*self.binned_weight[int(row.numGenJets)]
        else:
            weight = self.event_weight(row)

        self.histograms[sign+'all/jetN_30'].Fill(row.jetVeto30,weight)
        self.histograms[sign+'all/NUP'].Fill(row.NUP,weight)
        self.histograms[sign+'all/numGenJets'].Fill(row.numGenJets,weight)



        

    def process(self):
        
        cut_flow_histo = self.cut_flow_histo
        cut_flow_trk   = cut_flow_tracker(cut_flow_histo)
        myevent=()

        frw = []
        for row in self.tree:
            
            sign = 'ss' if row.e_m_SS else 'os'
            
            self.fill_jet_histos_all(row,sign)
            processtype ='gg'##changed from 20

            cut_flow_trk.new_row(row.run,row.lumi,row.evt)
           
            cut_flow_trk.Fill('allEvents')

            cut_flow_trk.Fill('HLTIsoPasstrg')


            jn = row.jetVeto30


            #take care of ecal gap
            if row.eEta > 1.4442 and row.eEta < 1.566 : continue             


            if not selections.muSelection(row, 'm'): continue
            cut_flow_trk.Fill('musel')
            if not selections.lepton_id_iso(row, 'm', 'MuIDTight_idiso025'): continue
            cut_flow_trk.Fill('muiso')

            frw=1. ## add the correct fakerate weight once we have it


            isMuonTight=False
            if selections.lepton_id_iso(row, 'm', 'MuIDTight_mutauiso015'): isMuonTight=True
            #event cleaning
            


            if row.bjetCISVVeto30Loose : continue
            cut_flow_trk.Fill('bjetveto')

            if deltaR(row.ePhi,row.mPhi,row.eEta,row.mEta)<0.1:continue
            cut_flow_trk.Fill('DR_e_mu')
            
            if row.muVetoPt5IsoIdVtx :continue
            cut_flow_trk.Fill('surplus_mu_veto')

            if row.eVetoMVAIso :continue
            cut_flow_trk.Fill('surplus_e_veto')


            if not selections.eSelection(row, 'e'): continue
            cut_flow_trk.Fill('esel')
           
            if not selections.lepton_id_iso(row, 'e', 'eid15Loose_etauiso01'): continue
            cut_flow_trk.Fill('eiso')

            #need to add Surpluse Veto
            if self.is_WJet or self.is_DYJet:
                self.fill_jet_histos(row,sign)
        cut_flow_trk.flush()        
            
    def finish(self):
        self.write_histos()

