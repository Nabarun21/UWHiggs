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
import FinalStateAnalysis.TagAndProbe.EGammaPOGCorrections as EGammaPOGCorrections
import FinalStateAnalysis.TagAndProbe.FakeRate2D as FakeRate2D
import FinalStateAnalysis.TagAndProbe.HetauCorrection as HetauCorrection
import bTagSF as bTagSF
from inspect import currentframe




cut_flow_step=['allEvents','HLTIsoPasstrg','esel','eiso','musel','muiso','bjetveto','DR_e_mu','surplus_mu_veto','surplus_e_veto','jet0sel','jet1sel','jet2loosesel','jet2tightsel']

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



def topPtreweight(pt1,pt2):
    #pt1=pt of top quark
    #pt2=pt of antitop quark
    #13 Tev parameters: a=0.0615,b=-0.0005
    #for toPt >400, apply SF at 400

    if pt1>400:pt1=400
    if pt2>400:pt2=400
    a=0.0615
    b=-0.0005 

    wt1=math.exp(a+b*pt1)
    wt2=math.exp(a+b*pt2)

    wt=sqrt(wt1*wt2)

    return wt

pu_distributions = glob.glob(os.path.join('inputs', os.environ['jobid'], 'data_SingleMu*pu.root'))


pu_corrector = PileupWeight.PileupWeight('MC_Spring16', *pu_distributions)
id_corrector  = MuonPOGCorrections.make_muon_pog_PFMedium_2016BCD()
iso_corrector = MuonPOGCorrections.make_muon_pog_TightIso_2016BCD()
trg_corrector  = MuonPOGCorrections.make_muon_pog_IsoMu22oIsoTkMu22_2016BCD()
mtrk_corrector = MuonPOGCorrections.mu_trackingEta_2016
#trk_corrector =  MuonPOGCorrections.make_muonptabove10_pog_tracking_corrections_2016()
#eId_corrector = EGammaPOGCorrections.make_egamma_pog_electronID_ICHEP2016( 'nontrigWP80')
etrk_corrector=EGammaPOGCorrections.make_egamma_pog_tracking_ICHEP2016()
eiso_corr0p10 =HetauCorrection.iso0p10_ele_2016
eiso_corr0p15 =HetauCorrection.iso0p15_ele_2016



fakerateWeight =FakeRate2D.make_fakerate2D()
class LFVHEMuAnalyzerMVAmakeBDTtrees(MegaBase):
    tree = 'em/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel='EMu'
        target = os.path.basename(os.environ['megatarget'])
        self.target=target
        self.is_WJet=('WJetsToLNu' in target or 'W1JetsToLNu' in target or 'W2JetsToLNu' in target or 'W3JetsToLNu' in target or 'W4JetsToLNu' in target)
        self.is_DYJet= ('DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target or  'DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target or 'DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target or 'DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target or 'DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target) 
        self.isDYlowmass=('DYJetsToLL_M-10to50_' in target)
        self.isData=('data' in target)
        self.isWGToLNuG=( 'WGToLNuG' in target)
        self.isWGstarToLNuEE=('WGstarToLNuEE' in target)
        self.isWGstarToLNuMuMu=('WGstarToLNuMuMu' in target)
        self.isST_tW_antitop=('ST_tW_antitop' in target)
        self.isST_tW_top=('ST_tW_top' in target)
        self.isWW=('WW_Tune' in target)
        self.isWZ=('WZ_Tune' in target)
        self.isZZ=('ZZ_Tune' in target)
        self.isTT=('TT_Tune' in target)
        self.isGluGluHTo=('GluGluHTo' in target)
        self.isGluGlu_LFV=('GluGlu_LFV_HToMuTau' in target)
        self.isVBFHTo=('VBFHTo' in target)
        self.isVBF_LFV=('VBF_LFV_HToMuTau' in target)
        self.isGluGluEtauSig=('GluGlu_LFV_HToETau' in target)
        self.isVBFEtauSig=('VBF_LFV_HToETau' in target)

        self.DYlowmass_weight=1.99619334706e-08
        self.WGToLNuG_weight=1.00735804765e-07
        self.WGstarToLNuEE_weight=1.58376883982e-06#0.00000564158
        self.WGstarToLNuMuMu_weight=1.25857015075e-06
        self.ST_tW_antitop_weight=3.61421319798e-05
        self.ST_tW_top_weight= 3.56570512821e-05
        self.WW_weight= 0.000119511001657
        self.WZ_weight=4.713e-05
        self.ZZ_weight=1.67015056929e-05
        self.TT_weight=8.70585890101e-06
        self.GluGluHTo_weight=2.04805444356e-06
        self.GluGlu_LFV_HToMuTau_weight=1.9428e-06
        self.VBFHTo_weight= 4.27406682083e-08
        self.VBF_LFV_HToMuTau_weight=4.88411526098e-08  
        self.GluGlu_LFV_HToETau_weight=1.9428e-06
        self.VBF_LFV_HToETau_weight=4.05239613191e-08  


        self.tree = EMTree(tree)
        self.output=outfile
        self.histograms = {}
        self.mym1='m'
        self.mye1='e'
        self.sysdir=['nosys','jetup','jetdown','tup','tdown','uup','udown']
#        self.sysdir=['nosys']
        if self.is_WJet:
            self.binned_weight=[0.618332066,0.199958214,0.106098513,0.053599448,0.058522049]
        elif self.is_DYJet:
            self.binned_weight=[0.064154079,0.014052138,0.01505218,0.015583224,0.012508924]
        else:
            self.binned_weight=[1,1,1,1,1]




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


    def mc_corrector_2015(self, row, region):
        
        pu = pu_corrector(row.nTruePU)
 #       pu=1
        muidcorr = id_corrector(getattr(row, self.mym1+'Pt'), abs(getattr(row, self.mym1+'Eta')))
        muisocorr = iso_corrector('Tight', getattr(row, self.mym1+'Pt'), abs(getattr(row, self.mym1+'Eta')))
#        print "id corr", muidcorr
 #       print "iso corr", muisocorr
        mutrcorr = trg_corrector(getattr(row, self.mym1+'Pt'), abs(getattr(row, self.mym1+'Eta'))) 
#        eidcorr = eId_corrector(getattr(row,self.mye1+'Eta'),getattr(row, self.mye1+'Pt'))
        mutrkcorr=mtrk_corrector(getattr(row,self.mym1+'Eta'))[0]
        eisocorr0p10= eiso_corr0p10(getattr(row, self.mye1+'Pt'),abs(getattr(row,self.mye1+'Eta')))[0]
        eisocorr0p15= eiso_corr0p15(getattr(row, self.mye1+'Pt'),abs(getattr(row,self.mye1+'Eta')))[0]
        etrkcorr=etrk_corrector(getattr(row,self.mye1+'Eta'),getattr(row, self.mye1+'Pt'))
#        mutrkcorr=trk_corrector(getattr(row,self.mym1+'Eta'))
 #       print "trk corr",mutrkcorr
  #      print "tr corr", mutrcorr
   #     print "eid corr", eidcorr
#       mutrcorr=1
     # if pu*muidcorr1*muisocorr1*muidcorr2*muisocorr2*mutrcorr==0: print pu, muidcorr1, muisocorr1, muidcorr2, muisocorr2, mutrcorr
    #    print "pileup--------   =",pu
   #     print pu*muidcorr*muisocorr*mutrcorr
#        print eisocorr

        topptreweight=1
        eidcorr=1
        if self.isTT:
            topptreweight=topPtreweight(row.topQuarkPt1,row.topQuarkPt2)

        return pu*muidcorr*muisocorr*mutrcorr*eidcorr*mutrkcorr*eisocorr0p10*etrkcorr*topptreweight
 
       # return pu*muidcorr*mutrcorr*eidcorr


    def correction(self,row,region):
	return self.mc_corrector_2015(row,region)
        
    def event_weight(self, row, region):
 
        if row.run > 2: #FIXME! add tight ID correction
            return 1.
       # if row.GenWeight*self.correction(row) == 0 : print 'weight==0', row.GenWeight*self.correction(row), row.GenWeight, self.correction(row), row.m1Pt, row.m2Pt, row.m1Eta, row.m2Eta
       # print row.GenWeight, "lkdfh"


        return row.GenWeight*self.correction(row,region) 
#        return self.correction(row) 


    def begin(self):
        self.weight_=array.array( 'f', [ 0 ] )
        self.mPt_=array.array( 'f', [ 0 ] )
        self.ePt_=array.array( 'f', [ 0 ] )
        self.deltaPhimue_=array.array( 'f', [ 0 ] )
        self.mEta_=array.array( 'f', [ 0 ] )
        self.eEta_=array.array( 'f', [ 0 ] )
        self.mTmuMet_=array.array( 'f', [ 0 ] )
        self.mTeMet_=array.array( 'f', [ 0 ] )
        self.eDphiPFMet_=array.array( 'f', [ 0 ] )
        self.mDphiPFMet_=array.array( 'f', [ 0 ] )
        self.metEt_=array.array( 'f', [ 0 ] )
        self.metPhi_=array.array( 'f', [ 0 ] )
        self.vbfMass_=array.array( 'f', [ 0 ] )
        self.vbfDeltaEta_=array.array( 'f', [ 0 ] )
        self.numjets_=array.array( 'f', [ 0 ] )
        self.mColl_=array.array( 'f', [ 0 ] )
        self.pZeta_=array.array( 'f', [ 0 ] )
        self.lepAsym_=array.array( 'f', [ 0 ] )        

        if self.isGluGlu_LFV or self.isVBF_LFV or self.isVBFEtauSig or self.isGluGluEtauSig:
            self.treeS=ROOT.TTree("treeS","treeS")
            self.treeS.Branch("weight_",self.weight_,"weight_/F")
            self.treeS.Branch("mPt_",self.mPt_,"mPt_/F")
            self.treeS.Branch("ePt_",self.ePt_,"ePt_/F")            
            self.treeS.Branch("deltaPhimue_",self.deltaPhimue_,"deltaPhimue_/F")
            self.treeS.Branch("mEta_",self.mEta_,"mEta_/F")
            self.treeS.Branch("eEta_",self.eEta_,"eEta_/F")
            self.treeS.Branch("mTmuMet_",self.mTmuMet_,"mTmuMet_/F")
            self.treeS.Branch("mTeMet_",self.mTeMet_,"mTeMet_/F")
            self.treeS.Branch("eDphiPFMet_",self.eDphiPFMet_,"eDphiPFMet_/F")
            self.treeS.Branch("mDphiPFMet_",self.mDphiPFMet_,"mDphiPFMet_/F")
            self.treeS.Branch("metEt_",self.metEt_,"metEt_/F")
            self.treeS.Branch("metPhi_",self.metPhi_,"metPhi_/F")
            self.treeS.Branch("vbfMass_",self.vbfMass_,"vbfMass_/F")
            self.treeS.Branch("vbfDeltaEta_",self.vbfDeltaEta_,"vbfDeltaEta_/F")
            self.treeS.Branch("numjets_",self.numjets_,"numjets_/F")
            self.treeS.Branch("mColl_",self.mColl_,"mColl_/F")
            self.treeS.Branch("pZeta_",self.pZeta_,"pZeta_/F")
            self.treeS.Branch("lepAsym_",self.lepAsym_,"lepAsym_/F")
        else:
            self.treeB=ROOT.TTree("treeB","treeB")
            self.treeB.Branch("weight_",self.weight_,"weight_/F")
            self.treeB.Branch("mPt_",self.mPt_,"mPt_/F")
            self.treeB.Branch("ePt_",self.ePt_,"ePt_/F")            
            self.treeB.Branch("deltaPhimue_",self.deltaPhimue_,"deltaPhimue_/F")
            self.treeB.Branch("mEta_",self.mEta_,"mEta_/F")
            self.treeB.Branch("eEta_",self.eEta_,"eEta_/F")
            self.treeB.Branch("mTmuMet_",self.mTmuMet_,"mTmuMet_/F")
            self.treeB.Branch("mTeMet_",self.mTeMet_,"mTeMet_/F")
            self.treeB.Branch("eDphiPFMet_",self.eDphiPFMet_,"eDphiPFMet_/F")
            self.treeB.Branch("mDphiPFMet_",self.mDphiPFMet_,"mDphiPFMet_/F")
            self.treeB.Branch("metEt_",self.metEt_,"metEt_/F")
            self.treeB.Branch("metPhi_",self.metPhi_,"metPhi_/F")
            self.treeB.Branch("vbfMass_",self.vbfMass_,"vbfMass_/F")
            self.treeB.Branch("vbfDeltaEta_",self.vbfDeltaEta_,"vbfDeltaEta_/F")
            self.treeB.Branch("numjets_",self.numjets_,"numjets_/F")
            self.treeB.Branch("mColl_",self.mColl_,"mColl_/F")
            self.treeB.Branch("pZeta_",self.pZeta_,"pZeta_/F")
            self.treeB.Branch("lepAsym_",self.lepAsym_,"lepAsym_/F")


        


    def fill_tree(self, row,btagweight=1,region='signal'):

        if self.is_WJet or self.is_DYJet:
            weight = self.event_weight(row,region) *self.binned_weight[int(row.numGenJets)]*0.001
        elif self.isWGToLNuG:
            weight=self.WGToLNuG_weight*self.event_weight(row,region) 
        elif self.isDYlowmass:
            weight=self.DYlowmass_weight*self.event_weight(row,region) 
        elif self.isWGstarToLNuEE:
            weight=self.WGstarToLNuEE_weight*self.event_weight(row,region) 
        elif self.isWGstarToLNuMuMu:
            weight=self.WGstarToLNuMuMu_weight*self.event_weight(row,region) 
        elif self.isST_tW_top:
            weight=self.ST_tW_top_weight*self.event_weight(row,region) 
        elif self.isST_tW_antitop:
            weight=self.ST_tW_antitop_weight*self.event_weight(row,region) 
        elif self.isWW:
            weight=self.WW_weight*self.event_weight(row,region) 
        elif self.isWZ:
            weight=self.WZ_weight*self.event_weight(row,region) 
        elif self.isZZ:
            weight=self.ZZ_weight*self.event_weight(row,region) 
        elif self.isTT:
            weight=self.TT_weight*self.event_weight(row,region) 
        elif self.isGluGluHTo:
            weight=self.GluGluHTo_weight*self.event_weight(row,region) 
        elif self.isGluGlu_LFV:
            weight=self.GluGlu_LFV_HToMuTau_weight*self.event_weight(row,region) 
        elif self.isVBFHTo:
            weight=self.VBFHTo_weight*self.event_weight(row,region) 
        elif self.isVBF_LFV:
            weight=self.VBF_LFV_HToMuTau_weight*self.event_weight(row,region) 
        elif self.isVBFEtauSig:
            weight=self.VBF_LFV_HToETau_weight*self.event_weight(row,region) 
        elif self.isGluGluEtauSig:
            weight=self.GluGlu_LFV_HToETau_weight*self.event_weight(row,region) 
        else:
            weight = self.event_weight(row,region) 
        
#        if btagweight<0:
 #           print "btagweight is negative:  ",btagweight
        weight=btagweight*weight

        self.weight_[0]=weight
        self.mPt_[0]=row.mPt
        self.ePt_[0]=row.ePt
        self.deltaPhimue_[0]=abs(row.e_m_DPhi)
        self.mEta_[0]=row.mEta
        self.eEta_[0]=row.eEta
        self.mTmuMet_[0]=row.mMtToPfMet_type1  
        self.mTeMet_[0]=row.eMtToPfMet_type1 
        self.eDphiPFMet_[0]=abs(row.eDPhiToPfMet_type1)
        self.mDphiPFMet_[0]=abs(row.mDPhiToPfMet_type1)
        self.metEt_[0]=row.type1_pfMetEt  
        self.metPhi_[0]=row.type1_pfMetPhi
        self.numjets_[0]=row.jetVeto30
        self.mColl_[0]=collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)
        self.pZeta_[0]=row.e_m_PZetaLess0p85PZetaVis
        self.lepAsym_[0]=(row.mPt-row.ePt)/(row.ePt+row.mPt)

        if row.jetVeto30==0 or row.jetVeto30==1:
            self.vbfMass_[0]=0.01
            self.vbfDeltaEta_[0]=0.01
        else:
            self.vbfMass_[0]=row.vbfMass
            self.vbfDeltaEta_[0]=row.vbfDeta
        
        if self.isGluGlu_LFV or self.isVBF_LFV or self.isVBFEtauSig or self.isGluGluEtauSig:
            self.treeS.Fill()
        else:
            self.treeB.Fill()

    def process(self):
        myevent=()
        frw = []
        curr_event=0
        for row in self.tree:
            sign = 'ss' if row.e_m_SS else 'os'

            if sign=='ss':continue
 #           ptthreshold = [30]
            repeatEvt=True
            if row.evt!=curr_event:
                curr_event=row.evt
                repeatEvt=False
            
            if repeatEvt:continue

#            print "non-repeat"
            processtype ='gg'##changed from 20

           


            #trigger
            if (not bool(row.singleIsoMu22Pass or row.singleIsoTkMu22Pass)) and self.isData: 
                continue   #notrigger in new MC; add later



            #vetoes and cleaning

            if deltaR(row.ePhi,row.mPhi,row.eEta,row.mEta)<0.3:continue

            
            if row.muVetoPt5IsoIdVtx :continue


            if row.eVetoMVAIsoVtx :continue


            if row.tauVetoPt20Loose3HitsVtx : continue
            

            #mu preselection
            if not selections.muSelection(row, 'm'): continue



            #E Preselection
            if not selections.eSelection(row, 'e'): continue
           
            if not selections.lepton_id_iso(row, 'e', 'eid15Loose_etauiso1','WP80'): continue
            

            #take care of ecal gap
            if row.eAbsEta > 1.4442 and row.eAbsEta < 1.566 : continue             


            nbtagged=row.bjetCISVVeto30Medium
            if nbtagged>2:
                nbtagged=2
            btagweight=1
            if (self.isData and nbtagged>0):
                continue
            if nbtagged>0:
                if nbtagged==1:
                    btagweight=bTagSF.bTagEventWeight(nbtagged,row.jb1pt,row.jb1hadronflavor,row.jb2pt,row.jb2hadronflavor,1,0,0) if (row.jb1pt>-990 and row.jb1hadronflavor>-990) else 0
                if nbtagged==2:
                    btagweight=bTagSF.bTagEventWeight(nbtagged,row.jb1pt,row.jb1hadronflavor,row.jb2pt,row.jb2hadronflavor,1,0,0) if (row.jb1pt>-990 and row.jb1hadronflavor>-990 and row.jb2pt>-990 and row.jb2hadronflavor>-990) else 0
#                print "btagweight,nbtagged,row.jb1pt,row.jb1hadronflavor,row.jb2pt,row.jb2hadronflavor"," ",btagweight," ",nbtagged," ",row.jb1pt," ",row.jb1hadronflavor," ",row.jb2pt," ",row.jb2hadronflavor

            if btagweight<0:btagweight=0

            if btagweight==0: continue

            if not selections.lepton_id_iso(row, 'm', 'MuIDTight_mutauiso0p15'):continue
 
            if not selections.lepton_id_iso(row, 'e', 'eid15Loose_etauiso0p1','WP80'):continue 
 
            self.fill_tree(row,btagweight)

            
    def finish(self):
        self.output.Write()


