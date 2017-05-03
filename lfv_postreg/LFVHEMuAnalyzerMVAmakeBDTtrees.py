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



pu_corrector = PileupWeight.PileupWeight('MC_Moriond17', *pu_distributions)
mid_corrector  = MuonPOGCorrections.make_muon_pog_PFMedium_2016ReReco()
miso_corrector = MuonPOGCorrections.make_muon_pog_TightIso_2016ReReco("Medium")
trg_corrector  = MuonPOGCorrections.make_muon_pog_IsoMu24oIsoTkMu24_2016ReReco()
mtrk_corrector = MuonPOGCorrections.mu_trackingEta_MORIOND2017
#trk_corrector =  MuonPOGCorrections.make_muonptabove10_pog_tracking_corrections_2016()
eId_corrector = EGammaPOGCorrections.make_egamma_pog_electronID_MORIOND2017( 'nontrigWP80')
erecon_corrector=EGammaPOGCorrections.make_egamma_pog_recon_MORIOND17()
etrk_corrector=EGammaPOGCorrections.make_egamma_pog_tracking_ICHEP2016()



class LFVHEMuAnalyzerMVAmakeBDTtrees(MegaBase):
    tree = 'em/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel='EMu'
        target = os.path.basename(os.environ['megatarget'])
        self.target=target

        self.is_WJet=('WJetsToLNu' in target or 'W1JetsToLNu' in target or 'W2JetsToLNu' in target or 'W3JetsToLNu' in target or 'W4JetsToLNu' in target)
        self.is_DYJet= ('DYJetsToLL_M-50' in target or  'DY1JetsToLL_M-50' in target or 'DY2JetsToLL_M-50' in target or 'DY3JetsToLL_M-50' in target or 'DY4JetsToLL_M-50' in target) 

        self.is_DYlowmass= ('DYJetsToLL_M-10to50' in target or  'DY1JetsToLL_M-10to50' in target or 'DY2JetsToLL_M-10to50' in target or 'DY3JetsToLL_M-10to50' in target or 'DY4JetsToLL_M-10to50' in target) 

        self.is_ZTauTau= ('ZTauTauJets_M-50' in target or  'ZTauTau1Jets_M-50' in target or 'ZTauTau2Jets_M-50' in target or 'ZTauTau3Jets_M-50' in target or 'ZTauTau4Jets_M-50' in target) 
        
        self.data_period="BCDEF" if ("Run2016B" in target or "Run2016C" in target or  "Run2016D" in target or  "Run2016E" in target or  "Run2016F" in target) else "GH"

        self.isData=('data' in target)

        self.generator=ROOT.TRandom3(0)
        #set systematics flag to true if you want shape syustematic histos
        self.syscalc=True
        
        self.isWGToLNuG=( 'WGToLNuG' in target)
        self.isWGstarToLNuEE=('WGstarToLNuEE' in target)
        self.isWGstarToLNuMuMu=('WGstarToLNuMuMu' in target)

        self.isST_tW_antitop=('ST_tW_antitop' in target)
        self.isST_tW_top=('ST_tW_top' in target)
        self.isST_t_antitop=('ST_t-channel_antitop' in target)
        self.isST_t_top=('ST_t-channel_top' in target)
        self.isTT=('TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_v6-v1' in target)
        self.isTTevtgen=('TT_TuneCUETP8M2T4_13TeV-powheg-pythia8-evtgen_v6-v1' in target)
        self.isWW=('WW_Tune' in target)
        self.isWZ=('WZ_Tune' in target)
        self.isZZ=('ZZ_Tune' in target)



        self.isGluGluHTo=('GluGluHTo' in target)
        self.isGluGluHToWW=('GluGluHToWW' in target)
        self.isGluGlu_LFV=('GluGlu_LFV_HToMuTau' in target)
        self.isGluGluEtauSig=('GluGlu_LFV_HToETau' in target)


        self.isVBFHTo=('VBFHTo' in target)
        self.isVBFHToWW=('VBFHToWW' in target)
        self.isVBF_LFV=('VBF_LFV_HToMuTau' in target)
        self.isVBFEtauSig=('VBF_LFV_HToETau' in target)

        self.isZHToTauTau=('ZHToTauTau' in target)
        self.isWplusHToTauTau=('WplusHToTauTau' in target)
        self.isWminusHToTauTau=('WminusHToTauTau' in target)
        

        self.isWZTo2L2Q=('WZTo2L2Q' in target)
        self.isVVTo2L2Nu=('VVTo2L2Nu' in target)
        self.isWWTo1L1Nu2Q=('WWTo1L1Nu2Q' in target)
        self.isWZJToLLLNu=('WZJToLLLNu' in target)
        self.isWZTo1L1Nu2Q=('WZTo1L1Nu2Q' in target)
        self.isWZTo1L3Nu=('WZTo1L3Nu' in target)
        self.isZZTo2L2Q=('ZZTo2L2Q' in target)
        self.isZZTo4L=('ZZTo4L' in target)

        self.WZTo2L2Q_weight=2.54725706507e-08 #2.66175641194e-08   #2.39150954356e-08 #2.40930214376e-08  #8.91328638844e-07
        self.VVTo2L2Nu_weight=5.53447020181e-08  #5.534702e-08
        self.WWTo1L1Nu2Q_weight=1.14858944695e-07 #1.14858944695e-07
        self.WZJToLLLNu_weight=3.33096092079e-07 #3.18878868498e-07 #7.77410062009e-07
        self.WZTo1L1Nu2Q_weight=2.54725706507e-08 #2.55418584842e-08  #2.54725706507e-08
        self.WZTo1L3Nu_weight=3.26075777445e-07# 3.26075777445e-07  #3.26075777445e-07
        self.ZZTo2L2Q_weight=4.13778925604e-08#4.17392810253e-08 #4.13778925604e-08 #4.13778925604e-08  #4.16766846027e-08
        self.ZZTo4L_weight=6.06324080221e-08 #6.10742507181e-08 #5.92779186525e-08 #5.92779186525e-08  #1.3202628761e-06

 #        self.DYlowmass_weight=1.99619334706e-08

        self.WGToLNuG_weight=3.36727421664e-08  #3.36727421664e-08 #1.02004951008e-07 #1.03976258739e-07
        self.WGstarToLNuEE_weight=1.56725042226e-06 #1.56725042226e-06   #1.56725042226e-06
        self.WGstarToLNuMuMu_weight=1.25890428e-06  #1.25890428e-06    #1.25890428e-06 #1.25890428e-06

        self.ST_tW_antitop_weight=5.1347926337e-06 #5.23465826064e-06
        self.ST_tW_top_weight=5.12021723529e-06   #5.16316171138e-06
        self.ST_t_antitop_weight=6.75839027873e-07  #6.77839939377e-07
        self.ST_t_top_weight=6.55405568597e-07 #6.57612514709e-07
        self.TT_weight=1.08602238426e-05 #1.07699999667e-05 #1.07699999667e-05  #1.08709111195e-05
        self.TTevtgen_weight=8.56505206262e-05 #8.41414729168e-05 #8.80724961084e-05  #1.08709111195e-05

        self.WW_weight= 1.48725639285e-05  #1.49334492783e-05
        self.WZ_weight= 1.17948019785e-05  #1.17948019785e-05
        self.ZZ_weight= 8.3109585141e-06  #4.14537072254e-06

        self.GluGluHTo_weight=2.05507004808e-06  #2.07059122633e-06 
        self.GluGluHToWW_weight=2.12615857113e-05
        self.VBFHToWW_weight=4.5508275481e-07

        self.ZHToTauTau_weight=1.28035908968e-07
        self.WplusHToTauTau_weight=2.33513589465e-07
        self.WminusHToTauTau_weight=3.59334302781e-07

        self.GluGlu_LFV_HToMuTau_weight=1.9432e-06            #1.9432e-06 
        self.VBFHTo_weight=4.23479177085e-08                                   #4.23479177085e-08 
        self.VBF_LFV_HToMuTau_weight=4.05154959453e-08   #4.05154959453e-08 
        self.GluGlu_LFV_HToETau_weight=1.9432e-06    #1.9432e-06 
        self.VBF_LFV_HToETau_weight=4.05239613191e-08  #1.83818883478e-07  

        self.tree = EMTree(tree)
        self.output=outfile
        self.histograms = {}
        self.mym1='m'
        self.mye1='e'

        if self.is_WJet:
            self.binned_weight=[0.709390278,0.190063899,0.060538503,0.01920696,0.019297825]
#0.709390278,0.190063899,0.058529964,0.019206445,0.01923548]
#0.709521921,0.190073347,0.059034569,0.019318685,0.019344044]

        elif self.is_DYJet:
            self.binned_weight=[0.040812575,0.012877836,0.013147445,0.013522627,0.011065415]
#0.117315804,0.016214144,0.016643878,0.017249744,0.01344205]
#0.119211763,0.016249863,0.016824004,0.017799422,0.014957552]
        elif self.is_ZTauTau:
            self.binned_weight=[0.040812575,0.012877836,0.013147445,0.013522627,0.011065415]
#            self.binned_weight=[0.117315804,0.016214144,0.016643878,0.017249744,0.01344205]
#0.119211763,0.016249863,0.016824004,0.017799422,0.014957552]

        elif self.is_DYlowmass:
            self.binned_weight=[0.527321457,0.011178183,0.009311641,0.527321457,0.527321457]
#0.527321457,0.011589863,0.009311641,0.527321457,0.527321457]

        else:
            self.binned_weight=[1,1,1,1,1]
        """need to think about this"""


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
        electron_Pt=self.my_elec.Pt()
        electron_Eta=self.my_elec.Eta()
        muon_Pt=self.my_muon.Pt()
        muon_Eta=self.my_muon.Eta()
#        print electron_Eta," ",muon_Eta
        muidcorr = mid_corrector(muon_Pt, abs(muon_Eta))
        muisocorr = miso_corrector(muon_Pt, abs(muon_Eta))
        mutrcorr = trg_corrector(muon_Pt, abs(muon_Eta))
        mutrkcorr=mtrk_corrector(muon_Eta)[0]
        eidcorr = eId_corrector(electron_Eta,electron_Pt)
        ereconcorr=erecon_corrector(electron_Eta,electron_Pt)
#        eidisocorr0p10= eidiso_corr0p10(getattr(row, self.mye1+'Pt'),abs(getattr(row,self.mye1+'Eta')))[0]
#        eidisocorr0p15= eidiso_corr0p15(getattr(row, self.mye1+'Pt'),abs(getattr(row,self.mye1+'Eta')))[0]
#        etrkcorr=etrk_corrector(getattr(row,self.mye1+'Eta'),getattr(row, self.mye1+'Pt'))
#        print "id corr", muidcorr
#        print "iso corr", muisocorr
#        print "pu  ",pu
#        print "trk corr",mutrkcorr
#        print "tr corr", mutrcorr
#        print "eidiso corr", eidcorr
#        print "etack ",ereconcorr
###       mutrcorr=1
     # if pu*muidcorr1*muisocorr1*muidcorr2*muisocorr2*mutrcorr==0: print pu, muidcorr1, muisocorr1, muidcorr2, muisocorr2, mutrcorr
 #       if pu>2:
#            print "pileup--------   =",pu
   #     print pu*muidcorr*muisocorr*mutrcorr
#        print eisocorr

        topptreweight=1

        if self.isTT:
            topptreweight=topPtreweight(row.topQuarkPt1,row.topQuarkPt2)

        return pu*muidcorr*muisocorr*mutrcorr*mutrkcorr*topptreweight*eidcorr*ereconcorr



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
        self.mVis_=array.array( 'f', [ 0 ] )
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
            self.treeS.Branch("mVis_",self.mVis_,"mVis_/F")
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
            self.treeB.Branch("mVis_",self.mVis_,"mVis_/F")
            self.treeB.Branch("pZeta_",self.pZeta_,"pZeta_/F")
            self.treeB.Branch("lepAsym_",self.lepAsym_,"lepAsym_/F")


        


    def fill_tree(self, row,btagweight=1,region='signal'):

        
        if self.is_WJet or self.is_DYJet or self.is_DYlowmass or self.is_ZTauTau:
            weight = self.event_weight(row,region) *self.binned_weight[int(row.numGenJets)]*0.001
        elif self.isWGToLNuG:
            weight=self.WGToLNuG_weight*self.event_weight(row,region) 
        elif self.isWGstarToLNuEE:
            weight=self.WGstarToLNuEE_weight*self.event_weight(row,region) 
        elif self.isWGstarToLNuMuMu:
            weight=self.WGstarToLNuMuMu_weight*self.event_weight(row,region) 
        elif self.isST_tW_top:
            weight=self.ST_tW_top_weight*self.event_weight(row,region) 
        elif self.isST_tW_antitop:
            weight=self.ST_tW_antitop_weight*self.event_weight(row,region) 
        elif self.isST_t_top:
            weight=self.ST_t_top_weight*self.event_weight(row,region) 
        elif self.isST_t_antitop:
            weight=self.ST_t_antitop_weight*self.event_weight(row,region) 
        elif self.isWW:
            weight=self.WW_weight*self.event_weight(row,region) 
        elif self.isWZ:
            weight=self.WZ_weight*self.event_weight(row,region) 
        elif self.isZZ:
            weight=self.ZZ_weight*self.event_weight(row,region) 
        elif self.isWZTo2L2Q:
            weight=self.WZTo2L2Q_weight*self.event_weight(row,region) 
        elif self.isVVTo2L2Nu:
            weight=self.VVTo2L2Nu_weight*self.event_weight(row,region) 
        elif self.isWWTo1L1Nu2Q:
            weight=self.WWTo1L1Nu2Q_weight*self.event_weight(row,region) 
        elif self.isWZJToLLLNu:
            weight=self.WZJToLLLNu_weight*self.event_weight(row,region) 
        elif self.isWZTo1L1Nu2Q:
            weight=self.WZTo1L1Nu2Q_weight*self.event_weight(row,region) 
        elif self.isWZTo1L3Nu:
            weight=self.WZTo1L3Nu_weight*self.event_weight(row,region) 
        elif self.isZZTo2L2Q:
            weight=self.ZZTo2L2Q_weight*self.event_weight(row,region) 
        elif self.isZZTo4L:
            weight=self.ZZTo4L_weight*self.event_weight(row,region) 
        elif self.isTT:
            weight=self.TT_weight*self.event_weight(row,region) 
        elif self.isTTevtgen:
            weight=self.TTevtgen_weight*self.event_weight(row,region) 
        elif self.isGluGluHTo:
            weight=self.GluGluHTo_weight*self.event_weight(row,region) 
        elif self.isWminusHToTauTau:
            weight=self.WminusHToTauTau_weight*self.event_weight(row,region) 
        elif self.isWplusHToTauTau:
            weight=self.WplusHToTauTau_weight*self.event_weight(row,region) 
        elif self.isVBFHToWW:
            weight=self.VBFHToWW_weight*self.event_weight(row,region) 
        elif self.isGluGluHToWW:
            weight=self.GluGluHToWW_weight*self.event_weight(row,region) 
        elif self.isZHToTauTau:
            weight=self.ZHToTauTau_weight*self.event_weight(row,region) 
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
        self.mVis_[0]=row.e_m_Mass
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
            if not bool(row.singleIsoMu24Pass or row.singleIsoTkMu24Pass): 
                continue   

            #vetoes and cleaning

            if deltaR(row.ePhi,row.mPhi,row.eEta,row.mEta)<0.3:continue

            
            if row.muVetoPt5IsoIdVtx :continue


            if row.eVetoMVAIsoVtx :continue


            if row.tauVetoPt20Loose3HitsVtx : continue
            
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

            self.my_muon=ROOT.TLorentzVector()
            self.my_muon.SetPtEtaPhiM(row.mPt,row.mEta,row.mPhi,row.mMass)
                
            self.my_elec=ROOT.TLorentzVector()
            self.my_elec.SetPtEtaPhiM(row.ePt,row.eEta,row.ePhi,row.eMass)


            #mu preselection
            if not selections.muSelection(row,self.my_muon, 'm'): continue



            #E Preselection
            if not selections.eSelection(row,self.my_elec, 'e'): continue
           
            



            if not selections.lepton_id_iso(row, 'm', 'MuIDTight_mutauiso0p15',dataperiod=self.data_period):continue
 
            if not selections.lepton_id_iso(row, 'e', 'eid15Loose_etauiso0p1',eIDwp='WP80'):continue 
 
            self.fill_tree(row,btagweight)

            
    def finish(self):
        self.output.Write()


