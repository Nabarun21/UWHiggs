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
import FinalStateAnalysis.TagAndProbe.HetauCorrection as HetauCorrection
import FinalStateAnalysis.TagAndProbe.FakeRate2D as FakeRate2D
import bTagSF as bTagSF
from inspect import currentframe
from FinalStateAnalysis.StatTools.RooFunctorFromWS import FunctorFromMVA




cut_flow_step=['allEvents','HLTIsoPasstrg','DR_e_mu','surplus_mu_veto','surplus_e_veto','surplus_tau_veto','musel','mulooseiso','esel','elooseiso','ecalgap','bjetveto','muiso','eiso','jet0sel','jet1sel','jet2loosesel','jet2tightsel']

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
mtrk_corrector = MuonPOGCorrections.mu_trackingEta_2016
#trk_corrector =  MuonPOGCorrections.make_muonptabove10_pog_tracking_corrections_2016()
#eId_corrector = EGammaPOGCorrections.make_egamma_pog_electronID_ICHEP2016( 'nontrigWP80')
eId_corrector = EGammaPOGCorrections.make_egamma_pog_electronID_MORIOND2017( 'nontrigWP80')
erecon_corrector=EGammaPOGCorrections.make_egamma_pog_recon_MORIOND17()

etrk_corrector=EGammaPOGCorrections.make_egamma_pog_tracking_ICHEP2016()
eidiso_corr0p10 =HetauCorrection.idiso0p10_ele_2016BtoHReReco
eidiso_corr0p15 =HetauCorrection.idiso0p15_ele_2016BtoHReReco



class LFVHEMuAnalyzerMVABDT1(MegaBase):
    tree = 'em/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):

        self.channel='EMu'
        super(LFVHEMuAnalyzerMVABDT1, self).__init__(tree, outfile, **kwargs)
        target = os.path.basename(os.environ['megatarget'])
        self.target=target
        self.var =['mTmuMet_','mTeMet_','deltaPhimue_','eDphiPFMet_','mDphiPFMet_','mColl_','pZeta_','mPt_','ePt_'] 
        self.xml_name = os.path.join(os.getcwd(),"BDTfitBDT/weights/TMVAClassification_BDT.weights.xml")  #weights from BDT
        self.functor = FunctorFromMVA('BDT method',self.xml_name, *self.var)

        self.is_WJet=('WJetsToLNu' in target or 'W1JetsToLNu' in target or 'W2JetsToLNu' in target or 'W3JetsToLNu' in target or 'W4JetsToLNu' in target)
        self.is_DYJet= ('DYJetsToLL_M-50' in target or  'DY1JetsToLL_M-50' in target or 'DY2JetsToLL_M-50' in target or 'DY3JetsToLL_M-50' in target or 'DY4JetsToLL_M-50' in target) 

        self.is_DYlowmass= ('DYJetsToLL_M-10to50' in target or  'DY1JetsToLL_M-10to50' in target or 'DY2JetsToLL_M-10to50' in target or 'DY3JetsToLL_M-10to50' in target or 'DY4JetsToLL_M-10to50' in target) 

        self.is_ZTauTau= ('ZTauTauJets_M-50' in target or  'ZTauTau1Jets_M-50' in target or 'ZTauTau2Jets_M-50' in target or 'ZTauTau3Jets_M-50' in target or 'ZTauTau4Jets_M-50' in target) 
        
        self.data_period="BCDEF" if ("Run2016B" in target or "Run2016C" in target or  "Run2016D" in target or  "Run2016E" in target or  "Run2016F" in target) else "GH"

        self.isData=('data' in target)

        self.isWGToLNuG=( 'WGToLNuG' in target)
        self.isWGstarToLNuEE=('WGstarToLNuEE' in target)
        self.isWGstarToLNuMuMu=('WGstarToLNuMuMu' in target)

        self.isST_tW_antitop=('ST_tW_antitop' in target)
        self.isST_tW_top=('ST_tW_top' in target)
        self.isST_t_antitop=('ST_t_antitop' in target)
        self.isST_t_top=('ST_t_top' in target)
        self.isTT=('TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_v6-v1' in target)
        self.isTTevtgen=('TT_TuneCUETP8M2T4_13TeV-powheg-pythia8-evtgen_v6-v1' in target)
        self.isWW=('WW_Tune' in target)
        self.isWZ=('WZ_Tune' in target)
        self.isZZ=('ZZ_Tune' in target)

        self.isGluGluHTo=('GluGluHTo' in target)
        self.isGluGlu_LFV=('GluGlu_LFV_HToMuTau' in target)
        self.isGluGluEtauSig=('GluGlu_LFV_HToETau' in target)
        self.isVBFHTo=('VBFHTo' in target)
        self.isVBF_LFV=('VBF_LFV_HToMuTau' in target)
        self.isVBFEtauSig=('VBF_LFV_HToETau' in target)

        self.isWZTo2L2Q=('WZTo2L2Q' in target)
        self.isVVTo2L2Nu=('VVTo2L2Nu' in target)
        self.isWWTo1L1Nu2Q=('WWTo1L1Nu2Q' in target)
        self.isWZJToLLLNu=('WZJToLLLNu' in target)
        self.isWZTo1L1Nu2Q=('WZTo1L1Nu2Q' in target)
        self.isWZTo1L3Nu=('WZTo1L3Nu' in target)
        self.isZZTo2L2Q=('ZZTo2L2Q' in target)
        self.isZZTo4L=('ZZTo4L' in target)

        self.WZTo2L2Q_weight=2.40930214376e-08 #8.91328638844e-07
        self.VVTo2L2Nu_weight=5.53447020181e-08 #5.534702e-08
        self.WWTo1L1Nu2Q_weight=1.14858944695e-07#1.14858944695e-07
        self.WZJToLLLNu_weight=3.18878868498e-07#7.77410062009e-07
        self.WZTo1L1Nu2Q_weight=2.55418584842e-08 #2.54725706507e-08
        self.WZTo1L3Nu_weight=3.26075777445e-07 #3.26075777445e-07
        self.ZZTo2L2Q_weight=4.13778925604e-08 #4.16766846027e-08
        self.ZZTo4L_weight=5.92779186525e-08 #1.3202628761e-06

#        self.DYlowmass_weight=1.99619334706e-08

        self.WGToLNuG_weight=1.02004951008e-07#1.03976258739e-07
        self.WGstarToLNuEE_weight=1.56725042226e-06
        self.WGstarToLNuMuMu_weight=1.25890428e-06#1.25890428e-06

        self.ST_tW_antitop_weight=5.23465826064e-06#5.1347926337e-06
        self.ST_tW_top_weight=5.16316171138e-06
        self.ST_t_antitop_weight=6.77839939377e-07
        self.ST_t_top_weight=6.57612514709e-07
        self.TT_weight=1.08709111195e-05
        self.TTevtgen_weight=8.54641838731e-05#1.08709111195e-05

        self.WW_weight= 1.49334492783e-05#0.000119415057363
        self.WZ_weight= 1.17948019785e-05#4.713e-05
        self.ZZ_weight=4.14537072254e-06#1.66888201167e-05

        self.GluGluHTo_weight=2.07059122633e-06#2.07059122633e-06
        self.GluGlu_LFV_HToMuTau_weight=1.9432e-06#1.9432e-06
        self.VBFHTo_weight=4.23479177085e-08#4.23479177085e-08
        self.VBF_LFV_HToMuTau_weight= 4.05154959453e-08#2.26629913711e-07 
        self.GluGlu_LFV_HToETau_weight=1.9432e-06#1.9432e-06
        self.VBF_LFV_HToETau_weight=4.05239613191e-08#1.83818883478e-07  


        self.tree = EMTree(tree)
        self.out=outfile
        self.histograms = {}
        self.mym1='m'
        self.mye1='e'
        self.sysdir=['nosys','jetup','jetdown','tup','tdown','uup','udown']
#        self.sysdir=['nosys']
        if self.is_WJet:
            self.binned_weight=[0.709521921,0.190073347,0.059034569,0.019318685,0.019344044]
#0.709763222,0.19009066,0.108822068,0.054793711,0.019344223]
        elif self.is_DYJet:
            self.binned_weight=[0.119211763,0.016249863,0.016824004,0.017799422,0.014957552]
#0.117315804,0.016233898,0.016725261,0.017545741,0.021544426]
        elif self.is_ZTauTau:
            self.binned_weight=[0.119211763,0.016249863,0.016824004,0.017799422,0.014957552]

        elif self.is_DYlowmass:
            self.binned_weight=[0.527321457,0.011589863,0.009311641,0.527321457,0.527321457]
#0.117315804,0.016233898,0.016725261,0.017545741,0.021544426]
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
        muidcorr = mid_corrector(getattr(row, self.mym1+'Pt'), abs(getattr(row, self.mym1+'Eta')))
        muisocorr = miso_corrector(getattr(row, self.mym1+'Pt'), abs(getattr(row, self.mym1+'Eta')))
        mutrcorr = trg_corrector(getattr(row, self.mym1+'Pt'), abs(getattr(row, self.mym1+'Eta'))) 
        mutrkcorr=mtrk_corrector(getattr(row,self.mym1+'Eta'))[0]
        eidcorr = eId_corrector(getattr(row,self.mye1+'Eta'),getattr(row, self.mye1+'Pt'))
        ereconcorr=erecon_corrector(getattr(row,self.mye1+'Eta'),getattr(row, self.mye1+'Pt'))
#        eidisocorr0p10= eidiso_corr0p10(getattr(row, self.mye1+'Pt'),abs(getattr(row,self.mye1+'Eta')))[0]
#        eidisocorr0p15= eidiso_corr0p15(getattr(row, self.mye1+'Pt'),abs(getattr(row,self.mye1+'Eta')))[0]
#        etrkcorr=etrk_corrector(getattr(row,self.mye1+'Eta'),getattr(row, self.mye1+'Pt'))
#        print "id corr", muidcorr
#        print "iso corr", muisocorr
#        print "pu  ",pu
#        print "trk corr",mutrkcorr
#        print "tr corr", mutrcorr
#        print "eidiso corr", eidisocorr0p10
#        print "etack ",etrkcorr
##       mutrcorr=1
     # if pu*muidcorr1*muisocorr1*muidcorr2*muisocorr2*mutrcorr==0: print pu, muidcorr1, muisocorr1, muidcorr2, muisocorr2, mutrcorr
 #       if pu>2:
#            print "pileup--------   =",pu
   #     print pu*muidcorr*muisocorr*mutrcorr
#        print eisocorr

        topptreweight=1

        if self.isTT:
            topptreweight=topPtreweight(row.topQuarkPt1,row.topQuarkPt2)

        return pu*muidcorr*muisocorr*mutrcorr*mutrkcorr*topptreweight*eidcorr*ereconcorr




 
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

        processtype=['gg']
#        threshold=[]
        sign=[ 'ss','os']
        jetN = [0,1,21,22,3]
        folder=[]
#        pudir = ['','p1s/', 'm1s/','trp1s/', 'trm1s/', 'eidp1s/','eidm1s/',  'eisop1s/','eisom1s/', 'mLoose/','mLooseUp/','mLooseDown/', ]

#        alldirs=['','antiIsolatedweighted/','antiIsolated/','antiIsolatedweightedelectron/','antiIsolatedweightedmuon/','antiIsolatedweightedmuonelectron/','fakeRateMethod/']
        alldirs=['','antiIsolatedweighted/','antiIsolated/','qcdshaperegion/']
#        alldirs=['']
        for d  in alldirs :
            for i in sign:
                for j in processtype:
 #                   for k in threshold:
                    for jn in jetN: 
                        folder.append(d+i+'/'+j+'/'+str(jn))
                        for s in self.sysdir:
                            folder.append(d+i+'/'+j+'/'+str(jn)+'/selected/'+s)
                        
                            

        for k in range(len(folder)):
            f=folder[k]
            if 'selected' not in f:
                self.book(f,"mPt", "mu p_{T}", 200, 0, 200)
                self.book(f,"mPhi", "mu phi", 100, -3.2, 3.2)
                self.book(f,"mEta", "mu eta",  50, -2.5, 2.5)
            
                self.book(f,"ePt", "e p_{T}", 200, 0, 200)
                self.book(f,"ePhi", "e phi",  100, -3.2, 3.2)
                self.book(f,"eEta", "e eta", 50, -2.5, 2.5)
            
                self.book(f, "em_DeltaPhi", "e-mu DeltaPhi" , 50, 0, 3.2)
                self.book(f, "em_DeltaR", "e-mu DeltaR" , 100, -7, 7)
#            
                self.book(f, "h_collmass_pfmet",  "h_collmass_pfmet",  300, 0, 300) 

                self.book(f, "BDT_value",  "BDT_value",  60, -0.6, 0.6)
#                self.book(f, "h_collmass_mvamet",  "h_collmass_mvamet",  300, 0, 300) 
#                
#                self.book(f,"scaledmPt","scaledmPt",60,0,1.2)
#                self.book(f,"scaledePt","scaledePt",60,0,1.2)
#
                self.book(f, "Met",  "Met",  30, 0, 300) 
                self.book(f, "pZeta",  "pZeta",  60,-300, 300) 
                self.book(f, "eMtToPfMet",  "eMtToPfMet",  300, 0, 300) 
                self.book(f, "mMtToPfMet",  "mMtToPfMet",  300, 0, 300) 
                self.book(f, "dPhiMetToE",  "dPhiMetToE",  50, 0, 3.2)
                
                self.book(f, "h_vismass",  "h_vismass",  300, 0, 300) 
                self.book(f, "mPFMET_Mt", "mu-PFMET M_{T}" , 200, 0, 200)
                self.book(f, "ePFMET_Mt", "e-PFMET M_{T}" , 200, 0, 200)           
                self.book(f, "mPFMET_DeltaPhi", "mu-PFMET DeltaPhi" , 50, 0, 3.2)
                self.book(f, "ePFMET_DeltaPhi", "e-PFMET DeltaPhi" , 50, 0, 3.2)
#                self.book(f, "mPFMETDeltaPhi_vs_ePFMETDeltaPhi",  "mPFMETDeltaPhi_vs_ePFMETDeltaPhi", 50, 0, 3.2, 50,0, 3.2, type=ROOT.TH2F)
                self.book(f, "vbfMass","vbf dijet mass",500,0,5000)
                self.book(f, "vbfDeta","vbf Delta Eta",50,0,5)

#                self.book(f, "jetN_20", "Number of jets, p_{T}>20", 10, -0.5, 9.5) 
                self.book(f, "jetN_30", "Number of jets, p_{T}>30", 10, -0.5, 9.5) 
            else:
#                self.book(f, "h_collmass_pfmet",  "h_collmass_pfmet",  300, 0, 300) 
                self.book(f, "Met",  "Met",  30, 0, 300) 
                self.book(f, "BDT_value",  "BDT_value",  95, -0.6, 0.35)
        for s in sign:
            self.book(s, "jetN_30", "Number of jets, p_{T}>30", 10, -0.5, 9.5) 
            self.book(s, "NUP", "Number of Partons", 12, -0.5, 11.5) 
            self.book(s, "numGenJets", "Number of Gen Level Jets", 12, -0.5, 11.5) 
            self.book(s, "numVertices", "Number of Vertices", 60, 0, 60) 
            self.book(s, "h_collmass_pfmet", "h_collmass_pfmet", 300, 0, 300) 
            self.book(s, "h_vismass",  "h_vismass",  300, 0, 300) 
            self.book(s+'/tNoCuts', "CUT_FLOW", "Cut Flow", len(cut_flow_step), 0, len(cut_flow_step))
            
            xaxis = self.histograms[s+'/tNoCuts/CUT_FLOW'].GetXaxis()
            self.cut_flow_histo = self.histograms[s+'/tNoCuts/CUT_FLOW']
            self.cut_flow_map   = {}
            for i, name in enumerate(cut_flow_step):
                xaxis.SetBinLabel(i+1, name)
                self.cut_flow_map[name] = i+0.5


    def fill_jet_histos(self,row,sign,region,btagweight):

        if self.is_WJet or self.is_DYJet:
            weight1 = self.event_weight(row,region)*self.binned_weight[int(row.numGenJets)]
        else:
            weight1 = self.event_weight(row,region)

        weight1=btagweight*weight1
        self.histograms[sign+'/jetN_30'].Fill(row.jetVeto30,weight1)
        self.histograms[sign+'/NUP'].Fill(row.NUP,weight1)
        self.histograms[sign+'/numGenJets'].Fill(row.numGenJets,weight1)


    def get_fakerate(self,row):
        if row.ePt<10:
            raise ValueError("electron Pt less than 10")
        elecPt=row.ePt if row.ePt<200 else 200
        if abs(row.eEta)<1.479: #par[0]+par[1]*TMath::Erf(par[2]*x - par[3])  Barrel
            fakerateFactor=0.449+0.289*ROOT.TMath.Erf(0.043*elecPt-2.686)
        else:#endcap par[0]+par[1]*x
            fakerateFactor=0.268+0.001*elecPt

        return fakerateFactor/(1-fakerateFactor)
            
            

    def fill_histos(self, row,MVAvalue,sign,f=None,fillCommonOnly=False,region=None,btagweight=1,sys='',qcdshape=None):


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
            
        weight=btagweight*weight
        #fill some histograms before dividing into categories (looking at some variables before dividing into njet categories makes more sense)
        if fillCommonOnly==True:
            if region=='signal':
                self.histograms[sign+'/jetN_30'].Fill(row.jetVeto30,weight)
                self.histograms[sign+'/NUP'].Fill(row.NUP,weight)
                self.histograms[sign+'/numGenJets'].Fill(row.numGenJets,weight)
                self.histograms[sign+'/numVertices'].Fill(row.nvtx,weight)
                self.histograms[sign+'/h_collmass_pfmet'].Fill(collmass(row,row.type1_pfMetEt, row.type1_pfMetPhi),weight)
                self.histograms[sign+'/h_vismass'].Fill(row.e_m_Mass,weight)
                return 1
            else:
                return 0

        histos = self.histograms
        pudir=['']
#        pudir =['','fakeRateMethod/'] uncoment if using old fr method
        elooseList = ['antiIsolatedweightedelectron/']
        mlooseList = ['antiIsolatedweightedmuon/']
        emlooseList= ['antiIsolatedweightedmuonelectron/']
        alllooselist=['antiIsolatedweighted/','antiIsolated/']
        qcdshapelist=['qcdshaperegion/']

        if region!='signal':
            for n, l in enumerate(qcdshapelist):
                qcdshapeweight= weight
                folder = l+f
                if sys=='presel':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,qcdshapeweight)
                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row,row.type1_pfMetEt, row.type1_pfMetPhi),qcdshapeweight)
                    histos[folder+'/mPt'].Fill(row.mPt, qcdshapeweight)
                    histos[folder+'/Met'].Fill(row.type1_pfMetEt, qcdshapeweight)
                    histos[folder+'/mEta'].Fill(row.mEta, qcdshapeweight)
                    histos[folder+'/mPhi'].Fill(row.mPhi, qcdshapeweight) 
                    histos[folder+'/ePt'].Fill(row.ePt, qcdshapeweight)
                    histos[folder+'/eEta'].Fill(row.eEta, qcdshapeweight)
                    histos[folder+'/ePhi'].Fill(row.ePhi, qcdshapeweight)
                    histos[folder+'/em_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.mPhi), qcdshapeweight)
                    histos[folder+'/em_DeltaR'].Fill(row.e_m_DR, qcdshapeweight)
                    histos[folder+'/h_vismass'].Fill(row.e_m_Mass, qcdshapeweight)
                    histos[folder+'/pZeta'].Fill(row.e_m_PZetaLess0p85PZetaVis, qcdshapeweight)
                    histos[folder+'/ePFMET_Mt'].Fill(row.eMtToPfMet_type1, qcdshapeweight)
                    histos[folder+'/mPFMET_Mt'].Fill(row.mMtToPfMet_type1, qcdshapeweight)
                    histos[folder+'/ePFMET_DeltaPhi'].Fill(abs(row.eDPhiToPfMet_type1), qcdshapeweight)
                    histos[folder+'/mPFMET_DeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1), qcdshapeweight)
#                    histos[folder+'/mPFMETDeltaPhi_vs_ePFMETDeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1),abs(row.eDPhiToPfMet_type1) , qcdshapeweight)
                    histos[folder+'/vbfMass'].Fill(row.vbfMass, qcdshapeweight)
                    histos[folder+'/vbfDeta'].Fill(row.vbfDeta, qcdshapeweight)
                    histos[folder+'/jetN_30'].Fill(row.jetVeto30, qcdshapeweight) 
                elif sys=='nosys':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,qcdshapeweight)
                    histos[folder+'/Met'].Fill(row.type1_pfMetEt, qcdshapeweight)
                elif sys=='jetup':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,qcdshapeweight)

                elif sys=='jetdown':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,qcdshapeweight)

                elif sys=='tup':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,qcdshapeweight)

                elif sys=='tdown':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,qcdshapeweight)

                elif sys=='uup':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,qcdshapeweight)

                elif sys=='udown':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,qcdshapeweight)
        
        """
        if region=='eLoosemTight':
            frarray=self.get_fakerate(row) 
            fakerateWeight=frarray
            antiIsolatedWeightList=[fakerateWeight,1]
            for n, l in enumerate(elooseList) :
                antiIsolatedWeight= weight*antiIsolatedWeightList[n]
                folder = l+f
                if sys=='presel':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,antiIsolatedWeight)
                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row,row.type1_pfMetEt, row.type1_pfMetPhi),antiIsolatedWeight)
                    histos[folder+'/mPt'].Fill(row.mPt, antiIsolatedWeight)
                    histos[folder+'/Met'].Fill(row.type1_pfMetEt, antiIsolatedWeight)
                    histos[folder+'/mEta'].Fill(row.mEta, antiIsolatedWeight)
                    histos[folder+'/mPhi'].Fill(row.mPhi, antiIsolatedWeight) 
                    histos[folder+'/ePt'].Fill(row.ePt, antiIsolatedWeight)
                    histos[folder+'/eEta'].Fill(row.eEta, antiIsolatedWeight)
                    histos[folder+'/ePhi'].Fill(row.ePhi, antiIsolatedWeight)
                    histos[folder+'/em_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.mPhi), antiIsolatedWeight)
                    histos[folder+'/em_DeltaR'].Fill(row.e_m_DR, antiIsolatedWeight)
                    histos[folder+'/h_vismass'].Fill(row.e_m_Mass, antiIsolatedWeight)
                    histos[folder+'/ePFMET_Mt'].Fill(row.eMtToPfMet_type1, antiIsolatedWeight)
                    histos[folder+'/mPFMET_Mt'].Fill(row.mMtToPfMet_type1, antiIsolatedWeight)
                    histos[folder+'/ePFMET_DeltaPhi'].Fill(abs(row.eDPhiToPfMet_type1), antiIsolatedWeight)
                    histos[folder+'/mPFMET_DeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1), antiIsolatedWeight)
#                    histos[folder+'/mPFMETDeltaPhi_vs_ePFMETDeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1),abs(row.eDPhiToPfMet_type1) , antiIsolatedWeight)
                    histos[folder+'/vbfMass'].Fill(row.vbfMass, antiIsolatedWeight)
                    histos[folder+'/vbfDeta'].Fill(row.vbfDeta, antiIsolatedWeight)
                    histos[folder+'/jetN_30'].Fill(row.jetVeto30, antiIsolatedWeight) 
                elif sys=='nosys':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,antiIsolatedWeight)

        if region=='eTightmLoose':
            fakerateWeight=4
            antiIsolatedWeightList=[fakerateWeight,1]
            for n, l in enumerate(mlooseList) :
                antiIsolatedWeight= weight*antiIsolatedWeightList[n]
                folder = l+f
                if sys=='presel':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,antiIsolatedWeight)
                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row,row.type1_pfMetEt, row.type1_pfMetPhi),antiIsolatedWeight)
                    histos[folder+'/mPt'].Fill(row.mPt, antiIsolatedWeight)
                    histos[folder+'/Met'].Fill(row.type1_pfMetEt, antiIsolatedWeight)
                    histos[folder+'/mEta'].Fill(row.mEta, antiIsolatedWeight)
                    histos[folder+'/mPhi'].Fill(row.mPhi, antiIsolatedWeight) 
                    histos[folder+'/ePt'].Fill(row.ePt, antiIsolatedWeight)
                    histos[folder+'/eEta'].Fill(row.eEta, antiIsolatedWeight)
                    histos[folder+'/ePhi'].Fill(row.ePhi, antiIsolatedWeight)
                    histos[folder+'/em_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.mPhi), antiIsolatedWeight)
                    histos[folder+'/em_DeltaR'].Fill(row.e_m_DR, antiIsolatedWeight)
                    histos[folder+'/h_vismass'].Fill(row.e_m_Mass, antiIsolatedWeight)
                    histos[folder+'/ePFMET_Mt'].Fill(row.eMtToPfMet_type1, antiIsolatedWeight)
                    histos[folder+'/mPFMET_Mt'].Fill(row.mMtToPfMet_type1, antiIsolatedWeight)
                    histos[folder+'/ePFMET_DeltaPhi'].Fill(abs(row.eDPhiToPfMet_type1), antiIsolatedWeight)
                    histos[folder+'/mPFMET_DeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1), antiIsolatedWeight)
#                    histos[folder+'/mPFMETDeltaPhi_vs_ePFMETDeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1),abs(row.eDPhiToPfMet_type1) , antiIsolatedWeight)
                    histos[folder+'/vbfMass'].Fill(row.vbfMass, antiIsolatedWeight)
                    histos[folder+'/vbfDeta'].Fill(row.vbfDeta, antiIsolatedWeight)
                    histos[folder+'/jetN_30'].Fill(row.jetVeto30, antiIsolatedWeight) 
                elif sys=='nosys':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,antiIsolatedWeight)

        if region=='eLoosemLoose':
            efrarray=self.get_fakerate(row) 
            efakerateWeight=efrarray
            mfakerateWeight=4
            antiIsolatedWeightList=[mfakerateWeight*efakerateWeight,1]
            for n, l in enumerate(emlooseList) :
                antiIsolatedWeight= weight*antiIsolatedWeightList[n]
                folder = l+f
                if sys=='presel':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,antiIsolatedWeight)
                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row,row.type1_pfMetEt, row.type1_pfMetPhi),antiIsolatedWeight)
                    histos[folder+'/mPt'].Fill(row.mPt, antiIsolatedWeight)
                    histos[folder+'/Met'].Fill(row.type1_pfMetEt, antiIsolatedWeight)
                    histos[folder+'/mEta'].Fill(row.mEta, antiIsolatedWeight)
                    histos[folder+'/mPhi'].Fill(row.mPhi, antiIsolatedWeight) 
                    histos[folder+'/ePt'].Fill(row.ePt, antiIsolatedWeight)
                    histos[folder+'/eEta'].Fill(row.eEta, antiIsolatedWeight)
                    histos[folder+'/ePhi'].Fill(row.ePhi, antiIsolatedWeight)
                    histos[folder+'/em_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.mPhi), antiIsolatedWeight)
                    histos[folder+'/em_DeltaR'].Fill(row.e_m_DR, antiIsolatedWeight)
                    histos[folder+'/h_vismass'].Fill(row.e_m_Mass, antiIsolatedWeight)
                    histos[folder+'/ePFMET_Mt'].Fill(row.eMtToPfMet_type1, antiIsolatedWeight)
                    histos[folder+'/mPFMET_Mt'].Fill(row.mMtToPfMet_type1, antiIsolatedWeight)
                    histos[folder+'/ePFMET_DeltaPhi'].Fill(abs(row.eDPhiToPfMet_type1), antiIsolatedWeight)
                    histos[folder+'/mPFMET_DeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1), antiIsolatedWeight)
#                    histos[folder+'/mPFMETDeltaPhi_vs_ePFMETDeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1),abs(row.eDPhiToPfMet_type1) , antiIsolatedWeight)
                    histos[folder+'/vbfMass'].Fill(row.vbfMass, antiIsolatedWeight)
                    histos[folder+'/vbfDeta'].Fill(row.vbfDeta, antiIsolatedWeight)
                    histos[folder+'/jetN_30'].Fill(row.jetVeto30, antiIsolatedWeight) 
                elif sys=='nosys':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,antiIsolatedWeight)
                    """
        if region!='signal':
            if region=='eLoosemTight':
                frarray=self.get_fakerate(row)
                fakerateWeight=frarray
            elif region=='eTightmLoose':
                fakerateWeight=4
            elif region=='eLoosemLoose':
                frarray=self.get_fakerate(row)
                efakerateWeight=frarray
                mfakerateWeight=4
                fakerateWeight=-efakerateWeight*mfakerateWeight
            antiIsolatedWeightList=[fakerateWeight,1]
            for n, l in enumerate(alllooselist):
                antiIsolatedWeight= weight*antiIsolatedWeightList[n]
                folder = l+f
                if sys=='presel':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,antiIsolatedWeight)
                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row,row.type1_pfMetEt, row.type1_pfMetPhi),antiIsolatedWeight)
                    histos[folder+'/mPt'].Fill(row.mPt, antiIsolatedWeight)
                    histos[folder+'/Met'].Fill(row.type1_pfMetEt, antiIsolatedWeight)
                    histos[folder+'/mEta'].Fill(row.mEta, antiIsolatedWeight)
                    histos[folder+'/mPhi'].Fill(row.mPhi, antiIsolatedWeight) 
                    histos[folder+'/ePt'].Fill(row.ePt, antiIsolatedWeight)
                    histos[folder+'/eEta'].Fill(row.eEta, antiIsolatedWeight)
                    histos[folder+'/ePhi'].Fill(row.ePhi, antiIsolatedWeight)
                    histos[folder+'/em_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.mPhi), antiIsolatedWeight)
                    histos[folder+'/em_DeltaR'].Fill(row.e_m_DR, antiIsolatedWeight)
                    histos[folder+'/h_vismass'].Fill(row.e_m_Mass, antiIsolatedWeight)
                    histos[folder+'/pZeta'].Fill(row.e_m_PZetaLess0p85PZetaVis, antiIsolatedWeight)
                    histos[folder+'/ePFMET_Mt'].Fill(row.eMtToPfMet_type1, antiIsolatedWeight)
                    histos[folder+'/mPFMET_Mt'].Fill(row.mMtToPfMet_type1, antiIsolatedWeight)
                    histos[folder+'/ePFMET_DeltaPhi'].Fill(abs(row.eDPhiToPfMet_type1), antiIsolatedWeight)
                    histos[folder+'/mPFMET_DeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1), antiIsolatedWeight)
#                    histos[folder+'/mPFMETDeltaPhi_vs_ePFMETDeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1),abs(row.eDPhiToPfMet_type1) , antiIsolatedWeight)
                    histos[folder+'/vbfMass'].Fill(row.vbfMass, antiIsolatedWeight)
                    histos[folder+'/vbfDeta'].Fill(row.vbfDeta, antiIsolatedWeight)
                    histos[folder+'/jetN_30'].Fill(row.jetVeto30, antiIsolatedWeight) 
                elif sys=='nosys':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,antiIsolatedWeight)
                    histos[folder+'/Met'].Fill(row.type1_pfMetEt, antiIsolatedWeight)

                elif sys=='jetup':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,antiIsolatedWeight)

                elif sys=='jetdown':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,antiIsolatedWeight)

                elif sys=='tup':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,antiIsolatedWeight)

                elif sys=='tdown':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,antiIsolatedWeight)

                elif sys=='uup':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,antiIsolatedWeight)

                elif sys=='udown':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,antiIsolatedWeight)
            """
            if not self.isData and not self.isGluGlu_LFV and not self.isVBF_LFV and not self.isGluGluEtauSig and not self.isVBFEtauSig:
                if region=='eLoosemTight':
                    frarray=self.get_fakerate(row)
                    efakerateWeight=frarray
                elif region=='eTightmLoose':
                    fakerateWeight=4
                elif region=='eLoosemLoose':
                    frarray=self.get_fakerate(row)
                    efakerateWeight=frarray
                    mfakerateWeight=4
                    fakerateweight=-efakerateWeight*mfakerateWeight
                antiIsolatedWeightList=[fakerateWeight,1]
                fakeRateMethodfolder=pudir[1]+f
                fakeRateMethodweight=fakerateWeight*weight*(-1)

                if sys=='presel':
                    histos[fakeRateMethodfolder+'/BDT_value'].Fill(MVAvalue,fakeRateMethodweight)
                    histos[fakeRateMethodfolder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi),fakeRateMethodweight)
                    histos[fakeRateMethodfolder+'/Met'].Fill(row.type1_pfMetEt, fakeRateMethodweight)
                    histos[fakeRateMethodfolder+'/mPt'].Fill(row.mPt, fakeRateMethodweight)
                    histos[fakeRateMethodfolder+'/mEta'].Fill(row.mEta, fakeRateMethodweight)
                    histos[fakeRateMethodfolder+'/mPhi'].Fill(row.mPhi, fakeRateMethodweight) 
                    histos[fakeRateMethodfolder+'/ePt'].Fill(row.ePt, fakeRateMethodweight)
                    histos[fakeRateMethodfolder+'/eEta'].Fill(row.eEta, fakeRateMethodweight)
                    histos[fakeRateMethodfolder+'/ePhi'].Fill(row.ePhi, fakeRateMethodweight)
                    histos[fakeRateMethodfolder+'/em_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.mPhi), fakeRateMethodweight)
                    histos[fakeRateMethodfolder+'/em_DeltaR'].Fill(row.e_m_DR, fakeRateMethodweight)
                    histos[fakeRateMethodfolder+'/h_vismass'].Fill(row.e_m_Mass, fakeRateMethodweight)
                    histos[fakeRateMethodfolder+'/ePFMET_Mt'].Fill(row.eMtToPfMet_type1, fakeRateMethodweight)
                    histos[fakeRateMethodfolder+'/mPFMET_Mt'].Fill(row.mMtToPfMet_type1, fakeRateMethodweight)
                    histos[fakeRateMethodfolder+'/ePFMET_DeltaPhi'].Fill(abs(row.eDPhiToPfMet_type1), fakeRateMethodweight)
                    histos[fakeRateMethodfolder+'/mPFMET_DeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1), fakeRateMethodweight)
                #               histos[fakeRateMethodfolder+'/mPFMETDeltaPhi_vs_ePFMETDeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1),abs(row.eDPhiToPfMet_type1) , fakeRateMethodweight)
                    histos[fakeRateMethodfolder+'/vbfMass'].Fill(row.vbfMass, fakeRateMethodweight)
                    histos[fakeRateMethodfolder+'/vbfDeta'].Fill(row.vbfDeta, fakeRateMethodweight)
                    histos[fakeRateMethodfolder+'/jetN_30'].Fill(row.jetVeto30, fakeRateMethodweight) 
                elif sys=='nosys':
                    histos[fakeRateMethodfolder+'/BDT_value'].Fill(MVAvalue,fakeRateMethodweight)

                elif sys=='jetup':
                    histos[fakeRateMethodfolder+'/BDT_value'].Fill(MVAvalue,fakeRateMethodweight)

                elif sys=='jetdown':
                    histos[fakeRateMethodfolder+'/BDT_value'].Fill(MVAvalue,fakeRateMethodweight)

                elif sys=='tup':
                    histos[fakeRateMethodfolder+'/BDT_value'].Fill(MVAvalue,fakeRateMethodweight)

                elif sys=='tdown':
                    histos[fakeRateMethodfolder+'/BDT_value'].Fill(MVAvalue,fakeRateMethodweight)

                elif sys=='uup':
                    histos[fakeRateMethodfolder+'/BDT_value'].Fill(MVAvalue,fakeRateMethodweight)

                elif sys=='udown':
                    histos[fakeRateMethodfolder+'/BDT_value'].Fill(MVAvalue,fakeRateMethodweight)
                    """
        if region=='signal' :
            for n,d  in enumerate(pudir) :
                folder = d+f                
                if sys=='presel':
                    histos[folder+'/mPt'].Fill(row.mPt, weight)
                    histos[folder+'/mEta'].Fill(row.mEta, weight)
                    histos[folder+'/mPhi'].Fill(row.mPhi, weight) 
                    histos[folder+'/ePt'].Fill(row.ePt, weight)
                    histos[folder+'/eEta'].Fill(row.eEta, weight)
                    histos[folder+'/ePhi'].Fill(row.ePhi, weight)
                    histos[folder+'/em_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.mPhi), weight)
                    histos[folder+'/em_DeltaR'].Fill(row.e_m_DR, weight)
                    histos[folder+'/h_vismass'].Fill(row.e_m_Mass, weight)

                    histos[folder+'/BDT_value'].Fill(MVAvalue,weight)

                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi),weight)
 #                   histos[folder+'/scaledmPt'].Fill(float(row.mPt)/float(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)),weight)
 #                   histos[folder+'/scaledePt'].Fill(float(row.ePt)/float(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)),weight)
 #
                    histos[folder+'/Met'].Fill(row.type1_pfMetEt, weight)

                    histos[folder+'/pZeta'].Fill(row.e_m_PZetaLess0p85PZetaVis, weight)
                    histos[folder+'/ePFMET_Mt'].Fill(row.eMtToPfMet_type1, weight)
                    histos[folder+'/mPFMET_Mt'].Fill(row.mMtToPfMet_type1, weight)
                    histos[folder+'/ePFMET_DeltaPhi'].Fill(abs(row.eDPhiToPfMet_type1), weight)
                    histos[folder+'/mPFMET_DeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1), weight)
#                    histos[folder+'/mPFMETDeltaPhi_vs_ePFMETDeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1),abs(row.eDPhiToPfMet_type1) , weight)
 
                    histos[folder+'/vbfMass'].Fill(row.vbfMass, weight)
                    histos[folder+'/vbfDeta'].Fill(row.vbfDeta, weight)
                    histos[folder+'/jetN_30'].Fill(row.jetVeto30, weight) 

                elif sys=='nosys':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,weight)
                    histos[folder+'/Met'].Fill(row.type1_pfMetEt, weight)
                elif sys=='jetup':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,weight)

                elif sys=='jetdown':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,weight)

                elif sys=='tup':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,weight)

                elif sys=='tdown':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,weight)

                elif sys=='uup':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,weight)

                elif sys=='udown':
                    histos[folder+'/BDT_value'].Fill(MVAvalue,weight)


    def process(self):
        
        cut_flow_histo = self.cut_flow_histo
        cut_flow_trk   = cut_flow_tracker(cut_flow_histo)
        myevent=()
        frw = []
        curr_event=0
        for row in self.tree:
            sign = 'ss' if row.e_m_SS else 'os'


            cut_flow_trk.new_row(row.run,row.lumi,row.evt)
           
            cut_flow_trk.Fill('allEvents')


            if (self.is_ZTauTau and not row.isZtautau and not self.is_DYlowmass):
                continue
            if (not self.is_ZTauTau and row.isZtautau and not self.is_DYlowmass):
                continue

 #           ptthreshold = [30]
            repeatEvt=True
            if row.evt!=curr_event:
                curr_event=row.evt
                repeatEvt=False
            
            if repeatEvt:continue

#            print "non-repeat"
            processtype ='gg'##changed from 20

            if not bool(row.singleIsoMu24Pass or row.singleIsoTkMu24Pass): 
                continue   

            cut_flow_trk.Fill('HLTIsoPasstrg')

            #vetoes and cleaning

            if deltaR(row.ePhi,row.mPhi,row.eEta,row.mEta)<0.3:continue
            cut_flow_trk.Fill('DR_e_mu')
            
            if row.muVetoPt5IsoIdVtx :continue
            cut_flow_trk.Fill('surplus_mu_veto')

            if row.eVetoMVAIsoVtx :continue
            cut_flow_trk.Fill('surplus_e_veto')

            if row.tauVetoPt20Loose3HitsVtx : continue
            cut_flow_trk.Fill('surplus_tau_veto')

            #mu preselection
            if not selections.muSelection(row, 'm'): continue
            cut_flow_trk.Fill('musel')
            if not selections.lepton_id_iso(row, 'm', 'MuIDTight_idiso0p25',dataperiod=self.data_period): continue
            cut_flow_trk.Fill('mulooseiso')


            #E Preselection
            if not selections.eSelection(row, 'e'): continue
            cut_flow_trk.Fill('esel')

            if not selections.lepton_id_iso(row, 'e', 'eid15Loose_etauiso1',eIDwp='WP80'): continue
            cut_flow_trk.Fill('elooseiso')

            #take care of ecal gap
            if row.eAbsEta > 1.4442 and row.eAbsEta < 1.566 : continue             

            cut_flow_trk.Fill('ecalgap')

            nbtagged=row.bjetCISVVeto30Medium
            if nbtagged>2:
                continue
            btagweight=1
            if (self.isData and nbtagged>0):
                continue
            if nbtagged>0:
                if nbtagged==1:
                    btagweight=bTagSF.bTagEventWeight(nbtagged,row.jb1pt,row.jb1hadronflavor,row.jb2pt,row.jb2hadronflavor,1,0,0) if (row.jb1pt>-990 and row.jb1hadronflavor>-990) else 0
                if nbtagged==2:
                    btagweight=bTagSF.bTagEventWeight(nbtagged,row.jb1pt,row.jb1hadronflavor,row.jb2pt,row.jb2hadronflavor,1,0,0) if (row.jb1pt>-990 and row.jb1hadronflavor>-990 and row.jb2pt>-990 and row.jb2hadronflavor>-990) else 0
#                print "btagweight,nbtagged,row.jb1pt,row.jb1hadronflavor,row.jb2pt,row.jb2hadronflavor"," ",btagweight," ",nbtagged," ",row.jb1pt," ",row.jb1hadronflavor," ",row.jb2pt," ",row.jb2hadronflavor

#            if btagweight<0:btagweight=0

            if btagweight==0: continue

            cut_flow_trk.Fill('bjetveto')
            ## All preselection passed


            
            jetN=row.jetVeto30
            
            if jetN>3:
                jetN=3
            if jetN==2:
                if row.vbfMass>=550:
                    jetN=22
                else:
                    jetN=21


            ##do stuff for qcd shape region
            qcdshaperegion=False
            if jetN==0:
                qcdshaperegion=(row.mRelPFIsoDBDefault>0.1) and (row.mRelPFIsoDBDefault<0.3) and (row.eIsoDB03<0.3)
            else:
                qcdshaperegion=(row.mRelPFIsoDBDefault>0.2) and (row.mRelPFIsoDBDefault<0.5) and (row.eIsoDB03<0.5)
          

            ## now divide by e-mu isolation regions, looseloose,loosetight,tightloose,tighttight
            isMuonTight=False
            if selections.lepton_id_iso(row, 'm', 'MuIDTight_mutauiso0p15',dataperiod=self.data_period):
                cut_flow_trk.Fill('muiso')
                isMuonTight=True

            isElecTight=False
            if selections.lepton_id_iso(row, 'e', 'eid15Loose_etauiso0p1',eIDwp='WP80'): 
                cut_flow_trk.Fill('eiso')
                isElecTight=True
 

            if not isMuonTight and not isElecTight: #double fakes, should be tiny
                region="eLoosemLoose"
            elif not isMuonTight and  isElecTight:   # mu fakes, should be small
                region="eTightmLoose"
            elif  isMuonTight and not isElecTight: #e fakes, most fakes should come from here
                region="eLoosemTight"
            elif isMuonTight and isElecTight: #signal region
                region="signal"
            
            MVAval=1    
            self.fill_histos(row,MVAval,sign,None,True,region,btagweight,'')
                

            self.var_2jet ={'mTmuMet_':row.mMtToPfMet_type1,'mTeMet_':row.eMtToPfMet_type1,'deltaPhimue_':abs(row.e_m_DPhi),'eDphiPFMet_':abs(row.eDPhiToPfMet_type1),'mDphiPFMet_':abs(row.mDPhiToPfMet_type1),'mColl_':collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi),'pZeta_':row.e_m_PZetaLess0p85PZetaVis,'mPt_':row.mPt,'ePt_':row.ePt}
                
            MVA2jet=self.functor(**self.var_2jet)

            MVAval=MVA2jet


            folder = sign+'/'+processtype+'/'+str(int(jetN))
            self.fill_histos(row,MVAval,sign,folder,False,region,btagweight,'presel',qcdshaperegion)


            for sys in self.sysdir:
                if sys =='nosys':
                    shifted_jetVeto30=row.jetVeto30
                    shifted_mDPhiToPfMet=row.mDPhiToPfMet_type1
                    shifted_mMtToPfMet=row.mMtToPfMet_type1
                    shifted_eDPhiToPfMet=row.eDPhiToPfMet_type1
                    shifted_eMtToPfMet=row.eMtToPfMet_type1
                    shifted_type1_pfMetPhi=row.type1_pfMetPhi
                    shifted_type1_pfMetEt=row.type1_pfMetEt
                    shifted_vbfMass=row.vbfMass
                    shifted_vbfDeta=row.vbfDeta
                elif sys =='jetup':
                    shifted_jetVeto30=row.jetVeto30_JetEnUp
                    shifted_mDPhiToPfMet=row.mDPhiToPfMet_JetEnUp
                    shifted_mMtToPfMet=row.mMtToPfMet_JetEnUp
                    shifted_eDPhiToPfMet=row.eDPhiToPfMet_JetEnUp
                    shifted_eMtToPfMet=row.eMtToPfMet_JetEnUp
                    shifted_type1_pfMetPhi=row.type1_pfMet_shiftedPhi_JetEnUp
                    shifted_type1_pfMetEt=row.type1_pfMet_shiftedPt_JetEnUp
                    shifted_vbfMass=row.vbfMass_JetEnUp
                    shifted_vbfDeta=row.vbfDeta_JetEnUp
                elif sys =='jetdown':
                    shifted_jetVeto30=row.jetVeto30_JetEnDown
                    shifted_mDPhiToPfMet=row.mDPhiToPfMet_JetEnDown
                    shifted_mMtToPfMet=row.mMtToPfMet_JetEnDown
                    shifted_eDPhiToPfMet=row.eDPhiToPfMet_JetEnDown
                    shifted_eMtToPfMet=row.eMtToPfMet_JetEnDown
                    shifted_type1_pfMetPhi=row.type1_pfMet_shiftedPhi_JetEnDown
                    shifted_type1_pfMetEt=row.type1_pfMet_shiftedPt_JetEnDown
                    shifted_vbfMass=row.vbfMass_JetEnDown
                    shifted_vbfDeta=row.vbfDeta_JetEnDown
                elif sys =='tup':
                    shifted_jetVeto30=row.jetVeto30
                    shifted_mDPhiToPfMet=row.mDPhiToPfMet_TauEnUp
                    shifted_mMtToPfMet=row.mMtToPfMet_TauEnUp
                    shifted_eDPhiToPfMet=row.eDPhiToPfMet_TauEnUp
                    shifted_eMtToPfMet=row.eMtToPfMet_TauEnUp
                    shifted_type1_pfMetPhi=row.type1_pfMet_shiftedPhi_TauEnUp
                    shifted_type1_pfMetEt=row.type1_pfMet_shiftedPt_TauEnUp
                    shifted_vbfMass=row.vbfMass
                    shifted_vbfDeta=row.vbfDeta
                elif sys =='tdown':
                    shifted_jetVeto30=row.jetVeto30
                    shifted_mDPhiToPfMet=row.mDPhiToPfMet_TauEnDown
                    shifted_mMtToPfMet=row.mMtToPfMet_TauEnDown
                    shifted_eDPhiToPfMet=row.eDPhiToPfMet_TauEnDown
                    shifted_eMtToPfMet=row.eMtToPfMet_TauEnDown
                    shifted_type1_pfMetPhi=row.type1_pfMet_shiftedPhi_TauEnDown
                    shifted_type1_pfMetEt=row.type1_pfMet_shiftedPt_TauEnDown
                    shifted_vbfMass=row.vbfMass
                    shifted_vbfDeta=row.vbfDeta
                elif sys =='uup':
                    shifted_jetVeto30=row.jetVeto30
                    shifted_mDPhiToPfMet=row.mDPhiToPfMet_UnclusteredEnUp
                    shifted_mMtToPfMet=row.mMtToPfMet_UnclusteredEnUp
                    shifted_eDPhiToPfMet=row.eDPhiToPfMet_UnclusteredEnUp
                    shifted_eMtToPfMet=row.eMtToPfMet_UnclusteredEnUp
                    shifted_type1_pfMetPhi=row.type1_pfMet_shiftedPhi_UnclusteredEnUp
                    shifted_type1_pfMetEt=row.type1_pfMet_shiftedPt_UnclusteredEnUp
                    shifted_vbfMass=row.vbfMass
                    shifted_vbfDeta=row.vbfDeta
                elif sys =='udown':
                    shifted_jetVeto30=row.jetVeto30
                    shifted_mDPhiToPfMet=row.mDPhiToPfMet_UnclusteredEnDown
                    shifted_mMtToPfMet=row.mMtToPfMet_UnclusteredEnDown
                    shifted_eDPhiToPfMet=row.eDPhiToPfMet_UnclusteredEnDown
                    shifted_eMtToPfMet=row.eMtToPfMet_UnclusteredEnDown
                    shifted_type1_pfMetPhi=row.type1_pfMet_shiftedPhi_UnclusteredEnDown
                    shifted_type1_pfMetEt=row.type1_pfMet_shiftedPt_UnclusteredEnDown
                    shifted_vbfMass=row.vbfMass
                    shifted_vbfDeta=row.vbfDeta

                jn = shifted_jetVeto30
                if jn > 3 : jn = 3
                if jn==2:
                    if row.vbfMass>=550:
                        jn=22
                    else:
                        jn=21
                self.var_d ={'mTmuMet_':shifted_mMtToPfMet,'mTeMet_':shifted_eMtToPfMet,'deltaPhimue_':abs(row.e_m_DPhi),'eDphiPFMet_':abs(shifted_eDPhiToPfMet),'mDphiPFMet_':abs(shifted_mDPhiToPfMet),'mColl_':collmass(row, shifted_type1_pfMetEt, shifted_type1_pfMetPhi),'pZeta_':row.e_m_PZetaLess0p85PZetaVis,'mPt_':row.mPt,'ePt_':row.ePt}
                
                MVAcutvalue=self.functor(**self.var_d)

                MVAval=MVAcutvalue

                folder = sign+'/'+processtype+'/'+str(int(jn))+'/selected/'+sys
                self.fill_histos(row,MVAval,sign,folder,False,region,btagweight,sys,qcdshaperegion)
        cut_flow_trk.flush()        
            
    def finish(self):
        self.write_histos()

