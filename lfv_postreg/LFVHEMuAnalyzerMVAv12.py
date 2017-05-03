# check in https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2015#MET when the mva met receipe is available.
from EMTree import EMTree
import os
import ROOT
import math
import glob
import array
#import mcCorrections
import baseSelections2 as selections
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


fakerateWeightGetter =FakeRate2D.make_fakerate2DVLoose()

class LFVHEMuAnalyzerMVAv12(MegaBase):
    tree = 'em/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):

        self.channel='EMu'
        super(LFVHEMuAnalyzerMVAv12, self).__init__(tree, outfile, **kwargs)
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
        self.out=outfile
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
 
        eidcorr=1
        return pu*muidcorr*muisocorr*mutrcorr*eidcorr*mutrkcorr*eisocorr0p10*etrkcorr
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
        alldirs=['','antiIsolatedweighted/','antiIsolated/','fakeRateMethod/']#,'antiIsolatedweightedUp/','antiIsolatedweightedDown/','antiIsolatedweighted_nofrweight/','subtracted/','subtractedup/','subtracteddown/','mLoose_nofrweight/','mLoose/','mLooseUp/','mLooseDown/','eLoose_nofrweight/','eLoose/','eLooseUp/','eLooseDown/','mLooseeLoose_nofrweight/','mLooseeLoose/','mLooseeLooseUp/','mLooseeLooseDown/']
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
            if k<5150:
                self.book(f,"mPt", "mu p_{T}", 200, 0, 200)
                self.book(f,"mPhi", "mu phi", 100, -3.2, 3.2)
                self.book(f,"mEta", "mu eta",  50, -2.5, 2.5)
            
                self.book(f,"ePt", "e p_{T}", 200, 0, 200)
                self.book(f,"ePhi", "e phi",  100, -3.2, 3.2)
                self.book(f,"eEta", "e eta", 50, -2.5, 2.5)
            
                self.book(f, "em_DeltaPhi", "e-mu DeltaPhi" , 50, 0, 3.2)
                self.book(f, "em_DeltaR", "e-mu DeltaR" , 100, -7, 7)
#            
                self.book(f, "h_collmass_pfmet",  "h_collmass_pfmet",  32, 0, 320)
#                self.book(f, "h_collmass_mvamet",  "h_collmass_mvamet",  32, 0, 320)
#                
#                self.book(f,"scaledmPt","scaledmPt",60,0,1.2)
#                self.book(f,"scaledePt","scaledePt",60,0,1.2)
#
                self.book(f, "Met",  "Met",  32, 0, 320)

                self.book(f, "eMtToPfMet",  "eMtToPfMet",  32, 0, 320)
                self.book(f, "mMtToPfMet",  "mMtToPfMet",  32, 0, 320)
                self.book(f, "dPhiMetToE",  "dPhiMetToE",  50, 0, 3.2)
                
                self.book(f, "h_vismass",  "h_vismass",  32, 0, 320)
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
                self.book(f, "h_collmass_pfmet",  "h_collmass_pfmet",  32, 0, 320)

        for s in sign:
            self.book(s, "jetN_30", "Number of jets, p_{T}>30", 10, -0.5, 9.5) 
            self.book(s, "NUP", "Number of Partons", 12, -0.5, 11.5) 
            self.book(s, "numGenJets", "Number of Gen Level Jets", 12, -0.5, 11.5) 
            self.book(s, "numVertices", "Number of Vertices", 60, 0, 60) 
            self.book(s, "h_collmass_pfmet", "h_collmass_pfmet", 32, 0, 320) 
            self.book(s, "h_vismass",  "h_vismass",  32, 0, 320)
            self.book(s+'/tNoCuts', "CUT_FLOW", "Cut Flow", len(cut_flow_step), 0, len(cut_flow_step))
            
            xaxis = self.histograms[s+'/tNoCuts/CUT_FLOW'].GetXaxis()
            self.cut_flow_histo = self.histograms[s+'/tNoCuts/CUT_FLOW']
            self.cut_flow_map   = {}
            for i, name in enumerate(cut_flow_step):
                xaxis.SetBinLabel(i+1, name)
                self.cut_flow_map[name] = i+0.5


    def geteFakeFactor(self,row):
        electronPt=row.ePt
        if electronPt<10:
            raise ValueError("Electron with Pt less than 10 got through!")
        if electronPt<=20:
            return 0.205
        elif electronPt<=30:
            return 0.1350
        elif electronPt<=40:
            return  0.315
        elif electronPt<=50:
            return  0.38
        elif electronPt<=70:
            return  0.49
        elif electronPt<=100:
            return  0.63
        elif electronPt<=150:
            return 0.70
        else:
            return 0.85

    def fill_jet_histos(self,row,sign,region,btagweight):

        if self.is_WJet or self.is_DYJet:
            weight1 = self.event_weight(row,region)*self.binned_weight[int(row.numGenJets)]
        else:
            weight1 = self.event_weight(row,region)

        weight1=btagweight*weight1
        self.histograms[sign+'/jetN_30'].Fill(row.jetVeto30,weight1)
        self.histograms[sign+'/NUP'].Fill(row.NUP,weight1)
        self.histograms[sign+'/numGenJets'].Fill(row.numGenJets,weight1)




    def fill_histos(self, row,sign,f=None,fillCommonOnly=False,region=None,btagweight=1,sys=''):
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
        pudir =['','fakeRateMethod/']#,'subtractedup/','subtracteddown/']
#        elooseList = ['eLoose_nofrweight/','eLoose/','eLooseUp/', 'eLooseDown/']
#        mlooseList = ['mLoose_nofrweight/','mLoose/','mLooseUp/', 'mLooseDown/']
#        emlooseList= ['mLooseeLoose_nofrweight/','mLooseeLoose/','mLooseeLooseUp/','mLooseeLooseDown/']
        alllooselist=['antiIsolatedweighted/','antiIsolated/']
#        alllooselist=['antiIsolatedweighted_nofrweight/','antiIsolatedweighted/','antiIsolatedweightedUp/','antiIsolatedweightedDown/']
#        if region!='signal':
#        fr_weights=None
#            fr_weights = [1,fakerate_weight,fakerate_weight+0.3*fakerate_weight,fakerate_weight-0.3*fakerate_weight]
#  elif region=='eTightmLoose':
#      fr_weights=[1,0.7,0.91,0.59]
#  elif region=='eLoosemLoose': 
#      frweight_array=fakerateWeight(row.ePt,abs(row.eEta))
#      efr_weights = [1,frweight_array[0],frweight_array[0]+frweight_array[1],frweight_array[0]-frweight_array[1]]
#      mfr_weights=[1,0.7,0.91,0.59]
#      comb_error=efr_weights[1]*mfr_weights[1]*sqrt(0.3*0.3+((frweight_array[1]/frweight_array[0])*(frweight_array[1]/frweight_array[0])))
#      fr_weights = [1,efr_weights[1]*mfr_weights[1],efr_weights[1]*mfr_weights[1]+comb_error,efr_weights[1]*mfr_weights[1]-comb_error]
        if region=='eLoosemTight':
            frarray=fakerateWeightGetter(row.ePt,abs(row.eEta)) 
            fakerateWeight=frarray[0]
#            fakerateFactor=self.geteFakeFactor(row)
           # fakerateWeight=fakerateFactor/(1-fakerateFactor)
            antiIsolatedWeightList=[fakerateWeight,1.0]
            for n, l in enumerate(alllooselist) :
                antiIsolatedWeight= weight*antiIsolatedWeightList[n]
                folder = l+f
                if sys=='presel' or sys=='nosys':
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
                elif sys=='jetup':
                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMet_shiftedPt_JetEnUp, row.type1_pfMet_shiftedPhi_JetEnUp),antiIsolatedWeight)
                elif sys=='jetdown':
                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMet_shiftedPt_JetEnDown, row.type1_pfMet_shiftedPhi_JetEnDown),antiIsolatedWeight)
                elif sys=='tup':
                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMet_shiftedPt_TauEnUp, row.type1_pfMet_shiftedPhi_TauEnUp),antiIsolatedWeight)
                elif sys=='tdown':
                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMet_shiftedPt_TauEnDown, row.type1_pfMet_shiftedPhi_TauEnDown),antiIsolatedWeight)
                elif sys=='uup':
                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMet_shiftedPt_UnclusteredEnUp, row.type1_pfMet_shiftedPhi_UnclusteredEnUp),antiIsolatedWeight)
                elif sys=='udown':
                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMet_shiftedPt_UnclusteredEnDown, row.type1_pfMet_shiftedPhi_UnclusteredEnDown),antiIsolatedWeight)
                    """
        if region=='eTightmLoose' and (sys=='presel' or sys=='nosys'):
            fr_weights=[1,0.7,0.91,0.59]
            if not self.isData :
                for n in range(1,4) :
                    subfolder=pudir[n]+f
                    subweight=fr_weights[n]*weight*(-1)
                    histos[subfolder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi),subweight)
                    histos[folder+'/mPt'].Fill(row.mPt, subweight)
                    histos[folder+'/mEta'].Fill(row.mEta, subweight)
                    histos[folder+'/mPhi'].Fill(row.mPhi, subweight) 
                    histos[folder+'/ePt'].Fill(row.ePt, subweight)
                    histos[folder+'/eEta'].Fill(row.eEta, subweight)
                    histos[folder+'/ePhi'].Fill(row.ePhi, subweight)
                    histos[folder+'/em_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.mPhi), subweight)
                    histos[folder+'/em_DeltaR'].Fill(row.e_m_DR, subweight)
                    histos[folder+'/h_vismass'].Fill(row.e_m_Mass, subweight)
                    histos[folder+'/ePFMET_Mt'].Fill(row.eMtToPfMet_type1, subweight)
                    histos[folder+'/mPFMET_Mt'].Fill(row.mMtToPfMet_type1, subweight)
                    histos[folder+'/ePFMET_DeltaPhi'].Fill(abs(row.eDPhiToPfMet_type1), subweight)
                    histos[folder+'/mPFMET_DeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1), subweight)
                    histos[folder+'/mPFMETDeltaPhi_vs_ePFMETDeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1),abs(row.eDPhiToPfMet_type1) , subweight)
                    histos[folder+'/vbfMass'].Fill(row.vbfMass, subweight)
                    histos[folder+'/vbfDeta'].Fill(row.vbfDeta, subweight)
                    histos[folder+'/jetN_30'].Fill(row.jetVeto30, subweight) 

            for n, l in enumerate(mlooseList) :
                frweight = weight*fr_weights[n]
                folder = l+f
                histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi),frweight)
        elif region=='eLoosemTight' and (sys=='presel' or sys=='nosys'):
            frweight_array=fakerateWeight(row.ePt,abs(row.eEta))
            fr_weights = [1,frweight_array[0],frweight_array[0]+frweight_array[1],frweight_array[0]-frweight_array[1]]"""
            if not self.isData and not self.isGluGlu_LFV and not self.isVBF_LFV and not self.isGluGluEtauSig and not self.isVBFEtauSig:
                fakeRateMethodfolder=pudir[1]+f
                fakeRateMethodweight=fakerateWeight*weight*(-1)
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
        """
            for n, l in enumerate(elooseList) :
                frweight = weight*fr_weights[n]
                folder = l+f
                histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi),frweight)

        elif region=='eLoosemLoose' and (sys=='presel' or sys=='nosys'):
            frweight_array=fakerateWeight(row.ePt,abs(row.eEta))
            efr_weights = [1,frweight_array[0],frweight_array[0]+frweight_array[1],frweight_array[0]-frweight_array[1]]
            mfr_weights=[1,0.7,0.91,0.59]
            comb_error=efr_weights[1]*mfr_weights[1]*sqrt(0.3*0.3+((frweight_array[1]/frweight_array[0])*(frweight_array[1]/frweight_array[0])))
            fr_weights = [1,efr_weights[1]*mfr_weights[1],efr_weights[1]*mfr_weights[1]+comb_error,efr_weights[1]*mfr_weights[1]-comb_error]
            if not self.isData and fr_weights[1]==123123:
                for n in range(1,4):
                    subfolder=pudir[n]+f
                    subweight=fr_weights[n]*weight*(-1)
                    histos[subfolder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi),subweight)
                    histos[folder+'/mPt'].Fill(row.mPt, subweight)
                    histos[folder+'/mEta'].Fill(row.mEta, subweight)
                    histos[folder+'/mPhi'].Fill(row.mPhi, subweight) 
                    histos[folder+'/ePt'].Fill(row.ePt, subweight)
                    histos[folder+'/eEta'].Fill(row.eEta, subweight)
                    histos[folder+'/ePhi'].Fill(row.ePhi, subweight)
                    histos[folder+'/em_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.mPhi), subweight)
                    histos[folder+'/em_DeltaR'].Fill(row.e_m_DR, subweight)
                    histos[folder+'/h_vismass'].Fill(row.e_m_Mass, subweight)
                    histos[folder+'/ePFMET_Mt'].Fill(row.eMtToPfMet_type1, subweight)
                    histos[folder+'/mPFMET_Mt'].Fill(row.mMtToPfMet_type1, subweight)
                    histos[folder+'/ePFMET_DeltaPhi'].Fill(abs(row.eDPhiToPfMet_type1), subweight)
                    histos[folder+'/mPFMET_DeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1), subweight)
                    histos[folder+'/mPFMETDeltaPhi_vs_ePFMETDeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1),abs(row.eDPhiToPfMet_type1) , subweight)
                    histos[folder+'/vbfMass'].Fill(row.vbfMass, subweight)
                    histos[folder+'/vbfDeta'].Fill(row.vbfDeta, subweight)
                    histos[folder+'/jetN_30'].Fill(row.jetVeto30, subweight) 

            for n, l in enumerate(emlooseList) :
                frweight = weight*fr_weights[n]
                folder = l+f
                histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi),frweight)
"""          
        if region=='signal' :
            for n,d  in enumerate(pudir) :
                folder = d+f                
                if sys=='presel' or sys=='nosys':
                    histos[folder+'/mPt'].Fill(row.mPt, weight)
                    histos[folder+'/mEta'].Fill(row.mEta, weight)
                    histos[folder+'/mPhi'].Fill(row.mPhi, weight) 
                    histos[folder+'/ePt'].Fill(row.ePt, weight)
                    histos[folder+'/eEta'].Fill(row.eEta, weight)
                    histos[folder+'/ePhi'].Fill(row.ePhi, weight)
                    histos[folder+'/em_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.mPhi), weight)
                    histos[folder+'/em_DeltaR'].Fill(row.e_m_DR, weight)
                    histos[folder+'/h_vismass'].Fill(row.e_m_Mass, weight)

                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi),weight)
 #                   histos[folder+'/scaledmPt'].Fill(float(row.mPt)/float(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)),weight)
 #                   histos[folder+'/scaledePt'].Fill(float(row.ePt)/float(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)),weight)
 #
                    histos[folder+'/Met'].Fill(row.type1_pfMetEt, weight)

                    histos[folder+'/ePFMET_Mt'].Fill(row.eMtToPfMet_type1, weight)
                    histos[folder+'/mPFMET_Mt'].Fill(row.mMtToPfMet_type1, weight)
                    histos[folder+'/ePFMET_DeltaPhi'].Fill(abs(row.eDPhiToPfMet_type1), weight)
                    histos[folder+'/mPFMET_DeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1), weight)
#                    histos[folder+'/mPFMETDeltaPhi_vs_ePFMETDeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1),abs(row.eDPhiToPfMet_type1) , weight)
 
                    histos[folder+'/vbfMass'].Fill(row.vbfMass, weight)
                    histos[folder+'/vbfDeta'].Fill(row.vbfDeta, weight)
                    histos[folder+'/jetN_30'].Fill(row.jetVeto30, weight) 
 #                   

                elif sys=='jetup':
                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMet_shiftedPt_JetEnUp, row.type1_pfMet_shiftedPhi_JetEnUp),weight)

                
                elif sys=='jetdown':
                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMet_shiftedPt_JetEnDown, row.type1_pfMet_shiftedPhi_JetEnDown),weight)


                elif sys=='tup':
                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMet_shiftedPt_TauEnUp, row.type1_pfMet_shiftedPhi_TauEnUp),weight)


                elif sys=='tdown':
                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMet_shiftedPt_TauEnDown, row.type1_pfMet_shiftedPhi_TauEnDown),weight)


                elif sys=='uup':
                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMet_shiftedPt_UnclusteredEnUp, row.type1_pfMet_shiftedPhi_UnclusteredEnUp),weight)

                elif sys=='udown':
                    folder = d+f
                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMet_shiftedPt_UnclusteredEnDown, row.type1_pfMet_shiftedPhi_UnclusteredEnDown),weight)

    def process(self):
        
        cut_flow_histo = self.cut_flow_histo
        cut_flow_trk   = cut_flow_tracker(cut_flow_histo)
        myevent=()
        frw = []
        curr_event=0
        for row in self.tree:
            sign = 'ss' if row.e_m_SS else 'os'

 #           ptthreshold = [30]
            repeatEvt=True
            if row.evt!=curr_event:
                curr_event=row.evt
                repeatEvt=False
            
            if repeatEvt:continue

#            print "non-repeat"
            processtype ='gg'##changed from 20

            cut_flow_trk.new_row(row.run,row.lumi,row.evt)
           
            cut_flow_trk.Fill('allEvents')

            #trigger
            if (not bool(row.singleIsoMu22Pass or row.singleIsoTkMu22Pass)) and self.isData: 
                continue   #notrigger in new MC; add later

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
            if not selections.lepton_id_iso(row, 'm', 'MuIDTight_idiso0p25'): continue
            cut_flow_trk.Fill('mulooseiso')


            #E Preselection
            if not selections.eSelection(row, 'e'): continue
            cut_flow_trk.Fill('esel')
           
            if not selections.lepton_id_iso(row, 'e', 'eid15Loose_etauiso0p5','WP80'): continue
            cut_flow_trk.Fill('elooseiso')

            #take care of ecal gap
            if row.eAbsEta > 1.4442 and row.eAbsEta < 1.566 : continue             

            cut_flow_trk.Fill('ecalgap')

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

            cut_flow_trk.Fill('bjetveto')
            ## All preselection passed


            ## now divide by e-mu isolation regions, looseloose,loosetight,tightloose,tighttight
            isMuonTight=False
            if selections.lepton_id_iso(row, 'm', 'MuIDTight_mutauiso0p15'):
                cut_flow_trk.Fill('muiso')
                isMuonTight=True

            isElecTight=False
            if selections.lepton_id_iso(row, 'e', 'eid15Loose_etauiso0p1','WP80'): 
                cut_flow_trk.Fill('eiso')
                isElecTight=True
 

            if not isMuonTight and not isElecTight: #double fakes, should be tiny
                region="eLoosemLoose"
                continue
            elif not isMuonTight and  isElecTight:   # mu fakes, should be small
                region="eTightmLoose"
                continue
            elif  isMuonTight and not isElecTight: #e fakes, most fakes should come from here
                region="eLoosemTight"
            elif isMuonTight and isElecTight: #signal region
                region="signal"

            self.fill_histos(row,sign,None,True,region,btagweight,'')
                
            jetN=row.jetVeto30
            if jetN>3:
                jetN=3
            if jetN==2:
                if row.vbfMass>=550:
                    jetN=22
                else:
                    jetN=21


            folder = sign+'/'+processtype+'/'+str(int(jetN))
            self.fill_histos(row,sign,folder,False,region,btagweight,'presel')


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
                if jn == 0 :
                    if row.mPt < 30: continue 
                    if row.ePt < 10: continue
                    if deltaPhi(row.ePhi,row.mPhi) < 2.5 : continue
                    if abs(shifted_eDPhiToPfMet) > 0.7 : continue
                    if shifted_mMtToPfMet < 60 : continue
 #                   if shifted_eMtToPfMet > 65 : continue
                    cut_flow_trk.Fill('jet0sel')
                
                if jn == 1 :
                    if row.mPt < 25: continue 
                    if row.ePt < 10 : continue
                    if abs(shifted_eDPhiToPfMet) > 0.7 : continue
                    if deltaPhi(row.ePhi,row.mPhi)<1:continue
                    if shifted_mMtToPfMet < 40 : continue
#                    if shifted_eMtToPfMet > 65 : continue
                    cut_flow_trk.Fill('jet1sel')
                    
                if jn == 21 :
                    if row.mPt < 25: continue 
                    if row.ePt < 10 : continue  #no cut as only electrons with pt>30 are in the ntuples
                    if abs(shifted_eDPhiToPfMet) > 0.5 : continue
                    if shifted_mMtToPfMet < 15 : continue
#                    if shifted_eMtToPfMet > 15 : continue
                    if shifted_vbfMass < 60 : continue
 #                   if shifted_vbfDeta < 0.5 : continue
                    cut_flow_trk.Fill('jet2loosesel')

                if jn == 22 :
                    if row.mPt < 25: continue 
                    if row.ePt < 10 : continue  #no cut as only electrons with pt>30 are in the ntuples
                    if abs(shifted_eDPhiToPfMet) > 0.3 : continue
                    if shifted_mMtToPfMet < 15 : continue
#                    if shifted_eMtToPfMet > 15 : continue
                    if shifted_vbfMass < 60 : continue
#                    if shifted_vbfDeta < 0.5 : continue
                    cut_flow_trk.Fill('jet2tightsel')


                folder = sign+'/'+processtype+'/'+str(int(jn))+'/selected/'+sys
                self.fill_histos(row,sign,folder,False,region,btagweight,sys)
        cut_flow_trk.flush()        
            
    def finish(self):
        self.write_histos()


