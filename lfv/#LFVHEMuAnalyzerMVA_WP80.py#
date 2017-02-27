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
import FinalStateAnalysis.TagAndProbe.FakeRate2D as FakeRate2D
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

pu_distributions = glob.glob(os.path.join('inputs', os.environ['jobid'], 'data_SingleMu*pu.root'))


pu_corrector = PileupWeight.PileupWeight('MC_Spring16', *pu_distributions)
id_corrector  = MuonPOGCorrections.make_muon_pog_PFMedium_2016BCD()
iso_corrector = MuonPOGCorrections.make_muon_pog_TightIso_2016BCD()
tr_corrector  = MuonPOGCorrections.make_muon_pog_IsoMu22oIsoTkMu22_2016BCD()
trk_corrector = MuonPOGCorrections.mu_trackingEta_2016
#trk_corrector =  MuonPOGCorrections.make_muonptabove10_pog_tracking_corrections_2016()
eId_corrector = EGammaPOGCorrections.make_egamma_pog_electronID_ICHEP2016( 'nontrigWP80')

fakerateWeight =FakeRate2D.make_fakerate2D()
class LFVHEMuAnalyzerMVA_WP80(object):
    tree = 'em/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel='EMu'
#        super(LFVHEMuAnalyzerMVA_WP80, self).__init__(tree, outfile, **kwargs)
        target = os.path.basename(os.environ['megatarget'])
        print outfile
        self.target=target
        self.output=outfile
        self.is_WJet=('WJetsToLNu' in target or 'W1JetsToLNu' in target or 'W2JetsToLNu' in target or 'W3JetsToLNu' in target or 'W4JetsToLNu' in target)
        self.is_DYJet= ('DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target or  'DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target or 'DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target or 'DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target or 'DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target) 
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
        self.isGluGlu_LFV=('GluGlu_LFV' in target)
        self.isVBFHTo=('VBFHTo' in target)
        self.isVBF_LFV=('VBF_LFV' in target)

        self.WGToLNuG_weight=0.00011736233
        self.WGstarToLNuEE_weight=0.00000564158
        self.WGstarToLNuMuMu_weight=0.00000180887
        self.ST_tW_antitop_weight=0.0000526974897638
        self.ST_tW_top_weight=0.000040713632205
        self.WW_weight=0.000132802570574
        self.WZ_weight=4.713e-05
        self.ZZ_weight=1.67015056929e-05
        self.TT_weight=8.71008159645e-06
        self.GluGluHTo_weight= 2.04805444356e-06
        self.GluGlu_LFV_weight=5.12622957741e-06
        self.VBFHTo_weight=4.67869114625e-08
        self.VBF_LFV_weight=7.57282991971e-08  
        self.tree = EMTree(tree)
        self.out=outfile
        self.histograms = {}
        self.mym1='m'
        self.mye1='e'
        self.sysdir=['nosys','jetup','jetdown','tup','tdown','uup','udown']
#        self.sysdir=['nosys']
        if self.is_WJet:
            self.binned_weight=[0.672454854,0.20530173,0.10758429,0.08283444,0.090256696]
        elif self.is_DYJet:
            self.binned_weight=[0.063407117,0.014028728,0.015010691,0.01553876,0.012480257]
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

    def mc_corrector_2015(self, row):
        
        pu = pu_corrector(row.nTruePU)
 #       pu=1
        muidcorr = id_corrector(getattr(row, self.mym1+'Pt'), abs(getattr(row, self.mym1+'Eta')))
        muisocorr = iso_corrector('Tight', getattr(row, self.mym1+'Pt'), abs(getattr(row, self.mym1+'Eta')))
#        print "id corr", muidcorr
 #       print "iso corr", muisocorr
        mutrcorr = tr_corrector(getattr(row, self.mym1+'Pt'), abs(getattr(row, self.mym1+'Eta'))) 
        eidcorr = eId_corrector(getattr(row,self.mye1+'Eta'),getattr(row, self.mye1+'Pt'))
        mutrkcorr=trk_corrector(getattr(row,self.mym1+'Eta'))[0]
#        mutrkcorr=trk_corrector(getattr(row,self.mym1+'Eta'))
 #       print "trk corr",mutrkcorr
  #      print "tr corr", mutrcorr
   #     print "eid corr", eidcorr
#       mutrcorr=1
     # if pu*muidcorr1*muisocorr1*muidcorr2*muisocorr2*mutrcorr==0: print pu, muidcorr1, muisocorr1, muidcorr2, muisocorr2, mutrcorr
    #    print "pileup--------   =",pu
   #     print pu*muidcorr*muisocorr*mutrcorr
        return pu*muidcorr*muisocorr*mutrcorr*eidcorr*mutrkcorr
       # return pu*muidcorr*mutrcorr*eidcorr


    def correction(self,row):
	return self.mc_corrector_2015(row)
        
    def event_weight(self, row):
 
        if row.run > 2: #FIXME! add tight ID correction
            return 1.
       # if row.GenWeight*self.correction(row) == 0 : print 'weight==0', row.GenWeight*self.correction(row), row.GenWeight, self.correction(row), row.m1Pt, row.m2Pt, row.m1Eta, row.m2Eta
       # print row.GenWeight, "lkdfh"
        return row.GenWeight*self.correction(row) 
#        return self.correction(row) 




    def begin(self):
        f=''
        self.weight_=array.array( 'f', [ 0 ] )
        self.treeS=ROOT.TTree("treeS","treeS")
        self.treeS.Branch("weight_",self.weight_,"weight_/F")


    def fill_histos(self, row, f=None,region=None,btagweight=1,sys=''):
        if self.is_WJet or self.is_DYJet:
            weight = self.event_weight(row)*self.binned_weight[int(row.numGenJets)]*0.001
        elif self.isWGToLNuG:
            weight=self.WGToLNuG_weight*self.event_weight(row)
        elif self.isWGstarToLNuEE:
            weight=self.WGstarToLNuEE_weight*self.event_weight(row)
        elif self.isWGstarToLNuMuMu:
            weight=self.WGstarToLNuMuMu_weight*self.event_weight(row)
        elif self.isST_tW_top:
            weight=self.ST_tW_top_weight*self.event_weight(row)
        elif self.isST_tW_antitop:
            weight=self.ST_tW_antitop_weight*self.event_weight(row)
        elif self.isWW:
            weight=self.WW_weight*self.event_weight(row)
        elif self.isWZ:
            weight=self.WZ_weight*self.event_weight(row)
        elif self.isZZ:
            weight=self.ZZ_weight*self.event_weight(row)
        elif self.isTT:
            weight=self.TT_weight*self.event_weight(row)
        elif self.isGluGluHTo:
            weight=self.GluGluHTo_weight*self.event_weight(row)
        elif self.isGluGlu_LFV:
            weight=self.GluGlu_LFV_weight*self.event_weight(row)
        elif self.isVBFHTo:
            weight=self.VBFHTo_weight*self.event_weight(row)
        elif self.isVBF_LFV:
            weight=self.VBF_LFV_weight*self.event_weight(row)
        else:
            weight = self.event_weight(row)
        weight=btagweight*weight
        self.weight_[0]=weight
        self.treeS.Fill()

        

    def process(self):
        myevent=()
        frw = []
        for row in self.tree:
            sign = 'ss' if row.e_m_SS else 'os'
#            print "row",row
 #           ptthreshold = [30]
            processtype ='gg'##changed from 20


            #trigger
            if (not bool(row.singleIsoMu22Pass or row.singleIsoTkMu22Pass)) and self.isData: 
                continue   #notrigger in new MC; add later

            #vetoes and cleaning
#            if row.bjetCISVVeto30Loose : continue
 #           cut_flow_trk.Fill('bjetveto')
            nbtagged=row.bjetCISVVeto30Medium
            if nbtagged>2:
                nbtagged=2
            btagweight=1
            if (self.isData and nbtagged>0):
                continue
            if nbtagged>0:
                if nbtagged==1:
                    btagweight=bTagSF.bTagEventWeight(nbtagged,row.jb1pt,row.jb1flavor,row.jb2pt,row.jb2flavor,1,0,0) if (row.jb1pt>-990 and row.jb1flavor>-990) else 0
                if nbtagged==2:
                    btagweight=bTagSF.bTagEventWeight(nbtagged,row.jb1pt,row.jb1flavor,row.jb2pt,row.jb2flavor,1,0,0) if (row.jb1pt>-990 and row.jb1flavor>-990 and row.jb2pt>-990 and row.jb2flavor>-990) else 0
#                print "btagweight,nbtagged,row.jb1pt,row.jb1flavor,row.jb2pt,row.jb2flavor"," ",btagweight," ",nbtagged," ",row.jb1pt," ",row.jb1flavor," ",row.jb2pt," ",row.jb2flavor
            if btagweight==0: continue

            if deltaR(row.ePhi,row.mPhi,row.eEta,row.mEta)<0.1:continue
            
            if row.muVetoPt5IsoIdVtx :continue

            if row.eVetoMVAIso :continue


            #E Preselection
            if not selections.eSelection(row, 'e'): continue

           
            if not selections.lepton_id_iso(row, 'e', 'eid15Loose_etauiso1'): continue


            #take care of ecal gap
            if row.eEta > 1.4442 and row.eEta < 1.566 : continue             

            #mu preselection
            if not selections.muSelection(row, 'm'): continue

            if not selections.lepton_id_iso(row, 'm', 'MuIDTight_idiso025'): continue


            ## All preselection passed

            ## now divide by e-mu isolation regions, looseloose,loosetight,tightloose,tighttight
            isMuonTight=False
            if selections.lepton_id_iso(row, 'm', 'MuIDTight_mutauiso015'):isMuonTight=True

            isElecTight=False
            if selections.lepton_id_iso(row, 'e', 'eid15Loose_etauiso01'): isElecTight=True
 

            if not isMuonTight and not isElecTight: #double fakes, should be tiny
                region="eLoosemLoose"
            elif not isMuonTight and  isElecTight:   # mu fakes, should be small
                region="eTightmLoose"
            elif  isMuonTight and not isElecTight: #e fakes, most fakes should come from here
                region="eLoosemTight"
            elif isMuonTight and isElecTight: #signal region
                region="signal"
            if self.is_WJet or self.is_DYJet:
                self.fill_jet_histos(row,sign)
                
            jetN=row.jetVeto30
            if jetN>3:
                jetN=3
            if jetN==2:
                if row.vbfMass>=550:
                    jetN=22
                else:
                    jetN=21

            folder = sign+'/'+processtype+'/'+str(int(jetN))
            self.fill_histos(row,folder,region,btagweight,'presel')


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
                    if deltaPhi(row.ePhi, row.mPhi) < 2.5 : continue
                    if abs(shifted_eDPhiToPfMet) > 0.7 : continue
                    if shifted_mMtToPfMet < 60 : continue
 #                   if shifted_eMtToPfMet > 65 : continue

                
                if jn == 1 :
                    if row.mPt < 25: continue 
                    if row.ePt < 10 : continue
                    if abs(shifted_eDPhiToPfMet) > 0.7 : continue
                    if deltaR(row.ePhi,row.mPhi,row.eEta,row.mEta)<1:continue
                    if shifted_mMtToPfMet < 40 : continue
#                    if shifted_eMtToPfMet > 65 : continue
                    
                if jn == 21 :
                    if row.mPt < 25: continue 
                    if row.ePt < 10 : continue  #no cut as only electrons with pt>30 are in the ntuples
                    if abs(shifted_eDPhiToPfMet) > 0.3 : continue
                    if shifted_mMtToPfMet < 15 : continue
#                    if shifted_eMtToPfMet > 15 : continue
                    if shifted_vbfMass < 50 : continue
                    if shifted_vbfDeta < 0.5 : continue

                if jn == 22 :
                    if row.mPt < 25: continue 
                    if row.ePt < 10 : continue  #no cut as only electrons with pt>30 are in the ntuples
                    if abs(shifted_eDPhiToPfMet) > 0.3 : continue
                    if shifted_mMtToPfMet < 15 : continue
#                    if shifted_eMtToPfMet > 15 : continue
                    if shifted_vbfMass < 50 : continue
                    if shifted_vbfDeta < 0.5 : continue


                folder = sign+'/'+processtype+'/'+str(int(jn))+'/selected/'+sys
#                self.fill_histos(row, folder, region,btagweight,sys)
            
    def finish(self):
        self.output.Write()


