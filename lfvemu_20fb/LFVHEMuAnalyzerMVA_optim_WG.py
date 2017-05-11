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
from inspect import currentframe

import optimizer as optimizer
import bTagSF as bTagSF

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
id_corrector  = MuonPOGCorrections.make_muon_pog_PFTight_2016BCD()
iso_corrector = MuonPOGCorrections.make_muon_pog_LooseIso_2016BCD()
tr_corrector  = MuonPOGCorrections.make_muon_pog_IsoMu22oIsoTkMu22_2016BCD()
trk_corrector = MuonPOGCorrections.mu_trackingEta_2016
#trk_corrector =  MuonPOGCorrections.make_muonptabove10_pog_tracking_corrections_2016()
eId_corrector = EGammaPOGCorrections.make_egamma_pog_electronID_ICHEP2016( 'nontrigWP80')

class LFVHEMuAnalyzerMVA_optim_WG(MegaBase):
    tree = 'em/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel='EMu'
        super(LFVHEMuAnalyzerMVA_optim_WG, self).__init__(tree, outfile, **kwargs)
        target = os.path.basename(os.environ['megatarget'])
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

#        self.dataC_weight=0.0003779335
#        self.dataD_weight=0.00022970294
#        self.dataB_weight=0.00017008874
#        self.WGToLNuG_weight=0.00011736233
#        self.WGstarToLNuEE_weight=0.00000564158
#        self.WGstarToLNuMuMu_weight=0.00000180887
#        self.ST_tW_antitop_weight=0.00006139587
#        self.ST_tW_top_weight=0.00012588401
#        self.WW_weight=0.00006415773
#        self.WZ_weight=0.0000405271
#        self.ZZ_weight=0.00001943092
#        self.TT_weight=0.00008897922
#        self.GluGluHTo_weight=0.00000185197
#        self.GluGlu_LFV_weight=0.00000235941
#        self.VBFHTo_weight=0.000000158182746
#        self.VBF_LFV_weight=0.00000024719694
#
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

        #self.pucorrector = mcCorrections.make_puCorrector('singlee')
        #self.pucorrectorUp = mcCorrections.make_puCorrectorUp('singlee')
        #self.pucorrectorDown = mcCorrections.make_puCorrectorDown('singlee')
     
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




    def begin(self):
        cuts={}
        cuts[0]=optimizer.compute_regions_0jet(100000,100000,100000,1000,-1000000,-100000,100000,100000,-10000)+['selected']
        cuts[1]=optimizer.compute_regions_1jet(100000,100000,100000,1000,-1000000,-100000,100000,1000000,-10000)+['selected']
        cuts[21]=optimizer.compute_regions_2jetgg(100000,100000,100000,-1000,-1000000,109000,100000,1000000,100000,-10000)+['selected']
        cuts[22]=optimizer.compute_regions_2jetvbf(100000,100000,100000,-1000,-1000000,109000,100000,1000000,100000,-10000)+['selected']

        cuts[3]=['selected']
        processtype=['gg']
        sign=[ 'ss','os']
        jetN = [0, 1, 21,22,3]
        folder=[]
        pudir=['','mLoose/']

        print cuts
        for d  in pudir :
            for i in sign:
                for j in processtype:
                    for jn in jetN: 
                        folder.append(d+i+'/'+j+'/'+str(jn))
                        for k in cuts[jn]:
                            for s in self.sysdir:
                                folder.append(d+i+'/'+j+'/'+str(jn)+'/'+k+'/'+s)
                            

        for f in folder:
            self.book(f, "h_collmass_pfmet",  "h_collmass_pfmet",  32, 0, 320)

        for s in sign:
            self.book(s+'/tNoCuts', "CUT_FLOW", "Cut Flow", len(cut_flow_step), 0, len(cut_flow_step))
            xaxis = self.histograms[s+'/tNoCuts/CUT_FLOW'].GetXaxis()
            self.cut_flow_histo = self.histograms[s+'/tNoCuts/CUT_FLOW']
            self.cut_flow_map   = {}
            for i, name in enumerate(cut_flow_step):
                xaxis.SetBinLabel(i+1, name)
                self.cut_flow_map[name] = i+0.5




    def fakerate_weights(self, tEta, central_weights, p1s_weights, m1s_weights):
        frweight=[1.,1.,1.]

        #central_weights = fakerate_central_histogram(25,0, 2.5)
        #p1s_weights = fakerate_central_histogram(25,0, 2.5)
        #m1s_weights = fakerate_central_histogram(25,0, 2.5)

        for n,w in enumerate( central_weights ):
            if abs(tEta) < w[1]:
                break
            ##frweight[0] = w[0]
            ##frweight[1] = p1s_weights[n][0]
            ##frweight[2] = m1s_weights[n][0]
            freight[0] = 1.
            freight[1] = 1.
            freight[2] = 1.
            
        
        return  frweight;

    
                    
    def fill_histos(self, row, f='os/gg/ept0/0',  isSignalRegion=False, frw=[1.,1.,1.],btagweight=1,sys=''):
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
        weight=weight*btagweight
        histos = self.histograms
        pudir =['']
        looseList = ['mLoose/']
        if not isSignalRegion :
            if not True:
                frweight_bv = 1.
                err = frweight_bv*0.05
                frweight_p1s = frweight_bv*(1+err)
                frweight_m1s = frweight_bv*(1-err)
        
                fr_weights = [frweight_bv, frweight_p1s, frweight_m1s]
        
                for n, l in enumerate(looseList) :
                    frweight = weight[0]*fr_weights[n]
                    folder = l+f
                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.raw_pfMetEt, row.raw_pfMetPhi), frweight)
                    histos[folder+'/h_vismass'].Fill(row.e_m_Mass, frweight)
        else:
            for n,d  in enumerate(pudir) :
                folder = d+f
                if sys=='presel' or sys=='nosys':
                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi),weight)

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
                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMet_shiftedPt_UnclusteredEnDown, row.type1_pfMet_shiftedPhi_UnclusteredEnDown),weight)

    def process(self):
        
        cut_flow_histo = self.cut_flow_histo
        cut_flow_trk   = cut_flow_tracker(cut_flow_histo)
        myevent=()


        frw = []

        for row in self.tree:
            sign = 'ss' if row.e_m_SS else 'os'

 #           ptthreshold = [30]
            processtype ='gg'##changed from 20

            cut_flow_trk.new_row(row.run,row.lumi,row.evt)
           
            cut_flow_trk.Fill('allEvents')
            if (not bool(row.singleIsoMu22Pass or row.singleIsoTkMu22Pass)) and self.isData: 
                continue   #notrigger in new MC; add later

#            if not bool(row.singleIsoMu20Pass) : continue   #notrigger in new MC; add later
            cut_flow_trk.Fill('HLTIsoPasstrg')


#            jn = row.jetVeto30
            #print 'number of jets', jn
 #           if jn > 3 : jn = 3

            #take care of ecal gap
            if row.eEta > 1.4442 and row.eEta < 1.566 : continue             


            if not selections.muSelection(row, 'm'): continue
            cut_flow_trk.Fill('musel')
            if not selections.lepton_id_iso(row, 'm', 'MuIDTight_idiso025'): continue
            cut_flow_trk.Fill('muiso')

            frw=1. ## add the correct fakerate weight once we have it


            

            if row.bjetCISVVeto30Medium : continue
            btagweight=1
 #           cut_flow_trk.Fill('bjetveto')

#            nbtagged=row.bjetCISVVeto30Medium
#            if nbtagged>2:
#                nbtagged=2
#            btagweight=1
#            if (self.isData and nbtagged>0):
#                continue
#            if nbtagged>0:
#                if nbtagged==1:
#                    btagweight=bTagSF.bTagEventWeight(nbtagged,row.jb1pt,row.jb1flavor,row.jb2pt,row.jb2flavor,1,0,0) if (row.jb1pt>-990 and row.jb1flavor>-990) else 0
#                if nbtagged==2:
#                    btagweight=bTagSF.bTagEventWeight(nbtagged,row.jb1pt,row.jb1flavor,row.jb2pt,row.jb2flavor,1,0,0) if (row.jb1pt>-990 and row.jb1flavor>-990 and row.jb2pt>-990 and row.jb2flavor>-990) else 0
##                print "btagweight,nbtagged,row.jb1pt,row.jb1flavor,row.jb2pt,row.jb2flavor"," ",btagweight," ",nbtagged," ",row.jb1pt," ",row.jb1flavor," ",row.jb2pt," ",row.jb2flavor
#            if btagweight<0:btagweight=0
#            
#            if btagweight==0: continue
#
#


            if deltaR(row.ePhi,row.mPhi,row.eEta,row.mEta)<0.1:continue
            cut_flow_trk.Fill('DR_e_mu')
  

            if row.muVetoPt5IsoIdVtx :continue
            cut_flow_trk.Fill('surplus_mu_veto')

            if row.eVetoMVAIso :continue
            cut_flow_trk.Fill('surplus_e_veto')

            #print "all else"
           
            #e Preselection
            if not selections.eSelection(row, 'e'): continue
            cut_flow_trk.Fill('esel')
            #print "ele sel--------------------------"
            
           
            if not selections.lepton_id_iso(row, 'e', 'eid15Loose_etauiso1'): continue
            cut_flow_trk.Fill('eiso')


            isMuonTight=False
            if selections.lepton_id_iso(row, 'm', 'MuIDTight_mutauiso015'):isMuonTight=True

            isElecTight=False
            if selections.lepton_id_iso(row, 'e', 'eid15Loose_etauiso01'): isElecTight=True


            isSignalRegion=False
            if isMuonTight and isElecTight:
                isSignalRegion=True
            

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

#                jn = shifted_jetVeto30
                jn = row.jetVeto30

                if jn > 3 : jn = 3

                if jn==2:
                    if row.vbfMass>=500:
                        jn=22
                    else:
                        jn=21

                selections_passed=[]
                for dummy in [1]:
                    if jn == 0 :
 #                       if row.mPt < 30: break  
                        if row.ePt < 10 : break
                        if abs(shifted_eDPhiToPfMet) > 0.7 : break
                        if deltaPhi(row.ePhi, row.mPhi) < 2.5 : break
                        if shifted_mMtToPfMet < 60 : break
#                        if shifted_eMtToPfMet > 65 : break
                        selections_passed.extend([('0',i) for i in optimizer.compute_regions_0jet(row.mPt, row.ePt,deltaPhi(row.ePhi, row.mPhi),shifted_mMtToPfMet,shifted_eMtToPfMet,abs(shifted_eDPhiToPfMet),abs(shifted_mDPhiToPfMet),abs(row.e_m_DR),float(row.ePt)/float(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)))])
    #                print "mPt  ",row.mPt
                        cut_flow_trk.Fill('jet0sel')
                        selections_passed.append(('0', 'selected'))
                    if jn == 1 :
                        if abs(shifted_eDPhiToPfMet) > 0.7 : break
#                        if deltaR(row.ePhi,row.mPhi,row.eEta,row.mEta)<1:break
                        if shifted_mMtToPfMet < 40 : break
                        if row.mPt < 25: break 
 #                       if shifted_eMtToPfMet > 65 : break
                        if row.ePt < 10 : break
                        selections_passed.extend([('1',i) for i in optimizer.compute_regions_1jet(row.mPt, row.ePt,deltaPhi(row.ePhi, row.mPhi),shifted_mMtToPfMet,shifted_eMtToPfMet,abs(shifted_eDPhiToPfMet),abs(shifted_mDPhiToPfMet),abs(row.e_m_DR),float(row.ePt)/float(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)))])
                    #                 print "mPt  ",row.mPt
                        cut_flow_trk.Fill('jet1sel')
                        selections_passed.append(('1', 'selected'))
                    if jn == 21 :
                        if row.mPt < 25: break 
                        if abs(shifted_eDPhiToPfMet) > 0.5 : break
  #                      if shifted_eMtToPfMet > 15 : break
#                        if shifted_mMtToPfMet < 15 : break
#                        if shifted_vbfDeta < 2.5 : break
                        if row.ePt < 10 : break # no cut as only electrons with pt>30 are in the ntuples
                        if shifted_vbfMass<100:break
                        selections_passed.extend([('21',i) for i in optimizer.compute_regions_2jetgg(row.mPt, row.ePt,shifted_mMtToPfMet,shifted_eMtToPfMet,abs(shifted_eDPhiToPfMet),abs(shifted_mDPhiToPfMet),shifted_vbfMass,shifted_vbfDeta,abs(row.e_m_DR),float(row.ePt)/float(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)))])
  #                  print "mPt  ",row.mPt
                        cut_flow_trk.Fill('jet2sel')
                        selections_passed.append(('21', 'selected'))
                    if jn == 22 :
                        if row.mPt < 25: break 
                        if abs(shifted_eDPhiToPfMet) > 0.3 : break
  #                      if shifted_eMtToPfMet > 15 : break
#                        if shifted_mMtToPfMet < 15 : break
#                        if shifted_vbfDeta < 2.5 : break
                        if row.ePt < 10 : break # no cut as only electrons with pt>30 are in the ntuples
                        selections_passed.extend([('22',i) for i in optimizer.compute_regions_2jetvbf(row.mPt, row.ePt,shifted_mMtToPfMet,shifted_eMtToPfMet,abs(shifted_eDPhiToPfMet),abs(shifted_mDPhiToPfMet),shifted_vbfMass,shifted_vbfDeta,abs(row.e_m_DR),float(row.ePt)/float(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)))])
  #                  print "mPt  ",row.mPt
                        cut_flow_trk.Fill('jet2sel')
                        selections_passed.append(('22', 'selected'))

#            print "continued on"
                for passed_selection in selections_passed:
                    folder = sign+'/'+processtype+'/'+str(passed_selection[0])+'/'+passed_selection[1]+'/'+sys
 #               print folder
                    self.fill_histos(row, folder, isSignalRegion,frw,btagweight,sys)
        cut_flow_trk.flush()        

   
             
            
    def finish(self):
        self.write_histos()

