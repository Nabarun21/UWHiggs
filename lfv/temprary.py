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

import optimizer as optimizer


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


class LFVHEMuAnalyzerMVA_optim(MegaBase):
    tree = 'em/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel='EMu'
        super(LFVHEMuAnalyzerMVA_optim, self).__init__(tree, outfile, **kwargs)
        target = os.path.basename(os.environ['megatarget'])
        self.is_WJet=('WJetsToLNu' in target or 'W1JetsToLNu' in target or 'W2JetsToLNu' in target or 'W3JetsToLNu' in target or 'W4JetsToLNu' in target)
        self.is_DYJet= ('DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target or  'DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target or 'DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target or 'DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target or 'DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target) 
        self.tree = EMTree(tree)
        self.out=outfile
        self.histograms = {}
        self.mym1='m'
        self.mye1='e'
        self.sysdir=['nosys','jetup','jetdown','tup','tdown','uup','udown']
#        self.sysdir=['nosys']
        if self.is_WJet:
            self.binned_weight=[0.003079413,0.00035568,0.000181734,0.000123275,0.000293479]
        elif self.is_DYJet:
            self.binned_weight=[0.000280512026,0.000027968676,0.000026799986,0.000033257606,0.000181650478]
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
        print "booking"

        cuts={}
        cuts[0]=optimizer.compute_regions_0jet(100000,100000,100000,1000,-1000000,-100000,100000,100000)+['selected']
        cuts[1]=optimizer.compute_regions_1jet(100000,100000,100000,1000,-1000000,-100000,100000,1000000)+['selected']
        cuts[2]=optimizer.compute_regions_2jet(100000,100000,100000,-1000,-1000000,109000,100000,1000000,100000)+['selected']

        cuts[3]=['selected']
        processtype=['gg']
        sign=[ 'ss','os']
        jetN = [0, 1, 2,3]
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

    
                    
    def fill_histos(self, row, f='os/gg/ept0/0',  isMuonTight=False, frw=[1.,1.,1.],sys=''):
        if self.is_WJet or self.is_DYJet:
            weight = [self.event_weight(row)*self.binned_weight[int(row.numGenJets)]]
        else:
            weight = [self.event_weight(row)]
        histos = self.histograms
        pudir =['']
        looseList = ['mLoose/']
        if not isMuonTight :
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
                if sys=='' or sys=='nosys':
                    histos[folder+'/h_collmass_pfmet'].Fill(row.e_m_collinearmass,weight[n])

                elif sys=='jetup':
                    histos[folder+'/h_collmass_pfmet'].Fill(row.e_m_collinearmass_JetEnUp,weight[n])
                
                elif sys=='jetdown':
                    histos[folder+'/h_collmass_pfmet'].Fill(row.e_m_collinearmass_JetEnDown,weight[n])

                elif sys=='tup':
                    histos[folder+'/h_collmass_pfmet'].Fill(row.e_m_collinearmass_TauEnUp,weight[n])


                elif sys=='tdown':
                    histos[folder+'/h_collmass_pfmet'].Fill(row.e_m_collinearmass_TauEnDown,weight[n])

                elif sys=='uup':
                    histos[folder+'/h_collmass_pfmet'].Fill(row.e_m_collinearmass_UnclusteredEnUp,weight[n])

                elif sys=='udown':
                    histos[folder+'/h_collmass_pfmet'].Fill(row.e_m_collinearmass_UnclusteredEnDown,weight[n])
        

    def process(self):
        
        cut_flow_histo = self.cut_flow_histo
        cut_flow_trk   = cut_flow_tracker(cut_flow_histo)
        myevent=()


        frw = []

        for row in self.tree:
            if row.jetVeto30!=0: continue
            sign = 'ss' if row.e_m_SS else 'os'

 #           ptthreshold = [30]
            processtype ='gg'##changed from 20

            cut_flow_trk.new_row(row.run,row.lumi,row.evt)
           
            cut_flow_trk.Fill('allEvents')

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

            #print "all else"
           
            #e Preselection
            if not selections.eSelection(row, 'e'): continue
            cut_flow_trk.Fill('esel')
            #print "ele sel--------------------------"
            
           
            if not selections.lepton_id_iso(row, 'e', 'eid15Loose_etauiso01'): continue
            cut_flow_trk.Fill('eiso')

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


                selections_passed=[]
                for dummy in [1]:
                    if jn == 0 :
                        if row.mPt < 30: break  
                        if row.ePt < 10 : break
                        if abs(shifted_eDPhiToPfMet) > 0.7 : break
                        if deltaPhi(row.ePhi, row.mPhi) < 2.5 : break
                        if shifted_mMtToPfMet < 60 : break
#                        if shifted_eMtToPfMet > 65 : break
                        selections_passed.extend([('0',i) for i in optimizer.compute_regions_0jet(row.mPt, row.ePt,deltaPhi(row.ePhi, row.mPhi),shifted_mMtToPfMet,shifted_eMtToPfMet,abs(shifted_eDPhiToPfMet),abs(shifted_mDPhiToPfMet),abs(row.e_m_DR))])
    #                print "mPt  ",row.mPt
                        cut_flow_trk.Fill('jet0sel')
                        selections_passed.append(('0', 'selected'))
                    if jn == 1 :
  #                      if abs(shifted_eDPhiToPfMet) > 0.7 : break
                        if deltaPhi(row.ePhi, row.mPhi) < 1.0 : break
                        if shifted_mMtToPfMet < 40 : break
                        if row.mPt < 25: break 
                        if shifted_eMtToPfMet > 65 : break
                        if row.ePt < 10 : break
                        selections_passed.extend([('1',i) for i in optimizer.compute_regions_1jet(row.mPt, row.ePt,deltaPhi(row.ePhi, row.mPhi),shifted_mMtToPfMet,shifted_eMtToPfMet,abs(shifted_eDPhiToPfMet),abs(shifted_mDPhiToPfMet),abs(row.e_m_DR))])
                    #                 print "mPt  ",row.mPt
#                        if deltaR(row.ePhi,row.mPhi,row.eEta,row.mEta)<1:break
                        cut_flow_trk.Fill('jet1sel')
                        selections_passed.append(('1', 'selected'))
                    if jn == 2 :
                        if row.mPt < 25: break 
                        if shifted_vbfMass < 200 : break
                        if abs(shifted_eDPhiToPfMet) > 0.3 : break
                        if shifted_eMtToPfMet > 15 : break
                        if shifted_mMtToPfMet < 15 : break
                       # if shifted_vbfDeta < 2.0 : break
                        if row.ePt < 10 : break # no cut as only electrons with pt>30 are in the ntuples
                        selections_passed.extend([('2',i) for i in optimizer.compute_regions_2jet(row.mPt, row.ePt,shifted_mMtToPfMet,shifted_eMtToPfMet,abs(shifted_eDPhiToPfMet),abs(shifted_mDPhiToPfMet),shifted_vbfMass,shifted_vbfDeta,abs(row.e_m_DR))])
  #                  print "mPt  ",row.mPt
                        cut_flow_trk.Fill('jet2sel')
                        selections_passed.append(('2', 'selected'))
#            print "continued on"
                for passed_selection in selections_passed:
                    folder = sign+'/'+processtype+'/'+str(passed_selection[0])+'/'+passed_selection[1]+'/'+sys
 #               print folder
                    self.fill_histos(row, folder, isMuonTight,frw,sys)
        cut_flow_trk.flush()        

   
             
            
    def finish(self):
        self.write_histos()


