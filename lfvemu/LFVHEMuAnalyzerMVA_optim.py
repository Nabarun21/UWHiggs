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
        cuts[0]=optimizer.compute_regions_0jet(100000,100000,100000,1000,-1000000,-100000)+['selected']
        cuts[1]=optimizer.compute_regions_1jet(100000,100000,100000,1000,-1000000,-100000)+['selected']
        cuts[2]=optimizer.compute_regions_2jet(100000,100000,100000,-1000,-1000000,100000,1000000)+['selected']
        cuts[3]=['selected']
        processtype=['gg']
        sign=[ 'ss','os']
        jetN = [0, 1, 2,3]
        folder=[]
        pudir=['','mLoose/']

        for d  in pudir :
            for i in sign:
                for j in processtype:
                    for jn in jetN: 
                        folder.append(d+i+'/'+j+'/'+str(jn))
                        for k in cuts[jn]:
                            folder.append(d+i+'/'+j+'/'+str(jn)+'/'+k)
                            
#        print folder
        for f in folder:
            self.book(f, "h_collmass_pfmet",  "h_collmass_pfmet",  32, 0, 320)
            self.book(f, "h_collmass_pfmet_jetdown",  "h_collmass_pfmet_jetdown",  32, 0, 320)
            self.book(f, "h_collmass_pfmet_jetup",  "h_collmass_pfmet_jetup",  32, 0, 320)
            self.book(f, "h_collmass_pfmet_tup",  "h_collmass_pfmet_tup",  32, 0, 320)
            self.book(f, "h_collmass_pfmet_uup",  "h_collmass_pfmet_uup",  32, 0, 320)
            self.book(f, "h_collmass_pfmet_tdown",  "h_collmass_pfmet_tdown",  32, 0, 320)
            self.book(f, "h_collmass_pfmet_udown",  "h_collmass_pfmet_udown",  32, 0, 320)

        for s in sign:
            self.book(s+'/tNoCuts', "CUT_FLOW", "Cut Flow", len(cut_flow_step), 0, len(cut_flow_step))
            xaxis = self.histograms[s+'/tNoCuts/CUT_FLOW'].GetXaxis()
            self.cut_flow_histo = self.histograms[s+'/tNoCuts/CUT_FLOW']
            self.cut_flow_map   = {}
            for i, name in enumerate(cut_flow_step):
                xaxis.SetBinLabel(i+1, name)
                self.cut_flow_map[name] = i+0.5

        print "booked"


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

    
                    
    def fill_histos(self, row, f='os/gg/ept0/0',  isMuonTight=False, frw=[1.,1.,1.]):
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
                histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi), weight[n])
                histos[folder+'/h_collmass_pfmet_jetdown'].Fill(row.e_m_collinearmass_JetEnDown,weight[n])
                histos[folder+'/h_collmass_pfmet_jetup'].Fill(row.e_m_collinearmass_JetEnUp,weight[n])
                histos[folder+'/h_collmass_pfmet_tdown'].Fill(row.e_m_collinearmass_TauEnDown,weight[n])
                histos[folder+'/h_collmass_pfmet_tup'].Fill(row.e_m_collinearmass_TauEnUp,weight[n])
                histos[folder+'/h_collmass_pfmet_uup'].Fill(row.e_m_collinearmass_UnclusteredEnDown,weight[n])
                histos[folder+'/h_collmass_pfmet_udown'].Fill(row.e_m_collinearmass_UnclusteredEnUp,weight[n])
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

#            if not bool(row.singleIsoMu20Pass) : continue   #notrigger in new MC; add later
            cut_flow_trk.Fill('HLTIsoPasstrg')


            jn = row.jetVeto30
            #print 'number of jets', jn
            if jn > 3 : jn = 3

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


            folder = sign+'/'+processtype+'/'+str(int(jn))
            self.fill_histos(row, folder,isMuonTight, frw)
            selections_passed=[]
            for i in [1]:
                if jn == 0 :
                    selections_passed.extend([('0',i) for i in optimizer.compute_regions_0jet(row.mPt, row.ePt,deltaPhi(row.ePhi, row.mPhi),row.mMtToPfMet_type1,row.eMtToPfMet_type1,deltaPhi(row.ePhi, row.raw_pfMetPhi))])
    #                print "mPt  ",row.mPt
                    if row.mPt < 50: break 
                    if row.ePt < 10 : break
                    if deltaPhi(row.ePhi, row.mPhi) < 2.7 : break
                    if deltaPhi(row.ePhi, row.raw_pfMetPhi) > 0.5 : break
                    if row.mMtToPfMet_type1 < 50 : break
                    if row.eMtToPfMet_type1 > 65 : break
                    cut_flow_trk.Fill('jet0sel')
                    selections_passed.append(('0', 'selected'))
                if jn == 1 :
                    selections_passed.extend([('1',i) for i in optimizer.compute_regions_1jet(row.mPt, row.ePt,deltaPhi(row.ePhi, row.mPhi),row.mMtToPfMet_type1,row.eMtToPfMet_type1,deltaPhi(row.ePhi, row.raw_pfMetPhi))])
   #                 print "mPt  ",row.mPt
                    if row.mPt < 45: break 
                    if row.ePt < 10 : break
                    if deltaPhi(row.ePhi, row.raw_pfMetPhi) > 0.5 : break
                    if deltaR(row.ePhi,row.mPhi,row.eEta,row.mEta)<1:break
                    if row.mMtToPfMet_type1 < 40 : break
                    if row.eMtToPfMet_type1 > 65 : break
                    cut_flow_trk.Fill('jet1sel')
                    selections_passed.append(('1', 'selected'))
                if jn == 2 :
                    selections_passed.extend([('2',i) for i in optimizer.compute_regions_2jet(row.mPt, row.ePt,row.mMtToPfMet_type1,row.eMtToPfMet_type1,deltaPhi(row.ePhi, row.raw_pfMetPhi),row.vbfMass,row.vbfDeta)])
  #                  print "mPt  ",row.mPt
                    if row.mPt < 25: break 
                    if row.ePt < 10 : break # no cut as only electrons with pt>30 are in the ntuples
                    if deltaPhi(row.ePhi, row.raw_pfMetPhi) > 0.3 : break
                    if row.mMtToPfMet_type1 < 15 : break
                    if row.eMtToPfMet_type1 > 25 : break
                    if row.vbfMass < 200 : break
                    if row.vbfDeta < 2.5 : break
                    cut_flow_trk.Fill('jet2sel')
                    selections_passed.append(('2', 'selected'))
#            print "continued on"
            for passed_selection in selections_passed:
                folder = sign+'/'+processtype+'/'+str(passed_selection[0])+'/'+passed_selection[1]
 #               print folder
                self.fill_histos(row, folder, isMuonTight,frw)
        cut_flow_trk.flush()        

   
             
            
    def finish(self):
        self.write_histos()


