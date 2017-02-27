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




cut_flow_step=['allEvents','HLTIsoPasstrg','esel','eiso','musel','muiso','bjetveto','DR_e_mu','surplus_mu_veto','surplus_e_veto','jet0sel','jet1sel','jet2sel']

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
iso_corrector = MuonPOGCorrections.make_muon_pog_TightIso_2016BCD()
tr_corrector  = MuonPOGCorrections.make_muon_pog_IsoMu22oIsoTkMu22_2016BCD()
eId_corrector = EGammaPOGCorrections.make_egamma_pog_electronID_ICHEP2016( 'nontrigWP90')

class LFVHEMuAnalyzerMVA(MegaBase):
    tree = 'em/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):

        self.channel='EMu'
        super(LFVHEMuAnalyzerMVA, self).__init__(tree, outfile, **kwargs)
        target = os.path.basename(os.environ['megatarget'])
        self.target=target
        self.is_WJet=('WJetsToLNu' in target or 'W1JetsToLNu' in target or 'W2JetsToLNu' in target or 'W3JetsToLNu' in target or 'W4JetsToLNu' in target)
        self.is_DYJet= ('DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target or  'DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target or 'DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target or 'DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target or 'DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM' in target) 
        self.isData=('data' in target)
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
#        print "iso corr", muisocorr
        mutrcorr = tr_corrector(getattr(row, self.mym1+'Pt'), abs(getattr(row, self.mym1+'Eta'))) 
        eidcorr = eId_corrector(getattr(row,self.mye1+'Eta'),getattr(row, self.mye1+'Pt'))
 #       print "tr corr", mutrcorr
#       mutrcorr=1
     # if pu*muidcorr1*muisocorr1*muidcorr2*muisocorr2*mutrcorr==0: print pu, muidcorr1, muisocorr1, muidcorr2, muisocorr2, mutrcorr
  #      print "pileup--------   =",pu
   #     print pu*muidcorr*muisocorr*mutrcorr
#        return pu*muidcorr*muisocorr*mutrcorr*eidcorr
        return pu*muidcorr*mutrcorr*eidcorr
      

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

        processtype=['gg']
#        threshold=[]
        sign=[ 'ss','os']
        jetN = [0, 1, 2, 3]
        folder=[]
#        pudir = ['','p1s/', 'm1s/','trp1s/', 'trm1s/', 'eidp1s/','eidm1s/',  'eisop1s/','eisom1s/', 'mLoose/','mLooseUp/','mLooseDown/', ]
        pudir=['','mLoose/','mLooseUp/','mLooseDown/']
        for d  in pudir :
            for i in sign:
                for j in processtype:
 #                   for k in threshold:
                    for jn in jetN: 
                        folder.append(d+i+'/'+j+'/'+str(jn))
                        for s in self.sysdir:
                            folder.append(d+i+'/'+j+'/'+str(jn)+'/selected/'+s)
                        


        for f in folder:
            self.book(f,"mPt", "mu p_{T}", 200, 0, 200)
            self.book(f,"mPhi", "mu phi", 100, -3.2, 3.2)
            self.book(f,"mEta", "mu eta",  50, -2.5, 2.5)
            
            self.book(f,"ePt", "e p_{T}", 200, 0, 200)
            self.book(f,"ePhi", "e phi",  100, -3.2, 3.2)
            self.book(f,"eEta", "e eta", 50, 0, 5)
            
            self.book(f, "em_DeltaPhi", "e-mu DeltaPhi" , 50, 0, 3.2)
            self.book(f, "em_DeltaR", "e-mu DeltaR" , 100, -7, 7)
            
            self.book(f, "h_collmass_pfmet",  "h_collmass_pfmet",  32, 0, 320)
            self.book(f, "h_collmass_mvamet",  "h_collmass_mvamet",  32, 0, 320)

            self.book(f, "eMtToPfMet",  "eMtToPfMet",  32, 0, 320)
            self.book(f, "mMtToPfMet",  "mMtToPfMet",  32, 0, 320)
            self.book(f, "dPhiMetToE",  "dPhiMetToE",  50, 0, 3.2)

            self.book(f, "h_vismass",  "h_vismass",  32, 0, 320)
            self.book(f, "mPFMET_Mt", "mu-PFMET M_{T}" , 200, 0, 200)
            self.book(f, "ePFMET_Mt", "e-PFMET M_{T}" , 200, 0, 200)           
            self.book(f, "mPFMET_DeltaPhi", "mu-PFMET DeltaPhi" , 50, 0, 3.2)
            self.book(f, "ePFMET_DeltaPhi", "e-PFMET DeltaPhi" , 50, 0, 3.2)
            self.book(f, "mPFMETDeltaPhi_vs_ePFMETDeltaPhi",  "mPFMETDeltaPhi_vs_ePFMETDeltaPhi", 50, 0, 3.2, 50,0, 3.2, type=ROOT.TH2F)
            self.book(f, "vbfMass","vbf dijet mass",500,0,5000)
            self.book(f, "vbfDeta","vbf Delta Eta",50,0,5)

            self.book(f, "jetN_20", "Number of jets, p_{T}>20", 10, -0.5, 9.5) 
            self.book(f, "jetN_30", "Number of jets, p_{T}>30", 10, -0.5, 9.5) 
            

        for s in sign:
            self.book(s, "jetN_30", "Number of jets, p_{T}>30", 10, -0.5, 9.5) 
            self.book(s, "NUP", "Number of Partons", 12, -0.5, 11.5) 
            self.book(s, "numGenJets", "Number of Gen Level Jets", 12, -0.5, 11.5) 
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

    
    def fill_jet_histos(self,row,sign):

        if self.is_WJet or self.is_DYJet:
            weight = self.event_weight(row)*self.binned_weight[int(row.numGenJets)]
        else:
            weight = self.event_weight(row)

        self.histograms[sign+'/jetN_30'].Fill(row.jetVeto30,weight)
        self.histograms[sign+'/NUP'].Fill(row.NUP,weight)
        self.histograms[sign+'/numGenJets'].Fill(row.numGenJets,weight)

    def fill_histos(self, row, f='os/gg/ept0/0',  isMuonTight=False, frw=[1.,1.,1.],sys=''):
        if self.is_WJet or self.is_DYJet:
            weight = [self.event_weight(row)*self.binned_weight[int(row.numGenJets)]]
        else:
            weight = [self.event_weight(row)]
        histos = self.histograms
        pudir =['']
#        if row.run < 2: pudir.extend( ['p1s/', 'm1s/', 'trp1s/', 'trm1s/', 'eidp1s/','eidm1s/',  'eisop1s/','eisom1s/'])
        looseList = ['mLoose/', 'mLooseUp/', 'mLooseDown/']
#       looseList = ['mLoose/']
        #print isElecTight

        if not isMuonTight and (sys=='' or sys=='nosys'):
            frweight_bv = 1.
            ####frweight_bv = frw[0]/(1.-frw[0])
            err = frweight_bv*0.05
            ####err = abs(frw[0] - frw[1])* abs(frw[0]/((1-frw[0])*(1-frw[0])) + 1/(1-frw[0]))
            frweight_p1s = frweight_bv*(1+err)
            frweight_m1s = frweight_bv*(1-err)
            
            fr_weights = [frweight_bv, frweight_p1s, frweight_m1s]
            
            for n, l in enumerate(looseList) :
                frweight = weight[0]*fr_weights[n]
                folder = l+f
                histos[folder+'/mPt'].Fill(row.mPt, frweight)
                histos[folder+'/mEta'].Fill(row.mEta, frweight)
                histos[folder+'/mPhi'].Fill(row.mPhi, frweight) 
                histos[folder+'/ePt'].Fill(row.ePt, frweight)
                histos[folder+'/eEta'].Fill(row.eEta, frweight)
                histos[folder+'/ePhi'].Fill(row.ePhi, frweight)
                histos[folder+'/em_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.mPhi), frweight)
                histos[folder+'/em_DeltaR'].Fill(row.e_m_DR, frweight)
                histos[folder+'/h_vismass'].Fill(row.e_m_Mass, frweight)
                
                histos[folder+'/h_collmass_pfmet'].Fill(row.e_m_collinearmass,frweight)
                
                histos[folder+'/ePFMET_Mt'].Fill(row.eMtToPfMet_type1, frweight)
                histos[folder+'/mPFMET_Mt'].Fill(row.mMtToPfMet_type1, frweight)
                histos[folder+'/ePFMET_DeltaPhi'].Fill(abs(row.eDPhiToPfMet_type1), frweight)
                histos[folder+'/mPFMET_DeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1), frweight)
                histos[folder+'/mPFMETDeltaPhi_vs_ePFMETDeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1),abs(row.eDPhiToPfMet_type1) , frweight)
                
                histos[folder+'/vbfMass'].Fill(row.vbfMass, frweight)
                histos[folder+'/vbfDeta'].Fill(row.vbfDeta, frweight)
                histos[folder+'/jetN_30'].Fill(row.jetVeto30, frweight) 
                    

        else:
            for n,d  in enumerate(pudir) :
                folder = d+f                
                if sys=='' or sys=='nosys':
                    histos[folder+'/mPt'].Fill(row.mPt, weight[n])
                    histos[folder+'/mEta'].Fill(row.mEta, weight[n])
                    histos[folder+'/mPhi'].Fill(row.mPhi, weight[n]) 
                    histos[folder+'/ePt'].Fill(row.ePt, weight[n])
                    histos[folder+'/eEta'].Fill(row.eEta, weight[n])
                    histos[folder+'/ePhi'].Fill(row.ePhi, weight[n])
                    histos[folder+'/em_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.mPhi), weight[n])
                    histos[folder+'/em_DeltaR'].Fill(row.e_m_DR, weight[n])
                    histos[folder+'/h_vismass'].Fill(row.e_m_Mass, weight[n])
                    
                    histos[folder+'/h_collmass_pfmet'].Fill(row.e_m_collinearmass,weight[n])

                    histos[folder+'/ePFMET_Mt'].Fill(row.eMtToPfMet_type1, weight[n])
                    histos[folder+'/mPFMET_Mt'].Fill(row.mMtToPfMet_type1, weight[n])
                    histos[folder+'/ePFMET_DeltaPhi'].Fill(abs(row.eDPhiToPfMet_type1), weight[n])
                    histos[folder+'/mPFMET_DeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1), weight[n])
                    histos[folder+'/mPFMETDeltaPhi_vs_ePFMETDeltaPhi'].Fill(abs(row.mDPhiToPfMet_type1),abs(row.eDPhiToPfMet_type1) , weight[n])

                    histos[folder+'/vbfMass'].Fill(row.vbfMass, weight[n])
                    histos[folder+'/vbfDeta'].Fill(row.vbfDeta, weight[n])
                    histos[folder+'/jetN_30'].Fill(row.jetVeto30, weight[n]) 
                    

                elif sys=='jetup':
                    histos[folder+'/h_collmass_pfmet'].Fill(row.e_m_collinearmass_JetEnUp,weight[n])
                    histos[folder+'/ePFMET_Mt'].Fill(row.eMtToPfMet_JetEnUp, weight[n])
                    histos[folder+'/mPFMET_Mt'].Fill(row.mMtToPfMet_JetEnUp, weight[n])
                    histos[folder+'/ePFMET_DeltaPhi'].Fill(abs(row.eDPhiToPfMet_JetEnUp), weight[n])
                    histos[folder+'/mPFMET_DeltaPhi'].Fill(abs(row.mDPhiToPfMet_JetEnUp), weight[n])
                    histos[folder+'/jetN_30'].Fill(row.jetVeto30_JetEnUp, weight[n]) 
                    histos[folder+'/vbfMass'].Fill(row.vbfMass_JetEnUp, weight[n])
                    histos[folder+'/vbfDeta'].Fill(row.vbfDeta_JetEnUp, weight[n])


                
                elif sys=='jetdown':
                    histos[folder+'/h_collmass_pfmet'].Fill(row.e_m_collinearmass_JetEnDown,weight[n])
                    histos[folder+'/ePFMET_Mt'].Fill(row.eMtToPfMet_JetEnDown, weight[n])
                    histos[folder+'/mPFMET_Mt'].Fill(row.mMtToPfMet_JetEnDown, weight[n])
                    histos[folder+'/ePFMET_DeltaPhi'].Fill(abs(row.eDPhiToPfMet_JetEnDown), weight[n])
                    histos[folder+'/mPFMET_DeltaPhi'].Fill(abs(row.mDPhiToPfMet_JetEnDown), weight[n])
                    histos[folder+'/jetN_30'].Fill(row.jetVeto30_JetEnDown, weight[n]) 
                    histos[folder+'/vbfMass'].Fill(row.vbfMass_JetEnDown, weight[n])
                    histos[folder+'/vbfDeta'].Fill(row.vbfDeta_JetEnDown, weight[n])


                elif sys=='tup':
                    histos[folder+'/h_collmass_pfmet'].Fill(row.e_m_collinearmass_TauEnUp,weight[n])
                    histos[folder+'/ePFMET_Mt'].Fill(row.eMtToPfMet_TauEnUp, weight[n])
                    histos[folder+'/mPFMET_Mt'].Fill(row.mMtToPfMet_TauEnUp, weight[n])
                    histos[folder+'/ePFMET_DeltaPhi'].Fill(abs(row.eDPhiToPfMet_TauEnUp), weight[n])
                    histos[folder+'/mPFMET_DeltaPhi'].Fill(abs(row.mDPhiToPfMet_TauEnUp), weight[n])
                    histos[folder+'/jetN_30'].Fill(row.jetVeto30, weight[n]) 
                    histos[folder+'/vbfMass'].Fill(row.vbfMass, weight[n])
                    histos[folder+'/vbfDeta'].Fill(row.vbfDeta, weight[n])



                elif sys=='tdown':
                    histos[folder+'/h_collmass_pfmet'].Fill(row.e_m_collinearmass_TauEnDown,weight[n])
                    histos[folder+'/ePFMET_Mt'].Fill(row.eMtToPfMet_TauEnDown, weight[n])
                    histos[folder+'/mPFMET_Mt'].Fill(row.mMtToPfMet_TauEnDown, weight[n])
                    histos[folder+'/ePFMET_DeltaPhi'].Fill(abs(row.eDPhiToPfMet_TauEnDown), weight[n])
                    histos[folder+'/mPFMET_DeltaPhi'].Fill(abs(row.mDPhiToPfMet_TauEnDown), weight[n])
                    histos[folder+'/jetN_30'].Fill(row.jetVeto30, weight[n]) 
                    histos[folder+'/vbfMass'].Fill(row.vbfMass, weight[n])
                    histos[folder+'/vbfDeta'].Fill(row.vbfDeta, weight[n])


                elif sys=='uup':
                    histos[folder+'/h_collmass_pfmet'].Fill(row.e_m_collinearmass_UnclusteredEnUp,weight[n])
                    histos[folder+'/ePFMET_Mt'].Fill(row.eMtToPfMet_UnclusteredEnUp, weight[n])
                    histos[folder+'/mPFMET_Mt'].Fill(row.mMtToPfMet_UnclusteredEnUp, weight[n])
                    histos[folder+'/ePFMET_DeltaPhi'].Fill(abs(row.eDPhiToPfMet_UnclusteredEnUp), weight[n])
                    histos[folder+'/mPFMET_DeltaPhi'].Fill(abs(row.mDPhiToPfMet_UnclusteredEnUp), weight[n])
                    histos[folder+'/jetN_30'].Fill(row.jetVeto30, weight[n]) 
                    histos[folder+'/vbfMass'].Fill(row.vbfMass, weight[n])
                    histos[folder+'/vbfDeta'].Fill(row.vbfDeta, weight[n])


                elif sys=='udown':
                    folder = d+f
                    histos[folder+'/h_collmass_pfmet'].Fill(row.e_m_collinearmass_UnclusteredEnDown,weight[n])
                    histos[folder+'/ePFMET_Mt'].Fill(row.eMtToPfMet_UnclusteredEnDown, weight[n])
                    histos[folder+'/mPFMET_Mt'].Fill(row.mMtToPfMet_UnclusteredEnDown, weight[n])
                    histos[folder+'/ePFMET_DeltaPhi'].Fill(abs(row.eDPhiToPfMet_UnclusteredEnDown), weight[n])
                    histos[folder+'/mPFMET_DeltaPhi'].Fill(abs(row.mDPhiToPfMet_UnclusteredEnDown), weight[n])
                    histos[folder+'/jetN_30'].Fill(row.jetVeto30, weight[n]) 
                    histos[folder+'/vbfMass'].Fill(row.vbfMass, weight[n])
                    histos[folder+'/vbfDeta'].Fill(row.vbfDeta, weight[n])
        

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

            
            if not bool(row.singleIsoMu22Pass or row.singleIsoTkMu22Pass) and self.isData: 
                continue   #notrigger in new MC; add later

            cut_flow_trk.Fill('HLTIsoPasstrg')

            #E Preselection
            if not selections.eSelection(row, 'e'): continue
            cut_flow_trk.Fill('esel')
            #print "ele sel--------------------------"
            
           
            if not selections.lepton_id_iso(row, 'e', 'eid15Loose_etauiso1'): continue
            cut_flow_trk.Fill('eiso')


            #take care of ecal gap
            if row.eEta > 1.4442 and row.eEta < 1.566 : continue             


            if not selections.muSelection(row, 'm'): continue
            cut_flow_trk.Fill('musel')
            if not selections.lepton_id_iso(row, 'm', 'MuIDTight_idiso025'): continue
            cut_flow_trk.Fill('muiso')

            frw=1. ## add the correct fakerate weight once we have it


            isMuonTight=False

            if selections.lepton_id_iso(row, 'm', 'MuIDTight_mutauiso015'):isMuonTight=True

#            print isMuonTight
            #event cleaning
            isMuonTight=True


            if row.bjetCISVVeto30Loose : continue
            cut_flow_trk.Fill('bjetveto')

            if deltaR(row.ePhi,row.mPhi,row.eEta,row.mEta)<0.1:continue
            cut_flow_trk.Fill('DR_e_mu')
            
            if row.muVetoPt5IsoIdVtx :continue
            cut_flow_trk.Fill('surplus_mu_veto')

            if row.eVetoMVAIso :continue
            cut_flow_trk.Fill('surplus_e_veto')

            #print "all else"
           

            #mu Preselection


            #need to add Surpluse Veto
            #            for j in ptthreshold:


            if self.is_WJet or self.is_DYJet:
                self.fill_jet_histos(row,sign)
                
            jetN=row.jetVeto30
            if jetN>3:
                jetN=3

            folder = sign+'/'+processtype+'/'+str(int(jetN))
            self.fill_histos(row, folder,isMuonTight, frw)


            
            continue


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
                if jn == 0 :
                    if row.mPt < 30: continue 
                    if row.ePt < 22 : continue
                    if deltaPhi(row.ePhi, row.mPhi) < 2.5 : continue
                    if abs(shifted_eDPhiToPfMet) > 0.7 : continue
                    if shifted_mMtToPfMet < 60 : continue
 #                   if shifted_eMtToPfMet > 65 : continue
                    cut_flow_trk.Fill('jet0sel')
                
                if jn == 1 :
                    if row.mPt < 25: continue 
                    if row.ePt < 10 : continue
                    if abs(shifted_eDPhiToPfMet) > 0.7 : continue
                    if deltaR(row.ePhi,row.mPhi,row.eEta,row.mEta)<1:continue
                    if shifted_mMtToPfMet < 40 : continue
#                    if shifted_eMtToPfMet > 65 : continue
                    cut_flow_trk.Fill('jet1sel')
                    
                if jn == 2 :
                    if row.mPt < 25: continue 
                    if row.ePt < 10 : continue  #no cut as only electrons with pt>30 are in the ntuples
                    if abs(shifted_eDPhiToPfMet) > 0.3 : continue
                    if shifted_mMtToPfMet < 15 : continue
#                    if shifted_eMtToPfMet > 15 : continue
                    if shifted_vbfMass < 200 : continue
                    if shifted_vbfDeta < 2.5 : continue
                    cut_flow_trk.Fill('jet2sel')

                folder = sign+'/'+processtype+'/'+str(int(jn))+'/selected/'+sys
                self.fill_histos(row, folder, isMuonTight,frw,sys)
        cut_flow_trk.flush()        
            
    def finish(self):
        self.write_histos()


