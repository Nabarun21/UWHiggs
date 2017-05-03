##Correction Factor still to add
from EEMTree import EEMTree
import os
import ROOT
import math
import optimizer
import glob
import array
import mcCorrections
import baseSelections as selections
import FinalStateAnalysis.PlotTools.pytree as pytree
from FinalStateAnalysis.PlotTools.decorators import  memo_last
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
from math import sqrt, pi, sin, cos, acos, sinh
from cutflowtracker import cut_flow_tracker
#Makes the cut flow histogram
cut_flow_step = ['allEvents', 'doubleEpass', 'bjetveto', 'esel','e2sel', 'msel', 'ZMass',  'mTightIso' ]

from inspect import currentframe

def get_linenumber():
    cf = currentframe()
    return cf.f_back.f_lineno

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
    

class MuFakeRateAnalyzerMVA_fromeem(MegaBase):
    tree = 'eem/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel='EEM'
        super(MuFakeRateAnalyzerMVA_fromeem, self).__init__(tree, outfile, **kwargs)
        self.tree = EEMTree(tree)
        self.out=outfile
        self.histograms = {}
        #self.pucorrector = mcCorrections.make_puCorrector('singlee')
        self.mye1 = 'e1'
        self.mye2 = 'e2'
        self.mym = 'm'
        #optimizer_keys   = [ i for i in optimizer.grid_search.keys() if i.startswith(self.channel) ]
        self.grid_search = {}
        #if len(optimizer_keys) > 1:
        #    for key in optimizer_keys:
        #        self.grid_search[key] = optimizer.grid_search[key]
        #else:
        #    self.grid_search[''] = optimizer.grid_search[optimizer_keys[0]]


    def event_weight(self, row):
        if row.run > 2: #FIXME! add tight ID correction
            return 1.

 
        #if bool(row.e1MatchesEle27WP80) and  not bool(row.e2MatchesEle27WP80) : etrig = 'e1'
        #if not bool(row.e1MatchesEle27WP80) and  bool(row.e2MatchesEle27WP80) :  etrig = 'e2'
        ##return self.pucorrector(row.nTruePU) * \
        ##    mcCorrections.eid_correction( row, self.mye1, self.mye2, self.mye3) * \
        ##    mcCorrections.eiso_correction(row, self.mye1, self.mye2, self.mye3) * \
        ##    mcCorrections.trig_correction(row, self.mye3   )
        return 1.


    def ee3DR(self, row):
        mye1_mye3_dr = 100.
        mye2_mye3_dr = 100.
        try:        
            mye1_mye3_dr = getattr(row, self.mye1+'_'+self.mye3+'_DR')
        except AttributeError:
            mye1_mye3_dr =getattr(row, self.mye3+'_'+self.mye1+'_DR')
        try :
            mye2_mye3_dr = getattr(row, self.mye2+'_'+self.mye3+'_DR')
        except AttributeError:
            mye2_mye3_dr =getattr(row, self.mye3+'_'+self.mye2+'_DR')

        return mye1_mye3_dr  if mye1_mye3_dr  < mye2_mye3_dr else mye1_mye3_dr 

    def ee3DPhi(self, row):
        e1e3DPhi=deltaPhi(getattr(row, self.mye1+'Phi'), getattr(row, self.mye3+'Phi'))
        e2e3DPhi=deltaPhi(getattr(row, self.mye2+'Phi'), getattr(row, self.mye3+'Phi'))
        return e1e3DPhi if e1e3DPhi < e2e3DPhi else e2e3DPhi

    def Z(self, row):
        e1p=ROOT.TVector3(getattr(row, self.mye1+'Pt')*cos(getattr(row, self.mye1+'Phi')),getattr(row, self.mye1+'Pt')*sin(getattr(row, self.mye1+'Phi')),getattr(row, self.mye1+'Pt')*sinh(getattr(row, self.mye1+'Eta')))
        e2p=ROOT.TVector3(getattr(row, self.mye2+'Pt')*cos(getattr(row, self.mye2+'Phi')),getattr(row, self.mye2+'Pt')*sin(getattr(row, self.mye2+'Phi')),getattr(row, self.mye2+'Pt')*sinh(getattr(row, self.mye2+'Eta')))
        e1FourVector= ROOT.TLorentzVector(e1p, sqrt(e1p.Mag2()+pow(getattr(row, self.mye1+'Mass'),2)))
        e2FourVector= ROOT.TLorentzVector(e2p, sqrt(e2p.Mag2()+pow(getattr(row, self.mye2+'Mass'),2)))
        zFourVector = e1FourVector+e2FourVector
        return zFourVector




##add the trigger correction 

    def begin(self):
        
        
        miso = ['mLoose', 'mTight']
        folder = []
        sign = ['ss','os']
        for iso in miso:
            for s in sign:
                folder.append(s+'/'+iso)
                j=0
                while j < 4 :
                    folder.append(s+'/'+iso+'/'+str(j))
                    j+=1
                    
        for f in folder: 
  
            self.book(f,"mPt", "m p_{T}", 200, 0, 200)
            ##self.book(f,"m3Phi", "m3 phi",  100, -3.2, 3.2)
            self.book(f,"mEta", "m eta", 46, -2.3, 2.3)
            self.book(f,"mAbsEta", "m abs eta", 23, 0, 2.3)
            self.book(f,"mPt_vs_mAbsEta", "m pt vs m abs eta", 23, 0, 2.3,  20, 0, 200.,  type=ROOT.TH2F)

            self.book(f, "ZMass",  " Inv Z Mass",  32, 0, 320)
  

        for s in sign:
            self.book(s+'/tNoCuts', "CUT_FLOW", "Cut Flow", len(cut_flow_step), 0, len(cut_flow_step))
            
            xaxis = self.histograms[s+'/tNoCuts/CUT_FLOW'].GetXaxis()
            self.cut_flow_histo = self.histograms[s+'/tNoCuts/CUT_FLOW']
            self.cut_flow_map   = {}
            for i, name in enumerate(cut_flow_step):
                xaxis.SetBinLabel(i+1, name)
                self.cut_flow_map[name] = i+0.5
                    
    def fill_histos(self, row, folder='os/tSuperLoose', fakeRate = False):
        weight = self.event_weight(row)
        histos = self.histograms
 
        histos[folder+'/mPt'].Fill( row.mPt , weight)
        histos[folder+'/mEta'].Fill(row.mEta, weight)
        ##histos[folder+'/m3Phi'].Fill(getattr(row, self.mye3+'Phi'), weight)
        histos[folder+'/mAbsEta'].Fill(abs(row.mEta), weight)
        histos[folder+'/mPt_vs_mAbsEta'].Fill(abs(row.mEta), row.mPt, weight)

        histos[folder+'/ZMass'].Fill(row.e1_e2_Mass, weight)
            

    def process(self):
        
        cut_flow_histo = self.cut_flow_histo
        cut_flow_trk   = cut_flow_tracker(cut_flow_histo)
        myevent =()
        #print self.tree.inputfilename
        for row in self.tree:
            jn = row.jetVeto30
            if jn > 3 : jn = 3
            cut_flow_trk.new_row(row.run,row.lumi,row.evt)
            #if row.run > 2: 
            cut_flow_trk.Fill('allEvents')

            if not bool(row.doubleEPass) : continue

            cut_flow_trk.Fill('doubleEpass')

            if row.bjetCISVVeto30Medium!=0 : continue 
                       
            cut_flow_trk.Fill('bjetveto')

            if not selections.eSelection(row, 'e1'): continue
            if not selections.lepton_id_iso(row, 'e1', 'eid15Loose_idiso06'): continue
            cut_flow_trk.Fill('esel')            
           
            
            if not selections.eSelection(row, 'e2'): continue
            if not selections.lepton_id_iso(row, 'e2', 'eid15Loose_idiso06'): continue
            cut_flow_trk.Fill('e2sel')
                       
            if not selections.muSelection(row, 'm'): continue
            if not selections.lepton_id_iso(row, 'm', 'muID_idiso06'): continue 
            cut_flow_trk.Fill('msel')
         
 
            if abs(row.e1_e2_Mass-91.2) > 25 : continue
            cut_flow_trk.Fill('ZMass')


            cut_flow_trk.Fill('MtToMet')


            miso = 'mLoose'
            sign = 'ss' if row.e1_e2_SS else 'os'
            folder = sign+'/'+miso
          
            self.fill_histos(row, folder)
            folder=folder+'/'+str(int(jn))
            self.fill_histos(row, folder)
            
            if selections.muTSelection(row, 'm'):
                miso = 'mTight' 
                folder = sign+'/'+miso
                self.fill_histos(row,  folder)
                cut_flow_trk.Fill('mTightIso')
                folder=folder+'/'+str(int(jn))
                self.fill_histos(row, folder)
                
                    
 
             
        cut_flow_trk.flush()
                                
             
            
    def finish(self):
        self.write_histos()
