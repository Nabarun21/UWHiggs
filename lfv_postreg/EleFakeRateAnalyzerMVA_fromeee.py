##Correction Factor still to add
from EEETree import EEETree
import os
import ROOT
import math
import optimizer
import glob
import array
#import mcCorrections
import baseSelections as selections
import FinalStateAnalysis.PlotTools.pytree as pytree
from FinalStateAnalysis.PlotTools.decorators import  memo_last
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
from math import sqrt, pi, sin, cos, acos, sinh
from cutflowtracker import cut_flow_tracker
#Makes the cut flow histogram
cut_flow_step = ['allEvents','triggerpass','bjetvetopass', 'e1sel', 'e1IDiso', 'e2sel', 'e2IDiso','e3sel','e3IDiso', 'ZMass', 'tauveto', 'eveto', 'eTightIso' ] 

import FinalStateAnalysis.TagAndProbe.PileupWeight as PileupWeight
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
    


pu_distributions = glob.glob(os.path.join('inputs', os.environ['jobid'], 'data_SingleE*pu.root'))
pu_corrector = PileupWeight.PileupWeight('25ns_matchData', *pu_distributions)

class EleFakeRateAnalyzerMVA_fromeee(MegaBase):
    tree = 'eee/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel='EEE'
        super(EleFakeRateAnalyzerMVA_fromeee, self).__init__(tree, outfile, **kwargs)
        self.tree = EEETree(tree)
        self.out=outfile
        self.histograms = {}
       # self.pucorrector = mcCorrections.make_puCorrector('singlemu')
        self.mye1 = 'e1'
        self.mye2 = 'e2'
        self.mye3 = 'e3'
        #optimizer_keys   = [ i for i in optimizer.grid_search.keys() if i.startswith(self.channel) ]
        self.grid_search = {}
        #if len(optimizer_keys) > 1:
        #    for key in optimizer_keys:
        #        self.grid_search[key] = optimizer.grid_search[key]
        #else:
        #    self.grid_search[''] = optimizer.grid_search[optimizer_keys[0]]

    def mc_corrector_2015(self, row):
        
        pu = pu_corrector(row.nTruePU)
#        muidcorr1 = id_corrector(getattr(row, self.mym1+'Pt'), abs(getattr(row, self.mym1+'Eta')))
 #       muisocorr1 = iso_corrector('Tight', getattr(row, self.mym1+'Pt'), abs(getattr(row, self.mym1+'Eta')))
  #      muidcorr2 = id_corrector(getattr(row, self.mym2+'Pt'), abs(getattr(row, self.mym2+'Eta')))
   #     muisocorr2 = iso_corrector('Tight',getattr(row, self.mym2+'Pt'), abs(getattr(row, self.mym2+'Eta')))
    #    mutrcorr = tr_corrector(getattr(row, self.mym1+'Pt'), abs(getattr(row, self.mym1+'Eta'))) if getattr(row, self.mym1+'Pt')>getattr(row, self.mym2+'Pt') else  tr_corrector(getattr(row, self.mym2+'Pt'), abs(getattr(row, self.mym2+'Eta'))) #match the electron instead
     #   if pu*muidcorr1*muisocorr1*muidcorr2*muisocorr2*mutrcorr==0: print pu, muidcorr1, muisocorr1, muidcorr2, muisocorr2, mutrcorr
        return pu

    

    def correction(self,row):
	return self.mc_corrector_2015(row)

    def event_weight(self, row):
 
        if row.run > 2: #FIXME! add tight ID correction
            return 1.
#        if row.GenWeight*self.correction(row) == 0 : print 'weight==0', row.GenWeight*self.correction(row), row.GenWeight, self.correction(row), row.m1Pt, row.m2Pt, row.m1Eta, row.m2Eta
        return row.GenWeight*self.correction(row) 


    def event_weight(self, row):
        if row.run > 2: #FIXME! add tight ID correction
            return 1.

        return 1
        #if bool(row.e1MatchesEle27WP80) and  not bool(row.e2MatchesEle27WP80) : etrig = 'e1'
        """   #if not bool(row.e1MatchesEle27WP80) and  bool(row.e2MatchesEle27WP80) :  etrig = 'e2'
        return self.pucorrector(row.nTruePU) * \
            mcCorrections.eid_correction( row, self.mye1, self.mye2, self.mye3) * \
            mcCorrections.eiso_correction(row, self.mye1, self.mye2, self.mye3) * \
            mcCorrections.trig_correction(row, self.mye3   )
            """
    def me3DR(self, row):
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

    def me3DPhi(self, row):
        e1e3DPhi=deltaPhi(getattr(row, self.mye1+'Phi'), getattr(row, self.mye3+'Phi'))
        e2e3DPhi=deltaPhi(getattr(row, self.mye2+'Phi'), getattr(row, self.mye3+'Phi'))
        return e1mDPhi if e1e3DPhi < e2e3DPhi else e2e3DPhi

    def Zbos(self, row):
        e1p=ROOT.TVector3(getattr(row, self.mye1+'Pt')*cos(getattr(row, self.mye1+'Phi')),getattr(row, self.mye1+'Pt')*sin(getattr(row, self.mye1+'Phi')),getattr(row, self.mye1+'Pt')*sinh(getattr(row, self.mye1+'Eta')))
        e2p=ROOT.TVector3(getattr(row, self.mye2+'Pt')*cos(getattr(row, self.mye2+'Phi')),getattr(row, self.mye2+'Pt')*sin(getattr(row, self.mye2+'Phi')),getattr(row, self.mye2+'Pt')*sinh(getattr(row, self.mye2+'Eta')))
        e1FourVector= ROOT.TLorentzVector(e1p, sqrt(e1p.Mag2()+pow(getattr(row, self.mye1+'Mass'),2)))
        e2FourVector= ROOT.TLorentzVector(e2p, sqrt(e2p.Mag2()+pow(getattr(row, self.mye2+'Mass'),2)))
        zFourVector = e1FourVector+e2FourVector
        return zFourVector




##add the trigger correction 

    def begin(self):
        
        eiso = ['eLoose', 'eTight']
        folder = []
        sign = ['ss','os']
        for iso in eiso:
            for s in sign:
                folder.append(s+'/'+iso)
                j=0
                while j < 4 :
                    folder.append(s+'/'+iso+'/'+str(j))
                    j+=1
                    
        for f in folder: 
            
            self.book(f,"ePt", "e p_{T}", 200, 0, 200)
            self.book(f,"eEta", "e eta", 46, -2.3, 2.3)
            self.book(f,"eAbsEta", "e abs eta", 23, 0, 2.3)
            self.book(f,"ePt_vs_eAbsEta", "e pt vs e abs eta", 23, 0, 2.3,  20, 0, 200.,  type=ROOT.TH2F)
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
        print weight
        histos = self.histograms
        histos[folder+'/ePt'].Fill( getattr(row, self.mye3+'Pt' ), weight)
        histos[folder+'/eEta'].Fill(getattr(row, self.mye3+'Eta'), weight)
        histos[folder+'/eAbsEta'].Fill(abs(getattr(row, self.mye3+'Eta')), weight)
        histos[folder+'/ePt_vs_eAbsEta'].Fill(abs(getattr(row, self.mye3+'Eta')), getattr(row, self.mye3+'Pt'), weight)
        histos[folder+'/ZMass'].Fill(abs(getattr(row, self.mye1+'_'+self.mye2+'_'+'Mass')), weight)
            

    def process(self):
        
        cut_flow_histo = self.cut_flow_histo
        cut_flow_trk   = cut_flow_tracker(cut_flow_histo)
        myevent =()
        #print self.tree.inputfilename
        for row in self.tree:
            jn = row.jetVeto30
            if jn > 3 : jn = 3
#            print row.run,row.lumi,row.evt
            cut_flow_trk.new_row(row.run,row.lumi,row.evt)

            cut_flow_trk.Fill('allEvents')

            if not bool(row.singleE23WPLoosePass): continue
            
            cut_flow_trk.Fill('triggerpass')

            if row.bjetCISVVeto30Medium!=0 : continue 

            cut_flow_trk.Fill('bjetvetopass')

            if not selections.eSelection(row, 'e1'): continue
            cut_flow_trk.Fill('e1sel')
            if not selections.lepton_id_iso(row, 'e1', 'eid15Loose_idiso05'): continue
            if abs(row.e1Eta) > 1.4442 and abs(row.e1Eta < 1.566) : continue
            cut_flow_trk.Fill('e1IDiso')            
            

            if not selections.eSelection(row, 'e2'): continue
            cut_flow_trk.Fill('e2sel')
            if not selections.lepton_id_iso(row, 'e2', 'eid15Loose_idiso05'): continue
            if abs(row.e2Eta) > 1.4442 and abs(row.e2Eta) < 1.566 : continue
            cut_flow_trk.Fill('e2IDiso')

            if not selections.eSelection(row, 'e3'): continue
            if not selections.lepton_id_iso(row, 'e3', 'eid15Loose_idiso05'): continue #very loose loose eid13Tight_mvaLoose
            if abs(row.e3Eta) > 1.4442 and abs(row.e3Eta) < 1.566 : continue
            cut_flow_trk.Fill('e3IDiso')



            Zs= [(abs(row.e1_e2_Mass-91.2), ['e1', 'e2', 'e3']) , (abs(row.e2_e3_Mass-91.2), ['e2', 'e3', 'e1']), (abs(row.e1_e3_Mass-91.2), ['e1', 'e3', 'e2'])]
                
            for ele in range(0, 2) :
                
                if Zs[ele][0] == min(Zs[z][0] for z in range (0,2)) :
                    self.mye1 = Zs[ele][1][0]
                    self.mye2 = Zs[ele][1][1]
                    self.mye3 = Zs[ele][1][2]
                    

#            triggerpass=bool(bool(getattr(row, firstmuonfromZ)) or bool(getattr(row, secondmuonfromZ)))
           
 #           if not triggerpass: continue

            myZ=self.Zbos(row)
            Z_mass=myZ.M()
            if abs(Z_mass-91.2)>25: continue
            
            cut_flow_trk.Fill('ZMass')


            if row.tauVetoPt20Loose3HitsNewDMVtx : continue 
            cut_flow_trk.Fill('tauveto')

            if row.eVetoMVAIso: continue
            cut_flow_trk.Fill('eveto')
            
            
            #if (row.run, row.lumi, row.evt, row.e1Pt, row.e2Pt)==myevent: continue
            #myevent=(row.run, row.lumi, row.evt, row.e1Pt, row.e2Pt)

            eleiso = 'eLoose'
            sign = 'ss' if row.e1_e2_SS else 'os'
            folder = sign+'/'+eleiso
          
            self.fill_histos(row, folder)
            folder=folder+'/'+str(int(jn))
            self.fill_histos(row, folder)
            print "passed"
            
            if selections.lepton_id_iso(row, self.mye3, 'eid15Loose_etauiso01'):
                print "passed tight cut"
                eleiso = 'eTight' 
                folder = sign+'/'+eleiso
                self.fill_histos(row,  folder)
                cut_flow_trk.Fill('eTightIso')
                folder=folder+'/'+str(int(jn))
                self.fill_histos(row, folder)
                

 
             
        cut_flow_trk.flush()
                                
             
            
    def finish(self):
        self.write_histos()
