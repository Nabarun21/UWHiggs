##Correction Factor still to add
from MMMTree import MMMTree
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
cut_flow_step = ['allEvents', 'doubleMuPass', 'bjetveto', 'm1sel', 'm2sel', 'm3sel', 'Zbosmass','tauveto','muveto','eveto', 'mTightIso' ]
import FinalStateAnalysis.TagAndProbe.PileupWeight as PileupWeight
import FinalStateAnalysis.TagAndProbe.MuonPOGCorrections as MuonPOGCorrections
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
    
pu_distributions = glob.glob(os.path.join('inputs', os.environ['jobid'], 'data_SingleMu*pu.root'))
pu_corrector = PileupWeight.PileupWeight('25ns_matchData', *pu_distributions)
id_corrector  = MuonPOGCorrections.make_muon_pog_PFTight_2015CD()
iso_corrector = MuonPOGCorrections.make_muon_pog_TightIso_2015CD()
tr_corrector  = MuonPOGCorrections.make_muon_pog_IsoMu20oIsoTkMu20_2015()
        
class MuFakeRateAnalyzerMVA(MegaBase):
    tree = 'mmm/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel='MMM'
        super(MuFakeRateAnalyzerMVA, self).__init__(tree, outfile, **kwargs)
        self.tree = MMMTree(tree)
        self.out=outfile
        self.histograms = {}
       # self.pucorrector = mcCorrections.make_puCorrector('singlemu')
        self.mym1 = 'm1'
        self.mym2 = 'm2'
        self.mym3 = 'm3'
        #optimizer_keys   = [ i for i in optimizer.grid_search.keys() if i.startswith(self.channel) ]
        self.grid_search = {}


    def mc_corrector_2015(self, row):
        
        pu = pu_corrector(row.nTruePU)
       # muidcorr1 = id_corrector(getattr(row, self.mym1+'Pt'), abs(getattr(row, self.mym1+'Eta')))
        muisocorr1 = iso_corrector('Tight', getattr(row, self.mym1+'Pt'), abs(getattr(row, self.mym1+'Eta')))
       # muidcorr2 = id_corrector(getattr(row, self.mym2+'Pt'), abs(getattr(row, self.mym2+'Eta')))
        muisocorr2 = iso_corrector('Tight',getattr(row, self.mym2+'Pt'), abs(getattr(row, self.mym2+'Eta')))
        mutrcorr = tr_corrector(getattr(row, self.mym1+'Pt'), abs(getattr(row, self.mym1+'Eta'))) if getattr(row, self.mym1+'Pt')>getattr(row, self.mym2+'Pt') else  tr_corrector(getattr(row, self.mym2+'Pt'), abs(getattr(row, self.mym2+'Eta'))) #match the electron instead
       # if pu*muidcorr1*muisocorr1*muidcorr2*muisocorr2*mutrcorr==0: print pu, muidcorr1, muisocorr1, muidcorr2, muisocorr2, mutrcorr
#        return pu*muidcorr1*muisocorr1*muidcorr2*muisocorr2*mutrcorr
        return pu*muisocorr1**muisocorr2*mutrcorr
    

    def correction(self,row):
	return self.mc_corrector_2015(row)
        
    def event_weight(self, row):
 
        if row.run > 2: #FIXME! add tight ID correction
            return 1.
       # if row.GenWeight*self.correction(row) == 0 : print 'weight==0', row.GenWeight*self.correction(row), row.GenWeight, self.correction(row), row.m1Pt, row.m2Pt, row.m1Eta, row.m2Eta
       # print row.GenWeight, "lkdfh"
        return row.GenWeight*self.correction(row) 


        #return 1
        #if bool(row.e1MatchesEle27WP80) and  not bool(row.e2MatchesEle27WP80) : etrig = 'e1'
        """   #if not bool(row.e1MatchesEle27WP80) and  bool(row.e2MatchesEle27WP80) :  etrig = 'e2'
        return self.pucorrector(row.nTruePU) * \
            mcCorrections.eid_correction( row, self.mym1, self.mym2, self.mym3) * \
            mcCorrections.eiso_correction(row, self.mym1, self.mym2, self.mym3) * \
            mcCorrections.trig_correction(row, self.mym3   )
            """
    def mm3DR(self, row):
        mym1_mym3_dr = 100.
        mym2_mym3_dr = 100.
        try:        
            mym1_mym3_dr = getattr(row, self.mym1+'_'+self.mym3+'_DR')
        except AttributeError:
            mym1_mym3_dr =getattr(row, self.mym3+'_'+self.mym1+'_DR')
        try :
            mym2_mym3_dr = getattr(row, self.mym2+'_'+self.mym3+'_DR')
        except AttributeError:
            mym2_mym3_dr =getattr(row, self.mym3+'_'+self.mym2+'_DR')

        return mym1_mym3_dr  if mym1_mym3_dr  < mym2_mym3_dr else mym1_mym3_dr 

    def mm3DPhi(self, row):
        m1m3DPhi=deltaPhi(getattr(row, self.mym1+'Phi'), getattr(row, self.mym3+'Phi'))
        m2m3DPhi=deltaPhi(getattr(row, self.mym2+'Phi'), getattr(row, self.mym3+'Phi'))
        return m1mDPhi if m1m3DPhi < m2m3DPhi else m2m3DPhi


    def Zbos(self, row):
        m1p=ROOT.TVector3(getattr(row, self.mym1+'Pt')*cos(getattr(row, self.mym1+'Phi')),getattr(row, self.mym1+'Pt')*sin(getattr(row, self.mym1+'Phi')),getattr(row, self.mym1+'Pt')*sinh(getattr(row, self.mym1+'Eta')))
        m2p=ROOT.TVector3(getattr(row, self.mym2+'Pt')*cos(getattr(row, self.mym2+'Phi')),getattr(row, self.mym2+'Pt')*sin(getattr(row, self.mym2+'Phi')),getattr(row, self.mym2+'Pt')*sinh(getattr(row, self.mym2+'Eta')))
        m1FourVector= ROOT.TLorentzVector(m1p, sqrt(m1p.Mag2()+pow(getattr(row, self.mym1+'Mass'),2)))
        m2FourVector= ROOT.TLorentzVector(m2p, sqrt(m2p.Mag2()+pow(getattr(row, self.mym2+'Mass'),2)))
        zFourVector = m1FourVector+m2FourVector
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
            

            self.book(f,"m3Pt", "m3 p_{T}", 200, 0, 200)

            self.book(f,"m3Eta", "m3 eta", 46, -2.3, 2.3)
            self.book(f,"m3AbsEta", "m3 abs eta", 23, 0, 2.3)
            self.book(f,"m3Pt_vs_m3AbsEta", "m3 pt vs m3 abs eta", 23, 0, 2.3,  20, 0, 200.,  type=ROOT.TH2F)
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
#        if weight==0 : 
        print "weight: ", weight

        histos[folder+'/m3Pt'].Fill( getattr(row, self.mym3+'Pt' ), weight)
        histos[folder+'/m3Eta'].Fill(getattr(row, self.mym3+'Eta'), weight)
        ##histos[folder+'/mPhi'].Fill(getattr(row, self.mym3+'Phi'), weight)
        histos[folder+'/m3AbsEta'].Fill(abs(getattr(row, self.mym3+'Eta')), weight)
        histos[folder+'/m3Pt_vs_m3AbsEta'].Fill(abs(getattr(row, self.mym3+'Eta')), getattr(row, self.mym3+'Pt'), weight)
        histos[folder+'/ZMass'].Fill(abs(getattr(row, self.mym1+'_'+self.mym2+'_'+'Mass')), weight)

    def process(self):
        
        cut_flow_histo = self.cut_flow_histo
        cut_flow_trk   = cut_flow_tracker(cut_flow_histo)
        myevent =()
        for row in self.tree:
            jn = row.jetVeto30
            if jn > 3 : jn = 3

            cut_flow_trk.Fill('allEvents')

            cut_flow_trk.new_row(row.run,row.lumi,row.evt)

            if not bool(row.singleIsoMu20Pass) and not bool(row.singleIsoTkMu20Pass) : continue
                
            cut_flow_trk.Fill('triggerselection')

            if row.bjetCISVVeto30Medium!=0 : continue 
          
            cut_flow_trk.Fill('bjetveto')
            #print row.m1JetPFCISVBtag, row.m1PixHits
            if not selections.muSelection(row, 'm1'): continue

            if not selections.lepton_id_iso(row, 'm1', 'muId_idiso025'): continue
            cut_flow_trk.Fill('m1sel')            
            
           
            if not selections.muSelection(row, 'm2'): continue
           
            if not selections.lepton_id_iso(row, 'm2', 'muId_idiso025'): continue
            cut_flow_trk.Fill('m2sel')
            
            #print row.m3JetPFCISVBtag, row.m3PixHits, row.m3Pt, row.m3AbsEta, row.m3PVDZ, row.m3PFIDTight, row.m3RelPFIsoDBDefault
            if not selections.muSelection(row, 'm3'): continue
            if not selections.lepton_id_iso(row, 'm3', 'muId_idiso025'): continue #very loose loose eid13Tight_mvaLoose
            cut_flow_trk.Fill('m3sel')
            

            Zs= [(abs(row.m1_m2_Mass-91.2), ['m1', 'm2', 'm3']) , (abs(row.m2_m3_Mass-91.2), ['m2', 'm3', 'm1']), (abs(row.m1_m3_Mass-91.2), ['m1', 'm3', 'm2'])]
                
            for ele in range(0, 2) :
                
                if Zs[ele][0] == min(Zs[z][0] for z in range (0,2)) :
                    self.mym1 = Zs[ele][1][0]
                    self.mym2 = Zs[ele][1][1]
                    self.mym3 = Zs[ele][1][2]
                    
  
            firstmuonfromZ=self.mym1+'MatchesSingleMu'        
            secondmuonfromZ=self.mym2+'MatchesSingleMu'        
            triggerpass=bool(bool(getattr(row, firstmuonfromZ)) or bool(getattr(row, secondmuonfromZ)))
           
#            if not triggerpass: continue

            myZ=self.Zbos(row)
            Z_mass=myZ.M()
            if abs(Z_mass-91.2)>25: continue
            
            cut_flow_trk.Fill('Zbosmass')

            if row.tauVetoPt20Loose3HitsNewDMVtx : continue 
            cut_flow_trk.Fill('tauveto')          
            if row.muVetoPt5IsoIdVtx : continue
            cut_flow_trk.Fill('muveto')          
            if row.eVetoMVAIso: continue
            cut_flow_trk.Fill('eveto')          
        



            miso = 'mLoose'
            sign = 'ss' if getattr(row, self.mym1+'_'+ self.mym2+'_SS') else 'os'
            folder = sign+'/'+miso
          
            self.fill_histos(row, folder)
            folder=folder+'/'+str(int(jn))
            self.fill_histos(row, folder)
            print "PASSED"
            if selections.muTSelection(row, self.mym3):
                miso = 'mTight' 
                folder = sign+'/'+miso
                self.fill_histos(row,  folder)
                cut_flow_trk.Fill('mTightIso')
                folder=folder+'/'+str(int(jn))
                self.fill_histos(row, folder)
                print "PASSED TIGHT"
             
        cut_flow_trk.flush()
                                
             
            
    def finish(self):
        self.write_histos()
