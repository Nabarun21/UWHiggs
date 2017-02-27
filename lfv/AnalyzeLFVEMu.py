'''

Run LFV H->EMu analysis in the e+mu channel.

Authors: Maria Cepeda, Aaron Levine, Evan K. Friis, UW

'''

import EMuTree
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import glob
import os
#import FinalStateAnalysis.TagAndProbe.MuonPOGCorrections as MuonPOGCorrections
#import FinalStateAnalysis.TagAndProbe.H2TauCorrections as H2TauCorrections
import FinalStateAnalysis.TagAndProbe.PileupWeight as PileupWeight
import ROOT
import math

from math import sqrt, pi

def deltaPhi(phi1, phi2):
  PHI = abs(phi1-phi2)
  if PHI<=pi:
      return PHI
  else:
      return 2*pi-PHI

def fullMT(met,ept,mupt, metphi, ephi, muphi):
	ex=ept*math.cos(ephi)
	ey=ept*math.sin(ephi)
        metx=met*math.cos(metphi)
        mety=met*math.sin(metphi)
        mux=mupt*math.cos(muphi)
        muy=mupt*math.sin(muphi)
	full_et=met+ept+mupt # for e and mu I am approximating pt~et (M<<P)
	full_x=metx+ex+mux
        full_y=mety+ey+muy
	full_mt_2 = full_et*full_et-full_x*full_x-full_y*full_y
	full_mt=0
	if (full_mt_2>0):
		full_mt= math.sqrt(full_mt_2)
	return full_mt

def fullPT(met,ept,mupt, metphi, ephi, muphi):
        ex=ept*math.cos(ephi)
        ey=ept*math.sin(ephi)
        metx=met*math.cos(metphi)
        mety=met*math.sin(metphi)
        mux=mupt*math.cos(muphi)
        muy=mupt*math.sin(muphi)
        full_x=metx+ex+mux
        full_y=mety+ey+muy
        full_pt_2 = full_x*full_x+full_y*full_y
        full_pt=0
        if (full_pt_2>0):
                full_pt= math.sqrt(full_pt_2)
        return full_pt

def collMass_type1(row):
        mupx=row.mPt*math.cos(row.mPhi)
        mupy=row.mPt*math.sin(row.mPhi)
        metpx = row.type1_pfMetEt*math.cos(row.type1_pfMetPhi)
        metpy = row.type1_pfMetEt*math.sin(row.type1_pfMetPhi)
        met = sqrt(metpx*metpx+metpy*metpy)

        METproj= abs(metpx*mupx+metpy*mupy)/row.mPt

        xth=row.mPt/(row.mPt+METproj)
        den=math.sqrt(xth)

        mass=row.e_m_Mass/den

        #print '%4.2f, %4.2f, %4.2f, %4.2f, %4.2f' %(scaleMass(row), den, xth, METproj,mass)


        return mass


################################################################################
#### MC-DATA and PU corrections ################################################
################################################################################

class AnalyzeLFVEMu(MegaBase):
    tree = 'em/final/Ntuple'
    #tree = 'New_Tree'

    def __init__(self, tree, outfile, **kwargs):
        super(AnalyzeLFVEMu, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = EMuTree.EMuTree(tree)
        self.out = outfile
        self.histograms = {}

    def begin(self):

        #names=["oldisotight","oldisoloose","oldisogg","oldisoboost","oldisovbf", "oldisoloosegg","oldisolooseboost","oldisoloosevbf",
         #     "newisotight","newisoloose","newisogg","newisoboost","newisovbf",  "newisoloosegg","newisolooseboost","newisoloosevbf",
          #    "noisogg","noisoboost","noisovbf",
          #    "noTauID","noiso"]
        names=["preselection"]
        namesize = len(names)
	for x in range(0,namesize):

            self.book(names[x], "weight", "Event weight", 100, 0, 5)
            self.book(names[x], "GenWeight", "Gen level weight", 200000 ,-1000000, 1000000)
    
            self.book(names[x], "rho", "Fastjet #rho", 100, 0, 25)
            self.book(names[x], "nvtx", "Number of vertices", 100, -0.5, 100.5)
            self.book(names[x], "prescale", "HLT prescale", 21, -0.5, 20.5)
   
            self.book(names[x], "ePt", "Electron  Pt", 300,0,300)
            self.book(names[x], "eEta", "Electron  eta", 100, -2.5, 2.5)
            self.book(names[x], "eMtToPfMet_Ty1", "Electron MT (PF Ty1)", 200, 0, 200)
            self.book(names[x], "eCharge", "Electron Charge", 5, -2, 2)

            self.book(names[x], "mPt", "Mu  Pt", 300,0,300)
            self.book(names[x], "mEta", "Mu  eta", 100, -2.5, 2.5)
            self.book(names[x], "mMtToPfMet_Ty1", "Mu MT (PF Ty1)", 200, 0, 200)
            self.book(names[x], "mCharge", "Mu  Charge", 5, -2, 2)
	    self.book(names[x], "mJetPt", "Mu Jet Pt" , 500, 0 ,500)	    
            self.book(names[x], "mMass", "Mu  Mass", 1000, 0, 10)

		       
            self.book(names[x], 'eJetBtag', 'E 1 JetBtag', 100, -5.5, 9.5)
            self.book(names[x], 'mJetBtag', 'M 1 JetBtag', 100, -5.5, 9.5)
    	   
            self.book(names[x],"collMass_type1","collMass_type1",500,0,500);
            self.book(names[x],"fullMT_type1","fullMT_type1",500,0,500);
            self.book(names[x],"fullPT_type1","fullPT_type1",500,0,500);	    
    	    self.book(names[x], "LT", "ht", 400, 0, 400)
            self.book(names[x], "type1_pfMetEt", "Type1 MET", 200, 0, 200)
    
            self.book(names[x], "e_m_Mass", "Electron + Mu Mass", 200, 0, 200)
            self.book(names[x], "e_m_Pt", "Electron + Mu Pt", 200, 0, 200)
            self.book(names[x], "e_m_DR", "Electron + Mu DR", 100, 0, 10)
            self.book(names[x], "e_m_DPhi", "Electron + Mu DPhi", 100, 0, 4)
            self.book(names[x], "e_m_SS", "Electron + Mu SS", 5, -2, 2)
            self.book(names[x], "e_m_ToMETDPhi_Ty1", "Electron Mu DPhi to MET", 100, 0, 4)
    
            # Vetoes
            self.book(names[x], 'muVetoPt5IsoIdVtx', 'Number of extra muons', 5, -0.5, 4.5)
	    self.book(names[x], 'muVetoPt15IsoIdVtx', 'Number of extra muons', 5, -0.5, 4.5)
            self.book(names[x], 'tauVetoPt20Loose3HitsVtx', 'Number of extra taus', 5, -0.5, 4.5)
            self.book(names[x], 'eVetoMVAIso', 'Number of extra CiC tight electrons', 5, -0.5, 4.5)
   
            #self.book(names[x], 'jetVeto30PUCleanedTight', 'Number of extra jets', 5, -0.5, 4.5)
            #self.book(names[x], 'jetVeto30PUCleanedLoose', 'Number of extra jets', 5, -0.5, 4.5)
            self.book(names[x], 'jetVeto30', 'Number of extra jets', 5, -0.5, 4.5)	
	    #Isolation
	    self.book(names[x], 'eRelPFIsoDB' ,'Electron Isolation', 100, 0.0,1.0)
            self.book(names[x], 'mRelPFIsoDBDefault' ,'Muon Isolation', 100, 0.0,1.0)
   
 
            self.book(names[x], "ePhiMtPhi", "", 100, 0,4)
            self.book(names[x], "ePhiMETPhiType1", "", 100, 0,4)
            self.book(names[x], "mPhiMtPhi", "", 100, 0,4)
            self.book(names[x], "mPhiMETPhiType1", "", 100, 0,4)

### vbf ###
            self.book(names[x], "vbfJetVeto30", "central jet veto for vbf", 5, -0.5, 4.5)
	    self.book(names[x], "vbfJetVeto20", "", 5, -0.5, 4.5)
	    self.book(names[x], "vbfMVA", "", 100, 0,0.5)
	    self.book(names[x], "vbfMass", "", 500,0,5000.0)
	    self.book(names[x], "vbfDeta", "", 100, -0.5,10.0)
            self.book(names[x], "vbfj1eta","",100,-2.5,2.5)
	    self.book(names[x], "vbfj2eta","",100,-2.5,2.5)
	    self.book(names[x], "vbfVispt","",100,0,200)
	    self.book(names[x], "vbfHrap","",100,0,5.0)
	    self.book(names[x], "vbfDijetrap","",100,0,5.0)
	    self.book(names[x], "vbfDphihj","",100,0,4)
            self.book(names[x], "vbfDphihjnomet","",100,0,4)
            self.book(names[x], "vbfNJets", "g", 5, -0.5, 4.5)
            #self.book(names[x], "vbfNJetsPULoose", "g", 5, -0.5, 4.5)
            #self.book(names[x], "vbfNJetsPUTight", "g", 5, -0.5, 4.5)

	     

    def fill_histos(self, row,name='gg'):
        histos = self.histograms
        weight=1
        if (row.GenWeight>=0):
          weight=1
        else:
          weight=-1
        histos[name+'/weight'].Fill(weight)
        histos[name+'/GenWeight'].Fill(row.GenWeight)
        histos[name+'/rho'].Fill(row.rho, weight)
        histos[name+'/nvtx'].Fill(row.nvtx, weight)
        histos[name+'/prescale'].Fill(row.doubleEPrescale, weight)

        histos[name+'/ePt'].Fill(row.ePt, weight)
        histos[name+'/eEta'].Fill(row.eEta, weight)
        histos[name+'/eMtToPfMet_Ty1'].Fill(row.eMtToPfMet_Ty1,weight)
        histos[name+'/eCharge'].Fill(row.eCharge, weight)
        histos[name+'/mPt'].Fill(row.mPt, weight)
        histos[name+'/mEta'].Fill(row.mEta, weight)
        histos[name+'/mMtToPfMet_Ty1'].Fill(row.mMtToPfMet_Ty1,weight)
        histos[name+'/mCharge'].Fill(row.mCharge, weight)
	histos[name+'/mJetPt'].Fill(row.mJetPt, weight)

        histos[name+'/mMass'].Fill(row.mMass,weight)

	histos[name+'/LT'].Fill(row.LT,weight)

        histos[name+'/collMass_type1'].Fill(collMass_type1(row),weight)
        histos[name+'/fullMT_type1'].Fill(fullMT(row.type1_pfMetEt,row.ePt,row.mPt,row.type1_pfMetPhi, row.ePhi, row.mPhi),weight)
        histos[name+'/fullPT_type1'].Fill(fullPT(row.type1_pfMetEt,row.ePt,row.mPt,row.type1_pfMetPhi, row.ePhi, row.mPhi),weight) 

	histos[name+'/type1_pfMetEt'].Fill(row.type1_pfMetEt,weight)

        histos[name+'/e_m_Mass'].Fill(row.e_m_Mass,weight)
        histos[name+'/e_m_Pt'].Fill(row.e_m_Pt,weight)
        histos[name+'/e_m_DR'].Fill(row.e_m_DR,weight)
        histos[name+'/e_m_DPhi'].Fill(row.e_m_DPhi,weight)
        histos[name+'/e_m_SS'].Fill(row.e_m_SS,weight)
	histos[name+'/e_m_ToMETDPhi_Ty1'].Fill(row.e_m_ToMETDPhi_Ty1,weight)

        histos[name+'/eJetBtag'].Fill(row.eJetBtag, weight)
        histos[name+'/mJetBtag'].Fill(row.eJetBtag, weight)

        histos[name+'/muVetoPt5IsoIdVtx'].Fill(row.muVetoPt5IsoIdVtx, weight)
        histos[name+'/muVetoPt15IsoIdVtx'].Fill(row.muVetoPt15IsoIdVtx, weight)
        histos[name+'/tauVetoPt20Loose3HitsVtx'].Fill(row.tauVetoPt20Loose3HitsVtx, weight)
        histos[name+'/eVetoMVAIso'].Fill(row.eVetoMVAIso, weight)
        histos[name+'/jetVeto30'].Fill(row.jetVeto30, weight)
        #histos[name+'/jetVeto30PUCleanedLoose'].Fill(row.jetVeto30PUCleanedLoose, weight)
        #histos[name+'/jetVeto30PUCleanedTight'].Fill(row.jetVeto30PUCleanedTight, weight)

	histos[name+'/eRelPFIsoDB'].Fill(row.eRelPFIsoDB, weight)
        histos[name+'/mRelPFIsoDBDefault'].Fill(row.mRelPFIsoDBDefault,weight)
        
	histos[name+'/ePhiMtPhi'].Fill(deltaPhi(row.ePhi,row.mPhi),weight)
        histos[name+'/ePhiMETPhiType1'].Fill(deltaPhi(row.ePhi,row.type1_pfMetPhi),weight)
        histos[name+'/mPhiMETPhiType1'].Fill(deltaPhi(row.mPhi,row.type1_pfMetPhi),weight)
	histos[name+'/vbfJetVeto30'].Fill(row.vbfJetVeto30, weight)
     	histos[name+'/vbfJetVeto20'].Fill(row.vbfJetVeto20, weight)
        histos[name+'/vbfMVA'].Fill(row.vbfMVA, weight)
        histos[name+'/vbfMass'].Fill(row.vbfMass, weight)
        histos[name+'/vbfDeta'].Fill(row.vbfDeta, weight)
        histos[name+'/vbfj1eta'].Fill(row.vbfj1eta, weight)
        histos[name+'/vbfj2eta'].Fill(row.vbfj2eta, weight)
        histos[name+'/vbfVispt'].Fill(row.vbfVispt, weight)
        histos[name+'/vbfHrap'].Fill(row.vbfHrap, weight)
        histos[name+'/vbfDijetrap'].Fill(row.vbfDijetrap, weight)
        histos[name+'/vbfDphihj'].Fill(row.vbfDphihj, weight)
        histos[name+'/vbfDphihjnomet'].Fill(row.vbfDphihjnomet, weight)
        histos[name+'/vbfNJets'].Fill(row.vbfNJets, weight)
        #histos[name+'/vbfNJetsPULoose'].Fill(row.vbfNJetsPULoose, weight)
        #histos[name+'/vbfNJetsPUTight'].Fill(row.vbfNJetsPUTight, weight)





    def presel(self, row):
        if not row.isoMu24eta2p1Pass:
            return False
        return True

    def kinematics(self, row):
        #if row.mPt < 30:
        #    return False
        if row.ePt < 10:
             return False
        if abs(row.eEta) >= 2.3:
            return False
        if row.mPt<25 :
            return False
        if abs(row.mEta)>=2.1 :
            return False
        return True


    def oppositesign(self,row):
	if row.eCharge*row.mCharge!=-1:
            return False
	return True

    #def obj1_id(self, row):
    #    return bool(row.mPFIDTight)  
 
    def obj1_id(self,row):
    	 return (row.eMVANonTrigWP90>0) and (row.ePVDXY < 0.02) and (row.ePVDZ < 0.5)

    def obj2_id(self, row):
        return row.mIsGlobal and row.mIsPFMuon and (row.mNormTrkChi2<10) and (row.mMuonHits > 0) and (row.mMatchedStations > 1) and (row.mPVDXY < 0.02) and (row.mPVDZ < 0.5) and (row.mPixHits > 0) and (row.mTkLayersWithMeasurement > 5)

    def vetos(self,row):
		return  (bool (row.muVetoPt5IsoIdVtx<1) and bool (row.eVetoMVAIso<1) and bool (row.tauVetoPt20Loose3HitsVtx<1) )

    #def obj1_iso(self, row):
    #    return bool(row.mRelPFIsoDBDefault <0.12)
   
    def obj1_iso(self,row):
         return bool(row.eRelPFIsoDB <0.1)

    def obj2_iso(self, row):
        return  bool(row.mRelPFIsoDBDefault < 0.1)

    def obj1_antiiso(self, row):
        return bool(row.eRelPFIsoDB >0.2) 

    def process(self):
        event =0
        sel=False
        for row in self.tree:
            if event!=row.evt:   # This is just to ensure we get the (E,Mu) with the highest Pt
                event=row.evt    # In principle the code saves all the E+Mu posibilities, if an event has several combinations
                sel = False      # it will save them all.
            if sel==True:
                continue

            if not self.oppositesign(row):
                continue   

            if not self.kinematics(row): 
                continue
 
            if not self.obj1_iso(row):
                continue

            if not self.obj1_id(row):
                continue

            if not self.vetos (row):
                continue


            #if not self.obj2_id (row):
            #    continue

            #self.fill_histos(row,'noiso')

            if self.obj2_iso(row):
              self.fill_histos(row,'preselection')
#              self.fill_histos(row,'oldisotight')
#
#              if self.gg(row):
#                  self.fill_histos(row,'oldisogg')
#
#              if self.boost(row):
#                  self.fill_histos(row,'oldisoboost')
#
#              if self.vbf(row):
#                  self.fill_histos(row,'oldisovbf')
# 
#            elif self.obj2_looseiso(row):
#
#              self.fill_histos(row,'oldisoloose')
#
#              if self.gg(row):
#                  self.fill_histos(row,'oldisoloosegg')
#
#              if self.boost(row):
#                  self.fill_histos(row,'oldisolooseboost')
#
#              if self.vbf(row):
#                  self.fill_histos(row,'oldisoloosevbf')
#
#
#            if self.obj2_newiso(row):
#
#              self.fill_histos(row,'newisotight')
#
#              if self.gg(row):
#                  self.fill_histos(row,'newisogg')
#
#              if self.boost(row):
#                  self.fill_histos(row,'newisoboost')
#
#              if self.vbf(row):
#                  self.fill_histos(row,'newisovbf')
#
#            elif self.obj2_newlooseiso(row):
#              
#              self.fill_histos(row,'newisoloose')
#
#              if self.gg(row):
#                  self.fill_histos(row,'newisoloosegg')
#
#              if self.boost(row):
#                  self.fill_histos(row,'newisolooseboost')
#
#              if self.vbf(row):
#                  self.fill_histos(row,'newisoloosevbf')
#
#
#            if self.gg(row):
#                self.fill_histos(row,'noisogg')
#
#            if self.boost(row):
#                self.fill_histos(row,'noisoboost')
#
#            if self.vbf(row):
#                self.fill_histos(row,'noisovbf')
#
            

            sel=True

    def finish(self):
        self.write_histos()
