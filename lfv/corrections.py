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
import FinalStateAnalysis.TagAndProbe.HetauCorrection as HetauCorrection
import FinalStateAnalysis.TagAndProbe.FakeRate2D as FakeRate2D
import bTagSF as bTagSF
from inspect import currentframe

#pu_distributions = glob.glob(os.path.join('inputs', os.environ['jobid'], 'data_SingleMu*pu.root'))


#pu_corrector = PileupWeight.PileupWeight('MC_Spring16', *pu_distributions)
id_corrector  = MuonPOGCorrections.make_muon_pog_PFMedium_2016BCD()
iso_corrector = MuonPOGCorrections.make_muon_pog_TightIso_2016BCD()
trg_corrector  = MuonPOGCorrections.make_muon_pog_IsoMu22oIsoTkMu22_2016BCD()
mtrk_corrector = MuonPOGCorrections.mu_trackingEta_2016
#trk_corrector =  MuonPOGCorrections.make_muonptabove10_pog_tracking_corrections_2016()
#eId_corrector = EGammaPOGCorrections.make_egamma_pog_electronID_ICHEP2016( 'nontrigWP80')
etrk_corrector=EGammaPOGCorrections.make_egamma_pog_tracking_ICHEP2016()
eiso_corr0p10 =HetauCorrection.iso0p10_ele_2016
eiso_corr0p15 =HetauCorrection.iso0p15_ele_2016
checktrg=MuonPOGCorrections.make_muon_pog_IsoMu24oIsoTkMu24_2016ReReco()
checkid1=MuonPOGCorrections.make_muon_pog_PFTight_2016ReReco()
checkid2=MuonPOGCorrections.make_muon_pog_PFMedium_2016ReReco()
checkid3=MuonPOGCorrections.make_muon_pog_PFLoose_2016ReReco()

checkIso1=MuonPOGCorrections.make_muon_pog_TightIso_2016ReReco("Tight") #Muon Id as string input for iso corr
checkIso2=MuonPOGCorrections.make_muon_pog_TightIso_2016ReReco("Medium")
#checkIso3=MuonPOGCorrections.make_muon_pog_TightIso_2016ReReco("Loose")
checkIso4=MuonPOGCorrections.make_muon_pog_LooseIso_2016ReReco("Tight")
checkIso5=MuonPOGCorrections.make_muon_pog_LooseIso_2016ReReco("Medium")
checkIso6=MuonPOGCorrections.make_muon_pog_LooseIso_2016ReReco("Loose")


print checktrg(40,0.6)       #all calls take argument pt and abseta
print ""
print checkid1(121,0.4)
print ""
print checkid2(55,0.4)
print ""
print checkid3(55,0.4)


print checkIso1(55,0.4)
print ""
print checkIso2(55,0.4)
print ""
#print checkIso3(55,0.4)
print ""
print checkIso4(55,0.4)
print ""
print checkIso5(55,0.4)
print ""
print checkIso6(55,0.4)
