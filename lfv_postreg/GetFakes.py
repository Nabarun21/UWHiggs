import math
import os
from sys import argv, stdout, stderr
import ROOT
import sys
import copy
ROOT.gROOT.SetStyle("Plain")
cat_now=['0','1','21','22']   #category names in analyzer                                                                                        
syst_names_now=['jetup','jetdown','tup','tdown','uup','udown']      #sysfolder names in analyzer                                             
vars = [
      ('h_collmass_pfmet', 'M_{coll}(e#mu) (GeV)', 1),
      ('mPt', 'p_{T}(mu) (GeV)', 4), 
      ('mEta', 'eta(mu)', 2),  
      ('mPhi', 'phi(mu)', 4), 
      ('ePt', 'p_{T}(e) (GeV)', 4), 
      ('eEta', 'eta(e)', 2),  
      ('ePhi', 'phi(e)', 4), 
      ('em_DeltaPhi', 'emu Deltaphi', 2), 
      ('em_DeltaR', 'emu Delta R', 2),
      ('h_vismass', 'M_{vis} (GeV)', 1),
      ('ePFMET_Mt', 'MT-e-MET (GeV)', 5),
      ('mPFMET_Mt', 'MT-mu-MET (GeV)', 5),
      ('ePFMET_DeltaPhi', 'Deltaphi-e-MET (GeV)', 2),
      ('mPFMET_DeltaPhi', 'Deltaphi-mu-MET (GeV)', 2),
      ('jetN_30', 'number of jets (p_{T} > 30 GeV)', 1),  
]

vars2 = [
      ('h_collmass_pfmet', 'M_{coll}(e#mu) (GeV)', 1),
]


histos={}
file=ROOT.TFile('LFVHEMuAnalyzerMVA_esffixed/data_obs.root')
for sign in ['os','ss']:
    for var in vars:
        for j in range(2):
            for i in range(4):
                if j==0:
                    hist_path="allfakes/"+sign+"/gg/"+cat_now[i]+"/"+var[0]
                else:
                    hist_path="allfakes/"+sign+"/gg/"+cat_now[i]+"/selected/nosys/"+var[0]
                histo=file.Get(hist_path)
                for bin in range(1,histo.GetNbinsX()+1):
                    binContent=histo.GetBinContent(bin)
                    binError=histo.GetBinError(bin)
                    histo.SetBinError(bin,math.sqrt(binError*binError+0.4*binContent*0.4*binContent))
                new_histo=copy.copy(histo)
                new_key=hist_path.split('/',1)[1]
                jojo= new_key.split('/')
                jojo1= '/'.join(new_key.split('/')[0:(len(jojo)-1)])
                new_key=new_key[:-17]
                new_key="subtracted/"+jojo1
                print new_key
                histos[(new_key,var[0])]=new_histo


for sign in ['os']:
    for var in vars2:
        for j in range(1):
            for i in range(4):
                for k in range(len(syst_names_now)):
                    hist_path="allfakes/os/gg/"+cat_now[i]+"/selected/"+syst_names_now[k]+"/h_collmass_pfmet"
                    histo_sys=file.Get(hist_path)
                    for bin in range(1,histo_sys.GetNbinsX()+1):
                        binContent=histo_sys.GetBinContent(bin)
                        binError=histo_sys.GetBinError(bin)
                        histo_sys.SetBinError(bin,math.sqrt(binError*binError+0.4*binContent*0.4*binContent))
                    new_histo_sys=copy.copy(histo_sys)
                    new_key=hist_path.split('/',1)[1]
                    jojo= new_key.split('/')
                    jojo1= '/'.join(new_key.split('/')[0:(len(jojo)-1)])
                    new_key=new_key[:-17]
                    new_key="subtracted/"+jojo1
                    print new_key
                    histos[(new_key,var[0])]=new_histo_sys


"""
fakeshapes
                if j==0:
                    hist_path="allfakesUp/"+sign+"/gg/"+cat_now[i]+"/h_collmass_pfmet"
                else:
                    hist_path="allfakesUp/"+sign+"/gg/"+cat_now[i]+"/selected/nosys/h_collmass_pfmet"
                histo=file.Get(hist_path)
                for bin in range(1,histo.GetNbinsX()+1):
                    binContent=histo.GetBinContent(bin)
                    binError=histo.GetBinError(bin)
                    histo.SetBinError(bin,math.sqrt(binError*binError+0.3*binContent*0.3*binContent))
                new_histo=copy.copy(histo)
                new_key=hist_path.split('/',1)[1]
                new_key=new_key[:-17]
                new_key="subtractedup/"+new_key
                histos[new_key]=new_histo



        for i in range(3):
            if j==0:
                hist_path="allfakesDown/"+sign+"/gg/"+cat_now[i]+"/h_collmass_pfmet"
            else:
                hist_path="allfakesDown/"+sign+"/gg/"+cat_now[i]+"/selected/nosys/h_collmass_pfmet"
            histo=file.Get(hist_path)
            new_histo=copy.copy(histo)
            new_key=hist_path.split('/',1)[1]
            new_key=new_key[:-17]
            new_key="subtracteddown/"+new_key
            histos[new_key]=new_histo
            


for i in range(3):
    for k in range(len(syst_names_now)):
        hist_path="allfakesUp/os/gg/"+cat_now[i]+"/selected/"+syst_names_now[k]+"/h_collmass_pfmet"
        histo_sys=file.Get(hist_path)
        new_histo_sys=copy.copy(histo_sys)
        new_key=hist_path.split('/',1)[1]
        new_key=new_key[:-17]
        new_key="subtractedup/"+new_key
        histos[new_key]=new_histo_sys
"""        

"""
        
for i in range(3):
    for k in range(len(syst_names_now)):
        hist_path="allfakesDown/os/gg/"+cat_now[i]+"/selected/"+syst_names_now[k]+"/h_collmass_pfmet"
        histo_sys=file.Get(hist_path)
        new_histo_sys=copy.copy(histo_sys)
        new_key=hist_path.split('/',1)[1]
        new_key=new_key[:-17]
        new_key="subtracteddown/"+new_key
        histos[new_key]=new_histo_sys
"""

outputfile=ROOT.TFile("FAKES.root","recreate")
outputfile.cd()
for key in histos.keys():
    print key
    dir0 = outputfile.mkdir(key[0]);
    dir0.Cd('FAKES.root:/'+key[0]);
    print dir0
    histos[key].Write()
outputfile.Close()
