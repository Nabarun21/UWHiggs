import array
import math
import os
from sys import argv, stdout, stderr
import ROOT
import sys
import copy
ROOT.gROOT.SetStyle("Plain")
cat_now=['0','1','21','22']   #category names in analyzer                                                                                        
syst_names_now=[]      #sysfolder names in analyzer                                             
if sys.argv[5]=="True":                                                                                                                                              
	syst_names_now=['jetup','jetdown','uup','udown']#,'mesup','mesdown','eesup','eesdown']      #sysfolder names in analyzer                                             
cutbasedvars = [
      ('h_collmass_pfmet', 'M_{coll}(e#mu) (GeV)', 1),
      ('mPt', 'p_{T}(mu) (GeV)', 4),
      ('mEta', 'eta(mu)', 1),
      ('mPhi', 'phi(mu)', 2),
      ('ePt', 'p_{T}(e) (GeV)', 4),
      ('eEta', 'eta(e)', 1),
      ('ePhi', 'phi(e)', 2),
      ('em_DeltaPhi', 'emu Deltaphi', 1),
      ('em_DeltaR', 'emu Delta R', 1),
      ('h_vismass', 'M_{vis} (GeV)', 1),
      ('Met', 'MET (GeV)', 1),
      ('ePFMET_Mt', 'MT-e-MET (GeV)', 5),
      ('mPFMET_Mt', 'MT-mu-MET (GeV)', 5),
      ('ePFMET_DeltaPhi', 'Deltaphi-e-MET (GeV)', 1),
      ('mPFMET_DeltaPhi', 'Deltaphi-mu-MET (GeV)', 1),
      ('jetN_30', 'number of jets (p_{T} > 30 GeV)', 1),
]

cutbasedvars2 = [
      ('h_collmass_pfmet', 'M_{coll}(e#mu) (GeV)', 1),
]

BDTvars = [
      ('BDT_value', 'BDT_value', 1),
      ('h_collmass_pfmet', 'M_{coll}(e#mu) (GeV)', 1),
      ('mPt', 'p_{T}(mu) (GeV)', 4),
      ('mEta', 'eta(mu)', 1),
      ('mPhi', 'phi(mu)', 2),
      ('ePt', 'p_{T}(e) (GeV)', 4),
      ('eEta', 'eta(e)', 1),
      ('ePhi', 'phi(e)', 2),
      ('em_DeltaPhi', 'emu Deltaphi', 1),
      ('em_DeltaR', 'emu Delta R', 1),
      ('h_vismass', 'M_{vis} (GeV)', 1),
      ('Met', 'MET (GeV)', 1),
      ('ePFMET_Mt', 'MT-e-MET (GeV)', 5),
      ('mPFMET_Mt', 'MT-mu-MET (GeV)', 5),
      ('ePFMET_DeltaPhi', 'Deltaphi-e-MET (GeV)', 1),
      ('mPFMET_DeltaPhi', 'Deltaphi-mu-MET (GeV)', 1),
      ('jetN_30', 'number of jets (p_{T} > 30 GeV)', 1),

]

BDTvars2 = [
      ('h_collmass_pfmet', 'M_{coll}(e#mu) (GeV)', 1),
]


BDT2vars = [
      ('BDT_value', 'BDT_value', 1),
      ('h_collmass_pfmet', 'M_{coll}(e#mu) (GeV)', 1),
      ('mPt', 'p_{T}(mu) (GeV)', 4),
      ('mEta', 'eta(mu)', 1),
      ('mPhi', 'phi(mu)', 2),
      ('ePt', 'p_{T}(e) (GeV)', 4),
      ('eEta', 'eta(e)', 1),
      ('ePhi', 'phi(e)', 2),
      ('em_DeltaPhi', 'emu Deltaphi', 1),
      ('em_DeltaR', 'emu Delta R', 1),
      ('h_vismass', 'M_{vis} (GeV)', 1),
      ('Met', 'MET (GeV)', 1),
      ('ePFMET_Mt', 'MT-e-MET (GeV)', 5),
      ('mPFMET_Mt', 'MT-mu-MET (GeV)', 5),
      ('ePFMET_DeltaPhi', 'Deltaphi-e-MET (GeV)', 1),
      ('mPFMET_DeltaPhi', 'Deltaphi-mu-MET (GeV)', 1),
      ('jetN_30', 'number of jets (p_{T} > 30 GeV)', 1),

]

BDT2vars2=[
	('BDT_value', 'BDT_value', 1)
]

if sys.argv[4]=='cut':
	vars=cutbasedvars
	vars2=cutbasedvars2
if sys.argv[4]=='BDT':
	vars=BDTvars
	vars2=BDTvars2
if sys.argv[4]=='BDT2':
	vars=BDT2vars
	vars2=BDT2vars2

binning=array.array( 'd', [-0.6,-0.5,-0.42,-0.34,-0.28,-0.22,-0.18,-0.14,-0.10,-0.06,-0.02,0.02,0.06,0.1,0.14,0.18,0.24,0.30,0.40,0.50,0.60])

binning1jet=array.array( 'd', [-0.6,-0.5,-0.42,-0.34,-0.26,-0.20,-0.16,-0.12,-0.08,-0.04,0.0,0.04,0.08,0.12,0.16,0.22,0.30,0.40,0.50,0.60])


binning2jet=array.array( 'd', [-0.6,-0.5,-0.42,-0.34,-0.28,-0.24,-0.20,-0.16,-0.12,-0.08,-0.04,0.0,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.34,0.42,0.50,0.60])



Analyzer=sys.argv[1]
Lumi=sys.argv[2]
histos={}
file=ROOT.TFile('LFVHEMuAnalyzerMVA'+Analyzer+Lumi+'/data_obs.root')
for sign in ['os','ss']:
    for var in vars:
        for j in range(2):
            for i in range(4):
                if j==0:
                    hist_path="antiIsolatedweighted/"+sign+"/gg/"+cat_now[i]+"/"+var[0]
                else:
                    hist_path="antiIsolatedweighted/"+sign+"/gg/"+cat_now[i]+"/selected/nosys/"+var[0]
                if j!=0 and ('collmass' not in var[0] and 'BDT' not in var[0]):continue
		if j!=0 and (sys.argv[4]=='BDT') and 'collmass' not in var[0]: continue
		if j!=0 and sys.argv[4]=='BDT2' and 'collmass' in var[0]:
			continue
                histo=file.Get(hist_path)


		rebin=var[2]
		if 'collmass' in var[0] or "MtToPfMet" in var[0] or "vismass" in var[0]:
			if (i==0 ):
				rebin = 5
			if ( i==1):
				rebin=10
			if ( i==2):
				rebin=25
			if ( i==3):
				rebin=25
		elif "BDT" in var[0]:
			if (i==0):
				histo=histo.Rebin(len(binning)-1,"",binning)
			elif (i==1):
				histo=histo.Rebin(len(binning1jet)-1,"",binning1jet)
			else:
				histo=histo.Rebin(len(binning2jet)-1,"",binning2jet)
		else:
			if (i==2):
				rebin = rebin*2
			if ( i==3 or  i==1 ):
				rebin=rebin*2
           
		if "BDT" not in var[0]:
			histo.Rebin(rebin)		
		lowBound=1
		highBound=1

		for bin in range(1,histo.GetNbinsX()):
			if histo.GetBinContent(bin) != 0:
				lowBound = bin
				break
		for bin in range(histo.GetNbinsX(),0,-1):
			if histo.GetBinContent(bin) != 0:
				highBound = bin
				break
		for bin in range(lowBound, highBound+1):
			if histo.GetBinContent(bin)<=0:
				histo.SetBinContent(bin,0.001)
				histo.SetBinError(bin,1.8)
		for bin in range(1,histo.GetNbinsX()+1):
			binContent=histo.GetBinContent(bin)
			binError=histo.GetBinError(bin)
			histo.SetBinError(bin,math.sqrt(binError*binError+0.4*binContent*0.4*binContent))

                new_histo=copy.copy(histo)
                new_key=hist_path.split('/',1)[1]
                jojo= new_key.split('/')
                jojo1= '/'.join(new_key.split('/')[0:(len(jojo)-1)])
                new_key=new_key[:-17]
                new_key="fakeRateMethod/"+jojo1
                #print new_key
                histos[(new_key,var[0])]=new_histo


for sign in ['os']:
    for var in vars2:
        for j in range(1):
            for i in range(4):
                for k in range(len(syst_names_now)):
                    hist_path="antiIsolatedweighted/os/gg/"+cat_now[i]+"/selected/"+syst_names_now[k]+"/"+var[0]
                    histo_sys=file.Get(hist_path)
		    for bin in range(1,histo_sys.GetNbinsX()):
			    if histo_sys.GetBinContent(bin) != 0:
				    lowBound = bin
				    break
		    for bin in range(histo_sys.GetNbinsX(),0,-1):
			    if histo_sys.GetBinContent(bin) != 0:
				    highBound = bin
				    break
		    for bin in range(lowBound, highBound+1):
			    if histo_sys.GetBinContent(bin)<=0:
				    histo_sys.SetBinContent(bin,0.001)
				    histo_sys.SetBinError(bin,1.8)
		    for bin in range(1,histo_sys.GetNbinsX()+1):
			    binContent=histo_sys.GetBinContent(bin)
			    binError=histo_sys.GetBinError(bin)
			    histo_sys.SetBinError(bin,math.sqrt(binError*binError+0.4*binContent*0.4*binContent))
                    new_histo_sys=copy.copy(histo_sys)
                    new_key=hist_path.split('/',1)[1]
                    jojo= new_key.split('/')
                    jojo1= '/'.join(new_key.split('/')[0:(len(jojo)-1)])
                    new_key=new_key[:-17]
                    new_key="fakeRateMethod/"+jojo1
                    #print new_key
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

outputfile=ROOT.TFile("FAKES"+sys.argv[1]+".root","recreate")
outputfile.cd()
for key in histos.keys():
    print key
    dir0 = outputfile.mkdir(key[0]);
    dir0.Cd('FAKES'+sys.argv[1]+'.root:/'+key[0]);
    print dir0
    histos[key].Write()
outputfile.Close()
