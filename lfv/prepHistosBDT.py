import array
import os 
from sys import argv, stdout, stderr
import ROOT
import sys
import copy
import getopt
import math
import array


#argv 1 =analyzer;argv2=lumi, argv3=jobid
ROOT.gROOT.SetStyle("Plain")
cat_now=['0','1','21','22']   #category names in analyzer                                                                                    
                                                                                                                                              
syst_names_now=['jetup','jetdown','tup','tdown','uup','udown']      #sysfolder names in analyzer                                             

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
      ('pZeta', 'pZeta', 1),
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
#
]

BDT2vars2=[
	('BDT_value', 'BDT_value', 1)]


if sys.argv[4]=='cut':
	vars=cutbasedvars
	vars2=cutbasedvars2
if sys.argv[4]=='BDT':
	vars=BDTvars
	vars2=BDTvars2
if sys.argv[4]=='BDT2':
	vars=BDT2vars
	vars2=BDT2vars2


commonvars=[('numVertices','number of vertices', 1),
	    ('numGenJets','Number of jets',1),
	    ('NUP','Number of Partons',1),
	    ('jetN_30','Number of jets with p_{T}>30',1),
	    ('h_collmass_pfmet', 'M_{coll}(e#mu) (GeV)', 1),
	    ('h_vismass', 'M_{vis}(e#mu) (GeV)', 1)
	    ]


binning=array.array( 'd', [-0.6,-0.5,-0.42,-0.34,-0.28,-0.22,-0.18,-0.14,-0.10,-0.06,-0.02,0.02,0.06,0.1,0.14,0.18,0.24,0.30])

binning1jet=array.array( 'd', [-0.6,-0.5,-0.42,-0.34,-0.26,-0.20,-0.16,-0.12,-0.08,-0.04,0.0,0.04,0.08,0.12,0.16,0.22,0.30])


binning2jet=array.array( 'd', [-0.6,-0.5,-0.42,-0.34,-0.28,-0.24,-0.20,-0.16,-0.12,-0.08,-0.04,0.0,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.30])


Analyzer='LFVHEMuAnalyzerMVA'+sys.argv[1]+sys.argv[2]+'plot'
JSONlumi=int(sys.argv[2])
def NoNegBins(histo):
        for i in range(1,histo.GetNbinsX()+1):
                if histo.GetBinContent(i) != 0:
                        lowBound = i
                        break
        for i in range(histo.GetNbinsX(),0,-1):
                if histo.GetBinContent(i) != 0:
                        highBound = i
                        break
        for i in range(lowBound, highBound+1):
                if histo.GetBinContent(i) <= 0:
                        histo.SetBinContent(i,0)

def set_poissonerrors(histo):
        histo.SetBinErrorOption(ROOT.TH1.kPoisson)

        for i in range(1,histo.GetNbinsX()+1):
                errorLow = histo.GetBinErrorLow(i)
                errorUp = histo.GetBinErrorUp(i)

def yieldHisto(histo,xmin,xmax):
        binmin = int(histo.FindBin(xmin))
        binwidth = histo.GetBinWidth(binmin)
        binmax = int(xmax/binwidth)
        signal = histo.Integral(binmin,binmax)
        return signal

def do_binbybin(histo,eff_lumi,lowBound,highBound,norm_uncert,rebin=None): #fill empty bins                                                  
	lowBound=histo.GetNbinsX()
	highBound=1
        for i in range(1,lowBound):
                if histo.GetBinContent(i) != 0:
                        lowBound = i
                        break
        for i in range(histo.GetNbinsX(),highBound,-1):
                if histo.GetBinContent(i) != 0:
                        highBound = i
                        break
	fillEmptyBins=True
	histo.Scale(JSONlumi)
        for i in range(lowBound, highBound+1):
                if fillEmptyBins: #fill empty bins                                                                                                                                 
                        if histo.GetBinContent(i) <=0:
                                histo.SetBinContent(i,0.001*eff_lumi*JSONlumi)
                                histo.SetBinError(i,1.8*eff_lumi*JSONlumi)
                else:
                        if histo.GetBinContent(i) <0:
                                histo.SetBinContent(i,0.001*eff_lumi*JSONlumi)
                                histo.SetBinError(i,1.8*eff_lumi*JSONlumi)
				
	for bin in range(1,histo.GetNbinsX()+1):
		binContent=histo.GetBinContent(bin)
		binError=histo.GetBinError(bin)
		histo.SetBinError(bin,math.sqrt(binError*binError+norm_uncert*binContent*norm_uncert*binContent))


regions= ['os','ss','antiIsolated/os','antiIsolated/ss']
#regions= ['os','ss','fakeRateMethod/os','fakeRateMethod/ss','antiIsolated/os','antiIsolated/ss','antiIsolatedweighted/os','antiIsolatedweighted/ss','antiIsolatedweightedmuon/ss','antiIsolatedweightedmuon/os','antiIsolatedweightedelectron/ss','antiIsolatedweightedelectron/os','antiIsolatedweightedmuonelectron/os','antiIsolatedweightedmuonelectron/ss']:
for sign in regions:
    for var in vars:
        for j in range(2):
            for i in range(4):
		     if j==0:
			     hist_path=sign+"/gg/"+cat_now[i]+"/"+var[0]
		     else:
			     hist_path=sign+"/gg/"+cat_now[i]+"/selected/nosys/"+var[0]
		     if j!=0 and 'collmass' not in var[0] and 'BDT' not in var[0]:
			     continue
		     if j!=0 and sys.argv[4]=='BDT' and 'collmass' not in var[0]:
			     continue
		     if j!=0 and sys.argv[4]=='BDT2' and 'collmass' in var[0]:
                             continue
		     jojo= hist_path.split('/')
		     folder= '/'.join(hist_path.split('/')[0:(len(jojo)-1)])
		     lowDataBin = 1
#		     highDataBin = data_histo.GetNbinsX()
		     highDataBin=1
#	     for i in range(1,data_histo.GetNbinsX()+1):
#		     if (data_histo.GetBinContent(i) > 0):
#			     lowDataBin = i
#			     break
#
#	     for i in range(data_histo.GetNbinsX(),0,-1):
#		     if (data_histo.GetBinContent(i) > 0):
#			     highDataBin = i
#			     break
#
#
#		     print hist_path
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
				     rebin=binning
			     elif (i==1):
				     rebin=binning1jet
			     else:
				     rebin=binning2jet
		     else:
			     if (i==2):
				     rebin = rebin*2
			     if ( i==3 or  i==1 ):
				     rebin=rebin*2

		     
		     
		     datafile=ROOT.TFile(Analyzer+"/data_obs.root","UPDATE")
		     data_histo=datafile.Get(hist_path)
		     
		     if isinstance(rebin,type(binning)):
			     data_histo=data_histo.Rebin(len(rebin)-1,"",rebin)
		     else:     
			     data_histo.Rebin(rebin)
		     datafile.cd(folder)
                     data_histo.Write()

                     datafile.Close()


		     GGfile=ROOT.TFile(Analyzer+"/LFVGG125.root","UPDATE")
		     GG_histo=GGfile.Get(hist_path)

		
		     if isinstance(rebin,type(binning)):
			     GG_histo=GG_histo.Rebin(len(rebin)-1,"",rebin)
		     else:     
			     GG_histo.Rebin(rebin)
		  		     
		     GGfile.cd(folder)
                     GG_histo.Write()

                     GGfile.Close()


		     VBFfile=ROOT.TFile(Analyzer+"/LFVVBF125.root","UPDATE")
		     VBF_histo=VBFfile.Get(hist_path)

		     if isinstance(rebin,type(binning)):
			     VBF_histo=VBF_histo.Rebin(len(rebin)-1,"",rebin)
		     else:     
			     VBF_histo.Rebin(rebin)

		     VBFfile.cd(folder)
                     VBF_histo.Write()

                     VBFfile.Close()


		     WGfile=ROOT.TFile(Analyzer+"/WG.root","UPDATE")


		     WGhisto=WGfile.Get(hist_path)
		     


		     if isinstance(rebin,type(binning)):

			     WGhisto=WGhisto.Rebin(len(rebin)-1,"",rebin)
		     else:     
			     WGhisto.Rebin(rebin)

		     do_binbybin(WGhisto,1.25890428e-06,lowDataBin,highDataBin,0.1,rebin)
		     
		     WGfile.cd(folder)
#		     print WGhisto.GetNbinsX()
		     WGhisto.Write()

		     WGfile.Close()

		     STfile=ROOT.TFile(Analyzer+"/T.root","UPDATE")


		     SThisto=STfile.Get(hist_path)

		     if isinstance(rebin,type(binning)):
			     SThisto=SThisto.Rebin(len(rebin)-1,"",rebin)
		     else:     
			     SThisto.Rebin(rebin)

		     do_binbybin(SThisto,5.23465826064e-06,lowDataBin,highDataBin,0.05,rebin)

		     STfile.cd(folder)
		     SThisto.Write()

		     STfile.Close()


		     
		     TTfile=ROOT.TFile(Analyzer+"/TT.root","UPDATE")


		     TThisto=TTfile.Get(hist_path)
		     if isinstance(rebin,type(binning)):
			     TThisto=TThisto.Rebin(len(rebin)-1,"",rebin)
		     else:     
			     TThisto.Rebin(rebin)

		     do_binbybin(TThisto,1.08709111195e-05,lowDataBin,highDataBin,0.06,rebin)

		     TTfile.cd(folder)
		     TThisto.Write()

		     TTfile.Close()

		     Dibosonsfile=ROOT.TFile(Analyzer+"/Dibosons.root","UPDATE")


		     Dibosonshisto=Dibosonsfile.Get(hist_path)
		     if isinstance(rebin,type(binning)):
			     Dibosonshisto=Dibosonshisto.Rebin(len(rebin)-1,"",rebin)
		     else:     
			     Dibosonshisto.Rebin(rebin)

		     do_binbybin(Dibosonshisto,1.17948019785e-05,lowDataBin,highDataBin,0.05,rebin)

		     Dibosonsfile.cd(folder)
		     Dibosonshisto.Write()

		     Dibosonsfile.Close()


		     DYfile=ROOT.TFile(Analyzer+"/DY.root","UPDATE")


		     DYhisto=DYfile.Get(hist_path)
		     if isinstance(rebin,type(binning)):
			     DYhisto=DYhisto.Rebin(len(rebin)-1,"",rebin)
		     else:     
			     DYhisto.Rebin(rebin)
 
		     do_binbybin(DYhisto,1.1e-05,lowDataBin,highDataBin,0.1,rebin)

		     DYfile.cd(folder)
		     DYhisto.Write()

		     DYfile.Close()

		     ZTauTaufile=ROOT.TFile(Analyzer+"/ZTauTau.root","UPDATE")

		     ZTauTauhisto=ZTauTaufile.Get(hist_path)
		     if isinstance(rebin,type(binning)):
			     ZTauTauhisto=ZTauTauhisto.Rebin(len(rebin)-1,"",rebin)
		     else:     
			     ZTauTauhisto.Rebin(rebin)


		     do_binbybin(ZTauTauhisto,1.1e-05,lowDataBin,highDataBin,0.1,rebin)

		     ZTauTaufile.cd(folder)
		     ZTauTauhisto.Write()

		     ZTauTaufile.Close()

		     WJETSMCfile=ROOT.TFile(Analyzer+"/WJETSMC.root","UPDATE")


		     WJETSMChisto=WJETSMCfile.Get(hist_path)

		     if isinstance(rebin,type(binning)):
			     WJETSMChisto=WJETSMChisto.Rebin(len(rebin)-1,"",rebin)
		     else:     
			     WJETSMChisto.Rebin(rebin)


		     do_binbybin(WJETSMChisto,1.5e-04,lowDataBin,highDataBin,0.1,rebin)

		     WJETSMCfile.cd(folder)
		     WJETSMChisto.Write()

		     WJETSMCfile.Close()

		     vbfHTauTaufile=ROOT.TFile(Analyzer+"/vbfHTauTau.root","UPDATE")


		     vbfHTauTauhisto=vbfHTauTaufile.Get(hist_path)

		     if isinstance(rebin,type(binning)):
			     vbfHTauTauhisto=vbfHTauTauhisto.Rebin(len(rebin)-1,"",rebin)
		     else:     
			     vbfHTauTauhisto.Rebin(rebin)


		     do_binbybin(vbfHTauTauhisto,4.2e-08,lowDataBin,highDataBin,0.1,rebin)

		     vbfHTauTaufile.cd(folder)
		     vbfHTauTauhisto.Write()

		     vbfHTauTaufile.Close()


		     ggHTauTaufile=ROOT.TFile(Analyzer+"/ggHTauTau.root","UPDATE")


		     ggHTauTauhisto=ggHTauTaufile.Get(hist_path)
		     if isinstance(rebin,type(binning)):
			     ggHTauTauhisto=ggHTauTauhisto.Rebin(len(rebin)-1,"",rebin)
		     else:     
			     ggHTauTauhisto.Rebin(rebin)


		     do_binbybin(ggHTauTauhisto,2.04e-06,lowDataBin,highDataBin,0.1,rebin)

		     ggHTauTaufile.cd(folder)
		     ggHTauTauhisto.Write()

		     ggHTauTaufile.Close()

		     
for sign in ['os','ss']:
    for var in commonvars:
	    hist_path=sign+"/"+var[0]
	    jojo= hist_path.split('/')
	    folder= '/'.join(hist_path.split('/')[0:(len(jojo)-1)])
	    datafile=ROOT.TFile(Analyzer+"/data_obs.root")
	    data_histo=datafile.Get(hist_path)
	    lowDataBin = 1
	    highDataBin = data_histo.GetNbinsX()

	    for i in range(1,data_histo.GetNbinsX()+1):
		    if (data_histo.GetBinContent(i) > 0):
			    lowDataBin = i
			    break
	
	    for i in range(data_histo.GetNbinsX(),0,-1):
		    if (data_histo.GetBinContent(i) > 0):
			    highDataBin = i
			    break


	    WGfile=ROOT.TFile(Analyzer+"/WG.root","UPDATE")


	    WGhisto=WGfile.Get(hist_path)

	    do_binbybin(WGhisto,1.25890428e-06,lowDataBin,highDataBin,0.1,1)

	    WGfile.cd(folder)
	    WGhisto.Write()

	    WGfile.Close()

	    STfile=ROOT.TFile(Analyzer+"/T.root","UPDATE")


	    SThisto=STfile.Get(hist_path)

	    do_binbybin(SThisto,5.23465826064e-06,lowDataBin,highDataBin,0.05,1)

	    STfile.cd(folder)
	    SThisto.Write()

	    STfile.Close()


	    
	    TTfile=ROOT.TFile(Analyzer+"/TT.root","UPDATE")


	    TThisto=TTfile.Get(hist_path)

	    do_binbybin(TThisto,1.08709111195e-05,lowDataBin,highDataBin,0.06,1)

	    TTfile.cd(folder)
	    TThisto.Write()

	    TTfile.Close()

	    Dibosonsfile=ROOT.TFile(Analyzer+"/Dibosons.root","UPDATE")


	    Dibosonshisto=Dibosonsfile.Get(hist_path)

	    do_binbybin(Dibosonshisto,1.17948019785e-05,lowDataBin,highDataBin,0.05,1)

	    Dibosonsfile.cd(folder)
	    Dibosonshisto.Write()

	    Dibosonsfile.Close()


	    DYfile=ROOT.TFile(Analyzer+"/DY.root","UPDATE")


	    DYhisto=DYfile.Get(hist_path)

	    do_binbybin(DYhisto,1.1e-05,lowDataBin,highDataBin,0.1,1)

	    DYfile.cd(folder)
	    DYhisto.Write()

	    DYfile.Close()

	    ZTauTaufile=ROOT.TFile(Analyzer+"/ZTauTau.root","UPDATE")


	    ZTauTauhisto=ZTauTaufile.Get(hist_path)

	    do_binbybin(ZTauTauhisto,1.1e-05,lowDataBin,highDataBin,0.1,1)

	    ZTauTaufile.cd(folder)
	    ZTauTauhisto.Write()

	    ZTauTaufile.Close()

	    WJETSMCfile=ROOT.TFile(Analyzer+"/WJETSMC.root","UPDATE")


	    WJETSMChisto=WJETSMCfile.Get(hist_path)

	    do_binbybin(WJETSMChisto,1.5e-04,lowDataBin,highDataBin,0.1,1)

	    WJETSMCfile.cd(folder)
	    WJETSMChisto.Write()

	    WJETSMCfile.Close()

	    vbfHTauTaufile=ROOT.TFile(Analyzer+"/vbfHTauTau.root","UPDATE")


	    vbfHTauTauhisto=vbfHTauTaufile.Get(hist_path)

	    do_binbybin(vbfHTauTauhisto,4.2e-08,lowDataBin,highDataBin,0.1,1)

	    vbfHTauTaufile.cd(folder)
	    vbfHTauTauhisto.Write()

	    vbfHTauTaufile.Close()


	    ggHTauTaufile=ROOT.TFile(Analyzer+"/ggHTauTau.root","UPDATE")


	    ggHTauTauhisto=ggHTauTaufile.Get(hist_path)

	    do_binbybin(ggHTauTauhisto,2.04e-06,lowDataBin,highDataBin,0.1,1)

	    ggHTauTaufile.cd(folder)
	    ggHTauTauhisto.Write()

	    ggHTauTaufile.Close()


