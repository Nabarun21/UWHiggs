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
      ('mPt', 'p_{T}(mu) (GeV)', 5),
      ('mEta', 'eta(mu)', 1),
      ('mPhi', 'phi(mu)', 2),
      ('ePt', 'p_{T}(e) (GeV)', 5),
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
      ('pZeta', 'pZeta', 2),
      ('h_collmass_pfmet', 'M_{coll}(e#mu) (GeV)', 1),
      ('mPt', 'p_{T}(mu) (GeV)', 5),
      ('mEta', 'eta(mu)', 1),
      ('mPhi', 'phi(mu)', 2),
      ('ePt', 'p_{T}(e) (GeV)', 5),
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
            ('h_collmass_pfmet','M_{coll}(e#mu) (GeV)', 1),
            ('h_vismass','M_{vis}(e#mu) (GeV)', 1)
            ]

if sys.argv[4]=='BDT2':
	commonvars=[
		('numVertices','number of vertices', 1),
		('numGenJets','Number of jets',1),
		('NUP','Number of Partons',1),
		('BDT_value', 'BDT_value', 1),
		('pZeta', 'pZeta', 2),
		('h_collmass_pfmet', 'M_{coll}(e#mu) (GeV)', 1),
		('mPt', 'p_{T}(mu) (GeV)', 5),
		('mEta', 'eta(mu)', 1),
		('mPhi', 'phi(mu)', 2),
		('ePt', 'p_{T}(e) (GeV)', 5),
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


binning=array.array( 'd', [-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,-0.05,0.00,0.05,0.10,0.15,0.20,0.25,0.30])

binning1jet=array.array( 'd', [-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,-0.05,0.00,0.05,0.10,0.15,0.20,0.25,0.30])


binning2jet=array.array( 'd', [-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,-0.05,0.00,0.05,0.10,0.15,0.20,0.25,0.30])


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


regions= ['os','ss']
regions_common=['os','ss']
try:
	if sys.argv[5]=='QCDCR':
		regions=['antiIsolated/ss','antiIsolated/os']
		regions_common=[]
except IndexError:
	pass
#regions= ['os','ss','fakeRateMethod/os','fakeRateMethod/ss','antiIsolated/os','antiIsolated/ss','antiIsolatedweighted/os','antiIsolatedweighted/ss','antiIsolatedweightedmuon/ss','antiIsolatedweightedmuon/os','antiIsolatedweightedelectron/ss','antiIsolatedweightedelectron/os','antiIsolatedweightedmuonelectron/os','antiIsolatedweightedmuonelectron/ss']:
for sign in regions:
    for var in vars:
        for j in range(2):
            for i in range(4):
		     if j==0:
			     hist_path=sign+"/gg/"+cat_now[i]+"/"+var[0]
		     else:
			     hist_path=sign+"/gg/"+cat_now[i]+"/selected/nosys/"+var[0]
		     if j!=0 and 'collmass' not in var[0] and 'BDT' not in var[0] and 'vismass' not in var[0]:
			     continue
		     if j!=0 and sys.argv[4]=='BDT' and 'collmass' not in var[0]:
			     continue
		     if j!=0 and sys.argv[4]=='BDT2' and 'collmass' in var[0]:
                             continue
		     if j!=0 and sys.argv[4]=='BDT2' and 'vismass' in var[0]:
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

		     rebin=var[2]
		     print hist_path,"     ",rebin
		     if 'collmass' in var[0] or "MtToPfMet" in var[0] or "vismass" in var[0]:
			     if (i==0 ):
				     rebin = 10
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


		     Wfile=ROOT.TFile(Analyzer+"/W.root","UPDATE")


		     Whisto=Wfile.Get(hist_path)
		     


		     if isinstance(rebin,type(binning)):

			     Whisto=Whisto.Rebin(len(rebin)-1,"",rebin)
		     else:     
			     Whisto.Rebin(rebin)

		     do_binbybin(Whisto,1.25890428e-06,lowDataBin,highDataBin,0.1,rebin)
		     
		     Wfile.cd(folder)
#		     print Whisto.GetNbinsX()
		     Whisto.Write()

		     Wfile.Close()

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

		     do_binbybin(TThisto,1.08709111195e-05,lowDataBin,highDataBin,0.1,rebin)

		     TTfile.cd(folder)
		     TThisto.Write()

		     TTfile.Close()

		     Dibosonfile=ROOT.TFile(Analyzer+"/Diboson.root","UPDATE")


		     Dibosonhisto=Dibosonfile.Get(hist_path)
		     if isinstance(rebin,type(binning)):
			     Dibosonhisto=Dibosonhisto.Rebin(len(rebin)-1,"",rebin)
		     else:     
			     Dibosonhisto.Rebin(rebin)

		     do_binbybin(Dibosonhisto,1.17948019785e-05,lowDataBin,highDataBin,0.05,rebin)

		     Dibosonfile.cd(folder)
		     Dibosonhisto.Write()

		     Dibosonfile.Close()


		     Zothersfile=ROOT.TFile(Analyzer+"/Zothers.root","UPDATE")


		     Zothershisto=Zothersfile.Get(hist_path)
		     if isinstance(rebin,type(binning)):
			     Zothershisto=Zothershisto.Rebin(len(rebin)-1,"",rebin)
		     else:     
			     Zothershisto.Rebin(rebin)
 
		     do_binbybin(Zothershisto,1.1e-05,lowDataBin,highDataBin,0.1,rebin)

		     Zothersfile.cd(folder)
		     Zothershisto.Write()

		     Zothersfile.Close()

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


		     qqH_httfile=ROOT.TFile(Analyzer+"/qqH_htt.root","UPDATE")


		     qqH_htthisto=qqH_httfile.Get(hist_path)

		     if isinstance(rebin,type(binning)):
			     qqH_htthisto=qqH_htthisto.Rebin(len(rebin)-1,"",rebin)
		     else:     
			     qqH_htthisto.Rebin(rebin)


		     do_binbybin(qqH_htthisto,4.2e-08,lowDataBin,highDataBin,0.1,rebin)

		     qqH_httfile.cd(folder)
		     qqH_htthisto.Write()

		     qqH_httfile.Close()


		     ggH_httfile=ROOT.TFile(Analyzer+"/ggH_htt.root","UPDATE")


		     ggH_htthisto=ggH_httfile.Get(hist_path)
		     if isinstance(rebin,type(binning)):
			     ggH_htthisto=ggH_htthisto.Rebin(len(rebin)-1,"",rebin)
		     else:     
			     ggH_htthisto.Rebin(rebin)


		     do_binbybin(ggH_htthisto,2.04e-06,lowDataBin,highDataBin,0.1,rebin)

		     ggH_httfile.cd(folder)
		     ggH_htthisto.Write()

		     ggH_httfile.Close()


		     qqH_hwwfile=ROOT.TFile(Analyzer+"/qqH_hww.root","UPDATE")


		     qqH_hwwhisto=qqH_hwwfile.Get(hist_path)

		     if isinstance(rebin,type(binning)):
			     qqH_hwwhisto=qqH_hwwhisto.Rebin(len(rebin)-1,"",rebin)
		     else:     
			     qqH_hwwhisto.Rebin(rebin)


		     do_binbybin(qqH_hwwhisto,4.2e-08,lowDataBin,highDataBin,0.1,rebin)

		     qqH_hwwfile.cd(folder)
		     qqH_hwwhisto.Write()

		     qqH_hwwfile.Close()


		     ggH_hwwfile=ROOT.TFile(Analyzer+"/ggH_hww.root","UPDATE")


		     ggH_hwwhisto=ggH_hwwfile.Get(hist_path)
		     if isinstance(rebin,type(binning)):
			     ggH_hwwhisto=ggH_hwwhisto.Rebin(len(rebin)-1,"",rebin)
		     else:     
			     ggH_hwwhisto.Rebin(rebin)


		     do_binbybin(ggH_hwwhisto,2.04e-06,lowDataBin,highDataBin,0.1,rebin)

		     ggH_hwwfile.cd(folder)
		     ggH_hwwhisto.Write()

		     ggH_hwwfile.Close()


		     
for sign in regions_common:
    for var in commonvars:
	    hist_path=sign+"/"+var[0]
	    jojo= hist_path.split('/')
	    folder= '/'.join(hist_path.split('/')[0:(len(jojo)-1)])

	    lowDataBin = 1
	    highDataBin = 1


	    rebin=var[2]
	    print hist_path,"     ",rebin
	    if 'collmass' in var[0] or "MtToPfMet" in var[0] or "vismass" in var[0]:
		    if (i==0 ):
			    rebin =10
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



	    Wfile=ROOT.TFile(Analyzer+"/W.root","UPDATE")


	    Whisto=Wfile.Get(hist_path)

	    if isinstance(rebin,type(binning)):
		    
		    Whisto=Whisto.Rebin(len(rebin)-1,"",rebin)
	    else:     
		    Whisto.Rebin(rebin)

	    do_binbybin(Whisto,1.25890428e-06,lowDataBin,highDataBin,0.1,1)

	    Wfile.cd(folder)
	    Whisto.Write()

	    Wfile.Close()

	    STfile=ROOT.TFile(Analyzer+"/T.root","UPDATE")


	    SThisto=STfile.Get(hist_path)
	    if isinstance(rebin,type(binning)):
		    SThisto=SThisto.Rebin(len(rebin)-1,"",rebin)
	    else:     
		    SThisto.Rebin(rebin)


	    do_binbybin(SThisto,5.23465826064e-06,lowDataBin,highDataBin,0.05,1)

	    STfile.cd(folder)
	    SThisto.Write()

	    STfile.Close()


	    
	    TTfile=ROOT.TFile(Analyzer+"/TT.root","UPDATE")


	    TThisto=TTfile.Get(hist_path)

	    if isinstance(rebin,type(binning)):
		    TThisto=TThisto.Rebin(len(rebin)-1,"",rebin)
	    else:     
		    TThisto.Rebin(rebin)


	    do_binbybin(TThisto,1.08709111195e-05,lowDataBin,highDataBin,0.1,1)

	    TTfile.cd(folder)
	    TThisto.Write()

	    TTfile.Close()

	    Dibosonfile=ROOT.TFile(Analyzer+"/Diboson.root","UPDATE")


	    Dibosonhisto=Dibosonfile.Get(hist_path)

	    if isinstance(rebin,type(binning)):
		    Dibosonhisto=Dibosonhisto.Rebin(len(rebin)-1,"",rebin)
	    else:     
		    Dibosonhisto.Rebin(rebin)

	    do_binbybin(Dibosonhisto,1.17948019785e-05,lowDataBin,highDataBin,0.05,1)

	    Dibosonfile.cd(folder)
	    Dibosonhisto.Write()

	    Dibosonfile.Close()


	    Zothersfile=ROOT.TFile(Analyzer+"/Zothers.root","UPDATE")


	    Zothershisto=Zothersfile.Get(hist_path)

	    if isinstance(rebin,type(binning)):
		    Zothershisto=Zothershisto.Rebin(len(rebin)-1,"",rebin)
	    else:     
		    Zothershisto.Rebin(rebin)

	    do_binbybin(Zothershisto,1.1e-05,lowDataBin,highDataBin,0.1,1)

	    Zothersfile.cd(folder)
	    Zothershisto.Write()

	    Zothersfile.Close()

	    ZTauTaufile=ROOT.TFile(Analyzer+"/ZTauTau.root","UPDATE")


	    ZTauTauhisto=ZTauTaufile.Get(hist_path)
	    if isinstance(rebin,type(binning)):
		    ZTauTauhisto=ZTauTauhisto.Rebin(len(rebin)-1,"",rebin)
	    else:     
		    ZTauTauhisto.Rebin(rebin)


	    do_binbybin(ZTauTauhisto,1.1e-05,lowDataBin,highDataBin,0.1,1)

	    ZTauTaufile.cd(folder)
	    ZTauTauhisto.Write()

	    ZTauTaufile.Close()

	    qqH_httfile=ROOT.TFile(Analyzer+"/qqH_htt.root","UPDATE")


	    qqH_htthisto=qqH_httfile.Get(hist_path)
	    if isinstance(rebin,type(binning)):
		    qqH_htthisto=qqH_htthisto.Rebin(len(rebin)-1,"",rebin)
	    else:     
		    qqH_htthisto.Rebin(rebin)


	    do_binbybin(qqH_htthisto,4.2e-08,lowDataBin,highDataBin,0.1,1)

	    qqH_httfile.cd(folder)
	    qqH_htthisto.Write()

	    qqH_httfile.Close()


	    ggH_httfile=ROOT.TFile(Analyzer+"/ggH_htt.root","UPDATE")


	    ggH_htthisto=ggH_httfile.Get(hist_path)

	    if isinstance(rebin,type(binning)):
		    ggH_htthisto=ggH_htthisto.Rebin(len(rebin)-1,"",rebin)
	    else:     
		    ggH_htthisto.Rebin(rebin)


	    do_binbybin(ggH_htthisto,2.04e-06,lowDataBin,highDataBin,0.1,1)

	    ggH_httfile.cd(folder)
	    ggH_htthisto.Write()

	    ggH_httfile.Close()



	    qqH_hwwfile=ROOT.TFile(Analyzer+"/qqH_hww.root","UPDATE")


	    qqH_hwwhisto=qqH_hwwfile.Get(hist_path)
	    if isinstance(rebin,type(binning)):
		    qqH_hwwhisto=qqH_hwwhisto.Rebin(len(rebin)-1,"",rebin)
	    else:     
		    qqH_hwwhisto.Rebin(rebin)


	    do_binbybin(qqH_hwwhisto,4.2e-08,lowDataBin,highDataBin,0.1,1)

	    qqH_hwwfile.cd(folder)
	    qqH_hwwhisto.Write()

	    qqH_hwwfile.Close()


	    ggH_hwwfile=ROOT.TFile(Analyzer+"/ggH_hww.root","UPDATE")


	    ggH_hwwhisto=ggH_hwwfile.Get(hist_path)

	    if isinstance(rebin,type(binning)):
		    ggH_hwwhisto=ggH_hwwhisto.Rebin(len(rebin)-1,"",rebin)
	    else:     
		    ggH_hwwhisto.Rebin(rebin)


	    do_binbybin(ggH_hwwhisto,2.04e-06,lowDataBin,highDataBin,0.1,1)

	    ggH_hwwfile.cd(folder)
	    ggH_hwwhisto.Write()

	    ggH_hwwfile.Close()


