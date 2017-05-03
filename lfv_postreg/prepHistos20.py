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
#
]

commonvars=[('numVertices','number of vertices', 1),
	    ('numGenJets','Number of jets',1),
	    ('NUP','Number of Partons',1),
	    ('jetN_30','Number of jets with p_{T}>30',1)
	    ]
Analyzer='LFVHEMuAnalyzerMVA'+sys.argv[1]+sys.argv[3]+'plot'
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

def do_binbybin(histo,eff_lumi,lowBound,highBound,norm_uncert): #fill empty bins                                                                                                    
        for i in range(1,lowBound):
                if histo.GetBinContent(i) != 0:
                        lowBound = i
                        break
        for i in range(histo.GetNbinsX(),highBound,-1):
                if histo.GetBinContent(i) != 0:
                        highBound = i
                        break
	fillEmptyBins=False
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


for sign in ['os','subtracted/os','subtracted/ss','allfakes/os','allfakes/ss']:
    for var in vars:
        for j in range(2):
            for i in range(4):
		     if j==0:
			     hist_path=sign+"/gg/"+cat_now[i]+"/"+var[0]
		     else:
			     hist_path=sign+"/gg/"+cat_now[i]+"/selected/nosys/"+var[0]

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

		     do_binbybin(WGhisto,8e-07,lowDataBin,highDataBin,0.1)

		     WGfile.cd(folder)
		     WGhisto.Write()

		     WGfile.Close()

		     STfile=ROOT.TFile(Analyzer+"/T.root","UPDATE")


		     SThisto=STfile.Get(hist_path)

		     do_binbybin(SThisto,3.5e-06,lowDataBin,highDataBin,0.1)

		     STfile.cd(folder)
		     SThisto.Write()

		     STfile.Close()


		     
		     TTfile=ROOT.TFile(Analyzer+"/TT.root","UPDATE")


		     TThisto=TTfile.Get(hist_path)

		     do_binbybin(TThisto,8.7e-06,lowDataBin,highDataBin,0.2)

		     TTfile.cd(folder)
		     TThisto.Write()

		     TTfile.Close()

		     Dibosonsfile=ROOT.TFile(Analyzer+"/Dibosons.root","UPDATE")


		     Dibosonshisto=Dibosonsfile.Get(hist_path)

		     do_binbybin(Dibosonshisto,5e-05,lowDataBin,highDataBin,0.1)

		     Dibosonsfile.cd(folder)
		     Dibosonshisto.Write()

		     Dibosonsfile.Close()


		     DYfile=ROOT.TFile(Analyzer+"/DY.root","UPDATE")


		     DYhisto=DYfile.Get(hist_path)

		     do_binbybin(DYhisto,1.1e-05,lowDataBin,highDataBin,0.1)

		     DYfile.cd(folder)
		     DYhisto.Write()

		     DYfile.Close()


		     WJETSMCfile=ROOT.TFile(Analyzer+"/WJETSMC.root","UPDATE")


		     WJETSMChisto=WJETSMCfile.Get(hist_path)

		     do_binbybin(WJETSMChisto,1.5e-04,lowDataBin,highDataBin,0.3)

		     WJETSMCfile.cd(folder)
		     WJETSMChisto.Write()

		     WJETSMCfile.Close()

		     vbfHTauTaufile=ROOT.TFile(Analyzer+"/vbfHTauTau.root","UPDATE")


		     vbfHTauTauhisto=vbfHTauTaufile.Get(hist_path)

		     do_binbybin(vbfHTauTauhisto,4.2e-08,lowDataBin,highDataBin,0.01)

		     vbfHTauTaufile.cd(folder)
		     vbfHTauTauhisto.Write()

		     vbfHTauTaufile.Close()


		     ggHTauTaufile=ROOT.TFile(Analyzer+"/ggHTauTau.root","UPDATE")


		     ggHTauTauhisto=ggHTauTaufile.Get(hist_path)

		     do_binbybin(ggHTauTauhisto,2.04e-06,lowDataBin,highDataBin,0.01)

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

	    do_binbybin(WGhisto,8e-07,lowDataBin,highDataBin,0.1)
		     
	    WGfile.cd(folder)
	    WGhisto.Write()
	    
	    WGfile.Close()

	    STfile=ROOT.TFile(Analyzer+"/T.root","UPDATE")


	    SThisto=STfile.Get(hist_path)
	    
	    do_binbybin(SThisto,3.5e-06,lowDataBin,highDataBin,0.1)

	    STfile.cd(folder)
	    SThisto.Write()

	    STfile.Close()


	    
	    TTfile=ROOT.TFile(Analyzer+"/TT.root","UPDATE")


	    TThisto=TTfile.Get(hist_path)

	    do_binbybin(TThisto,8.7e-06,lowDataBin,highDataBin,0.2)

	    TTfile.cd(folder)
	    TThisto.Write()

	    TTfile.Close()

	    Dibosonsfile=ROOT.TFile(Analyzer+"/Dibosons.root","UPDATE")


	    Dibosonshisto=Dibosonsfile.Get(hist_path)

	    do_binbybin(Dibosonshisto,5e-05,lowDataBin,highDataBin,0.1)

	    Dibosonsfile.cd(folder)
	    Dibosonshisto.Write()

	    Dibosonsfile.Close()


	    DYfile=ROOT.TFile(Analyzer+"/DY.root","UPDATE")


	    DYhisto=DYfile.Get(hist_path)

	    do_binbybin(DYhisto,1.1e-05,lowDataBin,highDataBin,0.1)

	    DYfile.cd(folder)
	    DYhisto.Write()

	    DYfile.Close()


	    WJETSMCfile=ROOT.TFile(Analyzer+"/WJETSMC.root","UPDATE")


	    WJETSMChisto=WJETSMCfile.Get(hist_path)

	    do_binbybin(WJETSMChisto,1.5e-04,lowDataBin,highDataBin,0.3)

	    WJETSMCfile.cd(folder)
	    WJETSMChisto.Write()

	    WJETSMCfile.Close()

	    vbfHTauTaufile=ROOT.TFile(Analyzer+"/vbfHTauTau.root","UPDATE")


	    vbfHTauTauhisto=vbfHTauTaufile.Get(hist_path)

	    do_binbybin(vbfHTauTauhisto,4.2e-08,lowDataBin,highDataBin,0.01)

	    vbfHTauTaufile.cd(folder)
	    vbfHTauTauhisto.Write()

	    vbfHTauTaufile.Close()


	    ggHTauTaufile=ROOT.TFile(Analyzer+"/ggHTauTau.root","UPDATE")


	    ggHTauTauhisto=ggHTauTaufile.Get(hist_path)

	    do_binbybin(ggHTauTauhisto,2.04e-06,lowDataBin,highDataBin,0.01)

	    ggHTauTaufile.cd(folder)
	    ggHTauTauhisto.Write()

	    ggHTauTaufile.Close()




"""
for var in vars:
	for sign in ['ss','os']:
		for j in range(2):
			for i in range(len(cat_now)):
				x=0
				y=0
				if j==0:
					hist_path=sign+"/gg/"+cat_now[i]+"/"+var[0]
				else:
					hist_path= sign+"/gg/"+cat_now[i]+"/selected/nosys/"+var[0]
					self.histomc=None
					self.histodata=None
					self.histoQCD=None
					for filename in os.listdir("results/Oct30/LFVHEMuAnalyzerMVA_highpt"):
						file=ROOT.TFile("results/Oct30/LFVHEMuAnalyzerMVA_highpt/"+filename)
						histo=file.Get(hist_path)
						if "data"  not in filename and "FAKES" not in filename and "LFV" not in filename and "QCD" not in filename:
							if x==0:
								self.histomc=histo.Clone()
								self.histomc.SetDirectory(0)
								x+=1
								else:
	        							self.histomc.Add(histo)
								
	        					elif "data" in filename:      		
	        						if y==0:
	        							y+=1
	        							self.histodata=histo.Clone()
									self.histodata.SetDirectory(0)
	        						else:
	        							self.histodata.Add(histo)
						self.histomc.Scale(12900)				
						print "data",self.histodata.Integral()
						print "MC",self.histomc.Integral()
						self.histoQCD=self.histodata.Clone()
						self.histoQCD.Add(self.histomc,-1)
						self.histoQCD.Scale(1.06)
						print "QCD",self.histoQCD.Integral()
						new_histo=copy.copy(self.histoQCD)
						jojo=hist_path.split('/')
						jojo1='/'.join(jojo[0:(len(jojo)-1)])
						jojo1=jojo1.replace('ss','os',1)
						print jojo1
						self.histos[(jojo1,var[0])]=new_histo
			    
		self.outputfile=ROOT.TFile("QCD.root","recreate")
		self.outputfile.cd()
		for key in self.histos.keys():

#			self.outputfile.cd()
			self.dir0 = self.outputfile.mkdir(key[0])
#			print self.dir0
			self.dir0.Cd("QCD.root:/"+key[0])
#    print dir0
#			print histos[key]
			self.histos[key].SetDirectory(self.dir0)
			self.histos[key].Write()
		self.outputfile.Close()




QCD=GetQCD()
"""
