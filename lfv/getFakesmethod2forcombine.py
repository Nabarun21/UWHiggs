import os 
from sys import argv, stdout, stderr
import ROOT
import sys
import math
import copy
ROOT.gROOT.SetStyle("Plain")
cat_now=['0','1','21','22']   #category names in analyzer                                                                                    
syst_names_now=['jetup','jetdown','tup','tdown','uup','udown','par0up','par1up','par2up','par3up','par0eeup','par1eeup','par0down','par1down','par2down','par3down','par0eedown','par1eedown'] 
      #sysfolder names in analyzer                                             
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

commonvars=[
	    ]
Analyzer=sys.argv[1]
Lumi=sys.argv[2]
jobid=sys.argv[3]
class GetFAKES(object):
	def __init__(self):
		self.histos={}
                self.histomc=None
	        self.histodata=None
	        self.histoFAKES=None
	        for var in vars:
	        	for sign in ['antiIsolatedweighted/ss','antiIsolatedweighted/os']:
	        		for j in range(2):
	        			for i in range(len(cat_now)):
	        				x=0
	        				y=0
	        				if j==0:
	        					hist_path=sign+"/gg/"+cat_now[i]+"/"+var[0]
	        				else:
	        					hist_path= sign+"/gg/"+cat_now[i]+"/selected/nosys/"+var[0]
						if j!=0 and 'collmass' not in var[0]:
								continue
						self.histomc=None
						self.histodata=None
						self.histoFAKES=None
	        				for filename in os.listdir('LFVHEMuAnalyzerMVA'+Analyzer+Lumi):
							if "FAKES" in filename or "QCD" in filename or "LFV" in filename or "WJETS" in filename: continue
	        					file=ROOT.TFile('LFVHEMuAnalyzerMVA'+Analyzer+Lumi+'/'+filename)
							histo=file.Get(hist_path)
	        					if "data"  not in filename:
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
						self.histomc.Scale(float(Lumi))				
#						print "data",self.histodata.Integral()
#						print "MC",self.histomc.Integral()
						self.histoFAKES=self.histodata.Clone()
						self.histoFAKES.Add(self.histomc,-1)
						for bin in range(1,self.histoFAKES.GetNbinsX()+1):
							if self.histoFAKES.GetBinContent(bin)<0:
                                                                self.histoFAKES.SetBinContent(bin,0.001)
								self.histoFAKES.SetBinError(bin,1.8)
						
#						print "FAKES",self.histoFAKES.Integral()
						new_histo=copy.copy(self.histoFAKES)
						hist_path=hist_path.split('/',1)[1]
						jojo=hist_path.split('/')
						jojo1='/'.join(jojo[0:(len(jojo)-1)])
						self.histos[(jojo1,var[0])]=new_histo
	        for var in vars2:
	        	for sign in ['antiIsolatedweighted/os']:
				for i in range(len(cat_now)):
					for k in range(len(syst_names_now)):
						x=0
	        				y=0
						hist_path= sign+"/gg/"+cat_now[i]+"/selected/"+syst_names_now[k]+"/"+var[0]
						self.histomc=None
						self.histodata=None
						self.histoFAKES=None
	        				for filename in os.listdir('LFVHEMuAnalyzerMVA'+Analyzer+Lumi):
							if "FAKES" in filename or "QCD" in filename or "LFV" in filename or "WJETS" in filename: continue
	        					file=ROOT.TFile('LFVHEMuAnalyzerMVA'+Analyzer+Lumi+'/'+filename)
							histo=file.Get(hist_path)
	        					if "data"  not in filename:
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
						self.histomc.Scale(float(Lumi))				
#						print "data",self.histodata.Integral()
#						print "MC",self.histomc.Integral()
						self.histoFAKES=self.histodata.Clone()
						self.histoFAKES.Add(self.histomc,-1)
						for bin in range(1,self.histoFAKES.GetNbinsX()+1):
							if self.histoFAKES.GetBinContent(bin)<0:
								self.histoFAKES.SetBinContent(bin,0.001)
 								self.histoFAKES.SetBinError(bin,1.8)
						new_histo=copy.copy(self.histoFAKES)
						hist_path=hist_path.split('/',1)[1]
						jojo=hist_path.split('/')
						jojo1='/'.join(jojo[0:(len(jojo)-1)])
						self.histos[(jojo1,var[0])]=new_histo
	        for var in commonvars:
	        	for sign in ['antiIsolatedweighted/os']:
	        		for j in range(1):
	        			for i in range(len(cat_now)):
	        				x=0
	        				y=0
	        				if j==0:
	        					hist_path=sign+"/"+var[0]
	        				else:
	        					hist_path= sign+"/"+var[0]
						self.histomc=None
						self.histodata=None
						self.histoFAKES=None
	        				for filename in os.listdir('LFVHEMuAnalyzerMVA'+Analyzer+Lumi):
							if "FAKES" in filename or "QCD" in filename or "LFV" in filename or "WJETS" in filename: continue
	        					file=ROOT.TFile('LFVHEMuAnalyzerMVA'+Analyzer+Lumi+'/'+filename)
							print filename
							histo=file.Get(hist_path)
	        					if "data"  not in filename:
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
						self.histomc.Scale(float(Lumi))				
#						print "data",self.histodata.Integral()
#						print "MC",self.histomc.Integral()
						self.histoFAKES=self.histodata.Clone()
						self.histoFAKES.Add(self.histomc,-1)
						for bin in range(1,self.histoFAKES.GetNbinsX()+1):
							if self.histoFAKES.GetBinContent(bin)<0:
                                                                self.histoFAKES.SetBinContent(bin,0.001)
								self.histoFAKES.SetBinError(bin,1.8)
						
#						print "FAKES",self.histoFAKES.Integral()
						new_histo=copy.copy(self.histoFAKES)
						hist_path=hist_path.split('/',1)[1]
						jojo=hist_path.split('/')
						jojo1='/'.join(jojo[0:(len(jojo)-1)])
						self.histos[(jojo1,var[0])]=new_histo
			    


			    
		self.outputfile=ROOT.TFile("FAKESmethod2forcombine.root","recreate")
		self.outputfile.cd()
		for key in self.histos.keys():

#			self.outputfile.cd()
			self.dir0 = self.outputfile.mkdir(key[0])
#			print self.dir0
			self.dir0.Cd("FAKESmethod2forcombine.root:/"+key[0])
#    print dir0
#			print histos[key]
			self.histos[key].SetDirectory(self.dir0)
			self.histos[key].Write()
		self.outputfile.Close()




FAKES=GetFAKES()
