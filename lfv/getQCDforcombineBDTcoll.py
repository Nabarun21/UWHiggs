import os 
from sys import argv, stdout, stderr
import ROOT
import sys
import math
import copy
ROOT.gROOT.SetStyle("Plain")
cat_now=['0','1','21','22']   #category names in analyzer                                                                                    
                                                                                                                                              
syst_names_now=['jetup','jetdown','tup','tdown','uup','udown']      #sysfolder names in analyzer                                             
vars = [
	("BDT_output", "BDT_output",1)
]

vars2 = [
	("BDT_output", "BDT_output",1)
]

Analyzer="LFVHEMuAnalyzerMVA"+sys.argv[1]
Lumi=int(sys.argv[2])
jobid=sys.argv[3]
class GetQCD(object):
	def __init__(self):
		self.histos={}
                self.histomc=None
	        self.histodata=None
	        self.histoQCD=None
	        for var in vars:
	        	for sign in ['ss']:
	        		for j in [1]:
	        			for i in range(len(cat_now)):
	        				x=0
	        				y=0
	        				if j==0:
							if "BDT" in var[0]:
								hist_path=sign+"/gg/"+cat_now[i]+"/"+var[0]
							else:
								continue
	        				else:
							if "BDT" in var[0]:
								hist_path= sign+"/gg/"+cat_now[i]+"/selected/nosys/"+var[0]
							else:
								continue
						self.histomc=None
						self.histodata=None
						self.histoQCD=None
	        				for filename in os.listdir("results/"+jobid+"/"+Analyzer):
							if "FAKES" in filename or "QCD" in filename: continue
	        					file=ROOT.TFile("results/"+jobid+"/"+Analyzer+"/"+filename)
							histo=file.Get(hist_path)
							print hist_path,"   ",filename,"   ",var[0],"  ",histo.Integral()
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
						self.histomc.Scale(Lumi)				
#						print "data",self.histodata.Integral()
#						print "MC",self.histomc.Integral()
						self.histoQCD=self.histodata.Clone()
						self.histoQCD.Add(self.histomc,-1)
						self.histoQCD.Scale(2.30)
						for bin in range(1,self.histomc.GetNbinsX()+1):
							if self.histomc.GetBinContent(bin)<0.0001:
								self.histoQCD.SetBinContent(bin,0)
						
#						print "QCD",self.histoQCD.Integral()
						new_histo=copy.copy(self.histoQCD)
						jojo=hist_path.split('/')
						jojo1='/'.join(jojo[0:(len(jojo)-1)])
						jojo1=jojo1.replace('ss','os',1)
						self.histos[(jojo1,var[0])]=new_histo
	        for var in vars2:
	        	for sign in ['ss']:
				for i in range(len(cat_now)):
					for k in range(len(syst_names_now)):
						x=0
	        				y=0
						hist_path= sign+"/gg/"+cat_now[i]+"/selected/"+syst_names_now[k]+"/"+var[0]
						self.histomc=None
						self.histodata=None
						self.histoQCD=None
	        				for filename in os.listdir("results/"+jobid+"/"+Analyzer):
							if "FAKES" in filename or "QCD" in filename: continue
	        					file=ROOT.TFile("results/"+jobid+"/"+Analyzer+"/"+filename)
							histo=file.Get(hist_path)
#							print hist_path,"   ",filename,"   ",var[0],"  ",histo.Integral()
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
						self.histomc.Scale(Lumi)				
#						print "data",self.histodata.Integral()
#						print "MC",self.histomc.Integral()
						self.histoQCD=self.histodata.Clone()
						self.histoQCD.Add(self.histomc,-1)
						self.histoQCD.Scale(2.30)
						for bin in range(1,self.histomc.GetNbinsX()+1):
							if self.histomc.GetBinContent(bin)<0.00001:
								self.histoQCD.SetBinContent(bin,0)
#						print "QCD",self.histoQCD.Integral()
						new_histo=copy.copy(self.histoQCD)
						jojo=hist_path.split('/')
						jojo1='/'.join(jojo[0:(len(jojo)-1)])
						jojo1=jojo1.replace('ss','os',1)
						self.histos[(jojo1,var[0])]=new_histo

			    
		self.outputfile=ROOT.TFile("QCDforcombineBDTcoll.root","recreate") 
		self.outputfile.cd()
		for key in self.histos.keys():

#			self.outputfile.cd()
			self.dir0 = self.outputfile.mkdir(key[0])
#			print self.dir0
			self.dir0.Cd("QCDforcombineBDTcoll.root:/"+key[0])
#    print dir0
#			print histos[key]
			self.histos[key].SetDirectory(self.dir0)
			self.histos[key].Write()
		self.outputfile.Close()




QCD=GetQCD()
