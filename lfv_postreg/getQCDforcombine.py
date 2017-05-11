import array
import os 
from sys import argv, stdout, stderr
import ROOT
import sys
import math
import copy
ROOT.gROOT.SetStyle("Plain")
cat_now=['0','1','21','22']   #category names in analyzer                                                                                    
syst_names_now=[]      #sysfolder names in analyzer                                             
if sys.argv[5]=="True":                                                                                                                                              
	syst_names_now=['uup','udown','mesup','mesdown','eesup','eesdown','eresrhoup','eresrhodown','eresphidown',
			'chargeduesdown','chargeduesup','ecaluesdown','ecaluesup','hcaluesdown','hcaluesup','hfuesdown','hfuesup','puup','pudown',
			'jes_JetAbsoluteFlavMapDown',
			'jes_JetAbsoluteMPFBiasDown',
			'jes_JetAbsoluteScaleDown',
			'jes_JetAbsoluteStatDown',
			'jes_JetFlavorQCDDown',
			'jes_JetFragmentationDown',
			'jes_JetPileUpDataMCDown',
			'jes_JetPileUpPtBBDown',
			'jes_JetPileUpPtEC1Down',
			'jes_JetPileUpPtEC2Down',
			'jes_JetPileUpPtHFDown',
			'jes_JetPileUpPtRefDown',
			'jes_JetRelativeBalDown',
			'jes_JetRelativeFSRDown',
			'jes_JetRelativeJEREC1Down',
			'jes_JetRelativeJEREC2Down',
			'jes_JetRelativeJERHFDown',
			'jes_JetRelativePtBBDown',
			'jes_JetRelativePtEC1Down',
			'jes_JetRelativePtEC2Down',
			'jes_JetRelativePtHFDown',
			'jes_JetRelativeStatECDown',
			'jes_JetRelativeStatFSRDown',
			'jes_JetRelativeStatHFDown',
			'jes_JetSinglePionECALDown',
			'jes_JetSinglePionHCALDown',
			'jes_JetTimePtEtaDown',
			'jes_JetAbsoluteFlavMapUp',
			'jes_JetAbsoluteMPFBiasUp',
			'jes_JetAbsoluteScaleUp',
			'jes_JetAbsoluteStatUp',
			'jes_JetFlavorQCDUp',
			'jes_JetFragmentationUp',
			'jes_JetPileUpDataMCUp',
			'jes_JetPileUpPtBBUp',
			'jes_JetPileUpPtEC1Up',
			'jes_JetPileUpPtEC2Up',
			'jes_JetPileUpPtHFUp',
			'jes_JetPileUpPtRefUp',
			'jes_JetRelativeBalUp',
			'jes_JetRelativeFSRUp',
			'jes_JetRelativeJEREC1Up',
			'jes_JetRelativeJEREC2Up',
			'jes_JetRelativeJERHFUp',
			'jes_JetRelativePtBBUp',
			'jes_JetRelativePtEC1Up',
			'jes_JetRelativePtEC2Up',
			'jes_JetRelativePtHFUp',
			'jes_JetRelativeStatECUp',
			'jes_JetRelativeStatFSRUp',
			'jes_JetRelativeStatHFUp',
			'jes_JetSinglePionECALUp',
			'jes_JetSinglePionHCALUp',
			'jes_JetTimePtEtaUp'
			]      #sysfolder names in analyzer                                             
cutbasedvars = [
	('h_vismass', 'M_{vis} (GeV)', 1),
	('h_collmass_pfmet', 'M_{coll}(e#mu) (GeV)', 1),
#     ('mPt', 'p_{T}(mu) (GeV)', 4),
#     ('mEta', 'eta(mu)', 2),
#     ('mPhi', 'phi(mu)', 4),
#     ('ePt', 'p_{T}(e) (GeV)', 4),
#     ('eEta', 'eta(e)', 2),
#     ('ePhi', 'phi(e)', 4),
#     ('em_DeltaPhi', 'emu Deltaphi', 2),
#     ('em_DeltaR', 'emu Delta R', 2),

#     ('Met', 'MET (GeV)', 5),
#     ('ePFMET_Mt', 'MT-e-MET (GeV)', 5),
#     ('mPFMET_Mt', 'MT-mu-MET (GeV)', 5),
#     ('ePFMET_DeltaPhi', 'Deltaphi-e-MET (GeV)', 2),
#     ('mPFMET_DeltaPhi', 'Deltaphi-mu-MET (GeV)', 2),
#     ('jetN_30', 'number of jets (p_{T} > 30 GeV)', 1),
]

cutbasedvars2 = [
	('h_vismass', 'M_{vis} (GeV)', 1),
	('h_collmass_pfmet', 'M_{coll}(e#mu) (GeV)', 1),
]

BDTvars = [
      ('BDT_value', 'BDT_value', 1),
      ('h_collmass_pfmet', 'M_{coll}(e#mu) (GeV)', 1),
#     ('mPt', 'p_{T}(mu) (GeV)', 4),
#     ('mEta', 'eta(mu)', 2),
#     ('mPhi', 'phi(mu)', 4),
#     ('ePt', 'p_{T}(e) (GeV)', 4),
#     ('eEta', 'eta(e)', 2),
#     ('ePhi', 'phi(e)', 4),
#     ('em_DeltaPhi', 'emu Deltaphi', 2),
#     ('em_DeltaR', 'emu Delta R', 2),
#     ('h_vismass', 'M_{vis} (GeV)', 1),
#     ('Met', 'MET (GeV)', 5),
#     ('ePFMET_Mt', 'MT-e-MET (GeV)', 5),
#     ('mPFMET_Mt', 'MT-mu-MET (GeV)', 5),
#     ('ePFMET_DeltaPhi', 'Deltaphi-e-MET (GeV)', 2),
#     ('mPFMET_DeltaPhi', 'Deltaphi-mu-MET (GeV)', 2),
#     ('jetN_30', 'number of jets (p_{T} > 30 GeV)', 1),
]

BDTvars2 = [
      ('h_collmass_pfmet', 'M_{coll}(e#mu) (GeV)', 1),
]


BDT2vars = [
      ('BDT_value', 'BDT_value', 1),
#     ('mPt', 'p_{T}(mu) (GeV)', 4),
#     ('mEta', 'eta(mu)', 2),
#     ('mPhi', 'phi(mu)', 4),
#     ('ePt', 'p_{T}(e) (GeV)', 4),
#     ('eEta', 'eta(e)', 2),
#     ('ePhi', 'phi(e)', 4),
#     ('em_DeltaPhi', 'emu Deltaphi', 2),
#     ('em_DeltaR', 'emu Delta R', 2),
#     ('h_vismass', 'M_{vis} (GeV)', 1),
#     ('Met', 'MET (GeV)', 5),
#     ('ePFMET_Mt', 'MT-e-MET (GeV)', 5),
#     ('mPFMET_Mt', 'MT-mu-MET (GeV)', 5),
#     ('ePFMET_DeltaPhi', 'Deltaphi-e-MET (GeV)', 2),
#     ('mPFMET_DeltaPhi', 'Deltaphi-mu-MET (GeV)', 2),
#     ('jetN_30', 'number of jets (p_{T} > 30 GeV)', 1),
]



BDT2vars2=[
	('BDT_value', 'BDT_value', 1),
	('h_vismass', 'M_{vis} (GeV)', 1),
	('h_collmass_pfmet', 'M_{coll}(e#mu) (GeV)', 1),
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
commonvars=[]
if sys.argv[4]=='BDT2':
	commonvars=[
		('BDT_value', 'BDT_value', 1),
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
		]
	


binning=array.array( 'd', [-0.6,-0.5,-0.42,-0.34,-0.28,-0.22,-0.18,-0.14,-0.10,-0.06,-0.02,0.02,0.06,0.1,0.14,0.18,0.24,0.30])

binning1jet=array.array( 'd', [-0.6,-0.5,-0.42,-0.34,-0.26,-0.20,-0.16,-0.12,-0.08,-0.04,0.0,0.04,0.08,0.12,0.16,0.22,0.30])


binning2jet=array.array( 'd', [-0.6,-0.5,-0.42,-0.34,-0.28,-0.24,-0.20,-0.16,-0.12,-0.08,-0.04,0.0,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.30])


regions=['ss']
regions_common=['ss']



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
	        	for sign in ['ss']:#,'antiIsolatedweighted/ss','antiIsolated/ss','antiIsolatedweightedmuonelectron/ss','antiIsolatedweightedelectron/ss','antiIsolatedweightedmuon/ss']:
	        		for j in range(2):
	        			for i in range(len(cat_now)):
	        				x=0
	        				y=0
	        				if j==0:
	        					hist_path=sign+"/gg/"+cat_now[i]+"/"+var[0]
	        				else:
	        					hist_path= sign+"/gg/"+cat_now[i]+"/selected/nosys/"+var[0]
						if j!=0 and 'collmass' not in var[0] and 'BDT' not in var[0] and 'vismass' not in var[0]:
							continue
						if j!=0 and sys.argv[4]=='BDT' and 'collmass' not in var[0]:
							continue
						if j!=0 and sys.argv[4]=='BDT2' and 'collmass' in var[0]:
							 continue
						self.histomc=None
						self.histodata=None
						self.histoQCD=None
	        				for filename in os.listdir(Analyzer+str(Lumi)):
							if "FAKES" in filename or "QCD" in filename: continue
	        					file=ROOT.TFile(Analyzer+str(Lumi)+"/"+filename)
							histo=file.Get(hist_path)
	#						print hist_path,"   ",filename,"   ",var[0],"  ",histo.Integral()
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
						if i==2:
							self.histoQCD.Scale(2.86)
						else:
							self.histoQCD.Scale(2.26)

						rebin=var[2]
						if 'collmass' in var[0] or "MtToPfMet" in var[0] or "vismass" in var[0]:
							if (i==0 ):
								rebin=10
							if ( i==1):
								rebin=10
							if ( i==2):
								rebin=25
							if ( i==3):
								rebin=25
						elif "BDT" in var[0]:
							  if (i==0):
								  self.histoQCD=self.histoQCD.Rebin(len(binning)-1,"",binning)
							  elif (i==1):
								  self.histoQCD=self.histoQCD.Rebin(len(binning1jet)-1,"",binning1jet)
							  else:
								  self.histoQCD=self.histoQCD.Rebin(len(binning2jet)-1,"",binning2jet)
						else:
							if (i==2):
								rebin = rebin*2
							if ( i==3 or  i==1 ):
								rebin=rebin*2

						if "BDT" not in var[0]:
							self.histoQCD.Rebin(rebin)		

						lowBound=0
						highBound=0
						for bin in range(1,self.histoQCD.GetNbinsX()):
							if self.histoQCD.GetBinContent(bin) != 0:
								lowBound = bin
								break
						for bin in range(self.histoQCD.GetNbinsX(),0,-1):
							if self.histoQCD.GetBinContent(bin) != 0:
								highBound = bin
								break
						for bin in range(lowBound, highBound+1):
							if lowBound==0:continue
                                                        if self.histoQCD.GetBinContent(bin)<=0:
                                                                self.histoQCD.SetBinContent(bin,0.001)
                                                                self.histoQCD.SetBinError(bin,1.8)
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
						self.histomc=None
						self.histodata=None
						self.histoQCD=None
	        				for filename in os.listdir(Analyzer+str(Lumi)):
							hist_path=sign+"/gg/"+cat_now[i]+"/selected/"+syst_names_now[k]+"/"+var[0]
							if 'data' in filename:
								hist_path=sign+"/gg/"+cat_now[i]+"/selected/nosys/"+var[0]

							if "FAKES" in filename or "QCD" in filename or "LFV" in filename: continue
	        					file=ROOT.TFile(Analyzer+str(Lumi)+"/"+filename)
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
						if i==2:
							self.histoQCD.Scale(2.86)
						else:
							self.histoQCD.Scale(2.26)

						rebin=var[2]
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
								  self.histoQCD=self.histoQCD.Rebin(len(binning)-1,"",binning)
							  elif (i==1):
								  self.histoQCD=self.histoQCD.Rebin(len(binning1jet)-1,"",binning1jet)
							  else:
								  self.histoQCD=self.histoQCD.Rebin(len(binning2jet)-1,"",binning2jet)
						else:
							if (i==2):
								rebin = rebin*2
							if ( i==3 or  i==1 ):
								rebin=rebin*2

						if "BDT" not in var[0]:
							self.histoQCD.Rebin(rebin)		


						lowBound=0
						highBound=0
						for bin in range(1,self.histoQCD.GetNbinsX()):
							if self.histoQCD.GetBinContent(bin) != 0:
								lowBound = bin
								break
						for bin in range(self.histoQCD.GetNbinsX(),0,-1):
							if self.histoQCD.GetBinContent(bin) != 0:
								highBound = bin
								break
						for bin in range(lowBound, highBound+1):
							if lowBound==0:continue
                                                        if self.histoQCD.GetBinContent(bin)<=0:
                                                                self.histoQCD.SetBinContent(bin,0.001)
								self.histoQCD.SetBinError(bin,1.8)
 						new_histo=copy.copy(self.histoQCD)
						jojo=hist_path.split('/')
						jojo1='/'.join(jojo[0:(len(jojo)-1)])
						jojo1=jojo1.replace('ss','os',1)
						self.histos[(jojo1,var[0])]=new_histo
	        for var in commonvars:
	        	for sign in regions_common:
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
						self.histoQCD=None
	        				for filename in os.listdir(Analyzer+str(Lumi)):
							if "FAKES" in filename or "QCD" in filename: continue
	        					file=ROOT.TFile(Analyzer+str(Lumi)+"/"+filename)
#							print filename
							histo=file.Get(hist_path)
		#					print hist_path,"   ",filename,"   ",var[0],"  ",histo.Integral()
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
						if i==2:
							self.histoQCD.Scale(2.86)
						else:
							self.histoQCD.Scale(2.26)

						rebin=var[2]
						if 'collmass' in var[0] or "MtToPfMet" in var[0] or "vismass" in var[0]:
							rebin=10
						elif "BDT" in var[0]:
							self.histoQCD=self.histoQCD.Rebin(len(binning)-1,"",binning)
						else:
							if (i==2 or i==0):
								rebin = rebin*2
							if ( i==3 or  i==1 ):
								rebin=rebin*2

						if "BDT" not in var[0]:
							self.histoQCD.Rebin(rebin)		

						lowBound=0
						highBound=0
						for bin in range(1,self.histoQCD.GetNbinsX()):
							if self.histoQCD.GetBinContent(bin) != 0:
								lowBound = bin
								break
						for bin in range(self.histoQCD.GetNbinsX(),0,-1):
							if self.histoQCD.GetBinContent(bin) != 0:
								highBound = bin
								break
						for bin in range(lowBound, highBound+1):
                                                        if self.histoQCD.GetBinContent(bin)<=0:
                                                                self.histoQCD.SetBinContent(bin,0.001)
                                                                self.histoQCD.SetBinError(bin,1.8)
						new_histo=copy.copy(self.histoQCD)
						jojo=hist_path.split('/')
						jojo1='/'.join(jojo[0:(len(jojo)-1)])
						jojo1=jojo1.replace('ss','os',1)
						self.histos[(jojo1,var[0])]=new_histo
			    
		self.outputfile=ROOT.TFile("QCDforcombine"+sys.argv[1]+".root","recreate")
		self.outputfile.cd()
		for key in self.histos.keys():

#			self.outputfile.cd()
			self.dir0 = self.outputfile.mkdir(key[0])
#			print self.dir0
			self.dir0.Cd("QCDforcombine"+sys.argv[1]+".root:/"+key[0])
#    print dir0
#			print histos[key]
			self.histos[key].SetDirectory(self.dir0)
			self.histos[key].Write()
		self.outputfile.Close()




QCD=GetQCD()