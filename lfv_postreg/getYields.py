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
      ('h_collmass_pfmet', 'M_{coll}(e#mu) (GeV)', 1),
]

regions=[('os','Signal(isolated OS) region using MC to estimated fakes'),('ss','Isolated SS region using MC to estimate Fakes'),('allfakes/ss','anti-isolated SS region'),('allfakes/os','anti-isolated OS region'),('subtracted/os','Signal(isolated OS) region using Fakes from data'),('subtracted/ss','Isolated SS region using fakes from data')]

class getYield(object):
	def __init__(self):
	        for var in vars:
	        		for j in range(2):
	        			for i in range(len(cat_now)):
						for sign in [regions[0]]:
							if j==0:
								print "PRINTING PRESELECTION YIELDS for ",i," jets"
								hist_path=sign[0]+"/gg/"+cat_now[i]+"/"+var[0]
							else:
								print "PRINTING Final selection Yields for",i, "jets"

								hist_path= sign[0]+"/gg/"+cat_now[i]+"/selected/nosys/"+var[0]
							print "printing yields for ",sign[1]
							for filename in os.listdir("results/Oct30/LFVHEMuAnalyzerMVAv1"):
								if "FAKES" in filename  and 'subtracted' not in sign[0]: continue
								if "QCD" in filename  and 'subtracted' in sign[0]: continue
								if "QCD" in filename  and 'ss' in sign[0]: continue
								file=ROOT.TFile("results/Oct30/LFVHEMuAnalyzerMVAv1/"+filename)
								histo=file.Get(hist_path)
								if "data" not in filename and "FAKES" not in filename and "QCD" not in filename:
									histo.Scale(18937)
								if "data" not in filename and "LFV" in filename:
									print filename,"  :  :",histo.Integral()
							print "                      "		
							print "                      "		




QCD=getYield()
