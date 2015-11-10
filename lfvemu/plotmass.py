#from lfvmutau plotter
import os
from sys import argv, stdout, stderr
import ROOT
import sys
from FinalStateAnalysis.PlotTools.MegaBase import make_dirs

ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

jobid = os.environ['jobid']
prod='vbf'

shape_norm = True
if shape_norm == False:
	ynormlabel = "Normalized to Data "
else:
	ynormlabel = "Normalized to 1 "

canvas = ROOT.TCanvas("canvas","canvas",800,800)
LFVStack = ROOT.THStack("stack","")

mydir = 'results/ntupleTest_12Oct/LFVHEMuAnalyzerMVA/'


ggfile125 =  ROOT.TFile(mydir+'VBF_LFV_HToMuTau_M125_13TeV_powheg_pythia8.root')
ggfile120 =  ROOT.TFile(mydir+'VBF_LFV_HToMuTau_M120_13TeV_powheg_pythia8.root')
ggfile130 =  ROOT.TFile(mydir+'VBF_LFV_HToMuTau_M130_13TeV_powheg_pythia8.root')
ggfile200 =  ROOT.TFile(mydir+'VBF_LFV_HToMuTau_M200_13TeV_powheg_pythia8.root')

gglumifile125 = open('inputs/'+jobid+'/VBF_LFV_HToMuTau_M125_13TeV_powheg_pythia8.lumicalc.sum', 'r')
gglumi125 = float(gglumifile125.readline())

gglumifile120 = open('inputs/'+jobid+'/VBF_LFV_HToMuTau_M120_13TeV_powheg_pythia8.lumicalc.sum', 'r')
gglumi120 = float(gglumifile120.readline())

gglumifile130 = open('inputs/'+jobid+'/VBF_LFV_HToMuTau_M130_13TeV_powheg_pythia8.lumicalc.sum', 'r')
gglumi130 = float(gglumifile130.readline())

gglumifile200 = open('inputs/'+jobid+'/VBF_LFV_HToMuTau_M200_13TeV_powheg_pythia8.lumicalc.sum', 'r')
gglumi200 = float(gglumifile200.readline())
jugu=0
juga=0

print "ggfile125=",ggfile125
print

wherehistosare='os/gg/ept30/0'
gendir = ggfile125.Get(wherehistosare)
print "gendir=",gendir
print
hlist = gendir.GetListOfKeys()
print "hlist=",hlist
print
iter = ROOT.TIter(hlist)
print "iter=",iter
print

print ggfile125.GetName()
filepath = 'plots/'+ggfile125.GetName()[8 : 44]+'massplots/'
print "filepath=",filepath
print
startingdir = os.getcwd()
print "startingdir=",startingdir
print
dirs = filepath.split('/')
del dirs[-1]
print "dirs=",dirs
print
thechannel = '#mu #tau_{e}'
for d in dirs:
        currentdir = os.getcwd()
        dirlist = os.listdir(currentdir)
        print "dirlist=",dirlist
        if d in dirlist:
                os.chdir(d)
        else: 
                os.makedirs(d)
                os.chdir(d)
print "dirs=",dirs		
                
os.chdir(startingdir)
        
for i in iter:
	if (i.GetName()=='h_vismass' and jugu==0) or (i.GetName()=='h_collmass_pfmet'and juga==0):
		print "name=",i.GetName()
		name=i.GetName()
		#120_histo = ggfile.Get('gen/'+i.GetName())
		histo_120 = ggfile120.Get('os/gg/ept30/0/'+i.GetName())
		histo_125 = ggfile125.Get(wherehistosare+'/'+i.GetName())
		histo_130 = ggfile130.Get(wherehistosare+'/'+i.GetName())
		histo_200 = ggfile200.Get(wherehistosare+'/'+i.GetName())        
		print
		print "histo120=",histo_120
		print "histo120=",histo_125
		print
		print "120_histo_int=",histo_120.Integral();
		histo_120.Scale(1./gglumi120)
		histo_125.Scale(1./gglumi125)
		histo_130.Scale(1./gglumi130)
		histo_200.Scale(1./gglumi200)
		print
		print "entered plotting loop"
		
		#p = gg_histo.GetName()[0:1]
	#	if p == 't' : p = '#tau'
	#	if p == 'm' : p = '#mu'
	#	if p == 'e' : p = 'e'
		histo_120.Draw("HIST")
		histo_125.Draw("HISTSAME") 
		histo_130.Draw("HISTSAME") 
		histo_200.Draw("HISTSAME") 
		
		histo_125.SetLineColor(2)
		histo_200.SetLineColor(3)
		histo_130.SetLineColor(4)

		histo_120.SetLineWidth(2)
		histo_125.SetLineWidth(2)
		histo_200.SetLineWidth(2)
		histo_130.SetLineWidth(2)

		
		canvas.SetLogy(0)
		variable=''
		
		if name == 'h_vismass' :
			nom='vissmass'
			jugu+=1
			axislabel = "visible mass (GeV)"
			histo_120.GetXaxis().SetTitle(axislabel)
		if name == 'h_collmass_pfmet' : 
			axislabel = "collmass (GeV)"
			histo_120.GetXaxis().SetTitle(axislabel)
			nom='collmass'
			juga+=1
		legend = ROOT.TLegend(0.6,0.75,0.8,0.85)                        
		legend.AddEntry(histo_120, "120 Gev VBF")                
		legend.AddEntry(histo_125, "125 Gev VBF")                
		legend.AddEntry(histo_130, "130 Gev VBF")                
		legend.AddEntry(histo_200, "200 Gev VBF")                
                        
		legend.Draw()
                        
                canvas.Update()
                canvas.SaveAs('/'.join(dirs)+'/'+prod+'_'+nom+'.png')
		
		
