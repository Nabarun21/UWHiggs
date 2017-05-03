#from mauro plotters
import rootpy.plotting.views as views
#Set logging before anything to override rootpy very verbose defaults
import sys
import logging
logging.basicConfig(stream=sys.stderr, level=logging.INFO)

import os
import ROOT
from pdb import set_trace
from FinalStateAnalysis.PlotTools.MegaBase import make_dirs
from FinalStateAnalysis.MetaData.data_styles import data_styles
from FinalStateAnalysis.PlotTools.BlindView import BlindView,  blind_in_range
from FinalStateAnalysis.PlotTools.SubtractionView      import SubtractionView, PositiveView
import itertools
import glob
import sys
from BasePlotter import BasePlotter
from argparse import ArgumentParser

def remove_name_entry(dictionary):
    return dict( [ i for i in dictionary.iteritems() if i[0] != 'name'] )



parser = ArgumentParser(description=__doc__)
parser.add_argument('--no-plots', dest='no_plots', action='store_true',
                    default=False, help='Does not print plots')
parser.add_argument('--no-shapes', dest='no_shapes', action='store_true',
                    default=True, help='Does not create shapes for limit computation')
args = parser.parse_args()

#jobid='MiniAODSIMv2-Spring15-25ns_LFV_October13'
jobid = os.environ['jobid']

#jobid = 'MCntuples_3March' 
channel = 'em'
import rootpy.plotting.views as views
        
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)



mc_samples = [
    'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
     'GluGluHToTauTau_M125_13TeV_powheg_pythia8',
    'TT_TuneCUETP8M1_13TeV-powheg-pythia8-evtgen', 
    'VBFHToTauTau_M125_13TeV_powheg_pythia8',
    'WW_TuneCUETP8M1_13TeV-pythia8',
    'WZ_TuneCUETP8M1_13TeV-pythia8',
    'ZZ_TuneCUETP8M1_13TeV-pythia8',
]

print "\nPlotting %s for %s\n" % (channel, jobid)

#check if blind
blind   = 'blind' not in os.environ or os.environ['blind'] == 'YES'
print 'blind?', blind
blind_region=[100, 150] if blind else None
#blind_region=[100, 200] if blind else None

embedded = False
print jobid
files=  glob.glob('results/%s/LFVHEMuAnalyzerMVA/*.root' % (jobid))
#print "files",files
outputdir = 'plots/%s/lfvemu/LFVHEMuAnalyzerMVA/' % (jobid)
plotter = BasePlotter(files, outputdir, blind_region, use_embedded=embedded)
EWKDiboson = views.StyleView(
    views.SumView( 
        *[ plotter.get_view(regex) for regex in 
          filter(lambda x : x.startswith('WW') or x.startswith('WZ') or x.startswith('ZZ') or x.startswith('WG'), mc_samples )]
    ), **remove_name_entry(data_styles['WW*'])
)


Wplus = views.StyleView(
    views.SumView( 
        *[ plotter.get_view(regex) for regex in 
          filter(lambda x : x.startswith('WJets') or x.startswith('W1Jets') or x.startswith('W2Jets') or x.startswith('W3Jets') or x.startswith('W4Jets'), mc_samples )]
    ), **remove_name_entry(data_styles['WJets*'])
)
DYLL = views.StyleView(
    views.SumView( 
        *[ plotter.get_view(regex) for regex in 
          filter(lambda x : x.startswith('DYJets') or x.startswith('DY1Jets') or x.startswith('DY2Jets') or x.startswith('DY3Jets') or x.startswith('DY4Jets'), mc_samples )]
    ), **remove_name_entry(data_styles['DYJets*'])
)
#
#Wplus = views.StyleView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x :  x.startswith('WJets'),mc_samples)]), **remove_name_entry(data_styles['WJets*']))
#DYLL = views.StyleView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x :  x.startswith('DYJets'),mc_samples)]), **remove_name_entry(data_styles['DYJets*']))
#DYLL = views.StyleView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x :  x.endswith('skimmedLL'), mc_samples )]), **remove_name_entry(data_styles['DY']))
#DYTT = views.StyleView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x :  x.endswith('jets_M50_skimmedTT'), mc_samples )]), **remove_name_entry(data_styles['Z*jets*TT']))
#singleT = views.StyleView(views.SumView(  *[ plotter.get_view(regex) for regex in  filter(lambda x : x.startswith('T_') or x.startswith('Tbar_'), mc_samples)]), **remove_name_entry(data_styles['T*_t*']))
SMH = views.StyleView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x : 'HToTauTau' in x , mc_samples)]), **remove_name_entry(data_styles['*HToTauTau*']))
TT = views.StyleView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x : x.startswith('TT_'), mc_samples)]), **remove_name_entry(data_styles['TT*']))


plotter.views['EWKDiboson']={'view' : EWKDiboson }
plotter.views['Wplus']={'view' : Wplus }
plotter.views['DYLL']={'view' : DYLL }
plotter.views['TT']={'view' : TT }
#plotter.views['DYTT']={'view' : DYTT }
#plotter.views['singleT']={'view' : singleT }
plotter.views['SMH']={'view' : SMH }


new_mc_samples = []



#print new_sigsamples 
new_mc_samples.extend(['EWKDiboson', 'SMH','DYLL', 'Wplus','TT'])
#new_mc_samples.extend(['EWKDiboson','DYLL', 'DYTT'])
#new_mc_samples.extend(['EWKDiboson'])


#rebins = [5, 5, 2, 5, 5, 2, 1, 5, 5, 2, 1]
#rebins = []
#for n in histoname :
#    rebins.append(1)


plotter.mc_samples = new_mc_samples





print "break 1"
if not args.no_plots:
   signs = ['os']
   jets = ['0',
      '1',
      '2'
   ]
   processtype = ['gg']
   threshold = []

   histo_info = [
#       ('mPFMET_DeltaPhi', 'Deltaphi-mu-MET (GeV)', 2),
      ('mPFMETDeltaPhi_vs_ePFMETDeltaPhi','mPFMETDeltaPhi_vs_ePFMETDeltaPhi',1)
   #   ('type1_pfMetEt', 'pfMet', 1)      
   ]
   
   logging.debug("Starting plotting")
   


   for sign,  njet in itertools.product(signs,  jets):
      path = os.path.join(sign,'gg',njet)
#      print path
#      plotter.set_subdir(os.path.join('embedded',path)) if embedded else plotter.set_subdir(path)
      if not os.path.exists(outputdir+path):
         os.makedirs(outputdir+path)
         os.makedirs(outputdir+path+'/selected')

      for var, xlabel, rebin in histo_info:
         logging.debug("Plotting %s/%s" % (path, var) )
         plotter.pad.SetLogy(False)
         print var
         if int(njet)==2: 
             rebin = rebin*2
         
         plotter.plot_mc_vs_data_witherrors_2D(path, var, rebin, xlabel,
                                                leftside=False, xrange=(0.,300.), show_ratio=True, ratio_range=1., 
                                                sort=True)
        
         print "**************************************************************************************************************************************************"
         plotter.save(path+"/"+var)
         print "wewedAdasda"
     # plotter.set_subdir(os.path.join('embedded', path+'/selected'))if embedded else plotter.set_subdir(path+'/selected')

