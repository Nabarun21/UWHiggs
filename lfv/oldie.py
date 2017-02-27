#from mauro plotters
import rootpy.plotting.views as views
#Set logging before anything to override rootpy very verbose defaults

import sys
import logging

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

import os
import ROOT
from pdb import set_trace
ROOT.gROOT.LoadMacro("tdrstyle.C")

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
    'WGstarToLNuMuMu_012Jets_13TeV-madgraph',
    'WGstarToLNuEE_012Jets_13TeV-madgraph',
    'WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
    'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1'
    
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
plotter = BasePlotter(files, outputdir, blind_region,use_embedded=embedded)
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


WGamma = views.StyleView(
    views.SumView( 
        *[ plotter.get_view(regex) for regex in 
          filter(lambda x : x.startswith('WG') , mc_samples )]
    ), **remove_name_entry(data_styles['WG*'])
)

SingleT = views.StyleView(
    views.SumView( 
        *[ plotter.get_view(regex) for regex in 
          filter(lambda x : x.startswith('ST') , mc_samples )]
    ), **remove_name_entry(data_styles['ST*'])
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
plotter.views['WGamma']={'view' : WGamma }
plotter.views['SingleT']={'view' :SingleT }


new_mc_samples = []



#print new_sigsamples 
new_mc_samples.extend(['DYLL', 'Wplus','TT','EWKDiboson', 'SMH','WGamma','SingleT'])
#new_mc_samples.extend(['EWKDiboson','DYLL', 'DYTT'])
#new_mc_samples.extend(['EWKDiboson'])


#rebins = [5, 5, 2, 5, 5, 2, 1, 5, 5, 2, 1]
#rebins = []
#for n in histoname :
#    rebins.append(1)


plotter.mc_samples = new_mc_samples





print "break 1"
if not args.no_plots:
   signs = ['os','ss']
   jets = ['0','1','2']
   processtype = ['gg']
   threshold = []
   
   histo_info1 = [
      ('jetN_30', 'number of jets (p_{T} > 30 GeV)', 1)  ,
      ('NUP', 'number of partons', 1)  ,
      ('numGenJets', 'number of gen level jets', 1)  

   ]
       
   histo_info = [
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
      ('scaledmPt', 'p_{T}(mu)/M_{coll}(e#mu)', 2), 
      ('scaledePt', 'p_{T}(e)/M_{coll}(e#mu)', 2)
#      ('mPFMETDeltaPhi_vs_ePFMETDeltaPhi','mPFMETDeltaPhi_vs_ePFMETDeltaPhi',1)
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
      if int(njet)==2:
         plotter.plot_mc_vs_data_witherrors(path,'vbfMass',njet, 1,'vbf dijet Mass',
                                      leftside=False, xrange=(0.,1500.), show_ratio=True, ratio_range=1., 
                                      sort=True,br=50)

         plotter.save(path+"/vbfMass",dotroot=True)
         plotter.plot_mc_vs_data_witherrors(path, 'vbfDeta', njet,2,'#Delta#eta Dijet',
                                      leftside=False, xrange=(0.,10.), show_ratio=True, ratio_range=1., 
                                      sort=True,br=50)
         plotter.save(path+"/vbfDeta",dotroot=True)
      for var, xlabel, rebin in histo_info:
         logging.debug("Plotting %s/%s" % (path, var) )
         plotter.pad.SetLogy(False)
         print var
         if int(njet)==2: 
             rebin = rebin*2
         
         plotter.plot_mc_vs_data_witherrors(path, var, njet,rebin, xlabel,
                                                leftside=False, xrange=(0.,300.), show_ratio=True, ratio_range=1., 
                                                sort=True,br=50)
        
         print "**************************************************************************************************************************************************"
         plotter.save(path+"/"+var,dotroot=True)
     # plotter.set_subdir(os.path.join('embedded', path+'/selected'))if embedded else plotter.set_subdir(path+'/selected')

      for var, xlabel, rebin in histo_info:
          if int(njet)==2:
              rebin = rebin*2
          print var
          logging.debug("Plotting %s/%s" % (path, var) )
          plotter.pad.SetLogy(False)
          plotter.plot_mc_vs_data_witherrors(path+'/selected/nosys', var,njet, rebin, xlabel,
                                             leftside=False, xrange=(0.,300.), show_ratio=True, ratio_range=1., 
                                      sort=True,br=5)
          
          plotter.save(path+"/selected/"+var,dotroot=True)

   for var, xlabel, rebin in histo_info1:
       plotter.plot_mc_vs_data_witherrors('os', var, njet,rebin, xlabel,
                                          leftside=False, xrange=(0.,300.), show_ratio=True, ratio_range=1.,
                                          sort=True,br=50)
       plotter.save("os/gg/"+var,dotroot=True)
 
print "break 2"

#make shapes for limit setting
if not args.no_shapes:
   signal_region = 'os/%s/selected'
   ##signal_region = 'os/gg/ept30/%s'
   jets_names = [
           ('0', 'gg0emu'  , 1),
           ('1', 'boostemu', 1),#was 2
           ('2', 'vbfemu'  , 1),#was 5
   ]
   pjoin = os.path.join
   for njets, cat_name, rebin in jets_names:
      output_path = plotter.base_out_dir
      tfile = ROOT.TFile(pjoin(output_path, 'shapes.%s.root' % njets), 'recreate')
      output_dir = tfile.mkdir(cat_name)
      unc_conf_lines, unc_vals_lines = plotter.write_shapes( 
         signal_region % njets, 'h_collmass_pfmet', output_dir, rebin=rebin,
         br_strenght=1, last=300)
      logging.warning('shape file %s created' % tfile.GetName()) 
      tfile.Close()
      with open(pjoin(output_path, 'unc.%s.conf' % njets), 'w') as conf:
         conf.write('\n'.join(unc_conf_lines))
      with open(pjoin(output_path, 'unc.%s.vals' % njets), 'w') as vals:
         vals.write('\n'.join(unc_vals_lines))

   with open(pjoin(output_path,'.shapes_timestamp'),'w') as stamp:
      stamp.write('no use')

