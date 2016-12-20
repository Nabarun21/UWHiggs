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



#parser = ArgumentParser(description=__doc__)
#parser.add_argument('--no-plots', dest='no_plots', action='store_true',
#                    default=False, help='Does not print plots')
#parser.add_argument('--no-shapes', dest='no_shapes', action='store_true',
#                    default=True, help='Does not create shapes for limit computation')
#args = parser.parse_args()

#jobid='MiniAODSIMv2-Spring15-25ns_LFV_October13'
jobid = os.environ['jobid']

#jobid = 'MCntuples_3March' 
channel = 'em'
import rootpy.plotting.views as views
        
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


mc_samples = [
    'DY',
    'WJETSMC',
    'ggHTauTau', 
    'vbfHTauTau',
    'Dibosons',
    'QCD',
    'WG',
    'T',
    'TT'
]

print "\nPlotting %s for %s\n" % (channel, jobid)

#check if blind
blind   = 'blind' not in os.environ or os.environ['blind'] == 'YES'
print 'blind?', blind
blind_region=[100, 150] if blind else None
#blind_region=[100, 200] if blind else None

embedded = False
print jobid

Analyzer='LFVHEMuAnalyzerMVA'+sys.argv[1]+sys.argv[2]+'plot'

print Analyzer
files=  glob.glob(Analyzer+'/*.root')
#print "files",files
outputdir = 'plots/'+sys.argv[3]+'/'+Analyzer+'/' 
plotter = BasePlotter(files, outputdir, blind_region,use_embedded=embedded,blind_path="os/.*ass*")
EWKDiboson = views.StyleView(
    views.SumView( 
        *[ plotter.get_view(regex) for regex in 
          filter(lambda x : x.startswith('Dibo'), mc_samples )]
    ), **remove_name_entry(data_styles['WW*'])
)


Wplus = views.StyleView(
    views.SumView( 
        *[ plotter.get_view(regex) for regex in 
          filter(lambda x : x.startswith('WJETSMC'), mc_samples )]
    ), **remove_name_entry(data_styles['WJets*'])
)


WGamma = views.StyleView(
    views.SumView( 
        *[ plotter.get_view(regex) for regex in 
          filter(lambda x : x.startswith('WG') , mc_samples )]
    ), **remove_name_entry(data_styles['WG*'])
)

QCD = views.StyleView(
    views.SumView( 
        *[ plotter.get_view(regex) for regex in 
          filter(lambda x : x.startswith('QCD') , mc_samples )]
    ), **remove_name_entry(data_styles['QCD*'])
)


SingleT = views.StyleView(
    views.SumView( 
        *[ plotter.get_view(regex) for regex in 
          filter(lambda x : x=='T', mc_samples )]
    ), **remove_name_entry(data_styles['ST*'])
)


DYLL = views.StyleView(
    views.SumView( 
        *[ plotter.get_view(regex) for regex in 
          filter(lambda x : x.startswith('DY') , mc_samples )]
    ), **remove_name_entry(data_styles['DYJets*'])
)


SMH = views.StyleView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x : 'HTauTau' in x , mc_samples)]), **remove_name_entry(data_styles['*HToTauTau*']))
TT = views.StyleView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x : x.startswith('TT'), mc_samples)]), **remove_name_entry(data_styles['TT*']))


plotter.views['EWKDiboson']={'view' : EWKDiboson }
plotter.views['Wplus']={'view' : Wplus }
plotter.views['DYLL']={'view' : DYLL }
plotter.views['TT']={'view' : TT }
plotter.views['QCD']={'view' : QCD }
plotter.views['SMH']={'view' : SMH }
plotter.views['WGamma']={'view' : WGamma }
plotter.views['SingleT']={'view' :SingleT }


new_mc_samples = []



#print new_sigsamples 
new_mc_samples.extend(['WGamma','SMH','SingleT','EWKDiboson','TT','Wplus','DYLL'])
#new_mc_samples.extend(['EWKDiboson','DYLL', 'DYTT'])
#new_mc_samples.extend(['EWKDiboson'])


#rebins = [5, 5, 2, 5, 5, 2, 1, 5, 5, 2, 1]
#rebins = []
#for n in histoname :
#    rebins.append(1)


plotter.mc_samples = new_mc_samples



no_plots=False
if not no_plots:
   signs = ['ss']
   jets = ['0','1','21','22']
#   jets = ['22']
   processtype = ['gg']
   threshold = []
   
   common_histo_info = [
      ('jetN_30', 'number of jets (p_{T} > 30 GeV)', 1)  ,
      ('NUP', 'number of partons', 1)  ,
      ('numGenJets', 'number of gen level jets', 1) ,
      ('numVertices', 'Number of vertices',1),
      ('h_collmass_pfmet', 'M_{coll}(emu) (GeV)', 1),
      ('h_vismass', 'M_{vis} (GeV)', 1),
   ]
       
   histo_info = [

      ('h_collmass_pfmet', 'M_{coll}(emu) (GeV)', 1),
      ('mPt', 'p_{T}(mu) (GeV)', 4), 
      ('mEta', 'eta(mu)', 2),  
      ('mPhi', 'phi(mu)', 4), 
      ('ePt', 'p_{T}(e) (GeV)', 4), 
      ('eEta', 'eta(e)', 2),  
      ('ePhi', 'phi(e)', 4), 
      ('em_DeltaPhi', 'emu Deltaphi', 2), 
      ('em_DeltaR', 'emu Delta R', 2),
      ('h_vismass', 'M_{vis} (GeV)', 1),
      ('Met', 'MET (GeV)', 1),
      ('ePFMET_Mt', 'MT-e-MET (GeV)', 5),
      ('mPFMET_Mt', 'MT-mu-MET (GeV)', 5),
      ('ePFMET_DeltaPhi', 'Deltaphi-e-MET (GeV)', 2),
      ('mPFMET_DeltaPhi', 'Deltaphi-mu-MET (GeV)', 2),
      ('jetN_30', 'number of jets (p_{T} > 30 GeV)', 1),  
####      ('scaledmPt', 'p_{T}(mu)/M_{coll}(emu)', 2), 
##      ('scaledePt', 'p_{T}(e)/M_{coll}(emu)', 2)
##      ('mPFMETDeltaPhi_vs_ePFMETDeltaPhi','mPFMETDeltaPhi_vs_ePFMETDeltaPhi',1)
   #   ('type1_pfMetEt', 'pfMet', 1)      
   ]
   
   logging.debug("Starting plotting")
   
   if not os.path.exists(outputdir+"/ss"):
         os.makedirs(outputdir+"/ss")


   for var, xlabel, rebin in common_histo_info:
         plotter.plot_mc_vs_data_new('ss', var, '0',rebin, xlabel,
                                                leftside=False, xrange=(-10.0,300.0), show_ratio=True, ratio_range=1., 
                                                sort=True,br=50)
         plotter.save("ss/"+var,dotroot=False)


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
         if (int(njet)==21 ):
             rebin = rebin*2
         if ( int(njet)==22 or int(njet)==112):
             rebin=rebin*2
         plotter.plot_mc_vs_data_new(path, var, njet,rebin, xlabel,
                                                leftside=False, xrange=(-10.0,300.0), show_ratio=True, ratio_range=1., 
                                                sort=True,br=50)
        
         print "**************************************************************************************************************************************************"
         plotter.save(path+"/"+var,dotroot=False)
     # plotter.set_subdir(os.path.join('embedded', path+'/selected'))if embedded else plotter.set_subdir(path+'/selected')

      for var, xlabel, rebin in histo_info:
          if "mass" not in var:continue
          if (int(njet)==21):
              rebin = rebin*2
          if ( int(njet)==22 or  int(njet)==1 ):
             rebin=rebin*2
          logging.debug("Plotting %s/%s" % (path, var) )
          plotter.pad.SetLogy(False)
          plotter.plot_mc_vs_data_new(path+'/selected/nosys', var,njet, rebin, xlabel,
                                             leftside=False, xrange=(-10.,300.0), show_ratio=True, ratio_range=1., 
                                      sort=True,br=5)
          
          plotter.save(path+"/selected/"+var,dotroot=False)



#make shapes for limit setting
no_shapes=True
if not no_shapes:
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


