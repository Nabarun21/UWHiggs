import rootpy.plotting.views as views
#from FinalStateAnalysis.PlotTools.SimplePlotter        import SimplePlotter
#from FinalStateAnalysis.PlotTools.Plotter        import Plotter
from BasePlotter import BasePlotter
from FinalStateAnalysis.PlotTools.BlindView      import BlindView
from FinalStateAnalysis.PlotTools.PoissonView    import PoissonView
from FinalStateAnalysis.PlotTools.MedianView     import MedianView
from FinalStateAnalysis.PlotTools.ProjectionView import ProjectionView
#from FinalStateAnalysis.PlotTools.FixedIntegralView import FixedIntegralView
from FinalStateAnalysis.PlotTools.RebinView  import RebinView
from FinalStateAnalysis.MetaData.data_styles import data_styles, colors
from FinalStateAnalysis.PlotTools.decorators import memo
from FinalStateAnalysis.MetaData.datacommon  import br_w_leptons, br_z_leptons
from optparse import OptionParser
import os
import ROOT
import glob
import math
import logging
from fnmatch import fnmatch
from yellowhiggs import xs, br, xsbr

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)

jobid = os.environ['jobid']

print jobid
mc_samples = [
    'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',
    #'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'GluGluHToTauTau_M125_13TeV_powheg_pythia8',
    'GluGlu_LFV_HToETau_M125_13TeV_powheg_pythia8',
    'TT_TuneCUETP8M1_13TeV-amcatnlo-pythia8',
    'VBFHToTauTau_M125_13TeV_powheg_pythia8',
    'VBF_LFV_HToETau_M125_13TeV_powheg_pythia8',
    'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',
    #'WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
    'WW_TuneCUETP8M1_13TeV-pythia8',
    'WZ_TuneCUETP8M1_13TeV-pythia8',
    'ZZ_TuneCUETP8M1_13TeV-pythia8',
]

files = []
lumifiles = []
channel = 'mmt'

files=    glob.glob('results/%s/TauFakeRateAnalyzerMVA/*.root' % (jobid))
lumifiles = glob.glob('inputs/%s/*.lumicalc.sum' % (jobid))

period = '13TeV'
sqrts = 13

def remove_name_entry(dictionary):
    return dict( [ i for i in dictionary.iteritems() if i[0] != 'name'] )

blind   = 'blind' not in os.environ or os.environ['blind'] == 'YES'
print 'blind?', blind
blind_region=[100, 150] if blind else None
blind=False
embedded = False


sign = ['os']
tauiso = ['tLoose', 'tTigh']
outputdir = 'plots/%s/TauFakeRateAnalyzerMVA/' % (jobid)
if not os.path.exists(outputdir):
    os.makedirs(outputdir)


plotter = BasePlotter(files, outputdir, blind_region,use_embedded=embedded)

EWKDiboson = views.StyleView(
    views.SumView( 
        *[ plotter.get_view(regex) for regex in \
          filter(lambda x : x.startswith('WW') or x.startswith('WZ') or x.startswith('ZZ') or x.startswith('WG'), mc_samples )]
#          filter(lambda x : x.startswith('WZ') , mc_samples )]
    ), **remove_name_entry(data_styles['WW*'#,'WZ*', 'WG*', 'ZZ*'
])
)

Wplus = views.StyleView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x :  x.startswith('WJets'),mc_samples)]), **remove_name_entry(data_styles['WJets*']))
DYLL = views.StyleView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x :  x.startswith('DYJets'),mc_samples)]), **remove_name_entry(data_styles['DYJets*']))
#DYLL = views.StyleView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x :  x.endswith('skimmedLL'), mc_samples )]), **remove_name_entry(data_styles['DY']))
#DYTT = views.StyleView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x :  x.endswith('jets_M50_skimmedTT'), mc_samples )]), **remove_name_entry(data_styles['Z*jets*TT']))
#singleT = views.StyleView(views.SumView(  *[ plotter.get_view(regex) for regex in  filter(lambda x : x.startswith('T_') or x.startswith('Tbar_'), mc_samples)]), **remove_name_entry(data_styles['T*_t*']))
SMH = views.StyleView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x : 'HToTauTau' in x , mc_samples)]), **remove_name_entry(data_styles['GluGluHToTauTau*']))
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
new_mc_samples.extend(['EWKDiboson', 'SMH', 'Wplus', 'DYLL', 'TT'
])
#new_mc_samples.extend(['EWKDiboson','DYLL', 'DYTT'])
print new_mc_samples

histoname = [('m1Pt','#mu_{1} p_{T} (GeV)', 2),('m1Phi','#mu_{1} #phi',4),('m1Eta','#mu_{1} #eta',2),('m2Pt','#mu_{2} p_{T} (GeV)',2),('m2Phi','#mu_{2} #phi',4),('m2Eta','#mu_{2} #eta',2),
             ('m1m2Mass', '#mu-#mu Inv Mass (GeV)',1), ('tPt','#tau p_{T} (GeV)',2), ('tPhi','#tau #phi',4),('tEta','#tau #eta',2),('tAbsEta','#tau |#eta|',2),
             ('mtDR', '#mu #tau dR', 2),  ('mtDPhi', '#mu #tau #Delta#phi', 2), ('ztDR', 'Z #tau dR', 2),  ('ztDPhi', 'Z #tau #Delta#phi', 2), ('type1_pfMetEt', 'type1_pfMet', 2), ( 'jetN_30', 'Number of jets, p_{T}>30', 1), ('bjetCSVVeto30', 'Number of b-jets',1) 
]

#rebins = [5, 5, 2, 5, 5, 2, 1, 5, 5, 2, 1]
#rebins = []
#for n in histoname :
#    rebins.append(1)


plotter.mc_samples = new_mc_samples

print plotter.mc_samples

for i in tauiso :
    for s in sign :        
        #for histo, axis in zip(histoname, axistitle):
#        for n,h in enumerate(histoname) :
        foldername = s+'/'+i
        if not os.path.exists(outputdir+foldername):
            os.makedirs(outputdir+foldername)
 
        if i == 'tNoCuts':
            plotter.plot_mc_vs_data(foldername,'CUT_FLOW', rebin=1, xaxis='CUT_FLOW', leftside=False, show_ratio=False, ratio_range=1.5, sort=True)
            plotter.save(foldername+'/CUT_FLOW')
            
        foldername = s+'/'+i
        for h in histoname:
            #plotter.canvas.SetLogy(True)
            #plotter.plot_data(foldername, h, rebin=rebins[n],  xaxis= axistitle[n], leftside=False)
            #plotter.plot_mc_vs_data(foldername,h, rebin=rebins[n], xaxis= axistitle[n], leftside=False, show_ratio=True, ratio_range=1.5, sort=True)
            print 'histoname',  h[0]
            plotter.plot_mc_vs_data_witherrors(foldername,h[0], rebin=h[2], xaxis=h[1], leftside=False, show_ratio=True, ratio_range=1.5, sort=True)
            plotter.save(foldername+'/'+h[0])
            
            
#        foldername = s+'/'+i+'/tptregion'
#        if not os.path.exists(outputdir+foldername):
#            os.makedirs(outputdir+foldername)
       
#        for h in histoname:
#            plotter.plot_mc_vs_data(foldername,h[0], rebin=h[2], xaxis=h[1], leftside=False, show_ratio=True, ratio_range=1.5, sorted=True)
#            plotter.save(foldername+'/'+h[0])


        jets=0
        while jets <  4 :
            foldername = s+'/'+i+'/'+str(int(jets))
            if not os.path.exists(outputdir+foldername):
                os.makedirs(outputdir+foldername)
            for h in histoname:
                plotter.plot_mc_vs_data_witherrors(foldername,h[0], rebin=h[2], xaxis=h[1], leftside=False, show_ratio=True, ratio_range=1.5, sort=True)
                plotter.save(foldername+'/'+h[0])
 #           foldername = s+'/'+i+'/'+str(int(jets))+'/tptregion'
 #           if not os.path.exists(outputdir+foldername):
 #               os.makedirs(outputdir+foldername)
 #           for h in histoname:
 #               plotter.plot_mc_vs_data(foldername,h[0], rebin=h[2], xaxis=h[1], leftside=False, show_ratio=True, ratio_range=1.5, sorted=True)
 #               plotter.save(foldername+'/'+h[0])
            jets+=1
                
