import rootpy.plotting.views as views
from FinalStateAnalysis.PlotTools.Plotter        import Plotter
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

class myPlotterReco(Plotter):
    def __init__(self):

        self.channel = 'et'
        jobid = os.environ['jobid']
        print "\nPlotting %s for %s\n" % (self.channel, jobid)
        #        period = '7TeV' if '7TeV' in jobid else '8TeV'
        period = '8TeV'
        self.period=period
        sqrts = 7 if '7TeV' in jobid else 8
        self.sqrts = sqrts
        mc_samples = [
            'ggHiggsToETau',
            'vbfHiggsToETau',
            'GluGluToHToTauTau_M-125_8TeV-powheg-pythia6',
            'VBF_HToTauTau_M-125_8TeV-powheg-pythia6',
            'Z*jets_M50',
            'TTJets*',
            'T_t*',
            'Tbar_t*', 
            'Wplus*',
            'WWJets*',
            'WZJets*',
            'ZZJets*',
# 
       ]

        self.files = []
        self.lumifiles = []

        for x in mc_samples:
            self.files.extend(glob.glob('results/%s/LFVHETauAnalyzer/%s.root' % (jobid, x)))
            self.lumifiles.extend(glob.glob('inputs/%s/%s.lumicalc.sum' % (jobid, x)))

        self.outputdir = 'plots/%s/LFVHETauAnalyzer/%s' % (jobid, self.channel.lower())
        if not os.path.exists(self.outputdir):
            os.makedirs(self.outputdir)
    
        blinder = None
            #blind   = 'blind' not in os.environ or os.environ['blind'] == 'YES'
            #print '\n\nRunning Blind: %s\n\n' % blind
        blind = False
        self.blind = blind
        if blind:
            # Don't look at the SS all pass region
            blinder = lambda x: BlindView(x, "ss/p1p2p3/.*")
    

        super(myPlotterReco, self).__init__(self.files, self.lumifiles, self.outputdir,
                                            blinder)
        self.mc_samples = mc_samples

if __name__ == "__main__":
    plotter = myPlotterReco()
    sqrts   = plotter.sqrts
    mc_samples = [
            'ggHiggsToETau',
            'vbfHiggsToETau',
            'GluGluToHToTauTau_M-125_8TeV-powheg-pythia6',
            'VBF_HToTauTau_M-125_8TeV-powheg-pythia6',
            'Z*jets_M50',
            'TTJets*',
            'T_t*',
            'Tbar_t*', 
            'Wplus*',
            'WWJets*',
            'WZJets*',
            'ZZJets*',
# 
        ]


    #plotter.plot('GluGluToHToTauTau_M-125_8TeV-powheg-pythia6', 'os/gg/ept0/ePt', ' ',rebin=1, xaxis='p_{T} (GeV)')
    #plotter.save('ePt')
 
    histoname = ['tPt','tPhi','tEta','ePt','ePhi','eEta','et_DeltaPhi','et_DeltaR','tPFMET_DeltaPhi','tPFMET_Mt','tMVAMET_DeltaPhi','tMVAMET_Mt','ePFMET_DeltaPhi','ePFMET_Mt','eMVAMET_DeltaPhi','eMVAMET_Mt','jetN_20','jetN_30']
    axistitle = ['#tau p_{T} (GeV)','#tau #phi','#tau #eta', 'e p_{T} (GeV)','e #phi','e #eta','e-#tau #Delta#phi','e-#tau #DeltaR','#tau-PFMET #Delta#phi','#tau-PFMET M_{T} (GeV) ','#tau-MVAMET #Delta#phi','#tau-MVAMET M_{T} (GeV)','e-PFMET #Delta#phi','e-PFMET M_{T} (GeV)','e-MVAMET #Delta#phi','e-MVAMET #M_{T} (GeV)','Number of jets','Number of jets']

    foldername = 'os/gg/ept0'

    for n,h in enumerate(histoname) :
        plotter.canvas.SetLogy(True)
        plotter.plot_mc(foldername, 'ggHiggsToETau',h, rebin=1, xaxis= axistitle[n], leftside=False)
        plotter.save('mc_'+h)
