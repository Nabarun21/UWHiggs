'''

Plot the MMT channel

Usage: python WHPlotterMMT.py

'''

import glob
import logging
import os
import ROOT
import sys
import WHPlotterBase
from WHPlotterBase import make_styler, parser, project_x
import rootpy.plotting.views as views
from FinalStateAnalysis.MetaData.data_styles import data_styles, colors

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

class WHPlotterMMT(WHPlotterBase.WHPlotterBase):
    def __init__(self):
        super(WHPlotterMMT, self).__init__('MMT')

rebin_slim = range(20, 81, 10)+[100, 130, 300]
categories = {
    'LTCut' : [80, 650],
    'Full'  : [0,  650],
    'LTLow' : [0,  130],
    'LTHigh': [130, 650],
}



if __name__ == "__main__":
    plotter = WHPlotterMMT()
    sqrts   = plotter.sqrts
    options,NOTUSED = parser.parse_args()
    if not options.dry_run:
        ###########################################################################
        ##  Zmm control plots #####################################################
        ###########################################################################
        plotter.set_subdir('mc_data/zmm')
        
        plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'm1_m2_Mass', rebin=10, 
                                xaxis='m_{#mu#mu} (GeV)', leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('mcdata-Zmm-m1m2Mass')

        plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'm1Pt#LT', rebin=10, 
                                xaxis='p_{T}^{#mu1} (GeV)', leftside=False, preprocess=project_x)
        plotter.add_cms_blurb(sqrts)
        plotter.save('mcdata-Zmm-m1Pt')

        plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'm2Pt#LT', rebin=10, 
                                xaxis='p_{T}^{#mu1} (GeV)', leftside=False, preprocess=project_x)
        plotter.add_cms_blurb(sqrts)
        plotter.save('mcdata-Zmm-m2Pt')

        plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'm1AbsEta#LT', rebin=10, 
                                xaxis='p_{T}^{#mu1} (GeV)', leftside=False, preprocess=project_x)
        plotter.add_cms_blurb(sqrts)
        plotter.save('mcdata-Zmm-m1AbsEta')

        plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'm2AbsEta#LT', rebin=10, 
                                xaxis='p_{T}^{#mu1} (GeV)', leftside=False, preprocess=project_x)
        plotter.add_cms_blurb(sqrts)
        plotter.save('mcdata-Zmm-m2AbsEta')

        # Check PU variables
        plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'rho')
        plotter.save('mcdata-Zmm-rho')
        
        plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'nvtx')
        plotter.save('mcdata-Zmm-nvtx')


        ###########################################################################
        ##  WZ control plots #####################################################
        ###########################################################################
        plotter.set_subdir('mc_data/wz_enhanced')
        
        plotter.plot_mc_vs_data('ss/tau_os/p1p2p3_enhance_wz', 'm2_t_Mass#LT', xaxis='m_{#mu#mu} (GeV)', 
                                xrange=(0, 120), rebin=10, leftside=False, preprocess=project_x)
        plotter.add_cms_blurb(sqrts)
        plotter.save('mcdata-ss-p1p2p3-enhance_wz-subMass')
        
        plotter.plot_mc_vs_data('ss/tau_os/p1p2p3_enhance_wz', 'm1_t_Mass', xaxis='m_{#mu#mu} (GeV)', 
                                xrange=(0, 120), rebin=10, leftside=False, preprocess=project_x)
        plotter.add_cms_blurb(sqrts)
        plotter.save('mcdata-ss-p1p2p3-enhance_wz-leadMass')
        
        #plotter.plot_mc_vs_data('ss/p1f2p3_enhance_wz', 'm1_t_Mass', xaxis='m_{#mu#mu} (GeV)', xrange=(0, 120), rebin=10)
        #plotter.add_cms_blurb(sqrts)
        #plotter.save('mcdata-ss-p1f2p3-enhance_wz-leadMass')
        #
        #plotter.plot_mc_vs_data('ss/p1f2p3_enhance_wz/w2', 'm1_t_Mass', xaxis='m_{#mu#mu} (GeV)', xrange=(0, 120), rebin=10)
        #plotter.add_cms_blurb(sqrts)
        #plotter.save('mcdata-ss-p1f2p3-w2-enhance_wz-leadMass')
        #
        #plotter.plot_mc_vs_data('ss/p1p2p3_enhance_wz', 'subMTMass', xaxis='m_{#mu#mu} (GeV)', xrange=(0, 120), rebin=10)
        #plotter.add_cms_blurb(sqrts)
        #plotter.save('mcdata-ss-p1p2p3-enhance_wz-subMTMass')
        
        plotter.set_subdir('WZ_enhanced')

        plotter.plot_final_wz('m1_t_Mass', 10, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-wz-leadMass')

        plotter.plot_final_wz('m2_t_Mass#LT', 10, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)',
                              project=[0,650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-wz-subMass')

        plotter.plot_final_wz('m2Pt#LT', 5, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)',
                              project=[0,650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-wz-m2Pt')
        
        plotter.plot_final_wz('m2JetPt#LT', 5, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)',
                              project=[0,650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-wz-m2JetPt')


        ###########################################################################
        ##  Signal region plots    ################################################
        ###########################################################################
        plotter.set_subdir('')
        #rebin_slim = range(20, 91, 10)+[110,200]
                
        for label, proj_range in categories.iteritems():
            for tau_charge in ['tau_os', 'tau_ss']:
                if tau_charge == 'tau_os':
                    plotter.set_subdir('%s' % label)
                else:
                    plotter.set_subdir('%s_charge3' % label)

                plotter.plot_final('m2_t_Mass#LT', rebin_slim, xaxis='m_{#mu_{2}#tau} (GeV)', maxy=None, 
                                   project=proj_range, project_axis='X', differential=True, 
                                   yaxis='Events / GeV', show_error=True, tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-subMass-%s' % label)

                plotter.plot_final('m2_t_Mass#LT', 20, xaxis='m_{#mu_{2}#tau} (GeV)', 
                                   project=proj_range, project_axis='X', 
                                   yaxis='Events', show_error=True, x_range=[0,199], 
                                   maxy=15.6, tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-subMass-const_bin-%s' % label, verbose=True)

                #data_view = plotter.make_signal_views(unblinded=(not plotter.blind), 
                #                                     qcd_weight_fraction=0.5,
                #                                     rebin=20, project=proj_range, project_axis='X')['data']
                #data_hist = data_view.Get('m2_t_Mass#LT')
                #data_hist.GetXaxis().SetRangeUser(0,200)
                #data_hist.GetXaxis().SetRangeUser(0,15.6)
                #data_hist.SetMarkerColor(2)
                #data_hist.Draw()
                #plotter.add_cms_blurb(sqrts)
                #plotter.save('bareData-subMass-const_bin-%s' % label)

                plotter.plot_final('m2_t_Mass#LT', 300, xaxis='m_{#mu_{2}#tau} (GeV)', maxy=None, 
                                   project=proj_range, project_axis='X', show_error=True, tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-subMass-%s-counting' % label) #, dotc=True, dotroot=True)
                
                #pt
                plotter.plot_final("m1Pt#LT"    , 10, xaxis='p_{T#mu_{1}} (GeV)', 
                                   maxy=None, project=proj_range, project_axis='X', 
                                   show_error=True, tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-m1Pt-%s' % label)
                
                plotter.plot_final("m2Pt#LT"    , 10, xaxis='p_{T#mu_{2}} (GeV)', maxy=None, 
                                   project=proj_range, project_axis='X', 
                                   show_error=True, tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-m2Pt-%s' % label)
                
                plotter.plot_final("m1JetPt#LT" , 10, xaxis='p_{T Jet #mu_{1}} (GeV)', 
                                   maxy=None, project=proj_range, project_axis='X', 
                                   show_error=True, tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-m1JetPt-%s' % label)
                
                plotter.plot_final("m2JetPt#LT" , 10, xaxis='p_{T Jet #mu_{2}} (GeV)', 
                                   maxy=None, project=proj_range, project_axis='X', 
                                   show_error=True, tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-m2JetPt-%s' % label)
                
                plotter.plot_final("tPt#LT"    , 10, xaxis='p_{T#tau} (GeV)', maxy=None, 
                                   project=proj_range, project_axis='X', show_error=True, 
                                   tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-tPt-%s' % label)
                
                #eta
                plotter.plot_final("m1AbsEta#LT", 10, xaxis='|#eta_{#mu_{1}}|', maxy=None, 
                                   project=proj_range, project_axis='X', show_error=True, 
                                   tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-m1AbsEta-%s' % label)
                
                plotter.plot_final("m2AbsEta#LT", 10, xaxis='|#eta_{#mu_{2}}|'  , maxy=None, 
                                   project=proj_range, project_axis='X', show_error=True, 
                                   tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-m2AbsEta-%s' % label)
                
                plotter.plot_final("tAbsEta#LT", 10, xaxis='|#eta_{#tau}|'  , maxy=None, 
                                   project=proj_range, project_axis='X', show_error=True, 
                                   tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-tAbsEta-%s' % label)
                
                #DR
                plotter.plot_final("m1_t_DR#LT", 10, xaxis='#DeltaR_{#mu_{1}#tau}', maxy=None, 
                                   project=proj_range, project_axis='X', show_error=True, 
                                   tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-m1_t_DR-%s' % label)
                
                plotter.plot_final("m2_t_DR#LT", 10, xaxis='#DeltaR_{#mu_{2}#tau}'  , maxy=None, 
                                   project=proj_range, project_axis='X', show_error=True, 
                                   tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-m2_t_DR-%s' % label)
                
                #Jet BTag
                
                #from pdb import set_trace; set_trace()
                plotter.plot_final("m2JetBtag#LT", [-100, -40, 20, 100], xaxis='#mu_{2} Jet Btag'  , 
                                   maxy=None, project=proj_range, project_axis='X', show_error=True, 
                                   tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-m2JetBtag-%s' % label, dotroot=True)
                
                plotter.plot_final("m1JetBtag#LT", 1, xaxis='#mu_{1} Jet Btag'  , maxy=None, 
                                   project=proj_range, project_axis='X', show_error=True, 
                                   tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-m1JetBtag-%s' % label)
                
                plotter.plot_final("m2JetCSVBtag#LT", 1, xaxis='#mu_{2} Jet CSV Btag'  , maxy=None, 
                                   project=proj_range, project_axis='X', x_range=[0,1], show_error=True, 
                                   tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-m2JetCSVBtag-%s' % label, dotroot=True)
                
                plotter.plot_final("m1JetCSVBtag#LT", 1, xaxis='#mu_{1} Jet CSV Btag'  , maxy=None, 
                                   project=proj_range, project_axis='X', x_range=[0,1], show_error=True, 
                                   tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-m1JetCSVBtag-%s' % label)

        #from pdb import set_trace; set_trace()
        plotter.plot_final('m2_t_Mass#LT', rebin_slim, xaxis='m_{#mu_{2}#tau} (GeV)', maxy=None, 
                           project=[0, 130], project_axis='X', differential=True, yaxis='Events / GeV')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-LTLow')

        plotter.plot_final('m2_t_Mass#LT', rebin_slim, xaxis='m_{#mu_{2}#tau} (GeV)', maxy=None, project=[130, 650], project_axis='X', differential=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-LTHigh')

        plotter.plot_final('m2_t_Mass#LT', rebin_slim, xaxis='m_{#mu_{2}#tau} (GeV)', maxy=None, project=[80, 650], project_axis='X', differential=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-LTCut')

        plotter.plot_final('m2_t_Mass#LT', 20, xaxis='m_{#mu_{2}#tau} (GeV)', maxy=None, project=[80, 650], project_axis='X', differential=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-LTCut-flatbin')

        plotter.plot_final('LT', 5, xaxis='LT (GeV)', maxy=15)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-LT')

        plotter.plot_final('m1Pt#LT', 10, project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-m1Pt')

        plotter.plot_final('m2Pt#LT', 10, project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-m2Pt')

        plotter.plot_final('m2Pt#LT', 10, project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-m2Pt')

        plotter.plot_final('m1AbsEta#LT', 10, project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-m1AbsEta')

        plotter.plot_final('m2AbsEta#LT', 10, project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-m2AbsEta')

        plotter.plot_final('m2AbsEta#LT', 10, project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-m2AbsEta')

        plotter.plot_final('m2_t_Mass#LT', 20, xaxis='m_{#mu_{2}#tau} (GeV)', project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass')

        plotter.plot_final('m2_t_Mass#LT', 20, xaxis='m_{#mu_{2}#tau} (GeV)', qcd_weight_fraction=1, project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-qweight-subMass')

        plotter.plot_final('m2_t_Mass#LT', 20, xaxis='m_{#mu_{2}#tau} (GeV)', qcd_weight_fraction=0.5, project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-qweight05-subMass')

        plotter.plot_final('m2_t_Mass#LT', 20, xaxis='m_{#mu_{2}#tau} (GeV)', show_error=True, project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-werror')

        plotter.plot_final('m2_t_Mass#LT', 20, xaxis='m_{#mu_{2}#tau} (GeV)',
                           show_error=True, fake_error=0, wz_error=0, zz_error=0
                           , project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-wshapeerror')

        plotter.plot_final('m2RelPFIsoDB', 10)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-m2Iso')

        ###########################################################################
        ##  F3 enhanced region plots    ###########################################
        ###########################################################################
                
        for label, proj_range in categories.iteritems():
            for tau_charge in ['tau_os', 'tau_ss']:
                if tau_charge == 'tau_os':
                    plotter.set_subdir('f3/%s' % label)
                else:
                    plotter.set_subdir('f3/%s_charge3' % label)
                plotter.plot_final_f3('m2_t_Mass#LT', rebin_slim, xaxis='m_{#mu_{2}#tau} (GeV)', 
                                      maxy=None, project=proj_range, project_axis='X', differential=True,
                                      show_error=True, tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-f3-subMass-%s' % label)

                plotter.plot_final_f3('m2_t_Mass#LT', 300, xaxis='m_{#mu_{2}#tau} (GeV)', 
                                      maxy=None, project=proj_range, project_axis='X',
                                      tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-f3-subMass-%s-counting' % label, dotc=True, dotroot=True)

                #pt
                plotter.plot_final_f3("m1Pt#LT"    , 10, xaxis='p_{T#mu_{1}} (GeV)', maxy=None, 
                                      project=proj_range, project_axis='X', tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-f3-m1Pt-%s' % label)

                plotter.plot_final_f3("m2Pt#LT"    , 10, xaxis='p_{T#mu_{2}} (GeV)', maxy=None, 
                                      project=proj_range, project_axis='X', tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-f3-m2Pt-%s' % label)

                plotter.plot_final_f3("m1JetPt#LT" , 10, xaxis='p_{T Jet #mu_{1}} (GeV)', 
                                      maxy=None, project=proj_range, project_axis='X', 
                                      tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-f3-m1JetPt-%s' % label)

                plotter.plot_final_f3("m2JetPt#LT" , 10, xaxis='p_{T Jet #mu_{2}} (GeV)', 
                                      maxy=None, project=proj_range, project_axis='X', 
                                      tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-f3-m2JetPt-%s' % label)

                plotter.plot_final_f3("tPt#LT", 10, xaxis='p_{T#tau} (GeV)', maxy=None, 
                                      project=proj_range, project_axis='X', tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-f3-tPt-%s' % label)

                #eta
                plotter.plot_final_f3("m1AbsEta#LT", 10, xaxis='|#eta_{#mu_{1}}|', maxy=None, 
                                      project=proj_range, project_axis='X', tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-f3-m1AbsEta-%s' % label)

                plotter.plot_final_f3("m2AbsEta#LT", 10, xaxis='|#eta_{#mu_{2}}|', maxy=None, 
                                      project=proj_range, project_axis='X', tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-f3-m2AbsEta-%s' % label)

                plotter.plot_final_f3("tAbsEta#LT", 10, xaxis='|#eta_{#tau}|'  , maxy=None, 
                                      project=proj_range, project_axis='X', tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-f3-tAbsEta-%s' % label)

                #DR
                plotter.plot_final_f3("m1_t_DR#LT", 5, xaxis='#DeltaR_{#mu_{1}#tau}', 
                                      maxy=None, project=proj_range, project_axis='X',
                                      tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-f3-m1_t_DR-%s' % label)

                plotter.plot_final_f3("m2_t_DR#LT", 5, xaxis='#DeltaR_{#mu_{2}#tau}'  , maxy=None, 
                                      project=proj_range, project_axis='X', tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-f3-m2_t_DR-%s' % label)

                #Jet BTag
            
                #from pdb import set_trace; set_trace()
                plotter.plot_final_f3("m2JetBtag#LT", [-100, -40, 20, 100], xaxis='#mu_{2} Jet Btag', 
                                      maxy=None, project=proj_range, project_axis='X', 
                                      tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-f3-m2JetBtag-%s' % label, dotroot=True)

                plotter.plot_final_f3("m1JetBtag#LT", 1, xaxis='#mu_{1} Jet Btag'  , maxy=None, 
                                      project=proj_range, project_axis='X', tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-f3-m1JetBtag-%s' % label)

                plotter.plot_final_f3("m2JetCSVBtag#LT", 1, xaxis='#mu_{2} Jet CSV Btag', maxy=None, 
                                      project=proj_range, project_axis='X', x_range=[0,1], 
                                      tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-f3-m2JetCSVBtag-%s' % label, dotroot=True)

                plotter.plot_final_f3("m1JetCSVBtag#LT", 1, xaxis='#mu_{1} Jet CSV Btag', maxy=None, 
                                      project=proj_range, project_axis='X', x_range=[0,1], 
                                      tau_charge=tau_charge)
                plotter.add_cms_blurb(sqrts)
                plotter.save('final-f3-m1JetCSVBtag-%s' % label)


        plotter.set_subdir('f3')
        plotter.plot_final_f3('LT', 5, xaxis='LT (GeV)', qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-LT')

        plotter.plot_final_f3('m1_t_Mass', 10, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)', qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-leadMass')

        plotter.plot_final_f3('m2_t_Mass#LT', 20, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)', qcd_correction=False, qcd_weight_fraction=0, show_error=True, project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-wjetfake-subMass')

        plotter.plot_final_f3('m2_t_Mass#LT', 20, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)', qcd_correction=False, qcd_weight_fraction=1, show_error=True, project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-qcdfake-subMass')

        plotter.plot_final_f3("m2JetBtag#LT", 4, xaxis='#mu_{2} Jet Btag'  , maxy=None, project=[0, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-m2JetBtag')

        plotter.plot_final_f3("m1JetBtag#LT", 4, xaxis='#mu_{1} Jet Btag'  , maxy=None, project=[0, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-m1JetBtag')

        plotter.plot_final_f3('m2_t_Mass#LT', 20, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)', qcd_correction=False, qcd_weight_fraction=0.5, show_error=True, project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-subMass')

        plotter.plot_final_f3('m2_t_Mass#LT', 300, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)', qcd_correction=False, qcd_weight_fraction=0.5, show_error=True, project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-subMass-counting-like', dotc=True, dotroot=True)

        #plotter.plot_final_f3_split('m2_t_Mass#LT', 10, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)', project=[130, 650], project_axis='X')
        #plotter.add_cms_blurb(sqrts)
        #plotter.save('final-f3-split-subMass')

        plotter.plot_final_f3('m1Pt#LT', 5, xaxis='p_{T#mu_{1}} (GeV)', project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-m1Pt')

        plotter.plot_final_f3('m1JetPt#LT', 5, xaxis='p_{TJet#mu_{1}} (GeV)', project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-m1JetPt')

        plotter.plot_final_f3('m1AbsEta#LT', 10, xaxis='|#eta_{#mu_{1}}| (GeV)', project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-m1AbsEta')

        plotter.plot_final_f3('m2Pt#LT', 5, xaxis='p_{T#mu_{2}} (GeV)', project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-m2Pt')

        plotter.plot_final_f3('m2JetPt#LT', 5, xaxis='p_{TJet#mu_{2}} (GeV)', project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-m2JetPt')

        plotter.plot_final_f3('m2AbsEta#LT', 10, xaxis='|#eta_{#mu_{2}}| (GeV)', project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-m2AbsEta')

        plotter.plot_final_f3('m1_t_DR#LT', 5, xaxis='#DeltaR_{#mu_{1}#tau}', project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-m1_t_DR')

        plotter.plot_final_f3('m2_t_DR#LT', 5, xaxis='#DeltaR_{#mu_{2}#tau}', project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-m2_t_DR')

        ###########################################################################
        ##  Check QCD contamination in control regions ############################
        ###########################################################################
        plotter.set_subdir('qcd_contamination')

        #plotter.plot_qcd_contamination('m2_t_Mass', 2, 10)
        #plotter.add_cms_blurb(sqrts)
        #plotter.save('final-qcd2-subMass')
        #
        #plotter.plot_qcd_contamination('m2_t_Mass', 1, 20)
        #plotter.add_cms_blurb(sqrts)
        #plotter.save('final-qcd1-subMass')
        #
        #plotter.plot_qcd_contamination('m2JetPt', 2, 10)
        #plotter.add_cms_blurb(sqrts)
        #plotter.save('final-qcd2-m2JetPt')
        #
        #plotter.plot_qcd_contamination('m1JetPt', 1, 20)
        #plotter.add_cms_blurb(sqrts)
        #plotter.save('final-qcd1-m1JetPt')

        #END OF if not options.dry_run:
    ###########################################################################
    ##  Making shape file     #################################################
    ###########################################################################
    plotter.set_subdir('')
    prefixes = [options.prefix+'$'] if options.prefix else ['']
    prefixes = [i+'$' for i in options.prefixes.split(',') if i] if options.prefixes else prefixes
    for prefix in prefixes:
        #plotter.plot_final(prefix+'m2_t_Mass', 20, xaxis='m_{#mu_{2}#tau} (GeV)', qcd_weight_fraction=0.5)
        #plotter.add_cms_blurb(sqrts)
        #plotter.save('final-%s-qweight05-subMass' % options.prefix)
        #
        #plotter.plot_final_f3(prefix+'m2_t_Mass', 20, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)', qcd_correction=False, qcd_weight_fraction=0.5, show_error=True)
        #plotter.add_cms_blurb(sqrts)
        #plotter.save('final-%s-f3-qweight05-subMass' % options.prefix)

        shape_prefix = prefix if len(prefixes) > 1 else ''
        shape_prefix = shape_prefix.replace(':','_').replace('$','_')
        
        shape_file = ROOT.TFile(
            os.path.join(plotter.outputdir, 'LTCut_mmt_shapes_%s.root' % ( plotter.period) ), 'RECREATE')
        shape_dir = shape_file.mkdir('mmtCatHigh')
        plotter.write_shapes(prefix+'m2_t_Mass#LT', 20, shape_dir, qcd_fraction=0.5, project=[80, 650], project_axis='X')
        shape_dir = shape_file.mkdir('mmtCatHigh_w')
        plotter.write_shapes(prefix+'m2_t_Mass#LT', 20, shape_dir, qcd_fraction=0.0, project=[80, 650], project_axis='X')
        shape_dir = shape_file.mkdir('mmtCatHigh_q')
        plotter.write_shapes(prefix+'m2_t_Mass#LT', 20, shape_dir, qcd_fraction=1.0, project=[80, 650], project_axis='X')
        shape_file.Close()

        #rebin_slim = range(20, 91, 10)+[110,200]
        shape_file = ROOT.TFile(
            os.path.join(plotter.outputdir, '%smmt_shapes_%s.root' % (shape_prefix, plotter.period) ), 'RECREATE')
        
        shape_dir = shape_file.mkdir('mmtCatLow')
        plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[0., 130], project_axis='X')
        shape_dir = shape_file.mkdir('mmtCatLow_w')
        plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[0., 130], project_axis='X')
        shape_dir = shape_file.mkdir('mmtCatLow_q')
        plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0., 130], project_axis='X')
        
        shape_dir = shape_file.mkdir('mmtCatHigh')
        plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[130, 650], project_axis='X')
        shape_dir = shape_file.mkdir('mmtCatHigh_w')
        plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[130, 650], project_axis='X')
        shape_dir = shape_file.mkdir('mmtCatHigh_q')
        plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[130, 650], project_axis='X')
        
        logging.warning('shape file %s created' % shape_file.GetName()) 
        shape_file.Close()



        shape_file = ROOT.TFile(
            os.path.join(plotter.outputdir, '%smmt_f3_shapes_%s.root' % (shape_prefix, plotter.period) ), 'RECREATE')
        
        shape_dir = shape_file.mkdir('mmtCatLow')
        plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[0., 130], project_axis='X')
        shape_dir = shape_file.mkdir('mmtCatLow_w')
        plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[0., 130], project_axis='X')
        shape_dir = shape_file.mkdir('mmtCatLow_q')
        plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0., 130], project_axis='X')
        
        shape_dir = shape_file.mkdir('mmtCatHigh')
        plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[130, 650], project_axis='X')
        shape_dir = shape_file.mkdir('mmtCatHigh_w')
        plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[130, 650], project_axis='X')
        shape_dir = shape_file.mkdir('mmtCatHigh_q')
        plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[130, 650], project_axis='X')
        
        logging.warning('shape file %s created' % shape_file.GetName()) 
        shape_file.Close()

        shape_file = ROOT.TFile(
            os.path.join(plotter.outputdir, '%smmt_f3_all_shapes_%s.root' % (shape_prefix, plotter.period) ), 'RECREATE')
        
        shape_dir = shape_file.mkdir('mmtCatLow')
        plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[0., 130], project_axis='X')
        shape_dir = shape_file.mkdir('mmtCatLow_w')
        plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[0., 130], project_axis='X')
        shape_dir = shape_file.mkdir('mmtCatLow_q')
        plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0., 130], project_axis='X')
        
        shape_dir = shape_file.mkdir('mmtCatHigh')
        plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[130, 650], project_axis='X')
        shape_dir = shape_file.mkdir('mmtCatHigh_w')
        plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[130, 650], project_axis='X')
        shape_dir = shape_file.mkdir('mmtCatHigh_q')
        plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[130, 650], project_axis='X')
        
        shape_dir = shape_file.mkdir('mmtCatLowf3')
        plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[0., 130], project_axis='X')
        shape_dir = shape_file.mkdir('mmtCatLowf3_w')
        plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[0., 130], project_axis='X')
        shape_dir = shape_file.mkdir('mmtCatLowf3_q')
        plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0., 130], project_axis='X')
        
        shape_dir = shape_file.mkdir('mmtCatHighf3')
        plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[130, 650], project_axis='X')
        shape_dir = shape_file.mkdir('mmtCatHighf3_w')
        plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[130, 650], project_axis='X')
        shape_dir = shape_file.mkdir('mmtCatHighf3_q')
        plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[130, 650], project_axis='X')
        
        logging.warning('shape file %s created' % shape_file.GetName()) 
        shape_file.Close()

        

