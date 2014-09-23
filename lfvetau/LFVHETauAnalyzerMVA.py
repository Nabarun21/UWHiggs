from ETauTree import ETauTree
import os
import ROOT
import math
import glob
import array
import mcCorrections
import baseSelections as selections
import FinalStateAnalysis.PlotTools.pytree as pytree
from FinalStateAnalysis.PlotTools.decorators import  memo_last
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
from math import sqrt, pi, cos
from fakerate_functions import fakerate_central_histogram, fakerate_p1s_histogram, fakerate_m1s_histogram

def collmass(row, met, metPhi):
    ptnu =abs(met*cos(deltaPhi(metPhi, row.tPhi)))
    visfrac = row.tPt/(row.tPt+ptnu)
    #print met, cos(deltaPhi(metPhi, row.tPhi)), ptnu, visfrac
    return (row.e_t_Mass / sqrt(visfrac))

def deltaPhi(phi1, phi2):
    PHI = abs(phi1-phi2)
    if PHI<=pi:
        return PHI
    else:
        return 2*pi-PHI
def deltaR(phi1, ph2, eta1, eta2):
    deta = eta1 - eta2
    dphi = abs(phi1-phi2)
    if (dphi>pi) : dphi = 2*pi-dphi
    return sqrt(deta*deta + dphi*dphi);

class LFVHETauAnalyzerMVA(MegaBase):
    tree = 'et/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel='ET'
        super(LFVHETauAnalyzerMVA, self).__init__(tree, outfile, **kwargs)
        self.tree = ETauTree(tree)
        self.out=outfile
        self.histograms = {}
        self.pucorrector = mcCorrections.make_puCorrector('singlee')
        self.pucorrectorUp = mcCorrections.make_puCorrectorUp('singlee')
        self.pucorrectorDown = mcCorrections.make_puCorrectorDown('singlee')
     

    @staticmethod 
    def tau_veto(row):
        if not row.tAntiMuonLoose2 or not row.tAntiElectronMVA5Tight or not row.tDecayFinding :
            return False

    @staticmethod
    def obj1_matches_gen(row):
        return row.eGenPdgId == -1*row.eCharge*11
    @staticmethod 
    def obj3_matches_gen(row):
        return t.genDecayMode != -2 

    
    def event_weight(self, row):
        if row.run > 2: #FIXME! add tight ID correction
            return [1.]
        allmcCorrections=    mcCorrections.get_electronId_corrections13_MVA(row, 'e') * \
                          mcCorrections.get_electronIso_corrections13_MVA(row, 'e') * mcCorrections.get_trigger_corrections_MVA(row,'e') 
                       

        trUp_mcCorrections =   mcCorrections.get_electronId_corrections13_MVA(row, 'e') * \
                          mcCorrections.get_electronIso_corrections13_MVA(row, 'e') *  mcCorrections.get_trigger_corrections_p1s_MVA(row,'e') 
        trDown_mcCorrections = mcCorrections.get_electronId_corrections13_MVA(row, 'e') * \
                               mcCorrections.get_electronIso_corrections13_MVA(row, 'e') *  mcCorrections.get_trigger_corrections_m1s_MVA(row,'e') 

        eidUp_mcCorrections=  mcCorrections.get_electronId_corrections13_p1s_MVA(row, 'e') *\
                              mcCorrections.get_electronIso_corrections13_MVA(row, 'e') *  mcCorrections.get_trigger_corrections_MVA(row,'e') 
        eidDown_mcCorrections= mcCorrections.get_electronId_corrections13_m1s_MVA(row, 'e') * \
                                mcCorrections.get_electronIso_corrections13_MVA(row, 'e') *  mcCorrections.get_trigger_corrections_MVA(row,'e') 
        eisoUp_mcCorrections=    mcCorrections.get_electronId_corrections13_MVA(row, 'e') * \
                               mcCorrections.get_electronIso_corrections13_p1s_MVA(row, 'e') * mcCorrections.get_trigger_corrections_MVA(row,'e') 
        eisoDown_mcCorrections= mcCorrections.get_electronId_corrections13_m1s_MVA(row, 'e') * \
                                 mcCorrections.get_electronIso_corrections13_p1s_MVA(row, 'e') * mcCorrections.get_trigger_corrections_MVA(row,'e') 
    
        #pucorrlist = self.pucorrector(row.nTruePU)
        
        weight =  self.pucorrector(row.nTruePU) *\
                 allmcCorrections
        weight_up =  self.pucorrectorUp(row.nTruePU) *\
                    allmcCorrections
        weight_down =  self.pucorrectorDown(row.nTruePU) *\
                      allmcCorrections
        
        weight_tr_up = self.pucorrector(row.nTruePU) *\
                       trUp_mcCorrections
        weight_tr_down = self.pucorrector(row.nTruePU) *\
                         trDown_mcCorrections

        
        weight_eid_up =  self.pucorrector(row.nTruePU) *\
                 eidUp_mcCorrections
        weight_eid_down =  self.pucorrector(row.nTruePU) *\
                 eidDown_mcCorrections
        weight_eiso_up =  self.pucorrector(row.nTruePU) *\
                 eisoUp_mcCorrections
        weight_eiso_down =  self.pucorrector(row.nTruePU) *\
                 eisoDown_mcCorrections
        
        return  [weight, weight_up, weight_down, weight_tr_up,  weight_tr_down, weight_eid_up, weight_eid_down, weight_eiso_up,  weight_eiso_down,]

## 
    def begin(self):

        processtype=['gg']
        threshold=['ept30']
        sign=['os', 'ss']
        jetN = [0, 1, 2, 3]
        folder=[]
        pudir = ['','p1s/', 'm1s/','trp1s/', 'trm1s/', 'eidp1s/','eidm1s/',  'eisop1s/','eisom1s/', 'tLoose/','tLooseUp/','tLooseDown/', ]

        for d  in pudir :
            for i in sign:
                for j in processtype:
                    for k in threshold:
                        for jn in jetN: 

                            folder.append(d+i+'/'+j+'/'+k +'/'+str(jn))
                            folder.append(d+i+'/'+j+'/'+k +'/'+str(jn)+'/selected')
                        
        for f in folder: 
            self.book(f,"tPt", "tau p_{T}", 200, 0, 200)
            self.book(f,"tPhi", "tau phi", 100, -3.2, 3.2)
            self.book(f,"tEta", "tau eta",  50, -2.5, 2.5)
            
            self.book(f,"ePt", "e p_{T}", 200, 0, 200)
            self.book(f,"ePhi", "e phi",  100, -3.2, 3.2)
            self.book(f,"eEta", "e eta", 50, -2.5, 2.5)
            
            self.book(f, "et_DeltaPhi", "e-tau DeltaPhi" , 50, 0, 3.2)
            self.book(f, "et_DeltaR", "e-tau DeltaR" , 50, 0, 3.2)
            
            self.book(f, "h_collmass_pfmet",  "h_collmass_pfmet",  32, 0, 320)
            self.book(f, "h_collmass_mvamet",  "h_collmass_mvamet",  32, 0, 320)
            self.book(f, "h_collmass_pfmet_Ty1",  "h_collmass_pfmet_Ty1",  32, 0, 320)
            self.book(f, "h_collmass_pfmet_jes", "h_collmass_pfmet_jes", 50, 0, 100)
            self.book(f, "h_collmass_pfmet_mes", "h_collmass_pfmet_mes", 50, 0, 100 )
            self.book(f, "h_collmass_pfmet_tes", "h_collmass_pfmet_tes", 50, 0, 100)
            self.book(f, "h_collmass_pfmet_ees", "h_collmass_pfmet_ees", 50, 0, 100)
            self.book(f, "h_collmass_pfmet_ues", "h_collmass_pfmet_ues", 50, 0, 100)



            self.book(f, "h_collmassSpread_pfmet",  "h_collmassSpread_pfmet",  40, -100, 100)
            self.book(f, "h_collmassSpread_mvamet",  "h_collmassSpread_mvamet",  40, -100, 100)
            self.book(f, "h_collmassSpread_lowPhi_pfmet",  "h_collmassSpread_lowPhi_pfmet",  40, -100, 100)
            self.book(f, "h_collmassSpread_lowPhi_mvamet",  "h_collmassSpread_lowPhi_mvamet", 40, -100, 100)
            self.book(f, "h_collmassSpread_highPhi_pfmet",  "h_collmassSpread_highPhi_pfmet", 40, -100, 100)
            self.book(f, "h_collmassSpread_highPhi_mvamet",  "h_collmassSpread_highPhi_mvamet", 40, -100, 100)
            self.book(f, "h_collmass_lowPhi_pfmet",  "h_collmass_lowPhi_pfmet",  32, 0, 320)
            self.book(f, "h_collmass_lowPhi_mvamet",  "h_collmass_lowPhi_mvamet",  32, 0, 320)
            self.book(f, "h_collmass_highPhi_pfmet",  "h_collmass_highPhi_pfmet",  32, 0, 320)
            self.book(f, "h_collmass_highPhi_mvamet", "h_collmass_highPhi_mvamet",  32, 0, 320)
            self.book(f, "h_collmass_vs_dPhi_pfmet",  "h_collmass_vs_dPhi_pfmet", 50, 0, 3.2, 32, 0, 320, type=ROOT.TH2F)
            self.book(f, "h_collmass_vs_dPhi_mvamet",  "h_collmass_vs_dPhi_mvamet", 50, 0, 3.2, 32, 0, 320, type=ROOT.TH2F)
            self.book(f, "h_collmassSpread_vs_dPhi_pfmet",  "h_collmassSpread_vs_dPhi_pfmet", 50, 0, 3.2, 20, -100, 100, type=ROOT.TH2F)
            self.book(f, "h_collmassSpread_vs_dPhi_mvamet",  "h_collmassSpread_vs_dPhi_mvamet", 50, 0, 3.2, 20, -100, 100, type=ROOT.TH2F)
            
            self.book(f, "h_vismass",  "h_vismass",  32, 0, 320)
            
            self.book(f, "type1_pfMetEt_vs_dPhi", "PFMet vs #Delta#phi(#tau,PFMet)", 50, 0, 3.2, 64, 0, 320, type=ROOT.TH2F)
            self.book(f, "mvaMetEt_vs_dPhi", "MVAMet vs #Delta#phi(#tau,MVAMet)", 50, 0, 3.2, 64, 0, 320, type=ROOT.TH2F)

            self.book(f, "tPFMET_DeltaPhi", "tau-PFMET DeltaPhi" , 50, 0, 3.2)
            self.book(f, "tPFMET_Mt", "tau-PFMET M_{T}" , 200, 0, 200)
            self.book(f, "tPFMET_DeltaPhi_Ty1", "tau-type1PFMET DeltaPhi" , 50, 0, 3.2)
            self.book(f, "tPFMET_Mt_Ty1", "tau-type1PFMET M_{T}" , 200, 0, 200)
            self.book(f, 'tPFMET_Mt_jes', "tau-MVAMET M_{T} JES" , 200, 0, 200)
            self.book(f, 'tPFMET_Mt_mes', "tau-MVAMET M_{T} JES" , 200, 0, 200)
            self.book(f, 'tPFMET_Mt_ees', "tau-MVAMET M_{T} JES" , 200, 0, 200)
            self.book(f, 'tPFMET_Mt_tes', "tau-MVAMET M_{T} JES" , 200, 0, 200)
            self.book(f, 'tPFMET_Mt_ues', "tau-MVAMET M_{T} JES" , 200, 0, 200)
            
            self.book(f, "tMVAMET_DeltaPhi", "tau-MVAMET DeltaPhi" , 50, 0, 3.2)
            self.book(f, "tMVAMET_Mt", "tau-MVAMET M_{T}" , 200, 0, 200)
               
            self.book(f, "ePFMET_DeltaPhi_Ty1", "e-type1PFMET DeltaPhi" , 50, 0, 3.2)
            self.book(f, "ePFMET_Mt_Ty1", "e-type1PFMET M_{T}" , 200, 0, 200)
            self.book(f, "ePFMET_DeltaPhi", "e-PFMET DeltaPhi" , 50, 0, 3.2)
            self.book(f, "ePFMET_Mt", "e-PFMET M_{T}" , 200, 0, 200)
            self.book(f, 'ePFMET_Mt_jes', "e-MVAMET M_{T} JES" , 200, 0, 200)
            self.book(f, 'ePFMET_Mt_mes', "e-MVAMET M_{T} JES" , 200, 0, 200)
            self.book(f, 'ePFMET_Mt_ees', "e-MVAMET M_{T} JES" , 200, 0, 200)
            self.book(f, 'ePFMET_Mt_tes', "e-MVAMET M_{T} JES" , 200, 0, 200)
            self.book(f, 'ePFMET_Mt_ues', "e-MVAMET M_{T} JES" , 200, 0, 200)

            self.book(f, "eMVAMET_DeltaPhi", "e-MVAMET DeltaPhi" , 50, 0, 3.2)
            self.book(f, "eMVAMET_Mt", "e-MVAMET M_{T}" , 200, 0, 200)
            
            self.book(f, "jetN_20", "Number of jets, p_{T}>20", 10, -0.5, 9.5) 
            self.book(f, "jetN_30", "Number of jets, p_{T}>30", 10, -0.5, 9.5) 

    def fakerate_weights(self, tEta, central_weights, p1s_weights, m1s_weights):
        frweight=[1.,1.,1.]

        #central_weights = fakerate_central_histogram(25,0, 2.5)
        #p1s_weights = fakerate_central_histogram(25,0, 2.5)
        #m1s_weights = fakerate_central_histogram(25,0, 2.5)

        for n,w in enumerate( central_weights ):
            if abs(tEta) < w[1]:
                break
            frweight[0] = w[0]
            frweight[1] = p1s_weights[n][0]
            frweight[2] = m1s_weights[n][0]
 
        
        return  frweight;

    
                    
    def fill_histos(self, row, f='os/gg/ept0/0',  isTauTight=False, frw=[1.,1.,1.]):
        weight = self.event_weight(row)
        histos = self.histograms
        pudir =['']
        if row.run < 2: pudir.extend( ['p1s/', 'm1s/', 'trp1s/', 'trm1s/', 'eidp1s/','eidm1s/',  'eisop1s/','eisom1s/'])
        looseList = ['tLoose/', 'tLooseUp/', 'tLooseDown/']

        #print isTauTight

        if not isTauTight:
            frweight_bv = frw[0]/(1.-frw[0])
            #print fakerate, frw[0], frweight
            err = abs(frw[0] - frw[1])* abs(frw[0]/((1-frw[0])*(1-frw[0])) + 1/(1-frw[0]))
            frweight_p1s = frweight_bv*(1+err)
            frweight_m1s = frweight_bv*(1-err)
        
            fr_weights = [frweight_bv, frweight_p1s, frweight_m1s]
        
            for n, l in enumerate(looseList) :
                frweight = weight[0]*fr_weights[n]
                folder = l+f
                #print folder , frweight, fr_weights[n], frw[0]
                histos[folder+'/tPt'].Fill(row.tPt, frweight)
                histos[folder+'/tEta'].Fill(row.tEta, frweight)
                histos[folder+'/tPhi'].Fill(row.tPhi, frweight) 
                histos[folder+'/ePt'].Fill(row.ePt, frweight)
                histos[folder+'/eEta'].Fill(row.eEta, frweight)
                histos[folder+'/ePhi'].Fill(row.ePhi, frweight)
                histos[folder+'/et_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.tPhi), frweight)
                histos[folder+'/et_DeltaR'].Fill(row.e_t_DR, frweight)
                histos[folder+'/h_collmass_vs_dPhi_pfmet'].Fill(deltaPhi(row.tPhi, row.type1_pfMetPhi), collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi), frweight)
                histos[folder+'/h_collmass_vs_dPhi_mvamet'].Fill(deltaPhi(row.tPhi, row.mva_metPhi), collmass(row, row.mva_metEt, row.mva_metPhi), frweight)
                histos[folder+'/h_collmassSpread_vs_dPhi_pfmet'].Fill(deltaPhi(row.tPhi, row.type1_pfMetPhi), collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)-125.0, frweight)
                histos[folder+'/h_collmassSpread_vs_dPhi_mvamet'].Fill(deltaPhi(row.tPhi, row.mva_metPhi), collmass(row, row.mva_metEt, row.mva_metPhi)-125.0, frweight)
                if deltaPhi(row.tPhi, row.pfMetPhi) > 1.57 :  
                    histos[folder+'/h_collmass_highPhi_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi), frweight)
                    histos[folder+'/h_collmassSpread_highPhi_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)-125.0, frweight)
                if deltaPhi(row.tPhi, row.pfMetPhi) < 1.57 :  
                    histos[folder+'/h_collmass_lowPhi_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi), frweight)
                    histos[folder+'/h_collmassSpread_lowPhi_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)-125.0, frweight)
                if deltaPhi(row.tPhi, row.mva_metPhi) > 1.57 :  
                    histos[folder+'/h_collmass_highPhi_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi), frweight)
                    histos[folder+'/h_collmassSpread_highPhi_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi)-125.0, frweight)
                if deltaPhi(row.tPhi, row.mva_metPhi) < 1.57 :  
                    histos[folder+'/h_collmass_lowPhi_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi), frweight)
                    histos[folder+'/h_collmassSpread_lowPhi_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi)-125.0, frweight)
                histos[folder+'/h_collmassSpread_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)-125.0, frweight)
                histos[folder+'/h_collmassSpread_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi)-125.0, frweight)
                histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.pfMetEt, row.pfMetPhi), frweight)
                histos[folder+'/h_collmass_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi), frweight)
                histos[folder+'/h_collmass_pfmet_Ty1'].Fill(collmass(row, row.type1_pfMetEt, row.pfMetPhi), frweight)
                histos[folder+'/h_collmass_pfmet_jes'].Fill(collmass(row, row.pfMet_jes_Et, row.pfMet_jes_Phi), frweight)
                histos[folder+'/h_collmass_pfmet_mes'].Fill(collmass(row, row.pfMet_mes_Et, row.pfMet_mes_Phi), frweight)
                histos[folder+'/h_collmass_pfmet_tes'].Fill(collmass(row, row.pfMet_tes_Et, row.pfMet_tes_Phi), frweight)
                histos[folder+'/h_collmass_pfmet_ees'].Fill(collmass(row, row.pfMet_ees_Et, row.pfMet_ees_Phi), frweight)
                histos[folder+'/h_collmass_pfmet_ues'].Fill(collmass(row, row.pfMet_ues_Et, row.pfMet_ues_Phi), frweight)


                histos[folder+'/h_vismass'].Fill(row.e_t_Mass, frweight)
                histos[folder+'/type1_pfMetEt_vs_dPhi'].Fill(deltaPhi(row.tPhi, row.type1_pfMetPhi), row.type1_pfMetEt, frweight)
                histos[folder+'/mvaMetEt_vs_dPhi'].Fill(deltaPhi(row.tPhi, row.mva_metPhi), row.mva_metEt, frweight)

                histos[folder+'/ePFMET_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.pfMetPhi), frweight)
                histos[folder+'/ePFMET_DeltaPhi_Ty1'].Fill(deltaPhi(row.ePhi, row.type1_pfMetPhi), frweight)
                histos[folder+'/eMVAMET_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.mva_metPhi), frweight)
                histos[folder+'/ePFMET_Mt'].Fill(row.eMtToPFMET, frweight)
                histos[folder+'/ePFMET_Mt_Ty1'].Fill(row.eMtToPfMet_Ty1, frweight)
                histos[folder+'/ePFMET_Mt_jes'].Fill(row.eMtToPfMet_jes, frweight)
                histos[folder+'/ePFMET_Mt_mes'].Fill(row.eMtToPfMet_mes, frweight)
                histos[folder+'/ePFMET_Mt_ees'].Fill(row.eMtToPfMet_ees, frweight)
                histos[folder+'/ePFMET_Mt_tes'].Fill(row.eMtToPfMet_tes, frweight)
                histos[folder+'/ePFMET_Mt_ues'].Fill(row.eMtToPfMet_ues, frweight)
                histos[folder+'/eMVAMET_Mt'].Fill(row.eMtToMVAMET, frweight)

                histos[folder+'/tPFMET_DeltaPhi'].Fill(deltaPhi(row.tPhi, row.pfMetPhi), frweight)
                histos[folder+'/tPFMET_DeltaPhi_Ty1'].Fill(deltaPhi(row.tPhi, row.type1_pfMetPhi), frweight)
                histos[folder+'/tMVAMET_DeltaPhi'].Fill(deltaPhi(row.tPhi, row.mva_metPhi), frweight)
                histos[folder+'/tPFMET_Mt'].Fill(row.tMtToPFMET, frweight)
                histos[folder+'/tMVAMET_Mt'].Fill(row.tMtToMVAMET, frweight)
                histos[folder+'/tPFMET_Mt_jes'].Fill(row.tMtToPfMet_jes, frweight)
                histos[folder+'/tPFMET_Mt_mes'].Fill(row.tMtToPfMet_mes, frweight)
                histos[folder+'/tPFMET_Mt_ees'].Fill(row.tMtToPfMet_ees, frweight)
                histos[folder+'/tPFMET_Mt_tes'].Fill(row.tMtToPfMet_tes, frweight)
                histos[folder+'/tPFMET_Mt_ues'].Fill(row.tMtToPfMet_ues, frweight)

                histos[folder+'/jetN_20'].Fill(row.jetVeto20, frweight) 
                histos[folder+'/jetN_30'].Fill(row.jetVeto30, frweight) 

        else:
            for n,d  in enumerate(pudir) :
        
                folder = d+f
                histos[folder+'/tPt'].Fill(row.tPt, weight[n])
                histos[folder+'/tEta'].Fill(row.tEta, weight[n])
                histos[folder+'/tPhi'].Fill(row.tPhi, weight[n]) 
                histos[folder+'/ePt'].Fill(row.ePt, weight[n])
                histos[folder+'/eEta'].Fill(row.eEta, weight[n])
                histos[folder+'/ePhi'].Fill(row.ePhi, weight[n])
                histos[folder+'/et_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.tPhi), weight[n])
                histos[folder+'/et_DeltaR'].Fill(row.e_t_DR, weight[n])
                histos[folder+'/h_collmass_vs_dPhi_pfmet'].Fill(deltaPhi(row.tPhi, row.type1_pfMetPhi), collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi), weight[n])
                histos[folder+'/h_collmass_vs_dPhi_mvamet'].Fill(deltaPhi(row.tPhi, row.mva_metPhi), collmass(row, row.mva_metEt, row.mva_metPhi), weight[n])
                histos[folder+'/h_collmassSpread_vs_dPhi_pfmet'].Fill(deltaPhi(row.tPhi, row.type1_pfMetPhi), collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)-125.0, weight[n])
                histos[folder+'/h_collmassSpread_vs_dPhi_mvamet'].Fill(deltaPhi(row.tPhi, row.mva_metPhi), collmass(row, row.mva_metEt, row.mva_metPhi)-125.0, weight[n])
                if deltaPhi(row.tPhi, row.pfMetPhi) > 1.57 :  
                    histos[folder+'/h_collmass_highPhi_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi), weight[n])
                    histos[folder+'/h_collmassSpread_highPhi_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)-125.0, weight[n])
                if deltaPhi(row.tPhi, row.pfMetPhi) < 1.57 :  
                    histos[folder+'/h_collmass_lowPhi_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi), weight[n])
                    histos[folder+'/h_collmassSpread_lowPhi_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)-125.0, weight[n])
                if deltaPhi(row.tPhi, row.mva_metPhi) > 1.57 :  
                    histos[folder+'/h_collmass_highPhi_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi), weight[n])
                    histos[folder+'/h_collmassSpread_highPhi_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi)-125.0, weight[n])
                if deltaPhi(row.tPhi, row.mva_metPhi) < 1.57 :  
                    histos[folder+'/h_collmass_lowPhi_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi), weight[n])
                    histos[folder+'/h_collmassSpread_lowPhi_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi)-125.0, weight[n])
                histos[folder+'/h_collmassSpread_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)-125.0, weight[n])
                histos[folder+'/h_collmassSpread_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi)-125.0, weight[n])
                histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi), weight[n])
                histos[folder+'/h_collmass_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi), weight[n])
                histos[folder+'/h_vismass'].Fill(row.e_t_Mass, weight[n])
                histos[folder+'/type1_pfMetEt_vs_dPhi'].Fill(deltaPhi(row.tPhi, row.type1_pfMetPhi), row.type1_pfMetEt, weight[n])
                histos[folder+'/mvaMetEt_vs_dPhi'].Fill(deltaPhi(row.tPhi, row.mva_metPhi), row.mva_metEt, weight[n])
                histos[folder+'/ePFMET_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.type1_pfMetPhi), weight[n])
                histos[folder+'/eMVAMET_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.mva_metPhi), weight[n])
                histos[folder+'/ePFMET_Mt'].Fill(row.eMtToPFMET, weight[n])
                histos[folder+'/eMVAMET_Mt'].Fill(row.eMtToMVAMET, weight[n])
                histos[folder+'/tPFMET_DeltaPhi'].Fill(deltaPhi(row.tPhi, row.type1_pfMetPhi), weight[n])
                histos[folder+'/tMVAMET_DeltaPhi'].Fill(deltaPhi(row.tPhi, row.mva_metPhi), weight[n])
                histos[folder+'/tPFMET_Mt'].Fill(row.tMtToPFMET, weight[n])
                histos[folder+'/tMVAMET_Mt'].Fill(row.tMtToMVAMET, weight[n])
                histos[folder+'/jetN_20'].Fill(row.jetVeto20, weight[n]) 
                histos[folder+'/jetN_30'].Fill(row.jetVeto30, weight[n]) 



        

    def process(self):
        
        central_weights = fakerate_central_histogram(25,0, 2.5)
        p1s_weights = fakerate_p1s_histogram(25,0, 2.5)#fakerate_p1s_histogram(25,0, 2.5)
        m1s_weights = fakerate_m1s_histogram(25,0, 2.5)#fakerate_m1s_histogram(25,0, 2.5)

        frw = []
        for row in self.tree:
 
            sign = 'ss' if row.e_t_SS else 'os'
            processtype = '' ## use a line as for sign when the vbf when selections are defined            
            ptthreshold = [30]
            processtype ='gg'##changed from 20
            jn = row.jetVeto30
            if jn > 3 : jn = 3
            #if row.run > 2 : #apply the trigger to data only (MC triggers enter in the scale factors)
            if not bool(row.singleE27WP80Pass) : continue
            if  not  bool(row.eMatchesSingleE27WP80): continue
            
            if not selections.eSelection(row, 'e'): continue
            if not selections.lepton_id_iso(row, 'e', 'eid13Tight_etauiso01'): continue
            if row.eEta > 1.4442 and row.eEta < 1.566 : continue
            if not selections.tauSelection(row, 't'): continue
 
            if not row.tAntiElectronMVA5Tight : continue
            if not row.tAntiMuon2Loose : continue
            
            if not row.tLooseIso3Hits : continue
            #isTauTight = False
            frw=self.fakerate_weights(row.tEta, central_weights, p1s_weights, m1s_weights )
            
            isTauTight = bool(row.tTightIso3Hits)
            
            #print bool(row.tTightIso3Hits)

            if row.tauVetoPt20EleTight3MuLoose : continue 
            #if row.tauHpsVetoPt20 : continue
            if row.muVetoPt5IsoIdVtx : continue
            if row.eVetoCicLooseIso : continue # change it with Loose
  
            for j in ptthreshold:
                folder = sign+'/'+processtype+'/ept'+str(j)+'/'+str(int(jn))
                
                if row.ePt < j : continue
                  
                self.fill_histos(row, folder,isTauTight, frw)
                print 'histo filled'
                #selections
                if jn == 0 :
                    if row.tPt < 35: continue 
                    if row.ePt < 40 : continue
                    if deltaPhi(row.ePhi, row.tPhi) < 2.7 : continue
                    if row.tMtToPFMET > 50 : continue
                if jn == 1 :
                    if row.tPt < 40: continue 
                    if row.ePt < 35 : continue
                    if row.tMtToPFMET > 35 : continue
                if jn == 2 :
                    if row.tPt < 40: continue 
                    if row.ePt < 30 : continue # no cut as only electrons with pt>30 are in the ntuples
                    if row.tMtToPFMET > 35 : continue
                    if row.vbfMass < 550 : continue
                    if row.vbfDeta < 3.5 : continue
                folder = sign+'/'+processtype+'/ept'+str(j)+'/'+str(int(jn))+'/selected'
                self.fill_histos(row, folder, isTauTight,frw)
                
                 
                    
             
            
    def finish(self):
        self.write_histos()


