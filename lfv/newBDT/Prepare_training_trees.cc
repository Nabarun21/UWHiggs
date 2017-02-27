#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include "TMultiGraph.h"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <utility>
#include <stdio.h>
#include <TF1.h>
#include <TDirectoryFile.h>
#include <TRandom3.h>
#include "TLorentzVector.h"
#include "TString.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TKey.h"
#include "THashList.h"
#include "THStack.h"
#include "TPaveLabel.h"
#include "TFile.h"
#include "myHelper.h"
#include "tr_Tree.h"
#include "ScaleFactor.h"
#include "LumiReweightingStandAlone.h"

using namespace std;

int main(int argc, char** argv) {

    std::string input = *(argv + 1);
    std::string output = *(argv + 2);
    std::string sample = *(argv + 3);
    float tes=0;
    if (argc > 1) {
        tes = atof(argv[4]);
    }

    TFile *f_Double = new TFile(input.c_str());
    cout<<"XXXXXXXXXXXXX "<<input.c_str()<<" XXXXXXXXXXXX"<<endl;
    TTree *arbre = (TTree*) f_Double->Get("RLE_tree");
    TH1F* nbevt = (TH1F*) f_Double->Get("nevents");
    float ngen = nbevt->GetBinContent(1);

    float xs=1.0; float weight=1.0; float luminosity=20100.0;
    if (sample=="Ztt" or sample=="Zll"){ xs=6025.0; weight=luminosity*xs/ngen;}
    else if (sample=="TT") {xs=809.0; weight=luminosity*xs/ngen;}
    else if (sample=="W") {xs=61526.7; weight=luminosity*xs/ngen;}
    else if (sample=="QCD") {xs=720648000*0.00042; weight=luminosity*xs/ngen;}
    else if (sample=="data_obs"){weight=1.0;}
    else if (sample=="WZ1L1Nu2Q") {xs=10.71; weight=luminosity*xs/ngen;}
    else if (sample=="WZ1L3Nu") {xs=3.05; weight=luminosity*xs/ngen;}
    else if (sample=="WZJets") {xs=5.26; weight=luminosity*xs/ngen;}
    else if (sample=="WZ2L2Q") {xs=5.595; weight=luminosity*xs/ngen;}
    else if (sample=="WW1L1Nu2Q") {xs=49.997; weight=luminosity*xs/ngen;}
    else if (sample=="ZZ4L") {xs=1.212; weight=luminosity*xs/ngen;}
    else if (sample=="ZZ2L2Q") {xs=3.22; weight=luminosity*xs/ngen;}
    else if (sample=="VV2L2Nu") {xs=11.95; weight=luminosity*xs/ngen;}
    else if (sample=="ST_tW_antitop") {xs=35.6; weight=luminosity*xs/ngen;}
    else if (sample=="ST_tW_top") {xs=35.6; weight=luminosity*xs/ngen;}
    else if (sample=="ST_t_antitop") {xs=80.95; weight=luminosity*xs/ngen;}
    else if (sample=="ST_t_top") {xs=136.02; weight=luminosity*xs/ngen;}
    else if (sample=="Wg") {xs=498.0; weight=luminosity*xs/ngen;}
    else if (sample=="ZZ") {xs=16.523; weight=luminosity*xs/ngen;}
    else if (sample=="WZ") {xs=47.13; weight=luminosity*xs/ngen;}
    else if (sample=="WW") {xs=118.7; weight=luminosity*xs/ngen;}
    else if (sample=="ggH") {xs=43.92; weight=luminosity*xs/ngen;}
    else if (sample=="VBF") {xs=3.748; weight=luminosity*xs/ngen;}

    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(10);
    arbre->SetBranchAddress("npu", &npu);
    arbre->SetBranchAddress("run", &run);
    arbre->SetBranchAddress("lumi", &lumi);
    arbre->SetBranchAddress("evt", &evt);
    arbre->SetBranchAddress("npv", &npv);
    arbre->SetBranchAddress("pt_1", &pt_1);
    arbre->SetBranchAddress("phi_1", &phi_1);
    arbre->SetBranchAddress("eta_1", &eta_1);
    arbre->SetBranchAddress("iso_1", &iso_1);
    arbre->SetBranchAddress("m_1", &m_1);
    arbre->SetBranchAddress("q_1", &q_1);
    arbre->SetBranchAddress("nbtag", &nbtag);
    arbre->SetBranchAddress("q_2", &q_2);
    arbre->SetBranchAddress("pt_2", &pt_2);
    arbre->SetBranchAddress("eta_2", &eta_2);
    arbre->SetBranchAddress("phi_2", &phi_2);
    arbre->SetBranchAddress("m_2", &m_2);
    arbre->SetBranchAddress("njets", &njets);
    arbre->SetBranchAddress("mt_2", &mt_2);
    arbre->SetBranchAddress("mt_1", &mt_1);
    arbre->SetBranchAddress("met", &met);
    arbre->SetBranchAddress("metphi", &metphi);
    arbre->SetBranchAddress("dphi_12", &dphi_12);
    arbre->SetBranchAddress("dphi_emet", &dphi_emet);
    arbre->SetBranchAddress("dphi_taumet", &dphi_taumet);
    arbre->SetBranchAddress("byTightIsolationMVArun2v1DBoldDMwLT_2",&byTightIsolationMVArun2v1DBoldDMwLT_2);
    arbre->SetBranchAddress("byLooseCombinedIsolationDeltaBetaCorr3Hits_2",&byLooseCombinedIsolationDeltaBetaCorr3Hits_2);
    arbre->SetBranchAddress("byMediumCombinedIsolationDeltaBetaCorr3Hits_2",&byMediumCombinedIsolationDeltaBetaCorr3Hits_2);
    arbre->SetBranchAddress("byTightCombinedIsolationDeltaBetaCorr3Hits_2",&byTightCombinedIsolationDeltaBetaCorr3Hits_2);
    arbre->SetBranchAddress("byVLooseIsolationMVArun2v1DBoldDMwLT_2",&byVLooseIsolationMVArun2v1DBoldDMwLT_2);
    arbre->SetBranchAddress("byLooseIsolationMVArun2v1DBoldDMwLT_2",&byLooseIsolationMVArun2v1DBoldDMwLT_2);
    arbre->SetBranchAddress("byMediumIsolationMVArun2v1DBoldDMwLT_2",&byMediumIsolationMVArun2v1DBoldDMwLT_2);
    arbre->SetBranchAddress("byTightIsolationMVArun2v1DBoldDMwLT_2",&byTightIsolationMVArun2v1DBoldDMwLT_2);
    arbre->SetBranchAddress("byVTightIsolationMVArun2v1DBoldDMwLT_2",&byVTightIsolationMVArun2v1DBoldDMwLT_2);
    arbre->SetBranchAddress("byIsolationMVA3oldDMwLTraw_2",&byIsolationMVA3oldDMwLTraw_2);
    arbre->SetBranchAddress("tau_id",&tau_id);
    arbre->SetBranchAddress("decayModeFinding_2",&decayModeFinding_2);
    arbre->SetBranchAddress("decayMode_2",&decayMode_2);
    arbre->SetBranchAddress("againstElectronVLooseMVA6_2",&againstElectronVLooseMVA6_2);
    arbre->SetBranchAddress("againstElectronLooseMVA6_2",&againstElectronLooseMVA6_2);
    arbre->SetBranchAddress("againstElectronMediumMVA6_2",&againstElectronMediumMVA6_2);
    arbre->SetBranchAddress("againstElectronTightMVA6_2",&againstElectronTightMVA6_2);
    arbre->SetBranchAddress("againstElectronVTightMVA6_2",&againstElectronVTightMVA6_2);
    arbre->SetBranchAddress("againstMuonTight3_2",&againstMuonTight3_2);
    arbre->SetBranchAddress("againstMuonLoose3_2",&againstMuonLoose3_2);
    arbre->SetBranchAddress("extratau_veto",&extratau_veto);
    arbre->SetBranchAddress("extramuon_veto",&extramuon_veto);
    arbre->SetBranchAddress("extraelec_veto",&extraelec_veto);
    arbre->SetBranchAddress("isZmt",&isZmt);
    arbre->SetBranchAddress("isZtt",&isZtt);
    arbre->SetBranchAddress("isZmm",&isZmm);
    arbre->SetBranchAddress("isZee",&isZee);
    arbre->SetBranchAddress("tightID_1",&tightID_1);
    arbre->SetBranchAddress("aMCatNLO_weight",&aMCatNLO_weight);
    arbre->SetBranchAddress("pu_weight",&pu_weight);
    arbre->SetBranchAddress("trg_weight",&trg_weight);
    arbre->SetBranchAddress("m_vis",&m_vis);
    arbre->SetBranchAddress("m_coll",&m_coll);
    arbre->SetBranchAddress("vbfmass",&vbfmass);
    arbre->SetBranchAddress("vbfeta",&vbfeta);
    arbre->SetBranchAddress("id90_1",&id90_1);
    arbre->SetBranchAddress("passEle25",&passEle25);
    arbre->SetBranchAddress("passEle27",&passEle27);
    arbre->SetBranchAddress("photonIso_2", &photonIso_2);
    arbre->SetBranchAddress("neutralIso_2", &neutralIso_2);
    arbre->SetBranchAddress("chargedIso_2", &chargedIso_2);
    arbre->SetBranchAddress("puIso_2", &puIso_2);
    arbre->SetBranchAddress("NUP", &NUP);
    arbre->SetBranchAddress("mjj", &mjj);
    arbre->SetBranchAddress("jdeta", &jdeta);

    TTree * arbre2 = new TTree("TreeS", "TreeS");
    arbre2->SetDirectory(0);
    arbre2->Branch("pt_1_", &pt_1_, "pt_1_/F");
    arbre2->Branch("pt_2_", &pt_2_, "pt_2_/F");
    arbre2->Branch("eta_1_", &eta_1_, "eta_1_/F");
    arbre2->Branch("eta_2_", &eta_2_, "eta_2_/F");
    arbre2->Branch("dphi_12_", &dphi_12_, "dphi_12_/F");
    arbre2->Branch("mt_1_", &mt_1_, "mt_1_/F");
    arbre2->Branch("mt_2_", &mt_2_, "mt_2_/F");
    arbre2->Branch("dphi_taumet_", &dphi_taumet_, "dphi_taumet_/F");
    arbre2->Branch("dphi_emet_", &dphi_emet_, "dphi_emet_/F");
    arbre2->Branch("met_", &met_, "met_/F");
    arbre2->Branch("njets_", &njets_, "njets_/F");
    arbre2->Branch("weight_", &weight_, "weight_/F");
    arbre2->Branch("vbfeta_", &vbfeta_, "vbfeta_/F");
    arbre2->Branch("vbfmass_", &vbfmass_, "vbfmass_/F");
    arbre2->Branch("m_coll_", &m_coll_, "m_coll_/F");

   reweight::LumiReWeighting* LumiWeights_12;
   LumiWeights_12 = new reweight::LumiReWeighting("MC_Spring16_PU25ns_V1.root", "MyDataPileupHistogram.root", "pileup", "pileup");

   ScaleFactor * myScaleFactor_trg = new ScaleFactor();
   myScaleFactor_trg->init_ScaleFactor("LeptonEfficiencies/Electron/Run2016BCD/Electron_Ele25eta2p1WPTight_eff.root");
   ScaleFactor * myScaleFactor_id = new ScaleFactor();
   myScaleFactor_id->init_ScaleFactor("LeptonEfficiencies/Electron/Run2016BCD/Electron_IdIso0p10_eff.root");

   float frw=0.5;
   Int_t nentries_wtn = (Int_t) arbre->GetEntries();
   for (Int_t i = 0; i < nentries_wtn; i++) {
        arbre->GetEntry(i);
        if (i % 20000 == 0) fprintf(stdout, "\r  Processed events: %8d of %8d ", i, nentries_wtn);
        fflush(stdout);
	if (extratau_veto>0) continue;
	if (sample=="Data" && !passEle25) continue;
   	if (!againstElectronVTightMVA6_2 or !againstMuonLoose3_2) continue;
	if (sample=="Ztt" && !isZtt) continue;
        if (sample=="Zll" && isZtt) continue;

      	frw=0.12392-0.000505*(pt_2-30);
	
	float sf_trg=1.0;
	float sf_id=1.0;

        if (sample=="W"){
            if (NUP==5)
                weight=13.516;
            if (NUP==6)
                weight=4.126;
            if (NUP==7)
                weight=2.162;
            if (NUP==8)
                weight=1.665;
            if (NUP==9)
                weight=1.814;
        }


        if (sample!="Data"){
           sf_trg=myScaleFactor_trg->get_EfficiencyData(pt_1,eta_1);
           sf_id=myScaleFactor_id->get_ScaleFactor(pt_1,eta_1);
        }
	float correction=sf_id*sf_trg*LumiWeights_12->weight(npu);//sf_trg*sf_id*pu_weight;//*LumiWeights_12->weight(npu);//sf_trg*sf_id*h_Trk->Eval(eta_1)*LumiWeights_12->weight(npu);
	float totalweight=weight*correction;//sf_trg*sf_id;//*pu_weight;
        if (sample=="Data") totalweight=1.0;

        byLooseCombinedIsolationDeltaBetaCorr3Hits_2=byVLooseIsolationMVArun2v1DBoldDMwLT_2;
        byTightCombinedIsolationDeltaBetaCorr3Hits_2=byTightIsolationMVArun2v1DBoldDMwLT_2;
	if (pt_1>26 && pt_2>20){
	   if (iso_1<0.1 && byTightCombinedIsolationDeltaBetaCorr3Hits_2 && q_1*q_2<0){
		fillTreeMVA(arbre2,pt_1,pt_2,eta_1, eta_2,mt_1,mt_2,dphi_12, dphi_emet, dphi_taumet,met,njets,jdeta,mjj,m_coll,totalweight);
	   }
           //if ((iso_1<0.1 && byLooseCombinedIsolationDeltaBetaCorr3Hits_2 && !byTightCombinedIsolationDeltaBetaCorr3Hits_2 && q_1*q_2<0)){
           //     fillTreeMVA(arbre2,pt_1,pt_2,eta_1, eta_2,mt_1,mt_2,dphi_12, dphi_emet, dphi_taumet,met,njets,vbfeta,vbfmass,m_coll,frw/(1-frw));
           //}
	}

    } // end of loop over events
    TFile *fout = TFile::Open(output.c_str(), "RECREATE");
    fout->cd();
    arbre2->Write();
    fout->Close();
} 


