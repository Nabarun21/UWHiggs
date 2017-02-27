#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

//#include "/Users/cecilecaillol/ROOT/root-v5-34-00-patches/tmva/test/TMVAGui.C"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


using namespace TMVA;

float mjj,jdeta,pt_1,pt_2,eta_1, eta_2, deltaphi_12,mt_1,mt_2,deltaphi_taumet,deltaphi_emet,met,njets,vbfeta,vbfmass, phi_1, phi_2, metphi,dphi_12, dphi_emet, dphi_taumet, m_coll, ptrat, deltaeta_12;
//float KNNforETau(TMVA::Reader *reader, float pt_1_ = 20, float pt_2_ = 20, float deltaphi_12_ = 20, float mt_1_ = 20, float mt_2_ = 20, float deltaphi_taumet_ = 20, float deltaphi_emet_ = 20, float met_, float njets_, float vbfeta_, float vbfmass_) {
float KNNforETau(TMVA::Reader *reader, float pt_1_ = 20, float pt_2_ = 20, float mt_1_ = 20, float mt_2_ = 20, float dphi_12_=20, float dphi_taumet_=20,float met_=20, float njets_=20, float vbfmass_=20, float vbfeta_=20, float eta_1_=20, float eta_2_=20, float m_coll_=20){
    gROOT->ProcessLine(".O0"); // turn off optimization in CINT
    pt_1=pt_1_;
    pt_2=pt_2_;
    mt_1=mt_1_;
    mt_2=mt_2_;
    met=met_;
    dphi_12=dphi_12_;
    dphi_taumet=dphi_taumet_;
    njets=1.0*njets_;
    if (njets<2){
       vbfmass=-0.1;
       vbfeta=-0.1;
    }
    else {
       vbfmass=vbfmass_/2600;
       vbfeta=vbfeta_/10;
    }
    ptrat=pt_1_/pt_2_;
    deltaeta_12=eta_1_-eta_2_;
    m_coll=m_coll_;
    return reader->EvaluateMVA("BDT method");
}

