#include "mvaXXX.h"

int main(int argc, char** argv) {

    std::string input = *(argv + 1);
    using namespace std;

    TFile *f_Double = new TFile(input.c_str(),"update");
    TTree *Run_Tree = (TTree*) f_Double->Get("RLE_tree");

    float bdt_3;
    TBranch *bbdt = Run_Tree->Branch("bdt_3",&bdt_3,"bdt_3/F");
    Run_Tree->SetBranchAddress("pt_1", &pt_1);
    Run_Tree->SetBranchAddress("pt_2", &pt_2);
    Run_Tree->SetBranchAddress("dphi_12", &dphi_12);
    Run_Tree->SetBranchAddress("dphi_taumet", &dphi_taumet);
    Run_Tree->SetBranchAddress("dphi_emet", &dphi_emet);
    Run_Tree->SetBranchAddress("m_coll", &m_coll);
    Run_Tree->SetBranchAddress("met", &met);
    Run_Tree->SetBranchAddress("mt_1", &mt_1);
    Run_Tree->SetBranchAddress("mt_2", &mt_2);
    Run_Tree->SetBranchAddress("mjj", &mjj);
    Run_Tree->SetBranchAddress("jdeta", &jdeta);
    Run_Tree->SetBranchAddress("njets", &njets);
    Run_Tree->SetBranchAddress("eta_2", &eta_2);
    Run_Tree->SetBranchAddress("eta_1", &eta_1);
    Run_Tree->SetBranchAddress("m_coll", &m_coll);

     //######################################################
    TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");
    //reader->AddVariable("pt_1_", &pt_1);
    //reader->AddVariable("pt_2_", &pt_2);
    reader->AddVariable("mt_1_", &mt_1);
    reader->AddVariable("mt_2_", &mt_2);
    reader->AddVariable("dphi_12_", &dphi_12);
    reader->AddVariable("dphi_taumet_", &dphi_taumet);
    reader->AddVariable("met_", &met);
    reader->AddVariable("njets_", &njets);
    reader->AddVariable("vbfmass_", &vbfmass);
    reader->AddVariable("vbfeta_", &vbfeta);
    //reader->AddVariable("pt_1_/pt_2_", &ptrat);
    reader->AddVariable("eta_1_-eta_2_", &deltaeta_12);
    //reader->AddVariable("m_coll_", &m_coll);
    TString weightfile = "weights/TMVAClassification_BDT.weights.xml";
    TString methodName= TString("BDT") + TString(" method");
    reader->BookMVA(methodName, weightfile);

    //#########################################################
    Int_t nentries_wtn = (Int_t) Run_Tree->GetEntries();
    for (Int_t i = 0; i < nentries_wtn; i++) {
        Run_Tree->GetEntry(i);
        if (i % 100000 == 0) fprintf(stdout, "\r  Processed events: %8d of %8d ", i, nentries_wtn);
        fflush(stdout);
        bdt_3 = KNNforETau(reader, pt_1, pt_2, mt_1, mt_2, dphi_12, dphi_taumet, met, njets, mjj, jdeta, eta_1, eta_2, m_coll);
        bbdt->Fill();
    }

    Run_Tree->Write();
    delete f_Double;
}

