//include ROOT classes
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "THistPainter.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TTree.h"
#include "THStack.h"
#include "TColor.h"

//include C++ library classes
#include <sstream>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using std::cout;
using std::endl;
using std::flush;

//include other parts of the code
#include "MultilepSUSYfunc.h"
//#include "tdrstyle.h"
//#include "tdrstyle_official.h"
#include "tdrstyle.h"
#include "plotCode.h"


extern const Color_t colors[]; 
int main(int argc, char* argv[]){
    setTDRStyle();
    //TFile* PUfile = TFile::Open("puw_2016_13fb_200.root");
    //TH1D* PUweights = (TH1D*) PUfile->Get("puw");

    const int nSamples = 3;
    const TString fileList[nSamples] = {"WZTo3LNu.root", "WGToLNuG_NLO80.root", "WJetsToLNu.root"};
    const double xSections[nSamples] = {4.42965 , 489, 61526.7};  //or is WZ xsection 5.26?
    const TString names[nSamples] = {"WZ", "W#gamma", "WJets"};

    TFile *hfile[nSamples];
    TTree *inputTree[nSamples];
    TApplication* rootapp = new TApplication("example",&argc, argv);

    //Declare branches
    //Reco level branches
    TBranch *b__nL;
    TBranch *b__nph;
    TBranch *b__lPt;
    TBranch *b__lEta;
    TBranch *b__lPhi;
    TBranch *b__lE;	
    TBranch *b__phPt;
    TBranch *b__phEta;
    TBranch *b__phPhi;
    TBranch *b__phE;
    TBranch *b__phmva_pass;
    TBranch *b__isloose;
    TBranch *b__charges;
    TBranch *b__flavors;
    TBranch *b__met;
    TBranch *b__met_phi;
    TBranch *b__pheveto_pass;
    TBranch *b__phpixelseed_pass;
    TBranch *b__nJets;
    TBranch *b__jetEta;
    TBranch *b__jetPhi;
    TBranch *b__jetPt;
    TBranch *b__jetE;
    TBranch *b__weight;
    TBranch *b__isolation;
    TBranch *b__n_bJets;
    TBranch *b__multiisolation;
    TBranch *b__mu_ismedium;
    TBranch *b__eventNb;
    TBranch *b__runNb;
    TBranch *b__lepMVA;
    TBranch *b__clusterpass;
    TBranch *b__Mu17_Photon22_CaloIdL_L1ISO;
    TBranch *b__IsoMu20;
    TBranch *b__n_trueInteractions;
    TBranch *b__origin;
    //Gen lepton branches
    TBranch *b__gen_nL;
    TBranch *b__gen_lPt;
    TBranch *b__gen_lEta;
    TBranch *b__gen_lPhi;
    TBranch *b__gen_lE;
    TBranch *b__gen_lmompdg;
    TBranch *b__gen_flavors;
    TBranch *b__gen_charges;
    TBranch *b__genmet;
    TBranch *b__genmet_phi;
    //Gen branches for cleaning Wjets sample
    TBranch *b__gen_nPh;
    TBranch *b__gen_phPt;
    TBranch *b__gen_phEta;
    TBranch *b__gen_phPhi;
    TBranch *b__gen_phE;
    TBranch *b__gen_isPromptPh;
    TBranch *b__gen_phmompdg;
    TBranch *b__gen_isPromptl;

    //Declare leaf types
    const int nL_max = 10;
    const int nPh_max = 10;
    const int nJets_max = 10;
    //Gen level leaf types for neutrinos
    Double_t _genmet;
    Double_t _genmet_phi;
    Double_t _weight;
    //Reco level leaf types for event selection
    Int_t _nL;
    Int_t _nEle;
    Int_t _nMu;
    Int_t _nph;
    Double_t _lPt[nL_max];
    Double_t _lE[nL_max];
    Double_t _lEta[nL_max];
    Double_t _lPhi[nL_max];
    Double_t _phPt[nPh_max];
    Double_t _phEta[nPh_max];
    Double_t _phPhi[nPh_max];
    Double_t _phE[nPh_max];
    Bool_t _phmva_pass[nPh_max];
    Bool_t _isloose[nL_max];
    Int_t _flavors[nL_max];
    Double_t _charges[nL_max];
    Double_t _met;
    Double_t _met_phi;
    Bool_t 	_pheveto_pass[nPh_max];
    Bool_t 	_phpixelseed_pass[nPh_max];
    Int_t 	_nJets;
    Double_t _HT;
    Double_t _jetEta[nJets_max];
    Double_t _jetPhi[nJets_max];
    Double_t _jetPt[nJets_max];
    Double_t _jetE[nJets_max];
    Double_t _isolation[nL_max];
    Int_t _n_bJets;
    ULong64_t _runNb;
    Double_t _lepMVA[nL_max];
    Bool_t _clusterpass[nL_max];
    Bool_t _Mu17_Photon22_CaloIdL_L1ISO;
    Bool_t _IsoMu20;
    //Gen level lepton branches
    Int_t _gen_nL;
    Double_t _gen_lPt[nL_max];
    Double_t _gen_lEta[nL_max];
    Double_t _gen_lPhi[nL_max];
    Double_t _gen_lE[nL_max];
    Int_t _gen_lmompdg[nL_max];
    Int_t _gen_flavors[nL_max];
    Double_t _gen_charges[nL_max];
    //Gen level branches for cleaning Wjets sample
    Int_t _gen_nPh;
    Bool_t _gen_isPromptPh[nPh_max];
    Double_t _gen_phPt[nPh_max];
    Double_t _gen_phEta[nPh_max];
    Double_t _gen_phPhi[nPh_max];
    Double_t _gen_phE[nPh_max];
    Bool_t _multiisolation[nL_max][5];
    Bool_t _mu_ismedium[nL_max];
    ULong64_t _eventNb;
    Int_t _gen_phmompdg[nPh_max];
    Int_t _n_trueInteractions;
    Int_t _origin[nL_max];
    Bool_t	_gen_isPromptl[nL_max];


    double hcounter[nSamples];
    for(int sam = 0; sam < nSamples; ++sam){
        cout << "name " << names[sam] << endl;
        hfile[sam] = new TFile("../data80/"+fileList[sam],"read");
        hfile[sam]->cd("FakeElectrons");
        inputTree[sam] = static_cast<TTree*>(hfile[sam]->Get("FakeElectrons/fakeTree"));
        //Set Reco level branch addresses
        inputTree[sam]->SetBranchAddress("_nL", &_nL, &b__nL);
        inputTree[sam]->SetBranchAddress("_nph",&_nph, &b__nph);
        inputTree[sam]->SetBranchAddress("_lPt", &_lPt, &b__lPt);
        inputTree[sam]->SetBranchAddress("_lEta", &_lEta, &b__lEta);
        inputTree[sam]->SetBranchAddress("_lPhi", &_lPhi, &b__lPhi);
        inputTree[sam]->SetBranchAddress("_lE", &_lE, &b__lE);
        inputTree[sam]->SetBranchAddress("_phPt", &_phPt, &b__phPt);
        inputTree[sam]->SetBranchAddress("_phEta", &_phEta, &b__phEta);
        inputTree[sam]->SetBranchAddress("_phPhi", &_phPhi, &b__phPhi);
        inputTree[sam]->SetBranchAddress("_phE", &_phE, &b__phE);
        inputTree[sam]->SetBranchAddress("_charges", &_charges, &b__charges);
        inputTree[sam]->SetBranchAddress("_flavors", &_flavors, &b__flavors);
        inputTree[sam]->SetBranchAddress("_phmva_pass", &_phmva_pass, &b__phmva_pass);
        inputTree[sam]->SetBranchAddress("_isloose", &_isloose, &b__isloose);
        inputTree[sam]->SetBranchAddress("_met", &_met, &b__met);
        inputTree[sam]->SetBranchAddress("_met_phi", &_met_phi, &b__met_phi);
        inputTree[sam]->SetBranchAddress("_pheveto_pass",&_pheveto_pass, &b__pheveto_pass);
        inputTree[sam]->SetBranchAddress("_phpixelseed_pass", &_phpixelseed_pass, &b__phpixelseed_pass);
        inputTree[sam]->SetBranchAddress("_nJets", &_nJets, &b__nJets);
        inputTree[sam]->SetBranchAddress("_jetEta", &_jetEta, &b__jetEta);
        inputTree[sam]->SetBranchAddress("_jetPhi", &_jetPhi, &b__jetPhi);
        inputTree[sam]->SetBranchAddress("_jetPt", &_jetPt, &b__jetPt);
        inputTree[sam]->SetBranchAddress("_jetE", &_jetE, &b__jetE);
        inputTree[sam]->SetBranchAddress("_isolation", &_isolation, &b__isolation);
        inputTree[sam]->SetBranchAddress("_n_bJets", &_n_bJets, &b__n_bJets);
        inputTree[sam]->SetBranchAddress("_multiisolation", &_multiisolation, &b__multiisolation);
        inputTree[sam]->SetBranchAddress("_mu_ismedium", &_mu_ismedium, &b__mu_ismedium);
        inputTree[sam]->SetBranchAddress("_eventNb", &_eventNb, &b__eventNb);
        inputTree[sam]->SetBranchAddress("_lepMVA", &_lepMVA, &b__lepMVA);
        inputTree[sam]->SetBranchAddress("_clusterpass", &_clusterpass, &b__clusterpass);
        inputTree[sam]->SetBranchAddress("_runNb", &_runNb, &b__runNb);
        inputTree[sam]->SetBranchAddress("_Mu17_Photon22_CaloIdL_L1ISO", &_Mu17_Photon22_CaloIdL_L1ISO, &b__Mu17_Photon22_CaloIdL_L1ISO);
        inputTree[sam]->SetBranchAddress("_IsoMu20", &_IsoMu20, &b__IsoMu20);
        inputTree[sam]->SetBranchAddress("_n_trueInteractions", &_n_trueInteractions, &b__n_trueInteractions);
        inputTree[sam]->SetBranchAddress("_origin", &_origin, &b__origin);
        //Gen level lepton addresses
        inputTree[sam]->SetBranchAddress("_gen_nL", &_gen_nL, &b__gen_nL);
        inputTree[sam]->SetBranchAddress("_gen_lmompdg", &_gen_lmompdg, &b__gen_lmompdg);
        inputTree[sam]->SetBranchAddress("_gen_charges", &_gen_charges, &b__gen_charges);
        inputTree[sam]->SetBranchAddress("_gen_flavors", &_gen_flavors, &b__gen_flavors);
        inputTree[sam]->SetBranchAddress("_gen_lPt", &_gen_lPt, &b__gen_lPt);
        inputTree[sam]->SetBranchAddress("_gen_lPhi", &_gen_lPhi, &b__gen_lPhi);
        inputTree[sam]->SetBranchAddress("_gen_lEta", &_gen_lEta, &b__gen_lEta);
        inputTree[sam]->SetBranchAddress("_gen_lE", &_gen_lE, &b__gen_lE);
        inputTree[sam]->SetBranchAddress("_gen_isPromptl", &_gen_isPromptl, &b__gen_isPromptl);
        inputTree[sam]->SetBranchAddress("_genmet", &_genmet, &b__genmet);
        inputTree[sam]->SetBranchAddress("_genmet_phi", &_genmet_phi, &b__genmet_phi);

        //Gen level addresses for cleaning the Wjets sample
        if(sam != 0){
            inputTree[sam]->SetBranchAddress("_gen_nPh", &_gen_nPh, &b__gen_nPh);
            inputTree[sam]->SetBranchAddress("_gen_phPt", &_gen_phPt, &b__gen_phPt);
            inputTree[sam]->SetBranchAddress("_gen_phEta", &_gen_phEta, &b__gen_phEta);
            inputTree[sam]->SetBranchAddress("_gen_phPhi", &_gen_phPhi, &b__gen_phPhi);
            inputTree[sam]->SetBranchAddress("_gen_phE", &_gen_phE, &b__gen_phE);
            inputTree[sam]->SetBranchAddress("_gen_isPromptPh", &_gen_isPromptPh, &b__gen_isPromptPh);
            inputTree[sam]->SetBranchAddress("_gen_phmompdg", &_gen_phmompdg, &b__gen_phmompdg);
        }		
        inputTree[sam]->SetBranchAddress("_weight", &_weight, &b__weight);
        TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
        _hCounter->Read("hCounter");
        hcounter[sam] = _hCounter->GetBinContent(1);
    }
    //Tweakable options////////////////////////////////////////////////////
    const bool DeltaRReweighing = false;
    const bool HTReweighing = false;
    const bool BosonPtReweighing = false;
    const bool leptonPtReweighing = false;
    const bool RecoleptonPtReweighing = false;
    const bool DeltaPhiReweighing = false;
    //const bool METreweighing = true;
    const bool TestRun = false;	//Break after a few events
    const bool SkipMuons = false;
    const bool SkipElectrons = true;
    const double MuMVACut = -0.2;
    const double EleMVACut = 0.5;
    const bool NeedPixSeed = true;
    const bool Pheveto = false;
    const bool bVeto = true;
    const bool SkipWZ = false;
    const bool SkipWgamma = false;
    const bool SkipWjets = true;
    const int nMETbins = 4;
    const int nPtBins = 4;
    const int nPhiBins = 4;
    const double MaxLeptonPtCut_Mu = 20; //25 for muons, 30 for electrons
    const double MaxLeptonPtCut_Ele = 35;
    const double MidLeptonPtCut = 10;
    //const double MidLeptonPtCutEle = 10;
    const double MinLeptonPtCut = 10;
    const double PhotonPtCut = 50;
    const double DeltaRCut = 0.3;				
    const double DeltaPhiCut = 1;
    const double DeltaEtaCut = 0;
    const double METCut = 50;
    const double MTCut = 0;
    const double MTMax = 300;
    const double DataLuminosity = 12.9;		//in units of fb^-1
    const TString extra = "_Mu_DeltaMET60";
    ///////////////////////////////////////////////////////////////////////

    TString lepton;
    if(SkipMuons) lepton = "e";
    else if(SkipElectrons) lepton = "#mu";
    else lepton = "lepton";
    double MaxLeptonPtCut;
    if(SkipElectrons) MaxLeptonPtCut = MaxLeptonPtCut_Mu;
    else if(SkipMuons) MaxLeptonPtCut = MaxLeptonPtCut_Ele;
    TH1D* Histos[nSamples][12];
    const TString Histnames[12] = {"MT", "BosonPt", "LeptonPt", "MET", "DeltaR", "DeltaPhi", "HT", "NJets", "NbJets", "DeltaPhi_lepMET", "DeltaEta", "Mlgamma"};
    const TString Xaxes[12] = { "M_{T}(" + lepton + " + MET) (GeV)", "P_{T}(#gamma) (GeV)", "P_{T}(" + lepton + ") (GeV)",  "MET (GeV)", "#DeltaR(" + lepton + ", #gamma)", "#Delta#Phi(" + lepton + ", #gamma)", "HT (GeV)", "number of jets", "number of b-jets", "#Delta#Phi(" + lepton + ", MET)", "#Delta#eta(" + lepton + ", #gamma)", "M_{" + lepton + "#gamma} (GeV)"};
    const TString Units[12] = {"GeV", "GeV", "GeV", "GeV", "", "", "GeV", "", "", "", "", "GeV"};
    const double HistMin[12] = {MTCut, PhotonPtCut, MaxLeptonPtCut, METCut, DeltaRCut, DeltaPhiCut, 30, 0, 0, 0, 0, 30};
    const double HistMax[12] = {MTMax, 300,150,250, 7, 3.2, 600, 10,10, 3.2, 5, 200};
    const int nBins[12] = {16, 20, 10, 20, 15, 15, 20, 10, 10, 15, 15, 30};
    for(int i = 0; i < 12; ++i){
        float BinWidth = (HistMax[i] - HistMin[i])/nBins[i];
        std::ostringstream strs; strs << BinWidth; std::string Yaxis = strs.str();
        for(int sam = 0; sam < nSamples; ++sam){
            Histos[sam][i] = new TH1D(Histnames[i] + names[sam], Histnames[i] + names[sam] + ";" + Xaxes[i] + "; events /" + Yaxis + Units[i], nBins[i], HistMin[i], HistMax[i]);
        }
    }
    TH1D* ReweighingHistos[11];
    TH1D* PTre;
    TH1D* Phire;

    const TString PTgbins[3] ={"$P_{T}(\\gamma) > 50$", "$P_{T}(\\gamma) > 70$", "$P_{T}(\\gamma) > 90$"};
    TH1D* MT_PTg[3];
    for(int i = 0; i < 3; ++i) MT_PTg[i] = (TH1D*) Histos[1][0]->Clone(PTgbins[i]);
    const TString Scalfluc[3] ={"+0 $sigma$", "+ 1 $sigma$", " -1 $sigma$"};
    TH1D* MT_Scalfluc[3];
    for(int i = 0; i < 3; ++i) MT_Scalfluc[i] = (TH1D*) Histos[1][0]->Clone(Scalfluc[i]);
    TH1D* PTre_up;
    TH1D* PTre_down;

    //2D MT MET histogram for SF's
    double MTbins[4] = {40, 120,160, 200};
    double METbins[5] = {50, 100, 150, 200, 250};
    TH2D *MTMET[2];
    for(int i = 0; i < 2; ++i){
        MTMET[i] = new TH2D("MTMET_" + names[i], "MTMET_" + names[i] +  Xaxes[0] + ";" + Xaxes[3], 3,MTbins , 4, METbins);
        MTMET[i]->Sumw2();
    }

    TH1D* WZPU[2];
    float BinWidth = (HistMax[0] - HistMin[0])/nBins[0];
    std::ostringstream strs; strs << BinWidth; std::string Yaxis = strs.str();
    WZPU[0] = new TH1D(Histnames[0] + "WZPU", Histnames[0] + "WZPU" + ";" + Xaxes[0] + "; events /" + Yaxis + Units[0], nBins[0], HistMin[0], HistMax[0]);
    BinWidth = (HistMax[3] - HistMin[3])/nBins[3];
    std::ostringstream strs2; strs2 << BinWidth; Yaxis = strs.str();
    WZPU[1] = new TH1D(Histnames[3] + "WZPU", Histnames[3] + "WZPU" + ";" + Xaxes[3] + "; events /" + Yaxis + Units[3], nBins[3], HistMin[3], HistMax[3]);

    TH2D* MTMETPU;
    MTMETPU = new TH2D("MTMET_PU", "MTMET_PU" +  Xaxes[0] + "(GeV)" + ";" + Xaxes[3] + "(GeV)", 3,MTbins , 4, METbins);

    //2 Histograms for splitting the WZ MT in events with correct and events with incorrect pairing	
    TH1D* WZ_split[6]; 
    const TString source[6] = {"correct_pairing", "wrong_pairing", "tau", "other_Z", "other", "fakes"};
    for(int i = 0; i < 6; ++i){
        float BinWidth = (HistMax[0] - HistMin[0])/nBins[0];
        std::ostringstream strs; strs << BinWidth; std::string Yaxis = strs.str();
        WZ_split[i] = new TH1D(source[i], source[i] + ";" + Xaxes[0] + "; events /" + Yaxis + Units[0], nBins[0], HistMin[0], HistMax[0]);
    }
	TH1D* WZ_METres[4];
	const TString metb[4] = {"20", "40", "100", "inf"};
	for(int i = 0; i < 4; ++i){
		float BinWidth = (HistMax[0] - HistMin[0])/nBins[0];
        std::ostringstream strs; strs << BinWidth; std::string Yaxis = strs.str();
        WZ_METres[i] = new TH1D(metb[i], metb[i] + ";" + Xaxes[0] + "; events /" + Yaxis + Units[0], nBins[0], HistMin[0], HistMax[0]);
    }
    TH1D* Mlltest[2];
    const TString pairing[2] = {"Z", "W + OS"};
    for(int i = 0; i < 2; ++i){
        float BinWidth = (0 - 300)/nBins[0];
        std::ostringstream strs; strs << BinWidth; std::string Yaxis = strs.str();
        Mlltest[i] = new TH1D(pairing[i], pairing[i] + ";" "M_{ll} ; events /" + Yaxis + Units[0], 60, HistMin[0], HistMax[0]);
    }	
    Double_t scale[nSamples];
    Double_t Total[nSamples];
    Double_t Pass[nSamples];
    //Text file for making dump
    std::ofstream dumpFile[nSamples];
    for(int i = 0;  i < nSamples; ++i){
        dumpFile[i].open(std::string(names[i]) + "_dump.txt");
    }
    unsigned dumpCount[nSamples];
    /////////////////////////////
    for(int r = 0; r < 2; ++r){ 
        if(!(BosonPtReweighing || HTReweighing ||DeltaRReweighing || leptonPtReweighing || RecoleptonPtReweighing || DeltaPhiReweighing) ) ++r;					//loop for reweighing
        for(int sam = 0; sam < nSamples; ++sam){
            //make sure to only loop over samples once if no reweighing is used
            if( (BosonPtReweighing || HTReweighing || DeltaRReweighing || leptonPtReweighing || RecoleptonPtReweighing || DeltaPhiReweighing) && r == 1 && sam == 0) continue;
            if(SkipWZ && names[sam] == "WZ") continue;
            if(SkipWjets && names[sam] == "WJets") continue;
            if(SkipWgamma &&  names[sam] == "WGToLNuG") continue;
            Total[sam] = 0;
            Pass[sam] = 0;
            scale[sam] = xSections[sam]*DataLuminosity*1000/(hcounter[sam]);
            Long64_t nEntries = inputTree[sam]->GetEntries();
            std::cout<<"Entries in "<< fileList[sam] <<" "<<nEntries<<std::endl;
				

            //Loop over all events in the sample
            for (Long64_t it = 0; it < nEntries; ++it){
                inputTree[sam]->GetEntry(it);
                if (it%10000 == 0) cout<<'.'<<flush;
                if(TestRun && it > 10000) break;

                //indices representing the position of the leptons from W and Z decay in the arrays
                int lz1 = -1, lz2 = -1, lw = -1;
                //index we use for the selected photon
                int ph = -1;

                //Clean Wjets sample of overlap with Wgamma
                if(sam == 2){
                    bool promptfail = false;
                    for(int j = 0; j < _gen_nPh; ++j){
                        if(fabs(_gen_phmompdg[j]) != 111 && (_gen_phPt[j] > 10) ){
                            promptfail = true;
                            break;
                        }
                    }
                    if(promptfail) continue;
                }
                //Kinematic cut on MET
                if(_met < METCut) continue;
                //Other kinematic cuts depend on the sample

                if(sam == 0){
                    const int nL = _nL;
                    int ind[nL];
                    int lepCounter = 0;
                    for(int j = 0; j < _nL; ++j){
                        if(_flavors[j] == 2) continue;
                        else if(_flavors[j] == 1){
                            if(!_isloose[j]) continue;
                            if(!_mu_ismedium[j]) continue;
                            if(_lepMVA[j] < MuMVACut) continue;
                            ind[lepCounter] = j;
                        }
                        else if(_flavors[j] == 0){
                            if(!_isloose[j]) continue;
                            if(!_clusterpass[j]) continue;
                            if(_lepMVA[j] < EleMVACut) continue;
                            ind[lepCounter] = j;
                        }
                        ++lepCounter;
                    }
                    if(lepCounter != 3) continue;
                    //look for OSSF pair(Z) and a different flavor lepton(W)
                    if(_flavors[ind[0]] == _flavors[ind[1]] && _flavors[ind[0]] != _flavors[ind[2]]){
                        lw = ind[2];
                        lz1 = ind[0];
                        lz2 = ind[1];
                    }
                    else if(_flavors[ind[0]] != _flavors[ind[1]] && _flavors[ind[0]] == _flavors[ind[2]]){
                        lw = ind[1];
                        lz1 = ind[0];
                        lz2 = ind[2];
                    }
                    else if(_flavors[ind[0]] != _flavors[ind[1]] && _flavors[ind[1]] == _flavors[ind[2]]){
                        lw = ind[0];
                        lz1 = ind[1];
                        lz2 = ind[2];
                    }
                    //Code to deal with 3mu or 3ele events with an OSSF pair
                    else if(_flavors[ind[0]] == _flavors[ind[1]] && _flavors[ind[0]] == _flavors[ind[2]]){
                        if(_charges[ind[0]] == _charges[ind[1]] && _charges[ind[0]] == _charges[ind[2]]) continue;
                        TLorentzVector leptons[3];
                        for(int n = 0; n < 3; ++n){
                            leptons[n].SetPtEtaPhiE(_lPt[ind[n]], _lEta[ind[n]], _lPhi[ind[n]], _lE[ind[n]]);
                        }
                        const double Zmass = 91.1876;
                        double Mindif = 99999;
                        for(int j = 0; j < 2 ; ++j){
                            for(int n = j + 1; n < 3; ++n){
                                if(_charges[ind[j]] == _charges[ind[n]]) continue;
                                if( fabs((leptons[j] + leptons[n]).Mag() - Zmass) < Mindif){
                                    Mindif = fabs((leptons[j] + leptons[n]).Mag() - Zmass);
                                    lz1 = ind[j];
                                    lz2 = ind[n];
                                }
                            }
                        }
                        if(lz1 == -1 || lz2 == -1) continue;  //No OSSF pair found!
                        for(int j = 0; j < 3; ++j){
                            if( ind[j] == lz1 || ind[j] == lz2) continue;
                            lw = ind[j];
                        }
                    }
                    else continue;
                    //WARNING REMOVE THIS
                    //if(_flavors[lw] != 1 || _flavors[lz1] != 1 || _flavors[lz2] != 1) continue;

                    //extra checks
                    if(_charges[lz1] == _charges[lz2]) continue;
                    if(_flavors[lz1] == 2 || _flavors[lw] == 2 || _flavors[lz2] == 2) continue;
                    if(SkipMuons && _flavors[lw] == 1) continue;
                    if(SkipElectrons && _flavors[lw] == 0) continue;

                    //Apply Pt cuts
                    int midind;
                    double MinPt = 0, MaxPt = 0, MidPt = 0;
                    if(_lPt[lw] > _lPt[lz1] && _lPt[lw] > _lPt[lz2]){
                        MaxPt = _lPt[lw];
                        if(_lPt[lz1] > _lPt[lz2]){
                            MidPt = _lPt[lz1];
                            MinPt = _lPt[lz2];
                            midind = lz1;
                        }
                        else{
                            MidPt = _lPt[lz2];
                            MinPt = _lPt[lz1];
                            midind = lz2;
                        }
                    }
                    else if(_lPt[lz1] > _lPt[lw] && _lPt[lz1] > _lPt[lz2]){
                        MaxPt = _lPt[lz1];
                        if(_lPt[lw] > _lPt[lz2]){
                            MidPt = _lPt[lw];
                            MinPt = _lPt[lz2];
                            midind = lw;
                        }
                        else{
                            MidPt = _lPt[lz2];
                            MinPt = _lPt[lw];
                            midind = lz2;
                        }
                    }
                    else{
                        MaxPt = _lPt[lz2];
                        if(_lPt[lw] > _lPt[lz1]){
                            MidPt = _lPt[lw];
                            MinPt = _lPt[lz1];
                            midind = lw;
                        }
                        else{
                            MidPt = _lPt[lz1];
                            MinPt = _lPt[lw];
                            midind = lz1;
                        }
                    }
                    if(MinPt < MinLeptonPtCut || (_flavors[midind] == 1 && MidPt < MidLeptonPtCut) || (_flavors[midind] == 0 && MidPt < MidLeptonPtCut)  || MaxPt < MaxLeptonPtCut) continue;

                    /*
                    //if( (_flavors[lw] == 1 && _lPt[lw] > MaxLeptonPtCut_Mu)  || _lPt[lw] < 10) continue;
                    if(_flavors[lw] == 1 && _lPt[lw] < MaxLeptonPtCut_Mu) continue;
                    //else if(_flavors[lw] == 0 && _lPt[lw] > MaxLeptonPtCut_Ele || _lPt[lw] < 10) continue;
                    else if(_flavors[lw] == 0 && _lPt[lw] < MaxLeptonPtCut_Ele) continue;
                    double MidPt, MinPt;
                    if(_lPt[lz1] > _lPt[lz2]){
                    MidPt = _lPt[lz1];
                    MinPt = _lPt[lz2];
                    } else{
                    MidPt = _lPt[lz2];	
                    MinPt = _lPt[lz1];
                    }
                    if(MidPt < MidLeptonPtCut || MinPt < MinLeptonPtCut) continue;
                     */

                }
                //Event selection for Wgamma and Wjets samples
                else{
                    //Cuts on the lepton
                    //Start lepton Selection
                    int lepCounter = 0;
                    for(int j = 0; j < _nL; ++j){
                        if(_flavors[j] == 2) continue;
                        else if(_flavors[j] == 1){
                            if( !_mu_ismedium[j]) continue;
                            if( !_isloose[j]) continue;
                            if( _lepMVA[j] <  -0.2) continue;
                            lw = j;
                        } else if(_flavors[j] == 0){
                            if( !_isloose[j]) continue;	
                            if( !_clusterpass[j]) continue;
                            if( _lepMVA[j] < 0.5) continue;
                            lw = j;
                        }
                        ++lepCounter;
                    }
                    if(lepCounter != 1) continue;
                    if(SkipMuons && _flavors[lw] == 1) continue;
                    if(SkipElectrons && _flavors[lw] == 0) continue;
                    if(_lPt[lw] < MaxLeptonPtCut) continue;

                    double MaxPhPt = 0;				
                    for(int j = 0; j <_nph; ++j){
                        if(!_phmva_pass[j]) continue;
                        if(Pheveto && !_pheveto_pass[j]) continue;
                        if(NeedPixSeed && _phpixelseed_pass[j]) continue;
                        if(fabs(_phEta[j]) > 2.5) continue;
                        if(_phPt[j] > MaxPhPt){
                            ph = j;
                            MaxPhPt = _phPt[j];
                        }
                    }					
                    if(ph == -1) continue;
                    if(MaxPhPt < PhotonPtCut) continue;
                    //if(_lPt[lw] < MaxLeptonPtCut) continue;
                    if(_flavors[lw] == 1 &&_lPt[lw] < MaxLeptonPtCut_Mu) continue;
                    if(_flavors[lw] == 0 &&_lPt[lw] < MaxLeptonPtCut_Ele) continue;
                    TLorentzVector lepton, photon;
                    lepton.SetPtEtaPhiE(_lPt[lw], _lEta[lw], _lPhi[lw], _lE[lw]);
                    photon.SetPtEtaPhiE(_phPt[ph], _phEta[ph], _phPhi[ph], _phE[ph]);
                    //if(SkipMuons && (lepton + photon).Mag() > 75 && (lepton + photon).Mag() < 105) continue;
                    if(lepton.DeltaR(photon) < DeltaRCut) continue;
                    if(fabs(lepton.DeltaPhi(photon)) < DeltaPhiCut) continue;
                    if(fabs(_lEta[lw] - _phEta[ph]) < DeltaEtaCut) continue;
                }

                int nJets = 0;
                double HT = 0;
                for(int i = 0; i < _nJets; ++i){
                    if(fabs(_jetEta[i]) > 2.5 || _jetPt[i] < 30) continue;
                    TLorentzVector Jet;
                    Jet.SetPtEtaPhiE(_jetPt[i], _jetEta[i], _jetPhi[i], _jetE[i]);
                    bool leptonoverlap = false;
                    for(int j = 0; j < _nL; ++j){
                        TLorentzVector lepton;
                        lepton.SetPtEtaPhiE(_lPt[j], _lEta[j], _lPhi[j], _lE[j]);
                        if(Jet.DeltaR(lepton) < -0.2){
                            leptonoverlap = true;
                            break;
                        }
                    }
                    if(leptonoverlap) continue;
                    if(sam != 0){
                        TLorentzVector photon;
                        photon.SetPtEtaPhiE(_phPt[ph], _phEta[ph], _phPhi[ph], _phE[ph]);
                        if(Jet.DeltaR(photon) < 0.4) continue;
                    }
                    ++nJets;
                    HT += fabs(_jetPt[i]);
                }

                TLorentzVector Wlepton;
                //WARNING WARNING TEMPORARAY
                //CHANGE LZ1 BACK TO LW !!!!!!!!!!!!!!!!!!!!
                Wlepton.SetPtEtaPhiE(_lPt[lw], _lEta[lw], _lPhi[lw], _lE[lw]);
                TLorentzVector METvec;
                METvec.SetPtEtaPhiE(_met, 0, _met_phi, _met);
                TLorentzVector Boson;
                if(sam == 0){
                    TLorentzVector Zlep1, Zlep2;
                    Zlep1.SetPtEtaPhiE(_lPt[lz1], _lEta[lz1], _lPhi[lz1], _lE[lz1]);
                    Zlep2.SetPtEtaPhiE(_lPt[lz2], _lEta[lz2], _lPhi[lz2], _lE[lz2]);
                    Boson = Zlep1 + Zlep2;
                    if(Boson.Mag() < 12 ) continue;  //values below are not simulated
                    //if( (Boson + Wlepton).Mag()  > 75 && (Boson + Wlepton).Mag() < 105) continue;
                }
                else{
                    Boson.SetPtEtaPhiE(_phPt[ph], _phEta[ph], _phPhi[ph], _phE[ph]);
                }

                //Cuts on MT and the number of bjets
                if(transmass(Wlepton,METvec) < MTCut) continue;
                if( bVeto && _n_bJets > 0) continue;


                //Generator matching
                //lepton generator matching

                const int nL = _nL;
                const int gen_nL = _gen_nL;
                TLorentzVector Lep[nL];
                TLorentzVector GenLep[gen_nL];
                for(int i = 0; i < _nL; ++i){
                    Lep[i].SetPtEtaPhiE(_lPt[i], _lEta[i], _lPhi[i], _lE[i]);
                }
                for(int i = 0; i < _gen_nL; ++i){
                    GenLep[i].SetPtEtaPhiE(_gen_lPt[i], _gen_lEta[i], _gen_lPhi[i], _gen_lE[i]);
                }
                int match[nL];
                for(int i = 0; i < _nL; ++i){
                    match[i] = -1;
                }

                std::vector<int> UsedGen;
                std::vector<int> UsedReco;
                for(int n = 0; n < _nL; ++n){
                    double MinDeltaR = 999;
                    int matched = -1;
                    for(int i = 0; i < _nL; ++i){
                        if(_lPt[i] == 0) continue;
                        bool recoused = false;
                        for(int j = 0; j < UsedReco.size(); ++j){
                            if(i == UsedReco[j]){
                                recoused = true;
                                break;
                            }
                        }
                        if(recoused) continue;
                        for(int j = 0; j < _gen_nL; ++j){
                            if(_gen_lPt[j] == 0) continue;
                            if(_charges[i] != _gen_charges[j] || _flavors[i] != _gen_flavors[j]) continue;
                            bool genused = false;
                            for(int k = 0; k < UsedGen.size(); ++k){
                                if(j == UsedGen[k]){
                                    genused = true;
                                    break;
                                }
                            }
                            if(genused) continue;
                            if(Lep[i].DeltaR(GenLep[j]) < MinDeltaR){
                                matched = i;
                                match[i] = j;
                                MinDeltaR = Lep[i].DeltaR(GenLep[j]);
                            }
                        }
                    }
                    UsedReco.push_back(matched);
                    UsedGen.push_back(match[matched]);                  
                }
                UsedGen.clear();
                UsedReco.clear();

                int ph_match = -1;
                if(sam != 0){
                    //photon generator matching
                    const int nPh = _nph;
                    const int gen_nph = _gen_nPh;
                    TLorentzVector photon;
                    photon.SetPtEtaPhiE(_phPt[ph], _phEta[ph], _phPhi[ph], _phE[ph]);
                    TLorentzVector genPhoton[gen_nph];
                    double MinDeltaR = 999;
                    for(int i = 0; i < _gen_nPh; ++i){
                        genPhoton[i].SetPtEtaPhiE(_gen_phPt[i], _gen_phEta[i], _gen_phPhi[i], _gen_phE[i]);
                        if(photon.DeltaR(genPhoton[i]) < MinDeltaR){
                            MinDeltaR = photon.DeltaR(genPhoton[i]);
                            ph_match = i;
                        }
                    }
                }
				//REMOVE this later
				//TLorentzVector genMETvec;
				//genMETvec.SetPtEtaPhiE(_genmet, 0, _genmet_phi, _genmet);
				//if(METvec.DeltaPhi(genMETvec) > 0.2) continue;
				if(fabs(_genmet -_met) < 40) continue;


                Pass[sam] +=scale[sam]*_weight;
                ++Total[sam];
                if(transmass(Wlepton,METvec) > 250 && dumpCount[sam] < 50){
                    TLorentzVector GenMETvec, GenWlepton;
                    GenMETvec.SetPtEtaPhiE(_genmet, 0, _genmet_phi, _met);
                    GenWlepton.SetPtEtaPhiE(_gen_lPt[match[lw]], _gen_lEta[match[lw]], _gen_lPhi[match[lw]], _gen_lE[match[lw]]);
                    dumpFile[sam]<< "===============" << names[sam] << "======================" << endl;
                    dumpFile[sam]<< "lepton Pt : " << _lPt[lw] << endl;
                    dumpFile[sam]<< "generator lepton Pt : " << _gen_lPt[match[lw]] << endl;
                    dumpFile[sam]<< "lepton mother : " << _gen_lmompdg[match[lw]] << endl;
                    dumpFile[sam]<< "lepton origin : " << _origin[lw] << endl;
                    if(sam !=0){
                        dumpFile[sam]<< "photon Pt : " << _phPt[ph] << endl;
                        dumpFile[sam]<< "generator photon Pt : " << _gen_phPt[ph_match] << endl;
                        dumpFile[sam]<< "photon mother : " << _gen_phmompdg[ph_match] << endl;
                    }

                    dumpFile[sam]<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
                    dumpFile[sam]<< "MET : " << _met << endl;
                    dumpFile[sam]<< "generator MET : " << _genmet << endl;
                    dumpFile[sam]<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
                    dumpFile[sam]<< "Transverse mass : " << transmass(Wlepton,METvec) << endl;
                    dumpFile[sam]<< "genetator transverse mass : " << transmass(GenWlepton,GenMETvec) << endl;
                    if(sam == 0){
                        dumpFile[sam]<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
                        dumpFile[sam]<< "Z lepton 1 mother : " << _gen_lmompdg[match[lz1]] << endl;
                        dumpFile[sam]<< "Z lepton 2 mother : " << _gen_lmompdg[match[lz2]] << endl;
                        TLorentzVector genZlep[2];
                        genZlep[0].SetPtEtaPhiE(_gen_lPt[match[lz1]], _gen_lEta[match[lz1]], _gen_lPhi[match[lz1]], _gen_lE[match[lz1]]);
                        genZlep[1].SetPtEtaPhiE(_gen_lPt[match[lz2]], _gen_lEta[match[lz2]], _gen_lPhi[match[lz2]], _gen_lE[match[lz2]]);
                        dumpFile[sam] << "generator Z Pt : " << (genZlep[0] + genZlep[1]).Pt() << endl;
                    }
                    dumpFile[sam]<< "=========================================================" << endl;
                    ++dumpCount[sam];
                }
	
				//if(sam == 0 && _gen_isPromptl[match[lw]] && fabs(_gen_lmompdg[lw]) != 24 && fabs(_gen_lmompdg[match[lw]]) != 15) continue;

                //Reweighing
                int effsam;
                if(sam != 2) effsam = sam;
                else effsam = 1;

                //if(RecoleptonPtReweighing && r == 1) _weight *= ReweighingHistos[1]->GetBinContent(ReweighingHistos[1]->FindBin(TMath::Min(_lPt[lw], ReweighingHistos[1]->GetBinCenter(ReweighingHistos[1]->GetNbinsX()) )));
                //if(_flavors[lw] == 1) _weight*= LeptonPtRatio_mu->GetBinContent(LeptonPtRatio_mu->FindBin(TMath::Min(_lPt[lw], LeptonPtRatio_mu->GetBinCenter(LeptonPtRatio_mu->GetNbinsX()) )));
                //else if(_flavors[lw] == 0) _weight*= LeptonPtRatio_ele->GetBinContent(LeptonPtRatio_ele->FindBin(TMath::Min(_lPt[lw], LeptonPtRatio_ele->GetBinCenter(LeptonPtRatio_ele->GetNbinsX()) )));

                if(r > 0 && sam ==1 && RecoleptonPtReweighing){
                    double tempweight[3] = {_weight, _weight, _weight};
                    tempweight[0] *= PTre->GetBinContent(PTre->FindBin(TMath::Min(_lPt[lw], PTre->GetBinCenter(PTre->GetNbinsX()) )));
                    tempweight[1] *= PTre_up->GetBinContent(PTre_up->FindBin(TMath::Min(_lPt[lw], PTre_up->GetBinCenter(PTre_up->GetNbinsX()) )));
                    tempweight[2] *= PTre_down->GetBinContent(PTre_down->FindBin(TMath::Min(_lPt[lw], PTre_down->GetBinCenter(PTre_down->GetNbinsX()) )));
                    for(int i = 0; i < 3; ++i){
                        MT_Scalfluc[i]->Fill(TMath::Min(transmass(Wlepton, METvec),Histos[effsam][0]->GetBinCenter(Histos[effsam][0]->GetNbinsX())), scale[sam]*tempweight[i]);
                    }
                }

                if(RecoleptonPtReweighing && r == 1) _weight *= PTre->GetBinContent(PTre->FindBin(TMath::Min(_lPt[lw], PTre->GetBinCenter(PTre->GetNbinsX()) )));
                if(DeltaPhiReweighing && r == 1) _weight *= Phire->GetBinContent(Phire->FindBin(TMath::Min(Wlepton.DeltaPhi(METvec), Phire->GetBinCenter(Phire->GetNbinsX()) )));

                Histos[effsam][1]->Fill(TMath::Min(_phPt[ph],Histos[effsam][1]->GetBinCenter(Histos[effsam][1]->GetNbinsX())), scale[sam]*_weight);
                Histos[effsam][2]->Fill(TMath::Min(_lPt[lw],Histos[effsam][2]->GetBinCenter(Histos[effsam][2]->GetNbinsX())), scale[sam]*_weight);
                Histos[effsam][3]->Fill(TMath::Min(_met,Histos[effsam][3]->GetBinCenter(Histos[effsam][3]->GetNbinsX())), scale[sam]*_weight);
                Histos[effsam][4]->Fill(TMath::Min( Boson.DeltaR(Wlepton),Histos[effsam][4]->GetBinCenter(Histos[effsam][4]->GetNbinsX())), scale[sam]*_weight);
                Histos[effsam][5]->Fill(TMath::Min( fabs(Boson.DeltaPhi(Wlepton)),Histos[effsam][5]->GetBinCenter(Histos[effsam][5]->GetNbinsX())), scale[sam]*_weight);
                Histos[effsam][6]->Fill(TMath::Min(HT,Histos[effsam][6]->GetBinCenter(Histos[effsam][6]->GetNbinsX())), scale[sam]*_weight);
                Histos[effsam][7]->Fill(nJets, scale[sam]*_weight);
                Histos[effsam][8]->Fill(_n_bJets, scale[sam]*_weight);
                Histos[effsam][9]->Fill(TMath::Min(Wlepton.DeltaPhi(METvec),Histos[effsam][9]->GetBinCenter(Histos[effsam][9]->GetNbinsX())), scale[sam]*_weight);
                Histos[effsam][10]->Fill(TMath::Min(fabs(_lEta[lw] - _phEta[ph]),Histos[effsam][10]->GetBinCenter(Histos[effsam][10]->GetNbinsX())), scale[sam]*_weight);
                Histos[effsam][11]->Fill(TMath::Min( (Wlepton + Boson).Mag(),Histos[effsam][11]->GetBinCenter(Histos[effsam][11]->GetNbinsX())), scale[sam]*_weight);
                if(r > 0 && sam == 1){
                    MT_PTg[0]->Fill(TMath::Min(transmass(Wlepton, METvec),MT_PTg[0]->GetBinCenter(MT_PTg[0]->GetNbinsX())), scale[sam]*_weight);
                    if(_phPt[ph] > 70) MT_PTg[1]->Fill(TMath::Min(transmass(Wlepton, METvec), MT_PTg[1]->GetBinCenter(MT_PTg[1]->GetNbinsX())), scale[sam]*_weight);
                    if(_phPt[ph] > 90) MT_PTg[2]->Fill(TMath::Min(transmass(Wlepton, METvec), MT_PTg[2]->GetBinCenter(MT_PTg[2]->GetNbinsX())), scale[sam]*_weight);	
                }
                if(r > 0 || effsam == 0 && transmass(Wlepton, METvec) > MTCut){
                    Histos[effsam][0]->Fill(TMath::Min(transmass(Wlepton, METvec),Histos[effsam][0]->GetBinCenter(Histos[effsam][0]->GetNbinsX())), scale[sam]*_weight);
                    MTMET[effsam]->Fill(TMath::Min(transmass(Wlepton, METvec),199.), TMath::Min(_met,249.), scale[sam]*_weight);
                    //_weight*= PUweights->GetBinContent(PUweights->FindBin(_n_trueInteractions));
                    if(effsam == 0){
                        WZPU[0]->Fill(TMath::Min(transmass(Wlepton, METvec), WZPU[0]->GetBinCenter( WZPU[0]->GetNbinsX())), scale[sam]*_weight);
                        WZPU[1]->Fill(TMath::Min(_met, WZPU[1]->GetBinCenter( WZPU[1]->GetNbinsX())), scale[sam]*_weight);
                        MTMETPU->Fill(TMath::Min(transmass(Wlepton, METvec),199.), TMath::Min(_met,249.), scale[sam]*_weight);
                    }
                    //Filling histograms with split WZ background
                    if(sam == 0){
                        if( fabs(_gen_lmompdg[match[lw]]) == 24){
                            WZ_split[0]->Fill(TMath::Min(transmass(Wlepton, METvec), WZ_split[0]->GetBinCenter(WZ_split[0]->GetNbinsX())), scale[sam]*_weight);
                        } else if(fabs(_gen_lmompdg[match[lw]]) == 23 ){
                            WZ_split[1]->Fill(TMath::Min(transmass(Wlepton, METvec), WZ_split[1]->GetBinCenter(WZ_split[1]->GetNbinsX())), scale[sam]*_weight);
                            TLorentzVector Zlep1, Zlep2, Wlep;
                            Zlep1.SetPtEtaPhiE(_lPt[lz1], _lEta[lz1], _lPhi[lz1], _lE[lz1]);
                            Zlep2.SetPtEtaPhiE(_lPt[lz2], _lEta[lz2], _lPhi[lz2], _lE[lz2]);
                            Wlep.SetPtEtaPhiE(_lPt[lw], _lEta[lw], _lPhi[lw], _lE[lw]);
                            Mlltest[0]->Fill(TMath::Min((Zlep1 + Zlep2).Mag(), Mlltest[0]->GetBinCenter(Mlltest[0]->GetNbinsX())), scale[sam]*_weight);
                            if(_charges[lw] != _charges[lz1]){
                                Mlltest[1]->Fill(TMath::Min((Zlep1 + Wlep).Mag(), Mlltest[1]->GetBinCenter(Mlltest[1]->GetNbinsX())), scale[sam]*_weight);
                            } else{
                                Mlltest[1]->Fill(TMath::Min((Zlep2 + Wlep).Mag(), Mlltest[1]->GetBinCenter(Mlltest[1]->GetNbinsX())), scale[sam]*_weight);
                            }
                        } else if( fabs(_gen_lmompdg[match[lw]]) == 15){
                            WZ_split[2]->Fill(TMath::Min(transmass(Wlepton, METvec), WZ_split[2]->GetBinCenter(WZ_split[2]->GetNbinsX())), scale[sam]*_weight); 
                        } else if(_gen_isPromptl[match[lw]] && ( fabs(_gen_lmompdg[lz1]) == 23 || fabs(_gen_lmompdg[match[lz2]]) == 23)){
                            WZ_split[3]->Fill(TMath::Min(transmass(Wlepton, METvec), WZ_split[3]->GetBinCenter(WZ_split[3]->GetNbinsX())), scale[sam]*_weight);
                        } else if(_gen_isPromptl[match[lw]]){
                            WZ_split[4]->Fill(TMath::Min(transmass(Wlepton, METvec), WZ_split[4]->GetBinCenter(WZ_split[4]->GetNbinsX())), scale[sam]*_weight);
                        } else{
                            WZ_split[5]->Fill(TMath::Min(transmass(Wlepton, METvec), WZ_split[5]->GetBinCenter(WZ_split[5]->GetNbinsX())), scale[sam]*_weight);
                        }
						//Filling MET res splitted histograms
						if( fabs(_met - _genmet) < 20){
						//TLorentzVector genMETvec;
						//genMETvec.SetPtEtaPhiE(_genmet, 0, _genmet_phi, _genmet);
						//if( fabs(METvec.DeltaPhi(genMETvec)) < 0.2){
							WZ_METres[0]->Fill( TMath::Min(transmass(Wlepton, METvec), WZ_METres[0]->GetBinCenter(WZ_METres[0]->GetNbinsX())), scale[sam]*_weight);
						} else if( fabs(_met - _genmet) < 40){
						//} else if( fabs(METvec.DeltaPhi(genMETvec)) < 0.4){
							WZ_METres[1]->Fill( TMath::Min(transmass(Wlepton, METvec), WZ_METres[1]->GetBinCenter(WZ_METres[1]->GetNbinsX())), scale[sam]*_weight);
						} else if( fabs(_met - _genmet) < 80){
						//} else if( fabs(METvec.DeltaPhi(genMETvec)) < 0.6){
							WZ_METres[2]->Fill( TMath::Min(transmass(Wlepton, METvec), WZ_METres[2]->GetBinCenter(WZ_METres[2]->GetNbinsX())), scale[sam]*_weight);
						} else{
							WZ_METres[3]->Fill( TMath::Min(transmass(Wlepton, METvec), WZ_METres[3]->GetBinCenter(WZ_METres[3]->GetNbinsX())), scale[sam]*_weight);
						}

                    }
                    /*
                       if(sam == 0){
                       TLorentzVector Zlep1, Zlep2, Wlep;
                       Zlep1.SetPtEtaPhiE(_lPt[lz1], _lEta[lz1], _lPhi[lz1], _lE[lz1]);
                       Zlep2.SetPtEtaPhiE(_lPt[lz2], _lEta[lz2], _lPhi[lz2], _lE[lz2]);
                       Wlep.SetPtEtaPhiE(_lPt[lw], _lEta[lw], _lPhi[lw], _lE[lw]);
                       Mlltest[0]->Fill(TMath::Min((Zlep1 + Zlep2).Mag(), Mlltest[0]->GetBinCenter(Mlltest[0]->GetNbinsX())), scale[sam]*_weight);
                       if(_charges[lw] != _charges[lz1]){
                       Mlltest[1]->Fill(TMath::Min((Zlep1 + Wlep).Mag(), Mlltest[1]->GetBinCenter(Mlltest[1]->GetNbinsX())), scale[sam]*_weight);
                       } else{
                       Mlltest[1]->Fill(TMath::Min((Zlep2 + Wlep).Mag(), Mlltest[1]->GetBinCenter(Mlltest[1]->GetNbinsX())), scale[sam]*_weight);
                       }
                       }
                     */
                }
            }
            if((BosonPtReweighing || HTReweighing || DeltaRReweighing || leptonPtReweighing || RecoleptonPtReweighing || DeltaPhiReweighing) && r == 0 && (sam == 2 || (sam == 1 && SkipWjets)) ){
                for(int i = 1; i < 12; ++i){
                    ReweighingHistos[i -1] = HistDiv(Histos[0][i], Histos[1][i]);
                    ReweighingHistos[i-1]->Sumw2();
                }
                PTre = HistDiv(Histos[0][2], Histos[1][2]);
                PTre->Sumw2();
                Phire= HistDiv(Histos[0][9], Histos[1][9]);
                PTre_up = (TH1D*) PTre->Clone("upwards");
                PTre_down = (TH1D*) PTre->Clone("downwards");	
                //plotHist(PTre, "Reweighingfactors" + extra);
                for(int i = 0; i < PTre_up->GetNbinsX() + 1; ++i){
                    PTre_up->SetBinContent(i, PTre_up->GetBinContent(i) + PTre_up->GetBinError(i));
                }
                for(int i = 0; i < PTre_down->GetNbinsX() + 1; ++i){
                    PTre_down->SetBinContent(i, PTre_down->GetBinContent(i) - PTre_down->GetBinError(i));
                }	
                TFile *outputfile = new TFile ("reweighing_hists.root" , "recreate");
                /*
                   for(int i = 0; i < 11; ++i){
                   ReweighingHistos[i]->Write();
                   }
                 */
                PTre->Write();
                Histos[0][0]->Write();
                Histos[0][3]->Write();
                WZPU[0]->Write();
                WZPU[1]->Write();
                MTMETPU->Write();
                for(int i = 0; i < 12; ++i){
                    Histos[1][i]->Reset();
                }
                outputfile->Close();
            }
            dumpFile[sam].close();
        }
    }

    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "N(Wgamma)/N(WZ) = " << Pass[1]/(Pass[0]) << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "N(Wjets)/N(Wgamma) = " << Pass[2]/Pass[1] << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;

    //cout << "integral ratio Wgamma/WZ : " << (MT[1]->Integral())/(MT[0]->Integral()) << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    for(int i = 0; i < nSamples; ++i){
        cout <<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
        cout << "Expected " << names[i] << " events : " << Pass[i] << endl;
        cout << names[i] << " MC events passing selection : " << Total[i] << endl;
        cout <<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    }

    //Histos[1][0]->Scale( 1/(Histos[1][0]->Integral()/Histos[0][0]->Integral()) );
    //plotHistRatio(Histos[0][0], MT_PTg[1], "WZ", "W#gamma", "MTcomp70", false, 0,0,true);
    plotHistRatio(Histos[0][0], Histos[1][0], "WZ", "W#gamma", "MTcomp" + extra, false, MTCut,MTMax,true);
    plotHistRatio(Histos[0][2], Histos[1][2], "WZ", "W#gamma", "LeptonPtcomp" + extra);
    plotHistRatio(Histos[0][3], Histos[1][3], "WZ", "W#gamma", "METcomp" + extra);
    //plotHistRatio(Histos[0][9], Histos[1][9], "WZ", "W#gamma", "DeltaPhi_lepMETcomp" + extra);

	THStack* WZ_Stack;
    WZ_Stack = new THStack("WZ_Stack", "WZ_Stack");
    for(int i = 0; i < 6; ++i){
        StackCol(WZ_split[i], colors[i]);
        WZ_Stack->Add(WZ_split[i], "f");
    }
    TLegend *WZ_Legend = new TLegend(0.65,0.40,0.95,0.9,NULL,"brNDC");
    WZ_Legend->SetFillStyle(0);
    const TString source_name[6] = {"W lepton", "mispaired", "#tau", "prompt + Z", "other prompt", "fake"};
    for(int i = 0; i < 6; ++i){
        WZ_Legend->AddEntry(WZ_split[i], source_name[i]);
    }
    plotDataVSMC(Histos[0][0], Histos[0][0], WZ_Stack, WZ_Legend, "WZ_split" + extra, true);
	
	
    THStack* WZ_metStack;
    WZ_metStack = new THStack("WZ_metStack", "WZ_metStack");
    for(int i = 0; i < 4; ++i){
        StackCol(WZ_METres[i], colors[i]);
        WZ_metStack->Add(WZ_METres[i], "f");
    }
    TLegend *WZ_metLegend = new TLegend(0.65,0.40,0.95,0.9,NULL,"brNDC");
    WZ_metLegend->SetFillStyle(0);
    const TString metb_name[4] = {"#DeltaMET < 20 GeV", "#DeltaMET < 40 GeV", "#DeltaMET < 80 GeV", "#DeltaMET > 80 GeV"};
	//const TString metb_name[4] = {"#Delta#PhiMET < 0.2", "#Delta#PhiMET < 0.4", "#DeltaPhiMET < 0.6", "#Delta#PhiMET > 0.6"};
    for(int i = 0; i < 4; ++i){
        WZ_metLegend->AddEntry(WZ_METres[i], metb_name[i]);
    }
    plotDataVSMC(Histos[0][0], Histos[0][0], WZ_metStack, WZ_metLegend, "WZ_METsplit" + extra, true);
	//plotDataVSMC(Histos[0][0], Histos[0][0], WZ_metStack, WZ_metLegend, "WZ_METPhisplit" + extra, true);
	
	
	WZ_METres[0]->Add(WZ_METres[1]);
	WZ_METres[2]->Add(WZ_METres[3]);
	cout << "MT > 200:  DeltaMET < 40 / DeltaMET > 40  = " << WZ_METres[0]->Integral()/ WZ_METres[2]->Integral() << endl;
	
	//plotHistRatio(WZ_METres[0], WZ_METres[3], "#DeltaMET < 80", "#DeltaMET > 80", "METratio" + extra, false, MTCut, MTMax, true);
	
    plotHistRatio(Mlltest[0], Mlltest[1], "Z pairing", "W + Z leptons", "Mllsplot" + extra);
	

    //plotHistRatio(Histos[0][0], WZPU[0], "WZ", "WZ PU rew", "MTPUcomp" + extra, true, MTCut,MTMax,true);	
    //plotHistRatio(Histos[0][3], WZPU[1], "WZ", "WZ PU rew", "METPUcomp" + extra);
    /*
       cout << "\\begin{tabular}{|c|c|c|c|}" << endl;
       cout << "\\hline" << endl;;
       cout << "& $M_{T} < 120$ & $120 < M_{T} < 160$ & $M_{T} >160$\\\\" << endl;	
       for(int i =0; i <3; ++i){
       MT_PTg[i]->Scale( 1/(MT_PTg[i]->Integral()/Histos[0][0]->Integral()) );
    //plotHistRatio(MT[0][0], MT[1][i],lepton +" + #gamma data", "WZ MC", "MTshape_dataMCcomp" + extra, false,MTMin,MTMax,true,true, true);

    //Calculate scale factors 
    double bins[4] = {40,120,160,200};
    TH1D* WZ_MT_reb =(TH1D*) Histos[0][0]->Rebin(3, "WZ_MT_reb", bins);
    TH1D* Wgamma_MT_reb =(TH1D*) MT_PTg[i]->Rebin(3, "Wgamma_MT_reb", bins);
    TH1D *ScalFacs = HistDiv(Wgamma_MT_reb, WZ_MT_reb);
    //plotHist(ScalFacs, "scaletest");
    cout << PTgbins[i] << " & ";
    for(int j = 1; j < ScalFacs->GetNbinsX() + 1; ++j){
    cout <<  std::setprecision(3) << ScalFacs->GetBinContent(j) << " $\\pm$ " <<   std::setprecision(3) << ScalFacs->GetBinError(j) << " & ";
    }
    cout << endl;
    cout << "\\hline" << endl;
    }
    cout << "\\end{tabular}" << endl;

    std::vector<TH1D*> histos;
    std::vector<TString> entries;

    histos.push_back(Histos[0][0]);
    for(int i = 0; i < 3; ++i) histos.push_back(MT_PTg[i]);
    entries.push_back("WZ");
    entries.push_back("W#gamma, P_{T} > 50 GeV");
    entries.push_back("W#gamma, P_{T} > 70 GeV");
    entries.push_back("W#gamma, P_{T} > 90 GeV");
    plotHist(histos,entries, "PTsyst" + extra, true);
     */
    /*
       cout << "\\begin{tabular}{|c|c|c|c|}" << endl;
       cout << "\\hline" << endl;
       cout << "& $M_{T} < 120$ & $120 < M_{T} < 160$ & $M_{T} >160$\\\\" << endl;
       for(int i =0; i <3; ++i){
       MT_Scalfluc[i]->Scale( 1/(MT_Scalfluc[i]->Integral()/Histos[0][0]->Integral()) );
    //plotHistRatio(MT[0][0], MT[1][i],lepton +" + #gamma data", "WZ MC", "MTshape_dataMCcomp" + extra, false,MTMin,MTMax,true,true, true);

    //Calculate scale factors 
    double bins[4] = {40,120,160,200};
    TH1D* WZ_MT_reb =(TH1D*) Histos[0][0]->Rebin(3, "WZ_MT_reb", bins);
    TH1D* Wgamma_MT_reb =(TH1D*) MT_Scalfluc[i]->Rebin(3, "Wgamma_MT_reb", bins);
    TH1D *ScalFacs = HistDiv(Wgamma_MT_reb, WZ_MT_reb);
    //plotHist(ScalFacs, "scaletest");
    cout <<Scalfluc[i] << " & ";
    for(int j = 1; j < ScalFacs->GetNbinsX() + 1; ++j){
    cout <<  std::setprecision(3) << ScalFacs->GetBinContent(j) << " $\\pm$ " <<   std::setprecision(3) << ScalFacs->GetBinError(j) << " & ";
    }
    cout << endl;
    cout << "\\hline" << endl;
    }
    cout << "\\end{tabular}" << endl;

    std::vector<TH1D*> histos;
    std::vector<TString> entries;
    histos.clear();
    entries.clear();
    histos.push_back(Histos[0][0]);
    for(int i = 0; i < 3; ++i) histos.push_back(MT_Scalfluc[i]);
    entries.push_back("WZ");
    entries.push_back("W#gamma, + 0 #sigma");
    entries.push_back("W#gamma, + 1 #sigma");
    entries.push_back("W#gamma, - 1 #sigma");
    plotHist(histos,entries, "Scalsyst" + extra, true);
     */

    TFile *MC_SF = new TFile("MC_SF.root", "recreate");
    MTMET[1]->Scale( 1/(MTMET[1]->Integral()/MTMET[0]->Integral()) );
    TH2D* MTMETC = (TH2D*) MTMET[0]->Clone("BinnedSF"); 
    const TString MTb[3] = {"$M_{T} < 120$",  "$120 < M_{T} < 160$", "$M_{T} > 160$"};
    cout << "\\textbf{MC binned Scale Factors :}" << endl;
    cout << "\\begin{tabular}{|c|c|c|c|c|}" << endl;
    cout << "\\hline" << endl;
    cout << "& 50 $<$ MET $<$ 100 & 100 $<$ MET $<$ 150 & 150 $<$ MET $<$ 200 & MET $>$ 200\\\\" << endl;
    cout << "\\hline" << endl;
    MTMETC->Divide(MTMET[1]);
    for(int i = 1; i < MTMETC->GetNbinsX() + 1; ++i){
        cout << MTb[i-1] << " & ";
        for(int j = 1; j < MTMETC->GetNbinsY() + 1; ++j){
            cout << std::setprecision(3) << MTMETC->GetBinContent(i,j) << " $\\pm$ " << std::setprecision(3) << MTMETC->GetBinError(i,j);
            if(j !=  MTMETC->GetNbinsY()) cout << " & ";
        }
        cout << "\\\\ \\hline" << endl;
    }
    cout << "\\end{tabular}" << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    MTMETC->Write();

    for( int r = 0; r < 2; ++r){
        for( int i = 0; i < MTMET[r]->GetNbinsX() + 1; ++i){
            for(int j = 0; j < MTMET[r]->GetNbinsY() + 1; ++j){
                double BinC = 0;
                double BinE2 = 0;
                for(int k = MTMET[r]->GetNbinsY(); k >= j; --k){
                    BinC += MTMET[r]->GetBinContent(i,k);
                    BinE2 += (MTMET[r]->GetBinError(i,k))*(MTMET[r]->GetBinError(i,k));
                }
                MTMET[r]->SetBinContent(i, j, BinC);
                MTMET[r]->SetBinError(i, j, sqrt(BinE2));
            }
        }
    }

    MTMET[0]->Divide(MTMET[1]);
    cout << "\\textbf{MC Inclusive Scale Factors :}" << endl;
    cout << "\\begin{tabular}{|c|c|c|c|c|}" << endl;
    cout << "\\hline" << endl;
    cout << "& MET $>$ 50 & MET $>$ 100 & MET $>$ 150 & MET $>$ 200\\\\" << endl;	
    cout << "\\hline" << endl;
    for(int i = 1; i < MTMET[0]->GetNbinsX() + 1; ++i){
        cout << MTb[i-1] << " & ";
        for(int j = 1; j < MTMET[0]->GetNbinsY() + 1; ++j){
            cout << std::setprecision(3) << MTMET[0]->GetBinContent(i,j) << " $\\pm$ " << std::setprecision(3) << MTMET[0]->GetBinError(i,j);
            if(j !=  MTMET[0]->GetNbinsY()) cout << " & ";
        }
        cout << "\\\\ \\hline" << endl;
    }
    cout << "\\end{tabular}" << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    MTMET[0]->Write("InclSF");	
    MC_SF->Close();

    rootapp->Run();
    return 0;	

};

