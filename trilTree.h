#ifndef trilTree_h
#define trilTree_h

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TLorentzVector.h"

#include <TF1.h>
#include <TH1.h>

#include "../bTag/BTagCalibrationStandalone.h"
//#include "TClonesArray.h"

class trilTree {
public :
   	TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   	Int_t           fCurrent; //!current Tree number in a TChain

    //Declare leaf types
	static const int nL_max = 10;
	static const int nPh_max = 10;
	static const int nJets_max = 10;
	Double_t _weight;
	//Reco level leaf types for event selection
	Int_t _nL;
	Int_t _nEle;
	Int_t _nMu;
	Double_t _lPt[nL_max];
	Double_t _lE[nL_max];
	Double_t _lEta[nL_max];
	Double_t _lPhi[nL_max];
	Int_t _flavors[nL_max];
	Double_t _charges[nL_max];
	Bool_t _isloose[nL_max];
	Bool_t _isFO[nL_max];
	Bool_t _istight[nL_max];
	Bool_t _islooseCut[nL_max];
	Bool_t _isFOCut[nL_max];
	Bool_t _istightCut[nL_max];
	Bool_t _lClean[nL_max];
	Bool_t _lCleanCut[nL_max];
	Bool_t _tightCharge[nL_max];
	Bool_t _mu_ismedium[nL_max];
	Int_t _hitsNumber[nL_max];
	Double_t _met;
	Double_t _met_phi;
	Int_t 	_nJets;
	Double_t _HT;
	Double_t _jetEta[nJets_max];
	Double_t _jetPhi[nJets_max];
	Double_t _jetPt[nJets_max];
	Double_t _jetE[nJets_max];
	Double_t _csv[nJets_max];
	Int_t  _jetFlavour[nJets_max];
	Double_t _closeJetCSVAll[nL_max];
	Double_t _closeJetPtAll[nL_max];
	Double_t _isolation[nL_max];
	Int_t _n_bJets;
	ULong64_t _runNb;
	ULong64_t _lumiBlock;
	ULong64_t _eventNb;
	Double_t _lepMVA[nL_max];
	Double_t _mvaValue[nL_max];
	Bool_t _clusterpass[nL_max];
	Bool_t _vtxFitConversion[nL_max];
	Bool_t _METfilters_pass;
	Int_t _n_trueInteractions;
	Int_t _n_PV;
	Bool_t _tril_trigger_all;
	Bool_t _tril_trigger_eee;
	Bool_t _tril_trigger_mmm;
	Bool_t _tril_trigger_eem;
	Bool_t _tril_trigger_emm;
	Bool_t _tril_trigger_emt; 
	Bool_t _tril_trigger_ett;
	Bool_t _tril_trigger_mtt;
	Bool_t _lowM_trigger_all;
    Bool_t _lowM_trigger_mmm;
    Bool_t _lowM_trigger_eee;
    Bool_t _lowM_trigger_mme;
    Bool_t _lowM_trigger_mee;
	Bool_t _Mu17_Photon22_CaloIdL_L1ISO;
	Bool_t _IsoMu24;
	Bool_t _IsoTkMu24;
	Double_t _ptrel[nL_max];
	Double_t _miniisolation[nL_max][2];
	Double_t _ptratio[nL_max];
	Double_t _3dIP[nL_max];
	Double_t _ipPV[nL_max];
	Double_t _ipZPV[nL_max];
	Double_t _3dIPsig[nL_max];
	//Electron trigger emulation
	/*
	Double_t _dEtaIn[nL_max];
    Double_t _dPhiIn[nL_max];
    Double_t _ele_hOverE[nL_max];
    Double_t _full5x5SigmaIetaIeta[nL_max];
    Double_t _ooEmooP[nL_max];
	*/
	//
	//Gen leptons
	Int_t _gen_nL;
	Double_t _gen_lPt[nL_max];
	Double_t _gen_lEta[nL_max];
	Double_t _gen_lPhi[nL_max];
	Double_t _gen_lE[nL_max];
	Bool_t _gen_isPromptl[nL_max];
	Int_t _gen_lmompdg[nL_max];
	Int_t _gen_flavors[nL_max];
	Double_t _gen_charges[nL_max];
	Int_t _pdgmc[nL_max];
	Bool_t _isPromptFinalState[nL_max];
	Bool_t _fromHardProcessFinalState[nL_max];
	Int_t _origin[nL_max];
	Int_t _originDetailed[nL_max];
	Int_t _originPhot[nL_max];
	Int_t _mompdg[nL_max];
	//MET uncertainties
	Double_t _metJECDown;
	Double_t _metJECUp;
	Double_t _metOtherUp;
	Double_t _metOtherDown;
	Double_t _met_phiJECDown;
	Double_t _met_phiJECUp;
	Double_t _met_phiOtherDown;
	Double_t _met_phiOtherUp;
	//JEC uncertainties	
	Double_t _jetPtDown[nJets_max];
	Double_t _jetPtUp[nJets_max];
	//pdf and scale uncertainty weights
	Double_t _scaleWeight[111];
	//Photons
	Int_t _nPh;
	Double_t _phPt[nPh_max];
    Double_t _phEta[nPh_max];
    Double_t _phPhi[nPh_max];
    Double_t _phE[nPh_max];
    Double_t _phmva[nPh_max];
    Bool_t _phmva_pass[nPh_max];
    Bool_t _pheveto_pass[nPh_max];
    Bool_t _phpixelseed_pass[nPh_max];
	//Generator photons
    Int_t  _gen_nPh;
	Double_t _gen_phPt[nPh_max];
	Double_t _gen_phEta[nPh_max];
	Double_t _gen_phPhi[nPh_max];
	Double_t _gen_phE[nPh_max];
    Bool_t _gen_isPromptPh[nPh_max];
	Int_t _gen_phmompdg[nPh_max];
	//Generator met
	Double_t _genmet;
	Double_t _genmet_phi;
	//Lifetime
	Double_t _ctau;
	

   	// List of branches
	TBranch *b__nL;
	TBranch *b__lPt;
	TBranch *b__lEta;
	TBranch *b__lPhi;
	TBranch *b__lE;
	TBranch *b__charges;
	TBranch *b__flavors;
	TBranch *b__isloose;
	TBranch *b__isFO;
	TBranch *b__istight;
	TBranch *b__islooseCut;
	TBranch *b__isFOCut;
	TBranch *b__lClean;
	TBranch *b__lCleanCut;
	TBranch *b__istightCut;
	TBranch *b__tightCharge;
	TBranch *b__mu_ismedium;
	TBranch *b__hitsNumber;
	TBranch *b__met;
	TBranch *b__met_phi;
	TBranch *b__nJets;
	TBranch *b__jetEta;
	TBranch *b__jetPhi;
	TBranch *b__jetPt;
	TBranch *b__jetE;
	TBranch *b__csv;
	TBranch *b__jetFlavour;
	TBranch *b__closeJetCSVAll;
	TBranch *b__closeJetPtAll;
	TBranch *b__weight;
	TBranch *b__isolation;
	TBranch *b__n_bJets;
	TBranch *b__eventNb;
	TBranch *b__runNb;
	TBranch *b__lumiBlock;
	TBranch *b__lepMVA;
	TBranch *b__mvaValue;
	TBranch *b__clusterpass;
	TBranch *b__vtxFitConversion;
	TBranch *b__METfilters_pass;
	TBranch *b__n_trueInteractions;
	TBranch *b__n_PV;
	TBranch *b__tril_trigger_all;
	TBranch *b__tril_trigger_eee;
	TBranch *b__tril_trigger_mmm;
	TBranch *b__tril_trigger_eem;
	TBranch *b__tril_trigger_emm;
	TBranch *b__tril_trigger_emt;
	TBranch *b__tril_trigger_ett;
	TBranch *b__tril_trigger_mtt;
	TBranch *b__lowM_trigger_all;
    TBranch *b__lowM_trigger_mmm;
    TBranch *b__lowM_trigger_eee;
    TBranch *b__lowM_trigger_mme;
    TBranch *b__lowM_trigger_mee;
	TBranch *b__Mu17_Photon22_CaloIdL_L1ISO;
	TBranch *b__IsoMu24;
	TBranch *b__IsoTkMu24;
	TBranch *b__ptrel;
	TBranch *b__miniisolation;
	TBranch *b__ptratio;
	TBranch *b__ipPV;
	TBranch *b__ipZPV;
	TBranch *b__3dIPsig;
	TBranch *b__3dIP;
	//Electron trigger emulation
	/*
	TBranch *b__dEtaIn;
    TBranch *b__dPhiIn;
    TBranch *b__ele_hOverE;
    TBranch *b__full5x5SigmaIetaIeta;
    TBranch *b__ooEmooP;
	*/
	//Generator lepton branches
	TBranch *b__gen_nL;
	TBranch *b__gen_lPt;
	TBranch *b__gen_lEta;
	TBranch *b__gen_lPhi;
	TBranch *b__gen_lE;
	TBranch *b__gen_isPromptl;
	TBranch *b__gen_lmompdg;
	TBranch *b__gen_flavors;
	TBranch *b__gen_charges;

	TBranch* b__pdgmc;
	TBranch* b__isPromptFinalState;
	TBranch* b__fromHardProcessFinalState;
	TBranch* b__origin;
	TBranch* b__originDetailed;
	TBranch* b__originPhot;
	TBranch* b__mompdg;
	//MET uncertainty branches
	TBranch* b__metJECDown;
	TBranch* b__metJECUp;
	TBranch* b__metOtherDown;
	TBranch* b__metOtherUp;
	TBranch* b__met_phiJECDown;
	TBranch* b__met_phiJECUp;
	TBranch* b__met_phiOtherDown;
	TBranch* b__met_phiOtherUp;
	//JEC uncertainty branches
	TBranch* b__jetPtDown;
	TBranch* b__jetPtUp;
	//pdf and scale uncertainty weight branches
	TBranch* b__scaleWeight;
	//Photon branches
	TBranch *b__nPh;
	TBranch *b__phPt; 
    TBranch *b__phEta;
    TBranch *b__phPhi;
    TBranch *b__phE;
    TBranch *b__phmva;
    TBranch *b__phmva_pass;
    TBranch *b__pheveto_pass;
    TBranch *b__phpixelseed_pass;
	//Generator photon branches
    TBranch *b__gen_nPh;
	TBranch *b__gen_phPt;
	TBranch *b__gen_phEta;
	TBranch *b__gen_phPhi;
	TBranch *b__gen_phE;
    TBranch *b__gen_isPromptPh;
	TBranch *b__gen_phmompdg;
	//Generator met branches
	TBranch *b__genmet;
	TBranch *b__genmet_phi;
	//lifetime branch
	TBranch *b__ctau;

   	trilTree(TTree *tree=0);
   	Int_t GetEntry(Long64_t entry);
   	//Long64_t LoadTree(Long64_t entry);
   	void Init(TTree *tree, const bool photonTree = false, const bool genTree = false);
   	void Loop();
	
	//New test functions
	void Loop(const TString& fileName, const double xSection, const bool isData, const bool isSignal);
	void hnlPlots();

	//Define Histrograms containing analysis scale factors and maps
	void readSF(const bool hnl = false);
	double getEventSF(const unsigned* ind, const unsigned lCount, const bool hnl = false);
	double bTagSF(const bool hnl = false, const unsigned bTagSyst = 0);
	//Read FR maps
	TH2D* frMap[3][3]; //First dimension concerns EWK contamination variation uncertainties, second the lepton flavor
	//Read SF maps
	TH2D* baseidSFMap[2];
	TH2D* idTightSFMap[2];
	TH1D* recSFMap_ele;
 	TGraph* recSFMap_mu_ptAbove10;
	TGraph* recSFMap_mu_ptBelow10;
	BTagCalibration calib;
	BTagCalibrationReader reader;
	TH2D* bTagEff[3]; // 0: udsg, 1: C, 2: b
	//Trigger SF
	/*
	TH2D* triggerEffMap3l[2];
	TH2D* triggerEffMap2l[2];
	TH2D* triggerEffEleLeg;
	*/
	//PU reweighing
	TH1D* PUweights[3]; //central, up, down
	//Functions for event selection
	bool baseline(const bool multilepton = true, const bool looseBtag = false, const bool jetCleaning = true, const bool lowMllVeto = true);
	bool vetoLowMll(const double mllCut = 12.);
	bool jetIsClean(unsigned jetI);
	void cleanJets();
	unsigned nCleanedJets();
	void selectJets();
	unsigned nBJets(const bool looseBtag = false, const bool jetCleaning = true, const unsigned jec = 0);
	unsigned lepOrder(unsigned*, const unsigned lepCut = 0, const bool skipTau = false, const bool cutBased = false);
	unsigned tightCount(const unsigned*, const unsigned);
	bool ptCuts(unsigned*, unsigned&);
	bool ptCuts_hnl(const unsigned*, const unsigned);
	bool promptMC(const unsigned*, const unsigned);
	bool trigger(const unsigned*, const unsigned);
	bool wzSel(const unsigned*, unsigned*, unsigned&,  const unsigned, const TLorentzVector*, const bool flavorDiff = true, const bool onZ = false, const bool newalgo = false);
	bool wgSel(unsigned&, const unsigned*, const unsigned, const double, const double, const double, const bool pixVeto = true, const bool eVeto = true, const bool onlyMu = true);
	//Match reco particles to generator particles
	void matchGenLep(unsigned*);
	void cutBased();
	//lifetime correction
	std::pair<double, double> findLifeTime(const double left, const double right, const std::string& sample, const bool reweighing = true);	
	void fitcTau(const std::string& sample);
	void comparecTau(const std::vector<std::string>&);
	void effVScTau(const std::string& sample);
};
#endif
