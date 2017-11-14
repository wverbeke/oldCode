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
#include "TROOT.h"
#include "TGraph.h"

//include C++ library classes
#include <sstream>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using std::cout;
using std::endl;
using std::flush;
using std::ofstream;

//include code to calculate btag SF
#include "../bTag/BTagCalibrationStandalone.h"
//include other parts of the code
#include "MultilepSUSYfunc.h"
#include "tdrstyle.h"
#include "plotCode.h"
#include "trilTree.h"
#include "hnlTools.h"


void trilTree::Loop(){
	//Set plotting style
	setTDRStyle();
	gROOT->SetBatch(kTRUE);
	const unsigned nSamples = 65;
	const unsigned nSamples_eff = 23;
	const unsigned nSig = 17;
	const double lowMCoupling = 0.00001;
	const double highMCoupling = 0.01;
	const double corrFactor = 200.;
	//~~~~~~~~ background normalizations~~~~~~~~~~~~~~~~~~~~
	const double ZZSF = 1.3043*1.00788;
	const double WZSF = 0.650269*0.983686;
	const double convSF = 0.897153*0.988728;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//TTGJets.root
	/*
	const TString fileList[nSamples] = {"MuonEG.root", "DoubleMuon.root", "DoubleEG.root", "SingleMuon.root","SingleElectron.root", "trilMaj5.root", "trilMaj15.root", "trilMaj20.root", "trilMaj30.root", "trilMaj40.root", "trilMaj50.root", "trilMaj60.root", "trilMaj80.root", "trilMaj100.root", "trilMaj130.root", "trilMaj150.root", "trilMaj200.root", "trilMaj300.root", "trilMaj400.root", "ZZTo4L.root",  "VHToNonbb.root", "GluGluHToZZTo4L_M125_trilepton.root", "VBF_HToZZTo4L_M125_trilepton.root", "WWG.root","WWW.root", "WWZ.root", "WWTo2L2Nu_DoubleScattering.root", "WWTo2L2Nu.root",  "ZZZ.root", "WZTo3LNu_mllmin01.root", "TTGJets.root","ZGTo2LG.root", "WGToLNuG.root", "TGJets.root", "ST_tW_antitop_NofullyHadronic.root", "ST_tW_top_NofullyHadronic.root", "ST_s-channel_leptonDecays.root", "ST_t-channel_top_inclusiveDecays.root", "ST_t-channel_antitop_inclusiveDecays.root", "ttHToNonbb.root", "TTWJetsToLNu.root", "TTZToLLNuNu.root",  "TTTT.root"};
	//const double xSections[nSamples - 5] = {0.215, 0.2043, 0.2529, 3.697, 1.256, 0.2147,  0.009103, 0.9561, 0.2086, 0.1651, 0.1729,   12.178, 0.01398, 4.4297, 123.9, 87.315, 182.175, 182.175, 61526.7, 350.674, 2.967, 38.09, 38.09, 10.11, 136.02, 80.95}; 
	//3.697
	*/
		const TString fileList[nSamples] = {"data_combined_trilepton.root", 
"MajoranaNeutrinoToSSSF_MuMuE_M5.root", "MajoranaNeutrinoToMuMuMu_M5.root",
"MajoranaNeutrinoToSSSF_MuMuE_M10.root", "MajoranaNeutrinoToMuMuMu_M10.root",
"MajoranaNeutrinoToSSSF_MuMuE_M20.root", "MajoranaNeutrinoToMuMuMu_M20.root",
"MajoranaNeutrinoToSSSF_MuMuE_M30.root", "MajoranaNeutrinoToMuMuMu_M30.root",
"MajoranaNeutrinoToSSSF_MuMuE_M40.root", "MajoranaNeutrinoToMuMuMu_M40.root",
"MajoranaNeutrinoToSSSF_MuMuE_M50.root", "MajoranaNeutrinoToMuMuMu_M50.root",
"MajoranaNeutrinoToSSSF_MuMuE_M60.root", "MajoranaNeutrinoToMuMuMu_M60.root",
"MajoranaNeutrinoToSSSF_MuMuE_M70.root", "MajoranaNeutrinoToMuMuMu_M70.root",
"MajoranaNeutrinoToSSSF_MuMuE_M90.root", "MajoranaNeutrinoToMuMuMu_M90.root",
"MajoranaNeutrinoToSSSF_MuMuE_M100.root", "MajoranaNeutrinoToMuMuMu_M100.root",
"MajoranaNeutrinoToSSSF_MuMuE_M150.root", "MajoranaNeutrinoToMuMuMu_M150.root",
"MajoranaNeutrinoToSSSF_MuMuE_M200.root", "MajoranaNeutrinoToMuMuMu_M200.root",
"MajoranaNeutrinoToSSSF_MuMuE_M300.root", "MajoranaNeutrinoToMuMuMu_M300.root",
"MajoranaNeutrinoToSSSF_MuMuE_M400.root", "MajoranaNeutrinoToMuMuMu_M400.root",
"MajoranaNeutrinoToSSSF_MuMuE_M500.root", "MajoranaNeutrinoToMuMuMu_M500.root",
"MajoranaNeutrinoToSSSF_MuMuE_M700.root", "MajoranaNeutrinoToMuMuMu_M700.root",
"MajoranaNeutrinoToSSSF_MuMuE_M1000.root", "MajoranaNeutrinoToMuMuMu_M1000.root",

"ZZTo4L_trilepton.root",  "VHToNonbb_trilepton.root", "GluGluHToZZTo4L_M125_trilepton.root", "VBF_HToZZTo4L_M125_trilepton.root", "WWG_trilepton.root","WWW_trilepton.root", "WWZ_trilepton.root", "WWTo2L2Nu_DoubleScattering_trilepton.root", "WWTo2L2Nu_trilepton.root",  "ZZZ_trilepton.root", "WZTo3LNu_mllmin01_trilepton.root", "TTGJets_trilepton.root","ZGTo2LG_trilepton.root", "WGToLNuG_trilepton.root", "TGJets_trilepton.root", "TTJets_DiLept_trilepton.root", "TTJets_SingleLeptFromTbar_trilepton.root", "TTJets_SingleLeptFromT_trilepton.root", "DYJetsToLL_M10to50_trilepton.root", "DYJetsToLL_M50_trilepton.root", "ST_tW_antitop_NofullyHadronic_trilepton.root", "ST_tW_top_NofullyHadronic_trilepton.root", "ST_s-channel_leptonDecays_trilepton.root", "ST_t-channel_top_inclusiveDecays_trilepton.root", "ST_t-channel_antitop_inclusiveDecays_trilepton.root", "ttHToNonbb_trilepton.root", "TTWJetsToLNu_trilepton.root", "TTZToLLNuNu_trilepton.root", "TTZToLL_M1to10_trilepton.root",  "TTTT_trilepton.root"};
	
	// 5.072e+03*lowMCoupling  = 15 GeV xsection
	const double xSections[nSamples - 1] = {
3.502*corrFactor*lowMCoupling, 4.046*corrFactor*lowMCoupling, //5
3.422*corrFactor*lowMCoupling, 3.982*corrFactor*lowMCoupling, //10
3.169*corrFactor*lowMCoupling, 3.711*corrFactor*lowMCoupling, //20
2.751*corrFactor*lowMCoupling, 3.241*corrFactor*lowMCoupling, //30
2.185*corrFactor*lowMCoupling, 2.619*corrFactor*lowMCoupling, //40
1.548*corrFactor*lowMCoupling, 1.875*corrFactor*lowMCoupling, //50
0.8685*corrFactor*lowMCoupling, 1.07*corrFactor*lowMCoupling, //60
0.2882*corrFactor*lowMCoupling*10, 0.3828*corrFactor*lowMCoupling*10, //70
0.01166*corrFactor*lowMCoupling*100, 0.02333*corrFactor*lowMCoupling*100, //90
0.005443*corrFactor*highMCoupling, 0.01082*corrFactor*highMCoupling, //100
0.0007544*corrFactor*highMCoupling, 0.001488*corrFactor*highMCoupling, //150
0.0002293*corrFactor*highMCoupling, 0.0004567*corrFactor*highMCoupling, //200
4.77e-05*corrFactor*highMCoupling, 9.52e-05*corrFactor*highMCoupling, //300
1.58e-05*corrFactor*highMCoupling, 3.14e-05*corrFactor*highMCoupling, //400
6.45e-06*corrFactor*highMCoupling, 1.29e-05*corrFactor*highMCoupling, //500
1.61e-06*corrFactor*highMCoupling, 3.21e-06*corrFactor*highMCoupling, //700
3.28e-07*corrFactor*highMCoupling, 6.48e-07*corrFactor*highMCoupling, //1000

 1.256*ZZSF,  0.9561, 0.01212, 0.001034,  0.2147, 0.2086, 0.1651, 0.1729,   12.178, 0.01398, 58.59*WZSF, 3.697, 123.9*convSF, 489*convSF,  2.967, 87.315, 182.175, 182.175, 18610, 1921.8*3, 38.09, 38.09, 10.11, 136.02, 80.95, 0.215, 0.2043, 0.2529, 0.0493, 0.009103}; 

	const TString names[nSamples] = {"data", 
"Majorana5", "Majorana5",
"Majorana10", "Majorana10",
"Majorana20", "Majorana20",
"Majorana30", "Majorana30",
"Majorana40", "Majorana40",
"Majorana50", "Majorana50",
"Majorana60", "Majorana60",
"Majorana70", "Majorana70",
"Majorana90", "Majorana90",
"Majorana100", "Majorana100",
"Majorana150", "Majorana150",
"Majorana200", "Majorana200",
"Majorana300", "Majorana300",
"Majorana400", "Majorana400",
"Majorana500", "Majorana500",
"Majorana700", "Majorana700",
"Majorana1000", "Majorana1000",

"ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "triboson", "triboson", "triboson", "triboson", "triboson", "triboson", "WZ", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X"};
	//Read Trees from ROOT files
	TFile* hfile[nSamples];
	TTree* inputTree[nSamples];
	double hcounter[nSamples];
	for(unsigned sam = 0; sam < nSamples; ++sam){
		cout << "name " << names[sam] << endl;
		hfile[sam] = new TFile("../data_reminiaod17/"+fileList[sam],"read");
       	hfile[sam]->cd("FakeElectrons");
		//Determine hcounter for cross section scaling
		TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
		_hCounter->Read("hCounter");
		hcounter[sam] = _hCounter->GetBinContent(1);
       	inputTree[sam] = (TTree*) (hfile[sam]->Get("FakeElectrons/fakeTree"));
		Init(inputTree[sam], false, sam > 0);
	}

	readSF(true);
	/*
	BTagCalibration calib("csvv2", "../bTag/CSVv2_Moriond17_B_H.csv");
	BTagCalibrationReader reader(BTagEntry::OP_LOOSE, "central", {"up", "down"});
	reader.load(calib, BTagEntry::FLAV_B, "comb");
	*/
	//Tweakable options////////////////////////////////////////////////////
	const bool TestRun = false;	//Break after a few events
	const double DataLuminosity = 35.9;    //units of fb^{-1}
	const TString extra = "_hnlKin";	//for plot file names
	//////////////////////////

	const TString eff_names[nSamples_eff + 1] = {"data", 
"Majorana5",
"Majorana10",
"Majorana20",
"Majorana30",
"Majorana40",
"Majorana50",
"Majorana60",
"Majorana70",
"Majorana90",
"Majorana100",
"Majorana150",
"Majorana200",
"Majorana300",
"Majorana400",
"Majorana500",
"Majorana700",
"Majorana1000",
"ZZ/H", "triboson", "WZ", "X + #gamma", "TT/T + X",  "non-prompt"};
	
	const unsigned nCat = 6;  //Number of categories
	const TString catNames[nCat] = {"lowM_3lOSSF_lowPt", "lowM_3lnoOSSF_lowPt", "lowM_3lOSSF_highPt", "lowM_3lnoOSSF_highPt", "highM_3lOSSF", "highM_3lnoOSSF"};
	
	//Define 1D histograms
	const unsigned nDist = 59;  //Number of distributions to plot
	TH1D* Histos[nDist][nCat][nSamples_eff + 1];
	const TString Histnames[nDist] = {"Mll", "M3l", "mt3l", "mnu", "minMos", "mt_minMos", "MosminDeltaR", "mt_minDeltaR", "mt", "mt2_ss", "mt2_maxPt", "LeptonPt_le","LeptonPt_sub", "LeptonPt_tr", "MET", "HT", "NJets", "NbJets", "DeltaPhi_lepMET_le", "DeltaPhi_lepMET_sub", "DeltaPhi_lepMET_tr", "Nvtx","Pt_trilep", "Pt_trilepscalarsum", "ConePt_le", "ConePt_sub", "ConePt_tr", "Eta_le", "Eta_sub", "Eta_tr", "MiniIso_le", "MiniIso_sub", "MiniIso_tr", "RelIso_le", "RelIso_sub", "RelIso_tr", "Ptrel_le", "Ptrel_sub", "Ptrel_tr", "Ptratio_le", "Ptratio_sub", "Ptratio_tr", "csv_le", "csv_sub", "csv_tr", "3dIP_le", "3dIP_sub", "3dIP_tr", "dxy_le", "dxy_sub", "dxy_tr", "leadPtOverM3l", "subPtOverM3l", "trailPtOverM3l", "leadPtOverM4", "subPtOverM4", "trailPtOverM4", "M4", "M2lmaxPtnu"};

	const TString Xaxes[nDist] = {"M_{ll}(GeV)", "M_{3l} (GeV)", "M_{T}(3l) (GeV)", "M_{2l + #nu}", "min(M_{OS}) (GeV)", "M_{T}(other min(M_{OS})) (GeV)", "M_{OS}(min #Delta R) (GeV)", "M_{T}(other min #Delta R) (GeV)",  "M_{T} (GeV)", "M_{T2}(SS) (GeV)",  "M_{T2}(max P_{T} 2l) (GeV)", "P_{T}(leading l) (GeV)", "P_{T}(subleading l) (GeV)", "P_{T}(trailing l) (GeV)",  "MET (GeV)", "HT (GeV)", "number of jets", "number of b-jets", "#Delta#Phi(leading, MET)", "#Delta#Phi(subleading, MET)", "#Delta#Phi(trailing, MET)", "Number of vertices", "P_{T}(3l) (GeV)", "#Sigma_{l}(|P_{T}(l)|) (GeV)", "P_{T}^{cone}(leading) (GeV)", "P_{T}^{cone}(subleading) (GeV)", "P_{T}^{cone}(trailing) (GeV)", "|#eta(leading)|","|#eta(subleading)|", "|#eta(trailing)|", "miniIso(leading)", "miniIso(subleading)", "miniIso(trailing)", "relIso(leading)", "relIso(subleading)", "relIso(trailing)", "P_{T}^{rel}(leading) (GeV)", "P_{T}^{rel}(subleading) (GeV)", "P_{T}^{rel}(trailing) (GeV)", "P_{T}^{ratio}(leading)", "P_{T}^{ratio}(subleading)", "P_{T}^{ratio}(trailing)", "closest Jet CSV(leading)", "closest Jet CSV(subleading)", "closest Jet CSV(trailing)", "|3DIP(leading)| (cm)", "|3DIP(subleading)| (cm)", "|3DIP(trailing)| (cm)", "|d_{xy}(leading)| (cm)", "|d_{xy}(subleading)| (cm)", "|d_{xy}(trailing)| (cm)", "P_{T}^{cone}(leading)/M_{3l}", "P_{T}^{cone}(subleading)/M_{3l}", "P_{T}^{cone}(tailing)/M_{3l}", "P_{T}^{cone}(leading)/M_{3l + #nu}", "P_{T}^{cone}(subleading)/M_{3l + #nu}", "P_{T}^{cone}(tailing)/M_{3l + #nu}", "M_{3l + #nu} (GeV)", "M_{max P_{T} 2l  + #nu} (GeV)"};

	const TString Units[nDist] = {"GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "", "", "", "", "", "",  "GeV", "GeV", "GeV", "GeV", "GeV", "", "", "", "", "", "", "", "", "",  "GeV", "GeV", "GeV",     "", "", "", "", "", "", "cm", "cm", "cm", "cm", "cm", "cm", "", "", "", "", "", "", "GeV", "GeV"};
	const double HistMin[nDist] = { 12, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10, 10, 0, 30,   0, 0, 0, 0, 0, 0, 0, 10, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
	const double HistMax[nDist] = {300, 600, 300, 300, 300, 300, 300, 300, 300, 300, 300, 200, 200, 200, 300, 600, 10,10, 3.2, 3.2, 3.2, 40, 300, 400, 200, 200, 200, 2.5, 2.5, 2.5, 0.4, 0.4, 0.4,  0.1, 0.1, 0.1, 100, 100, 100, 2, 2, 2, 0.8, 0.8, 0.8, 0.05, 0.05, 0.05, 0.015, 0.015, 0.015, 1.2, 0.5, 0.5, 1.2, 0.5, 0.5, 1200, 1200};
	unsigned nBins[nDist];
	for(unsigned dist = 0; dist < nDist; ++dist) nBins[dist] = 20;
	nBins[16] = 10;
	nBins[17] = 10;
	//nBins[4] = 10;
	
	for(unsigned dist = 0; dist < nDist; ++dist){
		float BinWidth = (HistMax[dist] - HistMin[dist])/nBins[dist];
		std::ostringstream strs; strs << BinWidth; std::string Yaxis = strs.str();
		for(unsigned cat = 0; cat < nCat; ++cat){
			for(unsigned effsam = 0; effsam < nSamples_eff + 1; ++effsam){
				Histos[dist][cat][effsam] = new TH1D(catNames[cat] +  eff_names[effsam] + Histnames[dist], catNames[cat] +  eff_names[effsam] + Histnames[dist] + ";" + Xaxes[dist] + "; events /" + Yaxis + Units[dist], nBins[dist], HistMin[dist], HistMax[dist]);
				Histos[dist][cat][effsam]->Sumw2();
			}
		}
	}

	double maxBinC[nDist];
	for(unsigned dist = 0; dist < nDist; ++dist){
		maxBinC[dist] = Histos[dist][0][0]->GetBinCenter(Histos[dist][0][0]->GetNbinsX());
	}
	/*
	//Define 2D histograms
	const unsigned nDist2D = 17;
	TH2D *hists2D[nDist2D][nCat][nSamples_eff + 1];
	const TString distNames2D[nDist] = {"mtminVStrailPt", "mtminDeltaRVStrailPt", "mtminVSm3l", "mtminDeltaRVSm3l", "trailPtVSm3l", "m3lVSmet", "mtminVSmet", "mtminDeltaRVSmet", "trailPtVSmet", "trailPtVSLeadPt", "trailPtVSsubPt", "subPtVSleadPt", "mtminVSleadPt", "mtminDeltaRVSleadPt", "mtminVSmll", "mtminDeltaRVSmll", "m3lVSmll"};
	const double xMin[nDist2D] = {0, 0, 0, 0, 20, 100, 0, 0, 20, 0, 0, 0, 0, 0, 0, 0, 0};
	const double xMax[nDist2D] = {300, 300, 300, 300, 200, 700, 300, 300, 200, 200, 200, 200, 300, 300, 300, 300, 600};
	const double yMin[nDist2D] = {20, 20, 100, 100, 100, 0, 0, 0, 0, 35, 0, 35, 35, 35, 0, 0, 0};
	const double yMax[nDist2D] = {200, 200, 700, 700, 700, 300, 300, 300, 300, 300, 200, 300, 300, 300, 300, 300, 300};
	const unsigned nBinsX[nDist2D] = {40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40};
	const unsigned nBinsY[nDist2D] = {40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40};
	const TString xTitle[nDist2D] = {"M_{T}(other min(M_{OS}) (GeV)", "M_{T}(other min #Delta R) (GeV)", "M_{T}(other min(M_{OS}) (GeV)", "M_{T}(other min #Delta R) (GeV)", "P_{T}^{cone}(trailing) (GeV)", "M_{3l} (GeV)",  "M_{T}(other min(M_{OS}) (GeV)", "M_{T}(other min #Delta R) (GeV)", "P_{T}^{cone}(trailing) (GeV)", "P_{T}^{cone}(trailing) (GeV)", "P_{T}^{cone}(trailing) (GeV)", "P_{T}^{cone}(subleading) (GeV)", "M_{T}(other min(M_{OS}) (GeV)", "M_{T}(other min #Delta R) (GeV)", "M_{T}(other min(M_{OS}) (GeV)", "M_{T}(other min #Delta R) (GeV)","M_{3l} (GeV)"}; 
	const TString yTitle[nDist2D] = {"P_{T}^{cone}(trailing) (GeV)", "P_{T}^{cone}(trailing) (GeV)", "M_{3l} (GeV)", "M_{3l} (GeV)", "M_{3l} (GeV)", "MET (GeV)", "MET (GeV)", "MET (GeV)","MET (GeV)", "P_{T}^{cone}(leading) (GeV)", "P_{T}^{cone}(subleading) (GeV)", "P_{T}^{cone}(leading) (GeV)", "P_{T}^{cone}(leading) (GeV)", "P_{T}^{cone}(leading) (GeV)", "M_{ll} (GeV)", "M_{ll} (GeV)", "M_{ll} (GeV)"};
	for(unsigned dist = 0; dist < nDist2D; ++dist){
		for(unsigned cat = 0; cat < nCat; ++cat){
			for(unsigned effsam = 0; effsam < nSamples_eff + 1; ++effsam){
				hists2D[dist][cat][effsam] = new TH2D(catNames[cat] +  eff_names[effsam] + distNames2D[dist], catNames[cat] +  eff_names[effsam] + distNames2D[dist] + ";" + xTitle[dist] + ";" + yTitle[dist], nBinsX[dist], xMin[dist], xMax[dist], nBinsY[dist], yMin[dist], yMax[dist]);
				hists2D[dist][cat][effsam]->Sumw2();
			}
		}
	}
	double maxBinC2D[nDist][2];
	for(unsigned dist = 0; dist < nDist2D; ++ dist){
		maxBinC2D[dist][0] = xMax[dist] - (0.5*xMin[dist]/( (double) nBinsX[dist]) );
		maxBinC2D[dist][1] = yMax[dist] - (0.5*yMin[dist]/( (double) nBinsY[dist]) );
	}
	*/
    Double_t scale[nSamples -1];
	//set to check which data events have already been processed
	std::set<std::tuple<unsigned long, unsigned long, unsigned long> > usedEvents; //runNb, lumiBlock, eventNb	
	//Loop over all samples
	for(unsigned sam = 0, effsam = 0; sam < nSamples; ++sam, ++effsam){
		if(sam != 0){
   			if(names[sam] == names[sam -1]) --effsam;
    	}		
    	Long64_t nEntries = inputTree[sam]->GetEntries();
    	if(sam > 0){
   			scale[sam -1] = xSections[sam -1]*DataLuminosity*1000/(hcounter[sam]);
    	}
		cout << eff_names[effsam] << endl;
    	std::cout<<"Entries in "<< fileList[sam] <<" "<<nEntries<<std::endl;
		cout << effsam << endl;
		
        for(Long64_t it = 0; it < nEntries; ++it){
        	inputTree[sam]->GetEntry(it);
        	if (it%10000 == 0) cout<<'.'<<flush;
        	if(TestRun && it > 10000) break;
        	double scal;
        	if(effsam == 0) scal = 1;
        	else{
				scal = scale[sam-1]*_weight;
			}		
			//Use HNL ID
			cutBased();
			//Baseline event selection
			unsigned _nJets_uncleaned = _nJets;
			if(!baseline(true, true, false, false)) continue;
			//Check if data events were used before
			if(effsam == 0){
				auto event = usedEvents.find(std::make_tuple(_eventNb, _lumiBlock, _runNb));
				if(event != usedEvents.end()) continue;
				usedEvents.insert(std::make_tuple(_eventNb, _lumiBlock, _runNb));
			}
			//Categorize according to the number of leptons and flavors
			unsigned* ind = new unsigned[_nL];
			unsigned lCount = lepOrder(ind, 3, true, true);
			if(lCount != 3) continue; //Veto 4th FO lepton considering signal model!
			//MC prompt matching
			if(effsam  > nSig){
				bool promptfail = false;
				for(unsigned l = 0; l < lCount; ++l){	
					//cout << _origin[ind[l]] << endl;
					if(_origin[ind[l]] != 0){
						promptfail = true;
						break;
					}
				}
				if(promptfail) continue;
			}
			//Analysis pt cuts 

			
			//Require 3 leptons to be tight in data and MC, and determine nonPrompt bkg in data
			unsigned nTight = tightCount(ind, lCount);
			bool tightFail = nTight < 3;
			//Apply FR maps to data control region
			double* conePt = new double[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				conePt[l] = _lPt[ind[l]]*(1 + std::max(_isolation[ind[l]] - 0.1, 0.));
			}
			//index used to fill events, needed to separate fakes from data
			unsigned fill = effsam;
			//Apply FR maps to data control region
			if(tightFail && (effsam == 0 || effsam > nSig) ){
				//fakes go in different histogram
				fill = nSamples_eff;
				//Apply FR maps
				if(effsam != 0) scal*=-1;
				scal*= fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap, lCount);
			} else if(tightFail) continue;
			//Apply triggers to data events;
			if(effsam == 0 || effsam > nSig){
				bool trigPass[4];
				trigPass[0] = _lowM_trigger_eee || _lowM_trigger_all;
				trigPass[1] = _lowM_trigger_mee || _lowM_trigger_all;
				trigPass[2] = _lowM_trigger_mme || _lowM_trigger_all;
				trigPass[3] = _lowM_trigger_mmm || _lowM_trigger_all;
				if(!trigPass[tril_flavorComb(ind, _flavors, lCount)]) continue;	
			}
			//Analysis Pt Cuts
			if(!ptCuts_hnl(ind,lCount)) continue;
			//determine search category
			unsigned cat = hnl::cat(ind, _flavors, _charges, lCount, conePt[0]);
			if(cat == 999){
				continue;	
			}
			//Apply ID and reco SF to simulation
			if(effsam != 0){
				//Apply ID and reco SF to simulation
				for(unsigned l = 0; l < lCount; ++l){  //CHANGE BACK BACK BACK
					if(_istight[ind[l]]){
						if(_flavors[ind[l]] == 2){
							scal*=0.83; //Twiki suggests flat 90% id SF for taus
						} else if(_flavors[ind[l]] == 0){
							scal*=idTightSFMap[0]->GetBinContent(idTightSFMap[0]->FindBin(TMath::Min(_lPt[ind[l]], 199.), fabs(_lEta[ind[l]])));
							scal*=baseidSFMap[0]->GetBinContent(baseidSFMap[0]->FindBin(TMath::Min(_lPt[ind[l]], 199.), fabs(_lEta[ind[l]])));
							scal*=recSFMap_ele->GetBinContent(recSFMap_ele->FindBin(_lEta[ind[l]]));
						} else if(_flavors[ind[l]] == 1){
							double mapPt = _lPt[ind[l]] > 10 ? _lPt[ind[l]] : 10.1;
							scal*=idTightSFMap[1]->GetBinContent(idTightSFMap[1]->FindBin(TMath::Min(mapPt, 119.), fabs(_lEta[ind[l]])));
							scal*=baseidSFMap[1]->GetBinContent(baseidSFMap[1]->FindBin(TMath::Min(mapPt, 119.), fabs(_lEta[ind[l]])));
							if(_lPt[l] > 10) scal*=recSFMap_mu_ptAbove10->Eval(_lEta[ind[l]]);
							else scal*=recSFMap_mu_ptBelow10->Eval(_lEta[ind[l]]);
						}
					} else if(_isFO[ind[l]]){
						;
					} else if(_isloose[ind[l]]){
						;
					}
				}			
				//Apply btag SF
				for(unsigned j = 0; j < _nJets_uncleaned; ++j){
					//scal*=reader.eval_auto_bounds("central", BTagEntry::FLAV_B, _jetEta[j], _jetPt[j], _csv[j]);
				}
				//Apply PU reweighing
				scal*= PUweights->GetBinContent(PUweights->FindBin(_n_trueInteractions));
			}
			//determine which leptons will be used for the calculation of mll
			unsigned mllI[2] = {99, 99};
			TLorentzVector* lepV = new TLorentzVector[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				lepV[l].SetPtEtaPhiE(conePt[l], _lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]*(conePt[l]/_lPt[ind[l]]));
			}
			mllIndices(mllI, ind, lepV, _charges, _flavors, lCount);				
			//determine mll
			double mll;
			if(mllI[0] == 99){
				if(mllI[1] != 99) std::cerr << "error one mll index is not -1 while the other is" << endl;
				mll = -1;
			} else{
				TLorentzVector lzV[2];
				for(unsigned l = 0; l < 2; ++l) lzV[l].SetPtEtaPhiE(_lPt[mllI[l]]*(1 + std::max(_isolation[mllI[l]] - 0.1, 0.)), _lEta[mllI[l]], _lPhi[mllI[l]], _lE[mllI[l]]*(1 + std::max(_isolation[mllI[l]] - 0.1, 0.)) );
				mll = (lzV[0] + lzV[1]).M();
			}
			//Determine W lepton indec
			unsigned lw = 9999;
			for(unsigned l = 0; l < lCount; ++l){
				if(ind[l] != mllI[0] && ind[l] != mllI[1]){
					lw = ind[l];
				}
			}

			//Calculate lepton system vector
			TLorentzVector lepSyst;
			for(int l = 0; l < 3; ++l) lepSyst += lepV[l];

			if(cat < 4){
				if(conePt[0] < 15) continue;
				if(conePt[1] < 10) continue;
				if(conePt[2] < 10 - 5*(_flavors[ind[2]])) continue;	//5(10) GeV cut on trailing muon(electron)
				//Calculate minDeltaR and maxDeltaR
				double minDeltaR = 9999.;
				double maxDeltaR = 0.;
				for(unsigned l = 0; l < lCount - 1; ++l){
					for(unsigned k = l + 1; k < lCount; ++k){
						if(lepV[l].DeltaR(lepV[k]) > maxDeltaR) maxDeltaR = lepV[l].DeltaR(lepV[k]);
						if(lepV[l].DeltaR(lepV[k]) < minDeltaR) minDeltaR = lepV[l].DeltaR(lepV[k]);
					}
				}
				if((cat == 0 || cat == 2) && minDeltaR < 0.05) continue;
				if((cat == 0 || cat == 2) && maxDeltaR < 2) continue;
				if(lepSyst.M() > 80) continue;
				if(_met > 75) continue;
			}
			//Analysis Pt thresholds
			else if(cat > 3){
				if(conePt[0] < 55) continue;
				if(conePt[1] < 15) continue;
				//if(conePt[2] < 10) continue;
				if(conePt[2] < 10 - 5*(_flavors[ind[2]])) continue;
				if(cat == 4){
					if(fabs(mll - 91) < 15) continue; //Veto onZ events (WZ and DY)
					if(fabs(lepSyst.M() - 91) < 15) continue; //veto conversions
				}
			}

			//Calculate MET vector
			TLorentzVector METvec;
			METvec.SetPtEtaPhiE(_met, 0, _met_phi, _met);
			//Calculate min(Mos) and Mos(min Delta R)
			double minMos = 0;
			unsigned minI[2] = {99, 99};
			double MosminDeltaR;
			unsigned minDeltaRI[2] = {99, 99};
			double minDeltaR = 99999.;
			for(unsigned l = 0; l < lCount -1 ; ++l){
				for(unsigned k = l + 1; k < lCount; ++k){
					if(_charges[ind[l]] != _charges[ind[k]]){
						if( (lepV[l] + lepV[k]).M() < minMos  || minMos == 0){
							minMos = (lepV[l] + lepV[k]).M();
							minI[0] = l;
							minI[1] = k;
						}
						if( lepV[l].DeltaR(lepV[k]) < minDeltaR){
							minDeltaR = lepV[l].DeltaR(lepV[k]);
							MosminDeltaR = (lepV[l] + lepV[k]).M();
							minDeltaRI[0] = l;
							minDeltaRI[1] = k;
						}
					}
				}
			}
			//Calculate alternete mt values
			unsigned lw_min, lw_minDeltaR;
			for(unsigned l = 0; l < lCount; ++l){
				if(l != minI[0] && l != minI[1]){
					lw_min = l;
				}
				if(l != minDeltaRI[0] && l != minDeltaRI[1]){
					lw_minDeltaR = l;
				}
			}
			double mt_min = transmass(lepV[lw_min], METvec);
			double mt_minDeltaR = transmass(lepV[lw_minDeltaR], METvec);
			//lepton used in MT
			TLorentzVector Wlep;
			Wlep.SetPtEtaPhiE(_lPt[lw]*(1 + std::max(_isolation[lw] - 0.1, 0.)), _lEta[lw], _lPhi[lw], _lE[lw]*(1 + std::max(_isolation[lw] - 0.1, 0.)));
			//Calculate mt3l
			double mt3l = transmass(lepV[0] + lepV[1] + lepV[2], METvec);
			//calculate mt2_ss
			double mt2_ss = transmass(Wlep, METvec);
			for(unsigned l = 0; l < lCount -1; ++l){
				for(unsigned k = l + 1; k < lCount; ++k){
					if(_charges[ind[k]] == _charges[ind[l]]){
						if(_flavors[ind[k]] == _flavors[ind[l]]){
							mt2_ss = mt2ll(lepV[k], lepV[l], METvec);
						}
					}
				}
			}
			//recontruct neutrino pt z
			/*
			const double mw = 80;
			TLorentzVector nu;
			double m2 = 0.5*mw*mw + _lPt[lw]*_met;
			double solplus = (m2/(_lPt[lw]*_lPt[lw]))*( Wlep.Pz() + fabs(Wlep.P())*sqrt(1 - (_met*_met*_lPt[lw]*_lPt[lw])/(m2*m2) ) );
			double solmin = (m2/(_lPt[lw]*_lPt[lw]))*( Wlep.Pz() - fabs(Wlep.P())*sqrt(1 - (_met*_met*_lPt[lw]*_lPt[lw])/(m2*m2) ) );
			double nupx = _met*cos(_met_phi);
			double nupy = _met*sin(_met_phi);
			TLorentzVector vecplus, vecmin;
			vecplus.SetPxPyPzE(nupx, nupy, solplus, sqrt(_met*_met + solplus*solplus));	
			vecmin.SetPxPyPzE(nupx, nupy, solmin, sqrt(_met*_met + solmin*solmin));
			double mindiff = fabs((vecmin + Wlep).M() - mw);
			if( fabs( (vecplus + Wlep).M() - mw) < mindiff ){
				nu = vecplus;
			} else{
				nu = vecmin;
			}
			TLorentzVector third;	
			if(_charges[mllI[0]] != _charges[lw]) third.SetPtEtaPhiE(_lPt[mllI[0]], _lEta[mllI[0]], _lPhi[mllI[0]], _lE[mllI[0]]);
			else third.SetPtEtaPhiE(_lPt[mllI[1]], _lEta[mllI[1]], _lPhi[mllI[1]], _lE[mllI[1]]);
			double mnu = (nu + Wlep + third).M();
			*/
			const double mw = 80;
			double minDiff = 99999.;
			double bestSol = 0.;
			unsigned nuI = 99;
			for(unsigned l = 0; l < lCount; ++l){
				double m2 = 0.5*mw*mw + conePt[l]*_met;
				double solplus = (m2/(conePt[l]*conePt[l]))*(lepV[l].Pz() + fabs(lepV[l].P())*sqrt(1 - (_met*_met*conePt[l]*conePt[l])/(m2*m2) ) );
				double solmin = (m2/(conePt[l]*conePt[l]))*(lepV[l].Pz() - fabs(lepV[l].P())*sqrt(1 - (_met*_met*conePt[l]*conePt[l])/(m2*m2) ) );
				double nupx = _met*cos(_met_phi);
				double nupy = _met*sin(_met_phi);
				TLorentzVector vecplus, vecmin;
				vecplus.SetPxPyPzE(nupx, nupy, solplus, sqrt(_met*_met + solplus*solplus));	
				vecmin.SetPxPyPzE(nupx, nupy, solmin, sqrt(_met*_met + solmin*solmin));
				if(fabs((vecplus + lepV[l]).M() - mw) < minDiff){
					minDiff = fabs((vecplus + lepV[l]).M() - mw);
					bestSol = solplus;
					nuI = l;
				}
				if(fabs((vecmin + lepV[l]).M() - mw) < minDiff){
					minDiff = fabs((vecmin + lepV[l]).M() - mw);
					bestSol = solmin;
					nuI = l;
				}
			}
			TLorentzVector nu;
			nu.SetPxPyPzE(_met*cos(_met_phi), _met*sin(_met_phi), bestSol, sqrt(_met*_met + bestSol*bestSol));
			double mnu = 0;
			for(unsigned l = 0; l < lCount; ++l){
				if(l == nuI) continue;
				if(_charges[ind[l]] != _charges[nuI]){
					mnu = (lepV[l] + lepV[nuI] + nu).M();
					break;
				}
			}
				
			//Calculate scalar sum of lepton Pt's
			double ptSum = 0;
			for(unsigned l = 0; l < lCount; ++l){
				ptSum += conePt[l];
			}
			//calculate mass of 3l + nu system:
			double m4 = (nu  + lepSyst).M();
						
			//fill 1D histograms
			double values[nDist] = { mll, lepSyst.M(), mt3l, mnu,  minMos, mt_min, MosminDeltaR, mt_minDeltaR, transmass(Wlep, METvec), mt2_ss,  mt2_maxPt(ind, _charges, lepV, METvec, lCount), _lPt[ind[0]], _lPt[ind[1]], _lPt[ind[2]], _met, _HT, static_cast<double>(_nJets), static_cast<double>(_n_bJets), fabs(METvec.DeltaPhi(lepV[0])), fabs(METvec.DeltaPhi(lepV[1])), fabs(METvec.DeltaPhi(lepV[2])), static_cast<double>(_n_PV), lepSyst.Pt(), ptSum,  conePt[0], conePt[1], conePt[2], fabs(_lEta[ind[0]]), fabs(_lEta[ind[1]]), fabs(_lEta[ind[2]]),  _miniisolation[ind[0]][0], _miniisolation[ind[1]][0], _miniisolation[ind[2]][0], _isolation[ind[0]],  _isolation[ind[1]],  _isolation[ind[2]], _ptrel[ind[0]], _ptrel[ind[1]], _ptrel[ind[2]], _ptratio[ind[0]], _ptratio[ind[1]], _ptratio[ind[2]], _closeJetCSVAll[ind[0]], _closeJetCSVAll[ind[1]], _closeJetCSVAll[ind[2]], fabs(_3dIP[ind[0]]), fabs(_3dIP[ind[1]]), fabs(_3dIP[ind[2]]), fabs(_ipPV[ind[0]]), fabs(_ipPV[ind[1]]), fabs(_ipPV[ind[2]]), conePt[0]/(lepSyst.M()), conePt[1]/(lepSyst.M()), conePt[2]/(lepSyst.M()), conePt[0]/(m4), conePt[1]/(m4), conePt[2]/(m4), m4, (lepV[0] + lepV[1] + nu).M()};

			for(unsigned dist = 0; dist < nDist; ++dist){
				Histos[dist][cat][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal);
			}
			/*
			//fill 2D histograms
			double valuesX[nDist2D] = {mt_min, mt_minDeltaR, mt_min, mt_minDeltaR, conePt[2], lepSyst.M(), mt_min, mt_minDeltaR, conePt[2], conePt[2], conePt[2], conePt[1], mt_min, mt_minDeltaR, mt_min, mt_minDeltaR, lepSyst.M()};
			double valuesY[nDist2D] = {conePt[2],conePt[2], lepSyst.M(), lepSyst.M(), lepSyst.M(), _met, _met, _met, _met, conePt[0], conePt[1], conePt[0], conePt[0], conePt[0], mll, mll, mll};
			for(unsigned dist = 0; dist < nDist2D; ++dist){
				hists2D[dist][cat][fill]->Fill(TMath::Min(valuesX[dist], maxBinC2D[dist][0]), TMath::Min(valuesY[dist], maxBinC2D[dist][1]), scal);
			}	
			*/		
    	}
	}

	//Split data and MC histograms for plotting and propagating uncertainties
	TH1D* dataYields[nDist][nCat];
	for(unsigned dist = 0; dist < nDist; ++dist){
		for(unsigned cat = 0; cat < nCat; ++cat){
			dataYields[dist][cat] = (TH1D*) Histos[dist][cat][nSig + 1]->Clone();
		}
	}
	TH1D* bkgYields[nDist][nCat][nSamples_eff -5]; //change to nSamples_eff if sig is removed
	for(unsigned dist = 0; dist < nDist; ++dist){
		for(unsigned cat = 0; cat < nCat; ++cat){
			for(unsigned effsam = nSig + 1; effsam < nSamples_eff + 1; ++effsam){
				bkgYields[dist][cat][effsam -nSig - 1] = (TH1D*) Histos[dist][cat][effsam]->Clone();
				if(effsam > nSig + 1){
					dataYields[dist][cat]->Add(bkgYields[dist][cat][effsam -nSig -1]);
				}
			}
		}
	}
	TString sigNames[nSig] = {"m_{N} = 5 GeV", "m_{N} = 10 GeV", "m_{N} = 20 GeV", "m_{N} = 30 GeV", "m_{N} = 40 GeV", "m_{N} = 50 GeV", "m_{N} = 60 GeV", "m_{N} = 70 GeV", "m_{N} = 90 GeV", "m_{N} = 100 GeV", "m_{N} = 150 GeV", "m_{N} = 200 GeV", "m_{N} = 300 GeV", "m_{N} = 400 GeV", "m_{N} = 500 GeV", "m_{N} = 700 GeV", "m_{N} = 1000 GeV"};
	/*
	for(unsigned sig = 0; sig < nSig; ++sig){
		sigNames[sig] += ", |V_{eff}|^{2} = ";
	}
	for(unsigned sig = 0; sig < 6; ++sig){
		sigNames[sig] += "10^{-5}";
	}
	sigNames[6] += "10^{-4}";
	for(unsigned sig = 7; sig < nSig; ++sig){
		sigNames[sig] += "10^{-2}";
	}
	*/
	
	for(unsigned sig = 0; sig < nSig; ++sig){
		sigNames[sig] += ", normalized to bkg.";
	}
	
	const TString procNames[nSamples_eff + 1 - nSig] = {"total BKG", "ZZ/H", "triboson", "WZ", "X + #gamma", "TT/T + X",  "non-prompt"};
	//Plot the yields as a function of the search region
	for(unsigned dist = 0; dist < nDist; ++dist){
		for(unsigned cat = 0; cat < nCat; ++cat){
			unsigned nSigToPlot;
			if(cat < 4) nSigToPlot = 7;
			else nSigToPlot = 7;
			TH1D* signals[nSigToPlot];
			TString sigPlotNames[nSigToPlot];
			if(cat < 4){
				for(unsigned sig = 0; sig < nSigToPlot; ++sig){
					signals[sig] = (TH1D*) Histos[dist][cat][sig + 1]->Clone();
					sigPlotNames[sig] = sigNames[sig];
				}
			} else{
				for(unsigned sig = nSig - nSigToPlot; sig < nSig; ++sig){
					signals[sig- nSig + nSigToPlot] = (TH1D*) Histos[dist][cat][sig + 1]->Clone();
					sigPlotNames[sig- nSig + nSigToPlot] = sigNames[sig];
				}
			}
			plotDataVSMC(dataYields[dist][cat], bkgYields[dist][cat], procNames, nSamples_eff - nSig, Histnames[dist] + "_" +  catNames[cat] + extra, false, 0, "HNL", true, signals,  sigPlotNames ,nSigToPlot, true);
		}
	}
	//Make slides
	/*
	ofstream slidedump;
	slidedump.open("./textdumps/slidedump.txt");
	unsigned plotCounter = 0;
	for(unsigned cat = 0; cat < nCat; ++cat){
		for(unsigned dist = 0; dist < nDist; ++dist){
			if(plotCounter%4 == 0){
				if(dist != 0){
					slidedump <<"\\end{figure} \n";
					slidedump <<"\\end{frame} \n";
					slidedump <<" \n";
				}
				slidedump << "\\begin{frame}{" << catNames[cat] << " plots} \n";
				slidedump << "\\begin{figure} \n";
			}
			slidedump << "\\subfloat{\\includegraphics[width = .4\\textwidth]{" + Histnames[dist] + "_" +  catNames[cat] + extra + ".pdf}}";
			if((plotCounter - 1)%2 == 0 && dist != 0) slidedump << "\\\\ \n";
			else slidedump << "\n";
			++plotCounter;
		}
		slidedump <<"\\end{figure} \n";
		slidedump <<"\\end{frame} \n";
		slidedump <<" \n";
		plotCounter -= plotCounter%4;
	}
	if(plotCounter%4 != 0){
		slidedump <<"\\end{figure} \n";
		slidedump <<"\\end{frame} \n";
		slidedump <<" \n";
	}
	*/
	//Plot 2D histograms
	/*
	TH2D* bkg2D[nDist][nCat];
	TH2D* signal2D[nDist][nCat][nSig];
	//ofstream slidedump2D;
	//slidedump2D.open("./textdumps/slidedump2D.txt");
	const TString signN[nSig] = {"mn5", "mn15", "mn20", "mn30", "mn40", "mn50", "mn60", "mn80", "mn100", "mn130", "mn150", "mm200", "mn300", "mn400"};
	//const TString masses[5] = {"$m_{N} = 100 \\mathrm{GeV}$", "$m_{N} = 150 \\mathrm{GeV}$", "$m_{N} = 200 \\mathrm{GeV}$", "$m_{N} = 300 \\mathrm{GeV}$", "$m_{N} = 400 \\mathrm{GeV}$"};
	for(unsigned dist = 0; dist < nDist2D; ++ dist){
		for(unsigned cat = 0; cat < nCat; ++cat){
			bkg2D[dist][cat] = (TH2D*) hists2D[dist][cat][nSig + 1]->Clone();
			for(unsigned effsam = nSig + 2; effsam < nSamples_eff + 1; ++ effsam){
				bkg2D[dist][cat]->Add(hists2D[dist][cat][effsam]);
			}
			*/
			/*
			slidedump2D << "\\begin{frame}{" << catNames[cat] << " plots} \n";
			slidedump2D << "\\begin{minipage}{.3\\textwidth} \n";
			slidedump2D << "\\scriptsize {\\color{blue} $\\Sigma(BKG)$}: \n";
			slidedump2D << "\\begin{figure} \n";
			slidedump2D << "\\includegraphics[width = \\textwidth]{" + distNames2D[dist] +catNames[cat] + "BKG" + extra + ".pdf} \n";
			slidedump2D << "\\end{figure} \n";
			slidedump2D << "\\end{minipage} \n";
			*/
			/*
			plotHist(bkg2D[dist][cat], distNames2D[dist] + catNames[cat] + "BKG" + extra, true);
			for(unsigned sig = 0; sig < nSig; ++sig){
				signal2D[dist][cat][sig] = (TH2D*) hists2D[dist][cat][sig + 1]->Clone();
				plotHist(signal2D[dist][cat][sig], distNames2D[dist] + "SIG_" + signN[sig] + catNames[cat] + extra, true);
				*/
				/*
				slidedump2D << "\\begin{minipage}{.3\\textwidth} \n";
				slidedump2D << "\\scriptsize {\\color{blue}" << masses[sig] << "}: \n";
				slidedump2D << "\\begin{figure} \n";
				slidedump2D << "\\includegraphics[width = \\textwidth]{" + distNames2D[dist] + "SIG_" + signN[sig] + catNames[cat] + extra + ".pdf}\n";
				slidedump2D << "\\end{figure} \n";
				slidedump2D << "\\end{minipage} \n";
				if(sig == 1) slidedump2D << "\\\\ \n";
				*/
			/*
			}
			//slidedump2D << "\\end{frame} \n";
		}
	}
	*/
	//slidedump2D.close();
}



int main(int argc, char* argv[]){
	TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;
	testtree.Loop();
	rootapp->Run();
    return 0;
}


