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

//include code to calculate btag SF
//#include "../bTag/BTagCalibrationStandalone.h"
//include other parts of the code
#include "MultilepSUSYfunc.h"
#include "tdrstyle.h"
#include "plotCode.h"
#include "trilTree.h"
#include "hnlTools.h"


void trilTree::hnlPlots(){
	//Set plotting style
	//	setTDRStyle();
	//Supress plot pop-up window
	//gROOT->SetBatch(kTRUE);
	//Define all samples and their cross_sections
	
	const unsigned nSamples = 50;
	const unsigned nSamples_eff = 24;
	const unsigned nSig = 18;
	const double lowMCoupling = 0.00001;
	const double highMCoupling = 0.01;
	const double couplingCorrection = 10000;
	//~~~~~~~~ background normalizations~~~~~~~~~~~~~~~~~~~~
	const double glugluToZZkFactor = 2.1;
	const double WZSF = 0.655;
	const double ZZSF = 1.032;
	const double XgammaSF = 0.945;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	const TString fileList[nSamples] = {"data_combined_trilepton.root", 
										"HeavyNeutrino_M1_2l.root", "HeavyNeutrino_M2_2l.root", "HeavyNeutrino_M5_2l.root", "HeavyNeutrino_M10_2l.root", "HeavyNeutrino_M20_2l.root", "HeavyNeutrino_M30_2l.root", "HeavyNeutrino_M40_2l.root", "HeavyNeutrino_M50_2l.root", "HeavyNeutrino_M60_2l.root", "HeavyNeutrino_M80_2l.root", "HeavyNeutrino_M100_2l.root", "HeavyNeutrino_M130_2l.root", "HeavyNeutrino_M150_2l.root", "HeavyNeutrino_M200_2l.root", "HeavyNeutrino_M400_2l.root", "HeavyNeutrino_M600_2l.root", "HeavyNeutrino_M800_2l.root", "HeavyNeutrino_M1000_2l.root",
 										"ZZTo4L.root",  "GluGluToZZTo4mu.root", "GluGluToZZTo4e.root", "GluGluToZZTo4tau.root", "GluGluToZZTo2e2mu.root", "GluGluToZZTo2e2tau.root", "GluGluToZZTo2mu2tau.root", "VHToNonbb.root", "GluGluHToZZTo4L_M125.root", "VBF_HToZZTo4L_M125.root", "WWG.root","WWW.root", "WWZ.root", "WWTo2L2Nu_DoubleScattering.root", "WWTo2L2Nu.root",  "ZZZ.root", "WZTo3LNu_mllmin01.root", "TTGJets.root","ZGTo2LG.root", "WGToLNuG.root", "TGJets.root", "TTJets_DiLept.root", "TTJets_SingleLeptFromTbar.root", "TTJets_SingleLeptFromT.root", "DYJetsToLL_M10to50.root", "DYJetsToLL_M50.root", "ttHToNonbb.root", "TTWJetsToLNu.root", "TTZToLLNuNu.root", "TTZToLL_M1to10.root",  "TTTT.root"};
	
	const double xSections[nSamples - 1] = { 5.237e-01*couplingCorrection*lowMCoupling, 5.238e-01*couplingCorrection*lowMCoupling, 5.215e-01*couplingCorrection*lowMCoupling, 5.132e-01*couplingCorrection*lowMCoupling, 4.764e-01*couplingCorrection*lowMCoupling, 4.160e-01*couplingCorrection*lowMCoupling, 3.335e-01*couplingCorrection*lowMCoupling, 2.364e-01*couplingCorrection*lowMCoupling, 1.336e-01*couplingCorrection*lowMCoupling, 2.473e-02*couplingCorrection*lowMCoupling*10, 1.033e-03*couplingCorrection*highMCoupling, 2.874e-04*couplingCorrection*highMCoupling, 1.561e-04*couplingCorrection*highMCoupling, 4.929e-05*couplingCorrection*highMCoupling, 3.505e-06*couplingCorrection*highMCoupling, 6.922e-07*couplingCorrection*highMCoupling, 2.013e-07*couplingCorrection*highMCoupling, 7.239e-08*couplingCorrection*highMCoupling, 
											1.256*ZZSF, 0.00159*glugluToZZkFactor*ZZSF, 0.00159*glugluToZZkFactor*ZZSF, 0.00159*glugluToZZkFactor*ZZSF, 0.00319*glugluToZZkFactor*ZZSF, 0.00319*glugluToZZkFactor*ZZSF, 0.00319*glugluToZZkFactor*ZZSF,  0.9561, 0.01212, 0.001034,  0.2147, 0.2086, 0.1651, 0.1729,   12.178, 0.01398, 58.59*WZSF, 3.697, 123.9*XgammaSF, 489*XgammaSF,  2.967, 87.315, 182.175, 182.175, 18610, 1921.8*3,  0.215, 0.2043, 0.2529, 0.0493, 0.009103}; 
	const TString names[nSamples] = {"data", "m_{N} = 1 GeV", "m_{N} = 2 GeV", "m_{N} = 5 GeV", "m_{N} = 10 GeV", "m_{N} = 20 GeV", "m_{N} = 30 GeV", "m_{N} = 40 GeV", "m_{N} = 50 GeV", "m_{N} = 60 GeV", "m_{N} = 80 GeV", "m_{N} = 100 GeV", "m_{N} = 130 GeV", "m_{N} = 150 GeV", "m_{N} = 200 GeV", "m_{N} = 400 GeV", "m_{N} = 600 GeV", "m_{N} = 800 GeV", "m_{N} = 1000 GeV",  									 
									 "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "triboson", "triboson", "triboson", "triboson", "triboson", "triboson", "WZ", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X"};
	
	//Read and initialize scale factors used in the analysis 
	readSF(true);
	//loop over all samples
	for(unsigned sam = 0; sam < nSamples; ++sam){
		//if(sam  < 1) continue;
		Loop(fileList[sam], xSections[sam], sam < 1, sam != 0 && sam <= nSig);
	}
	const TString eff_names[nSamples_eff + 1] = {"data", "m_{N} = 1 GeV", "m_{N} = 2 GeV", "m_{N} = 5 GeV", "m_{N} = 10 GeV", "m_{N} = 20 GeV", "m_{N} = 30 GeV", "m_{N} = 40 GeV", "m_{N} = 50 GeV", "m_{N} = 60 GeV", "m_{N} = 80 GeV", "m_{N} = 100 GeV", "m_{N} = 130 GeV", "m_{N} = 150 GeV", "m_{N} = 200 GeV", "m_{N} = 400 GeV", "m_{N} = 600 GeV", "m_{N} = 800 GeV", "m_{N} = 1000 GeV",
												 "ZZ/H", "triboson", "WZ", "X + #gamma", "TT/T + X",  "non-prompt"};
}


//loop over all events in an individual samples
void trilTree::Loop(const TString& fileName, const double xSection, const bool isData, const bool isSignal){
	cout << fileName << endl;
	//Read Trees from ROOT files
	TFile* hfile = new TFile("../data_april17/"+fileName,"read");
	hfile->cd("FakeElectrons");
	TTree* inputTree;
	double hcounter;
	//Determine hcounter for cross section scaling
	TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
	_hCounter->Read("hCounter");
	hcounter = _hCounter->GetBinContent(1);
    inputTree = (TTree*) (hfile->Get("FakeElectrons/fakeTree"));
	Init(inputTree, false, !isData);
	
	//Tweakable options////////////////////////////////////////////////////
	const bool TestRun = false;	//Break after a few events
	const double DataLuminosity = 35.867;    //units of fb^{-1}
	const TString extra = "";	//for plot file names
	//////////////////////////
	const unsigned nCat = 8;  //Number of categories
	const TString catNames[nCat] = {"lowM_3lOSSF_lowPt", "lowM_3lnoOSSF_lowPt", "lowM_3lOSSF_highPt", "lowM_3lnoOSSF_highPt", "highM_3lOSSF", "highM_3lnoOSSF", "baseline_OSSF", "baseline_noOSSF"};
	const unsigned nUnc = 6;//number of shape uncertainties
	const TString uncNames[nUnc] = {"jec", "metUnclustered", "pdf", "scale", "pu", "btagSF"};
	const unsigned nDist = 45;  //Number of distributions to plot	 	
	const TString suffix[2] = {"", "non-Prompt"};
	
	/*
	TH1D* histos[nCat][nDist][2];	//Kinematic distributions to plot
	//Kinematic shape uncertainty histograms
	TH1D* histosDown[nUnc][nCat][nDist][2];  //nSamples_eff +_ 1
	TH1D* histosUp[nUnc][nCat][nDist][2];
	TH1D* histosPdfVar[100][nCat][nDist][2];
	*/
	TH1D**** histos = new TH1D***[nCat]; //Kinematic distributions to plot
	for(unsigned cat = 0; cat < nCat; ++cat){
		histos[cat] = new TH1D**[nDist];
		for(unsigned dist = 0; dist < nDist; ++dist){
			histos[cat][dist] = new TH1D*[2]; //One entry for signal band one for sideband
		}
	}
	//Kinematic shape uncertainty histograms
	TH1D***** histosDown = new TH1D****[nUnc];
	TH1D***** histosUp = new TH1D****[nUnc];
	TH1D***** histosPdfVar = new TH1D****[100];
	if(!isData && !isSignal){
		for(unsigned unc = 0; unc < nUnc; ++unc){
			histosDown[unc] = new TH1D***[nCat];
			histosUp[unc] = new TH1D***[nCat];
			for(unsigned cat = 0; cat < nCat; ++cat){
				histosDown[unc][cat] = new TH1D**[nDist];
				histosUp[unc][cat] = new TH1D**[nDist];
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosDown[unc][cat][dist] = new TH1D*[2];
					histosUp[unc][cat][dist] = new TH1D*[2];
				}
			}
		}
		for(unsigned pdf = 0; pdf < 100; ++pdf){
			histosPdfVar[pdf] = new TH1D***[nCat];
			for(unsigned cat = 0; cat < nCat; ++cat){
				histosPdfVar[pdf][cat] = new TH1D**[nDist];
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosPdfVar[pdf][cat][dist] = new TH1D*[2];
				}
			}
		}
	}
	//Names of the distributions to plot
	const TString histNames[nDist] = {"Mll", "M3l", "minMos", "mt_minMos", "MosminDeltaR", "mt_minDeltaR", "mt", "mt2_ss", "mt2_maxPt", "MET", "HT", "NJets", "NbJets", "DeltaPhi_lepMET_le", "DeltaPhi_lepMET_sub", "DeltaPhi_lepMET_tr", "Nvtx","Pt_trilep", "ConePt_le", "ConePt_sub", "ConePt_tr", "Eta_le", "Eta_sub", "Eta_tr", "MiniIso_le", "MiniIso_sub", "MiniIso_tr", "RelIso_le", "RelIso_sub", "RelIso_tr", "Ptrel_le", "Ptrel_sub", "Ptrel_tr", "Ptratio_le", "Ptratio_sub", "Ptratio_tr", "csv_le", "csv_sub", "csv_tr", "dxy_le", "dxy_sub", "dxy_tr", "dz_le", "dz_sub", "dz_tr"};
  	//X-axis labels of distributions to plot
	const TString xAxes[nDist] = {"M_{ll}(GeV)", "M_{3l} (GeV)", "min(M_{OS}) (GeV)", "M_{T}(other min(M_{OS})) (GeV)", "M_{OS}(min #Delta R) (GeV)", "M_{T}(other min #Delta R) (GeV)",  "M_{T} (GeV)", "M_{T2}(SS) (GeV)",  "M_{T2}(max P_{T} 2l) (GeV)", "MET (GeV)", "HT (GeV)", "number of jets", "number of b-jets", "#Delta#Phi(leading, MET)", "#Delta#Phi(subleading, MET)", "#Delta#Phi(trailing, MET)", "Number of vertices", "P_{T}(3l) (GeV)", "P_{T}^{cone}(leading) (GeV)", "P_{T}^{cone}(subleading) (GeV)", "P_{T}^{cone}(trailing) (GeV)", "|#eta(leading)|","|#eta(subleading)|", "|#eta(trailing)|", "miniIso(leading)", "miniIso(subleading)", "miniIso(trailing)", "relIso(leading)", "relIso(subleading)", "relIso(trailing)", "P_{T}^{rel}(leading) (GeV)", "P_{T}^{rel}(subleading) (GeV)", "P_{T}^{rel}(trailing) (GeV)", "P_{T}^{ratio}(leading)", "P_{T}^{ratio}(subleading)", "P_{T}^{ratio}(trailing)", "closest Jet CSV(leading)", "closest Jet CSV(subleading)", "closest Jet CSV(trailing)", "|d_{xy}(leading)| (cm)", "|d_{xy}(subleading)| (cm)", "|d_{xy}(trailing)| (cm)", "|d_{z}(leading)| (cm)", "|d_{z}(subleading)| (cm)", "|d_{z}(trailing)| (cm)"};
	//Units of distributions to plot.
	const TString units[nDist] = {"GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "", "", "", "", "", "",  "GeV", "GeV", "GeV", "GeV", "", "", "", "", "", "", "", "", "",  "GeV", "GeV", "GeV",     "", "", "", "", "", "", "cm", "cm", "cm", "cm", "cm", "cm"};
	//Minimum x-range of histograms.
	const double histMin[nDist] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 30, 0, 0, 0, 0, 0, 0, 0, 15, 10, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
	//Minimum y-range of histograms.
	const double histMax[nDist] = {300, 600, 150, 300, 300, 300, 300, 300, 300, 300, 600, 10,10, 3.2, 3.2, 3.2, 40, 300, 200, 200, 200, 2.5, 2.5, 2.5, 0.5, 0.5, 0.5,  0.1, 0.1, 0.1, 100, 100, 100, 2, 2, 2, 0.8, 0.8, 0.8, 0.015, 0.015, 0.015, 0.15, 0.15, 0.15};
	unsigned nBins[nDist];
	for(unsigned dist = 0; dist < nDist; ++dist) nBins[dist] = 20;
	nBins[2] = 15;
	nBins[11] = 10;
	nBins[12] = 10;
	nBins[18] = 37;
	nBins[19] = 38;
	nBins[20] = 39;//39	
	//Initialize histograms of kinematic distributions
	for(unsigned dist = 0; dist < nDist; ++dist){
		float binWidth = (histMax[dist] - histMin[dist])/nBins[dist];
		std::ostringstream strs; strs << binWidth; std::string yAxis = strs.str();
		for(unsigned cat = 0; cat < nCat; ++cat){
			for(unsigned f = 0; f < 2; ++f){
				if(isSignal == true && f != 0) break;
				histos[cat][dist][f] = new TH1D(histNames[dist] + catNames[cat]+ suffix[f], histNames[dist] + catNames[cat] + suffix[f] + ";" + xAxes[dist] + "; Events /" + yAxis + units[dist], nBins[dist], histMin[dist], histMax[dist]);
				histos[cat][dist][f]->Sumw2();
				if(!isData && !isSignal){
					for(unsigned unc = 0; unc < nUnc; ++unc){
						histosDown[unc][cat][dist][f] = new TH1D(histNames[dist] + catNames[cat] + uncNames[unc] + "Down" + suffix[f], histNames[dist] + catNames[cat]  + uncNames[unc] + "Down" + suffix[f] + ";" + xAxes[dist] + "; Events /" + yAxis + units[dist], nBins[dist], histMin[dist], histMax[dist]);
						histosUp[unc][cat][dist][f] = new TH1D(histNames[dist] + catNames[cat] + uncNames[unc] + "Up" + suffix[f], histNames[dist] + catNames[cat]  + uncNames[unc] + "Up" + suffix[f] + ";" + xAxes[dist] + "; Events /" + yAxis + units[dist], nBins[dist], histMin[dist], histMax[dist]);
					}
					for(unsigned pdf = 0; pdf < 100; ++pdf){
						histosPdfVar[pdf][cat][dist][f] = new TH1D(histNames[dist] + catNames[cat] + "_pdf" + std::to_string(pdf) + suffix[f], histNames[dist] + catNames[cat] + "_pdf" + std::to_string(pdf) + suffix[f] + ";" + xAxes[dist] + "; Events /" + yAxis + units[dist], nBins[dist], histMin[dist], histMax[dist]);
					}
				}
			}
		}
	}
	double maxBinC[nDist];
	for(unsigned dist = 0; dist < nDist; ++dist){
		maxBinC[dist] = histos[0][dist][0]->GetBinCenter(histos[0][dist][0]->GetNbinsX());
	}

	//root file used for memory control, i.e. temporarily storing the histograms not to fill up the ram
	//TFile* storageFile = new TFile("tempStorage_" + fileName, "recreate");
	
	double scale;
	if(!isData) scale = xSection*DataLuminosity/hcounter;
	//set to check which data events have already been processed
	std::set<std::tuple<unsigned long, unsigned long, unsigned long> > usedEvents; //runNb, lumiBlock, eventNb	
	//Loop over all samples
   	Long64_t nEntries = inputTree->GetEntries();
    for(Long64_t it = 0; it < nEntries; ++it){
       	inputTree->GetEntry(it);
        if (it%10000 == 0) cout<<'.'<< std::flush;
        if(TestRun && it > 10000) break;
        double scal;
        if(isData) scal = 1;
        else scal = scale*_weight;	
		cutBased();
		//Baseline event selection
		if(!baseline(true, true, false, false)) continue;
		if(!isData &&  nBJets(true, false, 2) != 0) continue;
		else if(isData && nBJets(true, false, 0) != 0) continue;

		//Check if data events were used before
		if(isData){
			auto event = usedEvents.find(std::make_tuple(_eventNb, _lumiBlock, _runNb));
			if(event != usedEvents.end()) continue;
			usedEvents.insert(std::make_tuple(_eventNb, _lumiBlock, _runNb));
		}
		//Categorize according to the number of leptons and flavors
		unsigned* ind = new unsigned[_nL];
		unsigned lCount = lepOrder(ind, 3, true, true);	
		if(lCount != 3) continue; //Veto 4th FO lepton considering signal model!
		//Apply analysis Pt thresholds
		if(!ptCuts_hnl(ind,lCount)) continue;
		//MC prompt matching
		if(!isData){
			bool promptfail = false;
			for(unsigned l = 0; l < lCount; ++l){	
				if(_origin[ind[l]] != 0){
					promptfail = true;
					break;
				}
			}
			if(promptfail) continue;
		}
		//Require 3 leptons to be tight in data and MC, and determine nonPrompt bkg in data
		unsigned nTight = tightCount(ind, lCount);
		bool tightFail = nTight < 3;
		//index used to fill events, needed to separate fakes from data
		//Apply FR maps to data control region
		//Calculate and store conePt
		double* conePt = new double[lCount];
		for(unsigned l = 0; l < lCount; ++l){
			conePt[l] = _lPt[ind[l]]*(1 + std::max(_isolation[ind[l]] - 0.1, 0.));
		}
		unsigned fill = 0;
		if(tightFail && !isSignal){ //&& effsam == 0){
			//fakes go in different histogram
			fill = 1;
			//Apply FR maps
			if(!isData) scal *= -1;
			scal*= fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap, lCount);
		} else if(tightFail) continue;

		//Sample overlap removal
		if(fileName == "WGToLNuG.root"){
			bool promptfail =true;
			for(unsigned l = 0; l < lCount; ++l){	
				if(_pdgmc[ind[l]] == 22 && _originPhot[ind[l]] == 0){
					promptfail = false;
					break;
				}
			}
			if(promptfail) continue;
		}
		if(conePt[2] >= 10){
			if(fileName == "TTJets_DiLept.root" || fileName  == "DYJetsToLL_M10to50.root" || fileName  == "DYJetsToLL_M50.root"  || fileName  == "TTJets_SingleLeptFromTbar.root"  || fileName  == "TTJets_SingleLeptFromT.root" ) continue;
		} else if(conePt[2] < 10){
			if(fileName  == "TTJets_DiLept.root" || fileName  == "DYJetsToLL_M10to50.root" || fileName  == "DYJetsToLL_M50.root"  || fileName  == "TTJets_SingleLeptFromTbar.root"  || fileName  == "TTJets_SingleLeptFromT.root" ){
				bool promptfail = false;
				for(unsigned l = 0; l < lCount; ++l){	
					if(_pdgmc[ind[l]] == 22 && _originPhot[ind[l]] == 0){
						promptfail = true;
						break;
					}
				}
				if(promptfail) continue;	
			}
		}				
		if(fileName  == "TTGJets.root" || fileName  == "ZGTo2LG.root"){
			bool promptfail =true;
			for(unsigned l = 0; l < lCount; ++l){	
				if(_pdgmc[ind[l]] == 22 && _originPhot[ind[l]] == 0){
					promptfail = false;
					break;
				}
			}
			if(promptfail) continue;
		}
		//Apply triggers to data events;
		bool trigPass[4];
		trigPass[0] = _lowM_trigger_eee || _lowM_trigger_all;
		trigPass[1] = _lowM_trigger_mee || _lowM_trigger_all;
		trigPass[2] = _lowM_trigger_mme || _lowM_trigger_all;
		trigPass[3] = _lowM_trigger_mmm || _lowM_trigger_all;
		if(!trigPass[tril_flavorComb(ind, _flavors, lCount)]) continue;	

		//determine search category
		TLorentzVector* lepV = new TLorentzVector[lCount];
		for(unsigned l = 0; l < lCount; ++l){
			lepV[l].SetPtEtaPhiE(_lPt[ind[l]]*(1 + std::max(_isolation[ind[l]] - 0.1, 0.)), _lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]*(1 + std::max(_isolation[ind[l]] - 0.1, 0.)) );
		}
		//Calculate lepton system vector
		TLorentzVector lepSyst;
		for(int l = 0; l < 3; ++l) lepSyst += lepV[l];

		unsigned cat = hnl::cat(ind, _flavors, _charges, lCount, conePt[0]);
		if(cat == 999){
			continue;	
		}
		//Apply ID and reco SF to simulation
		if(!isData){
			scal*=getEventSF(ind, lCount, true);
		}
		//determine which leptons will be used for the calculation of mll
		unsigned mllI[2] = {99, 99};
		mllIndices(mllI, ind, lepV, _charges, _flavors, lCount);				
		//determine mll
		double mll;
		if(mllI[0] == 99){
			if(mllI[1] != 99) std::cerr << "error one mll index is not -1 while the other is" << endl;
			mll = -1;
		} else{
			TLorentzVector lzV[2];
			for(unsigned l = 0; l < 2; ++l) lzV[l].SetPtEtaPhiE(PtCone(_lPt[mllI[l]]*(1 + std::max(_isolation[mllI[l]] - 0.1, 0.)), _flavors[mllI[l]], _lepMVA[mllI[l]], _ptratio[mllI[l]]), _lEta[mllI[l]], _lPhi[mllI[l]], _lE[mllI[l]]*(1 + std::max(_isolation[mllI[l]] - 0.1, 0.)) );
			mll = (lzV[0] + lzV[1]).M();
		}
		//End of baseline seletion! fill total histograms here
		//determine the index of the W lepton
		unsigned lw = 9999;
		for(unsigned l = 0; l < lCount; ++l){
			if(ind[l] != mllI[0] && ind[l] != mllI[1]){
				lw = ind[l];
			}
		}
		TLorentzVector Wlep;
		Wlep.SetPtEtaPhiE(_lPt[lw]*std::max(1., 1 + (_isolation[lw] - 0.1)), _lEta[lw], _lPhi[lw], _lE[lw]*std::max(1., 1 + (_isolation[lw] - 0.1)));
		TLorentzVector METvec;
		METvec.SetPtEtaPhiE(_met, 0, _met_phi, _met);

		TLorentzVector lzV[2];
		for(unsigned l = 0; l < 2; ++l) lzV[l].SetPtEtaPhiE(_lPt[mllI[l]]*std::max(1., 1 + (_isolation[mllI[l]] - 0.1)), _lEta[mllI[l]], _lPhi[mllI[l]], _lE[mllI[l]]*std::max(1., 1 + (_isolation[mllI[l]] - 0.1)));
		mll = (lzV[0] + lzV[1]).M();
		///////////////////////////////////////////////////
			
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
		//Clean jets before plotting
		double HT = 0;
		unsigned nJets = 0;
		for(unsigned j = 0; j < _nJets; ++j){
			if(jetIsClean(j)){
				if(_jetPt[j] > 25){
					++nJets;
					if(_jetPt[j] > 30){
						HT += _jetPt[j];
					}
				}
			}
		}
		//Scale weights for glugluToZZ
		if(fileName.Contains("GluGluToZZTo")){
			for(unsigned pdf = 0; pdf < 110; ++pdf){
				_scaleWeight[pdf] = 1;
			}
		}
		double values[nDist] = { mll, lepSyst.M(), minMos, mt_min, MosminDeltaR, mt_minDeltaR, transmass(Wlep, METvec), mt2_ss,  mt2_maxPt(ind, _charges, lepV, METvec, lCount), _met, HT, static_cast<double>(nJets), static_cast<double>(nBJets(true, false, 0)), fabs(METvec.DeltaPhi(lepV[0])), fabs(METvec.DeltaPhi(lepV[1])), fabs(METvec.DeltaPhi(lepV[2])), static_cast<double>(_n_PV), lepSyst.Pt(),  _lPt[ind[0]]*std::max(1., 1 + (_isolation[ind[0]] - 0.1)), _lPt[ind[1]]*std::max(1., 1 + (_isolation[ind[1]] - 0.1)), _lPt[ind[2]]*std::max(1., 1 + (_isolation[ind[2]] - 0.1)), fabs(_lEta[ind[0]]), fabs(_lEta[ind[1]]), fabs(_lEta[ind[2]]),  _miniisolation[ind[0]][0], _miniisolation[ind[1]][0], _miniisolation[ind[2]][0], _isolation[ind[0]],  _isolation[ind[1]],  _isolation[ind[2]], _ptrel[ind[0]], _ptrel[ind[1]], _ptrel[ind[2]], _ptratio[ind[0]], _ptratio[ind[1]], _ptratio[ind[2]], _closeJetCSVAll[ind[0]], _closeJetCSVAll[ind[1]], _closeJetCSVAll[ind[2]], fabs(_ipPV[ind[0]]), fabs(_ipPV[ind[1]]), fabs(_ipPV[ind[2]]), fabs(_ipZPV[ind[0]]), fabs(_ipZPV[ind[1]]), fabs(_ipZPV[ind[2]])}; 

		//Fill total and per category kinematic distributions
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//Check if event passes category selection
		//Category spefific selection:
		//low mass general
		bool catPass = true;
		//OSSF or no OSSF
		bool ossfI = (hnl::flavorComposition(ind, _flavors, _charges, lCount) == 0);
		if(cat < 4){
			//Calculate minDeltaR and maxDeltaR
			double minDeltaR = 9999.;
			double maxDeltaR = 0.;
			for(unsigned l = 0; l < lCount - 1; ++l){
				for(unsigned k = l + 1; k < lCount; ++k){
					if(lepV[l].DeltaR(lepV[k]) > maxDeltaR) maxDeltaR = lepV[l].DeltaR(lepV[k]);
					if(lepV[l].DeltaR(lepV[k]) < minDeltaR) minDeltaR = lepV[l].DeltaR(lepV[k]);
				}
			}
			if((cat == 0 || cat == 2) && minDeltaR < 0.05) catPass = false;
			if((cat == 0 || cat == 2) && maxDeltaR < 2) catPass = false;
			if((cat == 0 || cat == 2) && !vetoLowMll(5)) catPass = false; //NEW
			if(lepSyst.M() > 80) catPass = false;
		}
		else if(cat > 3){
			if(conePt[1] < 15) catPass = false;
			if(cat == 4){
				if(fabs(mll - 91) < 15) catPass = false; //Veto onZ events (WZ and DY)
				if(fabs(lepSyst.M() - 91) < 15) catPass = false; //veto conversions
			}
			if(!vetoLowMll(5)) catPass = false;
		}
		//Nominal yields
		if(nBJets(true, false, 0) == 0){
			for(unsigned dist = 0; dist < nDist; ++dist){
				histos[6 + ossfI][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal);
			}
			if(catPass && (cat > 3 || _met < 75)){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histos[cat][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal);
				}
			}
		}
		if(isSignal || isData) continue;
		//vary JEC down
		if(nBJets(true, false, 1) == 0){
			METvec.SetPtEtaPhiE(_metJECDown, 0, _met_phiJECDown, _metJECDown);
			values[3] = transmass(lepV[lw_min], METvec);
			values[5] = transmass(lepV[lw_minDeltaR], METvec);
			values[6] = transmass(Wlep, METvec);
			values[7] = transmass(Wlep, METvec);
			for(unsigned l = 0; l < lCount -1; ++l){
				for(unsigned k = l + 1; k < lCount; ++k){
					if(_charges[ind[k]] == _charges[ind[l]]){
						if(_flavors[ind[k]] == _flavors[ind[l]]){
							values[7] = mt2ll(lepV[k], lepV[l], METvec);
						}
					}
				}
			}
			values[8] = mt2_maxPt(ind, _charges, lepV, METvec, lCount);
			values[9] = _metJECDown;
			HT = 0;
			nJets = 0;
			for(unsigned j = 0; j < _nJets; ++j){
				if(jetIsClean(j)){
						if(_jetPtDown[j] > 25){
							++nJets;
						if(_jetPtDown[j] > 30){
							HT += _jetPtDown[j];
						}
					}
				}
			}
			values[10] = HT;
			values[11] = static_cast<double>(nJets);
			values[12] = static_cast<double>(nBJets(true, false, 1));
			values[13] = fabs(METvec.DeltaPhi(lepV[0]));
			values[14] = fabs(METvec.DeltaPhi(lepV[1]));
			values[15] = fabs(METvec.DeltaPhi(lepV[2]));
			for(unsigned dist = 0; dist < nDist; ++dist){
				histosDown[0][6+ ossfI][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal);
			}
			if(catPass && (cat > 3 || _metJECDown < 75)){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosDown[0][cat][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal);
				}
			}
		}
		//vary JEC up
		if(nBJets(true, false, 2) == 0){
			METvec.SetPtEtaPhiE(_metJECUp, 0, _met_phiJECUp, _metJECUp);
			values[3] = transmass(lepV[lw_min], METvec);
			values[5] = transmass(lepV[lw_minDeltaR], METvec);
			values[6] = transmass(Wlep, METvec);
			values[7] = transmass(Wlep, METvec);
			for(unsigned l = 0; l < lCount -1; ++l){
				for(unsigned k = l + 1; k < lCount; ++k){
					if(_charges[ind[k]] == _charges[ind[l]]){
						if(_flavors[ind[k]] == _flavors[ind[l]]){
							values[7] = mt2ll(lepV[k], lepV[l], METvec);
						}
					}
				}
			}
			values[8] = mt2_maxPt(ind, _charges, lepV, METvec, lCount);
			values[9] = _metJECUp;
			HT = 0;
			nJets = 0;
			for(unsigned j = 0; j < _nJets; ++j){
				if(jetIsClean(j)){
					if(_jetPtUp[j] > 25){
						++nJets;
						if(_jetPtUp[j] > 30){
							HT += _jetPtUp[j];
						}
					}
				}
			}
			values[10] = HT;
			values[11] = static_cast<double>(nJets);
			values[12] = static_cast<double>(nBJets(true, false, 2));
			values[13] = fabs(METvec.DeltaPhi(lepV[0]));
			values[14] = fabs(METvec.DeltaPhi(lepV[1]));
			values[15] = fabs(METvec.DeltaPhi(lepV[2]));
			for(unsigned dist = 0; dist < nDist; ++dist){
				histosUp[0][6+ ossfI][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal);
			}
			if(catPass && (cat > 3 || _metJECUp < 75)){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosUp[0][cat][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal);
				}
			}
		}
		//nominal b-veto 
		if(nBJets(true, false, 0) != 0) continue;
		//nominal HT and nJets values
		HT = 0;
		nJets = 0;
		for(unsigned j = 0; j < _nJets; ++j){
			if(jetIsClean(j)){
				if(_jetPt[j] > 25){
					++nJets;
					if(_jetPt[j] > 30){
						HT += _jetPt[j];
					}
				}
			}
		}
		values[10] = HT;
		values[11] = static_cast<double>(nJets);
		values[12] = static_cast<double>(nBJets(true, false, 0));
					
		//vary unclustered met down
		METvec.SetPtEtaPhiE(_metOtherDown, 0, _met_phiOtherDown, _metOtherDown);
		values[3] = transmass(lepV[lw_min], METvec);
		values[5] = transmass(lepV[lw_minDeltaR], METvec);
		values[6] = transmass(Wlep, METvec);
		values[7] = transmass(Wlep, METvec);
		for(unsigned l = 0; l < lCount -1; ++l){
			for(unsigned k = l + 1; k < lCount; ++k){
				if(_charges[ind[k]] == _charges[ind[l]]){
				if(_flavors[ind[k]] == _flavors[ind[l]]){
						values[7] = mt2ll(lepV[k], lepV[l], METvec);
					}
				}
			}
		}
		values[8] = mt2_maxPt(ind, _charges, lepV, METvec, lCount);
		values[9] = _metOtherDown;
		values[13] = fabs(METvec.DeltaPhi(lepV[0]));
		values[14] = fabs(METvec.DeltaPhi(lepV[1]));
		values[15] = fabs(METvec.DeltaPhi(lepV[2]));
		for(unsigned dist = 0; dist < nDist; ++dist){
			histosDown[1][6+ ossfI][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal);
		}
		if(catPass && (cat > 3 || _metOtherDown < 75)){
			for(unsigned dist = 0; dist < nDist; ++dist){
				histosDown[1][cat][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal);
			}
		}
		
		//vary unclustered met up
		METvec.SetPtEtaPhiE(_metOtherUp, 0, _met_phiOtherUp, _metOtherUp);
		values[3] = transmass(lepV[lw_min], METvec);
		values[5] = transmass(lepV[lw_minDeltaR], METvec);
		values[6] = transmass(Wlep, METvec);
		values[7] = transmass(Wlep, METvec);
		for(unsigned l = 0; l < lCount -1; ++l){
			for(unsigned k = l + 1; k < lCount; ++k){
				if(_charges[ind[k]] == _charges[ind[l]]){
					if(_flavors[ind[k]] == _flavors[ind[l]]){
					values[7] = mt2ll(lepV[k], lepV[l], METvec);
					}
				}
			}
		}
		values[8] = mt2_maxPt(ind, _charges, lepV, METvec, lCount);
		values[9] = _metOtherDown;
		values[13] = fabs(METvec.DeltaPhi(lepV[0]));
		values[14] = fabs(METvec.DeltaPhi(lepV[1]));
		values[15] = fabs(METvec.DeltaPhi(lepV[2]));
		for(unsigned dist = 0; dist < nDist; ++dist){
			histosUp[1][6+ ossfI][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal);
		}
		if(catPass && (cat > 3 || _metOtherUp < 75)){
			for(unsigned dist = 0; dist < nDist; ++dist){
				histosUp[1][cat][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal);
			}
		}
		//Nominal values to fill
		METvec.SetPtEtaPhiE(_met, 0, _met_phi, _met);
		values[3] = transmass(lepV[lw_min], METvec);
		values[5] = transmass(lepV[lw_minDeltaR], METvec);
		values[6] = transmass(Wlep, METvec);
		values[7] = transmass(Wlep, METvec);
		for(unsigned l = 0; l < lCount -1; ++l){
			for(unsigned k = l + 1; k < lCount; ++k){
				if(_charges[ind[k]] == _charges[ind[l]]){
					if(_flavors[ind[k]] == _flavors[ind[l]]){
						values[7] = mt2ll(lepV[k], lepV[l], METvec);
					}
				}
			}
		}
		values[8] = mt2_maxPt(ind, _charges, lepV, METvec, lCount);
		values[9] = _met;
		values[13] = fabs(METvec.DeltaPhi(lepV[0]));
		values[14] = fabs(METvec.DeltaPhi(lepV[1]));
		values[15] = fabs(METvec.DeltaPhi(lepV[2]));	
			
		//vary scale down
		for(unsigned dist = 0; dist < nDist; ++dist){
			histosDown[3][6+ ossfI][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal*_scaleWeight[8]);
		}
		if(catPass && (cat > 3 || _met < 75)){
			for(unsigned dist = 0; dist < nDist; ++dist){
				histosDown[3][cat][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal*_scaleWeight[8]);
			}
		}
		//vary scale up
		for(unsigned dist = 0; dist < nDist; ++dist){
			histosUp[3][6+ ossfI][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal*_scaleWeight[4]);
		}
		if(catPass && (cat > 3 || _met < 75)){
			for(unsigned dist = 0; dist < nDist; ++dist){
				histosUp[3][cat][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal*_scaleWeight[4]);
			}
		}	
		//vary pu down
		for(unsigned dist = 0; dist < nDist; ++dist){
			histosDown[4][6+ ossfI][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), (scal/PUweights[0]->GetBinContent(PUweights[0]->FindBin( std::min(_n_trueInteractions, 49) )))*PUweights[1]->GetBinContent(PUweights[1]->FindBin(std::min(_n_trueInteractions, 49)  )) );
		}
		if(catPass && (cat > 3 || _met < 75)){
			for(unsigned dist = 0; dist < nDist; ++dist){
				histosDown[4][cat][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), (scal/PUweights[0]->GetBinContent(PUweights[0]->FindBin( std::min(_n_trueInteractions, 49) )))*PUweights[1]->GetBinContent(PUweights[1]->FindBin(std::min(_n_trueInteractions, 49)  )) );
			}
		}	
		//vary pu up
		for(unsigned dist = 0; dist < nDist; ++dist){
			histosUp[4][6+ ossfI][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), (scal/PUweights[0]->GetBinContent(PUweights[0]->FindBin( std::min(_n_trueInteractions, 49) )))*PUweights[2]->GetBinContent(PUweights[2]->FindBin(std::min(_n_trueInteractions, 49)  )) );
		}
		if(catPass && (cat > 3 || _met < 75)){
			for(unsigned dist = 0; dist < nDist; ++dist){
				histosUp[4][cat][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), (scal/PUweights[0]->GetBinContent(PUweights[0]->FindBin( std::min(_n_trueInteractions, 49) )))*PUweights[2]->GetBinContent(PUweights[2]->FindBin(std::min(_n_trueInteractions, 49)  )) );
			}
		}	
		//vary btag SF down	
		for(unsigned dist = 0; dist < nDist; ++dist){
			histosDown[5][6+ ossfI][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), (scal/bTagSF(true, 0))*bTagSF(true, 1) );
		}
		if(catPass && (cat > 3 || _met < 75)){
			for(unsigned dist = 0; dist < nDist; ++dist){
				histosDown[5][cat][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), (scal/bTagSF(true, 0))*bTagSF(true, 1) );
			}
		}	
		//vary btag SF up
		for(unsigned dist = 0; dist < nDist; ++dist){
			histosUp[5][6+ ossfI][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), (scal/bTagSF(true, 0))*bTagSF(true, 2) );
		}
		if(catPass && (cat > 3 || _met < 75)){
			for(unsigned dist = 0; dist < nDist; ++dist){
				histosUp[5][cat][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), (scal/bTagSF(true, 0))*bTagSF(true, 2) );
			}
		}
		//all pdf variations
		for(unsigned pdf = 0; pdf < 100; ++pdf){
			for(unsigned dist = 0; dist < nDist; ++dist){
				histosPdfVar[pdf][6+ ossfI][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal*_scaleWeight[9 + pdf] );
			}
			if(catPass && (cat > 3 || _met < 75)){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosPdfVar[pdf][cat][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal*_scaleWeight[9 + pdf] );
				}
			}
		}
		//Delete dynamical arrays
		delete[] conePt;
		delete[] ind;
		delete[] lepV;
	}		
	
	//Set negative bins to 0 before adding other processes
	if(!isData){
		for(unsigned cat = 0; cat < nCat; ++cat){
			for(unsigned dist = 0; dist < nDist; ++dist){
				for(unsigned bin = 1; bin < histos[cat][dist][0]->GetNbinsX() + 1; ++bin){
					if(histos[cat][dist][0]->GetBinContent(bin) < 0 ) histos[cat][dist][0]->SetBinContent(bin, 0.);
					if(!isSignal){
						for(unsigned unc = 0; unc < nUnc; ++unc){
							if(histosUp[unc][cat][dist][0]->GetBinContent(bin) < 0. ) histosUp[unc][cat][dist][0]->SetBinContent(bin, 0.);
							if(histosDown[unc][cat][dist][0]->GetBinContent(bin) < 0. ) histosDown[unc][cat][dist][0]->SetBinContent(bin, 0.);
						}	
						for(unsigned pdf = 0; pdf < 100; ++pdf){
							if(histosPdfVar[pdf][cat][dist][0]->GetBinContent(bin) < 0. ) histosPdfVar[pdf][cat][dist][0]->SetBinContent(bin, 0.);
						}
					}
				}
			}
		}
	}
	/*
	//Write histograms to temporary root file
	for(unsigned cat = 0; cat < nCat; ++cat){
		for(unsigned dist = 0; dist < nDist; ++dist){
			for(unsigned f = 0; f < 2; ++f){
				if(isSignal && f != 0) break;
				histos[cat][dist][f]->Write(catNames[cat] + histNames[dist] + fileName + suffix[f]);
				if(!isData && !isSignal){
					for(unsigned unc = 0; unc < nUnc; ++unc){
						histosDown[unc][cat][dist][f]->Write(catNames[cat] + histNames[dist] + fileName + uncNames[unc] + "Down" + suffix[f]);
						histosDown[unc][cat][dist][f]->Delete();
						histosUp[unc][cat][dist][f]->Write(catNames[cat] + histNames[dist] + fileName + uncNames[unc] + "Up" + suffix[f] );
						histosUp[unc][cat][dist][f]->Delete();
					}
					for(unsigned pdf = 0; pdf < 100; ++pdf){
						histosPdfVar[pdf][cat][dist][f]->Write(catNames[cat] + histNames[dist] + fileName + std::to_string(pdf) + suffix[f]);
						histosPdfVar[pdf][cat][dist][f]->Delete();
					}
				}
			}
		}
	}	
	*/
	//Delete histograms to free up ram
	if(!isData && !isSignal){
		cout << "DELETING UNCERTAINTY HISTOGRAMS" << endl;
		for(unsigned unc = 0; unc < nUnc; ++unc){
			for(unsigned cat = 0; cat < nCat; ++cat){
				for(unsigned dist = 0; dist < nDist; ++dist){					
					//delete[] histosDown[unc][cat][dist];
					//histosDown[unc][cat][dist] = nullptr;
					//free(histosDown[unc][cat][dist]);
					//delete[] histosUp[unc][cat][dist];
					//histosUp[unc][cat][dist] = nullptr;
					//free(histosDown[unc][cat][dist]);
				}
				//delete[] histosDown[unc][cat];	
				//histosDown[unc][cat] = nullptr;
				//free(histosDown[unc][cat]);
				//delete[] histosUp[unc][cat];
				//histosUp[unc][cat] = nullptr;
				//free(histosUp[unc][cat]);
			}
			//delete[] histosDown[unc];
			//histosDown[unc] = nullptr;
			//free(histosDown[unc]);
			//delete[] histosUp[unc];
			//histosUp[unc] = nullptr;
			//free(histosUp[unc]);
		}
		//delete[] histosDown;
		//histosDown = nullptr;
		//delete[] histosUp;
		//histosUp = nullptr;
		free(histosDown);
		free(histosUp);
		for(unsigned pdf = 0; pdf < 100; ++pdf){
			for(unsigned cat = 0; cat < nCat; ++cat){
				for(unsigned dist = 0; dist < nDist; ++dist){
					//delete[] histosPdfVar[pdf][cat][dist];
					//histosPdfVar[pdf][cat][dist] = nullptr;
					//free( histosPdfVar[pdf][cat][dist]);
				}
				//delete[] histosPdfVar[pdf][cat];
				//histosPdfVar[pdf][cat] = nullptr;
				//free(histosPdfVar[pdf][cat]);
			}
			//delete[] histosPdfVar[pdf];
			//histosPdfVar[pdf] = nullptr;
			//free(histosPdfVar[pdf]);
		}
		//delete[] histosPdfVar;
		//histosPdfVar = nullptr;
		free(histosPdfVar);
		cout << "DELETED UNCERTAINTY HISTOGRAMS" << endl;
	}
	
	for(unsigned cat = 0; cat < nCat; ++cat){
		for(unsigned dist = 0; dist < nDist; ++dist){		
			//delete[] histos[cat][dist];
		}
		//delete histos[cat];
	}
	//delete[] histos;
	free(histos);
	cout << "DELETED ALL NORMAL HISTOS" << endl;
	
	//storageFile->Close();
	//delete storageFile;
	
	//Delete file and tree
	/*
	delete _hCounter;
	delete inputTree;
	delete hfile;
	*/
	/*
	for(unsigned dist = 0; dist < nDist; ++dist){
		std::cout << "dist = " << dist << endl;
		for(unsigned cat = 0; cat < nCat; ++cat){
			std::cout << "cat = " << cat << endl;
			delete[] histos[cat][dist];
			if(!isData && !isSignal){
				for(unsigned unc = 0; unc < nUnc; ++unc){
					delete[] histosDown[unc][cat][dist];
					delete[] histosUp[unc][cat][dist];
				}
				for(unsigned pdf = 0; pdf < 100; ++pdf){
					delete[] histosPdfVar[pdf][cat][dist];
				}
			}
		}
	}
	*/
	//Delete inputTree
	//delete inputTree;
	//	cout << "DELETED inputTree" << endl;
	/*
	
	//Calculate rms of pdf shape variations
	for(unsigned effsam = 1; effsam < nSamples_eff; ++effsam){
		if(effsam <= nSig) continue;
		for(unsigned cat = 0; cat < nCat; ++cat){
			for(unsigned dist = 0; dist < nDist; ++dist){
				for(unsigned b = 1; b < histos[cat][dist][effsam]->GetNbinsX() + 1; ++b){
					double pdfVarRms = 0;
					for(unsigned pdf = 0; pdf < 100; ++pdf){
						pdfVarRms += (histos[cat][dist][effsam]->GetBinContent(b) - histosPdfVar[pdf][cat][dist][effsam]->GetBinContent(b))*(histos[cat][dist][effsam]->GetBinContent(b) - histosPdfVar[pdf][cat][dist][effsam]->GetBinContent(b));
					}
					pdfVarRms = 0.01*sqrt(pdfVarRms);
					histosDown[2][cat][dist][effsam]->SetBinContent(b, histos[cat][dist][effsam]->GetBinContent(b) - pdfVarRms);
					histosUp[2][cat][dist][effsam]->SetBinContent(b, histos[cat][dist][effsam]->GetBinContent(b) + pdfVarRms);
				}
			}
		}
	}	
	//Split data and MC histograms for plotting and propagating uncertainties
	//data and background kinematic distributions
	TH1D* dataHistos[nCat][nDist];
	for(unsigned cat = 0; cat < nCat; ++cat){
		for(unsigned dist = 0; dist < nDist; ++dist){
			dataHistos[cat][dist] = (TH1D*) histos[cat][dist][nSig + 1]->Clone();
		}
	}
		
	TH1D* bkgHistos[nCat][nDist][nSamples_eff - nSig]; //change to nSamples_eff if sig is removed	
	for(unsigned cat = 0; cat < nCat; ++cat){
		for(unsigned dist = 0; dist <nDist; ++dist){
			for(unsigned effsam = nSig + 1; effsam < nSamples_eff + 1; ++effsam){
				bkgHistos[cat][dist][effsam - nSig - 1] = (TH1D*) histos[cat][dist][effsam]->Clone();
				if(effsam > nSig + 1){
					dataHistos[cat][dist]->Add(bkgHistos[cat][dist][effsam - nSig - 1]);
				}
			}
		}
	}
	const double extraUnc[nSamples_eff] = {1.25, 1.5, 1.085, 1.15, 1.15, 1.3};
	double flatSyst[4] = {0.025, 0.04, 0.02};
	TH1D* bkgSystDist[nCat][nDist][nSamples_eff - nSig];
	for(unsigned bkg = 0; bkg < nSamples_eff - nSig; ++bkg){
		cout << "BACKGROUND = " << eff_names[nSig + 1 + bkg] << endl;
		for(unsigned cat = 0; cat < nCat; ++cat){
			for(unsigned dist = 0; dist < nDist; ++dist){
				bkgSystDist[cat][dist][bkg] = (TH1D*) bkgHistos[cat][dist][bkg]->Clone();
				for(unsigned b = 1; b < bkgHistos[cat][dist][bkg]->GetNbinsX() + 1; ++b){
					bkgSystDist[cat][dist][bkg]->SetBinContent(b,0);
					if(histos[cat][dist][bkg +nSig + 1]->GetBinContent(b) == 0){
						bkgSystDist[cat][dist][bkg]->SetBinContent(b,0);
					} else{
						if(bkg != nSamples_eff - nSig -1){
							double systbin = 0;
							//loop over shape uncertainties
							cout << "~~~~~~~~~~~~~~~~~~~~~~~" << endl;
							for(unsigned unc = 0; unc < nUnc; ++unc){
								double syst = std::max(fabs(histosUp[unc][cat][dist][bkg + nSig + 1]->GetBinContent(b) - histos[cat][dist][bkg + nSig + 1]->GetBinContent(b)), fabs(histosDown[unc][cat][dist][bkg + nSig + 1]->GetBinContent(b) - histos[cat][dist][bkg + nSig + 1]->GetBinContent(b)));
								systbin += syst*syst;	
								std::cout << "unc = " << unc << endl;
								std::cout << "histosUp[unc][cat][dist][bkg + nSig + 1]->GetBinContent(b) = " << histosUp[unc][cat][dist][bkg + nSig + 1]->GetBinContent(b) << endl;	
								std::cout << "histosDown[unc][cat][dist][bkg + nSig + 1]->GetBinContent(b)  = " << histosDown[unc][cat][dist][bkg + nSig + 1]->GetBinContent(b) << endl;
								std::cout << "histos[cat][dist][bkg + nSig + 1]->GetBinContent(b) = " << histos[cat][dist][bkg + nSig + 1]->GetBinContent(b) << endl;
								std::cout << "size of systematic " << uncNames[unc] << " = " << syst/histos[cat][dist][bkg + nSig + 1]->GetBinContent(b) << endl;
							}
							//loop over flat uncertainties
							for(unsigned unc = 0; unc < 4; ++unc){
								systbin += histos[cat][dist][bkg + nSig + 1]->GetBinContent(b)*histos[cat][dist][bkg + nSig + 1]->GetBinContent(b)*flatSyst[unc]*flatSyst[unc];						
							}
							bkgSystDist[cat][dist][bkg]->SetBinContent(b, sqrt(systbin));
						} else{
							bkgSystDist[cat][dist][bkg]->SetBinContent(b, histos[cat][dist][bkg + nSig + 1]->GetBinContent(b)*0.3);
						}
						bkgSystDist[cat][dist][bkg]->SetBinError(b, 0);
					}
				}
			}
		}
	}
	*/
	/*
	TString sigNames[nSig] = {"m_{N} = 1 GeV", "m_{N} = 2 GeV", "m_{N} = 5 GeV", "m_{N} = 10 GeV", "m_{N} = 20 GeV", "m_{N} = 30 GeV", "m_{N} = 40 GeV", "m_{N} = 50 GeV", "m_{N} = 60 GeV", "m_{N} = 80 GeV", "m_{N} = 100 GeV", "m_{N} = 130 GeV", "m_{N} = 150 GeV", "m_{N} = 200 GeV", "m_{N} = 400 GeV", "m_{N} = 600 GeV", "m_{N} = 800 GeV", "m_{N} = 1000 GeV"};
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
	/*
	const TString distNames[nSamples_eff + 1 - nSig] = {"total pred.", "ZZ/H", "triboson", "WZ", "X + #gamma", "TT/T + X",  "non-prompt"};
	
	//Plot the yields as a function of the search region
	for(unsigned cat = 0; cat < nCat; ++cat){
		const unsigned nSigToPlot = 6;
		TH1D* signals[nSigToPlot];
		TString sigPlotNames[nSigToPlot];
		if(cat < 4){
			sigPlotNames[0] = "m_{N} = 1 GeV, norm. to bkg.";
			sigPlotNames[1] = "m_{N} = 5 GeV, norm. to bkg.";
			sigPlotNames[2] = "m_{N} = 20 GeV, norm. to bkg.";
			sigPlotNames[3] = "m_{N} = 30 GeV, norm. to bkg.";
			sigPlotNames[4] = "m_{N} = 50 GeV, norm. to bkg.";
			sigPlotNames[5] = "m_{N} = 50 GeV, norm. to bkg.";
		} else if(cat < 6){
			sigPlotNames[0] = "m_{N} = 100 GeV, norm. to bkg.";
			sigPlotNames[1] = "m_{N} = 130 GeV, norm. to bkg.";
			sigPlotNames[2] = "m_{N} = 150 GeV, norm. to bkg.";
			sigPlotNames[3] = "m_{N} = 200 GeV, norm. to bkg.";
			sigPlotNames[4] = "m_{N} = 400 GeV, norm. to bkg.";
			sigPlotNames[5] = "m_{N} = 600 GeV, norm. to bkg.";
		} else{
			sigPlotNames[0] = "m_{N} = 5 GeV, norm. to bkg.";
			sigPlotNames[1] = "m_{N} = 20 GeV, norm. to bkg.";
			sigPlotNames[2] = "m_{N} = 40 GeV, norm. to bkg.";
			sigPlotNames[3] = "m_{N} = 100 GeV, norm. to bkg.";
			sigPlotNames[4] = "m_{N} = 200 GeV, norm. to bkg.";
			sigPlotNames[5] = "m_{N} = 400 GeV, norm. to bkg.";
		}
		//loop over all kinematic distributions
		for(unsigned dist = 0; dist < nDist; ++dist){
			if(cat < 4){
				signals[0] = histos[cat][dist][2];
				signals[1] = histos[cat][dist][3];
				signals[2] = histos[cat][dist][5];
				signals[3] = histos[cat][dist][6];
				signals[4] = histos[cat][dist][8];
				signals[5] = histos[cat][dist][9];
			} else if(cat < 6){
				signals[0] = histos[cat][dist][11];
				signals[1] = histos[cat][dist][12];
				signals[2] = histos[cat][dist][13];
				signals[3] = histos[cat][dist][14];
				signals[4] = histos[cat][dist][15];
				signals[5] = histos[cat][dist][16];
			} else{
				signals[0] = histos[cat][dist][2];
				signals[1] = histos[cat][dist][5];
				signals[2] = histos[cat][dist][7];
				signals[3] = histos[cat][dist][11];
				signals[4] = histos[cat][dist][14];
				signals[5] = histos[cat][dist][15];
			}
			plotDataVSMC(dataHistos[cat][dist], bkgHistos[cat][dist], distNames, nSamples_eff - nSig, "searchR/" + histNames[dist] + "_" + catNames[cat] + "_withSignal" + extra, false, 0, "HNL", bkgSystDist[cat][dist], true, signals,  sigPlotNames ,nSigToPlot, true);
			plotDataVSMC(dataHistos[cat][dist], bkgHistos[cat][dist], distNames, nSamples_eff - nSig, "searchR/" + histNames[dist] + "_" +  catNames[cat]  +  extra, false, 0, "HNL", bkgSystDist[cat][dist]);
		}
	}
	storageFile->Close();
	*/
}



int main(int argc, char* argv[]){
	//TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;	
	testtree.hnlPlots();	
	cout << "Function ended!" << endl;
	//rootapp->Run();
    return 0;
}


