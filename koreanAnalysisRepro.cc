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
	const TString extra = "";	//for plot file names
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

	TH1D* yields[nSamples_eff + 1]; //Seperate histogram for every category
	const unsigned nSR = nSig; //testing several SR options
	for(unsigned effsam = 0; effsam < nSamples_eff + 1; ++effsam){
		yields[effsam] = new TH1D(eff_names[effsam], eff_names[effsam] + "; search region ; events/search region", nSR, 0.5, nSR + 0.5);
	}
	Double_t scale[nSamples -1];
	//set to check which data events have already been processed
	std::set<std::tuple<unsigned long, unsigned long, unsigned long> > usedEvents; //runNb, lumiBlock, eventNb	

	//Cuts defining each search region
	const double m4Cuts[nSig] = {125, 125, 125, 125, 125, 125, 125, 125, 80, 110, 160, 240, 340, 480, 560, 770, 890};
	const double pt1Cuts[nSig] = {60, 55, 50, 40, 35, 30, 30, 35, 25, 25, 45, 65, 120, 120, 100, 200, 290};
	const double pt2Cuts[nSig] = {45, 45, 40, 40, 35, 30, 25, 30, 40, 15, 40, 55, 65, 65, 110, 100, 190};
	const double pt3Cuts[nSig] = {25, 35, 40, 30, 30, 30, 25, 25, 15, 15, 20, 30, 35, 50, 50, 45, 80};
	const double metCuts[nSig] = {0, 0, 0, 0, 0, 0, 0, 0, 20, 20, 20, 20, 20, 20, 20, 30, 20};

	double m200Count = 0.;
	unsigned long unweightedNEvents = 0;
	//Loop over all samples
	for(unsigned sam = 0, effsam = 0; sam < nSamples; ++sam, ++effsam){
		if(sam != 0){
   			if(names[sam] == names[sam -1]) --effsam;
    	}		
    	Long64_t nEntries = inputTree[sam]->GetEntries();

		//if(fileList[sam] != "MajoranaNeutrinoToSSSF_MuMuE_M150.root" && fileList[sam] != "MajoranaNeutrinoToMuMuMu_M150.root") continue;	
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
			cutBased();
			//Change analysis object selection
			for(unsigned l = 0; l < _nL; ++l){
				_isloose[l] == _isloose[l] && (_lPt[l] > 10.);
				_isFO[l] == _isFO[l] && (_lPt[l] > 10.);
				_istight[l] == _istight[l] && (_lPt[l] > 10.);
			}
			//Baseline event selection
			if(!baseline(true, false, true, false)) continue;
			if(!vetoLowMll(4)) continue;
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
			//Select exactly three muons
			if(tril_flavorComb(ind, _flavors, lCount) != 3) continue;
			
			if(_flavors[ind[0]] != 1) continue;
			if(_flavors[ind[1]] != 1) continue;
			if(_flavors[ind[2]] != 1) continue;
			
			//Apply analysis Pt thresholds
			if(_lPt[ind[0]]*(1 + std::max(_isolation[ind[0]] - 0.1, 0.)) < 20.) continue;
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
			//Require 3 leptons to be tight in data and MC, and determine nonPrompt bkg in data
			unsigned nTight = tightCount(ind, lCount);
			bool tightFail = nTight < 3;
			//index used to fill events, needed to separate fakes from data
			unsigned fill = effsam;
			//Apply FR maps to data control region
			double* conePt = new double[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				//conePt[l] =  PtCone(_lPt[ind[l]], _flavors[ind[l]], _lepMVA[ind[l]], _ptratio[ind[l]]);
				conePt[l] = _lPt[ind[l]]*(1 + std::max(_isolation[ind[l]] - 0.1, 0.));
			}
			if(tightFail && (effsam == 0 || effsam > nSig)){ //&& effsam == 0){
				//fakes go in different histogram
				fill = nSamples_eff;
				//Apply FR maps
				if(effsam != 0) scal *= -1;
				scal*= fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap, lCount);
			} else if(tightFail) continue;
			//determine search category
			unsigned cat = hnl::cat(ind, _flavors, _charges, lCount, conePt[0]);
			if(cat == 999) continue;

			//Sample overlap removal
			if(conePt[2] >= 10){
				if(fileList[sam] == "TTJets_DiLept_trilepton.root" || fileList[sam] == "DYJetsToLL_M10to50_trilepton.root" || fileList[sam] == "DYJetsToLL_M50_trilepton.root"  || fileList[sam] == "TTJets_SingleLeptFromTbar_trilepton.root"  || fileList[sam] == "TTJets_SingleLeptFromT_trilepton.root" ) continue;
			} else if(conePt[2] < 10){
				if(fileList[sam] == "TTJets_DiLept_trilepton.root" || fileList[sam] == "DYJetsToLL_M10to50_trilepton.root" || fileList[sam] == "DYJetsToLL_M50_trilepton.root"  || fileList[sam] == "TTJets_SingleLeptFromTbar_trilepton.root"  || fileList[sam] == "TTJets_SingleLeptFromT_trilepton.root" ){
					double maxMll = 0.;
					for(unsigned l = 0; l < _gen_nL -1; ++l){
						TLorentzVector lep1;
						lep1.SetPtEtaPhiE(_gen_lPt[l], _gen_lEta[l], _gen_lPhi[l], _gen_lE[l]);
						for(unsigned k = l + 1; k < _gen_nL; ++k){
							TLorentzVector lep2;
							lep2.SetPtEtaPhiE(_gen_lPt[k], _gen_lEta[k], _gen_lPhi[k], _gen_lE[k]);
							if( (lep1 + lep2).M() > maxMll) maxMll = (lep1 + lep2).M();
						}
					}
					if(maxMll > 30){
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
				if(fileList[sam] == "TTGJets_trilepton.root" || fileList[sam] == "ZGTo2LG_trilepton.root"){
					bool promptfail =true;
					for(unsigned l = 0; l < lCount; ++l){	
						//cout << _origin[ind[l]] << endl;
						if(_pdgmc[ind[l]] == 22 && _originPhot[ind[l]] == 0){
							promptfail = false;
							break;
						}
					}
					if(promptfail) continue;
				}
			}	
			//Apply triggers to data events;
			if(sam == 0 || effsam > nSig){
				bool trigPass[4];
				trigPass[0] = _lowM_trigger_eee || _lowM_trigger_all;
				trigPass[1] = _lowM_trigger_mee || _lowM_trigger_all;
				trigPass[2] = _lowM_trigger_mme || _lowM_trigger_all;
				trigPass[3] = _lowM_trigger_mmm || _lowM_trigger_all;
				if(!trigPass[tril_flavorComb(ind, _flavors, lCount)]) continue;	
			}				
			//determine search category
			TLorentzVector* lepV = new TLorentzVector[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				lepV[l].SetPtEtaPhiE(conePt[l], _lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]*(conePt[l]/_lPt[ind[l]]) );
			}
			//Calculate lepton system vector
			TLorentzVector lepSyst;
			for(int l = 0; l < 3; ++l) lepSyst += lepV[l];
			

			unsigned mllI[2] = {99, 99};
			mllIndices(mllI, ind, lepV, _charges, _flavors, lCount);				
			//determine mll
			double mll;
			if(mllI[0] == 99){
				if(mllI[1] != 99) std::cerr << "error one mll index is not -1 while the other is" << endl;
				mll = -1;
			} else{
				TLorentzVector lzV[2];
				//for(unsigned l = 0; l < 2; ++l) lzV[l].SetPtEtaPhiE(PtCone(_lPt[mllI[l]], _flavors[mllI[l]], _lepMVA[mllI[l]], _ptratio[mllI[l]]), _lEta[mllI[l]], _lPhi[mllI[l]], _lE[mllI[l]]);
				for(unsigned l = 0; l < 2; ++l) lzV[l].SetPtEtaPhiE(_lPt[mllI[l]]*std::max(1., 1 + (_isolation[mllI[l]] - 0.1)), _lEta[mllI[l]], _lPhi[mllI[l]], _lE[mllI[l]]*std::max(1., 1 + (_isolation[mllI[l]] - 0.1)));
				mll = (lzV[0] + lzV[1]).M();
			}
			if( fabs(mll - 91) < 15) continue; //consider onZ events for WZ CR	
			if( fabs(lepSyst.M() - 91) < 15) continue;


			//Reconstruct neutrino Z pt


			//Calculate m(3l + nu)

			/*			
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
			*/

			//Calculate MET vector
			TLorentzVector METvec;
			METvec.SetPtEtaPhiE(_met, 0, _met_phi, _met);

			//Check lepton charge combinations
			unsigned l1 = 99, l2 = 99, l3 = 99;
			unsigned SSI[2] = {99, 99};
			unsigned OSI = 99;
			if(_charges[ind[0]] == _charges[ind[1]] && _charges[ind[1]] != _charges[ind[2]]){
				SSI[0] = ind[0];
				SSI[1] = ind[1];
				OSI = ind[2];
			} else if(_charges[ind[0]] == _charges[ind[2]] && _charges[ind[0]] != _charges[ind[1]]){
				SSI[0] = ind[0];
				SSI[1] = ind[2];
				OSI = ind[1];
			} else if(_charges[ind[1]] == _charges[ind[2]] && _charges[ind[0]] != _charges[ind[1]]){
				SSI[0] = ind[1];
				SSI[1] = ind[2];
				OSI = ind[0];
			} else continue;
			//find lepton giving MT closest to mw
			const double mw = 80.385;
			unsigned bestMwI = 99;
			unsigned wlepI = 99;
			double minMwDiff = 99999.;
			for(unsigned l = 0; l < lCount; ++l){
				if(fabs(transmass(lepV[l], METvec) - mw) < minMwDiff){
					minMwDiff = fabs(transmass(lepV[l], METvec) - mw);
					bestMwI = ind[l];
					wlepI = l;
				}
			}
			//Calculate M4 in korean way
			//Low mass version
			double A = mw*mw - lepSyst.M()*lepSyst.M() + 2.*(lepSyst.Px()*_met*cos(_met_phi) + lepSyst.Py()*_met*sin(_met_phi));
			double a = lepSyst.E()*lepSyst.E() - lepSyst.Pz()*lepSyst.Pz();
			double b = -A*lepSyst.Pz();
			double c = lepSyst.E()*lepSyst.E()*_met*_met - A*A/4.;
			double root = (b*b - 4.*a*c);
			if(root < 0) root = 0;
			double pzPlus = (-b + root)/(2.*a);
			double pzMin = (-b - root)/(2.*a);
			double pz;
			if(fabs(pzMin) < fabs(pzPlus) ){
				pz = pzMin;
			} else{
				pz = pzPlus;
			}
			TLorentzVector nu;
			nu.SetPxPyPzE(_met*cos(_met_phi), _met*sin(_met_phi), pz, sqrt(_met*_met + pz*pz) );
			double m4LowM = (lepSyst + nu).M();
			//high mass version
			A = mw*mw - lepV[wlepI].M()*lepV[wlepI].M() + 2*(lepV[wlepI].Px()*_met*cos(_met_phi) + lepV[wlepI].Py()*_met*sin(_met_phi));
			a = lepV[wlepI].E()*lepV[wlepI].E() - lepV[wlepI].Pz()*lepV[wlepI].Pz();
			b = -A*lepV[wlepI].Pz();
			c = lepV[wlepI].E()*lepV[wlepI].E()*_met*_met - A*A/4.;
			root = (b*b - 4.*a*c);
			if(root < 0) root = 0;
			pzPlus = (-b + root)/(2.*a);
			pzMin = (-b - root)/(2.*a);
			if(fabs(pzMin) < fabs(pzPlus) ){
				pz = pzMin;
			} else{
				pz = pzPlus;
			}
			nu.SetPxPyPzE(_met*cos(_met_phi), _met*sin(_met_phi), pz, sqrt(_met*_met + pz*pz) );
			double m4HighM = (lepSyst + nu).M();

			//Fill search regions
			for(unsigned sig = 0; sig < nSig; ++sig){
				//only fill one signal in every bin
				if(effsam != 0 && effsam < nSig){
					if(effsam != sig + 1) continue;
				}
				if(sig < 8){
					//if(lepSyst.M() > 125.) continue;
					//Determine indices of l1, l2 and l3
					if(sig < 6){
						l1 = SSI[1];
					} else{
						l1 = SSI[0];
					}
					for(unsigned l = 0; l < lCount; ++l){
						if(ind[l] == l1) continue;
						l2 = ind[l];
						break;
					}
					for(unsigned l = 0; l < lCount; ++l){
						if(ind[l] == l1 || ind[l] == l2) continue;
						l3 = ind[l];
					}
					/*
					if(_lPt[l1]*(1 + std::max(_isolation[l1] - 0.1, 0.)) >  pt1Cuts[sig]) continue;
					if(_lPt[l2]*(1 + std::max(_isolation[l2] - 0.1, 0.)) >  pt2Cuts[sig]) continue;
					if(_lPt[l3]*(1 + std::max(_isolation[l3] - 0.1, 0.)) >  pt3Cuts[sig]) continue;
					*/
					if(_lPt[l1] >  pt1Cuts[sig]) continue;
					if(_lPt[l2] >  pt2Cuts[sig]) continue;
					if(_lPt[l3] >  pt3Cuts[sig]) continue;
					
					/*
					if(_lPt[ind[0]]*(1 + std::max(_isolation[ind[0]] - 0.1, 0.)) >  pt1Cuts[sig]) continue;
					if(_lPt[ind[1]]*(1 + std::max(_isolation[ind[1]] - 0.1, 0.)) >  pt2Cuts[sig]) continue;
					if(_lPt[ind[2]]*(1 + std::max(_isolation[ind[2]] - 0.1, 0.)) >  pt3Cuts[sig]) continue;
					*/
					if(m4LowM > m4Cuts[sig]) continue;	
				} else{
					l3 = bestMwI;
					if(l3 = SSI[0]){
						l1 = SSI[1];
					} else if(l3 = SSI[1]){
						l1 = SSI[0];
					} else{
						l1 = SSI[0];
					} 
					for(unsigned l = 0; l < lCount; ++l){
						if(ind[l] == l1 || ind[l] == l3) continue;
						l2 = ind[l];
					}
					/*
					if(_lPt[l1]*(1 + std::max(_isolation[l1] - 0.1, 0.)) <  pt1Cuts[sig]) continue;
					if(_lPt[l2]*(1 + std::max(_isolation[l2] - 0.1, 0.)) <  pt2Cuts[sig]) continue;
					if(_lPt[l3]*(1 + std::max(_isolation[l3] - 0.1, 0.)) <  pt3Cuts[sig]) continue;
					*/
					/*
					if(_lPt[ind[0]]*(1 + std::max(_isolation[ind[0]] - 0.1, 0.)) <  pt1Cuts[sig]) continue;
					if(_lPt[ind[1]]*(1 + std::max(_isolation[ind[1]] - 0.1, 0.)) <  pt2Cuts[sig]) continue;
					if(_lPt[ind[2]]*(1 + std::max(_isolation[ind[2]] - 0.1, 0.)) <  pt3Cuts[sig]) continue;
					*/
					if(_met < metCuts[sig]) continue;
					if(m4HighM < m4Cuts[sig]) continue;
					//Determine indices of l1, l2 and l3
				}
				yields[fill]->Fill(sig + 1, scal);
			}
				if(fileList[sam] == "MajoranaNeutrinoToSSSF_MuMuE_M150.root" || fileList[sam] == "MajoranaNeutrinoToMuMuMu_M150.root"){
					m200Count = m200Count + scal;
					++unweightedNEvents;
				}
			//}
    	}
	}
	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << "expected 150 GeV yields = " << m200Count << endl;
	cout << "unweighted 150 GeV yields = " << unweightedNEvents << endl;
	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;

	//CHANGE ALL NUMBERS BACK
	//Split data and MC histograms for plotting and propagating uncertainties
	TH1D* dataYields;
	dataYields = (TH1D*) yields[nSig + 1]->Clone();
		
	TH1D* bkgYields[nSamples_eff -nSig]; //change to nSamples_eff if sig is removed
	TH1D* bkgErros[nSamples_eff - nSig];
	for(unsigned effsam = nSig + 1; effsam < nSamples_eff + 1; ++effsam){
		bkgYields[effsam -nSig - 1] = (TH1D*) yields[effsam]->Clone();
		if(effsam > nSig + 1){
			dataYields->Add(bkgYields[effsam - nSig - 1]);
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//PRINT DATA CARDS FOR LIMIT CALCULATION
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//const unsigned nSig = 13;
	const unsigned nBkg = nSamples_eff - nSig; //CHANGED
	const unsigned nSyst = 5 + 1 + 2*nBkg; //5 general uncertainties, one stat unc for signal, stat + extra unc for every bkg;
	const TString sigN[nSig] = {"m5", "m10", "m20", "m30", "m40", "m50", "m60", "m70", "m90", "m100", "m150", "m200", "m300", "m400", "m500", "m700", "m1000"};
	const TString bkgNames[nBkg] = {"ZZH", "triboson", "WZ", "Xamma", "TTX",  "nonPrompt"}; //slightly rewrite bkg names not to confuse the combine tool
	std::vector<std::vector<double>> systUnc(nSyst, std::vector<double>(nBkg + 1, 0)); //2D array containing all systematics
	//initialize values non-bin-dependent systematics
	for(unsigned proc = 0; proc < nBkg; ++proc){ //last background is non-prompt and isn't susceptible to these uncertainty sources (else nBkg + 1 in loop_
		systUnc[0][proc] = 1.026;	//lumi
		systUnc[1][proc] = 1.04;	//id eff
		systUnc[2][proc] = 1.02;	//trig eff
		systUnc[3][proc] = 1.05;	//JEC	
		systUnc[4][proc] = 1.05;	//PU
	}
	for(unsigned syst = 0; syst < 4; ++syst) systUnc[syst][nBkg] = 0; //assign 0 for non-prompt
	//names of uncertainties;
	TString systNames[nSyst] = {"lumi", "id_eff", "trigeff", "JEC", "PU"};
	TString systDist[nSyst];
	for(unsigned syst = 0; syst < nSyst; ++syst) systDist[syst] = "lnN";
	const double extraUnc[nBkg] = {1.25, 1.5, 1.10, 1.15, 1.15, 1.3}; //extra unc assigned to backgrounds
	for(unsigned syst = 6 + nBkg; syst < nSyst; ++syst){
		unsigned bkg = syst - 6 - nBkg;
		systNames[syst] = "extra" + bkgNames[bkg];
		systUnc[syst][bkg + 1] = extraUnc[bkg];
	}
	for(unsigned sig = 0; sig < nSig; ++sig){
		for(unsigned bin = 1; bin < dataYields->GetNbinsX() + 1; ++bin){
			if(bin != (sig + 1)) continue;
			//signal stat unc
			systUnc[5][0] = 1 + std::max(0., yields[sig + 1]->GetBinError(bin)/yields[sig + 1]->GetBinContent(bin));
			systNames[5] = "statSig" + std::to_string(bin);
			//bkg stat uncertainties
			double bkgYieldVal[nBkg];
			for(unsigned bkg = 0; bkg < nBkg; ++bkg){
				bkgYieldVal[bkg] = bkgYields[bkg]->GetBinContent(bin);
				systUnc[6 + bkg][bkg + 1]  = 1 + std::max(0., bkgYields[bkg]->GetBinError(bin)/bkgYields[bkg]->GetBinContent(bin));
				systNames[6 + bkg] = "stat"  + bkgNames[bkg] + std::to_string(bin);
			} 
			hnl::printDataCard(dataYields->GetBinContent(bin), yields[sig + 1]->GetBinContent(bin), sigN[sig], bkgYieldVal, nBkg, bkgNames, systUnc, nSyst, systNames, systDist, "datacards/datacard" + sigN[sig] + "_korean_bin" + (TString) std::to_string(bin) + extra);
		}
	}
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	TString sigNames[nSig] = {"m_{N} = 5 GeV", "m_{N} = 10 GeV", "m_{N} = 20 GeV", "m_{N} = 30 GeV", "m_{N} = 40 GeV", "m_{N} = 50 GeV", "m_{N} = 60 GeV", "m_{N} = 70 GeV", "m_{N} = 90 GeV", "m_{N} = 100 GeV", "m_{N} = 150 GeV", "m_{N} = 200 GeV", "m_{N} = 300 GeV", "m_{N} = 400 GeV", "m_{N} = 500 GeV", "m_{N} = 700 GeV", "m_{N} = 1000 GeV"};
	for(unsigned sig = 0; sig < nSig; ++sig){
		sigNames[sig] += ", |V_{eff}|^{2} = ";
	}
	for(unsigned sig = 0; sig < 7; ++sig){
		sigNames[sig] += "10^{-5}";
	}
	sigNames[7] += "10^{-4}";
	sigNames[8] += "10^{-3}";
	for(unsigned sig = 9; sig < nSig; ++sig){
		sigNames[sig] += "10^{-2}";
	}
	const TString distNames[nSamples_eff + 1 - nSig] = {"total BKG", "ZZ/H", "triboson", "WZ", "X + #gamma", "TT/T + X",  "non-prompt"};
	
	//Plot the yields as a function of the search region
	unsigned nSigToPlot = 7;
	TH1D* signals[nSigToPlot];
	TString sigPlotNames[nSigToPlot];
	for(unsigned sig = nSig - nSigToPlot; sig < nSig; ++sig){
		cout << "sig = " << sig << endl;
		signals[sig- nSig + nSigToPlot] = (TH1D*) yields[sig + 1]->Clone();
		sigPlotNames[sig- nSig + nSigToPlot] = sigNames[sig];
	}	
	plotDataVSMC(dataYields, bkgYields, distNames, nSamples_eff - nSig, + "_koreanAnalysis_withsyst" + extra, true, 0, "HNL", true, signals,  sigPlotNames ,nSigToPlot);
}



int main(int argc, char* argv[]){
	TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;
	testtree.Loop();
	rootapp->Run();
    return 0;
}


