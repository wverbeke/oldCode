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


void trilTree::Loop(){
	//Set plotting style
	setTDRStyle();
	//gROOT->SetBatch(kTRUE);
	//Define list of samples
	const unsigned nSamples = 3;
	const TString fileList[nSamples] = {"WZTo3LNu.root", "WGToLNuG.root", "WJetsToLNu.root"};
	//489 350.674
	const double xSections[nSamples] = {4.42965, 489, 61526.7};
	const TString names[nSamples] = {"WZ", "W#gamma", "WJets"};

	//Read Trees from ROOT files
	TFile* hfile[nSamples];
	TTree* inputTree[nSamples];
	double hcounter[nSamples];

	for(unsigned sam = 0; sam < nSamples; ++sam){   //CHANGE BACK TO NSAMPLES
		cout << "name " << names[sam] << endl;
		hfile[sam] = new TFile("../data_EWKmoriond/"+fileList[sam],"read");
       	hfile[sam]->cd("FakeElectrons");
		//Determine hcounter for cross section scaling
		TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
		_hCounter->Read("hCounter");
		hcounter[sam] = _hCounter->GetBinContent(1);
       	inputTree[sam] = (TTree*) (hfile[sam]->Get("FakeElectrons/fakeTree"));
		bool photTree = (sam > 0);
		Init(inputTree[sam], photTree, true);
	}
	
	//Tweakable options////////////////////////////////////////////////////
	const bool TestRun = false;	//Break after a few events
	const double DataLuminosity = 35.867;    //units of fb^{-1}
	const bool ptReweighing = false;
	const bool onZ = false;
	const bool flavorDiff = true;	
	const bool pixVeto = true;
	const bool eVeto = true;
	const bool onlyMu = true;
	const double mtMin = 40;
	const double mtMax = 200;
	const double photonPtCut = 40;
	const double deltaRCut = 0.3;
	const double deltaPhiCut = 0;
	const TString extra = "_WgammaVsWjets";	//for plot file names
	//////////////////////////
	
	//Make histograms containing kinematic distributions
	const unsigned nDist = 9;
	const TString distNames[nDist] = {"mt", "met", "metPhi","lepPt", "bosonPt", "DeltaPhilepmet", "DeltaRlepgamma", "DeltaPhilepgamma", "ht"};
	const TString xAxes[nDist] = {"M_{T}", "MET", "#Phi(MET)", "P_{T}(lepton)", "P_{T}(Z/#gamma)", "#Delta#Phi(lepton, MET)", "#DeltaR(lepton, Z/#gamma)", "#Delta#Phi(lepton, Z/#gamma)", "H_{T}"};
	const TString units[nDist] = {"GeV", "GeV", "", "GeV", "GeV", "", "", "",  "GeV"};
	const double histMin[nDist] = {mtMin, 50, 0, 20, photonPtCut, 0, deltaRCut, deltaPhiCut, 30};
	const double histMax[nDist] = {mtMax, 200, 3.2, 150, 300, 3.2, 7, 3.2, 600};
	const int nBins[nDist] = {20, 20, 20, 20, 20, 20, 20, 20, 20};
	TH1D* hists[nDist][nSamples];
	for(unsigned dist = 0; dist < nDist; ++dist){
		float binWidth = (histMax[dist] - histMin[dist])/nBins[dist];
		std::ostringstream strs; strs << binWidth; std::string yAxis = strs.str();
		for(unsigned sam = 0; sam < nSamples; ++sam){
			hists[dist][sam] = new TH1D(distNames[dist] + names[sam], distNames[dist] + names[sam] + ";" + xAxes[dist] + " (" + units[dist] +  "); events /" + yAxis + units[dist], nBins[dist], histMin[dist], histMax[dist]);
		}
	}
	
	TH1D* ptWeights;
	
	double maxBinC[nDist];
	for(unsigned dist = 0; dist < nDist; ++dist){
		maxBinC[dist] = hists[dist][0]->GetBinCenter(hists[dist][0]->GetNbinsX());
	}

	
	TH2D *MTMET[nSamples];
	double MTbins[4] = {40, 100,160, 220};
	double METbins[6] = {50, 100, 150, 200, 250};
	for(int i = 0; i < nSamples; ++i){
		MTMET[i] = new TH2D("MTMET" + names[i], "MTMET" + names[i] + ";" + xAxes[0] + ";" + xAxes[1], 3, MTbins , 4, METbins);
		MTMET[i]->Sumw2();
	}

    Double_t scale[nSamples];
	//Loop over all samples
	for(unsigned re = 0; re < 2; ++re){
		if(!ptReweighing && re == 0) continue;
		for(unsigned sam = 0; sam < 3; ++sam){
			if(ptReweighing && re == 1 && sam == 0) continue;
			Long64_t nEntries = inputTree[sam]->GetEntries();
			scale[sam] = xSections[sam]*DataLuminosity*1000/(hcounter[sam]);
			std::cout<<"Entries in "<< fileList[sam] <<" "<<nEntries<<std::endl;
		    for(Long64_t it = 0; it < nEntries; ++it){
				if (it%10000 == 0) cout<<'.'<<flush;
		    	inputTree[sam]->GetEntry(it);
		    	if(TestRun && it > 10000) break;
		    	double scal;
		    	scal = scale[sam]*_weight;
				//Baseline event selection:
				if(!baseline(!sam, false, true, false)) continue;
				//Categorize according to the number of leptons and flavors
				unsigned* ind = new unsigned[_nL];
				unsigned lCount;
				unsigned lw = 99;
				TLorentzVector boson;
				//WZ event selection
				if(_met < 50) continue;
				if(sam == 0){
					if(!(lCount = lepOrder(ind, 3)) ) continue;
					if(lCount != 3) continue;
					//if(!ptCuts(ind, lCount) ) continue;
					unsigned mllI[2] = {99, 99};
					TLorentzVector* lepV = new TLorentzVector[lCount];
					for(unsigned l = 0; l < lCount; ++l){
						lepV[l].SetPtEtaPhiE(_lPt[ind[l]], _lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]);
					}
					if(!wzSel(ind, mllI, lw, lCount, lepV, flavorDiff, onZ) )continue;
					if(onlyMu && _flavors[lw] != 1) continue;
					TLorentzVector lzV[2];
					for(unsigned l = 0; l < 2; ++l) lzV[l].SetPtEtaPhiE(_lPt[mllI[l]], _lEta[mllI[l]], _lPhi[mllI[l]], _lE[mllI[l]]);
					boson = lzV[0] + lzV[1];
				//Wgamma and WJets selection
				} else{
					if(names[sam] == "WJets"){
						bool promptfail = false;
						for(unsigned ph = 0; ph < _gen_nPh; ++ph){
							if(_gen_lPt[ph] > 10 && fabs(_gen_phmompdg[ph] != 111)){
							//if(_gen_lPt[ph] > 10 && fabs(_gen_isPromptPh[ph])){
								promptfail = true;
								break;
							}
						}
						if(promptfail) continue;
					}
					unsigned ph = 99;
					if(!(lCount = lepOrder(ind, 1)) ) continue;
					//Does it work better with single lepton veto
					if(lCount != 1) continue;
					if(tightCount(ind, lCount) != 1) continue;
					if(!wgSel(ph, ind, lCount, photonPtCut, deltaRCut, deltaPhiCut, pixVeto, eVeto, onlyMu) ) continue;
					lw = ind[0];
					boson.SetPtEtaPhiE(_phPt[ph], _phEta[ph], _phPhi[ph], _phE[ph]);			
				}
				TLorentzVector lep, metV;
				lep.SetPtEtaPhiE(_lPt[lw], _lEta[lw], _lPhi[lw], _lE[lw]);
				metV.SetPtEtaPhiE(_met, 0, _met_phi, _met);
				if(transmass(lep,metV) < mtMin) continue;
				double fill[nDist] = {transmass(lep,metV),  _met, _met_phi, _lPt[lw], boson.Pt() , fabs(lep.DeltaPhi(metV)), lep.DeltaR(boson), fabs(lep.DeltaPhi(boson)), _HT};
				if(ptReweighing && re == 1 && sam == 1){
					scal*= ptWeights->GetBinContent(ptWeights->FindBin(std::min(_lPt[lw], ptWeights->GetBinCenter(ptWeights->GetNbinsX()) )));
					if(ptWeights->GetBinContent(ptWeights->FindBin(std::min(_lPt[lw], ptWeights->GetBinCenter(ptWeights->GetNbinsX()) ))) > 10){
						cout  << "~~~~~~~~~~~~~~~~~~~~~~~" << endl;
						cout  << ptWeights->GetBinContent(ptWeights->FindBin(TMath::Min(_lPt[lw], ptWeights->GetBinCenter(ptWeights->GetNbinsX()) ))) << endl;
						cout  << "WZ Ptbin : " << hists[3][0]->FindBin(TMath::Min(_lPt[lw], hists[3][0]->GetBinCenter(hists[3][0]->GetNbinsX()) )) << endl;
						cout  << "WG Ptbin : " << hists[3][1]->FindBin(TMath::Min(_lPt[lw], hists[3][1]->GetBinCenter(hists[3][1]->GetNbinsX()) )) << endl;
					}
				}
				//sam > 0
				MTMET[sam > 0]->Fill(TMath::Min(transmass(lep,metV),219.), TMath::Min(_met,249.), scal);
				for(unsigned dist = 0; dist < nDist; ++dist){
					hists[dist][sam]->Fill(std::min(fill[dist], maxBinC[dist]), scal);
				}			
			}
		}
		if(ptReweighing && re == 0){
			cout << "hists[3][0]->GetSumOfWeights() = " << hists[3][0]->GetSumOfWeights() << endl;
			cout << "hists[3][1]->GetSumOfWeights() = " << hists[3][1]->GetSumOfWeights() << endl;
			hists[3][0]->Scale(1/hists[3][0]->GetSumOfWeights());
			hists[3][1]->Scale(1/hists[3][1]->GetSumOfWeights());
			ptWeights = (TH1D*) hists[3][0]->Clone();
			ptWeights->Divide(hists[3][1]);
			//ptWeights = HistDiv(hists[3][0], hists[3][1]);
			cout << "ptWeights->GetSumOfWeights() : " <<  ptWeights->GetSumOfWeights() << endl;
			for(unsigned h = 1; h < hists[3][1]->GetNbinsX() + 1; ++h){
				cout <<"bin(" << h <<") Content = " << hists[3][1]->GetBinContent(h) << endl;
			}
		}	
	}
	plotHistRatio(hists[0][0], hists[0][1], "WZ", "W#gamma", "MTcomp" + extra, false, mtMin,mtMax,true);
	plotHistRatio(hists[1][0], hists[1][1], "WZ", "W#gamma", "METcomp" + extra);
	plotHistRatio(hists[2][0], hists[2][1], "WZ", "W#gamma", "METPhicomp" + extra);
	plotHistRatio(hists[3][0], hists[3][1], "WZ", "W#gamma", "leptonPtcomp" + extra);

	plotHistRatio(hists[0][1], hists[0][2], "W#gamma", "WJets", "WgammaVSWjetsMT", false, 0, 0, true);
	plotHistRatio(hists[1][1], hists[1][2], "W#gamma", "WJets", "WgammaVSWjetsMET", false, 0, 0, false);

	TFile *MC_SF = new TFile("../weights/WZoverWgammaMC_SF.root", "recreate");
	MTMET[1]->Scale(MTMET[0]->Integral()/MTMET[1]->Integral());
	MTMET[0]->Write("WZ_MTvsMET");
	MTMET[0]->Divide(MTMET[1]);
	MTMET[0]->Write("WZoverWgamma_MTvsMET");
	plotHist(MTMET[0], "MTMET_WZoverWgamma");
	//WZ and WZ/Wgamma MT distr
	hists[0][0]->Write("WZ_MT");
	hists[0][1]->Scale(hists[0][0]->Integral()/hists[0][1]->Integral());
	hists[0][0]->Divide(hists[0][1]);
	hists[0][0]->Write("WZoverWgamma_MT");
	MC_SF->Close();
}



int main(int argc, char* argv[]){
	TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;
	testtree.Loop();
	rootapp->Run();
    return 0;
}


