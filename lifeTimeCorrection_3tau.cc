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
#include "TGraphErrors.h"

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
//#include "drawLumi.h"
const double tau_x = 0.82046213139960145;
const double tau_p = 5.6857483814612042*2;

std::vector<double> trilTree::findLifeTime(const double left, const double right){
	//5 GeV
	//Initialize all samples and cross sections
	const unsigned nSamples = 1;	
	const TString fileList[nSamples] = {"HeavyNeutrino_M5_V00p00336_mu.root"};

	//Read Trees from ROOT files
	TFile* hfile[nSamples];
	TTree* inputTree[nSamples];
	double hcounter[nSamples];
	for(unsigned sam = 0; sam < nSamples; ++sam){
		hfile[sam] = new TFile("../unskimmedSignal/"+fileList[sam],"read");
       	hfile[sam]->cd("FakeElectrons");
		//Determine hcounter for cross section scaling
		TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
		_hCounter->Read("hCounter");
		hcounter[sam] = _hCounter->GetBinContent(1);
       	inputTree[sam] = (TTree*) (hfile[sam]->Get("FakeElectrons/fakeTree"));
		Init(inputTree[sam], false, true);
	}

	//Tweakable options////////////////////////////////////////////////////
	const bool TestRun = false;	//Break after a few events
	const TString extra = "";
	//////////////////////////

	TH1D* eff_left[2];
	TH1D* eff_right[2];
	eff_left[0] = new TH1D("efficiency num left", "", 1, 0, 1);
	eff_left[1] = new TH1D("efficiency denom left", "", 1, 0, 1);
	eff_left[0]->Sumw2();
	eff_left[1]->Sumw2();

	eff_right[0] = new TH1D("efficiency num right", "", 1, 0, 1);
	eff_right[1] = new TH1D("efficiency denom right", "", 1, 0, 1);
	eff_right[0]->Sumw2();
	eff_right[1]->Sumw2();
	
	//Loop over all samples
	for(unsigned sam = 0; sam < nSamples; ++sam){
    	Long64_t nEntries = inputTree[sam]->GetEntries();
    	std::cout<<"Entries in "<< fileList[sam] <<" "<<nEntries<<std::endl;

        for(Long64_t it = 0; it < nEntries; ++it){
        	inputTree[sam]->GetEntry(it);
        	if (it%10000 == 0) cout<<'.'<<flush;
        	if(TestRun && it > 10000) break;

			cutBased();
			//Order particles by gen lepton Pt and then store them
			if(_gen_nL < 3) continue;
			//Correct weight for lifetime
			double leftWeight = left*left*exp(-_ctau/tau_x + _ctau/(left*sqrt(tau_x*tau_p) ));
			double rightWeight =  right*right*exp(-_ctau/tau_x + _ctau/(left*sqrt(tau_x*tau_p) ));
			
			eff_left[1]->Fill(sam +0.1,_weight);
			eff_right[1]->Fill(sam +0.1,_weight);	

			//Apply HNL SELECTION
			//Baseline event selection:
			if(_nL < 3) continue;
			if(!baseline(true, true, false, false)) continue;
			if(nBJets(true, false, 0) != 0) continue;
			//uncleaned bveto;	
			//Categorize according to the number of leptons and flavors
			unsigned* ind = new unsigned[_nL];
			unsigned lCount = lepOrder(ind, 3, true, true);	
			if(lCount != 3) continue; //Veto 4th FO lepton considering signal model!
			//Require 3 leptons to be tight in data and MC, and determine nonPrompt bkg in data
			unsigned nTight = tightCount(ind, lCount);
			bool tightFail = nTight < 3;
	
			//Calculate conePt for every lepton
			double* conePt = new double[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				conePt[l] = _lPt[ind[l]]*std::max(1., 1 + (_isolation[ind[l]] - 0.1));
			}
			//require 3 tight leptons
			if(tightFail) continue;
			//HNL pt Cuts
			if(!ptCuts_hnl(ind,lCount)) continue;
			//veto events outside the HNL selection categories (+++, ---)
			unsigned cat = hnl::cat(ind, _flavors, _charges, lCount, conePt[0]);
			if(cat == 999) continue;
			eff_left[0]->Fill(sam +0.1,_weight*leftWeight);
			eff_right[0]->Fill(sam +0.1,_weight*rightWeight);
    	}
	}
	eff_left[0]->Divide(eff_left[1]);
	std::cout << "eff_left[0]->GetBinError(1)/eff_left[0]->GetBinContent(1) = " << eff_left[0]->GetBinError(1)/eff_left[0]->GetBinContent(1) << std::endl;
	eff_right[0]->Divide(eff_right	[1]);
	std::vector<double> acc = {eff_left[0]->GetBinContent(1), eff_right[0]->GetBinContent(1)};
	return acc;
}



int main(int argc, char* argv[]){
	trilTree testtree;
	//displaced closure 
	const double promptAcc =  0.000621369;
	const double corrRatio = 0.56466;
	
	double left = 1;
	double right = corrRatio;
	for(unsigned i = 0; i < 15; ++i){
		double left_old = left;
		double right_old = right;
		left = left - (left - right)/2.;
		right = right + (left - right)/2.;
		std::vector<double> test = testtree.findLifeTime(left, right);
		int j = 0;
		while(test[0] < corrRatio*promptAcc*corrRatio){
			if(j == 0) right = left; //NEW for faster convergence
			left = left + (left_old - left)/(2.);
			test = testtree.findLifeTime(left, right);
			++j;
		}
		j = 0;
		while(test[1] > corrRatio*promptAcc*corrRatio){
			if(j == 0) left = right; //NEW for faster convergence
			right = right - (right - right_old)/(2.);
			test = testtree.findLifeTime(left, right);
			++j;
		}
		std::cout << "left = " << (left*sqrt(tau_p * tau_x))/tau_x << std::endl;
		std::cout << "right = " << (right*sqrt(tau_p * tau_x))/tau_x << std::endl;
	}
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
	std::cout << "Final result: left = " << (left*sqrt(tau_p * tau_x))/tau_x << " 		right = " << (right*sqrt(tau_p * tau_x))/tau_x << std::endl;
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    return 0;
}


