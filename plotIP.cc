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
	//Define list of samples
	const unsigned nSamples = 2;
	const TString fileList[nSamples] = {"HeavyNeutrino_trilepton_M-5_V-0.00316_mu_displaced_unskimmedLep.root", "HeavyNeutrino_trilepton_M-6_V-0.00316_mu_prompt_unskimmedLep.root"};
	const double xSections[nSamples] = {1, 2};
	const TString names[nSamples] = {"m5", "m6"};

	//Read Trees from ROOT files
	TFile* hfile[nSamples];
	TTree* inputTree[nSamples];
	double hcounter[nSamples];

	for(unsigned sam = 0; sam < nSamples; ++sam){   //CHANGE BACK TO NSAMPLES
		cout << "name " << names[sam] << endl;
		//hfile[sam] = new TFile("../../public/ntuples/"+fileList[sam],"read");
		hfile[sam] = new TFile("../unskimmedSignal/"+fileList[sam],"read");
       	hfile[sam]->cd("FakeElectrons");
		//Determine hcounter for cross section scaling
		TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
		_hCounter->Read("hCounter");
		hcounter[sam] = _hCounter->GetBinContent(1);
       	inputTree[sam] = (TTree*) (hfile[sam]->Get("FakeElectrons/fakeTree"));
		bool photTree = (sam > 0);
		Init(inputTree[sam], true, true);
	}
	
	//Tweakable options////////////////////////////////////////////////////
	const bool TestRun = false;	//Break after a few events
	const double DataLuminosity = 35.9;    //units of fb^{-1}
	const TString extra = "";	//for plot file names
	//////////////////////////
	
	//Make histograms containing kinematic distributions
	const unsigned nDist = 3;
	const TString distNames[nDist] = {"dxy", "dz", "SIP"};
	const TString xAxes[nDist] = {"|d_{xy}| (cm)", "|d_{z}| (cm)", "SIP_{3D}"};
	const TString units[nDist] = {"", "", ""};
	const double histMin[nDist] = {0, 0, 0};
	const double histMax[nDist] = {1, 1, 20};
	const int nBins[nDist] = {100, 100, 100};
	TH1D* hists[nDist][nSamples];
	for(unsigned dist = 0; dist < nDist; ++dist){
		float binWidth = (histMax[dist] - histMin[dist])/nBins[dist];
		std::ostringstream strs; strs << binWidth; std::string yAxis = strs.str();
		for(unsigned sam = 0; sam < nSamples; ++sam){
			hists[dist][sam] = new TH1D(distNames[dist] + names[sam], distNames[dist] + names[sam] + ";" + xAxes[dist] + " (" + units[dist] +  "); normalized events /" + yAxis + units[dist], nBins[dist], histMin[dist], histMax[dist]);
		}
	}
	
	double maxBinC[nDist];
	for(unsigned dist = 0; dist < nDist; ++dist){
		maxBinC[dist] = hists[dist][0]->GetBinCenter(hists[dist][0]->GetNbinsX());
	}
	

    Double_t scale[nSamples];
	//Loop over all samples
	for(unsigned sam = 0; sam < nSamples; ++sam){
		Long64_t nEntries = inputTree[sam]->GetEntries();
		scale[sam] = xSections[sam]*DataLuminosity*1000/(hcounter[sam]);
		std::cout<<"Entries in "<< fileList[sam] <<" "<<nEntries<<std::endl;
	    for(Long64_t it = 0; it < nEntries; ++it){
			inputTree[sam]->GetEntry(it);
        	if (it%10000 == 0) cout<<'.'<<flush;
        	if(TestRun && it > 10000) break;
			double scal = scale[sam]*_weight;
						
			for(unsigned l = 0; l < _nL; ++l){
				if(_mompdg[l] != 9900012) continue;
				hists[0][sam]->Fill(std::min(_ipPV[l],  histMax[0] - 0.000000001), scal);
				hists[1][sam]->Fill(std::min(_ipZPV[l],  histMax[1]- 0.000000001), scal);
				hists[2][sam]->Fill(std::min(_3dIPsig[l],  histMax[2]- 0.000000001), scal);
			}
		}
		
	}
	//plotHistRatio(hists[0][0], hists[0][1], "new", "old", "mllCompWZ" + extra,  false, 0,0, true, true);
	//std::cout << "ratio WZ_minmll01/WZ_old = " << hists[0][0]->GetSumOfWeights()/hists[0][1]->GetSumOfWeights() << std::endl;

	//, SF = 0.7
	//plotHistRatio(hists[0][0], hists[1][0], "mt old", "mt new", "MT_newalgo_neutrinoZ" + extra, false, 0,0, true, true);
	
	//plotHistRatio(hists[1][0], hists[1][1], "mllmin01, SF = 0.7", "old powheg", "Mmumu_WZsampcomp" + extra,  false, 0,0, false, true);
	//plotHistRatio(hists[1][0], hists[1][1], "mllmin01, SF = 0.7", "old powheg", "Mmumu_WZsampcomp_log" + extra,  false, 0,0, true, true);
	
	//plotHist(hists[2][0], "lowGenMll");
	std::vector<TH1D*> distVec;
	std::vector<TString> distnames = {"5 GeV, c#tau = 5.7 mm", "6 GeV, prompt"};
	for(unsigned i = 0; i < 2; ++i){
		distVec.push_back(hists[0][i]);
	}
	plotHist(distVec, distnames, "signaldxyplot", true);
	distVec.clear();
	for(unsigned i = 0; i < 2; ++i){
		distVec.push_back(hists[1][i]);
	}
	plotHist(distVec, distnames, "signaldzplot", true);
	distVec.clear();
	for(unsigned i = 0; i < 2; ++i){
		distVec.push_back(hists[2][i]);
	}
	plotHist(distVec, distnames, "signalSIPplot", true);
}



int main(int argc, char* argv[]){
	TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;
	testtree.Loop();
	rootapp->Run();
    return 0;
}


