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
	//0.696758
	const TString fileList[nSamples] = {"WZTo3LNu_mllmin01.root", "WZTo3LNu_powheg.root"};
	//const TString fileList[nSamples] = {"WZTo3LNu_mllmin01.root", "WZTo3LNu.root"};
	//const TString fileList[nSamples] = {"WZTo3LNu_mllmin01.root", "test.root"};
	const double xSections[nSamples] = {58.59*0.652, 4.42965};
	//const double xSections[nSamples] = {58.59, 4.693};
	const TString names[nSamples] = {"WZTo3LNu_mllmin01", "WZTo3LNu"};

	//Read Trees from ROOT files
	TFile* hfile[nSamples];
	TTree* inputTree[nSamples];
	double hcounter[nSamples];

	for(unsigned sam = 0; sam < nSamples; ++sam){   //CHANGE BACK TO NSAMPLES
		cout << "name " << names[sam] << endl;
		//hfile[sam] = new TFile("../../public/ntuples/"+fileList[sam],"read");
		hfile[sam] = new TFile("../data_april17/"+fileList[sam],"read");
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
	const unsigned nDist = 1;
	const TString distNames[nDist] = {"mll"};
	const TString xAxes[nDist] = {"M_{#mu#mu}"};
	const TString units[nDist] = {"GeV"};
	const double histMin[nDist] = {10};
	const double histMax[nDist] = {200};
	const int nBins[nDist] = {200};
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
						
			if(_gen_nL < 3) continue;
			unsigned genLCount = 0;
			unsigned* genInd = new unsigned[3];
			for(unsigned l = 0; l < _gen_nL; ++l){
				if(!_gen_isPromptl[l]) continue;
				if(genLCount < 3) genInd[genLCount] = l;
				++genLCount;
			}
			if(genLCount != 3) continue;	
			/*
			if(tril_flavorComb(genInd, _gen_flavors, genLCount) != 2) continue; //mme

			TLorentzVector* lepV = new TLorentzVector[genLCount];
			for(unsigned l = 0; l < genLCount; ++l){
				lepV[l].SetPtEtaPhiE(_gen_lPt[genInd[l]], _lEta[genInd[l]], _lPhi[genInd[l]], _lE[genInd[l]]);
			}
			unsigned mllI[2] = {99, 99};
			mllIndices(mllI, genInd, lepV, _gen_charges, _gen_flavors, genLCount);				
			//determine mll
			double mll;
			if(mllI[0] == 99){
				if(mllI[1] != 99) std::cerr << "error one mll index is not -1 while the other is" << endl;
				mll = -1;
			} else{
				TLorentzVector lzV[2];
				for(unsigned l = 0; l < 2; ++l) lzV[l].SetPtEtaPhiE(_lPt[mllI[l]], _lEta[mllI[l]], _lPhi[mllI[l]], _lE[mllI[l]]);
				mll = (lzV[0] + lzV[1]).M();
			}
			*/
			unsigned muI[2];
			unsigned eI;
			unsigned muCount = 0;
			unsigned eCount = 0;
			for(unsigned l = 0; l < genLCount; ++l){
				if(_gen_flavors[l] == 0){
					++eCount;
					eI = l;
				} else if(_gen_flavors[l] == 1){ 
					if(muCount < 2) muI[muCount] = l;
					++muCount;
				}
			}
			if(muCount != 2) continue;
			if(eCount !=1) continue;
			if(_gen_charges[muI[0]] == _gen_charges[muI[1]]) continue;
			TLorentzVector l1,l2;
			l1.SetPtEtaPhiE(_gen_lPt[muI[0]], _gen_lEta[muI[0]], _gen_lPhi[muI[0]], _gen_lE[muI[0]]);
			l2.SetPtEtaPhiE(_gen_lPt[muI[1]], _gen_lEta[muI[1]], _gen_lPhi[muI[1]], _gen_lE[muI[1]]);
			double mll = (l1 + l2).M();	
			if(mll < 4) continue;
			hists[0][sam]->Fill(std::min(mll, 199.99), scal);
		}
		
	}
	plotHistRatio(hists[0][0], hists[0][1], "new", "old", "mllCompWZ" + extra,  false, 0,0, true, true);
	std::cout << "ratio WZ_minmll01/WZ_old = " << hists[0][0]->GetSumOfWeights()/hists[0][1]->GetSumOfWeights() << std::endl;

	//, SF = 0.7
	//plotHistRatio(hists[0][0], hists[1][0], "mt old", "mt new", "MT_newalgo_neutrinoZ" + extra, false, 0,0, true, true);
	
	//plotHistRatio(hists[1][0], hists[1][1], "mllmin01, SF = 0.7", "old powheg", "Mmumu_WZsampcomp" + extra,  false, 0,0, false, true);
	//plotHistRatio(hists[1][0], hists[1][1], "mllmin01, SF = 0.7", "old powheg", "Mmumu_WZsampcomp_log" + extra,  false, 0,0, true, true);
	
	//plotHist(hists[2][0], "lowGenMll");
	/*
	std::vector<TH1D*> distVec;
	std::vector<TString> distnames = {"old algorithm", "new algorithm", "unambiguous"};
	for(unsigned i = 0; i < 3; ++i){
		distVec.push_back(hists[i][0]);
	}
	plotHist(distVec, distnames, "WZ_newalgoComp", true);
	*/
}



int main(int argc, char* argv[]){
	TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;
	testtree.Loop();
	rootapp->Run();
    return 0;
}


