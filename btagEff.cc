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
#include <string>
using std::cout;
using std::endl;
using std::flush;

//include other parts of the code
#include "trilTree.h"

void trilTree::Loop(){
	//Set plotting style
	//setTDRStyle();
	gROOT->SetBatch(kTRUE);
	//Initialize all samples and cross sections
	const unsigned nSamples = 2;
	const TString fileList[nSamples] = {"TTJets_SingleLeptFromT_madgraph.root", "TTJets_SingleLeptFromTbar_madgraph.root"};
	const TString names[nSamples] = {"TT", "TT"};
	const double xSections[nSamples] = {182.175, 182.175};
	//Read Trees from ROOT files
	TFile* hfile[nSamples];
	TTree* inputTree[nSamples];
	double hcounter[nSamples];
	for(unsigned sam = 0; sam < nSamples; ++sam){
		cout << "name " << fileList[sam] << endl;
		hfile[sam] = new TFile("./"+fileList[sam],"read");
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
	const double DataLuminosity = 35.9;
	//const double bTagWP = 0.5426;
	const double bTagWP = 0.8484;
	const bool cleanJets = true;
	const bool hnlSel = false;
	const TString extra = "";
	//////////////////////////

	//const unsigned nTest = 18;  //number of WP to test
	const unsigned nTest = 80;
	TH2D* bTaggEff[2][3]; //2: num and denom 3: udsg, c, b
	const double ptBins[17] = {25,30,35,40,45,50,60,70,80,90,100,120,150,200,300,400,600};
	const double etaBins[7] = {0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
	const TString denNum[2] = {"numerator", "denominator"};
	const TString quarkFlav [3] = {"udsg", "charm", "beauty"};
	for(unsigned i = 0; i < 2; ++i){
		for(unsigned flav = 0; flav < 3; ++flav){
			bTaggEff[i][flav] = new TH2D("bTagEff" + quarkFlav[flav] + denNum[i], "bTagEff" + quarkFlav[flav] + denNum[i] + ";" + "P_{T}(jet) (GeV)" + "; |#eta(jet)|", 16, ptBins, 6, etaBins);
			bTaggEff[i][flav]->Sumw2();
		}
	}
	double scale[nSamples];
	//Loop over all samples
	for(unsigned sam = 0; sam < nSamples; ++sam){
    	Long64_t nEntries = inputTree[sam]->GetEntries();
    	std::cout<<"Entries in "<< fileList[sam] <<" "<<nEntries<<std::endl;
		scale[sam] = xSections[sam]*DataLuminosity*1000/(hcounter[sam]);
        for(Long64_t it = 0; it < nEntries; ++it){
        	inputTree[sam]->GetEntry(it);
        	if (it%10000 == 0) cout<<'.'<<flush;
        	if(TestRun && it > 10000) break;
			double scal= scale[sam]*_weight;
			for(unsigned j = 0; j < _nJets; ++j){
				if(fabs(_jetEta[j]) > 2.4) continue;
				if(_jetPt[j] < 25) continue;
				if(_jetFlavour[j] != 0 && _jetFlavour[j] != 4 && _jetFlavour[j] != 5) continue;
				if(cleanJets){
					TLorentzVector jet;
					jet.SetPtEtaPhiE(_jetPt[j], _jetEta[j], _jetPhi[j], _jetE[j]);
					bool clean = true;
					for(unsigned l = 0; l < _nL; ++l){
						if(_isFO[l]){
							TLorentzVector lep;
							lep.SetPtEtaPhiE(_lPt[l], _lEta[l], _lPhi[l], _lE[l]);
							if(jet.DeltaR(lep) < 0.4){
								clean = false;
								break;
							}
						}
					}	
					if(!clean) continue;
				}
				unsigned flav = 0 + (_jetFlavour[j] == 4) + 2*(_jetFlavour[j] == 5);
				if(_csv[j] > bTagWP){
					bTaggEff[0][flav]->Fill(std::min(_jetPt[j], 599.), std::min( fabs(_jetEta[j]), 2.4), scal);
				}
				bTaggEff[1][flav]->Fill(std::min(_jetPt[j], 599.), std::min( fabs(_jetEta[j]), 2.4), scal);	
			}			
    	}
	}
	
	TFile* frFile;
	if(hnlSel) frFile = new TFile("../weights/bTagEff.root", "recreate");
	else frFile = new TFile("../weights/bTagEff_ewkino.root", "recreate");
	for(unsigned flav = 0; flav < 3; ++flav){
		bTaggEff[0][flav]->Divide(bTaggEff[1][flav]);
		bTaggEff[0][flav]->Write("btagEff_" + quarkFlav[flav]);
	}
	frFile->Close();
}



int main(int argc, char* argv[]){
	TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;
	testtree.Loop();
	rootapp->Run();
    return 0;
}


