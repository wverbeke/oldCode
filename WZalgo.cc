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


#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"


void trilTree::Loop(){
	//Set plotting style
	setTDRStyle();
	//gROOT->SetBatch(kTRUE);
	//Define list of samples
	const unsigned nSamples = 1;
	//0.696758
	const TString fileList[nSamples] = {"WZTo3LNu_mllmin01.root"};
	const double xSections[nSamples] = {58.59*0.696758};
	const TString names[nSamples] = {"WZTo3LNu_mllmin01"};

	//Read Trees from ROOT files
	TFile* hfile[nSamples];
	TTree* inputTree[nSamples];
	double hcounter[nSamples];

	for(unsigned sam = 0; sam < nSamples; ++sam){   //CHANGE BACK TO NSAMPLES
		cout << "name " << names[sam] << endl;
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
	
	/*
	TMVA::Tools::Instance();

	TString outfileName("TMVA.root");
   	TFile* outputFile = TFile::Open(outfileName, "RECREATE" );
	TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
	                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
	dataloader->AddVariable( "mof", 'F' );
   	dataloader->AddVariable( "msf", 'F' );
   	//dataloader->AddVariable( "dphiOS", "units", 'F' );
   	//dataloader->AddVariable( "dphiSS",                "Variable 4", "units", 'F' );
	*/
	TFile treeFile("treetestWZ_ewkinoCuts.root","RECREATE");

	double mos, mother, dPhiMet, mt, wpt, deltaRos, lpt ,deltaPtos, motherAndMet, deltaEos, deltaPtMet, deltaBosonPt;
	TTree *tree = new TTree("signalTree","wzalgoTree correct");
    tree->Branch("mos",&mos,"mos/D");
    tree->Branch("mother", &mother, "mother/D");
	tree->Branch("dPhiMet", &dPhiMet, "dPhiMet/D");
	tree->Branch("mt", &mt, "mt/D");
	tree->Branch("wpt", &wpt, "wpt/D");
	tree->Branch("deltaBosonPt", &deltaBosonPt, "deltaBosonPt/D");
	tree->Branch("deltaRos", &deltaRos, "deltaRos/D");
	tree->Branch("lpt", &lpt, "lpt/D");
	tree->Branch("deltaPtos", &deltaPtos, "deltaPtos/D");
	tree->Branch("deltaEos", &deltaEos, "deltaEos/D");
	tree->Branch("motherAndMet", &motherAndMet, "motherAndMet/D");
	tree->Branch("deltaPtMet", &deltaPtMet, "deltaPtMet/D");

	TTree *tree2 = new TTree("bkgTree","wzalgoTree mispaired");
    tree2->Branch("mos",&mos,"mos/D");
    tree2->Branch("mother", &mother, "mother/D");
	tree2->Branch("dPhiMet", &dPhiMet, "dPhiMet/D");
	tree2->Branch("mt", &mt, "mt/D");
	tree2->Branch("wpt", &wpt, "wpt/D");
	tree2->Branch("deltaBosonPt", &deltaBosonPt, "deltaBosonPt/D");
	tree2->Branch("deltaRos", &deltaRos, "deltaRos/D");
	tree2->Branch("lpt", &lpt, "lpt/D");
	tree2->Branch("deltaPtos", &deltaPtos, "deltaPtos/D");
	tree2->Branch("deltaEos", &deltaEos, "deltaEos/D");
	tree2->Branch("motherAndMet", &motherAndMet, "motherAndMet/D");
	tree2->Branch("deltaPtMet", &deltaPtMet, "deltaPtMet/D");
	/*
	TMVA::Reader *readerWZ;
    Float_t mos, mother, dPhiMet;


    readerWZ = new TMVA::Reader( "!Color:!Silent" );
	readerWZ->AddVariable("mos", &mos);
	readerWZ->AddVariable("mother", &mother);
	readerWZ->AddVariable("dPhiMet", &dPhiMet);
    readerWZ->BookMVA( "BDT method", "/home/willem/Work/AnalysisCode/dataset/weights/TMVAClassification_BDT.weights.xml"); 
	*/
    Double_t scale[nSamples];
	//Loop over all samples
	for(unsigned sam = 0; sam < nSamples; ++sam){
		if(sam !=0) continue;
		Long64_t nEntries = inputTree[sam]->GetEntries();
		scale[sam] = xSections[sam]*DataLuminosity*1000/(hcounter[sam]);
		std::cout<<"Entries in "<< fileList[sam] <<" "<<nEntries<<std::endl;

		float progress = 0.0; //Progress bar 
	    for(Long64_t it = 0; it < nEntries; ++it){
			unsigned barWidth = 100;
			std::cout << "[";
			unsigned pos = barWidth*progress;
			for (unsigned i = 0; i < barWidth; ++i) {
				if (i < pos) std::cout << "=";
				else if (i == pos) std::cout << ">";
				else std::cout << " ";
			}
			std::cout << "] " << int(progress * 100.0) << " %\r";
			std::cout.flush();
			double a = (double) it;
			double b = (double) nEntries;
			progress = a/b; 

		    	inputTree[sam]->GetEntry(it);
	    	if(TestRun && it > 10000) break;
	    	double scal;
	    	scal = scale[sam]*_weight;

			//Apply HNL SELECTION
			//cutBased();
			//Baseline event selection:
			if(!baseline(true, false,true, false)) continue;	
			//Categorize according to the number of leptons and flavors
			unsigned* ind = new unsigned[_nL];
			unsigned lCount = lepOrder(ind, 3, true, true);	
			if(lCount != 3) continue; //Veto 4th FO lepton
			//Require 3 leptons to be tight in data and MC, and determine nonPrompt bkg in data
			unsigned nTight = tightCount(ind, lCount);
			bool tightFail = nTight < 3;
				
			//determine search category
			unsigned cat = hnl::cat(ind, _flavors, _charges, lCount, _lPt[ind[0]]);
			if(cat == 999) continue;
			if(cat != 0 && cat != 2 && cat != 4) continue; //This means there has to be an OSSF pair in the event
			//determine which leptons will be used for the calculation of mll
			TLorentzVector* lepV = new TLorentzVector[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				lepV[l].SetPtEtaPhiE(_lPt[ind[l]], _lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]);
			}
			unsigned mllI[2] = {99, 99};
			mllIndices(mllI, ind, lepV, _charges, _flavors, lCount);				
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
			if(_flavors[mllI[0]] != _flavors[mllI[1]]) continue;
			unsigned lw = 99;
			for(unsigned l = 0; l < lCount; ++l){
				if(ind[l] == mllI[0] || ind[l] == mllI[1]) continue;
				lw = ind[l];
			}
		
			//event selection!
			
			if(_met < 50) continue;
			if(!ptCuts_hnl(ind,lCount)) continue;
			if(!vetoLowMll(12)) continue;
			

			//only use unambiguous events 
			if(_flavors[mllI[0]] == _flavors[lw]) continue;
		    std::vector<Double_t> vars(2); // vector has size of number of input variables
			//find the most common charge
			unsigned posCount = 0;
			unsigned negCount = 0;
			for(unsigned l = 0; l < lCount; ++l){
				if(_charges[ind[l]] > 0) ++posCount;
				else ++negCount;
			}
			unsigned ssi;
			for(unsigned l = 0; l < lCount; ++l){
				if(posCount > negCount && _charges[ind[l]] < 0) ssi = l; 
				else if(negCount > posCount && _charges[ind[l]] > 0) ssi = l;
			}
			unsigned trueWi;
			unsigned fakeWi;
			for(unsigned l = 0; l < lCount; ++l){
				if(l == ssi) continue;
				if(_flavors[ind[l]] == _flavors[ind[ssi]]) fakeWi = l;
				else if(_flavors[ind[l]] != _flavors[ind[ssi]]) trueWi = l;
			}
			mos = (lepV[trueWi] + lepV[ssi]).M();
			mother = (lepV[fakeWi] + lepV[ssi]).M();
			TLorentzVector metVec;
			metVec.SetPtEtaPhiE(_met, 0, _met_phi, _met);
			mt = transmass(metVec, lepV[trueWi]);
			wpt = (metVec + lepV[trueWi]).Pt();
			deltaBosonPt = fabs( (lepV[fakeWi] + lepV[trueWi]).Pt() - wpt); 
			dPhiMet = fabs(lepV[trueWi].DeltaPhi(metVec));
			deltaRos = lepV[trueWi].DeltaR(lepV[ssi]);
			lpt = lepV[trueWi].Pt();
			deltaPtos = fabs(lepV[trueWi].Pt() - lepV[ssi].Pt());
			deltaEos = fabs(lepV[trueWi].E() - lepV[ssi].E());
			motherAndMet = transmass(lepV[fakeWi] + lepV[ssi], metVec);
			deltaPtMet = fabs(lepV[trueWi].E() - _met);
			//double mvaGood = readerWZ->EvaluateMVA( "BDT method");
			tree->Fill();

			mos = (lepV[fakeWi] + lepV[ssi]).M();
			mother = (lepV[trueWi] + lepV[ssi]).M();	
			mt = transmass(metVec, lepV[fakeWi]);
			wpt = (metVec + lepV[fakeWi]).Pt();
			deltaBosonPt = fabs( (lepV[trueWi] + lepV[trueWi]).Pt() - wpt); 
			dPhiMet = fabs(lepV[fakeWi].DeltaPhi(metVec));
			deltaRos = lepV[fakeWi].DeltaR(lepV[ssi]);
			lpt = lepV[fakeWi].Pt();
			deltaPtos = fabs(lepV[fakeWi].Pt() - lepV[ssi].Pt());
			deltaEos = fabs(lepV[fakeWi].E() - lepV[ssi].E());
			motherAndMet = transmass(lepV[trueWi] + lepV[ssi], metVec);
			deltaPtMet = fabs(lepV[fakeWi].E() - _met);
			//double mvaBad = readerWZ->EvaluateMVA( "BDT method");
			//std::cout << "mvaGood = " << mvaGood << "		 	mvaBad = " << mvaBad << std::endl;
			
			tree2->Fill();

			
			//std::cout << "signal : mof = " << mof << std::endl;
			//std::cout << "signal : msf = " << msf << std::endl;
			//dataloader->AddSignalTrainingEvent( vars, 1. );
			//mispairing! turn around the masses
			/*
			vars[1] = mof;
			vars[0] = msf;
		    dataloader->AddBackgroundTrainingEvent( vars, 1. );
			*/
		    			
		}
	}
	
	tree->Print();
	tree2->Print();
	treeFile.Write();
	treeFile.Close();
	/*
	TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   	TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";
    //dataloader->PrepareTrainingAndTestTree( mycuts, mycutb, "nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V" );
	//factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG","!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );
	//factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT");//, "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");
	factory->BookMethod( dataloader, TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );
	factory->TrainAllMethods();
	outputFile->Close();
	delete factory;
	delete dataloader;	
	*/
}



int main(int argc, char* argv[]){
	trilTree testtree;
	testtree.Loop();
    return 0;
}


