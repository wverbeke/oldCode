
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"


void Loop(){
	const TString fileName = "WZTo3LNu_mllmin01.root";
	//Read Trees from ROOT files
	/*
	TFile* file =  new TFile("../data_april17/" + fileName);
	file->cd("FakeElectrons");
	TTree* tree = (TTree*) (file->Get("FakeElectrons/fakeTree"));
	Double_t _lPt[20];
	Double_t _lEta[20];
	TBranch *b__lPt;
	TBranch *b__lEta;
	tree->SetBranchAddress("_lPt", _lPt, &b__lPt);
    tree->SetBranchAddress("_lEta", _lEta, &b__lEta);
	*/	

	TMVA::Tools::Instance();

	TString outfileName("TMVA.root");
   	TFile* outputFile = TFile::Open(outfileName, "RECREATE" );
	TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
	                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

    TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
	//dataloader->AddVariable( "a", 'F' , 0, 100);
   	//dataloader->AddVariable( "b", 'F' , 0, 2.5 );

	TFile* file =  new TFile("./treetestWZ.root");
	TTree* signalTree = (TTree*) (file->Get("signalTree"));
	TTree* bkgTree = (TTree*) (file->Get("bkgTree"));

	dataloader->AddSignalTree( signalTree,     1. );
  	dataloader->AddBackgroundTree(bkgTree, 1. );
	dataloader->AddVariable( "mos", 'F' );
   	dataloader->AddVariable( "mother", 'F' );
   	//dataloader->AddVariable( "dPhiMet", 'F' );
	dataloader->AddVariable( "mt", 'F' );
	dataloader->AddVariable( "deltaRos", 'F' );
	dataloader->AddVariable( "wpt", 'F' );
	//dataloader->AddVariable("deltaBosonPt", &deltaBosonPt);
	dataloader->AddVariable( "lpt", 'F' );
	dataloader->AddVariable( "deltaPtos", 'F' );
	dataloader->AddVariable( "deltaEos", 'F');
	dataloader->AddVariable( "motherAndMet", 'F');

	/*
	//TTree *signalTree     = (TTree*)input->Get("TreeS");
   	//TTree *background     = (TTree*)input->Get("TreeB");
	//Loop over all samples
	for(Long64_t it = 0; it < tree->GetEntries()/1000; ++it){
		tree->GetEntry(it);
		std::vector<Double_t> vars = {TMath::Min(_lPt[0], 99.9), TMath::Min(TMath::Abs(_lEta[0]), 2.49) };
		std::cout << "vars[0] = " << vars[0] << std::endl;
		std::cout << "vars[1] = " << vars[1] << std::endl;
		dataloader->AddSignalTrainingEvent( vars, 1. );
		vars[1] = TMath::Min(_lPt[0], 99.9);
		vars[0] =  TMath::Min(TMath::Abs(_lEta[0]), 2.49);
	    dataloader->AddBackgroundTrainingEvent( vars, 1. );
	}
	*/
	factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");
	//factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );
	//factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG", "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );
	//factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTB", "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );
	//factory->BookMethod( dataloader, TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );
	factory->TrainAllMethods();
	outputFile->Close();
	delete factory;
	delete dataloader;
}

