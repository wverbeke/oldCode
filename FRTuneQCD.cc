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


void trilTree::Loop(){
	//Set plotting style
	setTDRStyle();
	gROOT->SetBatch(kTRUE);

	//Initialize all samples and cross sections
	//Muon enriched QCD 
	/*
	const unsigned nSamples = 12;
	const TString fileList[nSamples] = {"QCD_MuEnriched_Pt15to20.root", "QCD_MuEnriched_Pt20to30.root", "QCD_MuEnriched_Pt30to50.root", "QCD_MuEnriched_Pt50to80.root", "QCD_MuEnriched_Pt80to120.root", "QCD_MuEnriched_Pt120to170.root", "QCD_MuEnriched_Pt170to300.root", "QCD_MuEnriched_Pt300to470.root", "QCD_MuEnriched_Pt470to600.root", "QCD_MuEnriched_Pt600to800.root", "QCD_MuEnriched_Pt800to1000.root", "QCD_MuEnriched_Pt1000toInf.root"};
	const TString names[nSamples] = {"QCD_MuEnriched", "QCD_MuEnriched", "QCD_MuEnriched", "QCD_MuEnriched", "QCD_MuEnriched", "QCD_MuEnriched", "QCD_MuEnriched", "QCD_MuEnriched", "QCD_MuEnriched", "QCD_MuEnriched", "QCD_MuEnriched", "QCD_MuEnriched"};
	const double xSections[nSamples] = {1273190000, 558528000, 139803000, 19222500, 2758420, 469797, 117989, 7820.25, 645.528, 187.109, 32.3486, 10.4305};
	*/
	//EMEnriched and bcToE QCD
	const unsigned nSamples = 13;
	const TString fileList[nSamples] = {"QCD_EMEnriched_Pt20to30.root", "QCD_EMEnriched_Pt30to50.root", "QCD_EMEnriched_Pt50to80.root", "QCD_EMEnriched_Pt80to120.root", "QCD_EMEnriched_Pt120to170.root", "QCD_EMEnriched_Pt170to300.root", "QCD_EMEnriched_Pt300toInf.root", "QCD_bcToE_Pt15to20.root", "QCD_bcToE_Pt20to30.root", "QCD_bcToE_Pt30to80.root", "QCD_bcToE_Pt80to170.root", "QCD_bcToE_Pt170to250.root", "QCD_bcToE_Pt250toInf.root"};
	const TString names[nSamples] = {"QCD_EMEnriched", "QCD_EMEnriched", "QCD_EMEnriched", "QCD_EMEnriched", "QCD_EMEnriched", "QCD_EMEnriched", "QCD_EMEnriched", "QCD_EMEnriched", "QCD_EMEnriched", "QCD_EMEnriched", "QCD_EMEnriched", "QCD_EMEnriched", "QCD_EMEnriched"};
	const double xSections[nSamples] = {557600000, 136000000, 19800000, 2800000, 477000, 114000, 9000, 1272980000, 557627000, 159068000, 3221000, 105771, 21094.1};
	
	//Read Trees from ROOT files
	TFile* hfile[nSamples];
	TTree* inputTree[nSamples];
	double hcounter[nSamples];
	for(unsigned sam = 0; sam < nSamples; ++sam){
		cout << "name " << fileList[sam] << endl;
		hfile[sam] = new TFile("../data_april17/"+fileList[sam],"read");
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
	const TString extra = "";
	//////////////////////////

	//const unsigned nTest = 18;  //number of WP to test
	const unsigned nTest = 160;
	TH1D* udsHist[2][3][nTest];	//[2]: FR numerator and denominator [3]: eta bins  [nTest]: entry for every FO WP to test
	TH1D* bcHist[2][3][nTest];
	//const unsigned nBins = 6;
	//const double ptBins[nBins + 1] = {5,10,20,30, 50, 70, 100};
	//const unsigned nBins = 5;
	//const double ptBins[nBins + 1] = {10,20,30, 50, 70, 100}; //100};
	//const unsigned nBins = 6;
	//const double ptBins[nBins + 1] = {5,10, 15, 25, 35, 50, 70};
	const unsigned nBins = 5;
	const double ptBins[nBins + 1] = {10, 15, 25, 35, 50, 70};

	/*
	const TString denNum[2] = {"numerator", "denominator"};
	for(unsigned i = 0; i < 2; ++i){
		for(unsigned et = 0; et < 3; ++et){
			for(unsigned test = 0; test < nTest; ++test){
				udsHist[i][et][test] = new TH1D("udsPt_etaBin" + std::to_string(et) + "_" + denNum[i] + "test_" + std::to_string(test) , "udsPt_etaBin" + std::to_string(et) + "_" + denNum[i] + "test_" + std::to_string(test) + ";" + "P_{T}(lepton) (GeV)" + "; events", nBins, ptBins);
				udsHist[i][et][test]->Sumw2();
				bcHist[i][et][test] = new TH1D("bcPt_etaBin" + std::to_string(et) + "_" + denNum[i] + "test_" + std::to_string(test) , "bcPt_etaBin" + std::to_string(et) + "_" + denNum[i] + "test_" + std::to_string(test) + ";" + "P_{T}(lepton) (GeV)" + "; events", nBins, ptBins);
				bcHist[i][et][test]->Sumw2();
			}
		}
	}
	*/

	TH1D* fr[2][3]; //num, denom, 3 delta R's
	const TString denNum[2] = {"numerator", "denominator"};
	const TString deltaRNames[3] = {"#DeltaR > 1", "#DeltaR > 1.5", "#DeltaR > 2"};
	for(unsigned i = 0; i < 2; ++i){
		for(unsigned r = 0; r < 3; ++r){
			fr[i][r] = new TH1D("fr" + deltaRNames[r] + denNum[i], "fr" + deltaRNames[r] + denNum[i] + ";" + "P_{T}(lepton) (GeV)" + "; Fake-rate", nBins, ptBins);
			fr[i][r]->Sumw2();
		}
	}


	double scale[nSamples];

	const double etaMVACut15[3] = {0.77, 0.56, 0.48};
	const double etaMVACut25[3] = {0.52, 0.11, -0.01};

	const double etaBins[4] = {0., 0.8, 1.479, 2.5};

	
	TH2D* FRMapEle[2];
	FRMapEle[0] = new TH2D("FRMapEle_num", "FRMapEle_num; P_{T}(GeV) ; |#eta(GeV)|", nBins, ptBins, 3, etaBins);
	FRMapEle[1] = new TH2D("FRMapEle_denom", "FRMapEle_denom; P_{T}(GeV) ; |#eta(GeV)|", nBins, ptBins, 3, etaBins);
	FRMapEle[0]->Sumw2();
	FRMapEle[1]->Sumw2();

	//const double mvaWPFO[3] = {0.48, 0.10, -0.36};
	const double mvaWPFO[3] = {-0.02, -0.52, -0.52};
	//const double mvaWPFO[3] = {0.5, 0.1, -0.02};
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

			cutBased();
			/*
			//Select exactly one very loose lepton
			unsigned ind = 0;
			unsigned lCount = 0;
			for(unsigned l = 0; l < _nL; ++l){
				if(_islooseCut[l] && _flavors[l] == 0 && _isolation[l] < 0.6 && _clusterpass[l] && fabs(_lEta[l]) < 2.5 && _3dIPsig[l] < 4){
				//if(_islooseCut[l] && _flavors[l] == 1 && _isolation[l] < 0.6 && _isFOCut[l] &&_3dIPsig[l] < 4){
					ind = l;
					++lCount;
				}
			}
			if(lCount != 1) continue;
			//Require one jet back to back to the loose lepton
			bool btbJet = false;
			TLorentzVector lep;
			lep.SetPtEtaPhiE(_lPt[ind], _lEta[ind], _lPhi[ind], _lE[ind]);
			for(unsigned j = 0; j < _nJets; ++j){
				TLorentzVector jet;
				jet.SetPtEtaPhiE(_jetPt[j], _jetEta[j], _jetPhi[j], _jetE[j]);
				if(jet.DeltaR(lep) > 1){
					btbJet = true;
					break;
				}
			}
			if(!btbJet) continue;
			//End baseline event selection
			//Determine eta bin;
			unsigned eta = (fabs(_lEta[ind]) > 0.8) + (fabs(_lEta[ind]) > 1.479);
			*/
			//Apply extra criteria to the lepton to tune FO and store FO denominator
			/*
			for(unsigned test = 0; test < nTest; ++test){
				double mvaToPass = std::min( (-1. + ((1. + etaMVACut15[eta])/nTest)*test), std::max( (-1. + ((1. + etaMVACut25[eta])/nTest)*test) , (-1. + ((1. + etaMVACut15[eta])/nTest)*test) + ( (-1. + ((1. + etaMVACut25[eta])/nTest)*test) - (-1 + ((1. + etaMVACut15[eta])/nTest)*test) )*0.1*(_lPt[ind]-15.) ) );
	
				//if(_isolation[ind] < (1 - test*0.05)){
				if(_mvaValue[ind] > mvaToPass){
					if(_origin[ind] == 3 || _origin[ind] == 4){
						//udsHist[1][0][test]->Fill(_lPt[ind], scal);	//CHANGE THE INDEX TO FILL TO ETA FOR ELECTRONS!	
						udsHist[1][eta][test]->Fill(std::min(_lPt[ind]*std::max(1., 1 + (_isolation[ind] - 0.1)), ptBins[nBins]), scal);
					} else if(_origin[ind] == 1 || _origin[ind] == 2){
						//bcHist[1][0][test]->Fill(_lPt[ind], scal);
						bcHist[1][eta][test]->Fill(std::min(_lPt[ind]*std::max(1., 1 + (_isolation[ind] - 0.1)), ptBins[nBins]), scal);
					}
				}
			}
			//Store tight numerator
			if(_istightCut[ind]){
				if(_origin[ind] == 3 || _origin[ind] == 4){
					//udsHist[0][0][0]->Fill(_lPt[ind], scal);	//CHANGE THE INDEX TO FILL TO ETA FOR ELECTRONS!
					udsHist[0][eta][0]->Fill(std::min(_lPt[ind]*std::max(1., 1 + (_isolation[ind] - 0.1)), ptBins[nBins]), scal);
				} else if(_origin[ind] == 1 || _origin[ind] == 2){
					//bcHist[0][0][0]->Fill(_lPt[ind], scal);	
					bcHist[0][eta][0]->Fill(std::min(_lPt[ind]*std::max(1., 1 + (_isolation[ind] - 0.1)), ptBins[nBins]), scal);
				}
			}
			*/
			//if(match[ind] >= 20 || fabs(_gen_lmompdg[match[ind]]) <= 3 || fabs(_gen_lmompdg[match[ind]]) <= 5){
			/*
				unsigned eta = (fabs(_lEta[ind]) > 0.8) + (fabs(_lEta[ind]) > 1.479);
				if(_istightCut[ind]){
					FRMapEle[0]->Fill(std::min(_lPt[ind]*std::max(1., 1 + (_isolation[ind] - 0.1)), ptBins[nBins]), fabs(_lEta[ind]), scal);
				}
				if(_mvaValue[ind] > mvaWPFO[eta]){
					FRMapEle[1]->Fill(std::min(_lPt[ind]*std::max(1., 1 + (_isolation[ind] - 0.1)), ptBins[nBins]), fabs(_lEta[ind]), scal);
				}
			*/
			//}
			/*
			if(_istightCut[ind]){
				FRMapEle[0]->Fill(std::min(_lPt[ind]*std::max(1., 1 + (_isolation[ind] - 0.1)), ptBins[nBins]), fabs(_lEta[ind]), scal);
			}
			FRMapEle[1]->Fill(std::min(_lPt[ind]*std::max(1., 1 + (_isolation[ind] - 0.1)), ptBins[nBins]), fabs(_lEta[ind]), scal);
			*/

			//Select exactly one very loose lepton
			unsigned ind = 0;
			unsigned lCount = 0;
			for(unsigned l = 0; l < _nL; ++l){
				//if(_islooseCut[l] && _flavors[l] == 0 && _isolation[l] < 0.6 && _clusterpass[l] && fabs(_lEta[l]) < 2.5 && _3dIPsig[l] < 4){
				if(_isFO[l] && _flavors[l] == 0){
				//if(_islooseCut[l] && _flavors[l] == 1 && _isolation[l] < 0.6 && _isFOCut[l] &&_3dIPsig[l] < 4){
					ind = l;
					++lCount;
				}
			}
			if(lCount != 1) continue;
			//Require one jet back to back to the loose lepton
			double maxDeltaR = 0;
			TLorentzVector lep;
			lep.SetPtEtaPhiE(_lPt[ind]*std::max(1., 1 + (_isolation[ind] - 0.1)), _lEta[ind], _lPhi[ind], _lE[ind]*std::max(1., 1 + (_isolation[ind] - 0.1)) );
			for(unsigned j = 0; j < _nJets; ++j){
				TLorentzVector jet;
				jet.SetPtEtaPhiE(_jetPt[j], _jetEta[j], _jetPhi[j], _jetE[j]);
				if(jet.DeltaR(lep) > maxDeltaR){
				 	maxDeltaR  = jet.DeltaR(lep);
				}
			}
			if(maxDeltaR < 1) continue;
			unsigned eta = (fabs(_lEta[ind]) > 0.8) + (fabs(_lEta[ind]) > 1.479);
			if(_istight[ind]){
				fr[0][0]->Fill(std::min(_lPt[ind]*std::max(1., 1 + (_isolation[ind] - 0.1)), ptBins[nBins]), scal);
			}
			if(_flavors[ind] == 1 || _mvaValue[ind] > mvaWPFO[eta]){
				fr[1][0]->Fill(std::min(_lPt[ind]*std::max(1., 1 + (_isolation[ind] - 0.1)), ptBins[nBins]), scal);
			}
			if(maxDeltaR > 1.5){
				if(_istight[ind]){
					fr[0][1]->Fill(std::min(_lPt[ind]*std::max(1., 1 + (_isolation[ind] - 0.1)), ptBins[nBins]), scal);
				}
				if(_flavors[ind] == 1 || _mvaValue[ind] > mvaWPFO[eta]){
					fr[1][1]->Fill(std::min(_lPt[ind]*std::max(1., 1 + (_isolation[ind] - 0.1)), ptBins[nBins]), scal);
				}
				if(maxDeltaR > 2.){
					if(_istight[ind]){
						fr[0][2]->Fill(std::min(_lPt[ind]*std::max(1., 1 + (_isolation[ind] - 0.1)), ptBins[nBins]), scal);
					}
					if(_flavors[ind] == 1 || _mvaValue[ind] > mvaWPFO[eta]){
						fr[1][2]->Fill(std::min(_lPt[ind]*std::max(1., 1 + (_isolation[ind] - 0.1)), ptBins[nBins]), scal);
					}
				}
			}
			
    	}
	}
	
	//TFile* frFile = new TFile("../QCDFRMap_electrons.root", "recreate");
	//Calculate all Fake-rates
	/*
	for(unsigned et = 0; et < 3; ++et){
		for(unsigned test = 0; test < nTest; ++test){
		*/
			/*
			if((et > 0 && (-1 + (1./nTest)*test) == -0.52) || (et == 0 && (-1 + (1./nTest)*test) == -0.02)){
				TH1D* udsnum = (TH1D*) udsHist[1][et][test]->Clone();
				TH1D* bcnum = (TH1D*) bcHist[1][et][test]->Clone();
				TH1D* udsden = (TH1D*)udsHist[0][et][0]->Clone();
				TH1D* bcden = (TH1D*) bcHist[0][et][0]->Clone();
				udsnum->Add(bcnum);
				udsden->Add(bcden);
				udsden->Divide(udsnum);
				TString name = "QCDFR_electron_eta" + (TString) std::to_string(et);
				udsnum->Write(name);
			}
			*/
			/*
			TH1D* numlight = (TH1D*) udsHist[0][et][0]->Clone();
			TH1D* numheavy = (TH1D*) bcHist[0][et][0]->Clone();

			numlight->Divide(udsHist[1][et][test]);
			numheavy->Divide(bcHist[1][et][test]);
			udsHist[1][et][test] = (TH1D*) numlight->Clone();
			bcHist[1][et][test] = (TH1D*) numheavy->Clone();
			*/
			//cout << "udsHist[1][et][test] entries = " << udsHist[1][et][test]->GetEntries() << endl;
			//cout << "bcHist[1][et][test] entries = " << bcHist[1][et][test]->GetEntries() << endl;
			//udsHist[1][et][test]->Divide(udsHist[0][et][0]);
			//bcHist[1][et][test]->Divide(bcHist[0][et][0]);
		/*
		}
	}
	//Determine the ratio of light to heavy fakerates
	for(unsigned et = 0; et < 3; ++et){
		double mindiff = 999;
		unsigned mindiffTest = 0;
		for(unsigned test = 0; test < nTest; ++test){
		*/
			/*
			cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
			//cout << "ISOLATION THRESHOLD: " << 1 - test*0.05 << endl;
			//cout << "MVA THRESHOLD: " << -1 + (1./nTest)*test << endl;
			//cout << "MVA THRESHOLD: " << -1 + (1./nTest)*test << endl; (1. + etaMVACut[eta])
			//TF1 *f = new TF1("f", "[0]", 5, 100);
			TF1 *f = new TF1("f", "[0]", 10, 70);
			//cout << "udsHist[1][et][test] entries = " << udsHist[1][et][test]->GetEntries() << endl;
			//cout << "bcHist[1][et][test] entries = " << bcHist[1][et][test]->GetEntries() << endl;
			//TF1 *f = new TF1("f", "[0]", 10, 100);
			*/
			/*
			udsHist[1][0][test]->Divide(bcHist[1][0][test]);
			udsHist[1][0][test]->Fit(f,"RSq");
			udsHist[1][0][test]->Fit(f,"RSq");
			*/
			/*
			udsHist[1][et][test]->Divide(bcHist[1][et][test]);
			udsHist[1][et][test]->Fit(f,"RSq");
			udsHist[1][et][test]->Fit(f,"RSq");
			if(fabs(f->GetParameter(0) - 1.) < mindiff){
				mindiff = fabs(f->GetParameter(0) - 1.);
				mindiffTest = test;
			}
			plotHist(udsHist[1][et	][test], "frFits/udsPt_electron _test" + std::to_string(test) + "_eta_" + std::to_string(et) );
		}
		cout << "######################################" << endl;
		cout << "######################################" << endl;
		//cout << "best closure got at ISO = " << 1 - mindiffTest*0.05 << endl;
		//cout << "best closure got at MVA = " << -1 + (1./nTest)*mindiffTest << endl;
		cout << "best closure got at MVA 15 = " << -1 + ((1. + etaMVACut15[et])/nTest)*mindiffTest << endl;
		cout << "best closure got at MVA 25	 = " << -1 + ((1. + etaMVACut25[et])/nTest)*mindiffTest << endl;
		cout << "mindiffTest = " << mindiffTest << endl;
		cout << "######################################" << endl;
	}

	TFile* frFile = new TFile("../QCDFRMap_muon.root", "recreate");
	FRMapEle[0]->Divide(FRMapEle[1]);
	FRMapEle[0]->Write("QCDFRMap_muon");
	frFile->Close();
	*/
	
	/*
	std::vector<TH1D*> distVec[nSamples];
	std::vector<TString> distnamesIso = {"iso leading", "iso subleading", "iso trailing"};
	std::vector<TString> distnamesPt = {"leading gen P_{T}", "subleading gen P_{T}", "trailing gen P_{T}" };
	const TString samNames[nSamples] = {"m5", "m15", "m20", "m30", "m40", "m50", "m60", "m80", "m100", "m130", "m150", "m200", "m300", "m400"};
	for(unsigned sam = 0; sam < nSamples; ++sam){
		for(unsigned i = 0; i < 3; ++i){
			distVec[sam].push_back(Histos[i][sam]);
		}
		plotHist(distVec[sam], distnamesIso, "signalIso" + samNames[sam] + extra);
	} 
	*/
	/*
	for(unsigned sam = 0; sam < nSamples; ++sam){
		distVec[sam].clear();
		for(unsigned i = 3; i < 6; ++i){
			distVec[sam].push_back(Histos[i][sam]);
		}
		plotHist(distVec[sam], distnamesPt, "genPt" + samNames[sam] + extra);
	}  	
	*/
		

	//Calculate fakeRate;
	for(unsigned r = 0; r < 3; ++r){
		fr[0][r]->Divide(fr[1][r]);
	}
	//plot fakeRate:

	std::vector<TH1D*> hists = {fr[0][0], fr[0][1], fr[0][2]};
	std::vector<TString> histNames = {"#DeltaR(lepton, jet) > 1", "#DeltaR(lepton, jet) > 1.5", "#DeltaR(lepton, jet) > 2"};
	plotHist(hists, histNames, "frDeltaRVariation_electron", true, true);
		

}



int main(int argc, char* argv[]){
	TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;
	testtree.Loop();
	rootapp->Run();
    return 0;
}


