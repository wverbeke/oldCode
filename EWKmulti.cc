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
	gROOT->SetBatch(kTRUE);

	//Initialize all samples and cross sections
	const unsigned nSamples = 19;
	const unsigned nSamples_eff = 6;
	const TString fileList[nSamples] = {"MuonEG.root", "DoubleMuon.root", "DoubleEG.root", "SingleElectron.root", "SingleMuon.root", "TTWJetsToLNu.root", "TTZToLLNuNu.root", "TTGJets.root", "ZZTo4L.root", "TTTT.root", "VHToNonbb.root", "WWW_4F.root", "WWZ.root", "WW_DoubleScattering.root", "WZZ.root", "WpWpJJ.root", "ZZZ.root", "ZGTo2LG.root", "WZTo3LNu.root"};
	const double xSections[nSamples - 5] = {0.2043, 0.2529, 3.697, 1.256, 0.009103, 0.9561, 0.2086, 0.1651, 1.61704, 0.05565, 0.03711, 0.01398, 123.9, 4.4297}; 
	const TString names[nSamples] = {"data", "data", "data", "data", "data", "TT + X", "TT + X", "TT + X", "ZZ", "rare SM", "rare SM", "rare SM", "rare SM", "rare SM", "rare SM", "rare SM", "rare SM", "conversion", "WZ"};

	//Read Trees from ROOT files
	TFile* hfile[nSamples];
	TTree* inputTree[nSamples];
	double hcounter[nSamples];
	for(unsigned sam = 0; sam < nSamples; ++sam){
		cout << "name " << names[sam] << endl;
		hfile[sam] = new TFile("../data_EWKmoriond/"+fileList[sam],"read");
       	hfile[sam]->cd("FakeElectrons");
		//Determine hcounter for cross section scaling
		TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
		_hCounter->Read("hCounter");
		hcounter[sam] = _hCounter->GetBinContent(1);
       	inputTree[sam] = (TTree*) (hfile[sam]->Get("FakeElectrons/fakeTree"));
		Init(inputTree[sam]);
	}

	//Read analysis scale factors:
	readSF();
	//Tweakable options////////////////////////////////////////////////////
	const bool TestRun = false;	//Break after a few events
	const double DataLuminosity = 36.22;    //units of fb^{-1}
	const TString extra = "";	//for plot file names
	const unsigned nCat = 12;
	//////////////////////////

	
	const TString eff_names[nSamples_eff + 1] = {"data", "TT + X", "ZZ", "rare SM", "conversion", "WZ", "non-prompt"};
	
	TH1D* yields[nCat][nSamples_eff + 1]; //Seperate histogram for every category
	const TString catNames[nCat] = {"3lOSSF", "3lnoOSSF", "OSSFtau", "OSOFtau", "SStau", "l2tau", "4l2OSSF", "4l1OSSF", "3ltau", "2l2tau2OSSF", "2l2tau1OSSF", "2lSS"};
	const unsigned nSR[nCat] = {44, 6, 18, 16, 12, 12, 5, 4, 4, 4, 3, 30};
	for(unsigned cat = 0; cat < nCat; ++cat){
		for(unsigned effsam = 0; effsam < nSamples_eff + 1; ++effsam){
			yields[cat][effsam] = new TH1D(catNames[cat] + eff_names[effsam], catNames[cat] + eff_names[effsam] + "; search region ; events/search region", nSR[cat], 0.5, nSR[cat] + 0.5);
		}
	}
	
	Double_t scale[nSamples -5];
	//set to check which data events have already been processed
	std::set<long> usedEvents[3]; //runNb, lumiBlock, eventNb
	
	//Loop over all samples
	for(unsigned sam = 0, effsam = 0; sam < nSamples; ++sam, ++effsam){
		if(sam != 0){
   			if(names[sam] == names[sam -1]) --effsam;
    	}
		if(fileList[sam] == "WWW_4F.root" || fileList[sam] == "WZZ.root" || fileList[sam] == "ZZZ.root") continue;

    	Long64_t nEntries = inputTree[sam]->GetEntries();
    	if(sam > 4){
   			scale[sam -5] = xSections[sam -5]*DataLuminosity*1000/(hcounter[sam]);
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
				scal = scale[sam-5]*_weight;
			}
			//Baseline event selection:
			if(_met < 50) continue;
			if(!baseline()) continue;
			//Check if data events were already used before

			if(effsam == 0){
				auto event = usedEvents[0].find(_eventNb);
				auto lumi = usedEvents[1].find(_lumiBlock);
				auto run = usedEvents[2].find(_runNb);
            	if (event != usedEvents[0].end() && lumi != usedEvents[1].end() && run != usedEvents[2].end() ) continue;
            	if(event == usedEvents[0].end()) usedEvents[0].insert(_eventNb); 
				if(lumi == usedEvents[1].end()) usedEvents[1].insert(_lumiBlock);
				if(run == usedEvents[2].end()) usedEvents[2].insert(_runNb);
			}
			//fill lepton indices and check there are enough FO leptons		 		
			unsigned* ind = new unsigned[_nL];
			unsigned lCount = lepOrder(ind, 3);
			if(lCount < 3) continue;
			//Apply analysis Pt and eta cuts
			if(!ptCuts(ind,lCount)) continue;
			//Require 3 leptons to be tight in data and MC, and determine nonPrompt bkg in data
			unsigned nTight = tightCount(ind, lCount);
			bool tightFail = nTight < 3;
			//Apply FR maps to data control region
			unsigned fill = effsam;
			if(tightFail && effsam == 0){
				//fakes go in different histogram
				fill = nSamples_eff;
				//Apply FR maps
				double* conePt = new double[lCount];
				for(unsigned l = 0; l < lCount; ++l){
					conePt[l] = PtCone(_lPt[ind[l]], _flavors[ind[l]], _lepMVA[ind[l]], _ptratio[ind[l]]);
				}
				scal*= fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap, 3);
			} else if(tightFail) continue;
			
			lCount = 3;
			//Apply triggers to data events;
			if(effsam == 0){
				bool trigPass[9];
				trigPass[0] = _tril_trigger_eee || _tril_trigger_all;
				trigPass[1] = _tril_trigger_eem || _tril_trigger_all;
				trigPass[2] = _tril_trigger_emm || _tril_trigger_all;
		 		trigPass[3] = _tril_trigger_mmm || _tril_trigger_all;
				trigPass[4] = _tril_trigger_eee || _tril_trigger_all; //eet and eee events have same trigger
				trigPass[5] = _tril_trigger_emt || _tril_trigger_all;
				trigPass[6] = _tril_trigger_mmm || _tril_trigger_all; //mmt and mmm events have same trigger
				trigPass[7] = _tril_trigger_ett || _tril_trigger_all;
				trigPass[8] = _tril_trigger_mtt || _tril_trigger_all;
				trigPass[9] = false;   //we do not consider 3tau events at the time
				if(!trigPass[tril_flavorComb(_flavors, lCount)]) continue;
			}						
			//determine search category
			unsigned cat = SR_EWK_cat(ind, _flavors, _charges, lCount);
			if(cat == 999){
				continue;	
			}
			if(cat > 5) continue; //Only considering 3l events here
			//Apply ID and reco SF to simulation
			/*
			if(effsam != 0){
				for(unsigned l = 0; l < lCount; ++l){
					if(_flavors[ind[l]] == 2){
						scal*=0.9; //Twiki suggests flat 90% id SF for taus
					} else if(_flavors[ind[l]] == 0){
						scal*=idSFMap[0]->GetBinContent(idSFMap[0]->FindBin(TMath::Min(_lPt[ind[l]], 199.), fabs(_lEta[ind[l]])));
						scal*=recSFMap_ele->GetBinContent(recSFMap_ele->FindBin(_lEta[ind[l]], TMath::Min(_lPt[ind[l]], 199.)));
					} else if(_flavors[ind[l]] == 1){
						scal*=idSFMap[1]->GetBinContent(idSFMap[1]->FindBin(TMath::Min(_lPt[ind[l]], 119.), fabs(_lEta[ind[l]])));
						scal*=recSFMap_mu->Eval(_lEta[ind[l]]);
					}
				}
			}
			//Apply trigger SF to simulation
			if(effsam != 0){
				if(cat < 2){ //3 light leptons
					scal*= triggerEffMap3l[_flavors[ind[2]]]->GetBinContent(triggerEffMap3l[_flavors[ind[2]]]->FindBin(TMath::Min(_lPt[ind[1]], 99.), TMath::Min(_lPt[ind[2]], 99.) ));
				} else if(cat < 5){ //2 light leptons, 1 tau
					unsigned lightInd[2];
					for(unsigned l = 0, k = 0; l < lCount; ++l){
						if(_flavors[ind[l]] != 2){
							lightInd[k] = ind[l];
							++k;
						}
					}
					scal*= triggerEffMap2l[_flavors[lightInd[0]]]->GetBinContent(triggerEffMap2l[_flavors[lightInd[0]]]->FindBin(TMath::Min(_lPt[lightInd[0]], 99.), TMath::Min(_lPt[lightInd[1]], 99.) ));
				} else{ //1 light lepton, 2 taus
					unsigned lightInd;
					for(unsigned l = 0, k = 0; l < lCount; ++l){
						if(_flavors[ind[l]] != 2){
							lightInd = ind[l];
							break;
						}
					}
					if(_flavors[lightInd] == 0){
						scal *=  triggerEffEleLeg->GetBinContent( triggerEffEleLeg->FindBin( TMath::Min(_lPt[lightInd], 499.), TMath::Min(_lEta[lightInd], 2.1)) );//Apply flat 86% efficiency for mtt events
					} else if(_flavors[lightInd] == 1){
						scal *= 0.86; 
					}
				}
			}
			//Apply btag SF
			if(effsam != 0){
				for(unsigned j = 0; j < nJets; ++j){
					//scal*=reader.eval_auto_bounds("central", BTagEntry::FLAV_B, _jetEta[jetInd[j]], _jetPt[jetInd[j]], _csv[jetInd[j]]);
				}
			}	
			//Apply PU reweighing
			if(effsam != 0){
				scal*= PUweights->GetBinContent(PUweights->FindBin(_n_trueInteractions));
			}
			*/
			//determine which leptons will be used for the calculation of mll
			unsigned mllI[2] = {99, 99};
			TLorentzVector* lepV = new TLorentzVector[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				lepV[l].SetPtEtaPhiE(_lPt[ind[l]], _lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]);
			}
			mllIndices(mllI, ind, lepV, _charges, _flavors, lCount);
			//In 3 lepton evenmts with ossf and ossf + tau
			if(cat == 0 || cat == 2){
				TLorentzVector tot;
				for(unsigned l = 0; l < 3; ++l) tot += lepV[l];
				if( fabs(tot.M() -  91) < 15) continue;
			}
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
			//determine MT and MT2ll
			TLorentzVector metV;
			metV.SetPtEtaPhiE(_met, 0, _met_phi, _met);
			double mt;
			if(cat < 2){
				//find index which belongs to the lepton from the W decay
				unsigned lw = 99;
				for(unsigned l = 0; l < lCount; ++l){
					if(ind[l] != mllI[0] && ind[l] != mllI[1]){
						lw  = l;
					}
				}
				mt = transmass(lepV[lw], metV); 
			}
			else{
				mt = find_mt2(ind, _flavors, _charges, lepV, metV, lCount);
			}
			unsigned searchR = SR_EWK_3lep(mt, _met, mll, cat);	
			/*if(names[sam] == "WZ"){
				cout << "category: " << cat << endl;
				cout << "search region: " << searchR << endl;
				cout << "effnames : " << eff_names[effsam] << endl;
			}*/
			yields[cat][fill]->Fill(searchR + 1, scal);
    	}
	}
	//Split data and MC histograms for plotting and propagating uncertainties
	TH1D* dataYields[6];
	for(unsigned cat = 0; cat < 6; ++cat) dataYields[cat] = (TH1D*) yields[cat][0]->Clone();
	TH1D* bkgYields[6][nSamples_eff];
	for(unsigned effsam = 1; effsam < nSamples_eff  + 1; ++ effsam){
		for(unsigned cat = 0; cat < 6; ++cat){
			bkgYields[cat][effsam -1] = (TH1D*) yields[cat][effsam]->Clone();
		}
	}
	//Propagate systematic uncertainties to the background
	
	double sampUnc[nSamples_eff] = {0.15, 0.25, 0.5, 0.2, 0.25, 0.36}; //25% dummy systematics added to WZ!
	for(unsigned cat = 0; cat < 6; ++cat){
		for(unsigned bkg = 0; bkg < nSamples_eff; ++bkg){
			if(cat < 2){
				double systUnc[7] = {0.05, 0.062, 0.04, 0.03, 0.03, 0.03, sampUnc[bkg]};
				addSyst(bkgYields[cat][bkg], systUnc, 7);
			} else if(cat < 5){
				double systUnc[7] = {0.05, 0.062, 0.04, 0.03, 0.03, 0.1166, sampUnc[bkg]};
				addSyst(bkgYields[cat][bkg], systUnc, 7);
			} else{
				double systUnc[7] = {0.05, 0.062, 0.03, 0.03, 0.1166, 0.1166, sampUnc[bkg]};
				addSyst(bkgYields[cat][bkg], systUnc, 7);
			}
		}
	}
	long full = 0;
	for(int i = 0; i < 36; ++i){
		cout << dataYields[0]->GetBinContent(i + 1) << endl;
		full += dataYields[0]->GetBinContent(i + 1);
	}
	cout << "total yield in category: " << full << endl;
	//Plot the yields as a function of the search region
	for(unsigned cat = 0; cat < 6; ++cat){
		plotDataVSMC(dataYields[cat], bkgYields[cat], eff_names, nSamples_eff, catNames[cat] + extra, true, 1);
		//printEWKTables(yields[cat], nSamples_eff, cat, catNames[cat]);
	}
}



int main(int argc, char* argv[]){
	TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;
	testtree.Loop();
	rootapp->Run();
    return 0;
}


