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

	//Initialize all samples and cross sections
	const unsigned nSamples = 1;
	const unsigned nSamples_eff = 1;
	const TString fileList[nSamples] = {"test.root"};
	const TString names[nSamples] = {"test"};

	const bool isData = (fileList[0] == "MuonEGsync.root");

	//Read Trees from ROOT files
	TFile* hfile[nSamples];
	TTree* inputTree[nSamples];
	double hcounter[nSamples];
	for(unsigned sam = 0; sam < nSamples; ++sam){
		cout << "name " << names[sam] << endl;
		hfile[sam] = new TFile("../data_april17/"+fileList[sam],"read");
		//hfile[sam] = new TFile("../data_april17/"+fileList[sam],"read");	
       	hfile[sam]->cd("FakeElectrons");
		//Determine hcounter for cross section scaling
		TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
		_hCounter->Read("hCounter");
		hcounter[sam] = _hCounter->GetBinContent(1);
       	inputTree[sam] = (TTree*) (hfile[sam]->Get("FakeElectrons/fakeTree"));
		Init(inputTree[sam], false, !isData);
	}
	//Read analysis scale factors
	//this->readSF(true);
	

	
	//Tweakable options////////////////////////////////////////////////////
	const bool TestRun = false;	//Break after a few events
	const double DataLuminosity = 35.9;    //units of fb^{-1}
	const TString extra = "";	//for plot file names
	//////////////////////////

	//Events for printing
	
	const unsigned nMiss = 1;
	const long missingEvents[nMiss] = {32416474};

	const unsigned nMany = 1;
	const long manyEvents[nMany] = {1};

	const TString eff_names[nSamples_eff] = {"test"};

	ofstream dump;
    dump.open ("sync_hnl" + extra +".txt");

	//Loop over all samples
	for(unsigned sam = 0; sam < nSamples; ++sam){
		
    	Long64_t nEntries = inputTree[sam]->GetEntries();
    	
    	std::cout<<"Entries in "<< fileList[sam] <<" "<<nEntries<<std::endl;
		
        for(Long64_t it = 0; it < nEntries; ++it){
        	inputTree[sam]->GetEntry(it);
        	if (it%10000 == 0) cout<<'.'<<flush;
        	if(TestRun && it > 10000) break;
        	double scal;
        	scal = 1;
			
			bool found = false;
			for(unsigned m = 0; m < nMiss; ++m){
				if(_eventNb == missingEvents[m]){
					found = true;
					break;
				}
			}
			if(found){
				cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
				cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
				cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
				cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
				cout << _eventNb << "	crash 1" << endl;
				cout << "_nL = " << _nL << endl;
				cout << "nBJets(true, false, 0) = " << nBJets(true, false, 0) << endl;
			}
			
			//Apply HNL SELECTION
			cutBased();
			//Baseline event selection
			//if(!baseline(true, true, false, false)) continue;
			if(found) cout << _eventNb << "	crash 2" << endl;
			//if(nBJets(true, false, 0) != 0) continue;

			if(found) cout << _eventNb << "	crash 3" << endl;

			//Categorize according to the number of leptons and flavors
			unsigned* ind = new unsigned[_nL];
			//unsigned lCount = lepOrder(ind, 3, true, true);
			unsigned lCount = lepOrder(ind, 1, true, true);


			if(found){
				unsigned counter = 0;
				for(unsigned l = 0; l < _nL; ++l){
					//if(_isFO[l]) ++counter;
					//if(!_isFO[l]){
						cout << "##############################" << endl;
						if(_flavors[l] == 0) cout << "electron : " << endl;
						else if(_flavors[l] == 1) cout << "muon : " << endl;
						else continue;
						if(_islooseCut[l]) cout << "pass loose" << endl;
						else cout << "fail loose" << endl;
						cout << "eta =  " << _lEta[l] << endl;
						cout << "pt = " << _lPt[l] << endl;
						cout << "iso = " << _isolation[l] << endl;
						if(_flavors[l] == 0)  cout << "mva = " << _mvaValue[l] << endl;
						cout << "3dipsig = " << _3dIPsig[l] << endl;
						if(_flavors[l] == 0) cout << "hitsNumber = " << _hitsNumber[l] << endl;
						if(_flavors[l] == 0 && _clusterpass[l] ) cout << "pass cluster" << endl;
						if(_flavors[l] == 0 && !_clusterpass[l] ) cout << "fail cluster" << endl;
					}
				//}
				cout << "lCount = " << counter << endl;
				for(unsigned l = 0; l < lCount; ++l){
					if(_flavors[ind[l]] == 0) cout << " electron ";
					else if(_flavors[ind[l]] == 1) cout << " muon ";
					else cout << " tau ";
					cout << "	 cone Pt = " << _lPt[ind[l]]*(1 + std::max(_isolation[ind[l]] - 0.1, 0.)) << endl;
				}
			}
			if(found) cout << _eventNb << "	before selection 1 lepton" << endl;
			unsigned nl = 0;
			unsigned indl = 0;
			for(unsigned l = 0; l < _nL; ++l){
				if(fabs(_ipPV[l]) > 0.05) continue;
				if(fabs(_ipZPV[l]) > 0.1) continue;
				if(_isolation[l] > 0.6) continue;
				if(_flavors[l] == 0){
					if(fabs(_lEta[l]) > 2.5) continue;
					if(_lPt[l] < 10) continue;
					if(_hitsNumber[l] != 0) continue;
					if(_vtxFitConversion[l]) continue;
				} else if(_flavors[l] == 1){
					if(fabs(_lEta[l]) > 2.4) continue;
					if(_lPt[l] < 5) continue;
				} else continue;
				indl = l;
				++nl;
			}
			if(nl != 1) continue;
			if(found) cout << _eventNb << "	crash 4" << endl;
			if(lCount != 1) continue; //Veto 4th FO lepton considering signal model!
			//Apply analysis Pt thresholds
			//if(!ptCuts_hnl(ind,lCount)) continue;

			if(found) cout << _eventNb << "	crash 5" << endl;
			//MC prompt matching
			/*
			bool promptfail = false;
			for(unsigned l = 0; l < lCount; ++l){	
				if(_origin[ind[l]] != 0){
					promptfail = true;
					break;
				}
			}
			*/
			//if(promptfail) continue;


			if(found) cout << _eventNb << "	crash 6" << endl;

			if(_eventNb == 248233){
				cout << "FAKE SYNC PRINT : " << endl;
				for(unsigned l = 0; l < lCount; ++l){
					if(_isFO[ind[l]] && !_istight[ind[l]]){
						if(_flavors[ind[l]] == 0) cout << "electron : ";
						else cout << "muon :";
						cout << "	eta = " << _lEta[ind[l]] << "		" << "cone Pt = " << _lPt[ind[l]]*(1 + std::max(_isolation[ind[l]] - 0.1, 0.)) << endl;
						cout << "_isolation = " << _isolation[ind[l]] << endl;
						if(_clusterpass[ind[l]]) cout << "passed cluster" << endl;
						cout << "mva value = " << _mvaValue[ind[l]] << endl;
					}
				}
			}
			
			double* conePt = new double[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				conePt[l] = _lPt[ind[l]]*(1 + std::max(_isolation[ind[l]] - 0.1, 0.));
			}
			
			if(_lPt[indl]*(1 + std::max(_isolation[indl] - 0.1,0.)) < 10 - 5*_flavors[ind[0]]) continue;
			//Require 3 leptons to be tight in data and MC, and determine nonPrompt bkg in data
			/*
			unsigned nTight = tightCount(ind, lCount);
			bool tightFail = nTight < 3;
			if(tightFail) continue;
			//index used to fill events, needed to separate fakes from data
			//unsigned fill = effsam;
			//Apply FR maps to data control region
			double* conePt = new double[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				conePt[l] = _lPt[ind[l]]*(1 + std::max(_isolation[ind[l]] - 0.1, 0.));
			}
			if(conePt[0] < 15) continue;
			if(conePt[1] < 10) continue;
			if(conePt[2] < 10 - 5*_flavors[ind[2]]) continue;
			if(tightFail){
				//Apply FR maps
				//scal *= -1; //MC fake subtraction
				scal*= fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap, lCount);
			} //else if(tightFail) continue;
			//Apply triggers to data events;
			bool trigPass[4];
			trigPass[0] = _lowM_trigger_eee || _lowM_trigger_all;
			trigPass[1] = _lowM_trigger_mee || _lowM_trigger_all;
			trigPass[2] = _lowM_trigger_mme || _lowM_trigger_all;
			trigPass[3] = _lowM_trigger_mmm || _lowM_trigger_all;
			//if(!trigPass[tril_flavorComb(ind, _flavors, lCount)]) continue;	
			//determine search category
			TLorentzVector* lepV = new TLorentzVector[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				lepV[l].SetPtEtaPhiE(conePt[l], _lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]*(conePt[l]/_lPt[ind[l]]) );
			}
			//Calculate lepton system vector
			TLorentzVector lepSyst;
			for(int l = 0; l < 3; ++l) lepSyst += lepV[l];
			*/
			/*
			unsigned cat = hnl::cat(ind, _flavors, _charges, lCount, conePt[0]);
			if(cat == 999){
				continue;	
			}
			*/
			/*
			//Apply ID and reco SF to simulation
			if(effsam != 0){
				scal*=getEventSF(ind, lCount, true);
			}
			*/
			//determine which leptons will be used for the calculation of mll
			/*
			unsigned mllI[2] = {99, 99};
			mllIndices(mllI, ind, lepV, _charges, _flavors, lCount);				
			//determine mll
			double mll;
			if(mllI[0] == 99){
				if(mllI[1] != 99) std::cerr << "error one mll index is not -1 while the other is" << endl;
				mll = -1;
			} else{
				TLorentzVector lzV[2];
				for(unsigned l = 0; l < 2; ++l) lzV[l].SetPtEtaPhiE(PtCone(_lPt[mllI[l]]*(1 + std::max(_isolation[mllI[l]] - 0.1, 0.)), _flavors[mllI[l]], _lepMVA[mllI[l]], _ptratio[mllI[l]]), _lEta[mllI[l]], _lPhi[mllI[l]], _lE[mllI[l]]*(1 + std::max(_isolation[mllI[l]] - 0.1, 0.)) );
				mll = (lzV[0] + lzV[1]).M();
			}
			
			//determine the index of the W lepton
			unsigned lw = 9999;
			for(unsigned l = 0; l < lCount; ++l){
				if(ind[l] != mllI[0] && ind[l] != mllI[1]){
					lw = ind[l];
				}
			}
			//Category spefific selection:
			//low mass general
			if(cat < 4){
				//Calculate minDeltaR and maxDeltaR
				double minDeltaR = 9999.;
				double maxDeltaR = 0.;
				for(unsigned l = 0; l < lCount - 1; ++l){
					for(unsigned k = l + 1; k < lCount; ++k){
						if(lepV[l].DeltaR(lepV[k]) > maxDeltaR) maxDeltaR = lepV[l].DeltaR(lepV[k]);
						if(lepV[l].DeltaR(lepV[k]) < minDeltaR) minDeltaR = lepV[l].DeltaR(lepV[k]);
					}
				}
				if((cat == 0 || cat == 2) && minDeltaR < 0.05) continue;
				if((cat == 0 || cat == 2) && maxDeltaR < 2) continue;
				if(lepSyst.M() > 80) continue;
				//if(_met > 75) continue;
			}
			else if(cat > 3){
				if(conePt[1] < 15) continue;
				if(conePt[2] < 10 - 5*(_flavors[ind[2]])) continue;
				if(cat == 4){
					if(fabs(mll - 91) < 15) continue; //Veto onZ events (WZ and DY)
					if(fabs(lepSyst.M() - 91) < 15) continue; //veto conversions
				}
				if(!vetoLowMll(5)) continue;
			}
			*/
			/*
			//Calculate MET vector
			TLorentzVector METvec;
			METvec.SetPtEtaPhiE(_met, 0, _met_phi, _met);
			//Calculate min(Mos)
			double minMos = 0;
			unsigned minI[2] = {99, 99};
			for(unsigned l = 0; l < lCount -1 ; ++l){
				for(unsigned k = l + 1; k < lCount; ++k){
					if(_charges[ind[l]] != _charges[ind[k]]){
						if( (lepV[l] + lepV[k]).M() < minMos  || minMos == 0){
							minMos = (lepV[l] + lepV[k]).M();
							minI[0] = l;
							minI[1] = k;
						}
					}
				}
			}
			//find lepton for MT calculation
			unsigned lw_min;
			for(unsigned l = 0; l < lCount; ++l){
				if(l != minI[0] && l != minI[1]){
					lw_min = l;
				}
			}
			//fill Nominal yields
			double mt_min = transmass(lepV[lw_min], METvec);			
			*/
			//Calculate MET vector
			if(found) cout << _eventNb << " crash 7 " << endl;
			double deltaR = 0;
			TLorentzVector lepton;
			lepton.SetPtEtaPhiE(_lPt[indl]*(1 + std::max(_isolation[indl] - 0.1,0.))  , _lEta[indl], _lPhi[indl], _lE[indl]*(1 + std::max(_isolation[indl] - 0.1,0.)));
			for(unsigned j = 0; j < _nJets; ++j){
				if(_jetPt[j] < 30) continue;
				TLorentzVector jet;
				jet.SetPtEtaPhiE(_jetPt[j], _jetEta[j], _jetPhi[j], _jetE[j]);
				if(found) cout << "deltaR = " << jet.DeltaR(lepton) << endl;
				if(jet.DeltaR(lepton)> deltaR){
					deltaR = jet.DeltaR(lepton);
				}
			}
			if(deltaR < 1) continue;
			/*
			int* lepID = new int[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				if(_flavors[ind[l]] == 0){
					lepID[l] = 11;
				} else if(_flavors[ind[l]] == 1){
					lepID[l] = 13;
				} else if(_flavors[ind[l]] == 2){
					lepID[l] = 15;
				}
				if(_charges[ind[l]] > 0) lepID[l]*=-1;
			}
			*/
			if(found) cout << _eventNb << " crash end " << endl;
			//dump << Form("%1d %9d %12d\t%+2d %5.1f %d\t%+2d %5.1f %d\t%+2d %5.1f %d\t%d\t%2d\t%5.1f\t%6.1f\t%6.1f\t%6.1f \t %1.5f", _runNb, _lumiBlock, _eventNb, lepID[0], conePt[0], 0, lepID[1],  conePt[1], 0, lepID[2],  conePt[2], 0, 0, 0, _met, 0., minMos, mt_min, scal) << std::endl;
			dump << _eventNb << endl;
			//dump << Form("%1d %9d %12d\t%+2d %5.1f %d\t%+2d %5.1f %d\t%+2d %5.1f %d\t%d\t%2d\t%5.1f\t%6.1f\t%6.1f\t%6.1f \t %1.5f", _runNb, _lumiBlock, _eventNb, lepID[0], conePt[0], 0, 0,  0, 0, 0,  0, 0, 0, 0, _met, 0., 0., 0., scal) << std::endl;
    	}
	}
	dump.close();
	//cout <<"number of entries : " << count << endl;
	
}



int main(int argc, char* argv[]){
	TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;
	testtree.Loop();
	rootapp->Run();
    return 0;
}


