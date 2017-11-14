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
	const unsigned nSamples = 1;
	const unsigned nSamples_eff = 1;
	const TString fileList[nSamples] = {"WZsync.root"};
	const TString names[nSamples] = {"test"};

	const bool isData = (fileList[0] == "MuonEGsync.root");

	//Read Trees from ROOT files
	TFile* hfile[nSamples];
	TTree* inputTree[nSamples];
	double hcounter[nSamples];
	for(unsigned sam = 0; sam < nSamples; ++sam){
		cout << "name " << names[sam] << endl;
		hfile[sam] = new TFile("../data_trilEWK/"+fileList[sam],"read");
       	hfile[sam]->cd("FakeElectrons");
		//Determine hcounter for cross section scaling
		TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
		_hCounter->Read("hCounter");
		hcounter[sam] = _hCounter->GetBinContent(1);
       	inputTree[sam] = (TTree*) (hfile[sam]->Get("FakeElectrons/fakeTree"));
		Init(inputTree[sam], false, !isData);
	}
	//Read analysis scale factors
	this->readSF();
	

	
	//Tweakable options////////////////////////////////////////////////////
	const bool TestRun = false;	//Break after a few events
	const double DataLuminosity = 12.9;    //units of fb^{-1}
	const double CutBased = false;
	const TString extra = "_MuonEG";	//for plot file names
	//////////////////////////
	const unsigned nMiss = 2;
	const long missingEvents[nMiss] = {7043811, 17557284};

	const unsigned nMany = 0;
	const long manyEvents[1] = {938010};


	const TString eff_names[nSamples_eff] = {"test"};

	long unsigned count = 0;
	ofstream dump;
    dump.open ("sync_EWK" + extra +".txt");
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
			//Baseline event selection
			for(unsigned e = 0; e < nMiss; ++e){
				if(_eventNb == missingEvents[e]){
					cout << " crash 1 : "  << _eventNb << endl;
				}
			}
			if(_nL < 3) continue;
			if(!_METfilters_pass) continue;
			for(unsigned e = 0; e < nMiss; ++e){
				if(_eventNb == missingEvents[e]){
					cout << " crash 2 : "  << _eventNb << endl;
				}
			}
			//Clean taus
			for(unsigned l = 0; l < _nL; ++l){
				if(_flavors[l] == 2){
					_isloose[l] = _isloose[l] && _lClean[l];
					_isFO[l] = _isFO[l] && _lClean[l];
					_istight[l] = _istight[l] && _istight[l];
				}
			}
			//Calculate number of bjets 
			unsigned nJets = 0;
			unsigned nbJets = 0;
			double HT = 0;
			unsigned* jetInd = new unsigned[_nJets];
			for(unsigned j = 0; j < _nJets; ++j){
				TLorentzVector jet;
				jet.SetPtEtaPhiE(_jetPt[j], _jetEta[j], _jetPhi[j], _jetE[j]);
				bool overlap = false;
				for(unsigned l = 0; l < _nL; ++l){
					if(_isFO[l]){
						TLorentzVector lep;
						lep.SetPtEtaPhiE(_lPt[l], _lEta[l], _lPhi[l], _lE[l]);
						if(lep.DeltaR(jet) < 0.4){
							overlap = true;
							break;
						}
					}
				}
				if(!overlap){
					jetInd[nJets] = j;
					++nJets;
					if(_csv[j] > 0.8484){
						++nbJets;
						break;
					}
					if(_jetPt[j] > 30){
						HT+= _jetPt[j];
					}
				}
			}
			if(nbJets > 0) continue;
			for(unsigned e = 0; e < nMiss; ++e){
				if(_eventNb == missingEvents[e]){
					cout << " crash 3 : "  << _eventNb << endl;
				}
			}	
			//veto events with low mass mll since they aren't covered by the simulation
			
			bool lowM_pair = false;
			for(unsigned l1 = 0; l1 < _nL -1; ++l1){
				if(_isloose[l1]){
					if(_flavors[l1] == 0){
						if(!_clusterpass[l1]) continue;
					}
					for(unsigned l2 = l1 + 1; l2 < _nL; ++l2){
						if(_isloose[l2]){
							if(_flavors[l1] == 0){
								if(!_clusterpass[l2]) continue;
							}
							if(_flavors[l1] == _flavors[l2]){
								if(_charges[l1] != _charges[l2]){
									TLorentzVector lep1, lep2;
									lep1.SetPtEtaPhiE(_lPt[l1], _lEta[l1], _lPhi[l1], _lE[l1]);
									lep2.SetPtEtaPhiE(_lPt[l2], _lEta[l2], _lPhi[l2], _lE[l2]);
									if(_eventNb ==  17557284){
									if(_flavors[l1] == 0) cout << "ele";
									else if(_flavors[l1] == 1) cout <<"mu";
									else cout << "tau ";
									if(_charges[l1] > 0) cout << "+ AND ";
									else cout << "- AND ";
									if(_flavors[l2] == 0) cout << "ele";
									else if(_flavors[l2] == 1) cout <<"mu";
									else cout << "tau ";
									if(_charges[l2] > 0) cout << "+";
									else cout << "-";
									cout << "  make mass of " << (lep1 + lep2).M()  << endl;
									}
									if( (lep1 + lep2).M() < 12){
										lowM_pair = true;
										break;
									}
								}
							}
						}
					}
				}
				if(lowM_pair) break;
			}
			if(lowM_pair) continue;
			for(unsigned e = 0; e < nMiss; ++e){
				if(_eventNb == missingEvents[e]){
					cout << " crash 4 : "  << _eventNb << endl;
				}
			}
			//Categorize according to the number of leptons and flavors
			unsigned* ind = new unsigned[_nL];
			unsigned lCount = lepOrder(ind, 3);
			if(lCount < 3) continue;
			
			for(unsigned e = 0; e < nMiss; ++e){
				if(_eventNb == missingEvents[e]){
					cout << " crash 5 : "  << _eventNb << endl;
				}
			}

			//MuonEG cut, code this in cleaner way
			if(_flavors[ind[0]] == 1){
				bool muoneg = false;
				if(lCount == 3){
					if(_flavors[ind[1]] == 0 && _flavors[ind[2]] != 2) muoneg = true;
				} else if(lCount > 3){
					if(_flavors[ind[2]] == 0 && _flavors[ind[3]] != 1) muoneg = true;
				}
				if(muoneg){
					if(PtCone(_lPt[ind[0]], _flavors[ind[0]], _lepMVA[ind[0]], _ptratio[ind[0]]) < 25){
						for(unsigned l = 1; l < lCount; ++l){
							ind[l-1] = ind[l];
						}
						lCount -= 1;
					}
				}
			}
			//////////////////////////////////////

			for(unsigned l = 0; l < lCount; ++l){
				for(unsigned e = 0; e < nMiss; ++e){
						if(_eventNb == missingEvents[e]){
							if(_flavors[ind[l]] == 0) cout << "ele";
							else if(_flavors[ind[l]] == 1) cout <<"mu";
							else cout << "tau ";
							if(_charges[ind[l]] > 0) cout << "+" << endl;
							else cout << "-" << endl;
							cout << "cone pt of lepton "  << l  <<  " = " <<  PtCone(_lPt[ind[l]], _flavors[ind[l]], _lepMVA[ind[l]], _ptratio[ind[l]]) << endl;
						}
					}
			}

			unsigned nTight = 0;
			for(unsigned l = 0; l < lCount; ++l){
				if(_istight[ind[l]]){
					++nTight;
				} else{
					break;
				}
			}
			if(nTight < 3) continue;

			for(unsigned e = 0; e < nMany; ++e){
				if(_eventNb == manyEvents[e]){
					cout << " event should have failed : "  << _eventNb << endl;
					cout << "~~~~~~~~~~~~~~~~~~~~~~" << endl;
					for(unsigned l = 0; l < lCount; ++l){
						if(_flavors[ind[l]] == 0) cout << "ele";
						else if(_flavors[ind[l]] == 1) cout <<"mu";
						else cout << "tau ";
						if(_charges[ind[l]] > 0) cout << "+" << endl;
						else cout << "-" << endl;
						cout << PtCone(_lPt[ind[l]], _flavors[ind[l]], _lepMVA[ind[l]], _ptratio[ind[l]]) << endl;
					}
					cout << "~~~~~~~~~~~~~~~~~~~~~" << endl;
				}
			}

			for(unsigned e = 0; e < nMiss; ++e){
				if(_eventNb == missingEvents[e]){
					cout << " crash 6 : "  << _eventNb << endl;
					cout << "~~~~~~~~~~~~~~~~~~~~~~" << endl;
					for(unsigned l = 0; l < lCount; ++l){
						if(_flavors[ind[l]] == 0) cout << "ele";
						else if(_flavors[ind[l]] == 1) cout <<"mu";
						else cout << "tau ";
						if(_charges[ind[l]] > 0) cout << "+" << endl;
						else cout << "-" << endl;
						cout << PtCone(_lPt[ind[l]], _flavors[ind[l]], _lepMVA[ind[l]], _ptratio[ind[l]]) << endl;
					}
					cout << "~~~~~~~~~~~~~~~~~~~~~" << endl;
				}
			}
			lCount = nTight;
			
			//Trigger data events
			if(isData && !trigger(lCount)) continue;
			//determine search category
			unsigned cat = SR_EWK_cat(ind, _flavors, _charges, lCount);

			for(unsigned e = 0; e < nMany; ++e){
				if(_eventNb == manyEvents[e]){
					cout << " event should have failed : "  << _eventNb << endl;
					cout << "category : " << cat << endl;
				}
			}
			//if(cat < 6 && _met < 50) continue;
			//if(cat != 6) continue;
			//if( _met < 50 || (cat != 0 && cat != 4) ) continue;
			if( _met < 50 || (cat != 0) ) continue;
			//if( _met < 50 || (cat != 4) ) continue;

			for(unsigned e = 0; e < nMiss; ++e){
				if(_eventNb == missingEvents[e]){
					cout << " crash 7 : "  << _eventNb << endl;
				}
			}

			if(_eventNb == 14197111){
				for(unsigned l = 0; l < lCount; ++l){
					if(_flavors[ind[l]] == 0) cout << "ele";
					else if(_flavors[ind[l]] == 1) cout <<"mu";
					else cout << "tau ";
					if(_charges[ind[l]] > 0) cout << "+ AND ";
					else cout << "- AND ";
					if(_flavors[ind[l]] == 0) cout << "ele";
					else if(_flavors[ind[l]] == 1) cout <<"mu";
					else cout << "tau ";
					if(_charges[ind[l]] > 0) cout << "+";
					else cout << "-";
					cout << PtCone(_lPt[ind[l]], _flavors[ind[l]], _lepMVA[ind[l]], _ptratio[ind[l]]) << endl;
				}
			}
			//Apply analysis Pt and eta cuts
			unsigned nTau = 0;
			for(unsigned l = 0; l < lCount; ++l){
				if(_flavors[ind[l]] == 2) ++nTau;
			}
		 	if(nTau < 2){
				if(_flavors[ind[1]] == 1 && _lPt[ind[1]] < 10) continue;
				else if(_flavors[ind[1]] == 0 && _lPt[ind[1]] < 10) continue;
				unsigned* lightInd = new unsigned[lCount - nTau];
				for(unsigned l = 0, count = 0; l < lCount; ++l){
					if(_flavors[ind[l]] != 2){
						lightInd[count] = ind[l];
						++count;
					}
				}
				if(_flavors[lightInd[1]] == 0 && _lPt[lightInd[1]] < 15) continue;
				//if(_lPt[ind[2]] < 10) continue;
				bool otherlepnotmu = true;
				for(unsigned l = 1; l < lCount; ++l){
					if(_flavors[ind[l]] == 1) otherlepnotmu = false;
				}
				if(_flavors[lightInd[0]] == 1 && otherlepnotmu){
					if(_lPt[lightInd[0]] < 20) continue; //CHANGE BACK TO 25
				} else if(_flavors[lightInd[0]] == 1){
					if(_lPt[lightInd[0]] < 20) continue;
				} else{
					if(_lPt[lightInd[0]] < 25) continue;
				}
				bool trailFail = false;
				for(unsigned l = 2; l < lCount; ++l){
					if(_lPt[ind[l]] < 10){
						trailFail = true;
						break;
					}
				}
				if(trailFail) continue;					
			} else{
				bool ptFail = false;
				for(unsigned l = 0; l < lCount; ++l){
					if(_flavors[ind[l]] == 2 && _lPt[ind[l]] < 20){
						ptFail = true;
					} else if (_flavors[ind[l]] == 1 && _lPt[ind[l]] < 25){
						ptFail = true;
					} else if (_flavors[ind[l]] == 0 && _lPt[ind[l]] < 30){
						ptFail = true;
					}
					if( fabs(_lEta[ind[l]]) > 2.1) ptFail = true;
					if(ptFail) break;
				}
				if(ptFail) continue;
			}
			for(unsigned e = 0; e < nMiss; ++e){
				if(_eventNb == missingEvents[e]){
					cout << " crash 8 : "  << _eventNb << endl;
				}
			}

			//if(!promptMC(ind, lCount)) continue;
			//Require 3 leptons to be tight in data and MC, and determine nonPrompt bkg in data
			/*
			bool nlepFail = false;
			bool tightFail = false;
			unsigned tightCount = 0;
			for(unsigned l = 0; l < lCount; ++l){
				if(l < 3){
					if(!_istight[ind[l]]) tightFail= true;
				} else{
					if(_istight[ind[l]]) nlepFail = true;
				}
				if(nlepFail) break;
			}
			if(nlepFail) continue;		
			if(tightFail) continue;
			lCount = 3;
			*/
			
			if(cat == 999){
				//		cout << "Skipping 3 tau event" << endl;
				continue;	
			}
			/*
			//Apply ID and reco SF to simulation
				for(unsigned l = 0; l < lCount; ++l){  //CHANGE BACK BACK BACK
					if(_istight[ind[l]]){
						if(_flavors[ind[l]] == 2){
							scal*=0.83; //Twiki suggests flat 90% id SF for taus
						} else if(_flavors[ind[l]] == 0){
							scal*=idTightSFMap[0]->GetBinContent(idTightSFMap[0]->FindBin(TMath::Min(_lPt[ind[l]], 199.), fabs(_lEta[ind[l]])));
							scal*=idLooseSFMap[0]->GetBinContent(idLooseSFMap[0]->FindBin(TMath::Min(_lPt[ind[l]], 199.), fabs(_lEta[ind[l]])));
							scal*=recSFMap_ele->GetBinContent(recSFMap_ele->FindBin(_lEta[ind[l]]));
						} else if(_flavors[ind[l]] == 1){
							scal*=idTightSFMap[1]->GetBinContent(idTightSFMap[1]->FindBin(TMath::Min(_lPt[ind[l]], 119.), fabs(_lEta[ind[l]])));
							scal*=idLooseSFMap[1]->GetBinContent(idLooseSFMap[1]->FindBin(TMath::Min(_lPt[ind[l]], 119.), fabs(_lEta[ind[l]])));
							scal*=recSFMap_mu->Eval(_lEta[ind[l]]);
						}
					} else if(_isFO[ind[l]]){
						;
					} else if(_isloose[ind[l]]){
						;
					}
				}
			*/
			//Apply trigger SF to simulation
			/*
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
						scal *=  triggerEffEleLeg->GetBinContent( triggerEffEleLeg->FindBin( TMath::Min(_lPt[lightInd], 499.), TMath::Min(fabs(_lEta[lightInd]), 2.1)) );//Apply flat 86% efficiency for mtt events
					} else if(_flavors[lightInd] == 1){
						scal *= 0.86; 
					}
				}
			//Apply btag SF
				for(unsigned j = 0; j < nJets; ++j){
					//scal*=reader.eval_auto_bounds("central", BTagEntry::FLAV_B, _jetEta[jetInd[j]], _jetPt[jetInd[j]], _csv[jetInd[j]]);
				}
			//Apply PU reweighing
				//scal*= PUweights->GetBinContent(PUweights->FindBin(_n_trueInteractions));
			*/
			//determine which leptons will be used for the calculation of mll
			unsigned mllI[2] = {99, 99};
			TLorentzVector* lepV = new TLorentzVector[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				lepV[l].SetPtEtaPhiE(_lPt[ind[l]], _lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]);
			}
			
			mllIndices(mllI, ind, lepV, _charges, _flavors, lCount);
			//In 3 lepton events with ossf, remove large conversion background by cutting on m3l
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
			if(cat < 6)dump << Form("%1d %9d %12d\t%+2d %5.1f %d\t%+2d %5.1f %d\t%+2d %5.1f %d\t%d\t%2d\t%5.1f\t%6.1f\t%6.1f\t%6.1f \t %1.5f", _runNb, _lumiBlock, _eventNb, lepID[0], _lPt[ind[0]], 0, lepID[1],  _lPt[ind[1]], 0, lepID[2],  _lPt[ind[2]], 0, nJets, nbJets, _met, HT, mll, mt, scal) << std::endl;
			else dump << Form("%1d %9d %12d\t %+2d %5.1f %d\t %+2d %5.1f %d\t %+2d %5.1f %d\t%+2d %5.1f %d \t%d\t%2d\t%5.1f\t%6.1f \t %1.5f", _runNb, _lumiBlock, _eventNb, lepID[0], _lPt[ind[0]], 0, lepID[1],  _lPt[ind[1]], 0, lepID[2],  _lPt[ind[2]], 0, lepID[3], _lPt[ind[3]], 0,  nJets, nbJets, _met, HT, scal) << std::endl;
			//cout << lepID[3] << std::endl;
			++count;
    	}
	}
	dump.close();
	cout <<"number of entries : " << count << endl;
	
}



int main(int argc, char* argv[]){
	TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;
	testtree.Loop();
	rootapp->Run();
    return 0;
}


