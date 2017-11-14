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
	//gROOT->SetBatch(kTRUE);
	//Define list of samples
	const unsigned nSamples = 1;
	//0.696758
	const TString fileList[nSamples] = {"WZTo3LNu_mllmin01.root"};
	const double xSections[nSamples] = {58.59*0.655};
	const TString names[nSamples] = {"WZ"};

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
	
	//Make histograms containing kinematic distributions
	const unsigned nDist = 3;
	const TString distNames[nDist] = {"mt", "mll", "genMll"};
	const TString xAxes[nDist] = {"M_{T}", "M_{#mu#mu}", "M_{#mu#mu}"};
	const TString units[nDist] = {"GeV", "GeV", "GeV"};
	const double histMin[nDist] = {0, 0, 0};
	const double histMax[nDist] = {300, 200, 2};
	const int nBins[nDist] = {100, 100, 100};
	TH1D* hists[nDist][nSamples];
	for(unsigned dist = 0; dist < nDist; ++dist){
		float binWidth = (histMax[dist] - histMin[dist])/nBins[dist];
		std::ostringstream strs; strs << binWidth; std::string yAxis = strs.str();
		for(unsigned sam = 0; sam < nSamples; ++sam){
			hists[dist][sam] = new TH1D(distNames[dist] + names[sam], distNames[dist] + names[sam] + ";" + xAxes[dist] + " (" + units[dist] +  "); events /" + yAxis + units[dist], nBins[dist], histMin[dist], histMax[dist]);
		}
	}
	
	double maxBinC[nDist];
	for(unsigned dist = 0; dist < nDist; ++dist){
		maxBinC[dist] = hists[dist][0]->GetBinCenter(hists[dist][0]->GetNbinsX());
	}
	

	const TString sourceName[3] = {"correct pairing", "mispaired", "#tau decay"};
	const TString mismName[3] = {"#Delta MET < 30 GeV", "30 GeV < #Delta MET < 60 GeV", "#Delta MET > 60 GeV"};
	TH1D* mtSplit[3][3];
	for(unsigned s = 0; s < 3; ++s){
		for(unsigned m = 0; m < 3; ++m){
			mtSplit[s][m] = new TH1D(sourceName[s] + mismName[m], sourceName[s] + mismName[m] + "; M_{T}(GeV); events /5 GeV", 60, 0, 300);
		}
	}

	const TString angmismName[3] = {"#Delta#Phi(MET) < 0.3", "0.3 < #Delta #Phi(MET) < 0.6", "#Delta #Phi(MET) > 0.6"};
	TH1D* mtAngSplit[3][3];
	for(unsigned s = 0; s < 3; ++s){
		for(unsigned m = 0; m < 3; ++m){
			mtAngSplit[s][m] = new TH1D(sourceName[s] + angmismName[m], sourceName[s] + angmismName[m] + "; M_{T}(GeV); events /5 GeV", 60, 0, 300);
		}
	}


	TH1D* mtResSplit[3][3];
	const TString resName[3] = {"MET + 0 GeV", "MET + 20 GeV", "MET - 20 GeV"};
	for(unsigned s = 0; s < 3; ++s){
		for(unsigned m = 0; m < 3; ++m){
			mtResSplit[s][m] = new TH1D(sourceName[s] + resName[m], sourceName[s] + resName[m] + "; M_{T}(GeV); events /5 GeV", 60, 0, 300);
		}
	}


	TH1D* mtResAngSplit[3][3];
	const TString angresName[3] = {"#PhiMET", "#PhiMET + #Delta#PhiMET", "#PhiMET - #Delta#PhiMET"};
	for(unsigned s = 0; s < 3; ++s){
		for(unsigned m = 0; m < 3; ++m){
			mtResAngSplit[s][m] = new TH1D(sourceName[s] + angresName[m], sourceName[s] + angresName[m] + "; M_{T}(GeV); events /5 GeV", 60, 0, 300);
		}
	}

	TH1D* leptonRes[2];
	leptonRes[0] = new TH1D("leptonres_mism", "leptonres_mism ; P_{T} - gen P_{T}/gen P_{T}; Events", 60, -2, 2);
	leptonRes[1] = new TH1D("leptonres_corr", "leptonres_corr ; P_{T} - gen P_{T}/gen P_{T}; Events", 60, -2, 2);

	TH1D* DeltaPhiGen[2];
	DeltaPhiGen[0] = new TH1D("DeltaPhiGen_corr", "DeltaPhiGen_corr ; #Delta#Phi(M_{T} lepton, MET); Events", 30, 0, 3.2);
	DeltaPhiGen[1] = new TH1D("DeltaPhiGen_mism", "DeltaPhiGen_mism ; #Delta#Phi(M_{T} lepton, MET); Events", 30, 0, 3.2);

    Double_t scale[nSamples];

	double mtRatio_corr[2][3];
	double mtRatio_misp[2][3];
	double mtRatioAng_corr[2][3];
	double mtRatioAng_misp[2][3];
	//Loop over all samples
	for(unsigned sam = 0; sam < nSamples; ++sam){
		Long64_t nEntries = inputTree[sam]->GetEntries();
		scale[sam] = xSections[sam]*DataLuminosity*1000/(hcounter[sam]);
		std::cout<<"Entries in "<< fileList[sam] <<" "<<nEntries<<std::endl;
	    for(Long64_t it = 0; it < nEntries; ++it){
			if (it%10000 == 0) cout<<'.'<<flush;
	    	inputTree[sam]->GetEntry(it);
	    	if(TestRun && it > 10000) break;
	    	double scal;
	    	scal = scale[sam]*_weight;
			//Apply HNL SELECTION
			cutBased();
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
			//if((cat == 0 || cat == 2 || cat == 4) && fabs(mll - 91) < 15) continue; //consider onZ events for WZ CR

			
			//determine the index of the W lepton
			unsigned lw = 9999;
			for(unsigned l = 0; l < lCount; ++l){
				if(ind[l] != mllI[0] && ind[l] != mllI[1]){
					lw = ind[l];
				}
			}
			if(_flavors[lw] !=  _flavors[mllI[0]]) continue;
			//if(_flavors[lw] == 1) continue;
			
			

			
			
			unsigned* match = new unsigned[_nL];
			matchGenLep(match);
			
			unsigned source;
			unsigned mism;
			unsigned angmism;
			if(_gen_nL == 0 || _nL == 0) continue;

			//cout << "match[lw] = " << match[lw] << endl;

			//Calculate min(Mos) and Mos(min Delta R)
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
			//Calculate alternete mt values
			unsigned lw_min;
			for(unsigned l = 0; l < lCount; ++l){
				if(l != minI[0] && l != minI[1]){
					lw_min = l;
				}
			}
			lw = lw_min;
			mllI[0] = minI[0];
			mllI[1] = minI[1];

			if(match[lw] > 20 || match[mllI[0]] > 20 ||  match[mllI[1]] > 20) continue;
			
			TLorentzVector Wlep, metV;
			Wlep.SetPtEtaPhiE(_lPt[lw], _lEta[lw], _lPhi[lw], _lE[lw]);
			metV.SetPtEtaPhiE(_met, 0, _met_phi, _met);
			//double fill[nDist] = {transmass(Wlep,metV), mll};
			/*
			//if(transmass(Wlep, metV) < 100) continue;
			for(unsigned dist = 0; dist < nDist; ++dist){
				hists[dist][sam]->Fill(std::min(fill[dist], maxBinC[dist]), scal);
			}	
			*/
			//cout << "lw = " << lw << endl;
			//cout << "match[lw] = " << match[lw] << endl;
			//cout << _gen_lmompdg[match[lw]] << endl;
			if( fabs(_gen_lmompdg[match[lw]]) == 24){
				source = 0;
			} else if(fabs(_gen_lmompdg[match[lw]]) == 23 ){
				source = 1;
			}  else if( fabs(_gen_lmompdg[match[lw]]) == 15){
				source = 2;
			} else if(_gen_isPromptl[match[lw]] && ( fabs(_gen_lmompdg[match[mllI[0]]]) == 23 && fabs(_gen_lmompdg[match[mllI[1]]]) == 23) ){
				source = 0;
			} else{
				source = 1;
			}	
			
			//cout << _genmet << endl;
			if(fabs(_met - _genmet) < 30){
				mism = 0;
			} else if(fabs(_met - _genmet) < 60){
				mism = 1;
			} else{
				mism = 2;
			}

			if(fabs(_met_phi - _genmet_phi) < 0.3){
				angmism = 0;
			} else if(fabs(_met_phi - _genmet_phi) < 0.6){
				angmism = 1;
			} else{
				angmism = 2;
			}
			mtSplit[source][mism]->Fill(TMath::Min(transmass(metV, Wlep), 299.), scal);	
			mtAngSplit[source][angmism]->Fill(TMath::Min(transmass(metV, Wlep), 299.), scal);	
			/*
			mtResAngSplit[source][0]->Fill(TMath::Min(transmass(metV, Wlep), 299.), scal);	
			
			mtResSplit[source][0]->Fill(TMath::Min(transmass(metV, Wlep), 299.), scal);	

			if(transmass(metV, Wlep) > 160.){
				if(source == 0) mtRatio_corr[0][0] = mtRatio_corr[0][0] +  transmass(metV, Wlep);
				if(source == 1) mtRatio_misp[0][0] =  mtRatio_misp[0][0] +  transmass(metV, Wlep);
			}
			else if(transmass(metV, Wlep) > 100.){
				if(source == 0) mtRatio_corr[1][0] = mtRatio_corr[1][0] +  transmass(metV, Wlep);
				if(source == 1) mtRatio_misp[1][0] =  mtRatio_misp[1][0] +  transmass(metV, Wlep);
			}
			double oldmet = _met;
			_met = oldmet + fabs(_met - _genmet);
			//cout << fabs(_met - _genmet) << endl;
			metV.SetPtEtaPhiE(_met, 0, _met_phi, _met);
			mtResSplit[source][1]->Fill(TMath::Min(transmass(metV, Wlep), 299.), scal);	
			if(transmass(metV, Wlep) > 160.){
				if(source == 0) mtRatio_corr[0][1] = mtRatio_corr[0][1] +  transmass(metV, Wlep);
				if(source == 1) mtRatio_misp[0][1] =  mtRatio_misp[0][1] +  transmass(metV, Wlep);
			}
			else if(transmass(metV, Wlep) > 100.){
				if(source == 0) mtRatio_corr[1][1] = mtRatio_corr[1][1] +  transmass(metV, Wlep);
				if(source == 1) mtRatio_misp[1][1] =  mtRatio_misp[1][1] +  transmass(metV, Wlep);
			}
			_met = oldmet - fabs(_met - _genmet);
			metV.SetPtEtaPhiE(_met, 0, _met_phi, _met);	
			mtResSplit[source][2]->Fill(TMath::Min(transmass(metV, Wlep), 299.), scal);
			if(transmass(metV, Wlep) > 160.){
				if(source == 0) mtRatio_corr[0][2] = mtRatio_corr[0][2] +  transmass(metV, Wlep);
				if(source == 1) mtRatio_misp[0][2] =  mtRatio_misp[0][2] +  transmass(metV, Wlep);
			}
			else if(transmass(metV, Wlep) > 100.){
				if(source == 0) mtRatio_corr[1][2] = mtRatio_corr[1][2] +  transmass(metV, Wlep);
				if(source == 1) mtRatio_misp[1][2] =  mtRatio_misp[1][2] +  transmass(metV, Wlep);
			}

			_met = oldmet;
			double oldmet_phi = _met_phi;
			_met_phi = oldmet_phi + fabs(_genmet_phi - _met_phi);
			metV.SetPtEtaPhiE(_met, 0, _met_phi, _met);	
			mtResAngSplit[source][1]->Fill(TMath::Min(transmass(metV, Wlep), 299.), scal);	
			if(transmass(metV, Wlep) > 160.){
				if(source == 0) mtRatioAng_corr[0][1] = mtRatioAng_corr[0][1] +  transmass(metV, Wlep);
				if(source == 1) mtRatioAng_misp[0][1] =  mtRatioAng_misp[0][1] +  transmass(metV, Wlep);
			}
			else if(transmass(metV, Wlep) > 100.){
				if(source == 0) mtRatioAng_corr[1][1] = mtRatioAng_corr[1][1] +  transmass(metV, Wlep);
				if(source == 1) mtRatioAng_misp[1][1] =  mtRatioAng_misp[1][1] +  transmass(metV, Wlep);
			}
			_met_phi = oldmet_phi - fabs(_met_phi - _genmet_phi);
			metV.SetPtEtaPhiE(_met, 0, _met_phi, _met);	
			mtResAngSplit[source][2]->Fill(TMath::Min(transmass(metV, Wlep), 299.), scal);
			if(transmass(metV, Wlep) > 160.){
				if(source == 0) mtRatioAng_corr[0][2] = mtRatioAng_corr[0][2] +  transmass(metV, Wlep);
				if(source == 1) mtRatioAng_misp[0][2] =  mtRatioAng_misp[0][2] +  transmass(metV, Wlep);
			}
			else if(transmass(metV, Wlep) > 100.){
				if(source == 0) mtRatioAng_corr[1][2] = mtRatioAng_corr[1][2] +  transmass(metV, Wlep);
				if(source == 1) mtRatioAng_misp[1][2] =  mtRatioAng_misp[1][2] +  transmass(metV, Wlep);
			}
			if(source == 0) leptonRes[1]->Fill( (_gen_lPt[match[lw]]- _lPt[lw])/_gen_lPt[match[lw]], scal);
			if(source == 1) leptonRes[0]->Fill(	(_gen_lPt[match[lw]]- _lPt[lw])/_gen_lPt[match[lw]], scal);

			if(source != 2){
				TLorentzVector genMetV, genWlep;
				genMetV.SetPtEtaPhiE(_genmet, 0, _genmet_phi, _genmet);
				genWlep.SetPtEtaPhiE(_gen_lPt[match[lw]], _gen_lEta[match[lw]], _gen_lPhi[match[lw]], _gen_lE[match[lw]]);
				double deltaphi = genMetV.DeltaPhi(genWlep);
				DeltaPhiGen[source]->Fill(deltaphi, scal);
			}
			*/
		}
	}
	
	TH1D* mtTot = (TH1D*) mtSplit[0][0]->Clone();
	for(unsigned source = 0; source < 3; ++source){
		for(unsigned i = 0; i < 3; ++i){
			if(i == 0 && source == 0) continue;
			mtTot->Add(mtSplit[source][i]);
		}

	}
	TH1D* mtSourceSplit[3];
	for(unsigned source = 0; source < 3; ++source){
		mtSourceSplit[source] = (TH1D*) mtSplit[source][0]->Clone();
		for(unsigned i = 1; i < 3; ++i){
			mtSourceSplit[source]->Add(mtSplit[source][i]);
		}
	}
	TH1D* mtMismSplit[3];
	for(unsigned m = 0; m < 3; ++m){
		mtMismSplit[m] = (TH1D*) mtSplit[0][m]->Clone();
		for(unsigned i = 1; i < 3; ++i){
			mtMismSplit[m]->Add(mtSplit[i][m]);
		}
	}
	
	TH1D* totalmisp = (TH1D*)  mtSourceSplit[1]->Clone();
	TH1D* totalCorr = (TH1D*)  mtSourceSplit[0]->Clone();
	const TString plotNames[4] = {"total WZ", "correct pairing", "mispaired", "#tau decay"};	
	TString plotNames2[4] = {"total WZ", "#Delta MET < 30 GeV", "30 GeV < #Delta MET < 60 GeV", "#Delta MET > 60 GeV"};
	
	plotDataVSMC(mtTot, mtSourceSplit, plotNames, 3, "WZsourcesplit" + extra, true);
	plotDataVSMC(mtTot, mtMismSplit, plotNames2, 3, "WZmismsplit" + extra, true);
	/*
	TH1D* mismairedMismSplit[3];
	for(unsigned m = 0; m < 3; ++m){
		mismairedMismSplit[m] = (TH1D*) mtSplit[1][m]->Clone();
	}
	plotNames2[0] = "total mispaired";
	plotDataVSMC( totalmisp , mismairedMismSplit, plotNames2, 3, "WZmispairingMETsplit" + extra, true);
	
	TH1D* correctpairMismSplit[3];
	for(unsigned m = 0; m < 3; ++m){
		correctpairMismSplit[m] = (TH1D*) mtSplit[0][m]->Clone();
	}
	plotNames2[0] = "total correct pairing";
	plotDataVSMC( totalCorr,  correctpairMismSplit, plotNames2, 3, "WZcorrectpairingMETsplit" + extra, true);

	std::vector<TH1D*> mtRes;
	mtRes.push_back(mtResSplit[0][0]);
	mtRes.push_back(mtResSplit[0][1]);
	mtRes.push_back(mtResSplit[0][2]);
	std::vector<TString> resnames = {"MET", "MET + #Delta MET", "MET - #Delta MET"};
	plotHist(mtRes, resnames, "correctpairing_metreseffect", true);

	mtRes.clear();
	mtRes.push_back(mtResSplit[1][0]);
	mtRes.push_back(mtResSplit[1][1]);
	mtRes.push_back(mtResSplit[1][2]);
	plotHist(mtRes, resnames, "mispairing_metreseffect", true);
	
	
	std::cout << "~~~~~~~~~~ Correct pairing ~~~~~~~~~~~~~ " << std::endl;
	cout << "mt > 160:  metup/met = " << mtRatio_corr[0][1]/mtRatio_corr[0][0] << endl;
	cout << "mt > 160:  metdown/met = " << mtRatio_corr[0][0]/mtRatio_corr[0][2] << endl;
	cout << "100 < mt < 160:  metup/met = " << mtRatio_corr[1][1]/mtRatio_corr[1][0] << endl;
	cout << "100 < mt < 160:  metdown/met = " << mtRatio_corr[1][0]/mtRatio_corr[1][2] << endl;
	
	std::cout << "~~~~~~~~~~ mispairing ~~~~~~~~~~~~~ " << std::endl;
	cout << "mt > 160:  metup/met = " << mtRatio_misp[0][1]/mtRatio_misp[0][0] << endl;
	cout << "mt > 160:  metdown/met = " << mtRatio_misp[0][0]/mtRatio_misp[0][2] << endl;
	cout << "100 < mt < 160:  metup/met = " << mtRatio_misp[1][1]/mtRatio_misp[1][0] << endl;
	cout << "100 < mt < 160:  metdown/met = " << mtRatio_misp[1][0]/mtRatio_misp[1][2] << endl;
	


	TH1D* mtAngMismSplit[3];
	for(unsigned m = 0; m < 3; ++m){
		mtAngMismSplit[m] = (TH1D*) mtAngSplit[0][m]->Clone();
		for(unsigned i = 1; i < 3; ++i){
			mtAngMismSplit[m]->Add(mtAngSplit[i][m]);
		}
	}
	TString plotNames3[4] = {"total WZ", "#Delta#Phi(MET) < 0.3", "0.3 < #Delta #Phi(MET) < 0.6", "#Delta #Phi(MET) > 0.6"};
	plotDataVSMC(mtTot, mtAngMismSplit, plotNames3, 3, "WZangmismsplit" + extra, true);

	TH1D* mismairedAngMismSplit[3];
	for(unsigned m = 0; m < 3; ++m){
		mismairedAngMismSplit[m] = (TH1D*) mtAngSplit[1][m]->Clone();
	}
	plotNames3[0] = "total mispaired";
	plotDataVSMC( totalmisp , mismairedAngMismSplit, plotNames3, 3, "WZmispairingMETPHIsplit" + extra, true);
	
	TH1D* correctpairAngMismSplit[3];
	for(unsigned m = 0; m < 3; ++m){
		correctpairAngMismSplit[m] = (TH1D*) mtAngSplit[0][m]->Clone();
	}
	plotNames3[0] = "total correct pairing";
	plotDataVSMC( totalCorr,  correctpairAngMismSplit, plotNames3, 3, "WZcorrectpairingMETPHIsplit" + extra, true);

	mtRes.clear();
	mtRes.push_back(mtResAngSplit[0][0]);
	mtRes.push_back(mtResAngSplit[0][1]);
	mtRes.push_back(mtResAngSplit[0][2]);
	resnames = {"#PhiMET", "#PhiMET + #Delta#PhiMET", "#PhiMET - #Delta#PhiMET"};
	plotHist(mtRes, resnames, "correctpairing_metangreseffect", true);
	

	mtRes.clear();
	mtRes.push_back(mtResAngSplit[1][0]);
	mtRes.push_back(mtResAngSplit[1][1]);
	mtRes.push_back(mtResAngSplit[1][2]);
	plotHist(mtRes, resnames, "mispairing_metangreseffect", true);
		

	mtRatioAng_corr[0][0] = mtRatio_corr[0][0];
	mtRatioAng_misp[0][0] = mtRatio_misp[0][0];
	mtRatioAng_corr[1][0] = mtRatio_corr[1][0];
	mtRatioAng_misp[1][0] = mtRatio_misp[1][0];
	std::cout << "~~~~~~~~~~ Correct pairing ~~~~~~~~~~~~~ " << std::endl;
	cout << "mt > 160:  metup/met = " << mtRatioAng_corr[0][1]/mtRatioAng_corr[0][0] << endl;
	cout << "mt > 160:  metdown/met = " << mtRatioAng_corr[0][0]/mtRatioAng_corr[0][2] << endl;
	cout << "100 < mt < 160:  metup/met = " << mtRatioAng_corr[1][1]/mtRatioAng_corr[1][0] << endl;
	cout << "100 < mt < 160:  metdown/met = " << mtRatioAng_corr[1][0]/mtRatioAng_corr[1][2] << endl;
	
	std::cout << "~~~~~~~~~~ mispairing ~~~~~~~~~~~~~ " << std::endl;
	cout << "mt > 160:  metup/met = " << mtRatioAng_misp[0][1]/mtRatioAng_misp[0][0] << endl;
	cout << "mt > 160:  metdown/met = " << mtRatioAng_misp[0][0]/mtRatioAng_misp[0][2] << endl;
	cout << "100 < mt < 160:  metup/met = " << mtRatioAng_misp[1][1]/mtRatioAng_misp[1][0] << endl;
	cout << "100 < mt < 160:  metdown/met = " << mtRatioAng_misp[1][0]/mtRatioAng_misp[1][2] << endl;

		plotHistRatio(leptonRes[1], leptonRes[0], "correctly paired", "mispaired", "LeptonResComp_inclusive",false, 0,0, true, false);

		plotHistRatio(DeltaPhiGen[0], DeltaPhiGen[1], "correctly paired", "mispaired", "DeltaPhiGenComp_inclusive",false, 0,0, false, false);
	*/
	/*
	//, SF = 0.7
	plotHistRatio(hists[0][0], hists[0][1], "mllmin01", "old", "MT_WZsampcomp" + extra, false, 0,0, true, true);
	plotHistRatio(hists[1][0], hists[1][1], "mllmin01, SF = 0.7", "old powheg", "Mmumu_WZsampcomp" + extra,  false, 0,0, false, true);
	plotHistRatio(hists[1][0], hists[1][1], "mllmin01, SF = 0.7", "old powheg", "Mmumu_WZsampcomp_log" + extra,  false, 0,0, true, true);
	*/
	//plotHist(hists[2][0], "lowGenMll");
}



int main(int argc, char* argv[]){
	TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;
	testtree.Loop();
	rootapp->Run();
    return 0;
}


