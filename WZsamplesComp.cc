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

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"


void trilTree::Loop(){
	//Set plotting style
	setTDRStyle();
	//gROOT->SetBatch(kTRUE);
	//Define list of samples
	const unsigned nSamples = 4;
	//0.696758
	//const TString fileList[nSamples] = {"WZTo3LNu_mllmin01.root", "WZTo3LNu_powheg.root"};
	const TString fileList[nSamples] = {"WZTo3LNu_mllmin01.root", "SMS-TChiSlepSnu_x0p5.root", "SMS-TChiSlepSnu_x0p05.root", "SMS-TChiWZ.root"};
	//const TString fileList[nSamples] = {"WZTo3LNu_mllmin01.root", "test.root"};
	const double xSections[nSamples] = {58.59*0.696758, 1 , 1, 1};
	const TString names[nSamples] = {"WZTo3LNu_mllmin01", "SMS-TChiSlepSnu_x0p5", "SMS-TChiSlepSnu_x0p05", "SMS-TChiWZ"};

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
		Init(inputTree[sam], true, true);
	}
	
	//Tweakable options////////////////////////////////////////////////////
	const bool TestRun = false;	//Break after a few events
	const double DataLuminosity = 35.9;    //units of fb^{-1}
	const TString extra = "_BDT_noCuts";	//for plot file names
	//////////////////////////
	
	//Make histograms containing kinematic distributions
	const unsigned nDist = 3;
	const TString distNames[nDist] = {"mt", "mtnew", "mtunamb"};
	const TString xAxes[nDist] = {"M_{T}", "M_{T}^{BDT}", "M_{T}^{unambig}"};
	const TString units[nDist] = {"GeV", "GeV", "GeV"};
	const double histMin[nDist] = {0, 0, 0};
	const double histMax[nDist] = {300, 300, 300};
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

	
	TH1D* mtHist[2][nSamples];
	const TString methods[2] = {"best Z", "BDT"};
	for(unsigned m = 0; m < 2; ++m){
		for(unsigned sam = 0; sam < nSamples; ++sam){
			mtHist[m][sam] = new TH1D("mtHist" + methods[m] + names[sam], "mtHist" + methods[m] + names[sam] + ";M_{T} (GeV); Events", 100, 0, 300);
			mtHist[m][sam]->Sumw2();
		}
	}


	TMVA::Reader *readerWZ;
    Float_t mos, mother, dPhiMet, mt, wpt, lpt, deltaRos, deltaPtos, deltaEos, motherAndMet;
    readerWZ = new TMVA::Reader( "!Color:!Silent" );
	readerWZ->AddVariable("mos", &mos);
	readerWZ->AddVariable("mother", &mother);
	//readerWZ->AddVariable("dPhiMet", &dPhiMet);
	readerWZ->AddVariable("mt", &mt);
	readerWZ->AddVariable("deltaRos", &deltaRos);
	readerWZ->AddVariable("wpt", &wpt);
	//readerWZ->AddVariable("deltaBosonPt", &deltaBosonPt);
	readerWZ->AddVariable("lpt", &lpt);
	readerWZ->AddVariable("deltaPtos", &deltaPtos);
	readerWZ->AddVariable("deltaEos", &deltaEos);
	readerWZ->AddVariable("motherAndMet", &motherAndMet);

   	readerWZ->BookMVA( "BDT method", "/home/willem/Work/AnalysisCode/dataset/weights/TMVAClassification_BDT.weights.xml"); 
	//readerWZ->BookMVA( "MLP method", "/home/willem/Work/AnalysisCode/dataset/weights/TMVAClassification_MLP.weights.xml"); 
   	//readerWZ->BookMVA( "BDTG method", "/home/willem/Work/AnalysisCode/dataset/weights/TMVAClassification_BDTG.weights.xml"); 
   	//readerWZ->BookMVA( "BDTB method", "/home/willem/Work/AnalysisCode/dataset/weights/TMVAClassification_BDTB.weights.xml"); 

    Double_t scale[nSamples];
	//Loop over all samples
	for(unsigned sam = 0; sam < nSamples; ++sam){
		Long64_t nEntries = inputTree[sam]->GetEntries();
		scale[sam] = xSections[sam]*DataLuminosity*1000/(hcounter[sam]);
		std::cout<<"Entries in "<< fileList[sam] <<" "<<nEntries<<std::endl;

		float progress = 0.0; //Progress bar 
	    for(Long64_t it = 0; it < nEntries/5; ++it){
			//if (it%10000 == 0) cout<<'.'<<flush;
			//print loading screen with progress bar
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
			TLorentzVector metV;
			metV.SetPtEtaPhiE(_met, 0, _met_phi, _met);
			/*
			if(_met < 50) continue;
			if(!ptCuts_hnl(ind,lCount)) continue;
			if(!vetoLowMll(12)) continue;
			*/		
			TLorentzVector Wlep;
			Wlep.SetPtEtaPhiE(_lPt[lw], _lEta[lw], _lPhi[lw], _lE[lw]);
			if(_flavors[lw] == _flavors[mllI[0]]){
				//old algo
				hists[0][sam]->Fill(transmass(Wlep,metV));
				for(unsigned c = 0; c < 100; ++c){
					if(transmass(Wlep,metV) > c*3){
						mtHist[0][sam]->Fill(c*3, scal);
					}
				}
				//new algo
				unsigned indOther;
				unsigned indZ;
				for(unsigned l = 0; l < 2; ++l){
					if(_charges[mllI[l]] == _charges[lw]) indOther = mllI[l];	
					else indZ = mllI[l];
				}
				double mvaBad, mvaGood;
				TLorentzVector Zlep, OtherLep;
				Zlep.SetPtEtaPhiE(_lPt[indZ], _lEta[indZ], _lPhi[indZ], _lE[indZ]);
				OtherLep.SetPtEtaPhiE(_lPt[indOther], _lEta[indOther], _lPhi[indOther], _lE[indOther]);
				mos = (Wlep + OtherLep).M();
				mother = mll;
				//dPhiMet = fabs(Wlep.DeltaPhi(metV));
				mt = transmass(Wlep,metV);
				wpt = (Wlep + metV).Pt();
				lpt = Wlep.Pt();
				deltaRos = Wlep.DeltaR(Zlep);

				deltaPtos = fabs(Wlep.Pt() - Zlep.Pt());
				deltaEos = fabs(Wlep.E() - Zlep.E());
				motherAndMet = transmass(OtherLep + Zlep, metV);
				mvaGood = readerWZ->EvaluateMVA( "BDT method");
				//mvaGood = readerWZ->EvaluateMVA( "MLP method");

				mos = mll;
				mother = (Wlep + OtherLep).M();
				//dPhiMet = fabs(OtherLep.DeltaPhi(metV));
				mt = transmass(OtherLep,metV);
				wpt = (OtherLep + metV).Pt();
				lpt = OtherLep.Pt();
				deltaRos = OtherLep.DeltaR(Zlep);
				deltaPtos = fabs(OtherLep.Pt() - Zlep.Pt());
				deltaEos = fabs(OtherLep.E() - Zlep.E());
				motherAndMet = transmass(Wlep + Zlep, metV);
				mvaBad = readerWZ->EvaluateMVA( "BDT method");
				//mvaBad = readerWZ->EvaluateMVA( "MLP method");
				//std::cout << "mvaGood = " << mvaGood << "		 	mvaBad = " << mvaBad << std::endl;
				if(mvaBad > mvaGood) {
					//std::cout << "change W lepton" << std::endl;
					lw = indOther;
				}
				Wlep.SetPtEtaPhiE(_lPt[lw], _lEta[lw], _lPhi[lw], _lE[lw]);
				hists[1][sam]->Fill(transmass(Wlep,metV), scal);


				for(unsigned c = 0; c < 100; ++c){
					if(transmass(Wlep,metV) > c*3){
						mtHist[1][sam]->Fill(c*3, scal);
					}
				}
			}
			else{
				hists[2][sam]->Fill(transmass(Wlep,metV), scal);

			}
			/*
			const double mw = 80.385;
			double minDiff = 99999.;
			double bestSol = 0.;
			unsigned nuI = 99;
			TLorentzVector* nu = new TLorentzVector[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				//if(lepV[l].DeltaPhi(metV) < 0.6) continue;
				double m2 = 0.5*mw*mw + _lPt[ind[l]]*_met;
				double solplus = (m2/(_lPt[ind[l]]*_lPt[ind[l]]))*(lepV[l].Pz() + fabs(lepV[l].P())*sqrt(1 - (_met*_met*_lPt[ind[l]]*_lPt[ind[l]])/(m2*m2) ) );
				double solmin = (m2/(_lPt[ind[l]]*_lPt[ind[l]]))*(lepV[l].Pz() - fabs(lepV[l].P())*sqrt(1 - (_met*_met*_lPt[ind[l]]*_lPt[ind[l]])/(m2*m2) ) );
				double nupx = _met*cos(_met_phi);
				double nupy = _met*sin(_met_phi);
				TLorentzVector vecplus, vecmin;
				vecplus.SetPxPyPzE(nupx, nupy, solplus, sqrt(_met*_met + solplus*solplus));	
				vecmin.SetPxPyPzE(nupx, nupy, solmin, sqrt(_met*_met + solmin*solmin));
				if(lepV[l].DeltaPhi(vecplus) > lepV[l].DeltaPhi(vecmin)){
					nu[l] = vecplus;
				} else{
					nu[l] = vecmin;
				}
			}	
			//Check which lepton has the "single" charge in the event, this one can't come from the W
			//for(unsigned l = 0; l < lCount; ++l){		



			double propagator1 = 1/(91.1876*2.4952*91.1876*2.4952 + ( (lepV[0] + lepV[1]).M()*(lepV[0] + lepV[1]).M() - 91.1876*91.1876)*( (lepV[0] + lepV[1]).M()*(lepV[0] + lepV[1]).M() - 91.1876*91.1876));
			double propagator2 = 1/(80.385*80.385*2.085*2.085 + ( (nu[2] + lepV[2]).M()*(nu[2] + lepV[2]).M() - 80.385*80.385)*( (nu[2] + lepV[2]).M()*(nu[2] + lepV[2]).M() - 80.385*80.385));
			double propagator = propagator1*propagator2;
			unsigned lwa = ind[2];
			propagator1 = 1/(91.1876*2.4952*91.1876*2.4952 + ( (lepV[0] + lepV[2]).M()*(lepV[0] + lepV[2]).M() - 91.1876*91.1876)*( (lepV[0] + lepV[2]).M()*(lepV[0] + lepV[2]).M() - 91.1876*91.1876));
			propagator2 = 1/(80.385*80.385*2.085*2.085 + ( (nu[1] + lepV[1]).M()*(nu[1] + lepV[1]).M() - 80.385*80.385)*( (nu[1] + lepV[1]).M()*(nu[1] + lepV[1]).M() - 80.385*80.385));
			double propagatortemp = propagator1*propagator2;
			if(propagatortemp > propagator){
				propagator = propagatortemp;
				lwa = ind[1];
			}
			propagator1 = 1/(91.1876*2.4952*91.1876*2.4952 + ( (lepV[1] + lepV[2]).M()*(lepV[1] + lepV[2]).M() - 91.1876*91.1876)*( (lepV[1] + lepV[2]).M()*(lepV[1] + lepV[2]).M() - 91.1876*91.1876));
			propagator2 = 1/(80.385*80.385*2.085*2.085 + ( (nu[0] + lepV[0]).M()*(nu[0] + lepV[0]).M() - 80.385*80.385)*( (nu[0] + lepV[0]).M()*(nu[0] + lepV[0]).M() - 80.385*80.385));
		    propagatortemp = propagator1*propagator2;
			if(propagatortemp > propagator){
				propagator = propagatortemp;
				lwa = ind[0];
			}
			*/
			/*
			TLorentzVector metV;
			metV.SetPtEtaPhiE(_met, 0, _met_phi, _met);
			double propagator1 = 1/(91.1876*2.4952*91.1876*2.4952 + ( (lepV[0] + lepV[1]).M()*(lepV[0] + lepV[1]).M() - 91.1876*91.1876)*( (lepV[0] + lepV[1]).M()*(lepV[0] + lepV[1]).M() - 91.1876*91.1876));
			double propagator2 = 1/(80.385*80.385*2.085*2.085 + ( transmass(lepV[2], metV)*transmass(lepV[2], metV) - 80.385*80.385)*( transmass(lepV[2], metV)*transmass(lepV[2], metV) - 80.385*80.385));
			double propagator = propagator1*propagator2;
			unsigned lwa = ind[2];
			propagator1 = 1/(91.1876*2.4952*91.1876*2.4952 + ( (lepV[0] + lepV[2]).M()*(lepV[0] + lepV[2]).M() - 91.1876*91.1876)*( (lepV[0] + lepV[2]).M()*(lepV[0] + lepV[2]).M() - 91.1876*91.1876));
			propagator2 = 1/(80.385*80.385*2.085*2.085 + ( transmass(lepV[1], metV)*transmass(lepV[1], metV) - 80.385*80.385)*( transmass(lepV[1], metV)*transmass(lepV[1], metV) - 80.385*80.385));
			double propagatortemp = propagator1*propagator2;
			if(propagatortemp > propagator){
				propagator = propagatortemp;
				lwa = ind[1];
			}
			propagator1 = 1/(91.1876*2.4952*91.1876*2.4952 + ( (lepV[1] + lepV[2]).M()*(lepV[1] + lepV[2]).M() - 91.1876*91.1876)*( (lepV[1] + lepV[2]).M()*(lepV[1] + lepV[2]).M() - 91.1876*91.1876));
			propagator2 = 1/(80.385*80.385*2.085*2.085 + ( transmass(lepV[0], metV)*transmass(lepV[0], metV) - 80.385*80.385)*( transmass(lepV[0], metV)*transmass(lepV[0], metV) - 80.385*80.385));
		    propagatortemp = propagator1*propagator2;
			if(propagatortemp > propagator){
				propagator = propagatortemp;
				lwa = ind[0];
			}
			*/
			//if(_flavors[lw] != _flavors[mllI[0]]) continue;
			/*
			if(!ptCuts_hnl(ind,lCount)) continue;
			if(!vetoLowMll(12)) continue;
			*/
			/*
			TLorentzVector Wlep, WlepNew;
			Wlep.SetPtEtaPhiE(_lPt[lw], _lEta[lw], _lPhi[lw], _lE[lw]);
			WlepNew.SetPtEtaPhiE(_lPt[lwa], _lEta[lwa], _lPhi[lwa], _lE[lwa]);
			double fill[nDist] = {transmass(Wlep,metV), transmass(WlepNew,metV), transmass(Wlep,metV)};
			for(unsigned dist = 0; dist < nDist; ++dist){
				if(_flavors[lw] == _flavors[mllI[0]] && dist == 2) continue;
				else if(_flavors[lw] != _flavors[mllI[0]] && dist != 2) continue;
				hists[dist][sam]->Fill(fill[dist], scal);
			}
			*/
			//if((cat == 0 || cat == 2 || cat == 4) && fabs(mll - 91) < 15) continue; //consider onZ events for WZ CR

			/*
			//determine the index of the W lepton
			unsigned lw = 9999;
			for(unsigned l = 0; l < lCount; ++l){
				if(ind[l] != mllI[0] && ind[l] != mllI[1]){
					lw = ind[l];
				}
			}
			//if(_flavors[lw] ==  _flavors[mllI[0]]) continue;
			if(_flavors[lw] == 0) continue;
			*/
			/*
			if(_gen_nL != 3) continue;
			unsigned nMu = 0;
			unsigned nEle = 0;
			unsigned* genI = new unsigned[_gen_nL];
			unsigned nLight;
			for(unsigned l = 0; l < _gen_nL; ++l){
				//if(_gen_lPt[l] < 0.01) continue;
				if(_gen_flavors[l] == 0) ++nEle;
				if(_gen_flavors[l] == 1) ++nMu;
			}
			if(nMu != 2) continue;
			if(nEle != 1) continue;
			unsigned mllI[2] = {99, 99};
			unsigned mllC = 0;
			unsigned lw = 0;
			for(unsigned l = 0; l < _gen_nL; ++l){
				//if(_gen_lPt[l] < 0.01) continue;
				if(_gen_flavors[l] == 1){
					mllI[mllC] = l;
					++mllC;
				} else if(_gen_flavors[l] == 0){
					lw = l;
				}
			}
			
			TLorentzVector lzV[2];
			for(unsigned l = 0; l < 2; ++l) lzV[l].SetPtEtaPhiE(_gen_lPt[mllI[l]], _gen_lEta[mllI[l]], _gen_lPhi[mllI[l]], _gen_lE[mllI[l]]);
			cout << "_gen_lPt[mllI[0]] = " << _gen_lPt[mllI[0]] << endl;
			double mll = (lzV[0] + lzV[1]).M();
			TLorentzVector Wlep, metV;
			Wlep.SetPtEtaPhiE(_gen_lPt[lw], _gen_lEta[lw], _gen_lPhi[lw], _gen_lE[lw]);
			metV.SetPtEtaPhiE(_met, 0, _met_phi, _met);
			double fill[nDist] = {transmass(Wlep,metV), mll, mll};
			for(unsigned dist = 0; dist < nDist; ++dist){
				hists[dist][sam]->Fill(fill[dist], scal);
			}
			*/
		}
	}
	//, SF = 0.7
	//plotHistRatio(hists[0][0], hists[1][0], "mt old", "mt new", "MT_newalgo_neutrinoZ" + extra, false, 0,0, true, true);
	
	//plotHistRatio(hists[1][0], hists[1][1], "mllmin01, SF = 0.7", "old powheg", "Mmumu_WZsampcomp" + extra,  false, 0,0, false, true);
	//plotHistRatio(hists[1][0], hists[1][1], "mllmin01, SF = 0.7", "old powheg", "Mmumu_WZsampcomp_log" + extra,  false, 0,0, true, true);
	
	//plotHist(hists[2][0], "lowGenMll");
	for(unsigned sam = 0; sam < nSamples; ++sam){
		std::vector<TH1D*> distVec;
		std::vector<TString> distnames = {"bestZ", "BDT", "unambiguous"};
		for(unsigned i = 0; i < 3; ++i){
			distVec.push_back(hists[i][sam]);
		}
		plotHist(distVec, distnames, "WZalgo/MT_newalgoComp" + names[sam] + extra, true);
	}


	//make S/sqrt(S + B) histograms as a function of MT cut 
	for(unsigned s = 0; s < nSamples - 1; ++s){
		TH1D* signif[2];
		for(unsigned m = 0; m < 2; ++m){
			signif[m] = (TH1D*) mtHist[m][0]->Clone();
			for(unsigned b = 1; b < signif[m]->GetNbinsX() + 1; ++b){
				signif[m]->SetBinContent(b, mtHist[m][s + 1]->GetBinContent(b)/( sqrt( mtHist[m][s + 1]->GetBinContent(b) + mtHist[m][0]->GetBinContent(b) ) )  );
				signif[m]->SetBinError(b, 0);
			}
			plotHist(signif[m], "WZalgo/S_over_sqrt_SplusB_" + methods[m] + names[s+1]);
		}
		signif[1]->Divide(signif[0]);
		plotHist(signif[1], "WZalgo/S_over_sqrt_SplusB_RATIO_" + names[s+1]);
	}
}



int main(int argc, char* argv[]){
	TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;
	testtree.Loop();
	rootapp->Run();
    return 0;
}


