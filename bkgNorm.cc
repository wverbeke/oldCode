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
#include <algorithm>

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
	//Define samples to loop over, their names and cross-sections
	const unsigned nSamples = 32;
	const unsigned nSamples_eff = 6;
	const TString fileList[nSamples] = {"data_combined_trilepton.root", "ZZTo4L.root",  "GluGluToZZTo4mu.root", "GluGluToZZTo4e.root", "GluGluToZZTo4tau.root", "GluGluToZZTo2e2mu.root", "GluGluToZZTo2e2tau.root", "GluGluToZZTo2mu2tau.root", 


"VHToNonbb.root", "GluGluHToZZTo4L_M125.root", "VBF_HToZZTo4L_M125.root", "WWG.root","WWW.root", "WWZ.root", "WWTo2L2Nu_DoubleScattering.root", "WWTo2L2Nu.root",  "ZZZ.root", "WZTo3LNu_mllmin01.root", "TTGJets.root","ZGTo2LG.root", "WGToLNuG.root", "TGJets.root", "TTJets_DiLept.root", "TTJets_SingleLeptFromTbar.root", "TTJets_SingleLeptFromT.root", "DYJetsToLL_M10to50.root", "DYJetsToLL_M50.root", "ttHToNonbb.root", "TTWJetsToLNu.root", "TTZToLLNuNu.root", "TTZToLL_M1to10.root",  "TTTT.root"};	
	/*
	const double glugluToZZkFactor = 2.1;
	const double WZSF = 0.652; //0.655
	const double ZZSF = 1.029; //1.032
	const double XgammaSF = 0.950; //0.948
	*/
	const double glugluToZZkFactor = 2.1; //1.7
	const double WZSF = 1;
	const double ZZSF = 1;
	const double XgammaSF = 1;
	const double xSections[nSamples - 1] = {1.256*ZZSF, 0.00159*glugluToZZkFactor*ZZSF, 0.00159*glugluToZZkFactor*ZZSF, 0.00159*glugluToZZkFactor*ZZSF, 0.00319*glugluToZZkFactor*ZZSF, 0.00319*glugluToZZkFactor*ZZSF, 0.00319*glugluToZZkFactor*ZZSF,

 0.9561, 0.01212, 0.001034,  0.2147, 0.2086, 0.1651, 0.1729,   12.178, 0.01398, 58.59*WZSF, 3.697, 123.9*XgammaSF, 489*XgammaSF,  2.967, 87.315, 182.175, 182.175, 18610, 1921.8*3, 0.215, 0.2043, 0.2529, 0.0493, 0.009103};
	const TString names[nSamples] = {"data", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H",  "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "triboson", "triboson", "triboson", "triboson", "triboson", "triboson", "WZ", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X"};
	
	//Read Trees from ROOT files
	TFile* hfile[nSamples];
	TTree* inputTree[nSamples];
	double hcounter[nSamples];
	for(unsigned sam = 0; sam < nSamples; ++sam){
		cout << "name " << names[sam] << endl;
		hfile[sam] = new TFile("../data_april17/"+fileList[sam],"read");
       	hfile[sam]->cd("FakeElectrons");
		//Determine hcounter for cross section scaling
		TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
		_hCounter->Read("hCounter");
		hcounter[sam] = _hCounter->GetBinContent(1);
       	inputTree[sam] = (TTree*) (hfile[sam]->Get("FakeElectrons/fakeTree"));
		Init(inputTree[sam], false, sam > 0);
	}
	//Read analysis scale factors
	readSF(true);
	//Tweakable options////////////////////////////////////////////////////
	const bool TestRun = false;	//Break after a few events
	const double DataLuminosity = 35.867;    //units of fb^{-1}
	const TString extra = "";	//for plot file names
	//const bool plotKinematics = false;
	//////////////////////////
	TH1D* yields[nSamples_eff + 1];	//Total yields in every CR, to be used in simultaneous fit
	const TString eff_names[nSamples_eff + 1] = {"data", "ZZ/H", "triboson", "WZ", "X + #gamma", "TT/T + X",  "non-prompt"}; //X + #gamma
	const unsigned nCont = 3;	//Number of control regions in simultaneous fit
	//const TString contNames[nCont] = {"WZ", "ZZ", "Xgamma", "DrellYan", "TT"};
	//const unsigned nDist = 57;  //Number of distributions to plot	
	//TH1D* histos[nCont][nDist][nSamples_eff + 1];	//Kinematic distributions to plot
	const unsigned nUnc = 9;//number of shape uncertainties
	const TString uncNames[nUnc] = {"jec", "metUnclustered", "pdf", "scale", "pu", "btagSF", "trigger", "ideff", "fakeEWK"};  //shape uncertainty names
	//Total yield shape uncertainty histograms
	TH1D* yieldsDown[nUnc][nSamples_eff + 1];
	TH1D* yieldsUp[nUnc][nSamples_eff + 1];
	TH1D* yieldsPdfVar[100][nSamples_eff + 1];
		
	//Initialize histogram containing total yields
	for(unsigned effsam = 0; effsam < nSamples_eff + 1; ++effsam){
		yields[effsam] = new TH1D("yields" + eff_names[effsam], "yields" + eff_names[effsam] + ";control region; Events/control region", 3, 0.5, 3 + 0.5);
		for(unsigned unc = 0; unc < nUnc; ++unc){
			yieldsDown[unc][effsam] = new TH1D("yields" + eff_names[effsam] + uncNames[unc] + "Down", "yields" + eff_names[effsam] + uncNames[unc] + "Down;control region; Events/control region", nCont, 0.5, nCont + 0.5);
			yieldsUp[unc][effsam] = new TH1D("yields" + eff_names[effsam] + uncNames[unc] + "Up", "yields" + eff_names[effsam] + uncNames[unc] + "Up;control region; Events/control region", nCont, 0.5, nCont + 0.5);
		}	
		for(unsigned pdf = 0; pdf < 100; ++pdf){
			yieldsPdfVar[pdf][effsam] = new TH1D("yields" + eff_names[effsam]  + "_pdf" + std::to_string(pdf), "yields" + eff_names[effsam]  + "_pdf" + std::to_string(pdf) + "; control region; Events/control region", nCont, 0.5, nCont + 0.5);
		}
	}
	
    Double_t scale[nSamples -1];
	//set to check which data events have already been processed
	std::set<std::tuple<unsigned long, unsigned long, unsigned long> > usedEvents; //runNb, lumiBlock, eventNb

	//root file used for memory control, i.e. temporarily storing the histograms not to fill up the ram
	
	
	//Loop over all samples
	for(unsigned sam = 0, effsam = 0; sam < nSamples; ++sam, ++effsam){
		if(sam != 0){
   			if(names[sam] == names[sam -1]) --effsam;
    	}		

		cout << eff_names[effsam] << endl;
		//if(fileList[sam] == "TTJets_DiLept.root" || fileList[sam] == "DYJetsToLL_M10to50.root" || fileList[sam] == "DYJetsToLL_M50.root"  || fileList[sam] == "TTJets_SingleLeptFromTbar.root"  || fileList[sam] == "TTJets_SingleLeptFromT.root" ) continue;
    	Long64_t nEntries = inputTree[sam]->GetEntries();
    	if(sam > 0){
   			scale[sam -1] = xSections[sam -1]*DataLuminosity*1000/(hcounter[sam]);
    	}

		std::cout<<"Entries in "<< fileList[sam] <<" "<<nEntries<<std::endl;
		cout << effsam << endl;
        for(Long64_t it = 0; it < nEntries; ++it){
        	inputTree[sam]->GetEntry(it);
        	if (it%10000 == 0) cout<<'.'<<flush;
        	if(TestRun && it > 10000) break;
        	double scal;
        	if(effsam == 0) scal = 1;
        	else{
				scal = scale[sam-1]*_weight;
			}		
			
			//Apply HNL SELECTION
			cutBased();
			//Baseline event selection:
			if(!baseline(true, true, false, false)) continue;
			if(effsam > 0 && nBJets(true, false, 1) != 0) continue;
			else if(effsam == 0 && nBJets(true, false, 0) != 0) continue;
			//Check if data events were used before
			if(effsam == 0){
				auto event = usedEvents.find(std::make_tuple(_eventNb, _lumiBlock, _runNb));
				if(event != usedEvents.end()) continue;
				usedEvents.insert(std::make_tuple(_eventNb, _lumiBlock, _runNb));
			}
			//Categorize according to the number of leptons and flavors
			unsigned* ind = new unsigned[_nL];
			//Order leptons by Pt
			unsigned lCount = lepOrder(ind, 3, true, true);
			if(lCount < 3 || lCount > 4) continue;
			//MC prompt matching
			if(effsam != 0){
				bool promptfail = false;
				for(unsigned l = 0; l < lCount; ++l){	
					if(_origin[ind[l]] != 0){
						promptfail = true;
						break;
					}
				}
				if(promptfail) continue;
			}	
			if(fileList[sam] == "WGToLNuG.root"){
				bool promptfail =true;
				for(unsigned l = 0; l < lCount; ++l){	
					if(_pdgmc[ind[l]] == 22 && _originPhot[ind[l]] == 0){
						promptfail = false;
						break;
					}
				}
				if(promptfail) continue;
			}
			
			//Require 3 leptons to be tight in data and MC, and determine nonPrompt bkg in data
			unsigned nTight = tightCount(ind, lCount);
			bool tightFail = nTight < lCount;
			//index used to fill events, needed to separate fakes from data
			unsigned fill = effsam;
			//Calculate conePt for every lepton
			double* conePt = new double[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				conePt[l] = _lPt[ind[l]]*std::max(1., 1 + (_isolation[ind[l]] - 0.1));
			}
			if(!ptCuts_hnl(ind,lCount)) continue;
			//Apply FR maps to data control region
			if(tightFail){
				//fakes go in different histogram
				fill = nSamples_eff;
				//MC fake subtraction
				if(effsam != 0) scal *= -1.;
				//Apply FR maps
				scal*= fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[0], lCount);
			} 
			if(conePt[2] >= 10){
				if(fileList[sam] == "TTJets_DiLept.root" || fileList[sam] == "DYJetsToLL_M10to50.root" || fileList[sam] == "DYJetsToLL_M50.root"  || fileList[sam] == "TTJets_SingleLeptFromTbar.root"  || fileList[sam] == "TTJets_SingleLeptFromT.root" ) continue;
			} else if(conePt[2] < 10){
				if(fileList[sam] == "TTJets_DiLept.root" || fileList[sam] == "DYJetsToLL_M10to50.root" || fileList[sam] == "DYJetsToLL_M50.root"  || fileList[sam] == "TTJets_SingleLeptFromTbar.root"  || fileList[sam] == "TTJets_SingleLeptFromT.root" ){
					double maxMll = 0.;
					for(unsigned l = 0; l < _gen_nL -1; ++l){
						TLorentzVector lep1;
						lep1.SetPtEtaPhiE(_gen_lPt[l], _gen_lEta[l], _gen_lPhi[l], _gen_lE[l]);
						for(unsigned k = l + 1; k < _gen_nL; ++k){
							TLorentzVector lep2;
							lep2.SetPtEtaPhiE(_gen_lPt[k], _gen_lEta[k], _gen_lPhi[k], _gen_lE[k]);
							if( (lep1 + lep2).M() > maxMll) maxMll = (lep1 + lep2).M();
						}
					}
					if(maxMll > 30){
						bool promptfail = false;
						for(unsigned l = 0; l < lCount; ++l){	
							if(_pdgmc[ind[l]] == 22 && _originPhot[ind[l]] == 0){
								promptfail = true;
								break;
							}
						}
						if(promptfail) continue;	
					}
				}
				if(fileList[sam] == "TTGJets.root" || fileList[sam] == "ZGTo2LG.root"){
					bool promptfail =true;
					for(unsigned l = 0; l < lCount; ++l){	
						//cout << _origin[ind[l]] << endl;
						if(_pdgmc[ind[l]] == 22 && _originPhot[ind[l]] == 0){
							promptfail = false;
							break;
						}
					}
					if(promptfail) continue;
				}
			}
			//Apply triggers
			bool trigPass[4];
			trigPass[0] = _lowM_trigger_eee || _lowM_trigger_all;
			trigPass[1] = _lowM_trigger_mee || _lowM_trigger_all;
			trigPass[2] = _lowM_trigger_mme || _lowM_trigger_all;
			trigPass[3] = _lowM_trigger_mmm || _lowM_trigger_all;
			if(!trigPass[tril_flavorComb(ind,_flavors, lCount)]) continue;
			//Apply all efficiency and reweighing SF to the simulation
			if(effsam != 0){
				scal*=getEventSF(ind, lCount, true);
			}
			//determine which leptons will be used for the calculation of mll
			TLorentzVector* lepV = new TLorentzVector[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				//lepV[l].SetPtEtaPhiE(conePt[l], _lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]);
				lepV[l].SetPtEtaPhiE(conePt[l], _lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]*(conePt[l]/_lPt[ind[l]]));
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
				for(unsigned l = 0; l < 2; ++l) lzV[l].SetPtEtaPhiE(_lPt[mllI[l]]*std::max(1., 1 + (_isolation[mllI[l]] - 0.1)), _lEta[mllI[l]], _lPhi[mllI[l]], _lE[mllI[l]]*std::max(1., 1 + (_isolation[mllI[l]] - 0.1)));
				mll = (lzV[0] + lzV[1]).M();
			}
			//Calculate lepton system vector
			TLorentzVector lepSyst;
			for(int l = 0; l < lCount; ++l) lepSyst += lepV[l];
			//Determine index of control region 0 = WZ, 1 = ZZ, 2 = conversions
			unsigned controlR = hnl::controlRegion(ind, _flavors, _charges, lCount, lepSyst.M(), mll);
			if(controlR == 999) continue;
			if(controlR == 0){ //WZ selection
				if(conePt[0] < 25.) continue;
				if(conePt[1] < 15.) continue;
				if(conePt[2] < 10.) continue;
			} else if(controlR == 1){ //ZZ selection
				if(mllI[1] == 99) continue;
				if(_flavors[mllI[1]] != _flavors[mllI[0]]) continue;
				unsigned unusedI[2] = {99, 99};
				unsigned unusedC = 0;
				for(unsigned l = 0; l < lCount; ++l){
					if(ind[l] == mllI[0] || ind[l] == mllI[1]) continue;
					unusedI[unusedC] = ind[l];
					++unusedC;
				}
				TLorentzVector lzV[2];
				for(unsigned l = 0; l < 2; ++l) lzV[l].SetPtEtaPhiE(_lPt[unusedI[l]]*std::max(1., 1 + (_isolation[unusedI[l]] - 0.1)), _lEta[unusedI[l]], _lPhi[unusedI[l]], _lE[unusedI[l]]*std::max(1., 1 + (_isolation[unusedI[l]] - 0.1)));
				if( fabs( (lzV[0] + lzV[1]).M() - 91.) > 15) continue;
				if(_flavors[unusedI[0]] != _flavors[unusedI[1]]) continue;
				if( fabs( (lzV[0] + lzV[1]).M() - 91.) > 15) continue;
			} else if(controlR == 2){ //Xgamma selection
				if(conePt[2] < 10.) continue;
			} else if(controlR == 3){
				//put MT cut later since it depends on met uncertainties
			}

			
			//Scale weights for glugluToZZ
			if(fileList[sam].Contains("GluGluToZZTo")){
				for(unsigned pdf = 0; pdf < 110; ++pdf){
					_scaleWeight[pdf] = 1;
				}
			}
			//Nominal yields
			if(nBJets(true, false, 0) == 0 && (controlR > 0 || _met > 50)) yields[fill]->Fill(controlR + 1, scal);
			//vary fake EWK subtraction
			if(fill == nSamples_eff && nBJets(true, false, 0) == 0 && (controlR > 0 || _met > 50)){
				yieldsDown[8][fill]->Fill(controlR + 1, scal*(fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[1], lCount)/fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[0], lCount) ) );
				yieldsUp[8][fill]->Fill(controlR + 1, scal*(fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[2], lCount)/fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[0], lCount) ) );
			}
			//vary JEC down
			if(nBJets(true, false, 1) == 0 && (controlR > 0 || _metJECDown > 50) )  yieldsDown[0][fill]->Fill(controlR + 1, scal);
			//vary JEC up
			if(nBJets(true, false, 2) == 0 && (controlR > 0 || _metJECUp > 50) )  yieldsUp[0][fill]->Fill(controlR + 1, scal);
			//nominal b-veto 
			if(nBJets(true, false, 0) != 0) continue;
			//vary unclustered met down
			if(controlR > 0 || _metOtherDown > 50) yieldsDown[1][fill]->Fill(controlR + 1, scal);
			//vary unclustered met up
			if(controlR > 0 || _metOtherUp > 50) yieldsUp[1][fill]->Fill(controlR + 1, scal);
			//nominal met cut
			if(controlR == 0 && _met < 50) continue;
			//vary scale down
			yieldsDown[3][fill]->Fill(controlR + 1, scal*_scaleWeight[8]);
			//vary scale up
			yieldsUp[3][fill]->Fill(controlR + 1, scal*_scaleWeight[4]); 	
			//vary pu down
			yieldsDown[4][fill]->Fill(controlR + 1, (scal/PUweights[0]->GetBinContent(PUweights[0]->FindBin( std::min(_n_trueInteractions, 49) )))*PUweights[1]->GetBinContent(PUweights[1]->FindBin(std::min(_n_trueInteractions, 49)  )) );
			//vary pu up
			yieldsUp[4][fill]->Fill(controlR + 1, (scal/PUweights[0]->GetBinContent(PUweights[0]->FindBin( std::min(_n_trueInteractions, 49) )))*PUweights[2]->GetBinContent(PUweights[2]->FindBin(std::min(_n_trueInteractions, 49)  )) );
			//vary btag SF down	
			yieldsDown[5][fill]->Fill(controlR + 1,  (scal/bTagSF(true, 0))*bTagSF(true, 1) );
			//vary btag SF up
			yieldsUp[5][fill]->Fill(controlR + 1,  (scal/bTagSF(true, 0))*bTagSF(true, 2) );
			//trigger unc 
			if(conePt[0] > 30){
				yieldsDown[6][fill]->Fill(controlR + 1, scal*0.98);
				yieldsUp[6][fill]->Fill(controlR + 1, scal*1.02);
			} else{
				yieldsDown[6][fill]->Fill(controlR + 1, scal*0.95);
				yieldsUp[6][fill]->Fill(controlR + 1, scal*1.05);
			}
			//id eff unc
			unsigned flavorC = tril_flavorComb(ind, _flavors, lCount);
			if(flavorC == 0 || flavorC == 3){
				yieldsDown[7][fill]->Fill(controlR + 1, scal*0.94);
				yieldsUp[7][fill]->Fill(controlR + 1, scal*1.06);
			} else{
				yieldsDown[7][fill]->Fill(controlR + 1, scal*0.9553);
				yieldsUp[7][fill]->Fill(controlR + 1, scal*1.0447);	
			}
			//all pdf variations
			for(unsigned pdf = 0; pdf < 100; ++pdf){
				yieldsPdfVar[pdf][fill]->Fill(controlR + 1, scal*_scaleWeight[pdf + 9]);
			}			
    	}
		//Set negative bins to 0 before adding other processes
		for(unsigned bin = 1; bin < yields[effsam]->GetNbinsX() + 1; ++bin){
			if(yields[effsam]->GetBinContent(bin) < 0.) yields[effsam]->SetBinContent(bin, 0.);
			for(unsigned unc = 0; unc < nUnc; ++unc){
				if(yieldsDown[unc][effsam]->GetBinContent(bin) < 0.) yieldsDown[unc][effsam]->SetBinContent(bin, 0.);
				if(yieldsUp[unc][effsam]->GetBinContent(bin) < 0.) yieldsUp[unc][effsam]->SetBinContent(bin, 0.);
			}
			for(unsigned pdf = 0; pdf < 100; ++pdf){
				if(yieldsPdfVar[pdf][effsam]->GetBinContent(bin) < 0.) yieldsPdfVar[pdf][effsam]->SetBinContent(bin, 0.);
			}
		}
	}
	//Calculate rms of pdf shape variations
	//for total yields
	for(unsigned effsam = 1; effsam < nSamples_eff; ++effsam){
		for(unsigned b = 1; b < yields[effsam]->GetNbinsX() + 1; ++b){
			double pdfVarRms = 0;
			for(unsigned pdf = 0; pdf < 100; ++pdf){
				pdfVarRms += (yields[effsam]->GetBinContent(b) - yieldsPdfVar[pdf][effsam]->GetBinContent(b))*(yields[effsam]->GetBinContent(b) - yieldsPdfVar[pdf][effsam]->GetBinContent(b));
			}
			pdfVarRms = 0.01*sqrt(pdfVarRms);
			yieldsDown[2][effsam]->SetBinContent(b, yields[effsam]->GetBinContent(b) - pdfVarRms);
			yieldsUp[2][effsam]->SetBinContent(b, yields[effsam]->GetBinContent(b) + pdfVarRms);
		}
	}
	
	//Split data and MC histograms for plotting and propagating uncertainties
	//Total data and background yields in each channel
	TH1D* dataYields;
	dataYields = (TH1D*) yields[0]->Clone();
	
	TH1D* bkgYields[nSamples_eff];	
	for(unsigned effsam = 1; effsam < nSamples_eff + 1; ++effsam){
		bkgYields[effsam - 1] = (TH1D*) yields[effsam]->Clone();
	}
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//Calculate histogram with systematic uncertainty for backgrounds and print the ranges of their size
	//const double extraUnc[nSamples_eff] = {1.0, 1.5, 1.0, 1.0, 1.5, 1.3}; //extra flat uncertainties assigned to each background
	//const double extraUnc[nSamples_eff] = {1.25, 1.5, 1.085, 1.15, 1.15, 1.3};
	const double extraUnc[nSamples_eff] = {1.103, 1.5, 1.094, 1.087, 1.5, 1.3};
	double flatSyst[2] = {0.025};
	TH1D* bkgSystYield[nSamples_eff];
	for(unsigned bkg = 0; bkg < nSamples_eff; ++bkg){
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		cout << "background name = " << eff_names[bkg + 1] << endl;	
		flatSyst[1] = fabs(extraUnc[bkg] -1);
		bkgSystYield[bkg] = (TH1D*) bkgYields[bkg]->Clone();
		for(unsigned b = 1;  b < dataYields->GetNbinsX() + 1; ++b){
			bkgSystYield[bkg]->SetBinContent(b,0);
			/*
			if(yields[bkg + 1]->GetBinContent(b) == 0){
				bkgSystYield[bkg]->SetBinContent(b,0);
			} else{
			*/
				if(bkg != nSamples_eff -1){
					double systbin = 0;
					//loop over shape uncertainties
					for(unsigned unc = 0; unc < nUnc; ++unc){
						if(unc == 8) continue;
						double syst = std::max(fabs(yieldsUp[unc][bkg + 1]->GetBinContent(b) - yields[bkg + 1]->GetBinContent(b)), fabs(yieldsDown[unc][bkg + 1]->GetBinContent(b) - yields[bkg + 1]->GetBinContent(b)));
						systbin += syst*syst;
						cout << "size of systematic : " << uncNames[unc] << "  =  " << syst/(yields[bkg + 1]->GetBinContent(b)) << endl;
					}
					//loop over flat uncertainties
					for(unsigned unc = 0; unc < 2; ++unc){
						systbin += bkgYields[bkg]->GetBinContent(b)*bkgYields[bkg]->GetBinContent(b)*flatSyst[unc]*flatSyst[unc];
					}
					bkgSystYield[bkg]->SetBinContent(b, sqrt(systbin));
				} else{
					double syst = std::max(fabs(yieldsUp[8][bkg + 1]->GetBinContent(b) - yields[bkg + 1]->GetBinContent(b)), fabs(yieldsDown[8][bkg + 1]->GetBinContent(b) - yields[bkg + 1]->GetBinContent(b)));
					bkgSystYield[bkg]->SetBinContent(b, sqrt(yields[bkg + 1]->GetBinContent(b)*0.3*yields[bkg + 1]->GetBinContent(b)*0.3 + syst*syst) );
				}
				bkgSystYield[bkg]->SetBinError(b, 0);
			//}
		}
	}

	plotDataVSMC(dataYields, bkgYields, eff_names, nSamples_eff, "crYields_preFit" + extra, true, 0, "HNL", bkgSystYield);
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//PRINT DATA CARDS FOR LIMIT CALCULATION
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	const unsigned nBkg = nSamples_eff;
	const unsigned nSyst = 10 + 2*nSamples_eff; //9 general uncertainties + stat signal + stat bkg + extra unc per bkg
	//TString systNames[nSyst] = {"lumi", "id_eff", "trigeff", "JEC", "metUncl", "pdf", "scale", "pu", "btagSF"};
	TString systNames[nSyst] = {"lumi", "JEC", "metUncl", "pdf", "scale", "pu", "btagSF", "trigeff", "ideff", "fakeEWK"};
	//const TString sigN[3] = {"WZ", "ZZ", "Xgamma"};
	//Reorder background yields
	/*
	TH1D* yieldsC[nSamples_eff + 1];
	TH1D* yieldsDownC[nUnc][nSamples_eff + 1];
	TH1D* yieldsUpC[nUnc][nSamples_eff + 1];
	unsigned samCounter = 0;
	for(unsigned effsam == 1; effsam < nSamples_eff; ++effsam){
		if(eff_names[effsam] != "ZZ/H" && eff_names[effsam] != "WZ" && eff_names[effsam] != "X + #gamma") continue;
		yieldsC[samCounter] = (TH1D*) yields[effsam]->Clone();
	*/
		


	const TString bkgNames[nBkg] = {"ZZH", "triboson", "WZ", "Xgamma", "TTX",  "nonPrompt"}; //slightly rewrite bkg names not to confuse the combine tool
	std::vector<std::vector<double>> systUnc(nSyst, std::vector<double>(nBkg + 1, 0)); //2D array containing all systematics
	//initialize flat and shape systematics
	for(unsigned proc = 0; proc < nBkg - 1; ++proc){ //last background is non-prompt and isn't susceptible to these uncertainty sources (note that entry 0 is signal so that the loop goes up to the second to last bkg)
		systUnc[0][proc] = 1.025;	//lumi
		systUnc[1][proc] = 1.00;	//id eff
		systUnc[2][proc] = 1.00;	//JEC shape
		systUnc[3][proc] = 1.00;	//PU shape
		systUnc[4][proc] = 1.00;	//MET unclustered shape
		systUnc[5][proc] = 1.00; 	//b-tag SF shape
		systUnc[6][proc] = 1.00;	//scale shape
		systUnc[7][proc] = 1.00;	//pdf shape
		systUnc[8][proc] = 1.00;	//trigger shape	
		systUnc[9][proc] = 0;		//fake EWK
	}
	systUnc[9][nBkg] = 1.00;
	TString systDist[nSyst]; //probability distribution of nuisances
	for(unsigned syst = 0; syst < nSyst; ++syst){
		if(syst > 1 && syst < 10) systDist[syst] = "shape";
		else systDist[syst] = "lnN";
	}
	for(unsigned syst = 10 + nBkg; syst < nSyst; ++syst){//loop over last nBkg uncertainties, being the exta uncertainties for each bkg
		unsigned bkg = syst - 10 - nBkg;//index of the background corresponding to the uncertainty index
		systNames[syst] = "extra" + bkgNames[bkg];
		systUnc[syst][bkg] = extraUnc[bkg];
	}
	const unsigned binC = 3;
	TH1D* histCent[nSamples_eff + 1][binC];// = new TH1D("yieldCent", "yieldCent; ; yield", 1, 0, 1);
	TH1D* histUp[nUnc][nSamples_eff + 1][binC];// = new TH1D("yieldUp", "yieldUp;  ; yield", 1, 0, 1);
	TH1D* histDown[nUnc][nSamples_eff + 1][binC];// = new TH1D("yieldDown", "yieldDown;  ; yield", 1, 0, 1);

	for(unsigned effsam = 0; effsam < nSamples_eff + 1; ++effsam){
		for(unsigned bin = 0; bin < binC; ++bin){
			histCent[effsam][bin] = new TH1D("yieldCent" + eff_names[effsam] + std::to_string(bin), "yieldCent"  + eff_names[effsam] + std::to_string(bin) + ";  ; yield", 1, 0, 1);
			for(unsigned unc = 0; unc < nUnc; ++unc){
				histDown[unc][effsam][bin]  = new  TH1D("yieldsDown"  +  eff_names[effsam]  + uncNames[unc] + std::to_string(bin), "yieldsDown"  +  eff_names[effsam]  + uncNames[unc] + std::to_string(bin) + "; ; yield", 1, 0, 1);
				histUp[unc][effsam][bin] = new  TH1D("yieldsUp" +  eff_names[effsam]  + uncNames[unc] + std::to_string(bin), "yieldsUp"  +  eff_names[effsam]  + uncNames[unc] + std::to_string(bin) + "; ; yield", 1, 0, 1);
			}
		}
	}
	//loop over all signals, categories and bins to fill final uncertainties and print datacards
	//loop over all categories and bins
	for(unsigned bin = 1;  bin < dataYields->GetNbinsX() + 1; ++bin){
		//change shape uncertainty names for uncorrelated uncertainty sources
		systNames[2] += "bin" + std::to_string(bin);
		systNames[4] += "bin" + std::to_string(bin);
		//make root file 
		TFile* shapeFile = new TFile("./datacards/shapes/shapeFile_bkgNorm_bin" + (TString) std::to_string(bin) + extra + ".root", "recreate");
		histCent[0][bin]->SetBinContent(1, dataYields->GetBinContent(bin));
		histCent[0][bin]->Write("data_obs");//"data_obs_ch"  + (TString) std::to_string(bin)
		for(unsigned bkg = 0; bkg < nBkg; ++bkg){
			histCent[bkg + 1][bin - 1]->SetBinContent(1, std::max(yields[bkg + 1]->GetBinContent(bin), 0.));
			histCent[bkg + 1][bin - 1]->Write(bkgNames[bkg]);//bkgNames[bkg] + "_ch" + std::to_string(bin)
			for(unsigned unc = 0; unc < nUnc; ++unc){
				if(unc != 8 && bkg == nBkg -1) continue;
				if(unc == 8 && bkg != nBkg -1) continue;
				histUp[unc][bkg + 1][bin - 1]->SetBinContent(1, std::max(yieldsUp[unc][bkg + 1]->GetBinContent(bin), 0.));
				histDown[unc][bkg + 1][bin - 1]->SetBinContent(1, std::max(yieldsDown[unc][bkg + 1]->GetBinContent(bin), 0.));
				histUp[unc][bkg + 1][bin - 1]->Write(bkgNames[bkg] + "_" + systNames[1 + unc] + "Up");//bkgNames[bkg] + "_ch"  + (TString) std::to_string(bin) + "_" + systNames[3 + unc] + "Up"
				histDown[unc][bkg + 1][bin -1]->Write(bkgNames[bkg] + "_" +  systNames[1 + unc] + "Down");//bkgNames[bkg] + "_ch"  + (TString) std::to_string(bin)+ "_" +  systNames[3 + unc] + "Down"
			}
		}
		shapeFile->Close();
		double bkgYieldVal[nBkg];//bkg rates
		//bkg stat unc
		for(unsigned bkg = 0; bkg < nBkg; ++bkg){
			bkgYieldVal[bkg] = bkgYields[bkg]->GetBinContent(bin);
			systUnc[10 + bkg][bkg]  = 1 + ((bkgYields[bkg]->GetBinContent(bin) == 0) ? 0. : std::max(0., bkgYields[bkg]->GetBinError(bin)/bkgYields[bkg]->GetBinContent(bin)) );
			systNames[10 + bkg] = "stat"  + bkgNames[bkg] + std::to_string(bin);
		}
		//print datacard
		std::ofstream card;
		card.open("datacards/datacard_bkgNorm_bin" + (TString) std::to_string(bin) + extra + ".txt");
		//define number of channels, background sources and systematics
		card << "imax 1 number of channels \n";
		card << "jmax " << nBkg  - 1<< " number of backgrounds \n";
		card << "kmax " << nSyst << " number of nuisance parameters (sources of systematical uncertainties) \n";
		card << "---------------------------------------------------------------------------------------- \n";
		//define the channels and the number of observed events
		card << "bin	" << bin << "\n";
		card << "observation " << dataYields->GetBinContent(bin) << "\n";
		//define all backgrounds and their yields
		card << "---------------------------------------------------------------------------------------- \n";
		card << "shapes * * " << "shapes/shapeFile_bkgNorm_bin" + (TString) std::to_string(bin) + extra + ".root  $PROCESS $PROCESS_$SYSTEMATIC\n";
		card << "---------------------------------------------------------------------------------------- \n";
		card << "bin	";
		for(unsigned proc = 0; proc < nBkg; ++proc){
			card << "	" << bin;
		}
		card << "\n";
		card << "process";
		for(unsigned bkg = 0; bkg < nBkg; ++bkg){
			//if(bkgNames[bkg] == "WZ" || bkgNames[bkg] == "ZZ/H" || bkgNames[bkg] == "Xgamma"){
				card << "	" << bkgNames[bkg];
			//}
		}
		/*
		for(unsigned bkg = 0; bkg < nBkg; ++bkg){
			if(bkgNames[bkg] != "WZ" && bkgNames[bkg] != "ZZ/H" && bkgNames[bkg] != "Xgamma"){
				card << "	" << bkgNames[bkg];
			}
		}
		*/
		card << "\n";
		card << "process";
		for(int bkg = -2; bkg < ((int) nBkg - 2); ++bkg){
			if(bkgNames[bkg + 2] == "ZZH") card << "	" << -2;
			else if(bkgNames[bkg + 2] == "WZ"  ) card << "	" << -1;
			else if(bkgNames[bkg + 2] == "Xgamma") card <<"	" << 0;
			else if(bkgNames[bkg + 2] == "triboson") card << "	" << 1;
			else if(bkgNames[bkg + 2] == "TTX") card << "	" << 2;
			else if(bkgNames[bkg + 2] == "nonPrompt") card << "	" << 3;
			
			//card << "	" << bkg;
		}
		card << "\n";
		card <<	"rate";
		for(unsigned bkg = 0; bkg < nBkg; ++bkg){
			//if(bkgNames[bkg] == "WZ" || bkgNames[bkg] == "ZZH" || bkgNames[bkg] == "Xgamma"){
				if(bkgYieldVal[bkg] <= 0) card << "	" << "0.00";
				else card << "	" << bkgYieldVal[bkg];
			//}
		}
		/*
		for(unsigned bkg = 0; bkg < nBkg; ++bkg){
			if(bkgNames[bkg] != "WZ" && bkgNames[bkg] != "ZZH" && bkgNames[bkg] != "Xgamma"){
				if(bkgYieldVal[bkg] <= 0) card << "	" << "0.00";
				else card << "	" << bkgYieldVal[bkg];
			}
		}
		*/
		card << "\n";
		card << "---------------------------------------------------------------------------------------- \n";
		//define sources of systematic uncertainty, what distibution they follow and how large their effect is
		for(unsigned syst = 0; syst < nSyst; ++syst){
			//if(syst > 8 && (bkgNames[(syst-9)%nBkg] != "WZ" && bkgNames[(syst-9)%nBkg] != "ZZH" && bkgNames[(syst-9)%nBkg] != "Xgamma") ) continue;
			card << systNames[syst] << "	" << systDist[syst];
			for(unsigned bkg = 0; bkg < nBkg; ++bkg){
				card << "	";
				if(systUnc[syst][bkg] == 0) card << "-";
				else card << systUnc[syst][bkg];
			}
			card << "\n";
		}
		/*
		for(unsigned syst = 9; syst < nSyst; ++syst){
			if(bkgNames[(syst-9)%nBkg] == "WZ" || bkgNames[(syst-9)%nBkg] == "ZZH" || bkgNames[(syst-9)%nBkg] == "Xgamma") continue;
			card << systNames[syst] << "	" << systDist[syst];
			for(unsigned bkg = 0; bkg < nBkg; ++bkg){
				card << "	";
				if(systUnc[syst][bkg] == 0) card << "-";
				else card << systUnc[syst][bkg];
			}
			card << "\n";
		}
		*/
		card.close();
		systNames[2] =  "metUncl";
		systNames[4] =  "scale";
		//increment bin Counter
	}
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

}



int main(int argc, char* argv[]){
	TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;
	testtree.Loop();
	rootapp->Run();
    return 0;
}


