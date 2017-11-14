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
#include <limits>

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

	const unsigned nSamples = 50;
	const unsigned nSamples_eff = 24;
	const unsigned nSig = 18;
	const double lowMCoupling = 0.00001;
	const double highMCoupling = 0.01;
	const double couplingCorrection = 10000;
	//~~~~~~~~ background normalizations~~~~~~~~~~~~~~~~~~~~
	const double glugluToZZkFactor = 2.1;
	const double WZSF = 0.655;
	const double ZZSF = 1.032;
	const double XgammaSF = 0.945;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	const TString fileList[nSamples] = {"data_combined_trilepton.root", 
										"HeavyNeutrino_M1_2l.root", "HeavyNeutrino_M2_2l.root", "HeavyNeutrino_M5_2l.root", "HeavyNeutrino_M10_2l.root", "HeavyNeutrino_M20_2l.root", "HeavyNeutrino_M30_2l.root", "HeavyNeutrino_M40_2l.root", "HeavyNeutrino_M50_2l.root", "HeavyNeutrino_M60_2l.root", "HeavyNeutrino_M80_2l.root", "HeavyNeutrino_M100_2l.root", "HeavyNeutrino_M130_2l.root", "HeavyNeutrino_M150_2l.root", "HeavyNeutrino_M200_2l.root", "HeavyNeutrino_M400_2l.root", "HeavyNeutrino_M600_2l.root", "HeavyNeutrino_M800_2l.root", "HeavyNeutrino_M1000_2l.root",
 										"ZZTo4L.root",  "GluGluToZZTo4mu.root", "GluGluToZZTo4e.root", "GluGluToZZTo4tau.root", "GluGluToZZTo2e2mu.root", "GluGluToZZTo2e2tau.root", "GluGluToZZTo2mu2tau.root", "VHToNonbb.root", "GluGluHToZZTo4L_M125.root", "VBF_HToZZTo4L_M125.root", "WWG.root","WWW.root", "WWZ.root", "WWTo2L2Nu_DoubleScattering.root", "WWTo2L2Nu.root",  "ZZZ.root", "WZTo3LNu_mllmin01.root", "TTGJets.root","ZGTo2LG.root", "WGToLNuG.root", "TGJets.root", "TTJets_DiLept.root", "TTJets_SingleLeptFromTbar.root", "TTJets_SingleLeptFromT.root", "DYJetsToLL_M10to50.root", "DYJetsToLL_M50.root", "ttHToNonbb.root", "TTWJetsToLNu.root", "TTZToLLNuNu.root", "TTZToLL_M1to10.root",  "TTTT.root"};
	
	const double xSections[nSamples - 1] = { 5.237e-01*couplingCorrection*lowMCoupling, 5.238e-01*couplingCorrection*lowMCoupling, 5.215e-01*couplingCorrection*lowMCoupling, 5.132e-01*couplingCorrection*lowMCoupling, 4.764e-01*couplingCorrection*lowMCoupling, 4.160e-01*couplingCorrection*lowMCoupling, 3.335e-01*couplingCorrection*lowMCoupling, 3.335e-01*couplingCorrection*lowMCoupling, 1.336e-01*couplingCorrection*lowMCoupling, 2.473e-02*couplingCorrection*lowMCoupling*10, 1.033e-03*couplingCorrection*highMCoupling, 2.874e-04*couplingCorrection*highMCoupling, 1.561e-04*couplingCorrection*highMCoupling, 4.929e-05*couplingCorrection*highMCoupling, 3.505e-06*couplingCorrection*highMCoupling, 6.922e-07*couplingCorrection*highMCoupling, 2.013e-07*couplingCorrection*highMCoupling, 7.239e-08*couplingCorrection*highMCoupling, 
											1.256*ZZSF, 0.00159*glugluToZZkFactor*ZZSF, 0.00159*glugluToZZkFactor*ZZSF, 0.00159*glugluToZZkFactor*ZZSF, 0.00319*glugluToZZkFactor*ZZSF, 0.00319*glugluToZZkFactor*ZZSF, 0.00319*glugluToZZkFactor*ZZSF,  0.9561, 0.01212, 0.001034,  0.2147, 0.2086, 0.1651, 0.1729,   12.178, 0.01398, 58.59*WZSF, 3.697, 123.9*XgammaSF, 489*XgammaSF,  2.967, 87.315, 182.175, 182.175, 18610, 1921.8*3,  0.215, 0.2043, 0.2529, 0.0493, 0.009103}; 
	const TString names[nSamples] = {"data", "m_{N} = 1 GeV", "m_{N} = 2 GeV", "m_{N} = 5 GeV", "m_{N} = 10 GeV", "m_{N} = 20 GeV", "m_{N} = 30 GeV", "m_{N} = 40 GeV", "m_{N} = 50 GeV", "m_{N} = 60 GeV", "m_{N} = 80 GeV", "m_{N} = 100 GeV", "m_{N} = 130 GeV", "m_{N} = 150 GeV", "m_{N} = 200 GeV", "m_{N} = 400 GeV", "m_{N} = 600 GeV", "m_{N} = 800 GeV", "m_{N} = 1000 GeV",  									 
									 "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "triboson", "triboson", "triboson", "triboson", "triboson", "triboson", "WZ", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X"};

	
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
	readSF(true);
	//Tweakable options////////////////////////////////////////////////////
	const bool TestRun = false;	//Break after a few events
	const double DataLuminosity = 35.867;    //units of fb^{-1}
	const TString extra = "";	//for plot file names
	//////////////////////////
	
	const TString eff_names[nSamples_eff + 1] = {"data", "m_{N} = 1 GeV", "m_{N} = 2 GeV", "m_{N} = 5 GeV", "m_{N} = 10 GeV", "m_{N} = 20 GeV", "m_{N} = 30 GeV", "m_{N} = 40 GeV", "m_{N} = 50 GeV", "m_{N} = 60 GeV", "m_{N} = 80 GeV", "m_{N} = 100 GeV", "m_{N} = 130 GeV", "m_{N} = 150 GeV", "m_{N} = 200 GeV", "m_{N} = 400 GeV", "m_{N} = 600 GeV", "m_{N} = 800 GeV", "m_{N} = 1000 GeV",
												 "ZZ/H", "triboson", "WZ", "X + #gamma", "TT/T + X",  "non-prompt"};
	const unsigned nCat = 3;  //Number of categories
	const TString catNames[nCat] = {"lowM", "highMOSSF", "highMnoOSSF"};
	const unsigned nCuts[nCat] = {9, 10, 7}; //numbers of search regions
	const unsigned nUnc = 6;//number of shape uncertainties
	const TString uncNames[nUnc] = {"jec", "metUnclustered",  "scale", "pdf", "pu", "btagSF"}; //"scaleAcc"};
	TH1D* yields[nCat][nSamples_eff + 1]; //nominal yields in every SR
	/*
	TH1D* yieldsDown[nUnc][nCat][nSamples_eff + 1]; //yields varied down by shape unc
	TH1D* yieldsUp[nUnc][nCat][nSamples_eff + 1]; //yields varied up by shape unc
	TH1D* yieldsPdfVar[100][nCat][nSamples_eff + 1]; //yields for all 100 possible pdf variations
	*/
	for(unsigned cat = 0; cat < nCat; ++cat){
		for(unsigned effsam = 0; effsam < nSamples_eff + 1; ++effsam){
			yields[cat][effsam] = new TH1D(catNames[cat] + eff_names[effsam], catNames[cat] + eff_names[effsam] + "; search region ; Events/search region", nCuts[cat], 0.5, nCuts[cat] + 0.5);	
			yields[cat][effsam]->Sumw2();
			/*
			for(unsigned unc = 0; unc < nUnc; ++unc){
				yieldsDown[unc][cat][effsam] = new TH1D(catNames[cat]  +  eff_names[effsam] + uncNames[unc] + "Up", catNames[cat] + eff_names[effsam] + uncNames[unc] + "Up; search region ; Events/search region", nSR[cat], 0.5, nSR[cat] + 0.5);
				yieldsDown[unc][cat][effsam]->Sumw2();
				yieldsUp[unc][cat][effsam] = new  TH1D(catNames[cat]  +  eff_names[effsam]  + uncNames[unc] + "Down", catNames[cat] + eff_names[effsam]  + uncNames[unc] + "Down; search region ; Events/search region", nSR[cat], 0.5, nSR[cat] + 0.5);
				yieldsUp[unc][cat][effsam]->Sumw2();
			}
			for(unsigned pdf = 0; pdf < 100; ++pdf){
				yieldsPdfVar[pdf][cat][effsam] = new  TH1D(catNames[cat]  +  eff_names[effsam]  + "_pdf" + std::to_string(pdf), catNames[cat] + eff_names[effsam]  + "_pdf" + std::to_string(pdf) + "; search region ; Events/search region", nSR[cat], 0.5, nSR[cat] + 0.5);
				yieldsPdfVar[pdf][cat][effsam]->Sumw2();
			}
			*/
		}
	}
    Double_t scale[nSamples -1];
	//set to check which data events have already been processed
	std::set<std::tuple<unsigned long, unsigned long, unsigned long> > usedEvents; //runNb, lumiBlock, eventNb	
	//Loop over all samples
	for(unsigned sam = 0, effsam = 0; sam < nSamples; ++sam, ++effsam){
		if(sam != 0){
   			if(names[sam] == names[sam -1]) --effsam;
    	}		
    	Long64_t nEntries = inputTree[sam]->GetEntries();
    	if(sam > 0){
   			scale[sam -1] = xSections[sam -1]*DataLuminosity*1000/(hcounter[sam]);
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
				scal = scale[sam-1]*_weight;
			}		
			cutBased();
	
			//Baseline event selection
			if(!baseline(true, true, false, false)) continue;
			//Check if data events were used before
			if(effsam == 0){
				auto event = usedEvents.find(std::make_tuple(_eventNb, _lumiBlock, _runNb));
				if(event != usedEvents.end()) continue;
				usedEvents.insert(std::make_tuple(_eventNb, _lumiBlock, _runNb));
			}
			//Categorize according to the number of leptons and flavors
			unsigned* ind = new unsigned[_nL];
			unsigned lCount = lepOrder(ind, 3, true, true);
			if(lCount < 3) continue;
			unsigned flavorComp = hnl::flavorComposition(ind, _flavors, _charges, lCount);
			if(flavorComp > 1) continue;
			double* conePt = new double[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				conePt[l] = _lPt[ind[l]]*(1 + std::max(_isolation[ind[l]] - 0.1, 0.));
			}

			//MC prompt matching
			if(effsam  > nSig){
				bool promptfail = false;
				for(unsigned l = 0; l < lCount; ++l){	
					if(_origin[ind[l]] != 0){
						promptfail = true;
						break;
					}
				}
				if(promptfail) continue;
			}
			//Sample overlap removal
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
						if(_pdgmc[ind[l]] == 22 && _originPhot[ind[l]] == 0){
							promptfail = false;
							break;
						}
					}
					if(promptfail) continue;
				}
			}
			//if(lCount != 3) continue; //Veto 4th FO lepton considering signal model!
			//unsigned cat = hnl::cat(ind, _flavors, _charges, lCount, conePt[0]);
			/*
			if(cat == 999){
				continue;	
			}
			*/

			//Fill first entry for signal
			if(effsam <= nSig && effsam != 0){
				if(flavorComp == 0){
					yields[1][effsam]->Fill(1, scal);
				} else if(flavorComp == 1){
					yields[0][effsam]->Fill(1, scal); 
					yields[2][effsam]->Fill(1, scal);
				}
			}
			//pt cuts
			if(!ptCuts_hnl(ind,lCount)) continue;
			//triggers
			bool trigPass[4];
			trigPass[0] = _lowM_trigger_eee || _lowM_trigger_all;
			trigPass[1] = _lowM_trigger_mee || _lowM_trigger_all;
			trigPass[2] = _lowM_trigger_mme || _lowM_trigger_all;
			trigPass[3] = _lowM_trigger_mmm || _lowM_trigger_all;
			if(!trigPass[tril_flavorComb(ind, _flavors, lCount)]) continue;	
			//Require 3 leptons to be tight in data and MC, and determine nonPrompt bkg in data
			unsigned nTight = tightCount(ind, lCount);
			bool tightFail = nTight < 3;
			//index used to fill events, needed to separate fakes from data
			unsigned fill = effsam;
			//Apply FR maps to data control region
			if(effsam != 0){
				scal*=getEventSF(ind, lCount, true);
			}

			if(tightFail && (effsam == 0 || effsam > nSig)){ //&& effsam == 0){
				//fakes go in different histogram
				fill = nSamples_eff;
				//Apply FR maps
				if(effsam != 0) scal *= -1;
				scal*= fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap, lCount);
			} else if(tightFail) continue;

			if(flavorComp == 0){
				yields[1][fill]->Fill(2, scal);
			} else if(flavorComp == 1){
				yields[0][fill]->Fill(2, scal); 
				yields[2][fill]->Fill(2, scal);
			}
		
			if(nBJets(true, false, 0) != 0) continue;
			
			if(flavorComp == 0){
				yields[1][fill]->Fill(3, scal);
			} else if(flavorComp == 1){
				yields[0][fill]->Fill(3, scal); 
				yields[2][fill]->Fill(3, scal);
			}
			
			if(lCount != 3) continue;
			
			if(flavorComp == 0){
				yields[1][fill]->Fill(4, scal);
			} else if(flavorComp == 1){
				yields[0][fill]->Fill(4, scal); 
				yields[2][fill]->Fill(4, scal);
			}


			
			
			//determine search category
			TLorentzVector* lepV = new TLorentzVector[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				lepV[l].SetPtEtaPhiE(conePt[l], _lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]*(conePt[l]/_lPt[ind[l]]) );
			}
			//Calculate lepton system vector
			TLorentzVector lepSyst;
			for(int l = 0; l < 3; ++l) lepSyst += lepV[l];
			//Apply ID and reco SF to simulation
			
			//determine which leptons will be used for the calculation of mll
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

			if(conePt[0] > 55){
				if(flavorComp == 0) yields[1][fill]->Fill(5, scal);
				else yields[2][fill]->Fill(5, scal);
				if(conePt[1] < 15) continue;
				if(flavorComp == 0) yields[1][fill]->Fill(6, scal);
				else yields[2][fill]->Fill(6, scal);
				if(conePt[2] < 10) continue;
				if(flavorComp == 0) yields[1][fill]->Fill(7, scal);
				else yields[2][fill]->Fill(7, scal);
				if(flavorComp == 0){
					if(fabs(mll - 91) < 15) continue;
					yields[1][fill]->Fill(8, scal);
					if(fabs(lepSyst.M() - 91) < 15) continue;
					yields[1][fill]->Fill(9, scal);
					if(!vetoLowMll(5)) continue;
					yields[1][fill]->Fill(10, scal);
				}
			} else{
				if(flavorComp == 0) continue;
				yields[0][fill]->Fill(5, scal);
				if(lepSyst.M() > 80) continue;
				yields[0][fill]->Fill(6, scal);
				if(_met > 75) continue;
				yields[0][fill]->Fill(7, scal);
				if(conePt[0] < 30){
					yields[0][fill]->Fill(8, scal);
				} else{
					yields[0][fill]->Fill(9, scal);
				}
			}
			

    	}
		//Set negative bins to 0 before adding other processes
		for(unsigned cat = 0; cat < nCat; ++cat){
			for(unsigned bin = 1; bin < yields[cat][effsam]->GetNbinsX() + 1; ++bin){
				if(yields[cat][effsam]->GetBinContent(bin) < 0.) yields[cat][effsam]->SetBinContent(bin, 0.);
			}
		} 
	}
	//Make total background histograms
	TH1D* bkgTot[nCat];
	for(unsigned cat = 0; cat < nCat; ++ cat){
		bkgTot[cat] = (TH1D*) yields[cat][nSig + 1]->Clone();
		for(unsigned bkg = 1; bkg < nSamples_eff - nSig; ++bkg){
			bkgTot[cat]->Add(yields[cat][nSig + 1 + bkg]);
		}
	}
	TH1D* sigPlot[nCat][5];
	TString sigPlotNames[nCat][5];
	for(unsigned cat = 0; cat < nCat; ++cat){
		if(cat == 0){	
			sigPlot[cat][0] = (TH1D*) yields[cat][3]->Clone();
			sigPlotNames[cat][0] = "$m_{N} =$ 5 GeV";
			sigPlot[cat][1] = (TH1D*) yields[cat][4]->Clone();
			sigPlotNames[cat][1] = "$m_{N} =$ 10 GeV";
			sigPlot[cat][2] = (TH1D*) yields[cat][6]->Clone();	
			sigPlotNames[cat][2] = "$m_{N} =$ 30 GeV";
			sigPlot[cat][3] = (TH1D*) yields[cat][8]->Clone();
			sigPlotNames[cat][3] = "$m_{N} =$ 50 GeV";
			sigPlot[cat][4] = (TH1D*) yields[cat][9]->Clone();
			sigPlotNames[cat][4] = "$m_{N} =$ 60 GeV";
		} else{
			sigPlot[cat][0] = (TH1D*) yields[cat][11]->Clone();
			sigPlotNames[cat][0] = "$m_{N} =$ 100 GeV";
			sigPlot[cat][1] = (TH1D*) yields[cat][13]->Clone();
			sigPlotNames[cat][1] = "$m_{N} =$ 150 GeV";
			sigPlot[cat][2] = (TH1D*) yields[cat][14]->Clone();	
			sigPlotNames[cat][2] = "$m_{N} =$ 200 GeV";
			sigPlot[cat][3] = (TH1D*) yields[cat][15]->Clone();
			sigPlotNames[cat][3] = "$m_{N} =$ 400 GeV";
			sigPlot[cat][4] = (TH1D*) yields[cat][16]->Clone();
			sigPlotNames[cat][4] = "$m_{N} =$ 600 GeV";
		}
	}
	//Make Cut flow tables
	std::string cutNames[nCat][11] = {{"3 FO", "3 tight + $P_{T}$ cuts", "b-veto",  "4th $\\ell$ veto", "$P_{T}(leading)$ $<$ 55 GeV", "$\\Mtril <$ 80 GeV", "$E_{T}^{miss} <$ 75 GeV", "$P_{T} <$ 30 GeV", "30 GeV $< P_{T} <$ 55 GeV"},
									 {"3 FO", "3 tight +  $P_{T}$ cuts", "b-veto",  "4th $\\ell$ veto", "$P_{T}(leading)$ $>$ 55 GeV", "$P_{T}(subleading)$ $>$ 15 GeV", "$P_{T}(trailing)$ $>$ 10 GeV", "$\\left| M_{\\ell\\ell}  - 91. \\right| > 15$ GeV", "$\\left| M_{3\\ell}  - 91. \\right| >$ 15 GeV", "min($M_{OSSF}$) $>$ 5"},
									 {"3 FO", "3 tight + \\pt cuts", "b-veto",  "4th $\\ell$ veto","$P_{T}(leading)$ $>$ 55 GeV", "$P_{T}(subleading)$ $>$ 15 GeV", "$P_{T}(trailing)$ $>$ 10 GeV"}};

	std::ofstream table[nCat];
	for(unsigned cat = 0; cat < nCat; ++ cat){
		table[cat].open("tables/Cutflow_" + catNames[cat] + ".txt");
		table[cat] << "\\resizebox{1\\textwidth}{!}{ \n";
		table[cat] << "\\begin{tabular}{|c|c|c|c|c|c|c|} \n";
		table[cat] << "\\hline \n";
		table[cat] << " & \\textbf{total bkg.}";
		for(unsigned sig = 0; sig < 5; ++sig){
			table[cat] << " & \\textbf{" << sigPlotNames[cat][sig] << "}";	
		}
		table[cat] << "\\\\ \\hline \n";
		for(unsigned cut = 0; cut < nCuts[cat]; ++ cut){
			table[cat] << "\\textbf{" << cutNames[cat][cut] << "}";
			if(cut != 0) table[cat] << " & " << std::fixed << std::setprecision(2) << bkgTot[cat]->GetBinContent(cut + 1);
			else table[cat] << " & / ";
			for(unsigned sig = 0; sig < 5; ++sig){
				table[cat] << " & " <<  std::fixed << std::setprecision(2) << sigPlot[cat][sig]->GetBinContent(cut + 1);
			}
			table[cat] << "\\\\ \\hline  \n";
		}
		table[cat] << "\\end{tabular}} \n";
		table[cat].close();
	}
}



int main(int argc, char* argv[]){
	TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;
	testtree.Loop();
	rootapp->Run();
    return 0;
}


