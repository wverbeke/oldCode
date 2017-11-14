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
	
	const TString eff_names[nSamples_eff + 1] = {"data", "m_{N} = 1 GeV", "m_{N} = 2 GeV", "m_{N} = 5 GeV", "m_{N} = 10 GeV", "m_{N} = 20 GeV", "m_{N} = 30 GeV", "m_{N} = 40 GeV", "m_{N} = 50 GeV", "m_{N} = 60 GeV", "m_{N} = 80 GeV", "m_{N} = 100 GeV", "m_{N} = 130 GeV", "m_{N} = 150 GeV", "m_{N} = 200 GeV", "m_{N} = 400 GeV", "m_{N} = 600 GeV", "m_{N} = 800 GeV", "m_{N} = 1000 GeV", "ZZ/H", "triboson", "WZ", "X + #gamma", "TT/T + X",  "non-prompt"};
	const unsigned nCat = 6;  //Number of categories
	const TString catNames[nCat] = {"lowM_3lOSSF_lowPt", "lowM_3lnoOSSF_lowPt", "lowM_3lOSSF_highPt", "lowM_3lnoOSSF_highPt", "highM_3lOSSF", "highM_3lnoOSSF"};
	const unsigned nSR[nCat] = {12, 4, 12, 4, 16, 9}; //numbers of search regions
	TH1D* sideBandYields[nCat][nSamples_eff + 1];
	TH1D* sideBandYields5GeVFake[nCat][nSamples_eff + 1];

	for(unsigned cat = 0; cat < nCat; ++cat){
		for(unsigned effsam = 0; effsam < nSamples_eff + 1; ++effsam){
			sideBandYields[cat][effsam] = new TH1D(catNames[cat] + eff_names[effsam], catNames[cat] + eff_names[effsam] + "; search region ; Events/search region", nSR[cat], 0.5, nSR[cat] + 0.5);	
			sideBandYields[cat][effsam]->Sumw2();
			sideBandYields5GeVFake[cat][effsam] = new TH1D(catNames[cat] + eff_names[effsam] + "_5GeVFake", catNames[cat] + eff_names[effsam] + "_5GeVFake; search region ; Events/search region", nSR[cat], 0.5, nSR[cat] + 0.5);	
			sideBandYields5GeVFake[cat][effsam]->Sumw2();
		}
	}

	double xGammaWeights[9][nCat];
	for(unsigned i = 0; i < 9; ++i){
		for(unsigned cat = 0; cat < nCat; ++cat){
			xGammaWeights[i][cat] = 0;
		}
	}
	double DYCounter = 0;
	double TTCounter = 0;

    Double_t scale[nSamples -1];

	//set to check which data events have already been processed
	std::set<std::tuple<unsigned long, unsigned long, unsigned long> > usedEvents; //runNb, lumiBlock, eventNb	
	//Loop over all samples
	for(unsigned sam = 0, effsam = 0; sam < nSamples; ++sam, ++effsam){
		if(sam != 0){
   			if(names[sam] == names[sam -1]) --effsam;
    	}		
		if(effsam > 0 && effsam <= nSig) continue;
		//if(effsam != 0) continue;
    	Long64_t nEntries = inputTree[sam]->GetEntries();
    	if(sam > 0){
   			scale[sam -1] = xSections[sam -1]*DataLuminosity*1000/(hcounter[sam]);
    	}
		cout << eff_names[effsam] << endl;
    	std::cout<<"Entries in "<< fileList[sam] <<" "<<nEntries<<std::endl;
		cout << effsam << endl;
		
        for(Long64_t it = 0; it < nEntries/5; ++it){
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
			if(nBJets(true, false, 0) != 0) continue;

			//Check if data events were used before
			if(effsam == 0){
				auto event = usedEvents.find(std::make_tuple(_eventNb, _lumiBlock, _runNb));
				if(event != usedEvents.end()) continue;
				usedEvents.insert(std::make_tuple(_eventNb, _lumiBlock, _runNb));
			}
			//Categorize according to the number of leptons and flavors
			unsigned* ind = new unsigned[_nL];
			unsigned lCount = lepOrder(ind, 3, true, true);
			if(lCount != 3) continue; //Veto 4th FO lepton considering signal model!
			//Apply analysis Pt thresholds
			if(!ptCuts_hnl(ind,lCount)) continue;
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
			//Require 3 leptons to be tight in data and MC, and determine nonPrompt bkg in data
			unsigned nTight = tightCount(ind, lCount);
			bool tightFail = nTight < 3;
			//index used to fill events, needed to separate fakes from data
			unsigned fill = effsam;
			//Apply FR maps to data control region
			double* conePt = new double[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				conePt[l] = _lPt[ind[l]]*(1 + std::max(_isolation[ind[l]] - 0.1, 0.));
			}
			if(!tightFail){
				if(fileList[sam] == "TTJets_DiLept.root" || fileList[sam] == "TTJets_SingleLeptFromTbar.root"  || fileList[sam] == "TTJets_SingleLeptFromT.root"){
					TTCounter += scal;
				} else if(fileList[sam] == "DYJetsToLL_M10to50.root" || fileList[sam] == "DYJetsToLL_M50.root"){
					DYCounter += scal;
				}
			}


			if(!tightFail) continue;
			fill = 0;
			if(effsam> nSig) fill = effsam;
			if(effsam != 0 && effsam <= nSig) continue;
			
			if(tightFail && (effsam == 0 || effsam > nSig)){ //&& effsam == 0){
				//fakes go in different histogram
				//fill = nSamples_eff;
				//Apply FR maps
				//if(effsam != 0) scal *= -1;
				//scal*= fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap, lCount);
				scal *= -1.;
				for(unsigned l = 0; l < lCount; ++l){
					if(!_istight[ind[l]]){
						//double fr = frMap[flavors[ind[l]]]->GetBinContent(frMap[flavors[ind[l]]]->FindBin(TMath::Min(conePt[l],49.), fabs(eta[ind[l]])));
						//double fr = frMap[flavors[ind[l]]]->GetBinContent(frMap[flavors[ind[l]]]->FindBin(TMath::Min(conePt[l], 65.), TMath::Min(fabs(eta[ind[l]]), 2.4) ) );
						//weight *= -fr/(1-fr);
						scal *= -1.;
					}
				}
			} else if(tightFail) continue;
			/*
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
			*/
			/*
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
			*/
			//Apply triggers to data events;
			if(sam == 0 || effsam > nSig){
				bool trigPass[4];
				trigPass[0] = _lowM_trigger_eee || _lowM_trigger_all;
				trigPass[1] = _lowM_trigger_mee || _lowM_trigger_all;
				trigPass[2] = _lowM_trigger_mme || _lowM_trigger_all;
				trigPass[3] = _lowM_trigger_mmm || _lowM_trigger_all;
				if(!trigPass[tril_flavorComb(ind, _flavors, lCount)]) continue;	
			}				
			//determine search category
			TLorentzVector* lepV = new TLorentzVector[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				lepV[l].SetPtEtaPhiE(conePt[l], _lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]*(conePt[l]/_lPt[ind[l]]) );
			}
			//Calculate lepton system vector
			TLorentzVector lepSyst;
			for(int l = 0; l < 3; ++l) lepSyst += lepV[l];


			unsigned cat = hnl::cat(ind, _flavors, _charges, lCount, conePt[0]);
			if(cat == 999){
				continue;	
			}
			
			//Apply ID and reco SF to simulation
			if(effsam != 0){
				scal*=getEventSF(ind, lCount, true);
			}
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
				if((cat == 0 || cat == 2) && !vetoLowMll(5)) continue; //NEW
				if(lepSyst.M() > 80) continue;
				if(_met > 75) continue;
			}
			else if(cat > 3){
				if(conePt[1] < 15) continue;
				if(conePt[2] < 10) continue;
				//if(conePt[2] < 10 - 5*(_flavors[ind[2]])) continue;
				if(cat == 4){
					if(fabs(mll - 91) < 15) continue; //Veto onZ events (WZ and DY)
					if(fabs(lepSyst.M() - 91) < 15) continue; //veto conversions
				}
				if(!vetoLowMll(5)) continue;
			}
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
			//TEMPORARY ZGAMMA REMOVAL
			//if(fileList[sam] == "ZGTo2LG.root") continue;
			//fill Nominal yields
			double mt_min = transmass(lepV[lw_min], METvec);
			unsigned searchR = hnl::sr(mt_min, minMos, lepSyst.M(), cat);
			sideBandYields[cat][fill]->Fill(searchR + 1, scal);

			if(sam >= (nSig + 18) && sam <= (nSig + 26)){
				cout << "scal = " << scal << endl;
				xGammaWeights[sam - nSig - 18][cat] = xGammaWeights[sam - nSig - 18][cat] +  scal;
			}

			//Fill 5GeV fake sideband
			bool lowFake = false;
			for(unsigned l = 0; l < lCount; ++l){
				if(_isFO[ind[l]] && !_istight[ind[l]]){
					if(_lPt[ind[l]]*(1 + std::max(_isolation[ind[l]] - 0.1, 0.)) < 10.){
						lowFake = true;
						break;
					}
				}
			}
			if(lowFake){
				sideBandYields5GeVFake[cat][fill]->Fill(searchR + 1, scal);
			}
    	}
		/*
		//Set negative bins to 0 before adding other processes
		for(unsigned cat = 0; cat < nCat; ++cat){
			for(unsigned bin = 1; bin < sideBandYields[cat][effsam]->GetNbinsX() + 1; ++bin){
				if(sideBandYields[cat][effsam]->GetBinContent(bin) < 0.) sideBandYields[cat][effsam]->SetBinContent(bin, 0.);
			}
		} 
		*/
	}
	/*
	//Set fakes to 0 if they are negative
	for(unsigned cat = 0; cat < nCat; ++cat){
		for(unsigned bin = 1; bin < yields[cat][nSamples_eff]->GetNbinsX() + 1; ++bin){
			if(yields[cat][nSamples_eff]->GetBinContent(bin) < 0) yields[cat][nSamples_eff]->SetBinContent(bin, 0.);
		}
	}
	*/
	for(unsigned cat = 0; cat < nCat; ++cat){
		if(cat == 0 || cat == 2) continue;
		plotHist(sideBandYields[cat][0], "sideBandYield_" + catNames[cat], false, true);
	}
	
	for(unsigned cat = 0; cat < nCat; ++cat){
		if(cat == 0 || cat == 2) continue;
		plotHist(sideBandYields5GeVFake[cat][0], "sideBandYield_5GeVFakes_" + catNames[cat], false, true);
	}

	TH1D* dataYields[nCat];
	for(unsigned cat = 0; cat < nCat; ++cat){
		dataYields[cat] = (TH1D*) sideBandYields[cat][nSig + 1]->Clone();
	}
	TH1D* bkgYields[nCat][nSamples_eff -nSig]; 
	TH1D* bkgErros[nCat][nSamples_eff - nSig];
	for(unsigned cat = 0; cat < nCat; ++cat){
		for(unsigned effsam = nSig + 1; effsam < nSamples_eff; ++effsam){
			bkgYields[cat][effsam -nSig - 1] = (TH1D*) sideBandYields[cat][effsam]->Clone();
			if(effsam > nSig + 1){
				dataYields[cat]->Add(bkgYields[cat][effsam - nSig - 1]);
			}
		}
	}

	TH1D* bkgSyst[nCat][nSamples_eff];
	for(unsigned cat = 0; cat < nCat; ++cat){
		for(unsigned bkg = 0; bkg < nSamples_eff - nSig - 1; ++bkg){
			bkgSyst[cat][bkg] = (TH1D*) bkgYields[cat][bkg]->Clone();
			for(unsigned b = 1; b < bkgSyst[cat][bkg]->GetNbinsX() + 1; ++b){
				bkgSyst[cat][bkg]->SetBinContent(b, 0);
			}
		}
	}

	const TString distNames[nSamples_eff + 1 - nSig] = {"total SB MC", "ZZ/H", "triboson", "WZ", "X + #gamma", "TT/T + X"};
	//Plot the yields as a function of the search region
	for(unsigned cat = 0; cat < nCat; ++cat){
		if(cat == 0 || cat == 2) continue;
		plotDataVSMC(dataYields[cat], bkgYields[cat], distNames, nSamples_eff - nSig - 1, "sideBandYieldMCPrompt_" + catNames[cat] + extra, true, 0, "HNL", bkgSyst[cat]);
	}
	
	std::cout << "Check of individual sample contribution for Xgamma" << std::endl;
	for(unsigned cat = 0; cat < nCat; ++ cat){
		if(cat == 0 || cat == 2) continue;
		std::cout << "####################################" << std::endl;
		std::cout << "category : " << catNames[cat] << std::endl;
		for(unsigned i = 0; i < 9; ++i){
			std::cout << "sample " << fileList[nSig + 18 + i] << " has SumOfWeights = " << xGammaWeights[i][cat] << std::endl;
			std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		}
	}
	
}



int main(int argc, char* argv[]){
	TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;
	testtree.Loop();
	rootapp->Run();
    return 0;
}


