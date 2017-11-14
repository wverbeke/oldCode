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
										"HeavyNeutrino_M1_e.root", "HeavyNeutrino_M2_e.root", "HeavyNeutrino_M5_e.root", "HeavyNeutrino_M10_e.root", "HeavyNeutrino_M20_e.root", "HeavyNeutrino_M30_e.root", "HeavyNeutrino_M40_e.root", "HeavyNeutrino_M50_e.root", "HeavyNeutrino_M60_e.root", "HeavyNeutrino_M80_e.root", "HeavyNeutrino_M100_e.root", "HeavyNeutrino_M130_e.root", "HeavyNeutrino_M150_e.root", "HeavyNeutrino_M200_e.root", "HeavyNeutrino_M400_e.root", "HeavyNeutrino_M600_e.root", "HeavyNeutrino_M800_e.root", "HeavyNeutrino_M1000_e.root",
 										"ZZTo4L.root",  "GluGluToZZTo4mu.root", "GluGluToZZTo4e.root", "GluGluToZZTo4tau.root", "GluGluToZZTo2e2mu.root", "GluGluToZZTo2e2tau.root", "GluGluToZZTo2mu2tau.root", 	"VHToNonbb.root", "GluGluHToZZTo4L_M125.root", "VBF_HToZZTo4L_M125.root", "WWG.root","WWW.root", "WWZ.root", "WWTo2L2Nu_DoubleScattering.root", "WWTo2L2Nu.root",  "ZZZ.root", "WZTo3LNu_mllmin01.root", "TTGJets.root","ZGTo2LG.root", "WGToLNuG.root", "TGJets.root", "TTJets_DiLept.root", "TTJets_SingleLeptFromTbar.root", "TTJets_SingleLeptFromT.root", "DYJetsToLL_M10to50.root", "DYJetsToLL_M50.root", "ttHToNonbb.root", "TTWJetsToLNu.root", "TTZToLLNuNu.root", "TTZToLL_M1to10.root",  "TTTT.root"};
	
	const double xSections[nSamples - 1] = { 2.617e-01*couplingCorrection*lowMCoupling, 2.607e-01*couplingCorrection*lowMCoupling, 2.600e-01*couplingCorrection*lowMCoupling, 2.571e-01*couplingCorrection*lowMCoupling, 2.381e-01*couplingCorrection*lowMCoupling, 2.088e-01*couplingCorrection*lowMCoupling, 1.677e-01*couplingCorrection*lowMCoupling, 1.190e-01*couplingCorrection*lowMCoupling, 6.715e-02*couplingCorrection*lowMCoupling, 1.242e-02*couplingCorrection*lowMCoupling*10, 5.139e-04*couplingCorrection*highMCoupling, 1.442e-04*couplingCorrection*highMCoupling, 7.833e-05*couplingCorrection*highMCoupling, 2.473e-05*couplingCorrection*highMCoupling, 1.748e-06*couplingCorrection*highMCoupling, 3.487e-07*couplingCorrection*highMCoupling, 1.027e-07*couplingCorrection*highMCoupling, 3.563e-08*couplingCorrection*highMCoupling, 
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
	const unsigned nCat = 6;  //Number of categories
	const TString catNames[nCat] = {"lowM_3lOSSF_lowPt", "lowM_3lnoOSSF_lowPt", "lowM_3lOSSF_highPt", "lowM_3lnoOSSF_highPt", "highM_3lOSSF", "highM_3lnoOSSF"};
	const unsigned nSR[nCat] = {12, 4, 12, 4, 16, 9}; //numbers of search regions
	const unsigned nUnc = 6;//number of shape uncertainties
	const TString uncNames[nUnc] = {"jec", "metUnclustered",  "scale", "pdf", "pu", "btagSF"}; //"scaleAcc"};
	TH1D* yields[nCat][nSamples_eff + 1]; //nominal yields in every SR
	TH1D* yieldsDown[nUnc][nCat][nSamples_eff + 1]; //yields varied down by shape unc
	TH1D* yieldsUp[nUnc][nCat][nSamples_eff + 1]; //yields varied up by shape unc
	TH1D* yieldsPdfVar[100][nCat][nSamples_eff + 1]; //yields for all 100 possible pdf variations
	for(unsigned cat = 0; cat < nCat; ++cat){
		for(unsigned effsam = 0; effsam < nSamples_eff + 1; ++effsam){
			yields[cat][effsam] = new TH1D(catNames[cat] + eff_names[effsam], catNames[cat] + eff_names[effsam] + "; search region ; Events/search region", nSR[cat], 0.5, nSR[cat] + 0.5);	
			yields[cat][effsam]->Sumw2();
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
			if(tightFail && (effsam == 0 || effsam > nSig)){ //&& effsam == 0){
				//fakes go in different histogram
				fill = nSamples_eff;
				//Apply FR maps
				if(effsam != 0) scal *= -1;
				scal*= fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap, lCount);
			} else if(tightFail) continue;

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
			//Apply triggers to data events;
			if(sam == 0 || effsam > nSig){
				bool trigPass[4];
				trigPass[0] = _lowM_trigger_eee || _lowM_trigger_all;
				trigPass[1] = _lowM_trigger_mee || _lowM_trigger_all;
				trigPass[2] = _lowM_trigger_mme || _lowM_trigger_all;
				trigPass[3] = _lowM_trigger_mmm || _lowM_trigger_all;
				if(!trigPass[tril_flavorComb(ind, _flavors, lCount)]) continue;	
			}				
			//Select only events with sensitivity to muons
			if(tril_flavorComb(ind, _flavors, lCount) != 0 && tril_flavorComb(ind, _flavors, lCount) != 1) continue;
			//determine search category
			unsigned cat = hnl::cat(ind, _flavors, _charges, lCount, conePt[0]);
			if(cat == 999){
				continue;	
			}

			//Set TLorentzVector for each lepton
			TLorentzVector* lepV = new TLorentzVector[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				lepV[l].SetPtEtaPhiE(conePt[l], _lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]*(conePt[l]/_lPt[ind[l]]) );
			}
			//Calculate lepton system vector
			TLorentzVector lepSyst;
			for(int l = 0; l < 3; ++l) lepSyst += lepV[l];

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
			
			//Category spefific selection:
			//low mass general
			if(cat < 4){
				if(cat == 0 || cat == 2){
					//Calculate minDeltaR and maxDeltaR
					double minDeltaR = 9999.;
					double maxDeltaR = 0.;
					for(unsigned l = 0; l < lCount - 1; ++l){
						for(unsigned k = l + 1; k < lCount; ++k){
							if(lepV[l].DeltaR(lepV[k]) > maxDeltaR) maxDeltaR = lepV[l].DeltaR(lepV[k]);
							if(lepV[l].DeltaR(lepV[k]) < minDeltaR) minDeltaR = lepV[l].DeltaR(lepV[k]);
						}
					}
					if(minDeltaR < 0.05) continue;
					if(maxDeltaR < 2) continue;
					if(!vetoLowMll(5)) continue;
				}				
				if(lepSyst.M() > 80) continue;
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
			//fill Nominal yields
			double mt_min = transmass(lepV[lw_min], METvec);
			unsigned searchR = hnl::sr(mt_min, minMos, lepSyst.M(), cat);
			if(nBJets(true, false, 0) == 0 && (cat > 3 || _met < 75)){
				yields[cat][fill]->Fill(searchR + 1, scal);
			}
			//Following uncertainties do not apply to data or data-driven backgrounds
			if(effsam == 0) continue;
			//yields with JEC varied down
			METvec.SetPtEtaPhiE(_metJECDown, 0, _met_phiJECDown, _metJECDown);	
			mt_min = transmass(lepV[lw_min], METvec);
			searchR = hnl::sr(mt_min, minMos, lepSyst.M(), cat);
			if(nBJets(true, false, 1) == 0 && (cat > 3 || _metJECDown < 75)){
				yieldsDown[0][cat][fill]->Fill(searchR + 1, scal);
			}
			//yields with JEC varied up
			METvec.SetPtEtaPhiE(_metJECUp, 0, _met_phiJECUp, _metJECUp);	
			mt_min = transmass(lepV[lw_min], METvec);
			searchR = hnl::sr(mt_min, minMos, lepSyst.M(), cat);
			if(nBJets(true, false, 2) == 0 && (cat > 3 || _metJECUp < 75)){
				yieldsUp[0][cat][fill]->Fill(searchR + 1, scal);
			}
			//nominal b-veto
			if(nBJets(true, false, 0) != 0) continue;
			//yields with unclustered met varied down
			METvec.SetPtEtaPhiE(_metOtherDown, 0, _met_phiOtherDown, _metOtherDown);			
			mt_min = transmass(lepV[lw_min], METvec);
			searchR = hnl::sr(mt_min, minMos, lepSyst.M(), cat);
			if(cat > 3 || _metOtherDown < 75){
				yieldsDown[1][cat][fill]->Fill(searchR + 1, scal);
			}
			//yields with unclustered met varied up
			METvec.SetPtEtaPhiE(_metOtherUp, 0, _met_phiOtherUp, _metOtherUp);	
			mt_min = transmass(lepV[lw_min], METvec);
			searchR = hnl::sr(mt_min, minMos, lepSyst.M(), cat);
			if(cat > 3 || _metOtherUp < 75){
				yieldsUp[1][cat][fill]->Fill(searchR + 1, scal);
			}	
			//nominal met cut for low-mass categories
			if(cat < 4 && _met > 75) continue;
			//nominal search region
			searchR = hnl::sr(mt_min, minMos, lepSyst.M(), cat);
			//Scale weights for glugluToZZ
			if(fileList[sam].Contains("GluGluToZZTo")){
				for(unsigned pdf = 0; pdf < 110; ++pdf){
					_scaleWeight[pdf] = 1;
				}
			}
			//yields with scale varied down
			yieldsDown[2][cat][fill]->Fill(searchR + 1, scal*_scaleWeight[8]);
			//yields with scale varied up
			yieldsUp[2][cat][fill]->Fill(searchR + 1, scal*_scaleWeight[4]);
			//yields with all pdf variations
			for(unsigned pdf = 0; pdf < 100; ++pdf){
				yieldsPdfVar[pdf][cat][fill]->Fill(searchR + 1, scal*_scaleWeight[pdf + 9]);
			}
			//yields with pileup varied down
			yieldsDown[4][cat][fill]->Fill(searchR + 1, (scal/PUweights[0]->GetBinContent(PUweights[0]->FindBin( std::min(_n_trueInteractions, 49) )))*PUweights[1]->GetBinContent(PUweights[1]->FindBin(std::min(_n_trueInteractions, 49)  )) );
			//yields with pileup varied up
			yieldsUp[4][cat][fill]->Fill(searchR + 1, (scal/PUweights[0]->GetBinContent(PUweights[0]->FindBin( std::min(_n_trueInteractions, 49) )))*PUweights[2]->GetBinContent(PUweights[2]->FindBin(std::min(_n_trueInteractions, 49) )) );
			//yields with b-tag SF varied down
			yieldsDown[5][cat][fill]->Fill(searchR + 1, (scal/bTagSF(true, 0))*bTagSF(true, 1));
			//yields with b-tag SF varied up
			yieldsUp[5][cat][fill]->Fill(searchR + 1, (scal/bTagSF(true, 0))*bTagSF(true, 2));

			delete[] ind;
			delete[] conePt;
			delete[] lepV;
    	}
		//Set negative bins to 0 before adding other processes
		for(unsigned cat = 0; cat < nCat; ++cat){
			for(unsigned bin = 1; bin < yields[cat][effsam]->GetNbinsX() + 1; ++bin){
				if(yields[cat][effsam]->GetBinContent(bin) < 0.) yields[cat][effsam]->SetBinContent(bin, 0.);
				for(unsigned unc = 0; unc < nUnc; ++unc){
					if(yieldsDown[unc][cat][effsam]->GetBinContent(bin) < 0.) yieldsDown[unc][cat][effsam]->SetBinContent(bin, 0.);
					if(yieldsUp[unc][cat][effsam]->GetBinContent(bin) < 0.) yieldsUp[unc][cat][effsam]->SetBinContent(bin, 0.);
				}
				for(unsigned pdf = 0; pdf < 100; ++pdf){
					if(yieldsPdfVar[pdf][cat][effsam]->GetBinContent(bin) < 0.) yieldsPdfVar[pdf][cat][effsam]->SetBinContent(bin, 0.);
				}
			}
		} 
	}
	//Set fakes to 0 if they are negative
	for(unsigned cat = 0; cat < nCat; ++cat){
		for(unsigned bin = 1; bin < yields[cat][nSamples_eff]->GetNbinsX() + 1; ++bin){
			if(yields[cat][nSamples_eff]->GetBinContent(bin) < 0) yields[cat][nSamples_eff]->SetBinContent(bin, 0.);
		}
	}
	//Calculate rms of pdf shape variations
	for(unsigned cat = 0; cat < nCat; ++cat){
		for(unsigned effsam = 1; effsam < nSamples_eff; ++effsam){
			for(unsigned b = 1; b < yields[cat][effsam]->GetNbinsX() + 1; ++b){
				double pdfVarRms = 0;
				for(unsigned pdf = 0; pdf < 100; ++pdf){
					pdfVarRms += (yields[cat][effsam]->GetBinContent(b) - yieldsPdfVar[pdf][cat][effsam]->GetBinContent(b))*(yields[cat][effsam]->GetBinContent(b) - yieldsPdfVar[pdf][cat][effsam]->GetBinContent(b));
				}
				pdfVarRms = 0.01*sqrt(pdfVarRms);
				yieldsDown[3][cat][effsam]->SetBinContent(b, yields[cat][effsam]->GetBinContent(b) - pdfVarRms);
				yieldsUp[3][cat][effsam]->SetBinContent(b, yields[cat][effsam]->GetBinContent(b) + pdfVarRms);
			}
		}
	}
	//Calculate flat acceptance uncertainties
	/*
	//scale variations
	for(unsigned cat = 0; cat < nCat; ++cat){
		for(unsigned effsam = 1; effsam < nSamples_eff; ++effsam){
			yieldsDown[6][cat][effsam] = (TH1D*) yields[cat][effsam]->Clone();
			yieldsDown[6][cat][effsam]->Scale(pdfCounter[effsam][8]/hcounter[effsam]);
			yieldsUp[6][cat][effsam] = (TH1D*) yields[cat][effsam]->Clone();
			yieldsUp[6][cat][effsam]->Scale(pdfCounter[effsam][4]/hcounter[effsam]);
		}
	}
	//rms of pdf acceptance variations
	double pdfAccUnc[nSamples_eff];
	for(unsigned effsam = 1; effsam < nSamples_eff; ++effsam){
		pdfAccUnc[effsam] = 0;
		for(unsigned pdf = 0; pdf < 100; ++pdf){
		 	pdfAccUnc[effsam] += (pdfCounter[effsam][pdf + 9]/hcounter[effsam])*(pdfCounter[effsam][pdf + 9]/hcounter[effsam]);
		}
		pdfAccUnc[effsam] = 0.01*sqrt(pdfAccUnc[effsam]);
	}
	*/
	//Split data and MC histograms for plotting and propagating uncertainties
	TH1D* dataYields[nCat];
	for(unsigned cat = 0; cat < nCat; ++cat){
		//dataYields[cat] = (TH1D*) yields[cat][nSig + 1]->Clone();
		dataYields[cat] = (TH1D*) yields[cat][0]->Clone();
	}
	TH1D* bkgYields[nCat][nSamples_eff -nSig];
	TH1D* bkgErros[nCat][nSamples_eff - nSig];
	for(unsigned cat = 0; cat < nCat; ++cat){
		for(unsigned effsam = nSig + 1; effsam < nSamples_eff + 1; ++effsam){
			bkgYields[cat][effsam -nSig - 1] = (TH1D*) yields[cat][effsam]->Clone();
			/*
			if(effsam > nSig + 1){
				dataYields[cat]->Add(bkgYields[cat][effsam - nSig - 1]);
			}
			*/
		}
	}
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//PRINT DATA CARDS FOR LIMIT CALCULATION
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	const unsigned nBkg = nSamples_eff - nSig;
	const unsigned nSyst = 9 + 1 + 2*nBkg; //9 general uncertainties + stat signal + stat bkg + extra unc per bkg
	TString systNames[nSyst] = {"lumi", "id_eff", "trigeff", "JEC", "metUncl", "scale", "pdf", "pu", "btagSF"};
	const TString sigN[nSig] = {"m1", "m2", "m5", "m10", "m20", "m30", "m40", "m50", "m60", "m80", "m100", "m130", "m150", "m200", "m400", "m600", "m800", "m1000"};
	const TString bkgNames[nBkg] = {"ZZH", "triboson", "WZ", "Xgamma", "TTX",  "nonPrompt"}; //slightly rewrite bkg names not to confuse the combine tool
	std::vector<std::vector<double>> systUnc(nSyst, std::vector<double>(nBkg + 1, 0)); //2D array containing all systematics
	//initialize flat and shape systematics
	for(unsigned proc = 0; proc < nBkg; ++proc){ //last background is non-prompt and isn't susceptible to these uncertainty sources (note that entry 0 is signal so that the loop goes up to the second to last bkg)
		systUnc[0][proc] = 1.025;	//lumi
		systUnc[1][proc] = 1.04;	//id eff
		systUnc[2][proc] = 1.02;	//trig eff
		systUnc[3][proc] = 1.00;	//JEC shape
		systUnc[4][proc] = 1.00;	//PU shape
		systUnc[5][proc] = 1.00;	//MET unclustered shape
		systUnc[6][proc] = 1.00; 	//b-tag SF shape
		systUnc[7][proc] = 1.00;	//scale shape
		systUnc[8][proc] = 1.00;	//pdf shape
		
	}
	TString systDist[nSyst]; //probability distribution of nuisances
	for(unsigned syst = 0; syst < nSyst; ++syst){
		if(syst > 2 && syst < 9) systDist[syst] = "shape";
		else systDist[syst] = "lnN";
	}
	const double extraUnc[nBkg] = {1.25, 1.5, 1.085, 1.15, 1.5, 1.3}; //extra flat uncertainties assigned to each background
	for(unsigned syst = 10 + nBkg; syst < nSyst; ++syst){//loop over last nBkg uncertainties, being the exta uncertainties for each bkg
		unsigned bkg = syst - 10 - nBkg;//index of the background corresponding to the uncertainty index
		systNames[syst] = "extra" + bkgNames[bkg];
		systUnc[syst][bkg + 1] = extraUnc[bkg];
	}
	unsigned binCount = 0;
	for(unsigned cat = 0; cat < nCat; ++cat){
		if(cat == 0 || cat == 2) continue;
		binCount += nSR[cat];
	}
	const unsigned binC = binCount;
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
	for(unsigned sig = 0; sig < nSig; ++sig){
		//Counter to determine the "combined bin index"
		unsigned binCounter = 1;
		//loop over all categories and bins
		for(unsigned cat = 0; cat < nCat; ++cat){
			if(cat == 0 || cat == 2) continue;
			for(unsigned bin = 1;  bin < dataYields[cat]->GetNbinsX() + 1; ++bin){
				//Define correlations between uncertainty sources
				//systNames[4] += catNames[cat];
				if(cat > 4) systNames[4] += catNames[cat]; //unclustered  met is fully correlated in low mass because of upper met cut
				if(cat == 4){
					if(bin < 4) systNames[4] += "part1";
					else if(bin < 9) systNames[4] += "part2";
					else if(bin < 13) systNames[4] += "part3";
					else systNames[4] += "part4";
				} else if(cat == 4){
					if(bin < 3) systNames[4] += "part1";
					else if(bin < 7) systNames[4] += "part2";
					else if(bin < 9) systNames[4] += "part3";
					else systNames[4] += "part4";
				}
				systNames[5] += "bin" + std::to_string(binCounter);	 //Assume scale uncorrelated between bins so the bulk does not constrain the tails
				//systNames[6] += "bin" + std::to_string(binCounter);
				
				//make root file 
				TFile* shapeFile = new TFile("./datacards/shapes/shapeFile_ElectronCoupling_"  + sigN[sig] + "_bin" + (TString) std::to_string(binCounter) + extra + ".root", "recreate");
				histCent[0][binCounter - 1]->SetBinContent(1, dataYields[cat]->GetBinContent(bin));
				histCent[0][binCounter - 1]->Write("data_obs");//_ch"  + (TString) std::to_string(binCounter));
				histCent[sig + 1][binCounter - 1]->SetBinContent(1, std::max(yields[cat][1+sig]->GetBinContent(bin), 0.) );
				histCent[sig + 1][binCounter - 1]->Write("signal");//_ch"  + (TString) std::to_string(binCounter));
				for(unsigned unc = 0; unc < nUnc; ++unc){
					double min = (yields[cat][1+sig]->GetBinContent(bin) == 0) ? 0. : (double) std::numeric_limits< float >::min();
					histUp[unc][sig + 1][binCounter - 1]->SetBinContent(1, std::max(yieldsUp[unc][cat][sig + 1]->GetBinContent(bin), min ));
					histUp[unc][sig + 1][binCounter - 1]->Write("signal_" + systNames[3 + unc] + "Up");//"signal_ch"  + (TString) std::to_string(binCounter) + "_" + systNames[3 + unc] + "Up"
					histDown[unc][sig + 1][binCounter - 1]->SetBinContent(1, std::max(yieldsDown[unc][cat][sig + 1]->GetBinContent(bin), min ));
					histDown[unc][sig + 1][binCounter - 1]->Write("signal_" + systNames[3 + unc] + "Down"); //"signal_ch"  + (TString) std::to_string(binCounter) + "_" + systNames[3 + unc] + "Down");
				}
				for(unsigned bkg = 0; bkg < nBkg; ++bkg){ //BUG IN BIN INDEX!
					histCent[bkg + nSig + 1][binCounter - 1]->SetBinContent(1, std::max(yields[cat][bkg + nSig + 1]->GetBinContent(bin), 0.)  );
					histCent[bkg + nSig + 1][binCounter - 1]->Write(bkgNames[bkg]); //+ "_Ch" + std::to_string(binCounter));
					if(bkg != nBkg -1){
						for(unsigned unc = 0; unc < nUnc; ++unc){						
							double min = (yields[cat][bkg + nSig + 1]->GetBinContent(bin) == 0) ? 0. : (double) std::numeric_limits< float >::min();
							histUp[unc][bkg + nSig + 1][binCounter - 1]->SetBinContent(1, std::max((double) yieldsUp[unc][cat][bkg + nSig + 1]->GetBinContent(bin), min ) );
							histDown[unc][bkg + nSig + 1][binCounter - 1]->SetBinContent(1, std::max((double) yieldsDown[unc][cat][bkg + nSig + 1]->GetBinContent(bin), min ) );
							histUp[unc][bkg + nSig + 1][binCounter - 1]->Write(bkgNames[bkg] + "_" + systNames[3 + unc] + "Up");//bkgNames[bkg] + "_ch"  + (TString) std::to_string(binCounter) + "_" + systNames[3 + unc] + "Up"
							histDown[unc][bkg + nSig + 1][binCounter - 1]->Write(bkgNames[bkg] + "_" + systNames[3 + unc] + "Down");//bkgNames[bkg] + "_ch"  + (TString) std::to_string(binCounter)+ "_" +  systNames[3 + unc] + "Down"
						}
					}
				}
				shapeFile->Close();
				//signal stat unc
				systUnc[9][0] = 1 + std::max(0., yields[cat][sig + 1]->GetBinError(bin)/yields[cat][sig + 1]->GetBinContent(bin));
				systNames[9] = "statSig" + std::to_string(binCounter);
				double bkgYieldVal[nBkg];//bkg rates
				//bkg stat unc
				for(unsigned bkg = 0; bkg < nBkg; ++bkg){
					bkgYieldVal[bkg] = bkgYields[cat][bkg]->GetBinContent(bin);
					systUnc[10 + bkg][bkg + 1]  = 1 + ((bkgYields[cat][bkg]->GetBinContent(bin) == 0) ? 0. : std::max(0., bkgYields[cat][bkg]->GetBinError(bin)/bkgYields[cat][bkg]->GetBinContent(bin)) );
					systNames[10 + bkg] = "stat"  + bkgNames[bkg] + std::to_string(binCounter);
				}
				//change shape uncertainty names for uncorrelated uncertainty sources
				//print datacard
				hnl::printDataCard(dataYields[cat]->GetBinContent(bin), yields[cat][sig + 1]->GetBinContent(bin), "signal", bkgYieldVal, nBkg, bkgNames, systUnc, nSyst, systNames, systDist, "datacards/datacard_ElectronCoupling" + sigN[sig] + "_bin" + (TString) std::to_string(binCounter) + extra, true, "shapes/shapeFile_ElectronCoupling_" + sigN[sig] + "_bin" + (TString) std::to_string(binCounter) + extra);
				//Reset uncertainty names that were changed for correlations
				systNames[4] = "metUncl";
				systNames[5] =  "scale";
				//systNames[6] =  "pdf";
				//increment bin Counter
				++binCounter;
			}
		}
	}
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//Calculate histogram with systematic uncertainty for backgrounds and print the ranges of their size
	double flatSyst[4] = {0.025, 0.06, 0.02};
	TH1D* bkgSyst[nCat][nSamples_eff -nSig];
	for(unsigned bkg = 0; bkg < nBkg; ++bkg){
		//flatSyst[4] = fabs(pdfAccUnc[bkg + nSig + 1] - 1);
		flatSyst[3] = fabs(extraUnc[bkg] -1);
		for(unsigned cat = 0; cat < nCat; ++cat){
			bkgSyst[cat][bkg] = (TH1D*) bkgYields[cat][bkg]->Clone();
			for(unsigned b = 1;  b < dataYields[cat]->GetNbinsX() + 1; ++b){
				bkgSyst[cat][bkg]->SetBinContent(b,0);
				/*
				if(yields[cat][bkg + nSig + 1]->GetBinContent(b) == 0){
					bkgSyst[cat][bkg]->SetBinContent(b,0);
				} else{
				*/
					if(bkg != nBkg -1){
						double systbin = 0;
						//loop over shape uncertainties
						for(unsigned unc = 0; unc < nUnc; ++unc){
							double syst = std::max(fabs(yieldsUp[unc][cat][bkg + nSig + 1]->GetBinContent(b) - yields[cat][bkg + nSig + 1]->GetBinContent(b)), fabs(yieldsDown[unc][cat][bkg + nSig + 1]->GetBinContent(b) - yields[cat][bkg + nSig + 1]->GetBinContent(b)));
							systbin += syst*syst;
							
						}
						//loop over flat uncertainties
						for(unsigned unc = 0; unc < 4; ++unc){
							systbin += bkgYields[cat][bkg]->GetBinContent(b)*bkgYields[cat][bkg]->GetBinContent(b)*flatSyst[unc]*flatSyst[unc];
						}
						bkgSyst[cat][bkg]->SetBinContent(b, sqrt(systbin));
					} else{
						bkgSyst[cat][bkg]->SetBinContent(b, yields[cat][bkg + nSig + 1]->GetBinContent(b)*0.3);
					}
				//}
				bkgSyst[cat][bkg]->SetBinError(b, 0);
			}
		}
	}

	const TString backGroundNames[nSamples_eff] = {"ZZ/H", "triboson", "WZ", "X + #gamma", "TT/T + X",  "non-prompt"};
	for(unsigned cat = 0; cat < nCat; ++cat){
		hnl::makeTable(cat, bkgYields[cat], bkgSyst[cat], dataYields[cat],  backGroundNames, nSamples_eff -nSig, catNames[cat], "_ElectronCoupling");
	}


	//const TString distNames[nSamples_eff + 1 - nSig] = {"total pred.", "ZZ/H", "triboson", "WZ", "X + #gamma", "TT/T + X",  "non-prompt"};
	const TString distNames[nSamples_eff + 1 - nSig] = {"obs.", "ZZ/H", "triboson", "WZ", "X + #gamma", "TT/T + X",  "non-prompt"};
	//Plot the yields as a function of the search region
	for(unsigned cat = 0; cat < nCat; ++cat){
		const unsigned nSigToPlot = 6;
		TH1D* signals[nSigToPlot];
		TString sigPlotNames[nSigToPlot];
		if(cat < 4){
			sigPlotNames[0] = "1 GeV, |V|^{2} = 10^{-5} ";
			sigPlotNames[1] = "5 GeV, |V|^{2} = 10^{-5} ";
			sigPlotNames[2] = "20 GeV, |V|^{2} = 10^{-5} ";
			sigPlotNames[3] = "30 GeV, |V|^{2} = 10^{-5} ";
			sigPlotNames[4] = "50 GeV, |V|^{2} = 10^{-5} ";
			sigPlotNames[5] = "60 GeV, |V|^{2} = 10^{-5} ";
		} else{
			sigPlotNames[0] = "100 GeV, |V|^{2} = 10^{-2} ";
			sigPlotNames[1] = "130 GeV, |V|^{2} = 10^{-2} ";
			sigPlotNames[2] = "150 GeV, |V|^{2} = 10^{-2} ";
			sigPlotNames[3] = "200 GeV, |V|^{2} = 10^{-2} ";
			sigPlotNames[4] = "400 GeV, |V|^{2} = 10^{-2} ";
			sigPlotNames[5] = "600 GeV, |V|^{2} = 10^{-2} ";
		}
		//loop over all kinematic distributions
		if(cat < 4){
			signals[0] = yields[cat][1];
			signals[1] = yields[cat][3];
			signals[2] = yields[cat][5];
			signals[3] = yields[cat][6];
			signals[4] = yields[cat][8];
			signals[5] = yields[cat][9];
		} else{
			signals[0] = yields[cat][11];
			signals[1] = yields[cat][12];
			signals[2] = yields[cat][13];
			signals[3] = yields[cat][14];
			signals[4] = yields[cat][15];
			signals[5] = yields[cat][16];
		}
		plotDataVSMC(dataYields[cat], bkgYields[cat], distNames, nSamples_eff - nSig, catNames[cat] + "_ElectronCoupling_withSignal" + extra, true, 0, "HNL", bkgSyst[cat], true, signals,  sigPlotNames ,nSigToPlot);
		plotDataVSMC(dataYields[cat], bkgYields[cat], distNames, nSamples_eff - nSig, catNames[cat] + "_ElectronCoupling" + extra, true, 0, "HNL", bkgSyst[cat]);
	}
}



int main(int argc, char* argv[]){
	TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;
	testtree.Loop();
	rootapp->Run();
    return 0;
}


