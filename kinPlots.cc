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
	const double couplingCorrection = 5000;
	//~~~~~~~~ background normalizations~~~~~~~~~~~~~~~~~~~~
	const double glugluToZZkFactor = 2.1;
	const double WZSF = 0.652; //0.655
	const double ZZSF = 1.029; //1.032
	const double XgammaSF = 0.950; //0.948
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	const TString fileList[nSamples] = {"data_combined_trilepton.root", 
										"HeavyNeutrino_M1_2l.root", "HeavyNeutrino_M2_2l.root", "HeavyNeutrino_M5_2l.root", "HeavyNeutrino_M10_2l.root", "HeavyNeutrino_M20_2l.root", "HeavyNeutrino_M30_2l.root", "HeavyNeutrino_M40_2l.root", "HeavyNeutrino_M50_2l.root", "HeavyNeutrino_M60_2l.root", "HeavyNeutrino_M80_2l.root", "HeavyNeutrino_M100_2l.root", "HeavyNeutrino_M130_2l.root", "HeavyNeutrino_M150_2l.root", "HeavyNeutrino_M200_2l.root", "HeavyNeutrino_M400_2l.root", "HeavyNeutrino_M600_2l.root", "HeavyNeutrino_M800_2l.root", "HeavyNeutrino_M1000_2l.root",
 										"ZZTo4L.root",  "GluGluToZZTo4mu.root", "GluGluToZZTo4e.root", "GluGluToZZTo4tau.root", "GluGluToZZTo2e2mu.root", "GluGluToZZTo2e2tau.root", "GluGluToZZTo2mu2tau.root", "VHToNonbb.root", "GluGluHToZZTo4L_M125.root", "VBF_HToZZTo4L_M125.root", "WWG.root","WWW.root", "WWZ.root", "WWTo2L2Nu_DoubleScattering.root", "WWTo2L2Nu.root",  "ZZZ.root", "WZTo3LNu_mllmin01.root", "TTGJets.root","ZGTo2LG.root", "WGToLNuG.root", "TGJets.root", "TTJets_DiLept.root", "TTJets_SingleLeptFromTbar.root", "TTJets_SingleLeptFromT.root", "DYJetsToLL_M10to50.root", "DYJetsToLL_M50.root", "ttHToNonbb.root", "TTWJetsToLNu.root", "TTZToLLNuNu.root", "TTZToLL_M1to10.root",  "TTTT.root"};
	
	const double xSections[nSamples - 1] = {5.237e-01*couplingCorrection*lowMCoupling, 5.238e-01*couplingCorrection*lowMCoupling, 5.215e-01*couplingCorrection*lowMCoupling, 5.132e-01*couplingCorrection*lowMCoupling, 4.764e-01*couplingCorrection*lowMCoupling, 4.160e-01*couplingCorrection*lowMCoupling, 3.335e-01*couplingCorrection*lowMCoupling, 2.364e-01*couplingCorrection*lowMCoupling, 1.336e-01*couplingCorrection*lowMCoupling, 2.473e-02*couplingCorrection*lowMCoupling*10, 1.033e-03*couplingCorrection*highMCoupling, 2.874e-04*couplingCorrection*highMCoupling, 1.561e-04*couplingCorrection*highMCoupling, 4.929e-05*couplingCorrection*highMCoupling, 3.505e-06*couplingCorrection*highMCoupling, 6.922e-07*couplingCorrection*highMCoupling, 2.013e-07*couplingCorrection*highMCoupling, 7.239e-08*couplingCorrection*highMCoupling, 
											1.256*ZZSF, 0.00159*glugluToZZkFactor*ZZSF, 0.00159*glugluToZZkFactor*ZZSF, 0.00159*glugluToZZkFactor*ZZSF, 0.00319*glugluToZZkFactor*ZZSF, 0.00319*glugluToZZkFactor*ZZSF, 0.00319*glugluToZZkFactor*ZZSF,  0.9561, 0.01212, 0.001034,  0.2147, 0.2086, 0.1651, 0.1729,   12.178, 0.01398, 58.59*WZSF, 3.697, 123.9*XgammaSF, 489*XgammaSF,  2.967, 87.315, 182.175, 182.175, 18610, 1921.8*3,  0.215, 0.2043, 0.2529, 0.0493, 0.009103};
	const TString names[nSamples] = {"data", "m_{N} = 1 GeV", "m_{N} = 2 GeV", "m_{N} = 5 GeV", "m_{N} = 10 GeV", "m_{N} = 20 GeV", "m_{N} = 30 GeV", "m_{N} = 40 GeV", "m_{N} = 50 GeV", "m_{N} = 60 GeV", "m_{N} = 80 GeV", "m_{N} = 100 GeV", "m_{N} = 130 GeV", "m_{N} = 150 GeV", "m_{N} = 200 GeV", "m_{N} = 400 GeV", "m_{N} = 600 GeV", "m_{N} = 800 GeV", "m_{N} = 1000 GeV",  									 
									 "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "triboson", "triboson", "triboson", "triboson", "triboson", "triboson", "WZ", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X"};

	
	//Read Trees from ROOT files
	TFile* hfile[nSamples];
	TTree* inputTree[nSamples];
	double hcounter[nSamples];
	double pdfCounter[nSig][110];
	for(unsigned sam = 0; sam < nSamples; ++sam){
		cout << "name " << names[sam] << endl;
		hfile[sam] = new TFile("../data_april17/"+fileList[sam],"read");
       	hfile[sam]->cd("FakeElectrons");
		//Determine hcounter for cross section scaling
		TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
		_hCounter->Read("hCounter");
		hcounter[sam] = _hCounter->GetBinContent(1);
		if(sam > 0 && sam <= nSig){
			TH1D* _pdfCounter = new TH1D("pdfCounter", "Events counter", 110, 0, 110);
			_pdfCounter->Read("pdfCounter");
			for(unsigned pdf = 0; pdf < 110; ++pdf){
				pdfCounter[sam -1][pdf] = _pdfCounter->GetBinContent(pdf + 1);
			}
		}
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
	const unsigned nCat = 8;  //Number of categories
	const TString catNames[nCat] = {"lowM_3lOSSF_lowPt", "lowM_3lnoOSSF_lowPt", "lowM_3lOSSF_highPt", "lowM_3lnoOSSF_highPt", "highM_3lOSSF", "highM_3lnoOSSF", "baseline_noOSSF", "baseline_OSSF"};
	const unsigned nUnc = 10;//number of shape uncertainties
	const TString uncNames[nUnc] = {"jec", "metUnclustered",  "scale", "pdf", "pu", "btagSF", "trigeff", "id_eff", "fakeEWK", "ZZmt"};
	const unsigned nDist = 5;
	const TString distNames[nDist] = {"M3l", "minMos", "mt_minMos", "MET", "ConePt_le"};
	const TString xAxes[nDist] = {"M_{3l} (GeV)", "M_{2lOS}^{min} (GeV)", "M_{T} (GeV)", "E_{T}^{miss} (GeV)", "P_{T}^{Cone}(leading) (GeV)"};
	double histMin[nDist] = {0, 0, 0, 0, 15};
	double histMax[nDist] = {600, 300, 300, 300, 215};
	unsigned nBins[nDist] = {20, 20, 20, 20, 20};

	TH1D* histos[nDist][nCat][nSamples_eff + 1]; //nominal yields in every SR
	TH1D* histosDown[nUnc][nDist][nCat][nSamples_eff + 1]; //yields varied down by shape unc
	TH1D* histosUp[nUnc][nDist][nCat][nSamples_eff + 1]; //yields varied up by shape unc
	TH1D* histosPdfVar[100][nDist][nCat][nSamples_eff + 1]; //yields for all 100 possible pdf variations

	for(unsigned cat = 0; cat < nCat; ++cat){
		if(cat < 4){
			histMax[3] = 75;
			nBins[3] = 10;
			histMax[0] = 80;
			nBins[0] = 10;
			histMax[1] = 80;
			nBins[1] = 10;
		}
		if(cat == 7){
			histMax[0] = 300;
			nBins[0] = 100;
		}
		for(unsigned effsam = 0; effsam < nSamples_eff + 1; ++effsam){
			for(unsigned dist = 0; dist < nDist; ++dist){
				histos[dist][cat][effsam] = new TH1D(distNames[dist] + catNames[cat] + eff_names[effsam], distNames[dist] + catNames[cat] + eff_names[effsam] + ";" + xAxes[dist] + "; Events", nBins[dist], histMin[dist],histMax[dist]);	
				histos[dist][cat][effsam]->Sumw2(); 
				for(unsigned unc = 0; unc < nUnc; ++unc){
					histosDown[unc][dist][cat][effsam] = new TH1D(distNames[dist] + catNames[cat]  +  eff_names[effsam] + uncNames[unc] + "Up", distNames[dist] + catNames[cat] + eff_names[effsam] + uncNames[unc] + "Up;"+ xAxes[dist] + ";Events", nBins[dist], histMin[dist],histMax[dist]);
					histosDown[unc][dist][cat][effsam]->Sumw2();
					histosUp[unc][dist][cat][effsam] = new  TH1D(distNames[dist] + catNames[cat]  +  eff_names[effsam]  + uncNames[unc] + "Down", distNames[dist] + catNames[cat] + eff_names[effsam]  + uncNames[unc] + "Down;" + xAxes[dist] + ";Events", nBins[dist], histMin[dist],histMax[dist]);
					histosUp[unc][dist][cat][effsam]->Sumw2();
				}
				for(unsigned pdf = 0; pdf < 100; ++pdf){
					histosPdfVar[pdf][dist][cat][effsam] = new  TH1D(distNames[dist] + catNames[cat]  +  eff_names[effsam]  + "_pdf" + std::to_string(pdf), distNames[dist] + catNames[cat] + eff_names[effsam]  + "_pdf" + std::to_string(pdf) + ";" + xAxes[dist] +  ";Events", nBins[dist], histMin[dist],histMax[dist]);
					histosPdfVar[pdf][dist][cat][effsam]->Sumw2();
				}
			}
		}
		if(cat < 4){
			histMax[3] = 300;
			nBins[3] = 20;
			histMax[0] = 300;
			nBins[0] = 20;
			histMax[1] = 300;
			nBins[1] = 20;
		}
		if(cat == 7){
			histMax[0] = 600;
			nBins[0] = 20;
		}
	}
    Double_t scale[nSamples -1];

	double maxBinC[nDist][nCat];
	for(unsigned dist = 0; dist < nDist; ++dist){
		for(unsigned cat = 0; cat < nCat; ++cat){
			maxBinC[dist][cat] = histos[dist][cat][0]->GetBinCenter(histos[dist][cat][0]->GetNbinsX());
		}
	}

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
				scal*= fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[0], lCount);
			} else if(tightFail) continue;

			//Sample overlap removal
			if(fileList[sam] == "WGToLNuG.root"){
				bool promptfail = true;
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

			//check if category specific selection is passed
			bool catPass = true;
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
					if(minDeltaR < 0.05) catPass = false;
					if(maxDeltaR < 2) catPass = false;
					//if(!vetoLowMll(5)) catPass = false;
					if(!vetoLowMll(5)) continue;
				}				
				if(lepSyst.M() > 80) catPass = false;
			}
			else if(cat > 3){
				if(conePt[1] < 15) catPass = false;
				if(conePt[2] < 10) catPass = false;
				if(cat == 4){
					if(fabs(mll - 91) < 15) catPass = false; //Veto onZ events (WZ and DY)
					if(fabs(lepSyst.M() - 91) < 15) catPass = false; //veto conversions
				}
				if(!vetoLowMll(5)) catPass = false;
			}



			bool ossfI = (hnl::flavorComposition(ind, _flavors, _charges, lCount) == 0);
			double values[nDist] = {lepSyst.M(), minMos, transmass(lepV[lw_min], METvec), _met, conePt[0]};
			if(nBJets(true, false, 0) == 0){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histos[dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), scal);
				}
				if(catPass && (cat > 3 || _met < 75)){
					for(unsigned dist = 0; dist < nDist; ++dist){
						histos[dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), scal);
					}
				}
			}
			
			//FR EWK contamination uncertainty
			if(fill == nSamples_eff && nBJets(true, false, 0) == 0){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosDown[8][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), scal*(fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[1], lCount)/fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[0], lCount) ));
					histosUp[8][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), scal*(fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[2], lCount)/fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[0], lCount) ));
				}
				if(catPass && (cat > 3 || _met< 75)){
					for(unsigned dist = 0; dist < nDist; ++dist){
						histosDown[8][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), scal*(fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[1], lCount)/fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[0], lCount) ));
						histosUp[8][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), scal*(fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[2], lCount)/fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[0], lCount) ));
					}
				}
			}
			
			//Following uncertainties do not apply to data or data-driven backgrounds
			if(effsam <= nSig) continue;
			//yields with JEC varied down
			METvec.SetPtEtaPhiE(_metJECDown, 0, _met_phiJECDown, _metJECDown);	
			values[2] = transmass(lepV[lw_min], METvec);
			values[3] = _metJECDown;
			if(nBJets(true, false, 1) == 0){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosDown[0][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), scal);
				}
				if(catPass && (cat > 3 || _metJECDown < 75)){
					for(unsigned dist = 0; dist < nDist; ++dist){
						histosDown[0][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), scal);
					}
				}
			}
			//yields with JEC varied up
			METvec.SetPtEtaPhiE(_metJECUp, 0, _met_phiJECUp, _metJECUp);	
			values[2] = transmass(lepV[lw_min], METvec);
			values[3] = _metJECUp;
			if(nBJets(true, false, 2) == 0){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosUp[0][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), scal);
				}
				if(catPass && (cat > 3 || _metJECUp < 75)){
					for(unsigned dist = 0; dist < nDist; ++dist){
						histosUp[0][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), scal);
					}
				}
			}
			//nominal b-veto
			if(nBJets(true, false, 0) != 0) continue;
			//yields with unclustered met varied down
			METvec.SetPtEtaPhiE(_metOtherDown, 0, _met_phiOtherDown, _metOtherDown);			
			values[2] = transmass(lepV[lw_min], METvec);
			values[3] = _metOtherDown;
			for(unsigned dist = 0; dist < nDist; ++dist){
				histosDown[1][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), scal);
			}
			if(catPass && (cat > 3 || _metOtherDown < 75)){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosDown[1][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), scal);
				}
			}
			//yields with unclustered met varied up
			METvec.SetPtEtaPhiE(_metOtherUp, 0, _met_phiOtherUp, _metOtherUp);	
			values[2] = transmass(lepV[lw_min], METvec);
			values[3] = _metOtherUp;
			for(unsigned dist = 0; dist < nDist; ++dist){
				histosUp[1][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), scal);
			}
			if(catPass && (cat > 3 || _metOtherUp < 75)){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosUp[1][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), scal);
				}
			}	
			//nominal met cut for low-mass categories
			if(cat < 4 && _met > 75) catPass = false;
			//nominal kinematics
			METvec.SetPtEtaPhiE(_met, 0, _met_phi, _met);	
			values[2] = transmass(lepV[lw_min], METvec);
			values[3] = _metOtherUp;
			//Scale weights for glugluToZZ
			if(fileList[sam].Contains("GluGluToZZTo")){
				for(unsigned pdf = 0; pdf < 110; ++pdf){
					_scaleWeight[pdf] = 1;
				}
			}
			//yields with scale varied down
			for(unsigned dist = 0; dist < nDist; ++dist){
				histosDown[2][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), scal*_scaleWeight[8]);
			}
			if(catPass){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosDown[2][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), scal*_scaleWeight[8]);
				}
			}	
			//yields with scale varied up
			for(unsigned dist = 0; dist < nDist; ++dist){
				histosUp[2][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), scal*_scaleWeight[4]);
			}
			if(catPass){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosUp[2][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), scal*_scaleWeight[4]);
				}
			}
			//yields with all pdf variations
			for(unsigned dist = 0; dist < nDist; ++dist){
				for(unsigned pdf = 0; pdf < 100; ++pdf){
					histosPdfVar[pdf][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), scal*_scaleWeight[pdf + 9]);
				}
			}
			if(catPass){
				for(unsigned dist = 0; dist < nDist; ++dist){
					for(unsigned pdf = 0; pdf < 100; ++pdf){
						histosPdfVar[pdf][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), scal*_scaleWeight[pdf + 9]);
					}
				}
			}
			//yields with pileup varied down
			for(unsigned dist = 0; dist < nDist; ++dist){
				histosDown[4][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), (scal/PUweights[0]->GetBinContent(PUweights[0]->FindBin( std::min(_n_trueInteractions, 49) )))*PUweights[1]->GetBinContent(PUweights[1]->FindBin(std::min(_n_trueInteractions, 49)  )) );
			}
			if(catPass){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosDown[4][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), (scal/PUweights[0]->GetBinContent(PUweights[0]->FindBin( std::min(_n_trueInteractions, 49) )))*PUweights[1]->GetBinContent(PUweights[1]->FindBin(std::min(_n_trueInteractions, 49)  )) );
				}
			}
			//yields with pileup varied up
			for(unsigned dist = 0; dist < nDist; ++dist){
				histosUp[4][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), (scal/PUweights[0]->GetBinContent(PUweights[0]->FindBin( std::min(_n_trueInteractions, 49) )))*PUweights[2]->GetBinContent(PUweights[2]->FindBin(std::min(_n_trueInteractions, 49) )) );
			}
			if(catPass){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosUp[4][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), (scal/PUweights[0]->GetBinContent(PUweights[0]->FindBin( std::min(_n_trueInteractions, 49) )))*PUweights[2]->GetBinContent(PUweights[2]->FindBin(std::min(_n_trueInteractions, 49) )) );
				}
			}
			//yields with b-tag SF varied down
			for(unsigned dist = 0; dist < nDist; ++dist){
				histosDown[5][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), (scal/bTagSF(true, 0))*bTagSF(true, 1));
			}
			if(catPass){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosDown[5][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), (scal/bTagSF(true, 0))*bTagSF(true, 1));
				}
			}
			//yields with b-tag SF varied up
			for(unsigned dist = 0; dist < nDist; ++dist){
				histosUp[5][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), (scal/bTagSF(true, 0))*bTagSF(true, 2));
			}
			if(catPass){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosUp[5][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), (scal/bTagSF(true, 0))*bTagSF(true, 2));
				}
			}
			//Trigger uncertainties
			if(conePt[0] > 30){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosDown[6][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), scal*0.98);
					histosUp[6][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), scal*1.02);
				}
				if(catPass){
					for(unsigned dist = 0; dist < nDist; ++dist){
						histosDown[6][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), scal*0.98);
						histosUp[6][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), scal*1.02);
					}
				}
			} else{
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosDown[6][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), scal*0.95);
					histosUp[6][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), scal*1.05);
				}
				if(catPass){
					for(unsigned dist = 0; dist < nDist; ++dist){
						histosDown[6][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), scal*0.95);
						histosUp[6][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), scal*1.05);
					}
				}
			}
			//Id efficiency uncertainties
			unsigned flavorC = tril_flavorComb(ind, _flavors, lCount);
			if(flavorC == 0 || flavorC == 3){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosDown[7][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), scal*0.94);
					histosUp[7][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), scal*1.06);
				}
				if(catPass){
					for(unsigned dist = 0; dist < nDist; ++dist){
						histosDown[7][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), scal*0.94);
						histosUp[7][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), scal*1.06);
					}
				}
			} else{
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosDown[7][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), scal*0.9553);
					histosUp[7][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), scal*1.0447);
				}
				if(catPass){
					for(unsigned dist = 0; dist < nDist; ++dist){
						histosDown[7][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), scal*0.9553);
						histosUp[7][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), scal*1.0447);
					}
				}
			}
			if(effsam - nSig -1 == 0){// only ZZ needs last nuisance!
				if(values[2] < 75){
					for(unsigned dist = 0; dist < nDist; ++dist){
						histosDown[9][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), scal);
						histosUp[9][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), scal);
					}
					if(catPass){
						for(unsigned dist = 0; dist < nDist; ++dist){
							histosDown[9][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), scal);
							histosUp[9][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), scal);
						}
					}
				} else{
					for(unsigned dist = 0; dist < nDist; ++dist){
						histosDown[9][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), scal*0.75);
						histosUp[9][dist][6 + ossfI][fill]->Fill(std::min(values[dist], maxBinC[dist][6 + ossfI]), scal*1.25);
					}
					if(catPass){
						for(unsigned dist = 0; dist < nDist; ++dist){
							histosDown[9][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), scal*0.75);
							histosUp[9][dist][cat][fill]->Fill(std::min(values[dist], maxBinC[dist][cat]), scal*1.25);
						}
					}
				}					
			}
			delete[] ind;
			delete[] conePt;
			delete[] lepV;
    	}
		//Set negative bins to 0 before adding other processes
		for(unsigned cat = 0; cat < nCat; ++cat){
			for(unsigned dist = 0; dist < nDist; ++dist){
				for(unsigned bin = 1; bin < histos[dist][cat][effsam]->GetNbinsX() + 1; ++bin){
					if(histos[dist][cat][effsam]->GetBinContent(bin) < 0.)histos[dist][cat][effsam]->SetBinContent(bin, 0.);
					for(unsigned unc = 0; unc < nUnc; ++unc){
						if(histosDown[unc][dist][cat][effsam]->GetBinContent(bin) < 0.) histosDown[unc][dist][cat][effsam]->SetBinContent(bin, 0.);
						if(histosUp[unc][dist][cat][effsam]->GetBinContent(bin) < 0.) histosUp[unc][dist][cat][effsam]->SetBinContent(bin, 0.);
					}
					for(unsigned pdf = 0; pdf < 100; ++pdf){
						if(histosPdfVar[pdf][dist][cat][effsam]->GetBinContent(bin) < 0.) histosPdfVar[pdf][dist][cat][effsam]->SetBinContent(bin, 0.);
					}
				}
			}
		} 
	}
	//Set fakes to 0 if they are negative
	for(unsigned cat = 0; cat < nCat; ++cat){	
		for(unsigned dist = 0; dist < nDist; ++dist){
			for(unsigned bin = 1; bin < histos[dist][cat][nSamples_eff]->GetNbinsX() + 1; ++bin){
				if(histos[dist][cat][nSamples_eff]->GetBinContent(bin) < 0) histos[dist][cat][nSamples_eff]->SetBinContent(bin, 0.);
				if(histosDown[8][dist][cat][nSamples_eff]->GetBinContent(bin) < 0) histosDown[8][dist][cat][nSamples_eff]->SetBinContent(bin, 0.);
				if(histosUp[8][dist][cat][nSamples_eff]->GetBinContent(bin) < 0) histosUp[8][dist][cat][nSamples_eff]->SetBinContent(bin, 0.);
			}
		}
	}
	
	//Calculate rms of pdf shape variations
	for(unsigned cat = 0; cat < nCat; ++cat){
		for(unsigned effsam = 1; effsam < nSamples_eff; ++effsam){
			for(unsigned dist = 0; dist < nDist; ++dist){
				for(unsigned b = 1; b < histos[dist][cat][effsam]->GetNbinsX() + 1; ++b){
					double pdfVarRms = 0;
					for(unsigned pdf = 0; pdf < 100; ++pdf){
						pdfVarRms += (histos[dist][cat][effsam]->GetBinContent(b) - histosPdfVar[pdf][dist][cat][effsam]->GetBinContent(b))*(histos[dist][cat][effsam]->GetBinContent(b) - histosPdfVar[pdf][dist][cat][effsam]->GetBinContent(b));
					}
					pdfVarRms = 0.01*sqrt(pdfVarRms);
					histosDown[3][dist][cat][effsam]->SetBinContent(b, histos[dist][cat][effsam]->GetBinContent(b) - pdfVarRms);
					histosUp[3][dist][cat][effsam]->SetBinContent(b, histos[dist][cat][effsam]->GetBinContent(b) + pdfVarRms);
				}
			}
		}
	}
	//Split data and MC histograms for plotting and propagating uncertainties
	TH1D* dataYields[nDist][nCat];
	for(unsigned cat = 0; cat < nCat; ++cat){
		for(unsigned dist = 0; dist < nDist; ++dist){
			dataYields[dist][cat] = (TH1D*) histos[dist][cat][0]->Clone();
		}
	}
	TH1D* bkgYields[nDist][nCat][nSamples_eff -nSig]; 
	for(unsigned cat = 0; cat < nCat; ++cat){
		for(unsigned dist = 0; dist < nDist; ++dist){
			for(unsigned effsam = nSig + 1; effsam < nSamples_eff + 1; ++effsam){
				bkgYields[dist][cat][effsam -nSig - 1] = (TH1D*) histos[dist][cat][effsam]->Clone();
			}
		}
	}

	plotHistRatio(histos[2][6][nSig + 1 + 2], histosDown[0][2][6][nSig + 1 + 2], "nominal", "JEC Down", "WZtest",  false, 0,0, false, true);
	plotHistRatio(histos[2][6][nSig + 1 + 2], histosUp[0][2][6][nSig + 1 + 2], "nominal", "JEC Down", "WZtest2",  false, 0,0, false, true);
	
	const unsigned nBkg = nSamples_eff - nSig;
	const double extraUnc[nBkg] = {1.1, 1.5, 1.094, 1.15, 1.5, 1.3}; //extra flat uncertainties assigned to each background
	//Calculate histogram with systematic uncertainty for backgrounds and print the ranges of their size
	double flatSyst[2] = {0.025};
	TH1D* bkgSyst[nDist][nCat][nSamples_eff -nSig];
	for(unsigned bkg = 0; bkg < nBkg; ++bkg){
		//flatSyst[4] = fabs(pdfAccUnc[bkg + nSig + 1] - 1);
		flatSyst[1] = fabs(extraUnc[bkg] -1);
		for(unsigned cat = 0; cat < nCat; ++cat){
			for(unsigned dist = 0; dist < nDist; ++dist){
				bkgSyst[dist][cat][bkg] = (TH1D*) bkgYields[dist][cat][bkg]->Clone();
				for(unsigned b = 1;  b < dataYields[dist][cat]->GetNbinsX() + 1; ++b){
					bkgSyst[dist][cat][bkg]->SetBinContent(b,0);
					if(bkg != nBkg -1){
						double systbin = 0;
						//loop over shape uncertainties
						for(unsigned unc = 0; unc < nUnc - 1; ++unc){
							if(unc == 8) continue; //EWK subtraction unc is not defined for MC backgrounds
							if(bkg != 0 && unc == 9) continue;
							double syst = std::max(fabs(histosUp[unc][dist][cat][bkg + nSig + 1]->GetBinContent(b) - histos[dist][cat][bkg + nSig + 1]->GetBinContent(b)), fabs(histosDown[unc][dist][cat][bkg + nSig + 1]->GetBinContent(b) - histos[dist][cat][bkg + nSig + 1]->GetBinContent(b)));
							systbin += syst*syst;
						}
						//loop over flat uncertainties
						for(unsigned unc = 0; unc < 2; ++unc){
							systbin += bkgYields[dist][cat][bkg]->GetBinContent(b)*bkgYields[dist][cat][bkg]->GetBinContent(b)*flatSyst[unc]*flatSyst[unc];
						}
						bkgSyst[dist][cat][bkg]->SetBinContent(b, sqrt(systbin));
					} else{
						double syst = std::max(fabs(histosUp[8][dist][cat][bkg + nSig + 1]->GetBinContent(b) - histos[dist][cat][bkg + nSig + 1]->GetBinContent(b)), fabs(histosDown[8][dist][cat][bkg + nSig + 1]->GetBinContent(b) - histos[dist][cat][bkg + nSig + 1]->GetBinContent(b)));
						bkgSyst[dist][cat][bkg]->SetBinContent(b, sqrt(histos[dist][cat][bkg + nSig + 1]->GetBinContent(b)*0.3*histos[dist][cat][bkg + nSig + 1]->GetBinContent(b)*0.3 + syst*syst));
					}
					bkgSyst[dist][cat][bkg]->SetBinError(b, 0);
				}
			}
		}
	}

	const TString procNames[nSamples_eff + 1 - nSig] = {"obs.", "ZZ/H", "triboson", "WZ", "X + #gamma", "TT/T + X",  "non-prompt"};
	
	//Plot the yields as a function of the search region
	for(unsigned cat = 0; cat < nCat; ++cat){
		const unsigned nSigToPlot = 6;
		TH1D* signals[nSigToPlot];
		TString sigPlotNames[nSigToPlot];
		if(cat < 4){
			sigPlotNames[0] = "1 GeV, norm. to bkg. ";
			sigPlotNames[1] = "5 GeV, norm. to bkg. ";
			sigPlotNames[2] = "20 GeV, norm. to bkg. ";
			sigPlotNames[3] = "30 GeV, norm. to bkg. ";
			sigPlotNames[4] = "50 GeV, norm. to bkg. ";
			sigPlotNames[5] = "60 GeV, norm. to bkg. ";
		} else if(cat < 6){
            sigPlotNames[0] = "100 GeV, norm. to bkg. ";
            sigPlotNames[1] = "130 GeV, norm. to bkg. ";
            sigPlotNames[2] = "150 GeV, norm. to bkg. ";
            sigPlotNames[3] = "200 GeV, norm. to bkg. ";
            sigPlotNames[4] = "400 GeV, norm. to bkg. ";
            sigPlotNames[5] = "600 GeV, norm. to bkg. ";
        } else{
            sigPlotNames[0] = "5 GeV, norm. to bkg. ";
            sigPlotNames[1] = "20 GeV, norm. to bkg. ";
            sigPlotNames[2] = "40 GeV, norm. to bkg. ";
            sigPlotNames[3] = "100 GeV, norm. to bkg. ";
            sigPlotNames[4] = "200 GeV, norm. to bkg. ";
            sigPlotNames[5] = "400 GeV, norm. to bkg. ";
        }
		//loop over all kinematic distributions
		for(unsigned dist = 0; dist < nDist; ++dist){
			if(cat < 4){
				signals[0] = histos[dist][cat][2];
				signals[1] = histos[dist][cat][3];
				signals[2] = histos[dist][cat][5];
				signals[3] = histos[dist][cat][6];
				signals[4] = histos[dist][cat][8];
				signals[5] = histos[dist][cat][9];
			} else if(cat < 6){
				signals[0] = histos[dist][cat][11];
				signals[1] = histos[dist][cat][12];
				signals[2] = histos[dist][cat][13];
				signals[3] = histos[dist][cat][14];
				signals[4] = histos[dist][cat][15];
				signals[5] = histos[dist][cat][16];
			} else{
				signals[0] = histos[dist][cat][2];
				signals[1] = histos[dist][cat][5];
				signals[2] = histos[dist][cat][7];
				signals[3] = histos[dist][cat][11];
				signals[4] = histos[dist][cat][14];
				signals[5] = histos[dist][cat][15];
			}
			plotDataVSMC(dataYields[dist][cat], bkgYields[dist][cat], procNames, nSamples_eff - nSig, "searchR/" + distNames[dist] + "_" + catNames[cat] + "_withSignal" + extra, false, 0, "HNL", bkgSyst[dist][cat], true, signals,  sigPlotNames ,nSigToPlot, true);
			plotDataVSMC(dataYields[dist][cat], bkgYields[dist][cat], procNames, nSamples_eff - nSig, "searchR/" + distNames[dist] + "_" +  catNames[cat]  +  extra, false, 0, "HNL", bkgSyst[dist][cat]);
			plotDataVSMC(dataYields[dist][cat], bkgYields[dist][cat], procNames, nSamples_eff - nSig, "searchR/" + distNames[dist] + "_" +  catNames[cat]  +  extra + "_log", true, 0, "HNL", bkgSyst[dist][cat]);
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


