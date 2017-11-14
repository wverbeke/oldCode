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
	const double glugluToZZkFactor = 2.1; //1.7
	const double WZSF = 0.652; //0.655
	const double ZZSF = 1.029; //1.032
	const double XgammaSF = 0.950; //0.948
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
	const TString extra = "_WZmispairing";	//for plot file names
	//////////////////////////
	TH1D* yields[nSamples_eff + 1];	//Total yields in every CR, to be used in simultaneous fit
	const TString eff_names[nSamples_eff + 1] = {"data", "ZZ/H", "triboson", "WZ", "X + #gamma", "TT/T + X",  "non-prompt"}; //X + #gamma
	const unsigned nDist = 2;  //Number of distributions to plot	
	TH1D* histos[nDist][nSamples_eff + 1];	//Kinematic distributions to plot
	const unsigned nUnc = 10;//number of shape uncertainties
	const TString uncNames[nUnc] = {"jec", "metUnclustered", "pdf", "scale", "pu", "btagSF", "triggeff", "ideff", "fakeEWK", "ZZ_mt"};  //shape uncertainty names
	//Kinematic shape uncertainty histograms
	TH1D* histosDown[nUnc][nDist][nSamples_eff + 1];
	TH1D* histosUp[nUnc][nDist][nSamples_eff + 1];
	TH1D* histosPdfVar[100][nDist][nSamples_eff + 1];
		
	//Names of the distributions to plot
	const TString histNames[nDist] = {"Mll", "mt"};
  	//X-axis labels of distributions to plot
	const TString xAxes[nDist] = {"M_{ll} (GeV)", "M_{T} (GeV)"};
	//Units of distributions to plot.
	const TString units[nDist] = {"GeV", "GeV"};
	//Minimum x-range of histograms.
	const double histMin[nDist] = { 0, 0}; 
	//Minimum y-range of histograms.
	const double histMax[nDist] = {300, 300};
	
	unsigned nBins[nDist];
	for(unsigned dist = 0; dist < nDist; ++dist) nBins[dist] = 20;
	
	//Initialize histograms of kinematic distributions
	for(unsigned dist = 0; dist < nDist; ++dist){
		float binWidth = (histMax[dist] - histMin[dist])/nBins[dist];
		std::ostringstream strs; strs << binWidth; std::string yAxis = strs.str();
		for(unsigned effsam = 0; effsam < nSamples_eff + 1; ++effsam){
			histos[dist][effsam] = new TH1D(eff_names[effsam] + histNames[dist], eff_names[effsam] + histNames[dist] + ";" + xAxes[dist] + "; Events /" + yAxis + units[dist], nBins[dist], histMin[dist], histMax[dist]);
			histos[dist][effsam]->Sumw2();
			for(unsigned unc = 0; unc < nUnc; ++unc){
				histosDown[unc][dist][effsam] = new TH1D(eff_names[effsam] + histNames[dist] + uncNames[unc] + "Down", eff_names[effsam] + histNames[dist]  + uncNames[unc] + "Down;" + xAxes[dist] + "; Events /" + yAxis + units[dist], nBins[dist], histMin[dist], histMax[dist]);
				histosUp[unc][dist][effsam] = new TH1D(eff_names[effsam] + histNames[dist] + uncNames[unc] + "Up", eff_names[effsam] + histNames[dist]  + uncNames[unc] + "Up;" + xAxes[dist] + "; Events /" + yAxis + units[dist], nBins[dist], histMin[dist], histMax[dist]);
			}
			for(unsigned pdf = 0; pdf < 100; ++pdf){
				histosPdfVar[pdf][dist][effsam] = new TH1D(eff_names[effsam] + histNames[dist] + "_pdf" + std::to_string(pdf), eff_names[effsam] + histNames[dist] + "_pdf" + std::to_string(pdf) + ";" + xAxes[dist] + "; Events /" + yAxis + units[dist], nBins[dist], histMin[dist], histMax[dist]);
			}
		}
	}
	double maxBinC[nDist];
	for(unsigned dist = 0; dist < nDist; ++dist){
		maxBinC[dist] = histos[dist][0]->GetBinCenter(histos[dist][0]->GetNbinsX());
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
			if(effsam > 0 && nBJets(true, false, 2) != 0) continue;
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
			if(lCount != 3) continue;
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
			if(_flavors[mllI[0]] != _flavors[mllI[1]]) continue;
			if( fabs(mll - 91.) > 15 && conePt[0] > 55) continue;
			if( fabs(lepSyst.M() - 91.) < 15) continue;
			if(conePt[0] < 25.) continue;
			if(conePt[1] < 15.) continue;
			if(conePt[2] < 10.) continue;

			//determine the index of the W lepton
			unsigned lw = 9999;
			for(unsigned l = 0; l < lCount; ++l){
				if(ind[l] != mllI[0] && ind[l] != mllI[1]){
					lw = ind[l];
				}
			}

			
			//MISPAIRING
			if(_flavors[mllI[0]] == _flavors[lw]) continue;
			unsigned mispI[2], mtI;
			if(_charges[mllI[0]] != _charges[lw]){
				mispI[0] = mllI[0];
				mtI = mllI[1];
			} else if( _charges[mllI[1]] != _charges[lw]){
				mispI[0] = mllI[1];
				mtI = mllI[0];
			}
			mispI[1] = lw;
			mllI[0] = mispI[0];
			mllI[1] = mispI[1];
			lw = mtI;


			TLorentzVector METvec;
			METvec.SetPtEtaPhiE(_met, 0, _met_phi, _met);
			//NEEDED FOR ZZ NUISANCE 
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
			double mt_min = transmass(lepV[lw_min], METvec);
			///////////////////////////////////////////////
			
			//new mll
			/*
			TLorentzVector lzV[2];
			for(unsigned l = 0; l < 2; ++l) lzV[l].SetPtEtaPhiE(PtCone(_lPt[mllI[l]], _flavors[mllI[l]], _lepMVA[mllI[l]], _ptratio[mllI[l]]), _lEta[mllI[l]], _lPhi[mllI[l]], _lE[mllI[l]]);
			mll = (lzV[0] + lzV[1]).M();
			*/
			TLorentzVector Wlep;
			Wlep.SetPtEtaPhiE(_lPt[lw]*std::max(1., 1 + (_isolation[lw] - 0.1)), _lEta[lw], _lPhi[lw], _lE[lw]*std::max(1., 1 + (_isolation[lw] - 0.1)));
			//Scale weights for glugluToZZ
			if(fileList[sam].Contains("GluGluToZZTo")){
				for(unsigned pdf = 0; pdf < 110; ++pdf){
					_scaleWeight[pdf] = 1;
				}
			}

			for(unsigned i = 0; i < 2; ++i){
				if(i == 1) continue;	
				///////////////////////////////////////////////////
				//if(i == 0) Wlep.SetPtEtaPhiE(_lPt[mllI[0]]*std::max(1., 1 + (_isolation[mllI[0]] - 0.1)), _lEta[mllI[0]], _lPhi[mllI[0]], _lE[mllI[0]]*std::max(1., 1 + (_isolation[mllI[0]] - 0.1)));
				//else Wlep.SetPtEtaPhiE(_lPt[mllI[1]]*std::max(1., 1 + (_isolation[mllI[1]] - 0.1)), _lEta[mllI[1]], _lPhi[mllI[1]], _lE[mllI[1]]*std::max(1., 1 + (_isolation[mllI[1]] - 0.1)));

				double values[nDist] = { mll, transmass(Wlep, METvec)}; 
				//Nominal yields
				if(_met > 50 && nBJets(true, true, 0) == 0){
					for(unsigned dist = 0; dist < nDist; ++dist){
						histos[dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal);
					}
				}

				//FR EWK contamination uncertainty
				if(fill == nSamples_eff && nBJets(true, false, 0) == 0 && _met > 50){
					for(unsigned dist = 0; dist < nDist; ++dist){			
						histosDown[8][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal*(fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[1], lCount)/fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[0], lCount) ) );
						histosUp[8][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal*(fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[2], lCount)/fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[0], lCount) ) );
					}
				}
			
				if(effsam == 0) continue;
				//vary JEC down
				METvec.SetPtEtaPhiE(_metJECDown, 0, _met_phiJECDown, _metJECDown);
				values[1] = transmass(Wlep, METvec);
				if(_metJECDown > 50 && nBJets(true, true, 1) == 0){
					for(unsigned dist = 0; dist < nDist; ++dist){
						histosDown[0][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal);
					}
				}
				//vary JEC up
				METvec.SetPtEtaPhiE(_metJECUp, 0, _met_phiJECUp, _metJECUp);
				values[1] = transmass(Wlep, METvec);
				if(_metJECUp > 50 && nBJets(true, true, 2) == 0){
					for(unsigned dist = 0; dist < nDist; ++dist){
						histosUp[0][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal);
					}
				}
				//nominal b-veto 
				if(nBJets(true, false, 0) != 0) continue;
				//nominal HT and nJets values
				//vary unclustered met down
				METvec.SetPtEtaPhiE(_metOtherDown, 0, _met_phiOtherDown, _metOtherDown);
				values[1] = transmass(Wlep, METvec);
				if(_metOtherDown > 50){
					for(unsigned dist = 0; dist < nDist; ++dist){
						histosDown[1][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal);
					}
				}
			
				//vary unclustered met up
				METvec.SetPtEtaPhiE(_metOtherUp, 0, _met_phiOtherUp, _metOtherUp);
				values[1] = transmass(Wlep, METvec);
				if(_metOtherUp > 50){
					for(unsigned dist = 0; dist < nDist; ++dist){
						histosUp[1][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal);
					}
				}
				//Nominal met cut 
				if(_met < 50) continue;
				METvec.SetPtEtaPhiE(_met, 0, _met_phi, _met);
				//Nominal values to fill
				values[2] = transmass(Wlep, METvec);
				//vary scale down
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosDown[3][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal*_scaleWeight[8]);
				}
				//vary scale up
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosUp[3][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal*_scaleWeight[4]);
				}
				//vary pu down
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosDown[4][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), (scal/PUweights[0]->GetBinContent(PUweights[0]->FindBin( std::min(_n_trueInteractions, 49) )))*PUweights[1]->GetBinContent(PUweights[1]->FindBin(std::min(_n_trueInteractions, 49)  )) );
				}
				//vary pu up
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosUp[4][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), (scal/PUweights[0]->GetBinContent(PUweights[0]->FindBin( std::min(_n_trueInteractions, 49) )))*PUweights[2]->GetBinContent(PUweights[2]->FindBin(std::min(_n_trueInteractions, 49)  )) );
				}
				//vary btag SF down	
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosDown[5][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]),  (scal/bTagSF(true, 0))*bTagSF(true, 1) );
				}
				//vary btag SF up
				for(unsigned dist = 0; dist < nDist; ++dist){
					histosUp[5][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]),  (scal/bTagSF(true, 0))*bTagSF(true, 2) );
				}
				//all pdf variations
				for(unsigned pdf = 0; pdf < 100; ++pdf){
					for(unsigned dist = 0; dist < nDist; ++dist){
						histosPdfVar[pdf][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal*_scaleWeight[pdf + 9]);
					}
				}
				//Trigger uncertainties
				if(conePt[0] > 30){
					for(unsigned dist = 0; dist < nDist; ++dist){
						histosDown[6][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal*0.95);
						histosUp[6][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal*1.05);
					}
				} else{
					for(unsigned dist = 0; dist < nDist; ++dist){
						histosDown[6][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal*0.98);
						histosUp[6][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal*1.02);
					}
				}
				//Id efficiency uncertainties
				unsigned flavorC = tril_flavorComb(ind, _flavors, lCount);
				if(flavorC == 0 || flavorC == 3){
					for(unsigned dist = 0; dist < nDist; ++dist){
						histosDown[7][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal*0.94);
						histosUp[7][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal*1.06);
					}
				} else{
					for(unsigned dist = 0; dist < nDist; ++dist){
						histosDown[7][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal*0.9553);
						histosUp[7][dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal*1.0447);
					}
				}	
				if(effsam - 1 == 0){// only ZZ needs last nuisance!
					if(mt_min < 75){
						for(unsigned dist = 0; dist < nDist; ++dist){
							histosDown[9][dist][fill]->Fill(std::min(values[dist], maxBinC[dist]), scal);
							histosUp[9][dist][fill]->Fill(std::min(values[dist], maxBinC[dist]), scal);
						}
					} else{
						for(unsigned dist = 0; dist < nDist; ++dist){
							histosDown[9][dist][fill]->Fill(std::min(values[dist], maxBinC[dist]), scal*0.75);
							histosUp[9][dist][fill]->Fill(std::min(values[dist], maxBinC[dist]), scal*1.25);
						}
					}					
				}
			}
    	}
		//Set negative bins to 0 before adding other processes
		for(unsigned dist = 0; dist < nDist; ++dist){
			for(unsigned bin = 1; bin < histos[dist][effsam]->GetNbinsX() + 1; ++bin){
				if(histos[dist][effsam]->GetBinContent(bin) < 0 ) histos[dist][effsam]->SetBinContent(bin, 0.);
				for(unsigned unc = 0; unc < nUnc; ++unc){
					if(histosUp[unc][dist][effsam]->GetBinContent(bin) < 0. ) histosUp[unc][dist][effsam]->SetBinContent(bin, 0.);
					if(histosDown[unc][dist][effsam]->GetBinContent(bin) < 0. ) histosDown[unc][dist][effsam]->SetBinContent(bin, 0.);
				}	
				for(unsigned pdf = 0; pdf < 100; ++pdf){
					if(histosPdfVar[pdf][dist][effsam]->GetBinContent(bin) < 0. ) histosPdfVar[pdf][dist][effsam]->SetBinContent(bin, 0.);
				}
			}
		}
	}
	//Calculate rms of pdf shape variations
	//for kinematic distributions
	for(unsigned effsam = 1; effsam < nSamples_eff; ++effsam){
		for(unsigned dist = 0; dist < nDist; ++dist){
			for(unsigned b = 1; b < histos[dist][effsam]->GetNbinsX() + 1; ++b){
				double pdfVarRms = 0;
				for(unsigned pdf = 0; pdf < 100; ++pdf){
					pdfVarRms += (histos[dist][effsam]->GetBinContent(b) - histosPdfVar[pdf][dist][effsam]->GetBinContent(b))*(histos[dist][effsam]->GetBinContent(b) - histosPdfVar[pdf][dist][effsam]->GetBinContent(b));
				}
				pdfVarRms = 0.01*sqrt(pdfVarRms);
				histosDown[2][dist][effsam]->SetBinContent(b, histos[dist][effsam]->GetBinContent(b) - pdfVarRms);
				histosUp[2][dist][effsam]->SetBinContent(b, histos[dist][effsam]->GetBinContent(b) + pdfVarRms);
			}
		}
	}
	//Split data and MC histograms for plotting and propagating uncertainties
	TH1D* dataHistos[nDist];
	for(unsigned dist = 0; dist < nDist; ++dist){
		dataHistos[dist] = (TH1D*) histos[dist][0]->Clone();
	}
	
	TH1D* bkgHistos[nDist][nSamples_eff]; //change to nSamples_eff if sig is removed	
	for(unsigned dist = 0; dist <nDist; ++dist){
		for(unsigned effsam = 1; effsam < nSamples_eff + 1; ++effsam){
			bkgHistos[dist][effsam -1] = (TH1D*) histos[dist][effsam]->Clone();
		}
	}
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//Calculate histogram with systematic uncertainty for backgrounds and print the ranges of their size
	//const double extraUnc[nSamples_eff] = {1.096, 1.5, 1.085, 1.079, 1.15, 1.3};
	const double extraUnc[nSamples_eff] = {1.1, 1.5, 1.094, 1.15, 1.15, 1.3}; //1.25
	double flatSyst[2] = {0.025};
	TH1D* bkgSystDist[nDist][nSamples_eff];
	for(unsigned bkg = 0; bkg < nSamples_eff; ++bkg){
		flatSyst[1] = extraUnc[bkg] - 1.;
		for(unsigned dist = 0; dist < nDist; ++dist){
			bkgSystDist[dist][bkg] = (TH1D*) bkgHistos[dist][bkg]->Clone();
			for(unsigned b = 1; b < bkgHistos[dist][bkg]->GetNbinsX() + 1; ++b){
				bkgSystDist[dist][bkg]->SetBinContent(b,0);
				/*
				if(histos[dist][bkg + 1]->GetBinContent(b) == 0){
					bkgSystDist[dist][bkg]->SetBinContent(b,0);
				} else{
				*/
					if(bkg != nSamples_eff -1){
						double systbin = 0;
						//loop over shape uncertainties
						for(unsigned unc = 0; unc < nUnc; ++unc){
							if(unc == 8) continue;
							if(bkg != 0 && unc == 9) continue;
							double syst = std::max(fabs(histosUp[unc][dist][bkg + 1]->GetBinContent(b) - histos[dist][bkg + 1]->GetBinContent(b)), fabs(histosDown[unc][dist][bkg + 1]->GetBinContent(b) - histos[dist][bkg + 1]->GetBinContent(b)));
							systbin += syst*syst;
						}
						//loop over flat uncertainties
						for(unsigned unc = 0; unc < 2; ++unc){
							systbin += histos[dist][bkg + 1]->GetBinContent(b)*histos[dist][bkg + 1]->GetBinContent(b)*flatSyst[unc]*flatSyst[unc];
						}
						bkgSystDist[dist][bkg]->SetBinContent(b, sqrt(systbin));
					} else{
						double syst = std::max(fabs(histosUp[8][dist][bkg + 1]->GetBinContent(b) - histos[dist][bkg + 1]->GetBinContent(b)), fabs(histosDown[8][dist][bkg + 1]->GetBinContent(b) - histos[dist][bkg + 1]->GetBinContent(b)));
						bkgSystDist[dist][bkg]->SetBinContent(b, sqrt(histos[dist][bkg + 1]->GetBinContent(b)*0.3*histos[dist][bkg + 1]->GetBinContent(b)*0.3 + syst*syst) );
					}
					bkgSystDist[dist][bkg]->SetBinError(b, 0);
				//}
			}
		}
	}

	for(unsigned dist = 0; dist < nDist; ++dist){
		plotDataVSMC(dataHistos[dist], bkgHistos[dist], eff_names, nSamples_eff, "controlR/" +  histNames[dist] + extra, false, 0, "HNL", bkgSystDist[dist]);
		plotDataVSMC(dataHistos[dist], bkgHistos[dist], eff_names, nSamples_eff, "controlR/" +  histNames[dist] + extra + "_log", true, 0, "HNL", bkgSystDist[dist]);
	}
}



int main(int argc, char* argv[]){
	TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;
	testtree.Loop();
	rootapp->Run();
    return 0;
}


