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
	/*
	const unsigned nSamples = 24;
	const unsigned nSamples_eff = 8;
	const TString fileList[nSamples] = {"MuonEG.root", "DoubleMuon.root", "DoubleEG.root", "SingleElectron.root", "SingleMuon.root", "TTWJetsToLNu.root", "TTZToLLNuNu.root", "TTGJets.root", "ZZTo4L.root", "TTTT.root", "VHToNonbb.root", "WWW_4F.root", "WWZ.root", "WW_DoubleScattering.root", "WZZ.root", "WpWpJJ.root", "ZZZ.root", "WZTo3LNu.root", "ZGTo2LG.root", "TTJets_DiLept.root", "TTJets_SingleLeptFromT.root", "TTJets_SingleLeptFromTbar.root", "WJetsToLNu.root", "WGToLNuG.root"};
	const double xSections[nSamples - 5] = {0.2043, 0.2529, 3.697, 1.256, 0.009103, 0.9561, 0.2086, 0.1651, 1.61704, 0.05565, 0.03711, 0.01398, 4.4297, 123.9, 87.315, 182.175, 182.175, 61526.7, 350.674}; 
	const TString names[nSamples] = {"data", "data", "data", "data", "data", "TT + X", "TT + X", "TT + X", "ZZ", "rare SM", "rare SM", "rare SM", "rare SM", "rare SM", "rare SM", "rare SM", "rare SM", "rare SM", "Z#gamma", "TT", "TT", "TT", "WJets", "W#gamma"};
	*/
	const unsigned nSamples = 31;
	const unsigned nSamples_eff = 9;
	const TString fileList[nSamples] = {"MuonEG.root", "DoubleMuon.root", "DoubleEG.root", "SingleElectron.root", "SingleMuon.root", "ttHToNonbb.root", "TTWJetsToLNu.root", "TTZToLLNuNu.root", "TTGJets.root", "ZZTo4L.root", "WWG.root", "TTTT.root", "VHToNonbb.root", "WWW.root", "WWZ.root", "WWTo2L2Nu_DoubleScattering.root", "WWTo2L2Nu.root",  "ZZZ.root",

 "WZTo3LNu.root", "ZGTo2LG.root", "TTJets_DiLept.root", "TTJets_SingleLeptFromT.root", "TTJets_SingleLeptFromTbar.root", "WJetsToLNu.root", "WGToLNuG.root", "TGJets.root", "ST_tW_antitop_NofullyHadronic.root", "ST_tW_top_NofullyHadronic.root", "ST_s-channel_leptonDecays.root", "ST_t-channel_top_inclusiveDecays.root", "ST_t-channel_antitop_inclusiveDecays.root"};
	//const double xSections[nSamples - 5] = {0.215, 0.2043, 0.2529, 3.697, 1.256, 0.2147,  0.009103, 0.9561, 0.2086, 0.1651, 0.1729,   12.178, 0.01398, 4.4297, 123.9, 87.315, 182.175, 182.175, 61526.7, 350.674, 2.967, 38.09, 38.09, 10.11, 136.02, 80.95}; 
	const double xSections[nSamples - 5] = {0.215, 0.2043, 0.2529, 3.697, 1.256, 0.2147,  0.009103, 0.9561, 0.2086, 0.1651, 0.1729,   12.178, 0.01398, 4.4297, 123.9, 87.315, 182.175, 182.175,50690, 489, 2.967, 38.09, 38.09, 10.11, 136.02, 80.95}; 
	//50690
	//350.674   489
	const TString names[nSamples] = {"data", "data", "data", "data", "data", "TT + X", "TT + X", "TT + X", "TT + X", "ZZ", "rare SM", "rare SM", "rare SM", "rare SM", "rare SM", "rare SM", "rare SM", "rare SM", "rare SM", "Z#gamma", "TT", "TT", "TT", "WJets", "W#gamma", "T + X", "T + X", "T + X", "T + X", "T + X", "T + X"};
	

	//Read Trees from ROOT files
	TFile* hfile[nSamples];
	TTree* inputTree[nSamples];
	double hcounter[nSamples];
	for(unsigned sam = 0; sam < nSamples; ++sam){
		if(fileList[sam] == "DoubleEG.root" || fileList[sam] == "SingleElectron.root"  || fileList[sam] == "DoubleMuon.root" ) continue;
		cout << "name " << names[sam] << endl;
		hfile[sam] = new TFile("../data_EWKmoriond/"+fileList[sam],"read");
       	hfile[sam]->cd("FakeElectrons");
		//Determine hcounter for cross section scaling
		TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
		_hCounter->Read("hCounter");
		hcounter[sam] = _hCounter->GetBinContent(1);
       	inputTree[sam] = (TTree*) (hfile[sam]->Get("FakeElectrons/fakeTree"));
		Init(inputTree[sam], true, sam > 5);
	}
	readSF();
	//Set up btag SF reader	
	//Btag SF
	BTagCalibration calib("csvv2", "../bTag/CSVv2_Moriond17_B_H.csv");
	BTagCalibrationReader reader(BTagEntry::OP_MEDIUM, "central", {"up", "down"});
	reader.load(calib, BTagEntry::FLAV_B, "comb");
	//photon SF
	TFile* photonSF_file =TFile::Open("../weights/photonEffSF.root");
	TH2D* photonSF = (TH2D*) photonSF_file->Get("EGamma_SF2D");
	TFile* eleVeto_file = TFile::Open("../weights/ScalingFactors_80X_Summer16.root");
	TH1D* eleVetoSF = (TH1D*) ((TH2D*) eleVeto_file->Get("Scaling_Factors_CSEV_R9 Inclusive"))->ProjectionX();
	TH1D* pixVetoSF = (TH1D*) ((TH2D*) eleVeto_file->Get("Scaling_Factors_HasPix_R9 Inclusive"))->ProjectionX();
	
	//Tweakable options////////////////////////////////////////////////////
	const bool TestRun = false;	//Break after a few events
	const double DataLuminosity = 35.867;    //units of fb^{-1}
	const bool ptReweighing = false;
	const bool pixVeto = true;
	const bool eVeto = true;
	const bool onlyMu = true;
	const double mtMin = 40;
	const double mtMax = 200;
	const double photonPtCut = 40;
	const double deltaRCut = 0.3;
	const double deltaPhiCut = 0;
	const TString extra = "_withpixveto_nomllveto_deltaR03_WjetsSubtr_cheat_wgammaCR";	//for plot file names
	//////////////////////////

	const TString eff_names[nSamples_eff + 1] = {"data", "TT + X", "ZZ", "rare SM", "Z#gamma", "TT", "WJets", "W#gamma", "T+X", "non-prompt"};
	//const TString eff_names[nSamples_eff + 1] = {"data", "TT + X", "ZZ", "rare SM", "TT", "WJets", "W#gamma", "non-prompt"};

	const unsigned nDist = 9;
	const TString distNames[nDist] = {"mt", "met", "metPhi","lepPt", "bosonPt", "DeltaPhilepmet", "DeltaRlepgamma", "DeltaPhilepgamma", "ht"};
	const TString xAxes[nDist] = {"M_{T}", "MET", "#Phi(MET)", "P_{T}(lepton)", "P_{T}(#gamma)", "#Delta#Phi(lepton, MET)", "#DeltaR(lepton, #gamma)", "#Delta#Phi(lepton, #gamma)", "H_{T}"};
	const TString units[nDist] = {"GeV", "GeV", "", "GeV", "GeV", "", "", "",  "GeV"};
	const double histMin[nDist] = {mtMin, 50, 0, 25, photonPtCut, 0, deltaRCut, deltaPhiCut, 30};
	const double histMax[nDist] = {mtMax, 200, 3.2, 150, 300, 3.2, 7, 3.2, 600};
	const int nBins[nDist] = {16, 15, 20, 15, 15, 20, 20, 20, 20};
	TH1D* hists[nDist][nSamples_eff + 1];
	for(unsigned dist = 0; dist < nDist; ++dist){
		float binWidth = (histMax[dist] - histMin[dist])/nBins[dist];
		std::ostringstream strs; strs << binWidth; std::string yAxis = strs.str();
		for(unsigned effsam = 0; effsam < nSamples_eff + 1; ++effsam){
			hists[dist][effsam] = new TH1D(distNames[dist] + eff_names[effsam], distNames[dist] + eff_names[effsam] + ";" + xAxes[dist] + " (" + units[dist] +  "); events /" + yAxis + units[dist], nBins[dist], histMin[dist], histMax[dist]);
		}
	}
	
	double maxBinC[nDist];
	for(unsigned dist = 0; dist < nDist; ++dist){
		maxBinC[dist] = hists[dist][0]->GetBinCenter(hists[dist][0]->GetNbinsX());
	}


	TH2D *MTMET[nSamples_eff + 1];
	double MTbins[4] = {40, 100,160, 220};
	double METbins[6] = {50, 100, 150, 200, 250};
	for(int i = 0; i < nSamples_eff + 1; ++i){
    	MTMET[i] = new TH2D("MTMET" + eff_names[i], "MTMET" + eff_names[i] + ";" + xAxes[0] + ";" + xAxes[1], 3, MTbins , 4, METbins);
    	MTMET[i]->Sumw2();
	}


    Double_t scale[nSamples -5];
	//set to check which data events have already been processed
	std::set<std::tuple<unsigned long, unsigned long, unsigned long> > usedEvents; //runNb, lumiBlock, eventNb
	
	//Loop over all samples
	for(unsigned sam = 0, effsam = 0; sam < nSamples; ++sam, ++effsam){
		if(sam != 0){
   			if(names[sam] == names[sam -1]) --effsam;
    	}
		if(fileList[sam] == "DoubleEG.root" || fileList[sam] == "SingleElectron.root"  || fileList[sam] == "DoubleMuon.root" ) continue;

		//if(fileList[sam] == "ttHToNonbb.root") continue;
		//if(fileList[sam] == "WW_4F.root" || fileList[sam] == "WZZ.root" || fileList[sam] == "ZZZ.root") continue;
		//if(fileList[sam] == "VHToNonbb.root" ||fileList[sam] == "WWW_4F.root" ) continue;
		cout << fileList[sam] << endl;
    	Long64_t nEntries = inputTree[sam]->GetEntries();
    	if(sam > 4){
   			scale[sam -5] = xSections[sam -5]*DataLuminosity*1000/(hcounter[sam]);
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
				scal = scale[sam-5]*_weight;
			}
			//Baseline event selection:
			if(_met < 50) continue;		//Control region with inverted MET!
			if(_nL < 1) continue;
			if(!_METfilters_pass) continue;
					
			//Clean Wjets sample
			if(names[sam] == "WJets" || fileList[sam] == "TTJets_DiLept.root" || fileList[sam] == "TTJets_SingleLeptFromTbar.root" || fileList[sam] == "TTJets_SingleLeptFromT.root" ||  fileList[sam] == "ST_s-channel_leptonDecays.root" || fileList[sam] == "ST_t-channel_top_inclusiveDecays.root" || fileList[sam] == "ST_t-channel_antitop_inclusiveDecays.root"){
				bool promptfail = false;
				for(unsigned ph = 0; ph < _gen_nPh; ++ph){
					if(_gen_lPt[ph] > 10 && fabs(_gen_phmompdg[ph] != 111)){
					//if(_gen_lPt[ph] > 10 && fabs(_gen_isPromptPh[ph])){
						promptfail = true;
						break;
					}
				}
				if(promptfail) continue;
			}
			//Calculate number of bjets 
			unsigned nJets = 0;
			unsigned nbJets = 0;
			unsigned* jetInd = new unsigned[_nJets];
			double _HT = 0;
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
						_HT += _jetPt[j];
					}
				}
			}
			if(nbJets > 0) continue;	
			//veto events with low mass mll since they aren't covered by the simulation
			/*
			bool lowM_pair = false;
			for(unsigned l1 = 0; l1 < _nL -1; ++l1){
				if(_isloose[l1]){
					for(unsigned l2 = l1 + 1; l2 < _nL; ++l2){
						if(_isloose[l2]){
							if(_flavors[l1] == _flavors[l2]){
								TLorentzVector lep1, lep2;
								lep1.SetPtEtaPhiE(_lPt[l1], _lEta[l1], _lPhi[l1], _lE[l1]);
								lep2.SetPtEtaPhiE(_lPt[l2], _lEta[l2], _lPhi[l2], _lE[l2]);
								if( (lep1 + lep2).M() < 12){
									lowM_pair = true;
									break;
								}
							}
						}
					}
				}
				if(lowM_pair) break;
			}
			if(lowM_pair) continue;
			*/
			//Check if data events were used before
			if(effsam == 0){
				auto event = usedEvents.find(std::make_tuple(_eventNb, _lumiBlock, _runNb));
				if(event != usedEvents.end()) continue;
				usedEvents.insert(std::make_tuple(_eventNb, _lumiBlock, _runNb));
			}
			unsigned* ind = new unsigned[_nL];
			unsigned lCount;
			if(!(lCount = lepOrder(ind, 1)) ) continue;
			if(lCount != 1) continue;
			//Apply Wgamma selection
			unsigned ph = 9999;
			if(!wgSel(ph, ind, lCount, photonPtCut, deltaRCut, deltaPhiCut, pixVeto, eVeto, onlyMu) ) continue;
			unsigned lw = ind[0];


			//Prompt matching
			if(effsam != 0){
				bool promptfail = false;
				for(unsigned l = 0; l < lCount; ++l){	
					//cout << _origin[ind[l]] << endl;
					if(_origin[ind[l]] != 0){
						promptfail = true;
						break;
					}
				}
				if(promptfail) continue;
			}

			//Require lepton to be tight in data and MC, and determine nonPrompt bkg in data
			bool tightFail = false;
			if(tightCount(ind, lCount) != 1) tightFail = true;
			//index used to fill events, needed to separate fakes from data
			unsigned fill = effsam;
			//Apply FR maps to data control region
			if(tightFail){ //&& effsam == 0){
				//fakes go in different histogram
				fill = nSamples_eff;
				if(effsam != 0) scal*= -1;
				//Apply FR maps
				double* conePt = new double[lCount];
				for(unsigned l = 0; l < 1; ++l){
					conePt[l] = PtCone(_lPt[ind[l]], _flavors[ind[l]], _lepMVA[ind[l]], _ptratio[ind[l]]);
				}
				scal*= fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap, lCount);
			}
			//Apply triggers to data events;
			if(effsam == 0){
				;
			}						
			//Apply ID and reco SF to simulation
			
			if(effsam != 0){
				//Apply ID and reco SF to simulation
				for(unsigned l = 0; l < lCount; ++l){  //CHANGE BACK BACK BACK
					if(_istight[ind[l]]){
						if(_flavors[ind[l]] == 2){
							scal*=0.83; //Twiki suggests flat 90% id SF for taus
						} else if(_flavors[ind[l]] == 0){
							scal*=idTightSFMap[0]->GetBinContent(idTightSFMap[0]->FindBin(TMath::Min(_lPt[ind[l]], 199.), fabs(_lEta[ind[l]])));
							//scal*=baseidSFMap[0]->GetBinContent(baseidSFMap[0]->FindBin(TMath::Min(_lPt[ind[l]], 199.), fabs(_lEta[ind[l]])));
							scal*=recSFMap_ele->GetBinContent(recSFMap_ele->FindBin(_lEta[ind[l]]));
						} else if(_flavors[ind[l]] == 1){
							scal*=idTightSFMap[1]->GetBinContent(idTightSFMap[1]->FindBin(TMath::Min(_lPt[ind[l]], 119.), fabs(_lEta[ind[l]])));
							scal*=baseidSFMap[1]->GetBinContent(baseidSFMap[1]->FindBin(TMath::Min(_lPt[ind[l]], 119.), fabs(_lEta[ind[l]])));
							//cout << recSFMap_mu_ptAbove10->GetNbinsX() << endl;
							scal*=recSFMap_mu_ptAbove10->Eval(_lEta[ind[l]]);	
						}
					} else if(_isFO[ind[l]]){
						;
					} else if(_isloose[ind[l]]){
						;
					}
				}
				scal*=photonSF->GetBinContent(photonSF->FindBin(TMath::Max(-2.5, TMath::Min(_phEta[ph], 2.5)), TMath::Min(_phPt[ph], 499.)));
				if(eVeto) scal *= eleVetoSF->GetBinContent(eleVetoSF->FindBin(TMath::Min(fabs(_phEta[ph]), 2.5) ) );
				if(pixVeto) scal *= pixVetoSF->GetBinContent(pixVetoSF->FindBin(TMath::Min(fabs(_phEta[ph]), 2.5) ) );
				//Apply btag SF
				for(unsigned j = 0; j < nJets; ++j){
					scal*=reader.eval_auto_bounds("central", BTagEntry::FLAV_B, _jetEta[jetInd[j]], _jetPt[jetInd[j]], _csv[jetInd[j]]);
				}
				//Apply PU reweighing
				scal*= PUweights->GetBinContent(PUweights->FindBin(_n_trueInteractions));
			}
			//Require triggers in data and MC
			//if(!(_IsoMu24 || _IsoTkMu24 || _Mu17_Photon22_CaloIdL_L1ISO)) continue;
			
			/////////////////////////////////
			TLorentzVector lep, metV;
			lep.SetPtEtaPhiE(PtCone(_lPt[lw], _flavors[lw], _lepMVA[lw], _ptratio[lw]), _lEta[lw], _lPhi[lw], _lE[lw]*(PtCone(_lPt[lw], _flavors[lw], _lepMVA[lw], _ptratio[lw])/_lPt[lw]) );
			metV.SetPtEtaPhiE(_met, 0, _met_phi, _met);
			if(transmass(lep,metV) < mtMin) continue;
			TLorentzVector boson;
			boson.SetPtEtaPhiE(_phPt[ph], _phEta[ph], _phPhi[ph], _phE[ph]);			
			double values[nDist] = {transmass(lep,metV),  _met, fabs(_met_phi), _lPt[lw], boson.Pt() , fabs(lep.DeltaPhi(metV)), lep.DeltaR(boson), fabs(lep.DeltaPhi(boson)), _HT};
			for(unsigned dist = 0; dist < nDist; ++dist){
				hists[dist][fill]->Fill(std::min(values[dist], maxBinC[dist]), scal);
			}
			MTMET[fill]->Fill(std::min(transmass(lep,metV) , 219.), std::min(_met, 249.), scal);
    	}
	}
	//Split data and MC histograms for plotting and propagating uncertainties
	TH1D* dataYields[nDist];
	for(unsigned dist = 0; dist < nDist; ++dist) dataYields[dist] = (TH1D*) hists[dist][0]->Clone();
	TH1D* bkgYields[nDist][nSamples_eff];
	for(unsigned effsam = 1; effsam < nSamples_eff  + 1; ++ effsam){
		for(unsigned dist = 0; dist < nDist; ++dist){
			bkgYields[dist][effsam -1] = (TH1D*)  hists[dist][effsam]->Clone();
		}
	}

	//Assign extra unc: 25% to Zgamma, 50% to other MC and 30% to non-prompt
	//Propagate systematic uncertainties to the background
	double sampUnc[nSamples_eff] = {0.15, 0.25, 0.5, 0.25, 0.5, 0 ,0, 0.5, 0.36};
	for(unsigned bkg = 0; bkg < nSamples_eff; ++bkg){
		double systUnc[1] = {sampUnc[bkg]};
		for(unsigned dist = 0; dist < nDist; ++dist){
			addSyst(bkgYields[dist][bkg], systUnc, 7);
		}
		addSyst(MTMET[bkg + 1], systUnc, 7);
	}


	
	//Calculate Wgamma Data MTvsMET histo
	for(unsigned effsam = 1; effsam < nSamples_eff + 1; ++effsam){
		if(eff_names[effsam] == "W#gamma") continue;	
		if(eff_names[effsam] == "WJets") continue;	
		MTMET[0]->Add(MTMET[effsam], -1);
	}
	MTMET[7]->Add(MTMET[6]);
	TH2D* pureSF = (TH2D*) MTMET[0]->Clone();
	pureSF->Divide(MTMET[7]);
	plotHist(pureSF, "pureSF" + extra, true);


	//Calculate Wgamma Data/WZ MC  MTvsMET histo
	TFile* WZoverWgammaMCfile = TFile::Open("../weights/WZoverWgammaMC_SF.root");

	TH2D* MTMET_WZMC = (TH2D*) WZoverWgammaMCfile->Get("WZ_MTvsMET");
	
	TH2D *MTMET_WGdata = (TH2D*) MTMET[0]->Clone();
	MTMET_WGdata->Scale( MTMET_WZMC->Integral()/MTMET_WGdata->Integral());
	MTMET_WGdata->Divide(MTMET_WZMC);
	plotHist(MTMET_WGdata, "WGammaSF_raw" + extra, true);
 	//Convolute with WZ/Wgamma MC MTvsMET histo	
	TH2D* MTMET_WZoverWgammaMC = (TH2D*) WZoverWgammaMCfile->Get("WZoverWgamma_MTvsMET");
    MTMET_WGdata->Multiply(MTMET_WZoverWgammaMC);
	plotHist(MTMET_WGdata, "WGammaSF_convoluted" + extra, true);
	
	//Plot the yields as a function of the search region
	for(unsigned dist = 0; dist < nDist; ++dist){
		plotDataVSMC(dataYields[dist], bkgYields[dist], eff_names, nSamples_eff, distNames[dist] + extra);
	}

	//Print inclusive SF
	TH1D* WZ_MT = (TH1D*) WZoverWgammaMCfile->Get("WZ_MT");
	TH1D* WZoverWgamma_MT = (TH1D*) WZoverWgammaMCfile->Get("WZoverWgamma_MT");
	//Calculate data mt distr:
	TH1D *Data_MT = (TH1D*) hists[0][0]->Clone();
	Data_MT->Sumw2();
	for(unsigned effsam = 1; effsam < nSamples_eff + 1; ++effsam){
		if(eff_names[effsam] == "W#gamma") continue;	
		Data_MT->Add(hists[0][effsam], -1);
	}

	
	const TString MTb[3] = {"$M_{T} < 100$",  "$100 < M_{T} < 160$", "$M_{T} > 160$"};
    double bins[4] = {40,100,160,220};
    TH1D* WZ_MT_reb =(TH1D*) WZ_MT->Rebin(3, "WZ_MT_reb", bins);
	WZ_MT_reb->Sumw2();
    TH1D* Data_MT_reb =(TH1D*) Data_MT->Rebin(3, "Data_MT_reb", bins);
	Data_MT_reb->Sumw2();
	Data_MT_reb->Scale(WZ_MT_reb->Integral()/ Data_MT_reb->Integral());
    TH1D *ScalFacs = HistDiv(Data_MT_reb, WZ_MT_reb);
    //plotHist(ScalFacs, "scaletest");
    for(int i = 1; i < ScalFacs->GetNbinsX() + 1; ++i){
        cout << MTb[i -1] << " : global $M_{T}$ SF = "<< std::setprecision(3) << ScalFacs->GetBinContent(i) << "$\\pm$" << std::setprecision(3) << ScalFacs->GetBinError(i) << " \\\\" << endl;
    }

	TH1D* WZoverWgamma_MT_rev = (TH1D*) WZoverWgamma_MT->Rebin(3, "WZoverWgamma_MT_reb", bins);
	//Data_MT_reb->Multiply(WZoverWgamma_MT_rev);
	//Data_MT_reb->Scale(WZ_MT_reb->Integral()/ Data_MT_reb->Integral());
	ScalFacs->Multiply(WZoverWgamma_MT_rev);
	 for(int i = 1; i < ScalFacs->GetNbinsX() + 1; ++i){
        cout << MTb[i -1] << " : global convoluted $M_{T}$ SF = "<< std::setprecision(3) << ScalFacs->GetBinContent(i) << "$\\pm$" << std::setprecision(3) << ScalFacs->GetBinError(i) << " \\\\" << endl;
    }
}



int main(int argc, char* argv[]){
	TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;
	testtree.Loop();
	rootapp->Run();
    return 0;
}


