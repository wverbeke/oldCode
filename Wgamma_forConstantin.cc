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

	const unsigned nSamples = 20; //28
	const unsigned nSamples_eff = 8;
	const TString fileList[nSamples] = {"MuonEG.root", "SingleMuon.root", "ttHToNonbb.root", "TTWJetsToLNu.root", "TTZToLLNuNu.root", "TTGJets.root","WWG.root", "TTTT.root", "WWW.root", "WWZ.root", "WWTo2L2Nu_DoubleScattering.root", "WWTo2L2Nu.root",  "ZZZ.root",
 "ZGTo2LG.root", "TTJets_DiLept.root", "TTJets_SingleLeptFromT.root", "TTJets_SingleLeptFromTbar.root", "WJetsToLNu.root", "WGToLNuG.root", "TGJets.root"};//, "ST_tW_antitop_NofullyHadronic.root", "ST_tW_top_NofullyHadronic.root", "ST_s-channel_leptonDecays.root", "ST_t-channel_top_inclusiveDecays.root", "ST_t-channel_antitop_inclusiveDecays.root"};
	const double xSections[nSamples - 2] = {0.215, 0.2043, 0.2529, 3.697, 0.2147,  0.009103, 0.2086, 0.1651, 0.1729,   12.178, 0.01398, 123.9, 87.315, 182.175, 182.175,50690, 489, 2.967};//, 38.09, 38.09, 10.11, 136.02, 80.95}; 
	const TString names[nSamples] = {"data", "data", "TT + X", "TT + X", "TT + X", "TT + X", "rare SM", "rare SM", "rare SM", "rare SM", "rare SM", "rare SM", "rare SM", "Z#gamma", "TT", "TT", "TT", "WJets", "W#gamma", "T + X"};//, "T + X", "T + X", "T + X", "T + X", "T + X"};

	
	//Read Trees from ROOT files
	TFile* hfile[nSamples];
	TTree* inputTree[nSamples];
	double hcounter[nSamples];
	for(unsigned sam = 0; sam < nSamples; ++sam){
		std::cout << "name " << names[sam] << ": \t " << fileList[sam] << std::endl;
		hfile[sam] = new TFile("../data_newWgamma/"+fileList[sam],"read");
       	hfile[sam]->cd("FakeElectrons");
		//Determine hcounter for cross section scaling
		TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
		_hCounter->Read("hCounter");
		hcounter[sam] = _hCounter->GetBinContent(1);
       	inputTree[sam] = (TTree*) (hfile[sam]->Get("FakeElectrons/fakeTree"));
		Init(inputTree[sam], true, sam > 0);
	}
	readSF(true);
	TFile* photonSF_file =TFile::Open("../weights/photonEffSF.root");
	TH2D* photonSF = (TH2D*) photonSF_file->Get("EGamma_SF2D");
	TFile* eleVeto_file = TFile::Open("../weights/ScalingFactors_80X_Summer16.root");
	TH1D* eleVetoSF = (TH1D*) ((TH2D*) eleVeto_file->Get("Scaling_Factors_CSEV_R9 Inclusive"))->ProjectionX();
	TH1D* pixVetoSF = (TH1D*) ((TH2D*) eleVeto_file->Get("Scaling_Factors_HasPix_R9 Inclusive"))->ProjectionX();
	//Tweakable options////////////////////////////////////////////////////
	const bool TestRun = false;	//Break after a few events
	const double DataLuminosity = 35.867;    //units of fb^{-1}
	const bool pixVeto = true;
	const bool eVeto = true;
	const bool onlyMu = true;
	const double mtMin = 40;
	const double mtMax = 200;
	const double photonPtCut = 40;
	const double deltaRCut = 0.3;
	const double deltaPhiCut = 0;
	const TString extra = "_withpixveto";	//for plot file names
	//////////////////////////
	
	const TString eff_names[nSamples_eff + 1] = {"data", "TT + X", "rare SM", "Z#gamma", "TT", "WJets", "W#gamma", "T+X", "non-prompt"};
	//Define histograms!
	const unsigned nDist = 9;
	const TString distNames[nDist] = {"mt", "met", "metPhi","lepPt", "bosonPt", "DeltaPhilepmet", "DeltaRlepgamma", "DeltaPhilepgamma", "ht"};
	const TString xAxes[nDist] = {"M_{T}", "MET", "#Phi(MET)", "P_{T}(lepton)", "P_{T}(#gamma)", "#Delta#Phi(lepton, MET)", "#DeltaR(lepton, #gamma)", "#Delta#Phi(lepton, #gamma)", "H_{T}"};
	const TString units[nDist] = {"GeV", "GeV", "", "GeV", "GeV", "", "", "",  "GeV"};
	const double histMin[nDist] = {mtMin, 50, 0, 25, photonPtCut, 0, deltaRCut, deltaPhiCut, 30};
	const double histMax[nDist] = {mtMax, 200, 3.2, 150, 300, 3.2, 7, 3.2, 600};
	const int nBins[nDist] = {16, 15, 20, 15, 15, 20, 20, 20, 20};
	const unsigned nUnc = 7;//number of shape uncertainties
	const TString uncNames[nUnc] = {"jec", "metUnclustered",  "scale", "pdf", "pu", "btagSF", "fakeEWK"};

	TH1D* hists[nDist][nSamples_eff + 1]; //Central values
	TH1D* histsDown[nUnc][nDist][nSamples_eff + 1]; //Downwards varied yields
	TH1D* histsUp[nUnc][nDist][nSamples_eff + 1]; //Upwards varies yields
	TH1D* histsPdfVar[100][nDist][nSamples_eff + 1];
	
	for(unsigned dist = 0; dist < nDist; ++dist){
		float binWidth = (histMax[dist] - histMin[dist])/nBins[dist];
		std::ostringstream strs; strs << binWidth; std::string yAxis = strs.str();
		for(unsigned effsam = 0; effsam < nSamples_eff + 1; ++effsam){
			hists[dist][effsam] = new TH1D(distNames[dist] + eff_names[effsam], distNames[dist] + eff_names[effsam] + ";" + xAxes[dist] + " (" + units[dist] +  "); events /" + yAxis + units[dist], nBins[dist], histMin[dist], histMax[dist]);
			for(unsigned unc = 0; unc < nUnc; ++unc){
				histsDown[unc][dist][effsam] = new TH1D(distNames[dist] + eff_names[effsam] + uncNames[unc] + "Down", distNames[dist] + eff_names[effsam]+ uncNames[unc] + "Down;" + xAxes[dist] + " (" + units[dist] +  "); events /" + yAxis + units[dist], nBins[dist], histMin[dist], histMax[dist]);
				histsUp[unc][dist][effsam] = new TH1D(distNames[dist] + eff_names[effsam] + uncNames[unc] + "Up", distNames[dist] + eff_names[effsam]+ uncNames[unc] + + "Up;" + xAxes[dist] + " (" + units[dist] +  "); events /" + yAxis + units[dist], nBins[dist], histMin[dist], histMax[dist]);
			}
			for(unsigned pdf = 0; pdf < 100; ++pdf){
				histsPdfVar[pdf][dist][effsam] = new TH1D(distNames[dist] + eff_names[effsam] + "_pdf" + std::to_string(pdf), distNames[dist] + eff_names[effsam] + "_pdf" + std::to_string(pdf) + ";" + xAxes[dist] + " (" + units[dist] +  "); events /" + yAxis + units[dist], nBins[dist], histMin[dist], histMax[dist]);
			}
		}
	}

	double maxBinC[nDist];
	for(unsigned dist = 0; dist < nDist; ++dist){
		maxBinC[dist] = hists[dist][0]->GetBinCenter(hists[dist][0]->GetNbinsX());
	}

	Double_t scale[nSamples - 2];

	//set to check which data events have already been processed
	std::set<std::tuple<unsigned long, unsigned long, unsigned long> > usedEvents; //runNb, lumiBlock, eventNb	
	//Loop over all samples
	for(unsigned sam = 0, effsam = 0; sam < nSamples; ++sam, ++effsam){
		if(sam != 0){
   			if(names[sam] == names[sam -1]) --effsam;
    	}		
    	Long64_t nEntries = inputTree[sam]->GetEntries();
    	if(sam > 1){
   			scale[sam - 2] = xSections[sam - 2]*DataLuminosity*1000/(hcounter[sam]);
    	}
		cout << eff_names[effsam] << endl;
    	std::cout<<"Entries in "<< fileList[sam] <<" "<<nEntries<<std::endl;
		cout << effsam << endl;
		
        for(Long64_t it = 0; it < nEntries; ++it){
        	inputTree[sam]->GetEntry(it);
        	if (it%10000 == 0) std::cout<<'.'<<std::flush;
        	if(TestRun && it > 10000) break;
        	double scal;
        	if(effsam == 0) scal = 1;
        	else{
				scal = scale[sam- 2]*_weight;
			}		
			cutBased();
			//Baseline event selection
			if(_nL < 1) continue;
			if(!_METfilters_pass) continue;
			//Preliminary b-veto
			cleanJets();
			if(effsam > 0 && nBJets(false, true, 1) != 0) continue;
			else if(effsam == 0 && nBJets(false, true, 0) != 0) continue;
			//Check if data events were used before
			if(effsam == 0){
				auto event = usedEvents.find(std::make_tuple(_eventNb, _lumiBlock, _runNb));
				if(event != usedEvents.end()) continue;
				usedEvents.insert(std::make_tuple(_eventNb, _lumiBlock, _runNb));
			}
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

			unsigned ind[(const unsigned) _nL]; // = new unsigned[_nL];
			const unsigned lCount = lepOrder(ind, 1);
			if(lCount != 1) continue;
			//Apply Wgamma selection
			unsigned ph = 9999;
			if(!wgSel(ph, ind, lCount, photonPtCut, deltaRCut, deltaPhiCut, pixVeto, eVeto, onlyMu) ) continue;
			unsigned lw = ind[0];


			//MC prompt matching
			if(effsam  > 0	){
				bool promptfail = false;
				for(unsigned l = 0; l < lCount; ++l){	
					if(_origin[ind[l]] != 0){
						promptfail = true;
						break;
					}
				}
				if(promptfail) continue;
			}
			//Require 1 leptons to be tight in data and MC, and determine nonPrompt bkg in data
			unsigned nTight = tightCount(ind, lCount);
			bool tightFail = nTight < 1;
			//index used to fill events, needed to separate fakes from data
			unsigned fill = effsam;
			//Apply FR maps to data control region
			double* conePt = new double[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				//EWKino cone pT
				conePt[l] = 	PtCone(_lPt[ind[l]], _flavors[ind[l]], _lepMVA[ind[l]], _ptratio[ind[l]]);
			}
			if(tightFail){ //&& effsam == 0){
				//fakes go in different histogram
				fill = nSamples_eff;
				//Apply FR maps
				if(effsam != 0) scal *= -1;
				scal*= fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[0], lCount);
			} else if(tightFail) continue;

			//Apply triggers to data events;
			if(effsam == 0){
				;
			}	
	
			//Apply ID and reco SF to simulation
			if(effsam != 0){
				scal*=getEventSF(ind, lCount, false);
				scal*=photonSF->GetBinContent(photonSF->FindBin(TMath::Max(-2.5, TMath::Min(_phEta[ph], 2.5)), TMath::Min(_phPt[ph], 499.)));
				if(eVeto) scal *= eleVetoSF->GetBinContent(eleVetoSF->FindBin(TMath::Min(fabs(_phEta[ph]), 2.5) ) );
				if(pixVeto) scal *= pixVetoSF->GetBinContent(pixVetoSF->FindBin(TMath::Min(fabs(_phEta[ph]), 2.5) ) );
			}					
			

			//Set TLorentzVector for each lepton
			TLorentzVector lepV[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				lepV[l].SetPtEtaPhiE(conePt[l], _lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]*(conePt[l]/_lPt[ind[l]]) );
			}
			//Calculate MET vector
			TLorentzVector METvec;
			METvec.SetPtEtaPhiE(_met, 0, _met_phi, _met);
			
			//photon lorentz vector
			TLorentzVector photVec;
			photVec.SetPtEtaPhiE(_phPt[ph], _phEta[ph], _phPhi[ph], _phE[ph]);
			//compute kinematic distributions
			//mt, met, metPhi, lepPt, bosonPt, DeltaPhilepmet, DeltaRlepgamma, DeltaPhilepgamma, ht
			double values[nDist] = {transmass(lepV[0], METvec), _met, _met_phi, conePt[0], _phPt[ph], lepV[0].DeltaPhi(METvec), photVec.DeltaR(lepV[0]), photVec.DeltaPhi(lepV[0]), _HT};

			//Fill nominal yields
			if(nBJets(false, true, 0) == 0 && _met > 50){
				for(unsigned dist = 0; dist < nDist; ++dist){
					hists[dist][fill]->Fill(std::min(values[dist], maxBinC[dist]), scal);
				}
			}
			
			//FR EWK contamination uncertainty
			if(fill == nSamples_eff && nBJets(false, true, 0) == 0 && _met > 50){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histsDown[6][dist][fill]->Fill(std::min(values[dist], maxBinC[dist]), scal*(fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[1], lCount)/fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[0], lCount) ));
					histsUp[6][dist][fill]->Fill(std::min(values[dist], maxBinC[dist]), scal*(fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[2], lCount)/fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap[0], lCount) ));

				}
			}		
			//Following uncertainties do not apply to data or data-driven backgrounds
			if(effsam != 0) continue;
			//yields with JEC varied down
			METvec.SetPtEtaPhiE(_metJECDown, 0, _met_phiJECDown, _metJECDown);	
			values[0] = transmass(lepV[0], METvec);
			values[1] = _metJECDown;
			values[2] = _met_phiJECDown;
			values[5] = lepV[0].DeltaPhi(METvec);	
			//Recompute HT 
			_HT = 0;
			for(unsigned j = 0; j < _nJets; ++j){
				if(_jetPtDown[j] > 30){
					_HT += _jetPtDown[j];
				}		
			}	
			values[8] = _HT;	
			if(nBJets(false, true, 1) == 0 && _metJECDown > 50){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histsDown[0][dist][fill]->Fill(std::min(values[dist], maxBinC[dist]), scal);
				}
			}
			//yields with JEC varied up
			METvec.SetPtEtaPhiE(_metJECUp, 0, _met_phiJECUp,_metJECUp);	
			values[0] = transmass(lepV[0], METvec);
			values[1] = _metJECUp;
			values[2] = _met_phiJECUp;
			values[5] = lepV[0].DeltaPhi(METvec);	
			//Recompute HT 
			_HT = 0;
			for(unsigned j = 0; j < _nJets; ++j){
				if(_jetPtUp[j] > 30){
					_HT += _jetPtUp[j];
				}		
			}	
			values[8] = _HT;	
			if(nBJets(false, true, 2) == 0 && _metJECUp > 50){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histsUp[0][dist][fill]->Fill(std::min(values[dist], maxBinC[dist]), scal);
				}
			}
			//nominal b-veto
			if(nBJets(false, true, 0) != 0) continue;	
			//recompute HT
			_HT = 0;
			for(unsigned j = 0; j < _nJets; ++j){
				if(_jetPtUp[j] > 30){
					_HT += _jetPt[j];
				}		
			}	
			values[8] = _HT;	
			//yields with unclustered met varied down
			METvec.SetPtEtaPhiE(_metOtherDown, 0, _met_phiOtherDown, _metOtherDown);			
			values[0] = transmass(lepV[0], METvec);
			values[1] = _metJECUp;
			values[2] = _met_phiJECUp;
			values[5] = lepV[0].DeltaPhi(METvec);	
			if(_metOtherDown > 50){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histsDown[1][dist][fill]->Fill(std::min(values[dist], maxBinC[dist]), scal);
				}
			}
			//yields with unclustered met varied up
			METvec.SetPtEtaPhiE(_metOtherUp, 0, _met_phiOtherUp, _metOtherUp);	
			values[0] = transmass(lepV[0], METvec);
			values[1] = _metJECUp;
			values[2] = _met_phiJECUp;
			values[5] = lepV[0].DeltaPhi(METvec);
			if(_metOtherUp > 50){
				for(unsigned dist = 0; dist < nDist; ++dist){
					histsUp[1][dist][fill]->Fill(std::min(values[dist], maxBinC[dist]), scal);
				}
			}
			//nominal met cut for low-mass categories
			if(_met < 50) continue;
			//nominal kinematics
			values[0] = transmass(lepV[0], METvec);
			values[1] = _metJECUp;
			values[2] = _met_phiJECUp;
			values[5] = lepV[0].DeltaPhi(METvec);	

			//yields with scale varied down
			for(unsigned dist = 0; dist < nDist; ++dist){
				histsDown[2][dist][fill]->Fill(std::min(values[dist], maxBinC[dist]), scal*_scaleWeight[8]);
			}
			//yields with scale varied up
			for(unsigned dist = 0; dist < nDist; ++dist){
				histsUp[2][dist][fill]->Fill(std::min(values[dist], maxBinC[dist]), scal*_scaleWeight[4]);
			}
			//yields with all pdf variations
			for(unsigned dist = 0; dist < nDist; ++dist){
				for(unsigned pdf = 0; pdf < 100; ++pdf){
					histsPdfVar[pdf][dist][fill]->Fill(std::min(values[dist], maxBinC[dist]), scal*_scaleWeight[pdf + 9]);
				}
			}
			//yields with pileup varied down
			for(unsigned dist = 0; dist < nDist; ++dist){
				histsDown[4][dist][fill]->Fill(std::min(values[dist], maxBinC[dist]), (scal/PUweights[0]->GetBinContent(PUweights[0]->FindBin( std::min(_n_trueInteractions, 49) )))*PUweights[1]->GetBinContent(PUweights[1]->FindBin(std::min(_n_trueInteractions, 49)  )) );
			}
			//yields with pileup varied up
			for(unsigned dist = 0; dist < nDist; ++dist){
				histsUp[4][dist][fill]->Fill(std::min(values[dist], maxBinC[dist]), (scal/PUweights[0]->GetBinContent(PUweights[0]->FindBin( std::min(_n_trueInteractions, 49) )))*PUweights[2]->GetBinContent(PUweights[2]->FindBin(std::min(_n_trueInteractions, 49) )) );
			}
			//yields with b-tag SF varied down
			for(unsigned dist = 0; dist < nDist; ++dist){
				histsDown[5][dist][fill]->Fill(std::min(values[dist], maxBinC[dist]), (scal/bTagSF(false, 0))*bTagSF(false, 1));
			}
			//yields with b-tag SF varied up
			for(unsigned dist = 0; dist < nDist; ++dist){
				histsUp[5][dist][fill]->Fill(std::min(values[dist], maxBinC[dist]), (scal/bTagSF(false, 0))*bTagSF(false, 2));
			}
	   	}
		//Set negative bins to 0 before adding other processes
		for(unsigned dist = 0; dist < nDist; ++dist){
			for(unsigned bin = 1; bin < hists[dist][effsam]->GetNbinsX() + 1; ++bin){
				if(hists[dist][effsam]->GetBinContent(bin) < 0.)hists[dist][effsam]->SetBinContent(bin, 0.);
				for(unsigned unc = 0; unc < nUnc; ++unc){
					if(histsDown[unc][dist][effsam]->GetBinContent(bin) < 0.) histsDown[unc][dist][effsam]->SetBinContent(bin, 0.);
					if(histsUp[unc][dist][effsam]->GetBinContent(bin) < 0.) histsUp[unc][dist][effsam]->SetBinContent(bin, 0.);
				}
				for(unsigned pdf = 0; pdf < 100; ++pdf){
					if(histsPdfVar[pdf][dist][effsam]->GetBinContent(bin) < 0.) histsPdfVar[pdf][dist][effsam]->SetBinContent(bin, 0.);
				}
			}
		}
	} 
	std::cout << "crash 1" << std::endl;
	//Set fakes to 0 if they are negative
	for(unsigned dist = 0; dist < nDist; ++dist){
		for(unsigned bin = 1; bin < hists[dist][nSamples_eff]->GetNbinsX() + 1; ++bin){
			if(hists[dist][nSamples_eff]->GetBinContent(bin) < 0) hists[dist][nSamples_eff]->SetBinContent(bin, 0.);
			if(histsDown[6][dist][nSamples_eff]->GetBinContent(bin) < 0) histsDown[6][dist][nSamples_eff]->SetBinContent(bin, 0.);
			if(histsUp[6][dist][nSamples_eff]->GetBinContent(bin) < 0) histsUp[6][dist][nSamples_eff]->SetBinContent(bin, 0.);
		}
	}
	std::cout << "crash 2" << std::endl;
	//Calculate rms of pdf shape variations
	for(unsigned effsam = 1; effsam < nSamples_eff; ++effsam){
		for(unsigned dist = 0; dist < nDist; ++dist){
			for(unsigned b = 1; b < hists[dist][effsam]->GetNbinsX() + 1; ++b){
				double pdfVarRms = 0;
				for(unsigned pdf = 0; pdf < 100; ++pdf){
					pdfVarRms += (hists[dist][effsam]->GetBinContent(b) - histsPdfVar[pdf][dist][effsam]->GetBinContent(b))*(hists[dist][effsam]->GetBinContent(b) - histsPdfVar[pdf][dist][effsam]->GetBinContent(b));
				}
				pdfVarRms = 0.01*sqrt(pdfVarRms);
				histsDown[3][dist][effsam]->SetBinContent(b, hists[dist][effsam]->GetBinContent(b) - pdfVarRms);
				histsUp[3][dist][effsam]->SetBinContent(b, hists[dist][effsam]->GetBinContent(b) + pdfVarRms);
			}
		}
	}
	std::cout << "crash 3" << std::endl;
	//Split data and MC histograms for plotting and propagating uncertainties
	TH1D* dataYields[nDist];
	for(unsigned dist = 0; dist < nDist; ++dist){
		dataYields[dist] = (TH1D*) hists[dist][0]->Clone();
	}
	TH1D* bkgYields[nDist][nSamples_eff]; 
	for(unsigned dist = 0; dist < nDist; ++dist){
		for(unsigned effsam = 1; effsam < nSamples_eff + 1; ++effsam){
			bkgYields[dist][effsam - 1] = (TH1D*) hists[dist][effsam]->Clone();
		}
	}
	std::cout << "crash 4" << std::endl;
	const unsigned nBkg = nSamples_eff;
	const double extraUnc[nBkg] = {1.5, 1.25, 1.5, 1.15, 1.3, 1.3, 1., 1.3}; //extra flat uncertainties assigned to each background
	//Calculate histogram with systematic uncertainty for backgrounds and print the ranges of their size
	double flatSyst[4] = {0.025, 0.02, 0.02}; //1 lepton syst
	TH1D* bkgSyst[nDist][nBkg];
	for(unsigned bkg = 0; bkg < nBkg; ++bkg){
		flatSyst[3] = fabs(extraUnc[bkg] -1);
		for(unsigned dist = 0; dist < nDist; ++dist){
			bkgSyst[dist][bkg] = (TH1D*) bkgYields[dist][bkg]->Clone();
			for(unsigned b = 1;  b < dataYields[dist]->GetNbinsX() + 1; ++b){
				bkgSyst[dist][bkg]->SetBinContent(b,0);
				if(bkg != nBkg -1){
					double systbin = 0;
					//loop over shape uncertainties
					for(unsigned unc = 0; unc < nUnc - 1; ++unc){
						if(unc == 6) continue; //EWK subtraction unc is not defined for MC backgrounds
						double syst = std::max(fabs(histsUp[unc][dist][bkg + 1]->GetBinContent(b) - hists[dist][bkg + 1]->GetBinContent(b)), fabs(histsDown[unc][dist][bkg + 1]->GetBinContent(b) - hists[dist][bkg + 1]->GetBinContent(b)));
						systbin += syst*syst;
					}
					//loop over flat uncertainties
					for(unsigned unc = 0; unc < 4; ++unc){
						systbin += bkgYields[dist][bkg]->GetBinContent(b)*bkgYields[dist][bkg]->GetBinContent(b)*flatSyst[unc]*flatSyst[unc];
					}
					bkgSyst[dist][bkg]->SetBinContent(b, sqrt(systbin));
				} else{
					double syst = std::max(fabs(histsUp[6][dist][bkg + 1]->GetBinContent(b) - hists[dist][bkg + 1]->GetBinContent(b)), fabs(histsDown[6][dist][bkg + 1]->GetBinContent(b) - hists[dist][bkg + 1]->GetBinContent(b)));
					bkgSyst[dist][bkg]->SetBinContent(b, sqrt(hists[dist][bkg + 1]->GetBinContent(b)*0.3*hists[dist][bkg + 1]->GetBinContent(b)*0.3 + syst*syst));
				}
				bkgSyst[dist][bkg]->SetBinError(b, 0);
			}
		}
	}
	std::cout << "crash 5" << std::endl;
	const TString procNames[nSamples_eff + 1] = {"obs.", "TT + X", "rare SM", "Z#gamma", "TT", "WJets", "W#gamma", "T+X", "non-prompt"};
	//Plot the yields as a function of the search region
	for(unsigned dist = 0; dist < nDist; ++dist){
		plotDataVSMC(dataYields[dist], bkgYields[dist], procNames, nSamples_eff, "controlR/" +  distNames[dist] + extra, false, 0, "", bkgSyst[dist]);
	}
}



int main(int argc, char* argv[]){
	TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;
	testtree.Loop();
	rootapp->Run();
    return 0;
}


