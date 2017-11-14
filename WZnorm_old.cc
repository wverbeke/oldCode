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

	
	const unsigned nSamples = 29;
	const unsigned nSamples_eff = 6;
	//TTGJets.root
	const TString fileList[nSamples] = {"data_combined_trilepton.root", "ZZTo4L.root",  "VHToNonbb.root", "WWG.root","WWW.root", "WWZ.root", "WWTo2L2Nu_DoubleScattering.root", "WWTo2L2Nu.root",  "ZZZ.root", "WZTo3LNu_mllmin01.root", "TTGJets.root","ZGTo2LG.root", "WGToLNuG.root", "TGJets.root", "TTJets_DiLept.root", "TTJets_SingleLeptFromTbar.root", "TTJets_SingleLeptFromT.root", "DYJetsToLL_M10to50.root", "DYJetsToLL_M50.root", "ST_tW_antitop_NofullyHadronic.root", "ST_tW_top_NofullyHadronic.root", "ST_s-channel_leptonDecays.root", "ST_t-channel_top_inclusiveDecays.root", "ST_t-channel_antitop_inclusiveDecays.root", "ttHToNonbb.root", "TTWJetsToLNu.root", "TTZToLLNuNu.root", "TTZToLL_M1to10.root",  "TTTT.root"};
	//const double xSections[nSamples - 5] = {0.215, 0.2043, 0.2529, 3.697, 1.256, 0.2147,  0.009103, 0.9561, 0.2086, 0.1651, 0.1729,   12.178, 0.01398, 4.4297, 123.9, 87.315, 182.175, 182.175, 61526.7, 350.674, 2.967, 38.09, 38.09, 10.11, 136.02, 80.95}; 
	//3.697 4.4297
	//*0.696758
	//0.672 ewkino WZ norm
	//0.696758 HNL WZ norm
	const double xSections[nSamples - 1] = {1.256,  0.9561, 0.2147, 0.2086, 0.1651, 0.1729,   12.178, 0.01398, 58.59*0.672, 3.697, 123.9, 489,  2.967, 87.315, 182.175, 182.175, 18610, 1921.8*3, 38.09, 38.09, 10.11, 136.02, 80.95, 0.215, 0.2043, 0.2529, 0.0493, 0.009103}; 



	const TString names[nSamples] = {"data", "ZZ/H", "ZZ/H", "triboson", "triboson", "triboson", "triboson", "triboson", "triboson", "WZ", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X"};
	
	//Read Trees from ROOT files
	TFile* hfile[nSamples];
	TTree* inputTree[nSamples];
	double hcounter[nSamples];
	for(unsigned sam = 0; sam < nSamples; ++sam){
		cout << "name " << names[sam] << endl;
		hfile[sam] = new TFile("../data_reminiaod17/"+fileList[sam],"read");
       	hfile[sam]->cd("FakeElectrons");
		//Determine hcounter for cross section scaling
		TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
		_hCounter->Read("hCounter");
		hcounter[sam] = _hCounter->GetBinContent(1);
       	inputTree[sam] = (TTree*) (hfile[sam]->Get("FakeElectrons/fakeTree"));
		Init(inputTree[sam], false, sam > 4);
	}

	//readSF(true);
	readSF();
	BTagCalibration calib("csvv2", "../bTag/CSVv2_Moriond17_B_H.csv");
	BTagCalibrationReader reader(BTagEntry::OP_LOOSE, "central", {"up", "down"});
	reader.load(calib, BTagEntry::FLAV_B, "comb");
	reader.load(calib, BTagEntry::FLAV_C, "comb");
	reader.load(calib, BTagEntry::FLAV_UDSG, "comb");
	//Tweakable options////////////////////////////////////////////////////
	const bool TestRun = false;	//Break after a few events
	const double DataLuminosity = 35.867;    //units of fb^{-1}
	const TString extra = "_hnl";	//for plot file names
	//////////////////////////

	const TString eff_names[nSamples_eff + 1] = {"data", "ZZ/H", "triboson", "WZ", "X + #gamma", "TT/T + X",  "non-prompt"}; //X + #gamma
	/*
	const unsigned nDist = 5;  //Number of distributions to plot
	TH1D* Histos[nDist][nSamples_eff + 1];
	const TString Histnames[nDist] = {"mt", "met", "mll", "mt_min", "minMos"};
	const TString Xaxes[nDist] = {"M_{T}(GeV)", "MET(GeV)", "M_{ll} (GeV)", "M_{T}(other min(M_{OS})) (GeV)", "min(M_{OS}) (GeV)"};
	const TString Units[nDist] = {"GeV", "GeV", "GeV", "GeV", "GeV"};
	const double HistMin[nDist] = {0, 0, 12, 0, 12}; 
	const double HistMax[nDist] = {300, 300, 300, 300, 300};
	unsigned nBins[nDist];
	for(unsigned dist = 0; dist < nDist; ++dist) nBins[dist] = 20;
	
	for(unsigned dist = 0; dist < nDist; ++dist){
		float BinWidth = (HistMax[dist] - HistMin[dist])/nBins[dist];
		std::ostringstream strs; strs << BinWidth; std::string Yaxis = strs.str();
		for(unsigned effsam = 0; effsam < nSamples_eff + 1; ++effsam){
			Histos[dist][effsam] = new TH1D(eff_names[effsam] + Histnames[dist], eff_names[effsam] + Histnames[dist] + ";" + Xaxes[dist] + "; events /" + Yaxis + Units[dist], nBins[dist], HistMin[dist], HistMax[dist]);
			Histos[dist][effsam]->Sumw2();
		}
	}
	*/
	
	const unsigned nDist = 57;  //Number of distributions to plot	
	TH1D* Histos[nDist][nSamples_eff + 1];
	const TString Histnames[nDist] = {"Mll", "M3l", "minMos", "mt_minMos", "MosminDeltaR", "mt_minDeltaR", "mt", "mt2_ss", "mt2_maxPt", "LeptonPt_le","LeptonPt_sub", "LeptonPt_tr", "MET", "HT", "NJets", "NbJets", "DeltaPhi_lepMET_le", "DeltaPhi_lepMET_sub", "DeltaPhi_lepMET_tr", "Nvtx","Pt_trilep", "Pt_trilepscalarsum", "ConePt_le", "ConePt_sub", "ConePt_tr", "Eta_le", "Eta_sub", "Eta_tr", "MiniIso_le", "MiniIso_sub", "MiniIso_tr", "RelIso_le", "RelIso_sub", "RelIso_tr", "Ptrel_le", "Ptrel_sub", "Ptrel_tr", "Ptratio_le", "Ptratio_sub", "Ptratio_tr", "csv_le", "csv_sub", "csv_tr", "3dIP_le", "3dIP_sub", "3dIP_tr", "dxy_le", "dxy_sub", "dxy_tr", "dz_le", "dz_sub", "dz_tr", "SIP3Dle", "SIP3Dsub", "SIP3Dtr", "ptw", "etaw"};

	const TString Xaxes[nDist] = {"M_{ll}(GeV)", "M_{3l} (GeV)", "min(M_{OS}) (GeV)", "M_{T}(other min(M_{OS})) (GeV)", "M_{OS}(min #Delta R) (GeV)", "M_{T}(other min #Delta R) (GeV)",  "M_{T} (GeV)", "M_{T2}(SS) (GeV)",  "M_{T2}(max P_{T} 2l) (GeV)", "P_{T}(leading l) (GeV)", "P_{T}(subleading l) (GeV)", "P_{T}(trailing l) (GeV)",  "MET (GeV)", "HT (GeV)", "number of jets", "number of b-jets", "#Delta#Phi(leading, MET)", "#Delta#Phi(subleading, MET)", "#Delta#Phi(trailing, MET)", "Number of vertices", "P_{T}(3l) (GeV)", "#Sigma_{l}(|P_{T}(l)|) (GeV)", "P_{T}^{cone}(leading) (GeV)", "P_{T}^{cone}(subleading) (GeV)", "P_{T}^{cone}(trailing) (GeV)", "|#eta(leading)|","|#eta(subleading)|", "|#eta(trailing)|", "miniIso(leading)", "miniIso(subleading)", "miniIso(trailing)", "relIso(leading)", "relIso(subleading)", "relIso(trailing)", "P_{T}^{rel}(leading) (GeV)", "P_{T}^{rel}(subleading) (GeV)", "P_{T}^{rel}(trailing) (GeV)", "P_{T}^{ratio}(leading)", "P_{T}^{ratio}(subleading)", "P_{T}^{ratio}(trailing)", "closest Jet CSV(leading)", "closest Jet CSV(subleading)", "closest Jet CSV(trailing)", "|3DIP(leading)| (cm)", "|3DIP(subleading)| (cm)", "|3DIP(trailing)| (cm)", "|d_{xy}(leading)| (cm)", "|d_{xy}(subleading)| (cm)", "|d_{xy}(trailing)| (cm)", "|d_{z}(leading)| (cm)", "|d_{z}(subleading)| (cm)", "|d_{z}(trailing)| (cm)", "SIP_{3D}(leading)", "SIP_{3D}(subleading)", "SIP_{3D}(trailing)", "P_{T}(W lepton) (GeV)", "|#eta(W lepton)|"};

	const TString Units[nDist] = {"GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "GeV", "", "", "", "", "", "",  "GeV", "GeV", "GeV", "GeV", "GeV", "", "", "", "", "", "", "", "", "",  "GeV", "GeV", "GeV",     "", "", "", "", "", "", "cm", "cm", "cm", "cm", "cm", "cm", "cm", "cm", "cm", "", "", "", "GeV"};
	const double HistMin[nDist] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10, 10, 0, 30, 0, 0, 0, 0, 0, 0, 0, 20, 15, 10, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0}; 
	const double HistMax[nDist] = {300, 600, 150, 300, 300, 300, 300, 300, 300, 200, 200, 200, 300, 600, 10,10, 3.2, 3.2, 3.2, 40, 300, 400, 200, 200, 200, 2.5, 2.5, 2.5, 0.5, 0.5, 0.5,  0.1, 0.1, 0.1, 100, 100, 100, 2, 2, 2, 0.8, 0.8, 0.8, 0.1, 0.1, 0.1, 0.015, 0.015, 0.015, 0.15, 0.15, 0.15, 10, 10, 10, 200, 2.5};
	
	unsigned nBins[nDist];
	for(unsigned dist = 0; dist < nDist; ++dist) nBins[dist] = 20;
	nBins[2] = 15;
	nBins[14] = 10;
	nBins[15] = 10;
	nBins[22] = 37;
	nBins[23] = 38;
	nBins[24] = 39;
	//nBins[1] = 250;


	for(unsigned dist = 0; dist < nDist; ++dist){
		float BinWidth = (HistMax[dist] - HistMin[dist])/nBins[dist];
		std::ostringstream strs; strs << BinWidth; std::string Yaxis = strs.str();
		for(unsigned effsam = 0; effsam < nSamples_eff + 1; ++effsam){
			Histos[dist][effsam] = new TH1D(eff_names[effsam] + Histnames[dist], eff_names[effsam] + Histnames[dist] + ";" + Xaxes[dist] + "; events /" + Yaxis + Units[dist], nBins[dist], HistMin[dist], HistMax[dist]);
			Histos[dist][effsam]->Sumw2();
		}
	}
	double maxBinC[nDist];
	for(unsigned dist = 0; dist < nDist; ++dist){
		maxBinC[dist] = Histos[dist][0]->GetBinCenter(Histos[dist][0]->GetNbinsX());
	}

	
    Double_t scale[nSamples -5];
	//set to check which data events have already been processed
	std::set<std::tuple<unsigned long, unsigned long, unsigned long> > usedEvents; //runNb, lumiBlock, eventNb
	
	//Loop over all samples
	for(unsigned sam = 0, effsam = 0; sam < nSamples; ++sam, ++effsam){
		if(sam != 0){
   			if(names[sam] == names[sam -1]) --effsam;
    	}		

		if(fileList[sam] == "TTJets_DiLept.root" || fileList[sam] == "DYJetsToLL_M10to50.root" || fileList[sam] == "DYJetsToLL_M50.root"  || fileList[sam] == "TTJets_SingleLeptFromTbar.root"  || fileList[sam] == "TTJets_SingleLeptFromT.root" ) continue;
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
			
			//Apply HNL SELECTION
			//cutBased();
			//Baseline event selection:
			unsigned nJets_unc = _nJets;
			//if(!baseline(true, true, false, false)) continue;	
			if(!baseline(true, true,false, false)) continue;	
			//Check if data events were used before
			if(effsam == 0){
				auto event = usedEvents.find(std::make_tuple(_eventNb, _lumiBlock, _runNb));
				if(event != usedEvents.end()) continue;
				usedEvents.insert(std::make_tuple(_eventNb, _lumiBlock, _runNb));
			}
			//Categorize according to the number of leptons and flavors
			unsigned* ind = new unsigned[_nL];
			//CHANGE BACK
			unsigned lCount = lepOrder(ind, 3, true, true);	
			//unsigned lCount = lepOrder(ind, 3, true, false);
			//Tau check: 
			if(lCount != 3) continue; //Veto 4th FO lepton considering signal model!
			//MC prompt matching
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

			if(fileList[sam] == "WGToLNuG.root"){
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

			if(fileList[sam] == "TTJets_DiLept.root" || fileList[sam] == "DYJetsToLL_M10to50.root" || fileList[sam] == "DYJetsToLL_M50.root"  || fileList[sam] == "TTJets_SingleLeptFromTbar.root"  || fileList[sam] == "TTJets_SingleLeptFromT.root" ){
				continue;
				
				bool lowMll = true;
				for(unsigned l = 0; l < _gen_nL -1; ++l){
					for(unsigned k = l + 1; k < _gen_nL; ++k){
						if(_gen_flavors[l] == _gen_flavors[k]){
							if(_gen_charges[l] != _gen_charges[k]){
								TLorentzVector lep1, lep2;
								lep1.SetPtEtaPhiE(_gen_lPt[l], _gen_lEta[l], _gen_lPhi[l], _gen_lE[l]);
								lep2.SetPtEtaPhiE(_gen_lPt[k], _gen_lEta[k], _gen_lPhi[k], _gen_lE[k]);
								if((lep1 + lep2).M() > 30){
									lowMll = false;
									break;
								}
							}
						}
					}
				}
				if(!lowMll){		
					bool promptfail =false;
					for(unsigned l = 0; l < lCount; ++l){	
						//cout << _origin[ind[l]] << endl;
						if(_pdgmc[ind[l]] == 22 && _originPhot[ind[l]] == 0){
							promptfail = true;
							break;
						}
					}
					if(promptfail) continue;
				}
			}
		
			if(fileList[sam] == "TTGJets.root" || fileList[sam] == "ZGToLNuG.root"){
				//if(fileList[sam] == "TTGJets.root") continue;
				/*
				continue;
				bool promptfail =true;
				bool photpass = false;
				for(unsigned l = 0; l < lCount; ++l){	
					//cout << _origin[ind[l]] << endl;
					if(_pdgmc[ind[l]] == 22 && _originPhot[ind[l]] == -2){
						promptfail = false;
						break;
					}
					else if(_pdgmc[ind[l]] == 22 && _originPhot[ind[l]] == 0){
						photpass = true;
					}
				}
				if(promptfail) continue;
				if(!photpass) continue;
				*/
			}
			//Require 3 leptons to be tight in data and MC, and determine nonPrompt bkg in data
			unsigned nTight = tightCount(ind, lCount);
			bool tightFail = nTight < 3;
			//cout << errcount - nTight << endl;
			
			//cout << "lCount = " << lCount << endl;
			//cout << "tightC = " << nTight << endl;
			//index used to fill events, needed to separate fakes from data
			unsigned fill = effsam;
			//Calculate conePt for every lepton
			double* conePt = new double[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				//conePt[l] =  PtCone(_lPt[ind[l]], _flavors[ind[l]], _lepMVA[ind[l]], _ptratio[ind[l]]);
				conePt[l] = _lPt[ind[l]]*std::max(1., 1 + (_isolation[ind[l]] - 0.1));
			}
			//Apply FR maps to data control region
			if(tightFail){ //&& (fileList[sam] != "TTJets_DiLept.root") ){
				//fakes go in different histogram
				fill = nSamples_eff;
				//MC fake subtraction
				if(effsam != 0) scal *= -1.;
				//Apply FR maps
				scal*= fakeWeight(ind, _flavors, conePt, _lEta, _istight, frMap, lCount);
			} 
			
			//Apply triggers
			bool trigPass[4];
			trigPass[0] = _lowM_trigger_eee || _lowM_trigger_all;
			trigPass[1] = _lowM_trigger_mee || _lowM_trigger_all;
			trigPass[2] = _lowM_trigger_mme || _lowM_trigger_all;
			trigPass[3] = _lowM_trigger_mmm || _lowM_trigger_all;
			if(!trigPass[tril_flavorComb(ind,_flavors, lCount)]) continue;
			//determine search category
			unsigned cat = hnl::cat(ind, _flavors, _charges, lCount, conePt[0]);
			if(cat == 999) continue;
			//if(cat == 0 || cat == 2 || cat == 4) continue;
			if(cat != 0 && cat != 2 && cat != 4) continue; //This means there has to be an OSSF pair in the event
			//Apply all efficiency and reweighing SF to the simulation
			if(effsam != 0){
				scal*=getEventSF(ind, lCount, true)
			}
			//determine which leptons will be used for the calculation of mll
			TLorentzVector* lepV = new TLorentzVector[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				lepV[l].SetPtEtaPhiE(conePt[l], _lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]);
				//lepV[l].SetPtEtaPhiE(conePt[l], _lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]*(conePt[l]/_lPt[ind[l]]));
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
				//for(unsigned l = 0; l < 2; ++l) lzV[l].SetPtEtaPhiE(PtCone(_lPt[mllI[l]], _flavors[mllI[l]], _lepMVA[mllI[l]], _ptratio[mllI[l]]), _lEta[mllI[l]], _lPhi[mllI[l]], _lE[mllI[l]]);
				for(unsigned l = 0; l < 2; ++l) lzV[l].SetPtEtaPhiE(_lPt[mllI[l]]*std::max(1., 1 + (_isolation[mllI[l]] - 0.1)), _lEta[mllI[l]], _lPhi[mllI[l]], _lE[mllI[l]]*std::max(1., 1 + (_isolation[mllI[l]] - 0.1)));
				mll = (lzV[0] + lzV[1]).M();
			}
			if((cat == 0 || cat == 2 || cat == 4) && fabs(mll - 91) > 15) continue; //consider onZ events for WZ CR

			
			//determine the index of the W lepton
			unsigned lw = 9999;
			for(unsigned l = 0; l < lCount; ++l){
				if(ind[l] != mllI[0] && ind[l] != mllI[1]){
					lw = ind[l];
				}
			}
			TLorentzVector Wlep;
			Wlep.SetPtEtaPhiE( PtCone(_lPt[lw], _flavors[lw], _lepMVA[lw], _ptratio[lw]), _lEta[lw], _lPhi[lw], _lE[lw]*(PtCone(_lPt[lw], _flavors[lw], _lepMVA[lw], _ptratio[lw])/_lPt[lw]) );
			//Wlep.SetPtEtaPhiE(_lPt[lw]*std::max(1., 1 + (_isolation[lw] - 0.1)), _lEta[lw], _lPhi[lw], _lE[lw]*std::max(1., 1 + (_isolation[lw] - 0.1)));
			TLorentzVector METvec;
			METvec.SetPtEtaPhiE(_met, 0, _met_phi, _met);
			/*
			//MISPAIRING, REMOVE LATER
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

			//new mll
			TLorentzVector lzV[2];
			for(unsigned l = 0; l < 2; ++l) lzV[l].SetPtEtaPhiE(PtCone(_lPt[mllI[l]], _flavors[mllI[l]], _lepMVA[mllI[l]], _ptratio[mllI[l]]), _lEta[mllI[l]], _lPhi[mllI[l]], _lE[mllI[l]]);
			mll = (lzV[0] + lzV[1]).M();
			Wlep.SetPtEtaPhiE( PtCone(_lPt[lw], _flavors[lw], _lepMVA[lw], _ptratio[lw]), _lEta[lw], _lPhi[lw], _lE[lw]*(PtCone(_lPt[lw], _flavors[lw], _lepMVA[lw], _ptratio[lw])/_lPt[lw]) );
			*/
			
			//Calculate lepton system vector
			TLorentzVector lepSyst;
			for(int l = 0; l < 3; ++l) lepSyst += lepV[l];
			//Category spefific selection:
			//inclusive pt thresholds
			if(conePt[0] < 25) continue;
			if(conePt[1] < 15) continue;
			//if(conePt[2] < 10 - 5*(_flavors[ind[2]])) continue;	//5(10) GeV cut on trailing muon(electron)
			if(conePt[2] < 10) continue;// - 5*(_flavors[ind[2]]) ) continue;
			if((cat == 0 || cat == 2 || cat == 4) && fabs(lepSyst.M() - 91) < 15) continue;



			//if(_flavors[ind[2]] != 0) continue;
			//truth matching
			/*
			if(effsam != 0){
				unsigned* match = new unsigned[_nL];
			    matchGenLep(match);
				bool promptfail = false;
				for(unsigned l = 0; l < 3; ++l){
					if(_flavors[ind[l]] == 2) continue;
					if(match[ind[l]] > 20){
						promptfail = true;
						break;
					}
					if( !_gen_isPromptl[match[ind[l]]] && fabs(_gen_lmompdg[match[ind[l]]]) != 15){
						promptfail = true;
						break;
					}
				}
				if(promptfail) continue;
			}
			*/
			/*
			if(fileList[sam] == "TTJets_DiLept.root"){
				unsigned* match = new unsigned[_nL];
			    matchGenLep(match);
				bool promptfail = false;
				unsigned photCount = 0;
				for(unsigned l = 0; l < 3; ++l){
					if(_flavors[ind[l]] == 2) continue;
					if(match[ind[l]] > 20){
						promptfail = true;
						break;
					}
					if(fabs(_gen_lmompdg[match[ind[l]]]) == 22 || fabs(_gen_lmompdg[match[ind[l]]]) == 23) ++photCount;
					if( !_gen_isPromptl[match[ind[l]]] && fabs(_gen_lmompdg[match[ind[l]]]) != 15){
						promptfail = true;
						break;
					}
				}
				if(promptfail) continue;
				//if(photCount < 1) continue;
			}
			*/
			//Calculate MET vector
			//TLorentzVector METvec;
			//METvec.SetPtEtaPhiE(_met, 0, _met_phi, _met);

			//REMOVE LATER
			//if(transmass(METvec, Wlep) < 50 || transmass(METvec, Wlep) > 100) continue;
			//Calculate min(Mos) and Mos(min Delta R)
			double minMos = 0;
			unsigned minI[2] = {99, 99};
			double MosminDeltaR;
			unsigned minDeltaRI[2] = {99, 99};
			double minDeltaR = 99999.;
			for(unsigned l = 0; l < lCount -1 ; ++l){
				for(unsigned k = l + 1; k < lCount; ++k){
					if(_charges[ind[l]] != _charges[ind[k]]){
						if( (lepV[l] + lepV[k]).M() < minMos  || minMos == 0){
							minMos = (lepV[l] + lepV[k]).M();
							minI[0] = l;
							minI[1] = k;
						}
						if( lepV[l].DeltaR(lepV[k]) < minDeltaR){
							minDeltaR = lepV[l].DeltaR(lepV[k]);
							MosminDeltaR = (lepV[l] + lepV[k]).M();
							minDeltaRI[0] = l;
							minDeltaRI[1] = k;
						}
					}
				}
			}
			//Calculate alternete mt values
			unsigned lw_min, lw_minDeltaR;
			for(unsigned l = 0; l < lCount; ++l){
				if(l != minI[0] && l != minI[1]){
					lw_min = l;
				}
				if(l != minDeltaRI[0] && l != minDeltaRI[1]){
					lw_minDeltaR = l;
				}
			}
			double mt_min = transmass(lepV[lw_min], METvec);
			double mt_minDeltaR = transmass(lepV[lw_minDeltaR], METvec);


			//Calculate mt3l
			double mt3l = transmass(lepV[0] + lepV[1] + lepV[2], METvec);
			//calculate mt2_ss
			double mt2_ss = transmass(Wlep, METvec);
			for(unsigned l = 0; l < lCount -1; ++l){
				for(unsigned k = l + 1; k < lCount; ++k){
					if(_charges[ind[k]] == _charges[ind[l]]){
						if(_flavors[ind[k]] == _flavors[ind[l]]){
							mt2_ss = mt2ll(lepV[k], lepV[l], METvec);
						}
					}
				}
			}
			//Calculate scalar sum of lepton Pt's
			double ptSum = 0;
			for(unsigned l = 0; l < lCount; ++l){
				ptSum += _lPt[ind[l]];
			}

			//fill 1D histograms
			/*
			double values[nDist] = { transmass(Wlep, METvec),  _met, mll, mt_min, minMos};

			for(unsigned dist = 0; dist < nDist; ++dist){
				Histos[dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal);
			} 
			*/
			//if(minMos > 10) continue;
			double values[nDist] = { mll, lepSyst.M(), minMos, mt_min, MosminDeltaR, mt_minDeltaR, transmass(Wlep, METvec), mt2_ss,  mt2_maxPt(ind, _charges, lepV, METvec, lCount), _lPt[ind[0]], _lPt[ind[1]], _lPt[ind[2]], _met, _HT, static_cast<double>(_nJets), static_cast<double>(_n_bJets), fabs(METvec.DeltaPhi(lepV[0])), fabs(METvec.DeltaPhi(lepV[1])), fabs(METvec.DeltaPhi(lepV[2])), static_cast<double>(_n_PV), lepSyst.Pt(), ptSum,  _lPt[ind[0]]*std::max(1., 1 + (_isolation[ind[0]] - 0.1)), _lPt[ind[1]]*std::max(1., 1 + (_isolation[ind[1]] - 0.1)), _lPt[ind[2]]*std::max(1., 1 + (_isolation[ind[2]] - 0.1)), fabs(_lEta[ind[0]]), fabs(_lEta[ind[1]]), fabs(_lEta[ind[2]]),  _miniisolation[ind[0]][0], _miniisolation[ind[1]][0], _miniisolation[ind[2]][0], _isolation[ind[0]],  _isolation[ind[1]],  _isolation[ind[2]], _ptrel[ind[0]], _ptrel[ind[1]], _ptrel[ind[2]], _ptratio[ind[0]], _ptratio[ind[1]], _ptratio[ind[2]], _closeJetCSVAll[ind[0]], _closeJetCSVAll[ind[1]], _closeJetCSVAll[ind[2]], fabs(_3dIP[ind[0]]), fabs(_3dIP[ind[1]]), fabs(_3dIP[ind[2]]), fabs(_ipPV[ind[0]]), fabs(_ipPV[ind[1]]), fabs(_ipPV[ind[2]]), fabs(_ipZPV[ind[0]]), fabs(_ipZPV[ind[1]]), fabs(_ipZPV[ind[2]]), _3dIPsig[ind[0]],  _3dIPsig[ind[1]], _3dIPsig[ind[2]], _lPt[lw]*std::max(1., 1 + (_isolation[lw] - 0.1)), fabs(_lEta[lw])}; 
	
			for(unsigned dist = 0; dist < nDist; ++dist){
				Histos[dist][fill]->Fill(TMath::Min(values[dist], maxBinC[dist]), scal);
			}
    	}
	}

	cout << "leading entries = " << Histos[22][3]->GetBinContent(1) << endl;
	cout << "subleading entries = " << Histos[23][3]->GetBinContent(1) << endl;
	cout << "trailing entries = " << Histos[24][3]->GetBinContent(1) << endl;
	cout << "Lower bin edge = " <<  Histos[24][3]->GetBinLowEdge(1) << endl;
	cout << "bin width = " << Histos[24][3]->GetBinWidth(1) << endl;
	//Split data and MC histograms for plotting and propagating uncertainties
	TH1D* dataYields[nDist];
	for(unsigned dist = 0; dist < nDist; ++dist){
		dataYields[dist] = (TH1D*) Histos[dist][0]->Clone();
	}
	
	TH1D* bkgYields[nDist][nSamples_eff]; //change to nSamples_eff if sig is removed
	for(unsigned dist = 0; dist <nDist; ++dist){
		for(unsigned effsam = 1; effsam < nSamples_eff + 1; ++effsam){
			bkgYields[dist][effsam -1] = (TH1D*) Histos[dist][effsam]->Clone();
		}
	}
	//Produce datacards for WZ normalization
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//PRINT DATA CARDS FOR LIMIT CALCULATION
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//const unsigned nSig = 13;
	const unsigned nBkg = nSamples_eff - 1; //WZ is our signal
	const unsigned nSyst = 5 + 1 + 2*nBkg; //5 general uncertainties, one stat unc for signal, stat + extra unc for every bkg;
	const TString bkgNames[nBkg] = {"ZZH", "triboson", "Xgamma", "TTX",  "nonPrompt"};//slightly rewrite bkg names not to confuse the combine tool
	std::vector<std::vector<double>> systUnc(nSyst, std::vector<double>(nBkg + 1, 0)); //2D array containing all systematics
	//initialize values non-bin-dependent systematics
	for(unsigned proc = 0; proc < nBkg; ++proc){ //last background is non-prompt and isn't susceptible to these uncertainty sources (else nBkg + 1 in loop_
		systUnc[0][proc] = 1.026;	//lumi
		systUnc[1][proc] = 1.04;	//id eff
		systUnc[2][proc] = 1.02;	//trig eff
		systUnc[3][proc] = 1.05;	//JEC	
		systUnc[4][proc] = 1.05;	//PU
	}
	for(unsigned syst = 0; syst < 4; ++syst) systUnc[syst][nBkg] = 0; //assign 0 for non-prompt
	//names of uncertainties;
	TString systNames[nSyst] = {"lumi", "id_eff", "trigeff", "JEC", "PU"};
	TString systDist[nSyst];
	for(unsigned syst = 0; syst < nSyst; ++syst) systDist[syst] = "lnN";
	const double extraUnc[nBkg] = {1.25, 1.5, 1.15, 1.15, 1.3}; //extra unc assigned to backgrounds
	for(unsigned syst = 6 + nBkg; syst < nSyst; ++syst){
		unsigned bkg = syst - 6 - nBkg;
		systNames[syst] = "extra" + bkgNames[bkg];
		systUnc[syst][bkg + 1] = extraUnc[bkg];
	}
	//Print one datacard with the total yields to be normalized
	//signal stat unc
	systUnc[5][0] = 1 + std::max(0., Histos[15][3]->GetBinError(1)/Histos[15][3]->GetBinContent(1));
	systNames[5] = "statSig";
	//bkg stat uncertainties
	double bkgYieldVal[nBkg];
	unsigned bkgCount = 0;
	for(unsigned bkg = 0; bkg < nBkg + 1; ++bkg){
		if(bkg == 2) continue;
		bkgYieldVal[bkgCount] = Histos[15][bkg + 1]->GetBinContent(1);
		systUnc[6 + bkgCount][bkgCount + 1]  = 1 + std::max(0., Histos[15][bkg + 1]->GetBinError(1)/Histos[15][bkg + 1]->GetBinContent(1));
		systNames[6 + bkgCount] = "stat" + bkgNames[bkgCount];
		++bkgCount;
	} 
	hnl::printDataCard(Histos[15][0]->GetBinContent(1), Histos[15][3]->GetBinContent(1), "WZ", bkgYieldVal, nBkg, bkgNames, systUnc, nSyst, systNames, systDist, "datacards/datacard_WZnormalization");
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


	
	for(unsigned dist = 0; dist < nDist; ++dist){
		plotDataVSMC(dataYields[dist], bkgYields[dist], eff_names, nSamples_eff, Histnames[dist] + extra, false, 0, "EWKino");
	}


}



int main(int argc, char* argv[]){
	TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;
	testtree.Loop();
	rootapp->Run();
    return 0;
}


