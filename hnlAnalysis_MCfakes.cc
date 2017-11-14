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

	const unsigned nSamples = 32;
	const unsigned nSamples_eff = 8;
	const unsigned nSig = 0;
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
 										"ZZTo4L.root",  "GluGluToZZTo4mu.root", "GluGluToZZTo4e.root", "GluGluToZZTo4tau.root", "GluGluToZZTo2e2mu.root", "GluGluToZZTo2e2tau.root", "GluGluToZZTo2mu2tau.root", "VHToNonbb.root", "GluGluHToZZTo4L_M125.root", "VBF_HToZZTo4L_M125.root", "WWG.root","WWW.root", "WWZ.root", "WWTo2L2Nu_DoubleScattering.root", "WWTo2L2Nu.root",  "ZZZ.root", "WZTo3LNu_mllmin01.root", "TTGJets.root","ZGTo2LG.root", "WGToLNuG.root", "TGJets.root", "TTJets_DiLept.root", "TTJets_SingleLeptFromTbar.root", "TTJets_SingleLeptFromT.root", "DYJetsToLL_M10to50.root", "DYJetsToLL_M50.root", "ttHToNonbb.root", "TTWJetsToLNu.root", "TTZToLLNuNu.root", "TTZToLL_M1to10.root",  "TTTT.root"};
	
	const double xSections[nSamples - 1] = {1.256*ZZSF, 0.00159*glugluToZZkFactor*ZZSF, 0.00159*glugluToZZkFactor*ZZSF, 0.00159*glugluToZZkFactor*ZZSF, 0.00319*glugluToZZkFactor*ZZSF, 0.00319*glugluToZZkFactor*ZZSF, 0.00319*glugluToZZkFactor*ZZSF,  0.9561, 0.01212, 0.001034,  0.2147, 0.2086, 0.1651, 0.1729,   12.178, 0.01398, 58.59*WZSF, 3.697, 123.9*XgammaSF, 489*XgammaSF,  2.967, 87.315, 182.175, 182.175, 18610, 1921.8*3,  0.215, 0.2043, 0.2529, 0.0493, 0.009103};
	const TString names[nSamples] = {"data", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "ZZ/H", "triboson", "triboson", "triboson", "triboson", "triboson", "triboson", "WZ", "X + #gamma", "X + #gamma", "X + #gamma", "X + #gamma", "TT", "TT", "TT", "Drell-Yan", "Drell-Yan", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X", "TT/T + X"};

	
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
	const TString extra = "_MCfakes";	//for plot file names
	//////////////////////////
	
	const TString eff_names[nSamples_eff + 1] = {"data", "ZZ/H", "triboson", "WZ", "X + #gamma", "TT", "Drell-Yan", "TT/T + X"};
	const unsigned nCat = 6;  //Number of categories
	const TString catNames[nCat] = {"lowM_3lOSSF_lowPt", "lowM_3lnoOSSF_lowPt", "lowM_3lOSSF_highPt", "lowM_3lnoOSSF_highPt", "highM_3lOSSF", "highM_3lnoOSSF"};
	const unsigned nSR[nCat] = {12, 4, 12, 4, 16, 9}; //numbers of search regions
	const unsigned nUnc = 8;//number of shape uncertainties
	const TString uncNames[nUnc] = {"jec", "metUnclustered",  "scale", "pdf", "pu", "btagSF", "trigeff", "id_eff"};
	TH1D* yields[nCat][nSamples_eff]; //nominal yields in every SR
	TH1D* yieldsDown[nUnc][nCat][nSamples_eff]; //yields varied down by shape unc
	TH1D* yieldsUp[nUnc][nCat][nSamples_eff]; //yields varied up by shape unc
	TH1D* yieldsPdfVar[100][nCat][nSamples_eff]; //yields for all 100 possible pdf variations

	for(unsigned cat = 0; cat < nCat; ++cat){
		for(unsigned effsam = 0; effsam < nSamples_eff; ++effsam){
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
			if(tightFail) continue;

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
				//if(fileList[sam] == "TTJets_DiLept.root" || fileList[sam] == "DYJetsToLL_M10to50.root" || fileList[sam] == "DYJetsToLL_M50.root"  || fileList[sam] == "TTJets_SingleLeptFromTbar.root"  || fileList[sam] == "TTJets_SingleLeptFromT.root" ) continue;
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

			//Trigger uncertainties
			if(conePt[0] > 30){
				yieldsDown[6][cat][fill]->Fill(searchR + 1, scal*0.98);
				yieldsUp[6][cat][fill]->Fill(searchR + 1, scal*1.02);
			} else{
				yieldsDown[6][cat][fill]->Fill(searchR + 1, scal*0.95);
				yieldsUp[6][cat][fill]->Fill(searchR + 1, scal*1.05);
			}
			//Id efficiency uncertainties
			unsigned flavorC = tril_flavorComb(ind, _flavors, lCount);
			if(flavorC == 0 || flavorC == 3){
				yieldsDown[7][cat][fill]->Fill(searchR + 1, scal*0.94);
				yieldsUp[7][cat][fill]->Fill(searchR + 1, scal*1.06);
			} else{
				yieldsDown[7][cat][fill]->Fill(searchR + 1, scal*0.9553);
				yieldsUp[7][cat][fill]->Fill(searchR + 1, scal*1.0447);	
			}
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
	//Calculate rms of pdf shape variations
	for(unsigned cat = 0; cat < nCat; ++cat){
		for(unsigned effsam = 1; effsam < nSamples_eff - 1; ++effsam){
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

	//Split data and MC histograms for plotting and propagating uncertainties
	TH1D* dataYields[nCat];
	for(unsigned cat = 0; cat < nCat; ++cat){
		dataYields[cat] = (TH1D*) yields[cat][0]->Clone();
	}
	TH1D* bkgYields[nCat][nSamples_eff - 1]; 
	TH1D* bkgErros[nCat][nSamples_eff - 1];
	for(unsigned cat = 0; cat < nCat; ++cat){
		for(unsigned effsam = 1; effsam < nSamples_eff; ++effsam){
			bkgYields[cat][effsam - 1] = (TH1D*) yields[cat][effsam]->Clone();
		}
	}
	
	const unsigned nBkg = nSamples_eff - 1;
	const double extraUnc[nBkg] = {1.25, 1.5, 1.094, 1.15, 1.3, 1.3, 1.5,}; //extra flat uncertainties assigned to each background
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	//Calculate histogram with systematic uncertainty for backgrounds and print the ranges of their size
	double flatSyst[2] = {0.025};
	TH1D* bkgSyst[nCat][nSamples_eff - 1];
	for(unsigned bkg = 0; bkg < nBkg; ++bkg){
		//flatSyst[4] = fabs(pdfAccUnc[bkg + nSig + 1] - 1);
		flatSyst[1] = fabs(extraUnc[bkg] -1);
		for(unsigned cat = 0; cat < nCat; ++cat){
			bkgSyst[cat][bkg] = (TH1D*) bkgYields[cat][bkg]->Clone();
			for(unsigned b = 1;  b < dataYields[cat]->GetNbinsX() + 1; ++b){
				bkgSyst[cat][bkg]->SetBinContent(b,0);
				double systbin = 0;
				//loop over shape uncertainties
				for(unsigned unc = 0; unc < nUnc - 1; ++unc){
					if(unc == 8) continue; //EWK subtraction unc is not defined for MC backgrounds
					double syst = std::max(fabs(yieldsUp[unc][cat][bkg + nSig + 1]->GetBinContent(b) - yields[cat][bkg + nSig + 1]->GetBinContent(b)), fabs(yieldsDown[unc][cat][bkg + nSig + 1]->GetBinContent(b) - yields[cat][bkg + nSig + 1]->GetBinContent(b)));
					systbin += syst*syst;
				}
				//loop over flat uncertainties
				for(unsigned unc = 0; unc < 2; ++unc){
					systbin += bkgYields[cat][bkg]->GetBinContent(b)*bkgYields[cat][bkg]->GetBinContent(b)*flatSyst[unc]*flatSyst[unc];
				}
				bkgSyst[cat][bkg]->SetBinContent(b, sqrt(systbin));
				bkgSyst[cat][bkg]->SetBinError(b, 0);
			}
		}
	}
	
	//const TString distNames[nSamples_eff + 1 - nSig] = {"total pred.", "ZZ/H", "triboson", "WZ", "X + #gamma", "TT/T + X",  "non-prompt"};
	const TString distNames[nSamples_eff] = {"obs.", "ZZ/H", "triboson", "WZ", "X + #gamma", "TT", "Drell-Yan", "TT/T + X"};
	//Plot the yields as a function of the search region
	for(unsigned cat = 0; cat < nCat; ++cat){
		plotDataVSMC(dataYields[cat], bkgYields[cat], distNames, nSamples_eff - 1, catNames[cat] + extra, true, 0, "HNL", bkgSyst[cat]);
	}
	/*
	//Make one histogram containing all search regions;
	unsigned totalSR = 0;
	for(unsigned cat = 0; cat < nCat; ++cat){
		if(cat == 0 || cat == 2) continue;
		totalSR += nSR[cat];
	}
	TH1D* totalDataYield = new TH1D("totaldata", "totaldata + ; search region ; Events/search region", totalSR, 0.5, totalSR + 0.5);
	unsigned binCounter = 1;
	for(unsigned cat = 0; cat < nCat; ++cat){
		if(cat == 0 || cat == 2 || cat == 4) continue;
		for(unsigned bin = 1; bin < dataYields[cat]->GetNbinsX() + 1; ++bin){
			totalDataYield->SetBinContent(binCounter, dataYields[cat]->GetBinContent(bin) );
			totalDataYield->SetBinError(binCounter, dataYields[cat]->GetBinError(bin) );
			++binCounter;
		}
	}
	for(unsigned bin = 1; bin < dataYields[4]->GetNbinsX() + 1; ++bin){
		totalDataYield->SetBinContent(binCounter, dataYields[4]->GetBinContent(bin) );
		totalDataYield->SetBinError(binCounter, dataYields[4]->GetBinError(bin) );
		++binCounter;
	}
	TH1D* totalBkgYields[nBkg];
	TH1D* totalBkgError[nBkg];
	for(unsigned bkg = 0; bkg < nBkg; ++bkg){
		binCounter = 1;
		totalBkgYields[bkg] = new TH1D(backGroundNames[bkg], backGroundNames[bkg] + "; search region ; Events/search region", totalSR, 0.5, totalSR + 0.5);
		totalBkgError[bkg] = new TH1D(backGroundNames[bkg] + "error",backGroundNames[bkg] + "error; search region ; Events/search region", totalSR, 0.5, totalSR + 0.5);
		for(unsigned cat = 0; cat < nCat; ++cat){
			if(cat == 0 || cat == 2 || cat == 4) continue;
			for(unsigned bin = 1; bin < bkgYields[cat][bkg]->GetNbinsX() + 1; ++bin){
				totalBkgYields[bkg]->SetBinContent(binCounter, bkgYields[cat][bkg]->GetBinContent(bin) );
				totalBkgYields[bkg]->SetBinError(binCounter, bkgYields[cat][bkg]->GetBinError(bin) );
				totalBkgError[bkg]->SetBinContent(binCounter, bkgSyst[cat][bkg]->GetBinContent(bin) );
				totalBkgError[bkg]->SetBinError(binCounter, 0);
				++binCounter;
			}
		}
		for(unsigned bin = 1; bin < bkgYields[4][bkg]->GetNbinsX() + 1; ++bin){
			totalBkgYields[bkg]->SetBinContent(binCounter, bkgYields[4][bkg]->GetBinContent(bin) );
			totalBkgYields[bkg]->SetBinError(binCounter, bkgYields[4][bkg]->GetBinError(bin) );
			totalBkgError[bkg]->SetBinContent(binCounter, bkgSyst[4][bkg]->GetBinContent(bin) );
			totalBkgError[bkg]->SetBinError(binCounter, 0);
			++binCounter;
		}
	}
	makePaperPlotHNL(totalDataYield, totalBkgYields, totalBkgError, backGroundNames, nBkg, "3 leptons, all events", "allSR");
	*/
}



int main(int argc, char* argv[]){
	TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;
	testtree.Loop();
	rootapp->Run();
    return 0;
}


