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
#include "TGraphErrors.h"

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
#include "hnlTools.h"
//#include "drawLumi.h"


void trilTree::Loop(){
	//Set plotting style
	setTDRStyle();
	gROOT->SetBatch(kTRUE);

	//Initialize all samples and cross sections

	const unsigned nSamples = 19; //18
	/*
	const TString fileList[nSamples] = {"HeavyNeutrino_M1_mu.root", "HeavyNeutrino_M2_mu.root", "HeavyNeutrino_M5_mu.root", "HeavyNeutrino_M10_mu.root", "HeavyNeutrino_M20_mu.root", "HeavyNeutrino_M30_mu.root", "HeavyNeutrino_M40_mu.root", "HeavyNeutrino_M50_mu.root", "HeavyNeutrino_M60_mu.root", "HeavyNeutrino_M80_mu.root", "HeavyNeutrino_M100_mu.root", "HeavyNeutrino_M130_mu.root", "HeavyNeutrino_M150_mu.root", "HeavyNeutrino_M200_mu.root", "HeavyNeutrino_M400_mu.root", "HeavyNeutrino_M600_mu.root", "HeavyNeutrino_M800_mu.root", "HeavyNeutrino_M1000_mu.root", "HeavyNeutrino_trilepton_M-5_V-0.00316_mu.root"};
	*/
	/*
	const TString fileList[nSamples] = {"HeavyNeutrino_M1_2l.root", "HeavyNeutrino_M2_2l.root", "HeavyNeutrino_M5_2l.root", "HeavyNeutrino_M10_2l.root", "HeavyNeutrino_M20_2l.root", "HeavyNeutrino_M30_2l.root", "HeavyNeutrino_M40_2l.root", "HeavyNeutrino_M50_2l.root", "HeavyNeutrino_M60_2l.root", "HeavyNeutrino_M80_2l.root", "HeavyNeutrino_M100_2l.root", "HeavyNeutrino_M130_2l.root", "HeavyNeutrino_M150_2l.root", "HeavyNeutrino_M200_2l.root", "HeavyNeutrino_M400_2l.root", "HeavyNeutrino_M600_2l.root", "HeavyNeutrino_M800_2l.root", "HeavyNeutrino_M1000_2l.root", "HeavyNeutrino_trilepton_M-5_V-0.00316_2l.root"};
	*/
	
	const TString fileList[nSamples] = {"HeavyNeutrino_M1_e.root", "HeavyNeutrino_M2_e.root", "HeavyNeutrino_M5_V0p00507_mu.root", "HeavyNeutrino_M10_e.root", "HeavyNeutrino_M20_e.root", "HeavyNeutrino_M30_e.root", "HeavyNeutrino_M40_e.root", "HeavyNeutrino_M50_e.root", "HeavyNeutrino_M60_e.root", "HeavyNeutrino_M80_e.root", "HeavyNeutrino_M100_e.root", "HeavyNeutrino_M130_e.root", "HeavyNeutrino_M150_e.root", "HeavyNeutrino_M200_e.root", "HeavyNeutrino_M400_e.root", "HeavyNeutrino_M600_e.root", "HeavyNeutrino_M800_e.root", "HeavyNeutrino_M1000_e.root", "HeavyNeutrino_M5_V00p00336_mu.root"}; //HeavyNeutrino_trilepton_M-5_V-0.00316_mu.root HeavyNeutrino_M5_V0p00507_mu.root
	
	const TString names[nSamples] = {"m_{N} = 1 GeV", "m_{N} = 2 GeV", "m_{N} = 5 GeV", "m_{N} = 10 GeV", "m_{N} = 20 GeV", "m_{N} = 30 GeV", "m_{N} = 40 GeV", "m_{N} = 50 GeV", "m_{N} = 60 GeV", "m_{N} = 80 GeV", "m_{N} = 100 GeV", "m_{N} = 130 GeV", "m_{N} = 150 GeV", "m_{N} = 200 GeV", "m_{N} = 400 GeV", "m_{N} = 600 GeV", "m_{N} = 800 GeV", "m_{N} = 1000 GeV", "m5 displaced"};
	
	//Read Trees from ROOT files
	TFile* hfile[nSamples];
	TTree* inputTree[nSamples];
	double hcounter[nSamples];
	for(unsigned sam = 0; sam < nSamples; ++sam){
		cout << "name " << names[sam] << endl;
		hfile[sam] = new TFile("../unskimmedSignal/"+fileList[sam],"read");
       	hfile[sam]->cd("FakeElectrons");
		//Determine hcounter for cross section scaling
		TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
		_hCounter->Read("hCounter");
		hcounter[sam] = _hCounter->GetBinContent(1);
       	inputTree[sam] = (TTree*) (hfile[sam]->Get("FakeElectrons/fakeTree"));
		Init(inputTree[sam], false, true);
	}

	//Tweakable options////////////////////////////////////////////////////
	const bool TestRun = false;	//Break after a few events
	const TString extra = "";
	//////////////////////////

	const unsigned nDist = 6;  //Number of distributions to plot
	TH1D* Histos[nDist][nSamples];
	const TString Histnames[nDist] = {"isole", "isosub", "isotr", "leadingGenPt", "subleadingGenPt", "trailingGenPt"};

	const TString Xaxes[nDist] = {"isolation", "isolation", "isolation", "P_{T}(GeV)", "P_{T}(GeV)", "P_{T}(GeV)"};

	const TString Units[nDist] = {"", "", "", "GeV", "GeV", "GeV"};
	const double HistMin[nDist] = {0,0,0, 0, 0, 0}; 
	const double HistMax[nDist] = {0.7,0.7,0.7, 100, 100, 100};
	unsigned nBins[nDist];
	for(unsigned dist = 0; dist < 3; ++dist) nBins[dist] = 7;
	for(unsigned dist = 3; dist < nDist; ++dist) nBins[dist] = 100;

	for(unsigned dist = 0; dist < nDist; ++dist){
		float BinWidth = (HistMax[dist] - HistMin[dist])/nBins[dist];
		std::ostringstream strs; strs << BinWidth; std::string Yaxis = strs.str();
		for(unsigned sam = 0; sam < nSamples; ++sam){
			Histos[dist][sam] = new TH1D(names[sam] + Histnames[dist],names[sam] + Histnames[dist] + ";" + Xaxes[dist] + "; Events /" + Yaxis + Units[dist], nBins[dist], HistMin[dist], HistMax[dist]);
			Histos[dist][sam]->Sumw2();
		}
	}

	
	const double mass[nSamples] = {1, 2, 5, 10, 20, 30, 40, 50, 60, 80, 100, 130, 150, 200, 400, 600, 800, 1000};
	
	TH1D* lowMeff[2];
	lowMeff[0] = new TH1D("lowM efficiency num", "", nSamples, 0, nSamples);
	lowMeff[1] = new TH1D("lowM efficiency denom", "", nSamples, 0, nSamples);
	lowMeff[0]->Sumw2();
	lowMeff[1]->Sumw2();
	TH1D* highMeff[2];
	highMeff[0] = new TH1D("highM efficiency num", "", nSamples, 0, nSamples);
	highMeff[1] = new TH1D("highM efficiency denom", "",nSamples, 0, nSamples);
	highMeff[0]->Sumw2();
	highMeff[1]->Sumw2();
	

	double tauCounter[nSamples];
	double totalCounter[nSamples];
	for(unsigned sam = 0; sam < nSamples; ++sam){
		tauCounter[sam] = 0.;
		totalCounter[sam] = 0.;
	}
	double denom;
	double numPrompt;	
	double numHard;
	
	int displacedCounter = 0;
	//Loop over all samples
	for(unsigned sam = 0; sam < nSamples; ++sam){
		if(sam != 2 && sam != 18) continue;	
    	Long64_t nEntries = inputTree[sam]->GetEntries();
    	std::cout<<"Entries in "<< fileList[sam] <<" "<<nEntries<<std::endl;

        for(Long64_t it = 0; it < nEntries; ++it){
        	inputTree[sam]->GetEntry(it);
        	if (it%10000 == 0) cout<<'.'<<flush;
        	if(TestRun && it > 10000) break;
			
			totalCounter[sam] += 1.;
			cutBased();
			//Order particles by gen lepton Pt and then store them
			if(_gen_nL < 3) continue;

			unsigned nMu = 0; 
			unsigned nEle = 0;
			
			//if(_gen_nL > 2){
				unsigned* genInd = new unsigned[_gen_nL];
				unsigned* ordind = new unsigned[_gen_nL];
				std::set<unsigned> usedLep;
				for(unsigned k =0; k < _gen_nL; ++k){
					double maxPt = 0;
					for(unsigned l = 0; l < _gen_nL; ++l){
						if(usedLep.find(l) == usedLep.end()){
							double conePt;
							if(_gen_lPt[l] > maxPt){
								maxPt = _gen_lPt[l];
								ordind[k] = l;
							}
						}
					}
					usedLep.insert(ordind[k]);
				}
				for(unsigned i = 0; i < _gen_nL; ++i){
					genInd[i] = ordind[i];
				}
				for(unsigned l = 0; l < 3; ++l){
					if(_gen_flavors[genInd[l]] == 0) ++nEle;
					else if(_gen_flavors[genInd[l]] == 1) ++nMu;
				}
				//if(nMu < 2) continue;
				//if(_gen_flavors[genInd[2]] == 0) continue;
				/*
				Histos[3][sam]->Fill(_gen_lPt[genInd[0]],1); 	
				Histos[4][sam]->Fill(_gen_lPt[genInd[1]],1); 			
				Histos[5][sam]->Fill(_gen_lPt[genInd[2]],1); 	
				*/
			//}	
			/*
			if(_gen_lPt[genInd[0]] < 15) continue;
			if(_gen_lPt[genInd[1]] < 10) continue;
			if(_gen_lPt[genInd[2]] < 10 - 5*_gen_flavors[genInd[2]]) continue;
			
			//if( tril_flavorComb(genInd, _gen_flavors, 3) != 0) continue;
			//Additional Pt cuts required by trigger path
			if(tril_flavorComb(genInd, _gen_flavors, 3) == 0){ //eee
				if(_gen_lPt[genInd[0]] < 19) continue;
				if( (_gen_lPt[genInd[0]] < 30) && (_gen_lPt[genInd[1]] < 15) ) continue;
			} else if(tril_flavorComb(genInd, _gen_flavors, 3) == 1){ //eem
				if(_flavors[genInd[2]] == 0){
					if(_gen_lPt[genInd[2]] < 15){
						if(_gen_lPt[genInd[0]] < 23) continue; //23
					}
				} else if(_flavors[genInd[2]] == 1){
					if(_gen_lPt[genInd[2]] < 8){//8
						if(_gen_lPt[genInd[0]] < 25) continue;
						if((_gen_lPt[genInd[0]] < 30) && (_gen_lPt[genInd[1]] < 15) ) continue;
					} else{
						if( (_gen_lPt[genInd[0]] < 23) && (_gen_lPt[genInd[1]] < 15)) continue;
					}
				}
			} else if(tril_flavorComb(genInd, _gen_flavors, 3) == 2){//emm
				if(_flavors[genInd[2]] == 1){
					if(_gen_lPt[genInd[2]] < 9){
						if(_gen_lPt[genInd[0]] < 23) continue; //23
					}
				}
			}


			bool etaFail = false;
			for(unsigned l = 0; l < 3; ++l){
				if(fabs(_gen_lEta[genInd[l]]) > (2.5 - 0.1*_gen_flavors[genInd[l]]) ){
					etaFail = true;
					break;
				}
			}
			if(etaFail) continue;
			*/
			lowMeff[1]->Fill(sam +0.1,_weight);
			highMeff[1]->Fill(sam + 0.1,_weight);


			//Apply HNL SELECTION
			//Baseline event selection:
			if(_nL < 3) continue;
			if(!baseline(true, true, false, false)) continue;
			if(nBJets(true, false, 0) != 0) continue;
			//uncleaned bveto;	
			//Categorize according to the number of leptons and flavors
			unsigned* ind = new unsigned[_nL];
			unsigned lCount = lepOrder(ind, 3, true, true);	
			//unsigned lCount = lepOrder(ind, 3, true, false);
			if(lCount != 3) continue; //Veto 4th FO lepton considering signal model!
		
			//lowMeff[1]->Fill(sam +0.1,1);
			//highMeff[1]->Fill(sam + 0.1,1);


			Histos[0][sam]->Fill(_isolation[ind[0]],1); 	
			Histos[1][sam]->Fill(_isolation[ind[1]],1); 			
			Histos[2][sam]->Fill(_isolation[ind[2]],1); 	
			//Require 3 leptons to be tight in data and MC, and determine nonPrompt bkg in data
			unsigned nTight = tightCount(ind, lCount);
			bool tightFail = nTight < 3;
	
			//Calculate conePt for every lepton
			double* conePt = new double[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				conePt[l] = _lPt[ind[l]]*std::max(1., 1 + (_isolation[ind[l]] - 0.1));
			}
			//require 3 tight leptons
			if(tightFail) continue;
			
			if(conePt[0] < 15) continue;
			if(conePt[1] < 10) continue;
			if(conePt[2] < 10 - 5*(_flavors[ind[2]])) continue;	//5(10) GeV cut on trailing muon(electron)
			if(!ptCuts_hnl(ind,lCount)) continue;
			
			//Apply triggers to data events;
			//if(effsam == 0){	
			/*
				bool trigPass[4];
				trigPass[0] = _lowM_trigger_eee || _lowM_trigger_all;
				trigPass[1] = _lowM_trigger_mee || _lowM_trigger_all;
				trigPass[2] = _lowM_trigger_mme || _lowM_trigger_all;
				trigPass[3] = _lowM_trigger_mmm || _lowM_trigger_all;
				if(!trigPass[tril_flavorComb(_flavors, ind, lCount)]) continue;
			*/
			//}						
			//determine search category
			unsigned cat = hnl::cat(ind, _flavors, _charges, lCount, conePt[0]);
			if(cat == 999) continue;	
			if(sam == 18) ++displacedCounter;
			//if(cat < 4){
				lowMeff[0]->Fill(sam + 0.1,_weight);
			//}
			//} else{
			//	if(conePt[2] < 10.) continue;
			//	highMeff[0]->Fill(sam + 0.1 ,_weight);
			//}				
    	}
	}
	/*
	cout << "Fraction of events with 3 prompt gen leptons = " <<  numPrompt/denom << endl;
	cout << "Fraction of events with 3 hard gen leptons = " << numHard/denom << endl;
	

	for(unsigned sam = 0; sam < nSamples; ++sam){
		cout << "Fraction of m = " << names[sam] << " samples having gen taus is : " << tauCounter[sam]/totalCounter[sam] << endl;
	}
	cout << "5 GeV displaced number of entries = " << lowMeff[0]->GetEntries() << endl;
	*/
	lowMeff[0]->Divide(lowMeff[1]);
//	highMeff[0]->Divide(highMeff[1]);	


	//plotRoc(lowMeff[0], highMeff[0], "hnl_sigeff"  + extra, names, true, "high mass selection efficiency", "low mass selection efficiency");
	
	//make TGraph as a function of mass
	unsigned nSigLow = 13; //9
	double yields[nSigLow];
	double lowMpoints[nSigLow];
	double lowMerrs[nSigLow];
	for(unsigned sam = 0; sam < nSigLow; ++ sam){ //9
		yields[sam] = lowMeff[0]->GetBinContent(sam + 1);
		lowMpoints[sam] = mass[sam];
		lowMerrs[sam] = lowMeff[0]->GetBinError(sam + 1);
	}
	cout << "1 GeV efficiency = " << yields[0] << endl;
	cout << "2 GeV efficiency = " << yields[1] << endl;	

	cout << "5 GeV prompt efficiency = " << lowMeff[0]->GetBinContent(3) << " pm " << lowMeff[0]->GetBinError(3) << endl;
	cout << "5 GeV displaced efficiency = " << lowMeff[0]->GetBinContent(19) << " pm " << lowMeff[0]->GetBinError(19) << endl;
	cout << "5 GeV efficiency ratio = " << lowMeff[0]->GetBinContent(19)/lowMeff[0]->GetBinContent(3) << endl;
	cout << "5 GeV displaced nEntries = " << displacedCounter << std::endl;
	/*
	TGraphErrors* lowMgraph = new TGraphErrors(nSigLow, lowMpoints, yields, 0, lowMerrs);
	TCanvas* c = new TCanvas("canv","canv",500,500);
	c->SetTopMargin(0.08);
	c->SetLeftMargin(0.16);
	lowMgraph->GetYaxis()->SetTitleOffset(1.2);
    //h->GetXaxis()->SetTitleOffset(1.4);
	lowMgraph->GetXaxis()->SetTitle("m_{N} (GeV)");
	lowMgraph->GetYaxis()->SetTitle("low mass selection efficiency");
	lowMgraph->GetYaxis()->SetRangeUser(0, 0.75);
	lowMgraph->Draw("CAPE");
	//drawLumi(c, "Simulation", false);
    c->SaveAs("plots/efficiency/hnl_lowMsigeff" + extra + ".pdf");
	c->SaveAs("plots/efficiency/hnl_lowMsigeff" + extra + ".png");

	double highMyields[18];
	double highMpoints[18];
	double highMerrs[18];
	for(unsigned sam = 0; sam < nSamples; ++ sam){
		highMyields[sam] = highMeff[0]->GetBinContent(sam + 1);
		highMpoints[sam] = mass[sam];
		highMerrs[sam] = highMeff[0]->GetBinError(sam + 1);
	}
	TGraphErrors* highMgraph = new TGraphErrors(18, highMpoints, highMyields, 0, highMerrs);
	TCanvas* c1 = new TCanvas("canv1","canv1",500,500);
	c->SetTopMargin(0.08);
	c->SetLeftMargin(0.16);
	highMgraph->GetYaxis()->SetTitleOffset(1.2);
    //h->GetXaxis()->SetTitleOffset(1.4);
	highMgraph->GetXaxis()->SetTitle("m_{N} (GeV)");
	highMgraph->GetYaxis()->SetTitle("high mass selection efficiency");
	highMgraph->GetYaxis()->SetRangeUser(0, 1);
	highMgraph->Draw("CAPE");
	//drawLumi(c1, "Simulation", false);
    c1->SaveAs("plots/efficiency/hnl_highMsigeff" + extra + ".pdf");
	c1->SaveAs("plots/efficiency/hnl_highMsigeff" + extra + ".png");
	*/
	/*
	std::vector<TH1D*> distVec[nSamples];
	std::vector<TString> distnamesIso = {"iso leading", "iso subleading", "iso trailing"};
	std::vector<TString> distnamesPt = {"leading gen P_{T}", "subleading gen P_{T}", "trailing gen P_{T}" };
	const TString samNames[nSamples] = {"m1", "m2", "m5", "m10", "m20", "m30", "m40", "m50", "m60", "m80", "m100", "m130", "m150", "m200", "m400", "m600", "m800", "m1000"};
	for(unsigned sam = 0; sam < nSamples; ++sam){
		for(unsigned i = 0; i < 3; ++i){
			distVec[sam].push_back(Histos[i][sam]);
		}
		plotHist(distVec[sam], distnamesIso, "efficiency/signalIso" + samNames[sam] + extra);
	} 
	
	for(unsigned sam = 0; sam < nSamples; ++sam){
		distVec[sam].clear();
		for(unsigned i = 3; i < 6; ++i){
			distVec[sam].push_back(Histos[i][sam]);
		}
		plotHist(distVec[sam], distnamesPt, "ptSpectrum/genPt" + samNames[sam] + extra);
	}  	
	*/
}



int main(int argc, char* argv[]){
	TApplication* rootapp = new TApplication("example",&argc, argv);
	trilTree testtree;
	testtree.Loop();
	rootapp->Run();
    return 0;
}


