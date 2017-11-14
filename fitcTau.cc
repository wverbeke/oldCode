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
#include <map>
#include <memory>

//include code to calculate btag SF
#include "../bTag/BTagCalibrationStandalone.h"
//include other parts of the code
#include "MultilepSUSYfunc.h"
#include "tdrstyle.h"
#include "plotCode.h"
#include "trilTree.h"
#include "hnlTools.h"
//#include "drawLumi.h"

void trilTree::fitcTau(const std::string& sample){
	//Read Trees from ROOT files
	//TFile* sampleFile;
	std::shared_ptr<TFile> sampleFile = std::make_shared<TFile>("../unskimmedSignal/"+(const TString&) sample,"read");
    sampleFile->cd("FakeElectrons");
    TTree* inputTree = (TTree*) (sampleFile->Get("FakeElectrons/fakeTree"));
	Init(inputTree, false, true);
	
	TH1D cTau;
	cTau = TH1D("ctau", "", 400, 0, 200);
	cTau.Sumw2();
	
   	Long64_t nEntries = inputTree->GetEntries();
	unsigned counter = 0;
    //std::cout<<"Entries in "<< sample <<" "<<nEntries<<std::endl;
    for(Long64_t it = 0; it < nEntries; ++it){
		if (it%10000 == 0) std::cout <<'.'<< std::flush;
       	inputTree->GetEntry(it);
		/*
		//Correct weight for lifetime
		cutBased();
		//Order particles by gen lepton Pt and then store them
		if(_gen_nL < 3) continue;
		//Apply HNL SELECTION
		//Baseline event selection:
		if(_nL < 3) continue;
		if(!baseline(true, true, false, false)) continue;
		if(nBJets(true, false, 0) != 0) continue;
		std::shared_ptr<unsigned> ind(new unsigned[_nL], [] (unsigned* ptr) -> void {delete[] ptr;}); 
		const unsigned lCount = lepOrder(ind.get(), 3, true, true);	
		if(lCount != 3) continue; //Veto 4th FO lepton considering signal model!
		//Require 3 leptons to be tight in data and MC, and determine nonPrompt bkg in data
		unsigned nTight = tightCount(ind.get(), lCount);
		bool tightFail = nTight < 3;
		double conePt[lCount];
		for(unsigned l = 0; l < lCount; ++l){
			conePt[l] = _lPt[ind.get()[l]]*std::max(1., 1 + (_isolation[ind.get()[l]] - 0.1));
		}
		//require 3 tight leptons
		if(tightFail) continue;
		//HNL pt Cuts
		if(!ptCuts_hnl(ind.get(),lCount)) continue;
		//veto events outside the HNL selection categories (+++, ---)
		unsigned cat = hnl::cat(ind.get(), _flavors, _charges, lCount, conePt[0]);
		if(cat == 999) continue;
		*/	
		cTau.Fill(_ctau, _weight);
	}
	TF1 *f1 = new TF1("fit","[1]/[0]*exp(-x/[0])",0,200);
	f1->SetParameter(0,2.);
	cTau.Fit("fit");
	cTau.Fit("fit");
	cTau.GetXaxis()->SetTitle("c#tau (mm)");
	cTau.GetYaxis()->SetTitle("Events");
	std::shared_ptr<TCanvas> c = std::make_shared<TCanvas>("canv","canv",500,500);
	c->SetLogy();
	c->SetTopMargin(0.08);
	c->SetLeftMargin(0.18);
	c->SetBottomMargin(0.14);
	cTau.Draw();
	c->SaveAs("plots/ctauFit.pdf");
	c->Close();
}

void trilTree::comparecTau(const std::vector< std::string>& samples){
	setTDRStyle();
	gStyle->SetOptFit(0);
	//Read Trees from ROOT files
	TH1D* cTau = new TH1D[samples.size()];	
	const double ctau_0 = 11.3714830365;
	const double ctau_1 = 0.92160604666;

	const unsigned nSamples = samples.size();
	const unsigned nDist = 9;
	TH1D* histos[nDist][nSamples];
	const TString distNames[nDist] = {"dxy_displaced", "dz_displaced", "sip3d_displaced", "mt", "met", "minMos", "leadingPt", "M3l", "SearchR"};
	const TString xAxes[nDist] = {"d_{xy} (cm)", "d_{z} (cm)", "#sigma(IP_{3D})/IP_{3D}", "M_{T} (GeV)", "MET (GeV)", "M^{min}_{2l OS} (GeV)", "P_{T}^{lead} (GeV)", "M_{3l} (GeV)", "Search Region"};
	const double histMin[nDist] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	const double histMax[nDist] = {0.3, 0.3, 30, 300, 300, 10, 150, 150, 12};
	unsigned nBins[nDist];
	for(unsigned dist = 0; dist < nDist; ++dist) nBins[dist] = 50;
	nBins[8] = 12;
	nBins[5] = 20;
	for(unsigned dist = 0; dist < nDist; ++dist){
		for(unsigned sam = 0; sam < nSamples; ++sam){
			histos[dist][sam] = new TH1D(distNames[dist] + samples[sam], distNames[dist] + samples[sam] + ";" + xAxes[dist] + "; normalized Events", nBins[dist], histMin[dist], histMax[dist]);
		}
	}
	readSF(true);
	//TFile* sampleFile;
	for(unsigned sam = 0; sam < samples.size(); ++sam){
		//initialize cTau histogram for given sample
		//cTau[sam] = TH1D("ctau" + (const TString&) samples[sam], "", 200, 0, 30);
		cTau[sam] = TH1D("ctau" + (const TString&) samples[sam], "", 200, 0, 6);
		cTau[sam].Sumw2();
		//read data from file
		std::shared_ptr<TFile> sampleFile = std::make_shared<TFile>("../unskimmedSignal/"+(const TString&) samples[sam],"read");
		sampleFile->cd("FakeElectrons");
		TTree* inputTree = (TTree*) (sampleFile->Get("FakeElectrons/fakeTree"));
		Init(inputTree, false, true);
		
		Long64_t nEntries = inputTree->GetEntries();
		unsigned counter = 0;
		TH1D averageSF = TH1D((const TString&) "averageSF" + samples[sam], (const TString&) "averageSF" + samples[sam] + "; ; normalized Events", 1, 0.,1.);
		//std::cout<<"Entries in "<< sample <<" "<<nEntries<<std::endl;
		for(Long64_t it = 0; it < nEntries; ++it){
			if (it%10000 == 0) std::cout <<'.'<< std::flush;
		   	inputTree->GetEntry(it);
			if(sam == 0) _weight*= ( ctau_0)/(ctau_1)*exp(- _ctau/ctau_1 + _ctau/ ctau_0);
			cTau[sam].Fill(_ctau, _weight);
			if(sam == 1){
				//std::cout << "weight = " << _weight << std::endl;
				//std::cout << "ctau = " << _ctau << std::endl;
			}

			//find leading lepton Pt
			/*
			if(_nL != 3) continue;
			unsigned maxI = 0;
			for(unsigned l = 1; l < _nL; ++l){
				if(_lPt[l] > _lPt[maxI]) maxI = l;
			}
			//compute minMos, MT and M3l
			//make lorentz vector for every lepton
			TLorentzVector lepV[(const unsigned) _nL];
			for(unsigned l = 0; l < _nL; ++l){
				lepV[l].SetPtEtaPhiE(_lPt[l], _lEta[l], _lPhi[l], _lE[l]);
			}
			//min(Mos)
			double minMos = 0;
			unsigned minI[2] = {99, 99};
			for(unsigned l = 0; l < _nL -1 ; ++l){
				for(unsigned k = l + 1; k < _nL; ++k){
					if(_charges[l] != _charges[k]){
						if( (lepV[l] + lepV[k]).M() < minMos  || minMos == 0){
							minMos = (lepV[l] + lepV[k]).M();
							minI[0] = l;
							minI[1] = k;
						}
					}
				}
			}
			//Calculate MET vector
			TLorentzVector METvec;
			METvec.SetPtEtaPhiE(_met, 0, _met_phi, _met);
			//find lepton for MT calculation
			unsigned lw_min;
			for(unsigned l = 0; l < _nL; ++l){
				if(l != minI[0] && l != minI[1]){
					lw_min = l;
				}
			}
			double mt_min = transmass(lepV[lw_min], METvec);
			double m3l = (lepV[0] + lepV[1] + lepV[2]).M();
			//find indices for displaced leptons
			unsigned minDisp = 0;
			for(unsigned l = 1; l < _nL; ++l){
				if(_3dIPsig[l] < _3dIPsig[minDisp]) minDisp = l;
			}
			for(unsigned l = 0; l < _nL; ++l){
				//if(_mompdg[l] == 9900012){
				if(l != minDisp){
					histos[0][sam]->Fill(_ipPV[l], _weight);
					histos[1][sam]->Fill(_ipZPV[l], _weight);
					histos[2][sam]->Fill(_3dIPsig[l], _weight);
				}
			}
			histos[3][sam]->Fill(mt_min, _weight);
			histos[4][sam]->Fill(_met, _weight);
			histos[5][sam]->Fill(minMos, _weight);
			histos[6][sam]->Fill(_lPt[maxI], _weight);
			histos[7][sam]->Fill(m3l, _weight);
			*/
			/*
			unsigned cat = hnl::cat(ind, _flavors, _charges, lCount, conePt[0]);
			if(cat == 999){
				continue;	
			}
			if(cat != 1 && cat != 3) continue;
			unsigned searchR = 4*(cat == 3) + hnl::sr(mt_min, minMos, lepSyst.M(), cat); 
			histos[8][sam]->Fill(searchR, _weight);
			*/
			
			for(unsigned l = 0; l < _nL; ++l){
				_islooseCut[l] = _islooseCut[l] && fabs(_ipPV[l]) < 0.05 && fabs(_ipZPV[l]) < 0.1;
			}
			cutBased();
			//Order particles by gen lepton Pt and then store them
			//if(_gen_nL < 3) continue;
		
			//Apply HNL SELECTION
			//Baseline event selection:
			if(_nL < 3) continue;
			if(!baseline(true, true, false, false)) continue;
			if(nBJets(true, false, 0) != 0) continue;
			//uncleaned bveto;	
			//Categorize according to the number of leptons and flavors
			//unsigned* ind = new unsigned[_nL];	
			const unsigned nLtemp = _nL;
			unsigned ind[nLtemp];
			//std::shared_ptr<unsigned> ind(new unsigned[_nL], [] (unsigned* ptr) -> void {delete[] ptr;}); 
			const unsigned lCount = lepOrder(ind, 3, true, true);	
			if(lCount != 3) continue; //Veto 4th FO lepton considering signal model!
			//Require 3 leptons to be tight in data and MC, and determine nonPrompt bkg in data
			unsigned nTight = tightCount(ind, lCount);
			bool tightFail = nTight < 3;

			//Calculate conePt for every lepton
			//double* conePt = new double[lCount];
			double conePt[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				conePt[l] = _lPt[ind[l]]*std::max(1., 1 + (_isolation[ind[l]] - 0.1));
			}
			//require 3 tight leptons
			if(tightFail) continue;
			//HNL pt Cuts
			if(!ptCuts_hnl(ind,lCount)) continue;
			//veto events outside the HNL selection categories (+++, ---)
			unsigned cat = hnl::cat(ind, _flavors, _charges, lCount, conePt[0]);
			if(cat == 999) continue;

			//compute minMos, MT and M3l
			//make lorentz vector for every lepton
			TLorentzVector lepV[lCount];
			for(unsigned l = 0; l < lCount; ++l){
				lepV[l].SetPtEtaPhiE(_lPt[ind[l]], _lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]);
			}
			//min(Mos)
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
			//Calculate MET vector
			TLorentzVector METvec;
			METvec.SetPtEtaPhiE(_met, 0, _met_phi, _met);
			//find lepton for MT calculation
			unsigned lw_min;
			for(unsigned l = 0; l < lCount; ++l){
				if(l != minI[0] && l != minI[1]){
					lw_min = l;
				}
			}
			double mt_min = transmass(lepV[lw_min], METvec);
			double m3l = (lepV[0] + lepV[1] + lepV[2]).M();



			//find indices for displaced leptons
			for(unsigned l = 0; l < lCount; ++l){
				//if(_mompdg[ind[l]] == 9900012){
					histos[0][sam]->Fill(_ipPV[ind[l]], _weight);
					histos[1][sam]->Fill(_ipZPV[ind[l]], _weight);
					histos[2][sam]->Fill(_3dIPsig[ind[l]], _weight);
				//}
			}
			
			histos[3][sam]->Fill(mt_min, _weight);
			histos[4][sam]->Fill(_met, _weight);
			histos[5][sam]->Fill(minMos, _weight);
			histos[6][sam]->Fill(conePt[0], _weight);
			histos[7][sam]->Fill(m3l, _weight);

			if(cat != 1 && cat != 3) continue;
			averageSF.Fill(0.5, bTagSF(true, 0));
			++counter;
			unsigned searchR = 4*(cat == 3) + hnl::sr(mt_min, minMos, m3l, cat); 
			histos[8][sam]->Fill(searchR, _weight);
		}
		std::cout << "Average b-tag SF = " << averageSF.GetBinContent(1)/averageSF.GetEntries() << " +- " <<  averageSF.GetBinError(1)/averageSF.GetEntries() << std::endl;
	}
	//Fit ctau to each sample 
	TF1** f = new TF1*[samples.size()];
	//Normalize histograms
	for(unsigned sam = 0; sam < samples.size(); ++sam){
		cTau[sam].Scale(1/cTau[sam].GetSumOfWeights());
	}
	const Color_t colors[8] = {kBlue, kRed};//kGreen - 7
	for(unsigned sam = 0; sam < samples.size(); ++sam){
		//f[sam] = TF1("fit" + (const TString&) samples[sam], "[1]/[0]*exp(-x/[0])",0,30);
		//f[sam] = new TF1("fit" + (const TString&) std::to_string(sam), "[1]/[0]*exp(-x/[0])",0,30);
		f[sam] = new TF1("fit" + (const TString&) std::to_string(sam), "[1]/[0]*exp(-x/[0])",0,7);
		f[sam]->SetParameter(0,11.);	
		f[sam]->SetLineWidth(3);
		if(sam == 1){
			f[1]->SetLineColor(kBlue);
			f[1]->SetMarkerColor(kBlack);
		} else if(sam == 0){
			f[0]->SetLineColor(kRed);
			f[0]->SetMarkerColor(kBlue);
		}
		cTau[sam].Fit("fit" + (const TString&) std::to_string(sam), "M");
		cTau[sam].Fit("fit"+ (const TString&) std::to_string(sam), "M");
	}
	cTau[1].SetMarkerColor(kBlack);
	const unsigned markerStyles[2] = {9, 22};
	const double markerSize[2] = {0.5, 0.7};
	for(unsigned sam = 0; sam < samples.size(); ++sam){
		cTau[sam].SetMarkerSize( markerSize[sam]);
		cTau[sam].SetMarkerStyle(markerStyles[sam]);
	
	}
	//Plot samples with fit
	std::shared_ptr<TCanvas> c = std::make_shared<TCanvas>("canv","canv",500,500);
	c->SetLogy();
	c->SetTopMargin(0.08);
	c->SetLeftMargin(0.18);
	c->SetBottomMargin(0.14);
	TLegend* legend = new TLegend(0.60,0.7,0.92,0.90,NULL,"brNDC");
	legend->AddEntry(&cTau[0], "V_{#muN} = 0.00316, c#tau = 11mm", "p");
	legend->AddEntry(&cTau[1], "V_{#muN} = 0.0111, c#tau = 0.9mm", "p");
	//legend->AddEntry(&cTau[1], "V_{#muN} = 0.0111, prompt", "p");
	cTau[0].GetXaxis()->SetTitle("c#tau (mm)");
	cTau[0].GetYaxis()->SetTitle("normalized Events");
	cTau[0].Draw("");
	//f[0]->Draw("same");
	for(unsigned sam = 1; sam < samples.size(); ++sam){
		cTau[sam].Draw("same");
		//f[sam]->Draw("same");
	}
	legend->Draw("same");
	c->SaveAs("plots/ctauComparison_reweighted_CHECK.pdf");
	c->Close();

	//plot kinematic comparison of displaced samples
	for(unsigned dist = 0; dist < nDist; ++dist){
		plotHistRatio(histos[dist][0], histos[dist][1], "1 mm", "prompt", "displacedComparison/displacedComp_" + distNames[dist] + "_Prompt_Selected_5GeV");// + "_Selected");

	}
}


void trilTree::effVScTau(const std::string& sample){
	//Read Trees from ROOT files
	//TFile* sampleFile;
	std::shared_ptr<TFile> sampleFile = std::make_shared<TFile>("../unskimmedSignal/"+(const TString&) sample,"read");
    sampleFile->cd("FakeElectrons");
    TTree* inputTree = (TTree*) (sampleFile->Get("FakeElectrons/fakeTree"));
	Init(inputTree, false, true);
	//efficiency histograms
	TH1D eff[2];
	eff[0] = TH1D("efficiency num left", "", 15, 0, 2);
	eff[1] = TH1D("efficiency denom left", "", 15, 0, 2);
	eff[0].Sumw2();
	eff[1].Sumw2(); 	
	//loop over events
	
   	Long64_t nEntries = inputTree->GetEntries();
	unsigned counter = 0;
    for(Long64_t it = 0; it < nEntries; ++it){
       	inputTree->GetEntry(it);
		//Correct weight for lifetime
		cutBased();
		//Order particles by gen lepton Pt and then store them
		//if(_gen_nL < 3) continue;
		eff[1].Fill(_ctau,_weight);
		//Apply HNL SELECTION
		//Baseline event selection:
		if(_nL < 3) continue;
		if(!baseline(true, true, false, false)) continue;
		if(nBJets(true, false, 0) != 0) continue;
		//uncleaned bveto;	
		//Categorize according to the number of leptons and flavors
		//unsigned* ind = new unsigned[_nL];	
		const unsigned nLtemp = _nL;
		//unsigned ind[nLtemp];
		std::shared_ptr<unsigned> ind(new unsigned[_nL], [] (unsigned* ptr) -> void {delete[] ptr;}); 
		const unsigned lCount = lepOrder(ind.get(), 3, true, true);	
		if(lCount != 3) continue; //Veto 4th FO lepton considering signal model!
		//Require 3 leptons to be tight in data and MC, and determine nonPrompt bkg in data
		unsigned nTight = tightCount(ind.get(), lCount);
		bool tightFail = nTight < 3;

		//Calculate conePt for every lepton
		//double* conePt = new double[lCount];
		double conePt[lCount];
		for(unsigned l = 0; l < lCount; ++l){
			conePt[l] = _lPt[ind.get()[l]]*std::max(1., 1 + (_isolation[ind.get()[l]] - 0.1));
		}
		//require 3 tight leptons
		if(tightFail) continue;
		//if(tril_flavorComb(ind.get(), _flavors,  lCount) != 1) continue;
		//HNL pt Cuts
		if(!ptCuts_hnl(ind.get(),lCount)) continue;
		//veto events outside the HNL selection categories (+++, ---)
		unsigned cat = hnl::cat(ind.get(), _flavors, _charges, lCount, conePt[0]);
		if(cat == 999) continue;
		eff[0].Fill(_ctau,_weight);
	}
	eff[0].Divide(&eff[1]);
	eff[0].GetXaxis()->SetTitle("c#tau (mm)");
	eff[0].GetYaxis()->SetTitle("Efficiency");
	std::shared_ptr<TCanvas> c = std::make_shared<TCanvas>("canv","canv",500,500);
	c->SetLogy();
	c->SetTopMargin(0.08);
	c->SetLeftMargin(0.18);
	c->SetBottomMargin(0.14);
	eff[0].Draw();
	c->SaveAs("plots/effVsCTau_1GeV.pdf");
	c->Close();	

	setTDRStyle();
	plotHist(&eff[0], "effVsCTau_1GeV_alt.pdf");	
}


int main(int argc, char* argv[]){
	trilTree testtree;
	//testtree.fitcTau("HeavyNeutrino_trilepton_M-1_V-0.11505_mu_displaced.root");
	testtree.comparecTau({"HeavyNeutrino_trilepton_M-5_V-0.00316_mu_displaced_unskimmedLep.root", "HeavyNeutrino_trilepton_M-5_V-0.0111_mu_prompt_unskimmedLep.root"});
	//testtree.comparecTau({"HeavyNeutrino_trilepton_M-5_V-0.0111_mu_displaced_unskimmedLep.root", "HeavyNeutrino_trilepton_M-5_V-0.0111_mu_prompt_unskimmedLep.root"});
	//testtree.comparecTau({"HeavyNeutrino_trilepton_M-12_V-0.0011934501_mu_displaced.root", "HeavyNeutrino_trilepton_M-12_V-0.0011934501_mu_prompt.root"});
	//testtree.effVScTau("HeavyNeutrino_trilepton_M-1_V-0.11505_mu_displaced.root");
	//testtree.comparecTau({"HeavyNeutrino_M1_mu.root", "HeavyNeutrino_trilepton_M-1_V-0.11505_mu_displaced.root"});
    return 0;
}


