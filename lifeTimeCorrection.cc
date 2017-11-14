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
#include "TLine.h"

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

std::pair<double, double> sampleTau(const std::string& sample){
	static bool first = true;
	static std::map<std::string, std::pair<double,double> > tauMap; //return {tau_prompt, tau_displaced} 
	if(first){
		tauMap["HeavyNeutrino_M5_V00p00336_mu.root"] = {5.68574838114612042*2, 0.82046213139960145};
		tauMap["HeavyNeutrino_trilepton_M-5_V-0.00316_mu.root"] = {10.0584706495, 11.3714830365};//{5.68574838114612042*2, 5.68574838114612042*2};
		tauMap["HeavyNeutrino_trilepton_M-5_V-0.0111_mu.root"] = {10.0584706495, 0.92160604666};//{5.68574838114612042*2, 5.68574838114612042*2};
		tauMap["HeavyNeutrino_trilepton_M-5_V-0.00316_e.root"] = {11.0529212344, 11.3714830365}; //{5.6857483814612042*2, 5.6857483814612042*2};
		tauMap["HeavyNeutrino_M5_V0p00316_e.root"]= {11.0529212344, 11.3714830365};
		tauMap["HeavyNeutrino_trilepton_M-6_V-0.00316_mu_displaced.root"] = {4.30719500732,	4.56612878826};//{4.56612878826, 4.56612878826};
		tauMap["HeavyNeutrino_trilepton_M-6_V-0.00316_e_displaced.root"] = {4.7642694571, 4.56612878826};//{4.56612878826, 4.56612878826};
		tauMap["HeavyNeutrino_trilepton_M-7_V-0.00316_mu_displaced.root"] = {1.88098144141,	2.32787};//2.11179046242};//{2.11179046242, 2.11179046242};
		tauMap["HeavyNeutrino_trilepton_M-7_V-0.00316_e_displaced.root"] = {2.0068420451, 2.11179046242};//{2.11179046242, 2.11179046242};
	
		tauMap["HeavyNeutrino_trilepton_M-1_V-0.5_mu_displaced.root"] = {7970.70862994, 1.41978088057};
 		tauMap["HeavyNeutrino_trilepton_M-1_V-0.5_e_displaced.root"] = {8142.12191536, 1.41978088057};
		
		tauMap["HeavyNeutrino_trilepton_M-1_V-0.11505_mu_displaced.root"] = {7970.70862994, 26.815644318875442};

		tauMap["HeavyNeutrino_trilepton_M-2_V-0.07071_mu_displaced.root"] = {621.567395371, 2.22408538668};
		tauMap["HeavyNeutrino_trilepton_M-2_V-0.07071_e_displaced.root"] = {761.168401415, 2.22408538668};

		tauMap["HeavyNeutrino_trilepton_M-3_V-0.0245_e_displaced.root"] = {111.702504399, 2.43519833788};
		tauMap["HeavyNeutrino_trilepton_M-3_V-0.0245_mu_displaced.root"] = {132.226817762, 2.43519833788};

		tauMap["HeavyNeutrino_trilepton_M-4_V-0.01225_mu_displaced.root"] = {27.4400433003, 2.30714878353};
		tauMap["HeavyNeutrino_trilepton_M-4_V-0.01225_e_displaced.root"] = {30.989940326, 2.30714878353};

	
		tauMap["HeavyNeutrino_trilepton_M-8_V-0.00328876176_e_displaced.root"] = {1.04213318774 ,9.97392e-01};
		tauMap["HeavyNeutrino_trilepton_M-8_V-0.00328876176_mu_displaced.root"] = {0.99185445394 ,9.99415e-01};

		tauMap["HeavyNeutrino_trilepton_M-9_V-0.00244991552_e_displaced.root"] = {0.63586194257, 9.97327e-01};
		tauMap["HeavyNeutrino_trilepton_M-9_V-0.00244991552_mu_displaced.root"] = {0.58818212698, 9.91448e-01};

		tauMap["HeavyNeutrino_trilepton_M-8_V-0.00328876176_mu_displaced_noIPCuts.root"] = {0.58818212698 ,9.99415e-01};
		tauMap["HeavyNeutrino_trilepton_M-9_V-0.00244991552_mu_displaced_noIPCuts.root"] = {0.58818212698, 9.91448e-01};

		tauMap["HeavyNeutrino_trilepton_M-10_V-0.00188259708_e_displaced.root"] = {0.32860775905, 9.95993e-01};
		tauMap["HeavyNeutrino_trilepton_M-10_V-0.00188259708_mu_displaced.root"] = {0.33887334453, 9.97228e-01};

		tauMap["HeavyNeutrino_trilepton_M-11_V-0.00148345941_e_displaced.root"] = {0.20123589399, 9.89449e-01};
		tauMap["HeavyNeutrino_trilepton_M-11_V-0.00148345941_mu_displaced.root"] = {0.20102824119, 9.88428e-01};

		tauMap["HeavyNeutrino_trilepton_M-12_V-0.0011934501_e_displaced.root"] = {0.13563783783, 9.93217e-01};
		tauMap["HeavyNeutrino_trilepton_M-12_V-0.0011934501_mu_displaced.root"] = {0.13014938142 ,9.88719e-01}; 
		first = false;
	}
	return tauMap[sample];
}



std::pair<double, double> trilTree::findLifeTime(const double left, const double right, const std::string& sample, const bool reweighing){
	//Read Trees from ROOT files
	//TFile* sampleFile;
	std::shared_ptr<TFile> sampleFile = std::make_shared<TFile>("../unskimmedSignal/"+(const TString&) sample,"read");
    sampleFile->cd("FakeElectrons");
    TTree* inputTree = (TTree*) (sampleFile->Get("FakeElectrons/fakeTree"));
	Init(inputTree, false, true);
	
	//mean ctau for sample
	std::pair<double, double> tau;
	if(reweighing) tau = sampleTau(sample);

	TH1D eff_left[2];
	TH1D eff_right[2];
	eff_left[0] = TH1D("efficiency num left", "", 1, 0, 1);
	eff_left[1] = TH1D("efficiency denom left", "", 1, 0, 1);
	eff_left[0].Sumw2();
	eff_left[1].Sumw2();

	eff_right[0] = TH1D("efficiency num right", "", 1, 0, 1);
	eff_right[1] = TH1D("efficiency denom right", "", 1, 0, 1);
	eff_right[0].Sumw2();
	eff_right[1].Sumw2();
	
   	Long64_t nEntries = inputTree->GetEntries();
	unsigned counter = 0;
    //std::cout<<"Entries in "<< sample <<" "<<nEntries<<std::endl;
    for(Long64_t it = 0; it < nEntries; ++it){
       	inputTree->GetEntry(it);
       	//if (it%10000 == 0) std::cout <<'.'<< std::flush;
		//Correct weight for lifetime
		cutBased();
		//Order particles by gen lepton Pt and then store them
		//if(_gen_nL < 3) continue;
		
		double leftWeight = 1, rightWeight = 1;
		if(reweighing){
			//leftWeight = left*left*exp(-_ctau/tau.second + _ctau/(left*sqrt(tau.second*tau.first) ));
			//rightWeight =  right*right*exp(-_ctau/tau.second + _ctau/(right*sqrt(tau.second*tau.first) ));
			leftWeight = left*left*exp(_ctau/tau.second - _ctau*left/sqrt(tau.second*tau.first) );
			//std::cout << " exp = " << exp(_ctau/tau.second - _ctau*left/sqrt(tau.second*tau.first) ) << std::endl;
			rightWeight =  right*right*exp(_ctau/tau.second - _ctau*right/sqrt(tau.second*tau.first) );
		}
		eff_left[1].Fill(0.1,_weight);
		eff_right[1].Fill(0.1,_weight);
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
		/*
		for(unsigned l = 0; l < lCount; ++l){
			if(fabs(_ipPV[ind.get()[l]]) > 0.05) std::cout << "dxy cut fail" << std::endl;
			if(fabs(_ipZPV[ind.get()[l]]) > 0.1) std::cout << "dz cut fail" << std::endl;
			if(fabs(_3dIPsig[ind.get()[l]]) > 4) std::cout << "3DIPSig cut fail" << std::endl;
		}
		*/
		//HNL pt Cuts
		if(!ptCuts_hnl(ind.get(),lCount)) continue;
		//veto events outside the HNL selection categories (+++, ---)
		unsigned cat = hnl::cat(ind.get(), _flavors, _charges, lCount, conePt[0]);
		if(cat == 999) continue;
		/*
		bool trigPass[4];
		trigPass[0] = _lowM_trigger_eee || _lowM_trigger_all;
		trigPass[1] = _lowM_trigger_mee || _lowM_trigger_all;
		trigPass[2] = _lowM_trigger_mme || _lowM_trigger_all;
		trigPass[3] = _lowM_trigger_mmm || _lowM_trigger_all;
		if(!trigPass[tril_flavorComb(_flavors,  lCount)]) continue;
		*/
		++counter;
		//std::cout << "leftWeight = " << leftWeight << std::endl;
		//std::cout << "rightWeight = " << rightWeight << std::endl;
		eff_left[0].Fill(0.1,_weight*leftWeight);
		eff_right[0].Fill(0.1,_weight*rightWeight);
	}
	std::cout << "number of pass events for " << sample << " = " << counter << std::endl;
	eff_left[0].Divide(&eff_left[1]);
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
	std::cout << "Stat unc on correction = " << eff_left[0].GetBinError(1)/eff_left[0].GetBinContent(1) << std::endl;
	eff_right[0].Divide(&eff_right[1]);
	std::pair<double, double> acc = {eff_left[0].GetBinContent(1), eff_right[0].GetBinContent(1)};
	return acc;
}

int main(int argc, char* argv[]){
	trilTree testtree;
	//loop over all samples
	const std::vector< std::pair < std::string , std::string > > samples = {
		//MUONS
		/*
		{"HeavyNeutrino_M1_mu.root", "HeavyNeutrino_trilepton_M-1_V-0.5_mu_displaced.root"}//"HeavyNeutrino_trilepton_M-1_V-0.11505_mu_displaced.root"}//"HeavyNeutrino_trilepton_M-1_V-0.5_mu_displaced.root"}
		{"HeavyNeutrino_M2_mu.root", "HeavyNeutrino_trilepton_M-2_V-0.07071_mu_displaced.root"},
		{"HeavyNeutrino_trilepton_M-3_V-0.0245_mu_prompt.root", "HeavyNeutrino_trilepton_M-3_V-0.0245_mu_displaced.root"},
		{"HeavyNeutrino_trilepton_M-4_V-0.01_mu_prompt.root", "HeavyNeutrino_trilepton_M-4_V-0.01225_mu_displaced.root"},
		{"HeavyNeutrino_M5_mu.root", "HeavyNeutrino_trilepton_M-5_V-0.0111_mu.root"},//00316  0.0111
		{"HeavyNeutrino_trilepton_M-6_V-0.00316_mu_prompt.root", "HeavyNeutrino_trilepton_M-6_V-0.00316_mu_displaced.root"},
		{"HeavyNeutrino_trilepton_M-7_V-0.00316_mu_prompt.root", "HeavyNeutrino_trilepton_M-7_V-0.00316_mu_displaced.root"},
		{"HeavyNeutrino_trilepton_M-8_V-0.00328876176_mu_prompt.root", "HeavyNeutrino_trilepton_M-8_V-0.00328876176_mu_displaced.root"},
		{"HeavyNeutrino_trilepton_M-9_V-0.00244991552_mu_prompt.root", "HeavyNeutrino_trilepton_M-9_V-0.00244991552_mu_displaced.root"},
		{"HeavyNeutrino_M10_mu.root", "HeavyNeutrino_trilepton_M-10_V-0.00188259708_mu_displaced.root"},
		{"HeavyNeutrino_trilepton_M-11_V-0.00148345941_mu_prompt.root", "HeavyNeutrino_trilepton_M-11_V-0.00148345941_mu_displaced.root"},
		{"HeavyNeutrino_trilepton_M-12_V-0.0011934501_mu_prompt.root", "HeavyNeutrino_trilepton_M-12_V-0.0011934501_mu_displaced.root"},
		//ELECTRONS
		{"HeavyNeutrino_M1_e.root", "HeavyNeutrino_trilepton_M-1_V-0.5_e_displaced.root"},
		{"HeavyNeutrino_M2_e.root", "HeavyNeutrino_trilepton_M-2_V-0.07071_e_displaced.root"},
		{"HeavyNeutrino_trilepton_M-3_V-0.0245_e_prompt.root", "HeavyNeutrino_trilepton_M-3_V-0.0245_e_displaced.root"},
		{"HeavyNeutrino_trilepton_M-4_V-0.01_e_prompt.root", "HeavyNeutrino_trilepton_M-4_V-0.01225_e_displaced.root"},
		{"HeavyNeutrino_M5_e.root", "HeavyNeutrino_M5_V0p00316_e.root"},
		{"HeavyNeutrino_trilepton_M-6_V-0.00316_e_prompt.root", "HeavyNeutrino_trilepton_M-6_V-0.00316_e_displaced.root"},
		{"HeavyNeutrino_trilepton_M-7_V-0.00316_e_prompt.root", "HeavyNeutrino_trilepton_M-7_V-0.00316_e_displaced.root"},
		{"HeavyNeutrino_trilepton_M-8_V-0.00328876176_e_prompt.root", "HeavyNeutrino_trilepton_M-8_V-0.00328876176_e_displaced.root"},
		{"HeavyNeutrino_trilepton_M-9_V-0.00244991552_e_prompt.root", "HeavyNeutrino_trilepton_M-9_V-0.00244991552_e_displaced.root"},
		{"HeavyNeutrino_M10_e.root", "HeavyNeutrino_trilepton_M-10_V-0.00188259708_e_displaced.root"},
		{"HeavyNeutrino_trilepton_M-11_V-0.00148345941_e_prompt.root", "HeavyNeutrino_trilepton_M-11_V-0.00148345941_e_displaced.root"},
		{"HeavyNeutrino_trilepton_M-12_V-0.0011934501_e_prompt.root", "HeavyNeutrino_trilepton_M-12_V-0.0011934501_e_displaced.root"}
		*/
};

	for(auto samIt = samples.begin(); samIt != samples.end(); ++samIt){	
		const double promptAcc = (testtree.findLifeTime(0, 0, samIt->first, false)).first;
		const double corrRatio = (testtree.findLifeTime(0, 0, samIt->second, false)).first/promptAcc;
		std::cout << "promptAcc = "  << promptAcc << std::endl;
		std::cout << "corrRatio = "  << corrRatio << std::endl;
		//double left = 1;
		//double right = corrRatio;
		double left = 1/corrRatio;
		double right = 1;
		//TGraph ygraph = TGraph();
		std::vector<double> xPoints;
		std::vector<double> yPoints;

		std::pair<double, double> tau = sampleTau(samIt->second);
		/*
		std::pair<double, double> test = testtree.findLifeTime(left, right, samIt->second);
		xPoints.push_back(left);
		xPoints.push_back(right);
		yPoints.push_back(test.first);
		yPoints.push_back(test.second);
		double step = (left - right)/50;
		*/
		for(unsigned i = 0; i < 7 && left != right; ++i){	
			double left_old = left;
			double right_old = right;
			left = left - (left - right)/2.;
			right = right + (left - right)/2.;
			//left -= step/2;
			//right += step/2;
			std::pair<double, double> test = testtree.findLifeTime(left, right, samIt->second);
			/*
			xPoints.push_back(left);
			xPoints.push_back(right);
			yPoints.push_back(test.first);
			yPoints.push_back(test.second);
			*/
			//if(test.second > corrRatio*promptAcc*corrRatio) std::cout << "Disaster in starting conditions" << std::endl;
			int j = 0;
			//while(test.first < corrRatio*promptAcc*corrRatio){
			while(test.first < promptAcc){
				//std::cout << "test.first = " << test.first << "			corrRatio*promptAcc*corrRatio = " << corrRatio*promptAcc*corrRatio << std::endl;
				if(j == 0) right = left; //NEW for faster convergence
				left = left + (left_old - left)/(2.);
				test = testtree.findLifeTime(left, right, samIt->second);
				xPoints.push_back(left);
				yPoints.push_back(test.first);
				++j;
			}
			j = 0;
			//while(test.second > corrRatio*promptAcc*corrRatio){
			while(test.second > promptAcc){
				//std::cout << "test.second = " << test.second << "			corrRatio*promptAcc*corrRatio = " << corrRatio*promptAcc*corrRatio << std::endl;
				if(j == 0) left = right; //NEW for faster convergence
				right = right - (right - right_old)/(2.);
				test = testtree.findLifeTime(left, right, samIt->second);
				xPoints.push_back(right);
				yPoints.push_back(test.second);
				++j;
			}
			//std::cout << "left = " << (left*sqrt(tau.second * tau.first))/tau.first << " test.first = " << test.first << std::endl;
			//std::cout << "right = " << (right*sqrt(tau.second * tau.first))/tau.first << " test.second = " << test.second  << std::endl;
			//std::cout << "left = " << left*tau.first/sqrt(tau.second * tau.first) << " test.first = " << test.first << std::endl;
			//std::cout << "right = " << right*tau.first/sqrt(tau.second * tau.first) << " test.second = " << test.second  << std::endl;
		}
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		//std::cout << "Final result for " << samIt->first << "   : left = " << (left*sqrt(tau.first * tau.second))/tau.first  << " 		right = " << (right*sqrt(tau.first * tau.second))/tau.first  << std::endl;
		std::cout << "Final result for " << samIt->first << "   : left = " << left*tau.first/sqrt(tau.second * tau.first) << " 		right = " << right*tau.first/sqrt(tau.second * tau.first) << std::endl;
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		/*
		TGraph testgraph(xPoints.size(), &xPoints[0], &yPoints[0]);
		testgraph.SetMinimum(0);	
		//testgraph.SetMaximum(0.008);	
		testgraph.GetXaxis()->SetTitle("y");
		testgraph.GetYaxis()->SetTitle("Acceptance ratio");
		//TLine line(0, corrRatio*promptAcc*corrRatio, 1, corrRatio*promptAcc*corrRatio);
		TLine line(0, promptAcc, 10, promptAcc);
		line.SetLineColor(kBlue);
		TCanvas* c = new TCanvas("", "");
		c->SetLeftMargin(0.12);
		testgraph.Draw("");
		//testgraph.SetMarkerSize(3);
		testgraph.Draw("a*same");
		line.Draw("lsame");
		c->SaveAs("testplot_m10_inverted.pdf");
		*/
	}
    return 0;
}


