#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <tuple>
#include <set>

#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH2.h"
#include "TROOT.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "drawLumi.h"
#include "tdrstyle.h"


void plotAllLimits(){
	setTDRStyle();
	//plot muon limits
	std::vector<std::string> fileNames = {"atlas_8mu", "belle_mu", "delphi_short", "delphi_long", "fmmf", "kmumupi", "kmunu", "lhcb", "nutevmu", "bebc", "charm2b", "cms_8mu", "e949", "l3_mu", "na3_mu", "ps191_mu"};
	std::vector<std::string> expNames = {"ATLAS",  "Belle", "Delphi (prompt)", "Delphi (displaced)", "FMMF", "K #rightarrow #mu#mu#pi", "K #rightarrow #mu#mu", "LHCb", "NuTeV", "BEBC", "Charm", "CMS (8 TeV)", "E949", "L3", "NA3", "PS191"};
	for(std::vector<std::string>::iterator it = fileNames.begin(); it != fileNames.end(); ++it){
		it->append(".dat");
		std::string temp = "../../Downloads/hnlLimits/known/muons/";
		it->insert(it->begin(), temp.cbegin(), temp.cend());
	}
	std::cout << "crash 1" << std::endl;
	std::vector< std::vector<double> > xlimits;
	std::vector< std::vector<double> > ylimits;
	for(std::vector<std::string>::const_iterator it = fileNames.cbegin(); it != fileNames.cend(); ++it){
		std::ifstream file(*it);
		std::vector<double> tempx, tempy;
		while(!file.eof()){
			double xval, yval;
			file >> xval >> yval;
			tempx.push_back(xval);
			tempy.push_back(yval);
		}
		xlimits.push_back(tempx);
		ylimits.push_back(tempy);
	}
	std::map<std::string, TGraphAsymmErrors*> limitGraphs;
	for(unsigned i = 0; i < expNames.size(); ++i){
		if(expNames[i] == "LHCb" || expNames[i] == "E949"){
			for(auto it = xlimits[i].begin(); it != xlimits[i].end(); ++it){
				(*it) *= 0.001;
			}
		} else if(expNames[i] ==  "Delphi (prompt)" || expNames[i] ==  "Delphi (displaced)"){
			for(unsigned j = 0; j < xlimits[i].size(); ++j){
				xlimits[i][j] = pow(10.0, xlimits[i][j]);
				ylimits[i][j] = pow(10.0, ylimits[i][j]);
			}
		}
		limitGraphs[expNames[i]] = new TGraphAsymmErrors(xlimits[i].size(), &xlimits[i][0], &ylimits[i][0]);		
	}
	//make canvas
	TCanvas* c = new TCanvas("canv","canv",1000,500);
	c->SetTopMargin(0.08);
	c->SetLeftMargin(0.09);
	c->SetBottomMargin(0.14);
	c->SetLogy();
	c->SetLogx();

	TPad* p1, *p2;
    p1 = new TPad("p1","",0,0,0.85,1);
    p1->Draw();
    p1->cd();
    //p1->SetBottomMargin(0);
	p1->SetTopMargin(0.08);
	p1->SetLeftMargin(0.12);
	p1->SetBottomMargin(0.16);
	p1->SetRightMargin(0.03);
	p1->SetLogy();
	p1->SetLogx();
	
	
	limitGraphs["ATLAS"]->SetTitle("");
	limitGraphs["ATLAS"]->GetYaxis()->SetTitleOffset(0.95);
	limitGraphs["ATLAS"]->GetXaxis()->SetTitleOffset(1.2);
	limitGraphs["ATLAS"]->GetXaxis()->SetTitle("m_{N} (GeV)");
	limitGraphs["ATLAS"]->GetYaxis()->SetTitle("|V_{#muN}|^{2}");
	limitGraphs["ATLAS"]->GetXaxis()->SetLimits(0.1, 1000);
	limitGraphs["ATLAS"]->GetYaxis()->SetRangeUser(10e-12, 1);
	limitGraphs["ATLAS"]->Draw("lA");


	const unsigned nPoints = 17; 
	const TString fileNamesHNL[nPoints] = {"higgsCombineMuonCouplingm5.Asymptotic.mH120.root", "higgsCombineMuonCouplingm6.Asymptotic.mH120.root",  "higgsCombineMuonCouplingm10.Asymptotic.mH120.root", "higgsCombineMuonCouplingm20.Asymptotic.mH120.root", "higgsCombineMuonCouplingm30.Asymptotic.mH120.root", "higgsCombineMuonCouplingm40.Asymptotic.mH120.root", "higgsCombineMuonCouplingm50.Asymptotic.mH120.root", "higgsCombineMuonCouplingm60.Asymptotic.mH120.root", "higgsCombineMuonCouplingm80.Asymptotic.mH120.root", "higgsCombineMuonCouplingm100.Asymptotic.mH120.root", "higgsCombineMuonCouplingm130.Asymptotic.mH120.root", "higgsCombineMuonCouplingm150.Asymptotic.mH120.root", "higgsCombineMuonCouplingm200.Asymptotic.mH120.root", "higgsCombineMuonCouplingm400.Asymptotic.mH120.root", "higgsCombineMuonCouplingm600.Asymptotic.mH120.root", "higgsCombineMuonCouplingm800.Asymptotic.mH120.root", "higgsCombineMuonCouplingm1000.Asymptotic.mH120.root"};
	const float massPoints[nPoints] = {5, 6, 10, 20, 30, 40, 50, 60, 80, 100, 130, 150, 200, 400, 600, 800, 1000};
	const float massCouplings[nPoints] = {0.00001, 0.00001, 0.00001, 0.00001,  0.00001, 0.00001, 0.00001, 0.00001, 0.0001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
	
	//Read the limits and uncertainties for every point
	TFile* limitFiles[nPoints];
	for(unsigned p = 0; p < nPoints; ++p){
		limitFiles[p] = TFile::Open("../limits/august2017/" + fileNamesHNL[p]);
	}
	float limits[6][nPoints]; //0: central, 1: 65% down 2: 65% up, 2: 95% down, 3: 95% up
	for(unsigned p = 0; p < nPoints; ++p){
		TTreeReader limitReader("limit",limitFiles[p]);
		TTreeReaderValue<double> limitVal(limitReader, "limit");
		TTreeReaderValue<float> quantVal(limitReader, "quantileExpected");
		float limitAndQuantile[5][2];
		unsigned i = 0;
		while(limitReader.Next()){
			if(i == 0) limits[3][p] = (float) *limitVal;
			else if(i == 1) limits[1][p] = (float) *limitVal;
			else if(i == 2) limits[0][p] = (float) *limitVal;
			else if(i == 3) limits[2][p] = (float) *limitVal;
			else if(i == 4) limits[4][p] = (float) *limitVal;
			else if(i == 5) limits[5][p] = (float) *limitVal;
			++i;
		}
	}
	//Correct for couplings
	for(unsigned p = 0; p < nPoints; ++p){
		for(unsigned i = 0; i < 6; ++i){
			limits[i][p] *= massCouplings[p];
			if(p == 0) limits[i][p] *= (1./0.44);
			if(p == 1) limits[i][p] *= (1./0.520445);
		}
	}
	//make graph for errors and expected limit
	TGraphAsymmErrors* limitGraph[6];
	for(unsigned i = 0; i < 6; ++i){
		limitGraph[i] = new TGraphAsymmErrors(nPoints, massPoints, limits[i]);
	}

	for(unsigned p = 0; p < nPoints; ++p){
		double yValDown;
		double yValUp;
		double yValCenter;
		double xVal;
		limitGraph[0]->GetPoint(p, xVal, yValCenter);
		limitGraph[3]->GetPoint(p, xVal, yValDown);	
		limitGraph[4]->GetPoint(p, xVal, yValUp);
		limitGraph[0]->SetPointError(p,0 , 0, yValCenter - yValDown, yValUp - yValDown);
	}
	TGraphAsymmErrors* twoSigmaBand = (TGraphAsymmErrors*) limitGraph[0]->Clone();
	twoSigmaBand->SetFillStyle(1001);
  	twoSigmaBand->SetFillColor(kOrange);
	twoSigmaBand->GetYaxis()->SetTitleOffset(0.95);
	twoSigmaBand->GetXaxis()->SetTitleOffset(1.2);
	twoSigmaBand->GetXaxis()->SetTitle("m_{N} (GeV)");
	twoSigmaBand->GetYaxis()->SetTitle("|V_{#muN}|^{2}");
	twoSigmaBand->GetXaxis()->SetLimits(0.1, 1000);
	twoSigmaBand->GetYaxis()->SetRangeUser(10e-12, 1);
	twoSigmaBand->Draw("lA");
	twoSigmaBand->Draw("CAE3same");
	for(unsigned p = 0; p < nPoints; ++p){
		double yValDown;
		double yValUp;
		double yValCenter;
		double xVal;
		limitGraph[0]->GetPoint(p, xVal, yValCenter);
		limitGraph[1]->GetPoint(p, xVal, yValDown);	
		limitGraph[2]->GetPoint(p, xVal, yValUp);
		std::cout << "limit for mass : " << massPoints[p]  << " is " <<  yValCenter << std::endl;
		limitGraph[0]->SetPointError(p,0 , 0, yValCenter - yValDown, yValUp - yValDown);
	}
	TGraphAsymmErrors* oneSigmaBand = (TGraphAsymmErrors*) limitGraph[0]->Clone();
	oneSigmaBand->SetFillStyle(1001);
  	oneSigmaBand->SetFillColor(kGreen + 1);
	oneSigmaBand->Draw("CE3same");
	//limitGraph[0]->Draw("E4same");
	for(unsigned p = 0; p < nPoints; ++p){
		limitGraph[0]->SetPointError(p,0 , 0, 0,0);
	}
	limitGraph[0]->SetLineStyle(2);
	limitGraph[0]->Draw("Lsame");
	
	limitGraph[5]->SetLineStyle(1);
	limitGraph[5]->SetLineWidth(2);
	limitGraph[5]->Draw("Lsame");




	TLegend* legend = new TLegend(0,0.14*0.15/(1- 0.15),1,1 - 0.08*0.15/(1- 0.15),NULL,"brNDC"); //xPad
	legend->SetFillStyle(0);
	legend->SetTextFont(42); //62
	legend->SetTextSize(0.085);
	//legend->SetTextSize(0.06);
	
	std::vector<Color_t> colors = {kRed, kBlue, kGreen, kYellow, kCyan, kMagenta, kGreen + 3, kAzure + 1, kGreen - 7, kMagenta -7, kRed - 7, kBlue -3, kOrange + 6, kCyan + 1, kMagenta +3, kBlack};
	unsigned i = 0; 
	for(auto it = expNames.begin(); it != expNames.end(); ++it){
		limitGraphs[*it]->SetLineColor(colors[i]);
		limitGraphs[*it]->SetLineWidth(2);
		limitGraphs[*it]->Draw("lsame");
    	legend->AddEntry(limitGraphs[*it], (TString) *it, "l");	
		++i;
	}


	drawLumi(p1);
	c->cd();
	p2 = new TPad("p2","",0.85, 0, 1, 1);
	p2->Draw();
	p2->cd();
	p2->SetTopMargin(0.08);
	p2->SetRightMargin(0.01);
	p2->SetBottomMargin(0.14);
	legend->Draw("");
	c->cd();
	//drawLumi(c);
	c->SaveAs("plots/testplot.pdf");
	c->Close();
	/*
	std::map<std::string, std::string> fileMap;
	for(std::vector<std::string>::const_iterator fit = fileNames.cbegin(), eit = expNames.cbegin(); fit != fileNames.cend(), eit != expNames.cend(); ++fit, ++eit){
		fileMap[*fit] = *eit;
	}
	std::set< std::tuple<std::string, std::vector<double>, std::vector<double> > > limits;
	
	//fill limits tuple with values from txt files
	for(std::vector<std::string>::iterator it = fileNames.begin(); it != fileNames.end(); ++it){
		std::ifstream file(*it);
		std::vector<double> x,y;
		while(!file.eof()){
			double xval, yval;
			file >> xval >> yval;
			x.push_back(xval);
			y.push_back(yval); 		
		}
		std::tuple<std::string, std::vector<double>, std::vector<double> > tempTuple(fileMap[*it], x, y);
		limits.insert(tempTuple);
		TGraphAsymmErrors test((unsigned) x.size(), &x[0], &y[0]);
	}
	//Create TGraphAsymmErrors
	std::map< std::name 
	for(std::vector<std::string>::iterator it = fileNames.begin(); it != fileNames.end(); ++it){
		
	}		
	*/
}
/*
int main(){
	plotAllLimits();
	return 0;

}
*/
