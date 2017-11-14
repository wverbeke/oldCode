#include "TFile.h"
#include "TH2D.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TROOT.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include <iostream>
#include <utility>      // std::pair
using std::cout;
using std::endl;
//Include other parts of the code
#include "drawLumi.h"
#include "tdrstyle.h"
void plotLimits(){
	setTDRStyle();
	gROOT->SetBatch(kTRUE);
	const unsigned nPoints = 18; 
	const TString fileNames[nPoints] = {"higgsCombinem1_significanceScan.Asymptotic.mH120.root", "higgsCombinem2_significanceScan.Asymptotic.mH120.root", "higgsCombinem5_significanceScan.Asymptotic.mH120.root", "higgsCombinem10_significanceScan.Asymptotic.mH120.root", "higgsCombinem20_significanceScan.Asymptotic.mH120.root", "higgsCombinem30_significanceScan.Asymptotic.mH120.root", "higgsCombinem40_significanceScan.Asymptotic.mH120.root", "higgsCombinem50_significanceScan.Asymptotic.mH120.root", "higgsCombinem60_significanceScan.Asymptotic.mH120.root", "higgsCombinem80_significanceScan.Asymptotic.mH120.root", "higgsCombinem100_significanceScan.Asymptotic.mH120.root", "higgsCombinem130_significanceScan.Asymptotic.mH120.root", "higgsCombinem150_significanceScan.Asymptotic.mH120.root", "higgsCombinem200_significanceScan.Asymptotic.mH120.root", "higgsCombinem400_significanceScan.Asymptotic.mH120.root", "higgsCombinem600_significanceScan.Asymptotic.mH120.root", "higgsCombinem800_significanceScan.Asymptotic.mH120.root", "higgsCombinem1000_significanceScan.Asymptotic.mH120.root"};
	const float massPoints[nPoints] = {1, 2, 5, 10, 20, 30, 40, 50, 60, 80, 100, 130, 150, 200, 400, 600, 800, 1000};
	const float massCouplings[nPoints] = {0.00001, 0.00001, 0.00001, 0.00001, 0.00001,  0.00001, 0.00001, 0.00001, 0.00001, 0.0001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
	//Read the limits and uncertainties for every point
	TFile* limitFiles[nPoints];
	for(unsigned p = 0; p < nPoints; ++p){
		limitFiles[p] = TFile::Open("../limits/significanceScan/" + fileNames[p]);
	}
	float limits[5][nPoints]; //0: central, 1: 65% down 2: 65% up, 2: 95% down, 3: 95% up
	float observedLimit[nPoints];
	for(unsigned p = 0; p < nPoints; ++p){
		TTreeReader limitReader("limit",limitFiles[p]);
		TTreeReaderValue<double> limitVal(limitReader, "limit");
		TTreeReaderValue<float> quantVal(limitReader, "quantileExpected");
		//std::pair<Double_t,Double_t> limitAndQuantile[6];
		//float limitAndQuantile[5][2];
		unsigned i = 0;
		while(limitReader.Next()){
			//limitAndQuantile[i] = (std::pair<Double_t,Double_t>) std::make_pair<(constexpr) *limitVal, (constexpr) *quantVal>;
			//limitAndQuantile[i][0] = (float) *limitVal;
			//limitAndQuantile[i][1] = *quantVal;
			if(i == 0) limits[3][p] = (float) *limitVal;
			else if(i == 1) limits[1][p] = (float) *limitVal;
			else if(i == 2) limits[0][p] = (float) *limitVal;
			else if(i == 3) limits[2][p] = (float) *limitVal;
			else if(i == 4) limits[4][p] = (float) *limitVal;
			else if(i == 5) observedLimit[p] = (float) *limitVal;
			++i;
		}
	}
	//Correct for couplings
	for(unsigned p = 0; p < nPoints; ++p){
		observedLimit[p] *= massCouplings[p]*2;
		for(unsigned i = 0; i < 5; ++i){
			limits[i][p] *= massCouplings[p]*2;
			//cout << limits[i][p] << endl;
		}
	}
	//make graph for errors and expected limit
	TGraphAsymmErrors* limitGraph[5];
	for(unsigned i = 0; i < 5; ++i){
		limitGraph[i] = new TGraphAsymmErrors(nPoints, massPoints, limits[i]);
	}
	TCanvas* c = new TCanvas("canv","canv",500,500);
	c->SetTopMargin(0.08);
	c->SetLeftMargin(0.18);
	c->SetBottomMargin(0.14);
    //h->GetXaxis()->SetTitleOffset(1.4);
	for(unsigned i = 0; i < 5; ++i){
		limitGraph[i]->GetYaxis()->SetTitleOffset(1.4);
		limitGraph[i]->GetXaxis()->SetTitleOffset(1);
		limitGraph[i]->GetXaxis()->SetTitle("m_{N} (GeV)");
		limitGraph[i]->GetYaxis()->SetTitle("|V_{e}|^{2} + |V_{#mu}|^{2}");
	}
	c->SetLogy();
	//c->SetLogx();
	//lowMgraph->GetYaxis()->SetRangeUser(0, 0.75);
	//limitGraph[0]->Draw("CAE");
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
	twoSigmaBand->Draw("CAE3");
	for(unsigned p = 0; p < nPoints; ++p){
		double yValDown;
		double yValUp;
		double yValCenter;
		double xVal;
		limitGraph[0]->GetPoint(p, xVal, yValCenter);
		std::cout << "limit for mass : " << massPoints[p]  << " is " <<  yValCenter << std::endl;
		limitGraph[1]->GetPoint(p, xVal, yValDown);	
		limitGraph[2]->GetPoint(p, xVal, yValUp);
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
	
	TGraphAsymmErrors* observed = new TGraphAsymmErrors(nPoints, massPoints, observedLimit);
	observed->SetLineStyle(1);
	observed->SetLineWidth(3);
	observed->Draw("same");
	
	TLegend* legend = new TLegend(0.6,0.3,0.9,0.6,NULL,"brNDC");
	//TLegend* legend = new TLegend(0.3,0.6,0.6,0.9,NULL,"brNDC");
	legend->SetFillStyle(0);
	//legend->AddEntry( (TObject*)nullptr, "95\% CL upper limits", "");
	legend->SetHeader("95\% CL upper limits");
	legend->AddEntry(twoSigmaBand, "\\pm 2 std. deviation", "f");	
	legend->AddEntry(oneSigmaBand, "\\pm 1 std. deviation", "f");	
	legend->AddEntry(limitGraph[0], "expected", "l");
	legend->AddEntry(observed, "pseudo-observed", "l");
	legend->Draw("same");

	drawLumi(c);
    c->SaveAs("plots/hnl_SignificanceScanm30and200.pdf");
	c->SaveAs("plots/hnl_SignificanceScanm30and200.png");
	c->Close();

	//Draw logarithmic plot
	TCanvas* clog = new TCanvas("canv","canv",500,500);
	clog->SetTopMargin(0.08);
	clog->SetLeftMargin(0.18);
	clog->SetBottomMargin(0.14);
	clog->SetLogy();
	clog->SetLogx();
	twoSigmaBand->Draw("CAE3");
	oneSigmaBand->Draw("CE3same");
	limitGraph[0]->Draw("Lsame");
	observed->Draw("Lsame");
	TLegend* legend2 = new TLegend(0.3,0.6,0.6,0.9,NULL,"brNDC");
	legend2->SetFillStyle(0);
	//legend->AddEntry( (TObject*)nullptr, "95\% CL upper limits", "");
	legend2->SetHeader("95\% CL upper limits");
	legend2->AddEntry(twoSigmaBand, "\\pm 2 std. deviation", "f");	
	legend2->AddEntry(oneSigmaBand, "\\pm 1 std. deviation", "f");	
	legend2->AddEntry(limitGraph[0], "expected", "l");
	legend2->AddEntry(observed, "pseudo-observed", "l");
	legend2->Draw("same");

	drawLumi(clog);
    clog->SaveAs("plots/hnl_SignificanceScanm30and200_logx.pdf");
	clog->SaveAs("plots/hnl_SignificanceScanm30and200_logx.png");
	//return 0;
}















