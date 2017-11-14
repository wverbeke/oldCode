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
	const unsigned nPoints = 17; 
	const TString fileNames[nPoints] = {"higgsCombinem5korean.Asymptotic.mH120.root", 
										"higgsCombinem10korean.Asymptotic.mH120.root",
										"higgsCombinem20korean.Asymptotic.mH120.root",
										"higgsCombinem30korean.Asymptotic.mH120.root",
										"higgsCombinem40korean.Asymptotic.mH120.root",
										"higgsCombinem50korean.Asymptotic.mH120.root",
										"higgsCombinem60korean.Asymptotic.mH120.root",
										"higgsCombinem70korean.Asymptotic.mH120.root",
										"higgsCombinem90korean.Asymptotic.mH120.root",
										"higgsCombinem100korean.Asymptotic.mH120.root",
										"higgsCombinem150korean.Asymptotic.mH120.root",
										"higgsCombinem200korean.Asymptotic.mH120.root",
										"higgsCombinem300korean.Asymptotic.mH120.root",
										"higgsCombinem400korean.Asymptotic.mH120.root",
										"higgsCombinem500korean.Asymptotic.mH120.root",
										"higgsCombinem700korean.Asymptotic.mH120.root",	
										"higgsCombinem1000korean.Asymptotic.mH120.root",
};
	const float massPoints[nPoints] = {5, 10, 20, 30, 40, 50, 60, 70, 90, 100, 150, 200, 300, 400, 500, 700, 1000};
	const float massCouplings[nPoints] = {0.00001, 0.00001, 0.00001, 0.00001, 0.00001,  0.00001, 0.00001, 0.0001, 0.001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
	//Read the limits and uncertainties for every point
	TFile* limitFiles[nPoints];
	for(unsigned p = 0; p < nPoints; ++p){
		limitFiles[p] = TFile::Open("../limits/" + fileNames[p]);
	}
	float limits[5][nPoints]; //0: central, 1: 65% down 2: 65% up, 2: 95% down, 3: 95% up
	for(unsigned p = 0; p < nPoints; ++p){
		TTreeReader limitReader("limit",limitFiles[p]);
		TTreeReaderValue<double> limitVal(limitReader, "limit");
		TTreeReaderValue<float> quantVal(limitReader, "quantileExpected");
		//std::pair<Double_t,Double_t> limitAndQuantile[6];
		float limitAndQuantile[5][2];
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
			else ;
			++i;
		}
		/*
		for(unsigned i = 0; i < 6; ++i){
			if(limitAndQuantile[i][1] == 0.0250000) limits[3][p] = limitAndQuantile[i][0];
			else if(limitAndQuantile[i][1] == 0.1599999) limits[1][p] = limitAndQuantile[i][0];
			else if(limitAndQuantile[i][1] == 0.5) limits[0][p] = limitAndQuantile[i][0];
			else if(limitAndQuantile[i][1] == 0.8399999) limits[2][p] = limitAndQuantile[i][0];
			else if(limitAndQuantile[i][1] ==  0.9750000) limits[4][p] =  limitAndQuantile[i][0];
			else ;
		}
		*/
	}
	//Correct for couplings
	for(unsigned p = 0; p < nPoints; ++p){
		for(unsigned i = 0; i < 5; ++i){
			limits[i][p] *= massCouplings[p];
			//cout << limits[i][p] << endl;
		}
	}
	//make graph for errors and expected limit
	/*
	TGraph* limitGraph[5];
	for(unsigned i = 0; i < 5; ++i){
		limitGraph[i] = new TGraphErrors(12, massPoints, limits[i]);
	}
	TCanvas* c = new TCanvas("canv","canv",500,500);
	c->SetTopMargin(0.08);
	c->SetLeftMargin(0.16);
    //h->GetXaxis()->SetTitleOffset(1.4);
	for(unsigned i = 0; i < 5; ++i){
		limitGraph[i]->GetYaxis()->SetTitleOffset(1.2);
		limitGraph[i]->GetXaxis()->SetTitle("m_{N} (GeV)");
		limitGraph[i]->GetYaxis()->SetTitle("upper limit: #Sigma_{l}(|V_{l#nu}|^{2}");
	}
	c->SetLogy();
	//lowMgraph->GetYaxis()->SetRangeUser(0, 0.75);
	limitGraph[0]->Draw("CAE");
	*/
	/*
	TGraph *shade = new TGraph(2*nPoints);
	TGraph *shade2 = new TGraph(2*nPoints);
	for(unsigned p = 0; p < nPoints; ++p){
		shade->SetPoint(p, massPoints[p], limits[1][p]);
		shade->SetPoint(nPoints + p, massPoints[nPoints - p - 1], limits[2][nPoints - p - 1]);
		shade2->SetPoint(p, massPoints[p], limits[3][p]);
		shade2->SetPoint(nPoints + p, massPoints[nPoints - p - 1], limits[4][nPoints - p - 1]);
	}
  	shade->SetFillStyle(1001);
   	shade->SetFillColor(kGreen + 1);
	shade2->SetFillStyle(1001);
   	shade2->SetFillColor(kOrange);
	shade2->Draw("fsame");
	shade->Draw("fsame");
	for(unsigned i = 1; i < 5; ++i){
		limitGraph[0]->Draw("same");
		//limitGraph[i]->Draw("same");
	}
	TLegend* legend = new TLegend(0.6,0.3,0.9,0.6,NULL,"brNDC");
	//Avoid legend box
	legend->SetFillStyle(0);
	legend->AddEntry( (TObject*)nullptr, "95\% CL upper limits", "");
	legend->AddEntry(shade2, "\\pm 2 std. deviation", "f");	
	legend->AddEntry(shade, "\\pm 1 std. deviation", "f");	
	legend->AddEntry(limitGraph[0], "expected limit", "l");
	legend->Draw("same");
	//Add data to legend
	*/
	TGraphAsymmErrors* limitGraph[5];
	for(unsigned i = 0; i < 5; ++i){
		limitGraph[i] = new TGraphAsymmErrors(nPoints, massPoints, limits[i]);
	}
	TCanvas* c = new TCanvas("canv","canv",500,500);
	c->SetTopMargin(0.08);
	c->SetLeftMargin(0.16);
    //h->GetXaxis()->SetTitleOffset(1.4);
	for(unsigned i = 0; i < 5; ++i){
		//limitGraph[i]->GetYaxis()->SetTitleOffset(1.2);
		limitGraph[i]->GetXaxis()->SetTitle("m_{N} (GeV)");
		limitGraph[i]->GetYaxis()->SetTitle("|V_{#muN}|^{2}");
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
	
	//TLegend* legend = new TLegend(0.6,0.3,0.9,0.6,NULL,"brNDC");
	TLegend* legend = new TLegend(0.3,0.6,0.6,0.9,NULL,"brNDC");
	legend->SetFillStyle(0);
	//legend->AddEntry( (TObject*)nullptr, "95\% CL upper limits", "");
	legend->SetHeader("95\% CL upper limits");
	legend->AddEntry(twoSigmaBand, "\\pm 2 std. deviation", "f");	
	legend->AddEntry(oneSigmaBand, "\\pm 1 std. deviation", "f");	
	legend->AddEntry(limitGraph[0], "expected", "l");
	legend->Draw("same");

	drawLumi(c);
    c->SaveAs("plots/hnl_Limits_koreanAnalysis.pdf");
	//return 0;
}


/*
void grshade() {
   TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);

   c1->SetGrid();
   c1->DrawFrame(0,0,2.2,12);
   
   const Int_t n = 20;
   Double_t x[n], y[n],ymin[n], ymax[n];
   Int_t i;
   for (i=0;i<n;i++) {
     x[i] = 0.1+i*0.1;
     ymax[i] = 10*sin(x[i]+0.2);
     ymin[i] = 8*sin(x[i]+0.1);
     y[i] = 9*sin(x[i]+0.15);
   }
   TGraph *grmin = new TGraph(n,x,ymin);
   TGraph *grmax = new TGraph(n,x,ymax);
   TGraph *gr    = new TGraph(n,x,y);
   TGraph *grshade = new TGraph(2*n);
   for (i=0;i<n;i++) {
      grshade->SetPoint(i,x[i],ymax[i]);
      grshade->SetPoint(n+i,x[n-i-1],ymin[n-i-1]);
   }
   grshade->SetFillStyle(3013);
   grshade->SetFillColor(16);
   grshade->Draw("f");
   grmin->Draw("l");
   grmax->Draw("l");
   gr->SetLineWidth(4);
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->Draw("CP");
}
*/













