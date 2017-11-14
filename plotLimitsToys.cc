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
	const unsigned nPoints = 16; 


	//REMOVED 30 GEV and 130 GEV POINT, PUT IT BACK
	const TString fileNames[nPoints] = {"higgsCombinem1.Asymptotic.mH120.root", "higgsCombinem2.Asymptotic.mH120.root", "higgsCombinem5.Asymptotic.mH120.root", "higgsCombinem10.Asymptotic.mH120.root", "higgsCombinem20.Asymptotic.mH120.root", "higgsCombinem40.Asymptotic.mH120.root", "higgsCombinem50.Asymptotic.mH120.root", "higgsCombinem60.Asymptotic.mH120.root", "higgsCombinem80.Asymptotic.mH120.root", "higgsCombinem100.Asymptotic.mH120.root", "higgsCombinem150.Asymptotic.mH120.root", "higgsCombinem200.Asymptotic.mH120.root", "higgsCombinem400.Asymptotic.mH120.root", "higgsCombinem600.Asymptotic.mH120.root", "higgsCombinem800.Asymptotic.mH120.root", "higgsCombinem1000.Asymptotic.mH120.root"};
	/*
	const TString fileNamesToys[nPoints] = {"higgsCombinem1.HybridNew.mH120.root", "higgsCombinem2.HybridNew.mH120.root", "higgsCombinem5.HybridNew.mH120.root", "higgsCombinem10.HybridNew.mH120.root", "higgsCombinem20.HybridNew.mH120.root", "higgsCombinem30.HybridNew.mH120.root", "higgsCombinem40.HybridNew.mH120.root", "higgsCombinem50.HybridNew.mH120.root", "higgsCombinem60.HybridNew.mH120.root", "higgsCombinem80.HybridNew.mH120.root", "higgsCombinem100.HybridNew.mH120.root", "higgsCombinem130.HybridNew.mH120.root", "higgsCombinem150.HybridNew.mH120.root", "higgsCombinem200.HybridNew.mH120.root", "higgsCombinem400.HybridNew.mH120.root", "higgsCombinem600.HybridNew.mH120.root", "higgsCombinem800.HybridNew.mH120.root", "higgsCombinem1000.HybridNew.mH120.root"};
	*/
	const TString fileNamesToysCentral[nPoints] = {"higgsCombinem1_0.5.HybridNew.mH120.quant0.500.root", "higgsCombinem2_0.5.HybridNew.mH120.quant0.500.root", "higgsCombinem5_0.5.HybridNew.mH120.quant0.500.root", "higgsCombinem10_0.5.HybridNew.mH120.quant0.500.root", "higgsCombinem20_0.5.HybridNew.mH120.quant0.500.root", "higgsCombinem40_0.5.HybridNew.mH120.quant0.500.root", "higgsCombinem50_0.5.HybridNew.mH120.quant0.500.root", "higgsCombinem60_0.5.HybridNew.mH120.quant0.500.root", "higgsCombinem80_0.5.HybridNew.mH120.quant0.500.root", "higgsCombinem100_0.5.HybridNew.mH120.quant0.500.root",  "higgsCombinem150_0.5.HybridNew.mH120.quant0.500.root", "higgsCombinem200_0.5.HybridNew.mH120.quant0.500.root", "higgsCombinem400_0.5.HybridNew.mH120.quant0.500.root", "higgsCombinem600_0.5.HybridNew.mH120.quant0.500.root", "higgsCombinem800_0.5.HybridNew.mH120.quant0.500.root", "higgsCombinem1000_0.5.HybridNew.mH120.quant0.500.root"};

	const TString fileNamesToys1SigmaDown[nPoints] = {"higgsCombinem1_0.16.HybridNew.mH120.quant0.160.root", "higgsCombinem2_0.16.HybridNew.mH120.quant0.160.root", "higgsCombinem5_0.16.HybridNew.mH120.quant0.160.root", "higgsCombinem10_0.16.HybridNew.mH120.quant0.160.root", "higgsCombinem20_0.16.HybridNew.mH120.quant0.160.root", "higgsCombinem40_0.16.HybridNew.mH120.quant0.160.root", "higgsCombinem50_0.16.HybridNew.mH120.quant0.160.root", "higgsCombinem60_0.16.HybridNew.mH120.quant0.160.root", "higgsCombinem80_0.16.HybridNew.mH120.quant0.160.root", "higgsCombinem100_0.16.HybridNew.mH120.quant0.160.root",  "higgsCombinem150_0.16.HybridNew.mH120.quant0.160.root", "higgsCombinem200_0.16.HybridNew.mH120.quant0.160.root", "higgsCombinem400_0.16.HybridNew.mH120.quant0.160.root", "higgsCombinem600_0.16.HybridNew.mH120.quant0.160.root", "higgsCombinem800_0.16.HybridNew.mH120.quant0.160.root", "higgsCombinem1000_0.16.HybridNew.mH120.quant0.160.root"};

	const TString fileNamesToys1SigmaUp[nPoints] =  {"higgsCombinem1_0.84.HybridNew.mH120.quant0.840.root", "higgsCombinem2_0.84.HybridNew.mH120.quant0.840.root", "higgsCombinem5_0.84.HybridNew.mH120.quant0.840.root", "higgsCombinem10_0.84.HybridNew.mH120.quant0.840.root", "higgsCombinem20_0.84.HybridNew.mH120.quant0.840.root", "higgsCombinem40_0.84.HybridNew.mH120.quant0.840.root", "higgsCombinem50_0.84.HybridNew.mH120.quant0.840.root", "higgsCombinem60_0.84.HybridNew.mH120.quant0.840.root", "higgsCombinem80_0.84.HybridNew.mH120.quant0.840.root", "higgsCombinem100_0.84.HybridNew.mH120.quant0.840.root",  "higgsCombinem150_0.84.HybridNew.mH120.quant0.840.root", "higgsCombinem200_0.84.HybridNew.mH120.quant0.840.root", "higgsCombinem400_0.84.HybridNew.mH120.quant0.840.root", "higgsCombinem600_0.84.HybridNew.mH120.quant0.840.root", "higgsCombinem800_0.84.HybridNew.mH120.quant0.840.root", "higgsCombinem1000_0.84.HybridNew.mH120.quant0.840.root"};

	const TString fileNamesToys2SigmaDown[nPoints] = {"higgsCombinem1_0.025.HybridNew.mH120.quant0.025.root", "higgsCombinem2_0.025.HybridNew.mH120.quant0.025.root", "higgsCombinem5_0.025.HybridNew.mH120.quant0.025.root", "higgsCombinem10_0.025.HybridNew.mH120.quant0.025.root", "higgsCombinem20_0.025.HybridNew.mH120.quant0.025.root", "higgsCombinem40_0.025.HybridNew.mH120.quant0.025.root", "higgsCombinem50_0.025.HybridNew.mH120.quant0.025.root", "higgsCombinem60_0.025.HybridNew.mH120.quant0.025.root", "higgsCombinem80_0.025.HybridNew.mH120.quant0.025.root", "higgsCombinem100_0.025.HybridNew.mH120.quant0.025.root",  "higgsCombinem150_0.025.HybridNew.mH120.quant0.025.root", "higgsCombinem200_0.025.HybridNew.mH120.quant0.025.root", "higgsCombinem400_0.025.HybridNew.mH120.quant0.025.root", "higgsCombinem600_0.025.HybridNew.mH120.quant0.025.root", "higgsCombinem800_0.025.HybridNew.mH120.quant0.025.root", "higgsCombinem1000_0.025.HybridNew.mH120.quant0.025.root"};

	const TString fileNamesToys2SigmaUp[nPoints] = {"higgsCombinem1_0.975.HybridNew.mH120.quant0.975.root", "higgsCombinem2_0.975.HybridNew.mH120.quant0.975.root", "higgsCombinem5_0.975.HybridNew.mH120.quant0.975.root", "higgsCombinem10_0.975.HybridNew.mH120.quant0.975.root", "higgsCombinem20_0.975.HybridNew.mH120.quant0.975.root", "higgsCombinem40_0.975.HybridNew.mH120.quant0.975.root", "higgsCombinem50_0.975.HybridNew.mH120.quant0.975.root", "higgsCombinem60_0.975.HybridNew.mH120.quant0.975.root", "higgsCombinem80_0.975.HybridNew.mH120.quant0.975.root", "higgsCombinem100_0.975.HybridNew.mH120.quant0.975.root",  "higgsCombinem150_0.975.HybridNew.mH120.quant0.975.root", "higgsCombinem200_0.975.HybridNew.mH120.quant0.975.root", "higgsCombinem400_0.975.HybridNew.mH120.quant0.975.root", "higgsCombinem600_0.975.HybridNew.mH120.quant0.975.root", "higgsCombinem800_0.975.HybridNew.mH120.quant0.975.root", "higgsCombinem1000_0.975.HybridNew.mH120.quant0.975.root"};

		
	const float massPoints[nPoints] = {1, 2, 5, 10, 20, 40, 50, 60, 80, 100, 150, 200, 400, 600, 800, 1000};
	const float massCouplings[nPoints] = {0.00001, 0.00001, 0.00001, 0.00001,  0.00001, 0.00001, 0.00001, 0.00001, 0.0001, 0.01,  0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
	//Read the limits and uncertainties for every point
	TFile* limitFiles[nPoints];
	for(unsigned p = 0; p < nPoints; ++p){
		limitFiles[p] = TFile::Open("../limits/FullToys/" + fileNames[p]);
	}
	TFile* limitFilesToys[5][nPoints];
	for(unsigned p =0; p < nPoints; ++p){
		limitFilesToys[0][p] = TFile::Open("../limits/FullToys/" + fileNamesToysCentral[p]);
		limitFilesToys[1][p] = TFile::Open("../limits/FullToys/" + fileNamesToys1SigmaDown[p]);
		limitFilesToys[2][p] = TFile::Open("../limits/FullToys/" + fileNamesToys1SigmaUp[p]);
		limitFilesToys[3][p] = TFile::Open("../limits/FullToys/" + fileNamesToys2SigmaDown[p]);
		limitFilesToys[4][p] = TFile::Open("../limits/FullToys/" + fileNamesToys2SigmaUp[p]);
	}
	float limits[5][nPoints]; //0: central, 1: 65% down 2: 65% up, 2: 95% down, 3: 95% up
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
			else ;
			++i;
		}
	}
	float toyLimits[5][nPoints];
	for(unsigned p = 0; p < nPoints; ++p){
		for(unsigned i = 0; i < 5; ++i){
			TTreeReader limitReader("limit",limitFilesToys[i][p]);
			TTreeReaderValue<double> limitVal(limitReader, "limit");
			TTreeReaderValue<float> quantVal(limitReader, "quantileExpected");	
			while(limitReader.Next()){
				toyLimits[i][p] = (float) *limitVal;
			}
		}
	}
	//Correct for couplings
	for(unsigned p = 0; p < nPoints; ++p){
		for(unsigned i = 0; i < 5; ++i){
			limits[i][p] *= massCouplings[p]*2;
			toyLimits[i][p] *= massCouplings[p]*2;
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

	//TGraph containing limit with toys
	/*
	TGraphAsymmErrors* toyLimitGraph = new TGraphAsymmErrors(nPoints, massPoints, toyLimits);
	toyLimitGraph->SetLineStyle(1);
	toyLimitGraph->SetLineWidth(3);
	toyLimitGraph->Draw("same");
	*/
	TGraphAsymmErrors* limitGraphToys[5];
	for(unsigned i = 0; i < 5; ++i){
		limitGraphToys[i] = new TGraphAsymmErrors(nPoints, massPoints, toyLimits[i]);
	}
	for(unsigned p = 0; p < nPoints; ++p){
		double yValDown;
		double yValUp;
		double yValCenter;
		double xVal;
		limitGraphToys[0]->GetPoint(p, xVal, yValCenter);
		limitGraphToys[3]->GetPoint(p, xVal, yValDown);	
		limitGraphToys[4]->GetPoint(p, xVal, yValUp);
		limitGraphToys[0]->SetPointError(p,0 , 0, yValCenter - yValDown, yValUp - yValDown);
	}
	TGraphAsymmErrors* twoSigmaBandToys = (TGraphAsymmErrors*) limitGraphToys[0]->Clone();
	//twoSigmaBandToys->SetFillStyle(3244);
  	//twoSigmaBandToys->SetFillColor(kBlue);
	twoSigmaBandToys->SetLineColor(kBlue);
	twoSigmaBandToys->SetLineWidth(2);
	twoSigmaBandToys->Draw("Esame");
	for(unsigned p = 0; p < nPoints; ++p){
		double yValDown;
		double yValUp;
		double yValCenter;
		double xVal;
		limitGraph[0]->GetPoint(p, xVal, yValCenter);
		std::cout << "limit for mass : " << massPoints[p]  << " is " <<  yValCenter << std::endl;
		limitGraphToys[1]->GetPoint(p, xVal, yValDown);	
		limitGraphToys[2]->GetPoint(p, xVal, yValUp);
		limitGraphToys[0]->SetPointError(p,0 , 0, yValCenter - yValDown, yValUp - yValDown);
	}
	TGraphAsymmErrors* oneSigmaBandToys = (TGraphAsymmErrors*) limitGraphToys[0]->Clone();
  	//oneSigmaBandToys->SetFillStyle(3244); //3005  3244
    //oneSigmaBandToys->SetFillColor(kRed + 1 );
	oneSigmaBandToys->Draw("Esame");
	oneSigmaBandToys->SetLineWidth(2);	
	limitGraphToys[0]->SetLineStyle(1);
	limitGraphToys[0]->SetLineWidth(3);	
	for(unsigned p = 0; p < nPoints; ++p){
		limitGraphToys[0]->SetPointError(p,0 , 0, 0, 0);
	}
	limitGraphToys[0]->Draw("Lsame");

	/*
	//Limit ratio print:
	std::cout << "printing limit ratios for Asymptotic vs full frequentist calculation:" <<std::endl;
	for(unsigned p = 0; p < nPoints; ++p){
		std::cout << "m_{N} = " << massPoints[p] << " has ratio : full frequensist/Asymptotic = " << limits[0][p]/toyLimits[p] << std::endl;
	}
	*/
	TLegend* legend = new TLegend(0.6,0.3,0.9,0.6,NULL,"brNDC");
	//TLegend* legend = new TLegend(0.3,0.6,0.6,0.9,NULL,"brNDC");
	legend->SetFillStyle(0);
	//legend->AddEntry( (TObject*)nullptr, "95\% CL upper limits", "");
	legend->SetHeader("95\% CL upper limits");
	legend->AddEntry(twoSigmaBand, "\\pm 2 std. deviation asymptotic", "f");	
	legend->AddEntry(oneSigmaBand, "\\pm 1 std. deviation asymptotic", "f");	
	legend->AddEntry(limitGraph[0], "asymptotic expected", "l");
	legend->AddEntry(twoSigmaBandToys, "\\pm 2 std. deviation frequentist", "e");	
	legend->AddEntry(oneSigmaBandToys, "\\pm 1 std. deviation frequentist", "e");	
	legend->AddEntry(limitGraphToys[0], "fully frequentist expected", "l");
	legend->Draw("same");

	drawLumi(c);
    c->SaveAs("plots/hnl_Limits_VSToys.pdf");
	c->SaveAs("plots/hnl_Limits_VSToys.png");
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
	twoSigmaBandToys->Draw("Esame");
	oneSigmaBandToys->Draw("Esame");
	limitGraphToys[0]->Draw("Lsame");
	TLegend* legend2 = new TLegend(0.3,0.6,0.6,0.9,NULL,"brNDC");
	legend2->SetFillStyle(0);
	//legend->AddEntry( (TObject*)nullptr, "95\% CL upper limits", "");
	legend2->SetHeader("95\% CL upper limits");
	legend2->AddEntry(twoSigmaBand, "\\pm 2 std. deviation asymptotic", "f");	
	legend2->AddEntry(oneSigmaBand, "\\pm 1 std. deviation asymptotic", "f");	
	legend2->AddEntry(limitGraph[0], "asymptotic expected", "l");
	legend2->AddEntry(twoSigmaBandToys, "\\pm 2 std. deviation frequentist", "e");	
	legend2->AddEntry(oneSigmaBandToys, "\\pm 1 std. deviation frequentist", "e");	
	legend2->AddEntry(limitGraphToys[0], "fully frequentist expected", "l");
	legend2->Draw("same");

	drawLumi(clog);
    clog->SaveAs("plots/hnl_Limits_VSToys_logx.pdf");
	clog->SaveAs("plots/hnl_Limits_VSToys_logx.png");
	//return 0;
}















