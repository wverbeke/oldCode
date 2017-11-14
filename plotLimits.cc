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
	const TString fileNames[nPoints] = {"higgsCombinem1.Asymptotic.mH120.root", "higgsCombinem2.Asymptotic.mH120.root", "higgsCombinem5.Asymptotic.mH120.root", "higgsCombinem10.Asymptotic.mH120.root", "higgsCombinem20.Asymptotic.mH120.root", "higgsCombinem30.Asymptotic.mH120.root", "higgsCombinem40.Asymptotic.mH120.root", "higgsCombinem50.Asymptotic.mH120.root", "higgsCombinem60.Asymptotic.mH120.root", "higgsCombinem80.Asymptotic.mH120.root", "higgsCombinem100.Asymptotic.mH120.root", "higgsCombinem130.Asymptotic.mH120.root", "higgsCombinem150.Asymptotic.mH120.root", "higgsCombinem200.Asymptotic.mH120.root", "higgsCombinem400.Asymptotic.mH120.root", "higgsCombinem600.Asymptotic.mH120.root", "higgsCombinem800.Asymptotic.mH120.root", "higgsCombinem1000.Asymptotic.mH120.root"};
	const float massPoints[nPoints] = {1, 2, 5, 10, 20, 30, 40, 50, 60, 80, 100, 130, 150, 200, 400, 600, 800, 1000};
	const float massCouplings[nPoints] = {0.00001, 0.00001, 0.00001, 0.00001, 0.00001,  0.00001, 0.00001, 0.00001, 0.00001, 0.0001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
	//Read the limits and uncertainties for every point
	/*
	const unsigned nPoints = 16; 
	const TString fileNames[nPoints] = {"higgsCombinem5.Asymptotic.mH120.root", "higgsCombinem10.Asymptotic.mH120.root", "higgsCombinem20.Asymptotic.mH120.root", "higgsCombinem30.Asymptotic.mH120.root", "higgsCombinem40.Asymptotic.mH120.root", "higgsCombinem50.Asymptotic.mH120.root", "higgsCombinem60.Asymptotic.mH120.root", "higgsCombinem80.Asymptotic.mH120.root", "higgsCombinem100.Asymptotic.mH120.root", "higgsCombinem130.Asymptotic.mH120.root", "higgsCombinem150.Asymptotic.mH120.root", "higgsCombinem200.Asymptotic.mH120.root", "higgsCombinem400.Asymptotic.mH120.root", "higgsCombinem600.Asymptotic.mH120.root", "higgsCombinem800.Asymptotic.mH120.root", "higgsCombinem1000.Asymptotic.mH120.root"};
	const float massPoints[nPoints] = {5, 10, 20, 30, 40, 50, 60, 80, 100, 130, 150, 200, 400, 600, 800, 1000};
	const float massCouplings[nPoints] = {0.00001, 0.00001, 0.00001,  0.00001, 0.00001, 0.00001, 0.00001, 0.0001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
	*/
	//Delphi limits
	double massPointsDelphi[38] = {0.260273972602739, 0.292237442922374, 0.342465753424657, 0.390410958904109, 0.447488584474885, 0.481735159817351, 0.534246575342466, 0.591324200913242, 0.625570776255707, 0.652968036529680, 0.694063926940639, 0.732876712328766, 0.776255707762557, 0.815068493150685, 0.860730593607306, 0.915525114155251, 0.972602739726027, 1.038812785388127, 1.098173515981735, 1.166666666666666, 1.230593607305935, 1.280821917808219, 1.331050228310502, 1.381278538812785, 1.456621004566209, 1.513698630136986, 1.586757990867580, 1.630136986301369, 1.666666666666666, 1.712328767123287, 1.751141552511415, 1.785388127853881, 1.815068493150684, 1.837899543378995, 1.856164383561643, 1.872146118721461, 1.888127853881278, 1.901826484018265};
	double limitsDelphi[38] = {-2.00000000000000, -2.19424460431654, -2.49640287769784, -2.78417266187050, -3.11510791366906, -3.33812949640287, -3.63309352517985, -3.95683453237410, -4.14388489208633, -4.30935251798561, -4.47482014388489, -4.59712230215827, -4.68345323741007, -4.72661870503597, -4.75539568345323, -4.74820143884892, -4.73381294964028, -4.72661870503597, -4.72661870503597,  -4.70503597122302, -4.68345323741007, -4.67625899280575, -4.66187050359712, -4.63309352517985, -4.60431654676259, -4.60431654676259, -4.60431654676259, -4.59712230215827, -4.56115107913669, -4.48920863309352, -4.38129496402877, -4.20143884892086, -3.98561151079136, -3.73381294964028, -3.50359712230215, -3.20863309352517, -2.92805755395683, -2.64748201438848};

	for(unsigned d = 0; d < 38; ++d){
		massPointsDelphi[d] = pow(10.0, massPointsDelphi[d]);
		limitsDelphi[d] = pow(10.0, limitsDelphi[d]);
	}
	//Make Delphi limit graph
	TGraphAsymmErrors* DelphiGraph;
	DelphiGraph = new TGraphAsymmErrors(38, massPointsDelphi, limitsDelphi);


	TFile* limitFiles[nPoints];
	for(unsigned p = 0; p < nPoints; ++p){
		limitFiles[p] = TFile::Open("../limits/june2017/doubleCoupling/" + fileNames[p]);
		//limitFiles[p] = TFile::Open("../limits/new/" + fileNames[p]);
	}
	float limits[6][nPoints]; //0: central, 1: 65% down 2: 65% up, 2: 95% down, 3: 95% up
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
			else if(i == 5) limits[5][p] = (float) *limitVal;
			++i;
		}
	}
	//Correct for couplings
	for(unsigned p = 0; p < nPoints; ++p){
		for(unsigned i = 0; i < 6; ++i){
			limits[i][p] *= massCouplings[p];
			//cout << lim its[i][p] << endl;
		}
	}
	//make graph for errors and expected limit
	TGraphAsymmErrors* limitGraph[6];
	for(unsigned i = 0; i < 6; ++i){
		limitGraph[i] = new TGraphAsymmErrors(nPoints, massPoints, limits[i]);
	}
	//h->GetXaxis()->SetTitleOffset(1.4);
	for(unsigned i = 0; i < 6; ++i){
		limitGraph[i]->GetYaxis()->SetTitleOffset(1.4);
		limitGraph[i]->GetXaxis()->SetTitleOffset(1);
		limitGraph[i]->GetXaxis()->SetTitle("m_{N} (GeV)");
		//limitGraph[i]->GetYaxis()->SetTitle("|V_{e} + V_{#mu}|^{2}");
		limitGraph[i]->GetYaxis()->SetTitle("|V_{eN}|^{2} + |V_{#muN}|^{2}");
	}
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
	//limitGraph[0]->Draw("E4same");
	for(unsigned p = 0; p < nPoints; ++p){
		limitGraph[0]->SetPointError(p,0 , 0, 0,0);
	}
	limitGraph[0]->SetLineStyle(2);

	limitGraph[5]->SetLineStyle(1);
	limitGraph[5]->SetLineWidth(2);

	//DelphiGraph->SetLineStyle(3);
	DelphiGraph->SetLineStyle(1);
	DelphiGraph->SetLineWidth(3);
	DelphiGraph->SetLineColor(kBlue);
	
	TLegend* legend = new TLegend(0.6,0.3,0.9,0.6,NULL,"brNDC");
	//TLegend* legend = new TLegend(0.3,0.6,0.6,0.9,NULL,"brNDC");
	legend->SetFillStyle(0);
	//legend->AddEntry( (TObject*)nullptr, "95\% CL upper limits", "");
	legend->SetHeader("95\% CL upper limits");
	legend->AddEntry(limitGraph[0], "Expected", "l");
	legend->AddEntry(twoSigmaBand, "Expected \\pm 2 s.d.", "f");	
	legend->AddEntry(oneSigmaBand, "Expected \\pm 1 s.d.", "f");	
	legend->AddEntry(limitGraph[5], "Observed", "l");
	legend->AddEntry(DelphiGraph, "DELPHI", "l");

	const TString extra = "_cutAt5";	
	for(unsigned i = 0; i < 2; ++i){
		std::shared_ptr<TCanvas> c = std::make_shared<TCanvas>("canv","canv",500,500);
		c->SetTopMargin(0.08);
		c->SetLeftMargin(0.18);
		c->SetBottomMargin(0.14);
		c->SetLogy();
		twoSigmaBand->Draw("CAE3");
		oneSigmaBand->Draw("CE3same");
		twoSigmaBand->GetYaxis()->SetRangeUser(4*0.000001, 2);
		if(i == 0) twoSigmaBand->GetXaxis()->SetLimits(1, 1000);
		if(i == 1) twoSigmaBand->GetXaxis()->SetLimits(5, 1000);
		limitGraph[0]->Draw("Lsame");
		limitGraph[5]->Draw("Lsame");
		DelphiGraph->Draw("Lsame");
		legend->Draw("same");
		drawLumi(c.get());
	    c->SaveAs("plots/hnl_Limits" + ((i == 1) ? extra : "" ) + ".pdf");
		c->Close();
	}	

	//Draw logarithmic plot
	twoSigmaBand->GetYaxis()->SetTitleOffset(1.4);
	twoSigmaBand->GetXaxis()->SetTitleOffset(1);
	twoSigmaBand->GetXaxis()->SetTitle("m_{N} (GeV)");
	twoSigmaBand->GetYaxis()->SetTitle("|V_{eN}|^{2} + |V_{#muN}|^{2}");

	TLegend* legend2 = new TLegend(0.3,0.6,0.6,0.9,NULL,"brNDC");
	legend2->SetFillStyle(0);
	//legend->AddEntry( (TObject*)nullptr, "95\% CL upper limits", "");
	legend2->SetHeader("95\% CL upper limits");
	legend2->AddEntry(limitGraph[0], "Expected", "l");
	//legend2->AddEntry(twoSigmaBand, "Expected \\pm 2 std. deviation", "f");	
	//legend2->AddEntry(oneSigmaBand, "Expected \\pm 1 std. deviation", "f");	
	legend2->AddEntry(twoSigmaBand, "Expected \\pm 2 s.d.", "f");	
	legend2->AddEntry(oneSigmaBand, "Expected \\pm 1 s.d.", "f");
	legend2->AddEntry(limitGraph[5], "Observed", "l");
	legend2->AddEntry(DelphiGraph, "DELPHI", "l");
	
	for(unsigned i = 0; i < 2; ++i){
		std::shared_ptr<TCanvas> clog = std::make_shared<TCanvas>("canv","canv",500,500);
		clog->SetTopMargin(0.08);
		clog->SetLeftMargin(0.18);
		clog->SetBottomMargin(0.14);
		clog->SetLogy();
		clog->SetLogx();
		twoSigmaBand->GetYaxis()->SetRangeUser(4*0.000001, 2);
		if(i == 0) twoSigmaBand->GetXaxis()->SetLimits(1, 1000);
		if(i == 1) twoSigmaBand->GetXaxis()->SetLimits(5, 1000);
		twoSigmaBand->Draw("CAE3");
		oneSigmaBand->Draw("CE3same");
		limitGraph[0]->Draw("Lsame");
		limitGraph[5]->Draw("Lsame");
		DelphiGraph->Draw("Lsame");
		legend2->Draw("same");
		drawLumi(clog.get());
    	clog->SaveAs("plots/hnl_Limits_logx" + ((i == 1) ? extra : "" ) + ".pdf");
	}
}















