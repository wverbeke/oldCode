#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TROOT.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include <iostream>
#include <utility>      // std::pair
#include <memory>
//Include other parts of the code
#include "drawLumi.h"
#include "tdrstyle.h"



void plotLimitsMuon(){
	setTDRStyle();
	gROOT->SetBatch(kTRUE);
	const unsigned nPoints = 26; 
	const TString fileNames[nPoints] = {"higgsCombineMuonCouplingm1.Asymptotic.mH120.root", "higgsCombineMuonCouplingm2.Asymptotic.mH120.root", "higgsCombineMuonCouplingm3.Asymptotic.mH120.root", "higgsCombineMuonCouplingm4.Asymptotic.mH120.root", "higgsCombineMuonCouplingm5.Asymptotic.mH120.root", "higgsCombineMuonCouplingm6.Asymptotic.mH120.root", "higgsCombineMuonCouplingm7.Asymptotic.mH120.root",  "higgsCombineMuonCouplingm8.Asymptotic.mH120.root",  "higgsCombineMuonCouplingm9.Asymptotic.mH120.root", "higgsCombineMuonCouplingm10.Asymptotic.mH120.root", "higgsCombineMuonCouplingm11.Asymptotic.mH120.root", "higgsCombineMuonCouplingm12.Asymptotic.mH120.root", "higgsCombineMuonCouplingm20.Asymptotic.mH120.root", "higgsCombineMuonCouplingm30.Asymptotic.mH120.root", "higgsCombineMuonCouplingm40.Asymptotic.mH120.root", "higgsCombineMuonCouplingm50.Asymptotic.mH120.root", "higgsCombineMuonCouplingm60.Asymptotic.mH120.root", "higgsCombineMuonCouplingm80.Asymptotic.mH120.root", "higgsCombineMuonCouplingm100.Asymptotic.mH120.root", "higgsCombineMuonCouplingm130.Asymptotic.mH120.root", "higgsCombineMuonCouplingm150.Asymptotic.mH120.root", "higgsCombineMuonCouplingm200.Asymptotic.mH120.root", "higgsCombineMuonCouplingm400.Asymptotic.mH120.root", "higgsCombineMuonCouplingm600.Asymptotic.mH120.root", "higgsCombineMuonCouplingm800.Asymptotic.mH120.root", "higgsCombineMuonCouplingm1000.Asymptotic.mH120.root"};
	const float massPoints[nPoints] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 20, 30, 40, 50, 60, 80, 100, 130, 150, 200, 400, 600, 800, 1000};
	const float massCouplings[nPoints] = {0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001,  0.00001, 0.00001, 0.00001, 0.00001, 0.0001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
	/*
	const unsigned nPoints = 17; 
	const TString fileNames[nPoints] = {"higgsCombineMuonCouplingm5.Asymptotic.mH120.root", "higgsCombineMuonCouplingm6.Asymptotic.mH120.root",  "higgsCombineMuonCouplingm10.Asymptotic.mH120.root", "higgsCombineMuonCouplingm20.Asymptotic.mH120.root", "higgsCombineMuonCouplingm30.Asymptotic.mH120.root", "higgsCombineMuonCouplingm40.Asymptotic.mH120.root", "higgsCombineMuonCouplingm50.Asymptotic.mH120.root", "higgsCombineMuonCouplingm60.Asymptotic.mH120.root", "higgsCombineMuonCouplingm80.Asymptotic.mH120.root", "higgsCombineMuonCouplingm100.Asymptotic.mH120.root", "higgsCombineMuonCouplingm130.Asymptotic.mH120.root", "higgsCombineMuonCouplingm150.Asymptotic.mH120.root", "higgsCombineMuonCouplingm200.Asymptotic.mH120.root", "higgsCombineMuonCouplingm400.Asymptotic.mH120.root", "higgsCombineMuonCouplingm600.Asymptotic.mH120.root", "higgsCombineMuonCouplingm800.Asymptotic.mH120.root", "higgsCombineMuonCouplingm1000.Asymptotic.mH120.root"};
	const float massPoints[nPoints] = {5, 6, 10, 20, 30, 40, 50, 60, 80, 100, 130, 150, 200, 400, 600, 800, 1000};
	const float massCouplings[nPoints] = {0.00001, 0.00001, 0.00001, 0.00001,  0.00001, 0.00001, 0.00001, 0.00001, 0.0001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
	*/
	//Delphi limits
	double massPointsDelphi[38] = {0.260273972602739, 0.292237442922374, 0.342465753424657, 0.390410958904109, 0.447488584474885, 0.481735159817351, 0.534246575342466, 0.591324200913242, 0.625570776255707, 0.652968036529680, 0.694063926940639, 0.732876712328766, 0.776255707762557, 0.815068493150685, 0.860730593607306, 0.915525114155251, 0.972602739726027, 1.038812785388127, 1.098173515981735, 1.166666666666666, 1.230593607305935, 1.280821917808219, 1.331050228310502, 1.381278538812785, 1.456621004566209, 1.513698630136986, 1.586757990867580, 1.630136986301369, 1.666666666666666, 1.712328767123287, 1.751141552511415, 1.785388127853881, 1.815068493150684, 1.837899543378995, 1.856164383561643, 1.872146118721461, 1.888127853881278, 1.901826484018265};
	double limitsDelphi[38] = {-2.00000000000000, -2.19424460431654, -2.49640287769784, -2.78417266187050, -3.11510791366906, -3.33812949640287, -3.63309352517985, -3.95683453237410, -4.14388489208633, -4.30935251798561, -4.47482014388489, -4.59712230215827, -4.68345323741007, -4.72661870503597, -4.75539568345323, -4.74820143884892, -4.73381294964028, -4.72661870503597, -4.72661870503597,  -4.70503597122302, -4.68345323741007, -4.67625899280575, -4.66187050359712, -4.63309352517985, -4.60431654676259, -4.60431654676259, -4.60431654676259, -4.59712230215827, -4.56115107913669, -4.48920863309352, -4.38129496402877, -4.20143884892086, -3.98561151079136, -3.73381294964028, -3.50359712230215, -3.20863309352517, -2.92805755395683, -2.64748201438848};

	double massPointsAtlas[25] = {100.237932 ,110.246164 ,119.409396 ,140.781448 ,160.210482 ,179.085100 ,199.904443 ,221.280479 ,239.880451 ,257.926575 ,282.358730 ,298.739614 ,324.006808 ,345.941813 ,365.100470 ,386.202713 ,409.248828 ,427.297514 ,442.847060
 ,457.563558 ,471.447010 ,483.942430 ,491.994837 ,499.769752 ,499.769752};
	double limitsAtlas[25] = {0.003474 ,0.002831 ,0.003084 ,0.004687 ,0.006876 ,0.009806 ,0.013886 ,0.019115 ,0.025217 ,0.031881 ,0.043279 ,0.052810 ,0.068711 ,0.086890 ,0.106042 ,0.132205 ,0.167192 ,0.198336 ,0.230308 ,0.265536 ,0.303980 ,0.340654 ,0.368393
 ,0.395574 ,0.395574};


	for(unsigned d = 0; d < 38; ++d){
		massPointsDelphi[d] = pow(10.0, massPointsDelphi[d]);
		limitsDelphi[d] = pow(10.0, limitsDelphi[d]);
	}
	//Make Delphi limit graph
	TGraphAsymmErrors* DelphiGraph;
	DelphiGraph = new TGraphAsymmErrors(38, massPointsDelphi, limitsDelphi);
	
	//Make ATLAS limit graph
	TGraphAsymmErrors* AtlasGraph;
	AtlasGraph = new TGraphAsymmErrors(25, massPointsAtlas, limitsAtlas);

	//Read the limits and uncertainties for every point
	TFile* limitFiles[nPoints];
	for(unsigned p = 0; p < nPoints; ++p){
		limitFiles[p] = TFile::Open("../limits/september2017/" + fileNames[p]);
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
		limitGraph[i]->GetYaxis()->SetTitle("|V_{#muN}|^{2}");
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
		limitGraph[1]->GetPoint(p, xVal, yValDown);	
		limitGraph[2]->GetPoint(p, xVal, yValUp);
		std::cout << "for mass : " << massPoints[p]  << " exp = " <<  yValCenter << " + " << fabs(yValCenter - yValUp) << " - " << fabs(yValCenter - yValDown) << "   	obs = " << limits[5][p] << std::endl;
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

	DelphiGraph->SetLineStyle(1);
	DelphiGraph->SetLineWidth(3);
	DelphiGraph->SetLineColor(kBlue);

	AtlasGraph->SetLineStyle(1);
	AtlasGraph->SetLineWidth(3);
	AtlasGraph->SetLineColor(kRed);
	
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
	legend->AddEntry(AtlasGraph, "ATLAS", "l");

	const TString extra = "_cutAt5";	
	const TString extra2 = "_corrected";
	for(unsigned i = 0; i < 2; ++i){
		std::shared_ptr<TCanvas> c = std::make_shared<TCanvas>("canv","canv",500,500);
		c->SetLogy();
		c->SetTopMargin(0.08);
		c->SetLeftMargin(0.18);
		c->SetBottomMargin(0.14);
		//if(i == 1) twoSigmaBand->SetMinimum(5.);
		twoSigmaBand->GetYaxis()->SetRangeUser(4*0.000001, 2);
		if(i == 0) twoSigmaBand->GetXaxis()->SetLimits(1, 1000);
		if(i == 1) twoSigmaBand->GetXaxis()->SetLimits(5, 1000);
		twoSigmaBand->Draw("CAE3");
		oneSigmaBand->Draw("CE3same");
		limitGraph[0]->Draw("Lsame");
		limitGraph[5]->Draw("Lsame");
		DelphiGraph->Draw("Lsame");
		AtlasGraph->Draw("Lsame");
		legend->Draw("same");
		drawLumi(c.get());
		gPad->RedrawAxis();
		c->SaveAs("plots/hnl_Limits_MuonCoupling" + ((i == 1) ? extra : "" ) + ".pdf");
		c->Close();
	}

	//Draw logarithmic plot

	TLegend* legend2 = new TLegend(0.3,0.6,0.6,0.9,NULL,"brNDC");
	legend2->SetFillStyle(0);
	legend2->SetHeader("95\% CL upper limits");
	legend2->AddEntry(limitGraph[0], "Expected", "l");
	legend2->AddEntry(twoSigmaBand, "Expected \\pm 2 s.d.", "f");	
	legend2->AddEntry(oneSigmaBand, "Expected \\pm 1 s.d.", "f");	
	legend2->AddEntry(limitGraph[5], "Observed", "l");
	legend2->AddEntry(DelphiGraph, "DELPHI", "l");
	legend2->AddEntry(AtlasGraph, "ATLAS", "l");
	legend2->Draw("same");


	for(unsigned i = 0; i < 2; ++i){
		std::shared_ptr<TCanvas> clog = std::make_shared<TCanvas>("canvlog","canvlog",500,500);
		clog->SetTopMargin(0.08);
		clog->SetLeftMargin(0.18);
		clog->SetBottomMargin(0.14);
		clog->SetLogy();
		clog->SetLogx();
		twoSigmaBand->GetYaxis()->SetRangeUser(4*0.000001, 2);
		if(i == 0) twoSigmaBand->GetXaxis()->SetLimits(1, 1000);
		if(i == 1) twoSigmaBand->GetXaxis()->SetLimits(5, 1000);
		twoSigmaBand->GetYaxis()->SetTitleOffset(1.4);
		twoSigmaBand->GetXaxis()->SetTitleOffset(1);
		twoSigmaBand->GetXaxis()->SetTitle("m_{N} (GeV)");
		twoSigmaBand->GetYaxis()->SetTitle("|V_{#muN}|^{2}");
		twoSigmaBand->Draw("CAE3");
		oneSigmaBand->Draw("CE3same");
		limitGraph[0]->Draw("Lsame");
		limitGraph[5]->Draw("Lsame");
		DelphiGraph->Draw("Lsame");
		AtlasGraph->Draw("Lsame");
		legend2->Draw("same");
		gPad->RedrawAxis();
		drawLumi(clog.get());
		clog->SaveAs("plots/hnl_Limits_MuonCoupling_logx" + ((i == 1) ? extra : "" ) + ".pdf");
	}
}


void plotLimitsMuonCorrected(){
	setTDRStyle();
	gROOT->SetBatch(kTRUE);
	const unsigned nPoints = 26; 
	const TString fileNames[nPoints] = {"higgsCombineMuonCouplingm1.Asymptotic.mH120.root", "higgsCombineMuonCouplingm2.Asymptotic.mH120.root", "higgsCombineMuonCouplingm3.Asymptotic.mH120.root", "higgsCombineMuonCouplingm4.Asymptotic.mH120.root", "higgsCombineMuonCouplingm5.Asymptotic.mH120.root", "higgsCombineMuonCouplingm6.Asymptotic.mH120.root", "higgsCombineMuonCouplingm7.Asymptotic.mH120.root",  "higgsCombineMuonCouplingm8.Asymptotic.mH120.root",  "higgsCombineMuonCouplingm9.Asymptotic.mH120.root", "higgsCombineMuonCouplingm10.Asymptotic.mH120.root", "higgsCombineMuonCouplingm11.Asymptotic.mH120.root", "higgsCombineMuonCouplingm12.Asymptotic.mH120.root", "higgsCombineMuonCouplingm20.Asymptotic.mH120.root", "higgsCombineMuonCouplingm30.Asymptotic.mH120.root", "higgsCombineMuonCouplingm40.Asymptotic.mH120.root", "higgsCombineMuonCouplingm50.Asymptotic.mH120.root", "higgsCombineMuonCouplingm60.Asymptotic.mH120.root", "higgsCombineMuonCouplingm80.Asymptotic.mH120.root", "higgsCombineMuonCouplingm100.Asymptotic.mH120.root", "higgsCombineMuonCouplingm130.Asymptotic.mH120.root", "higgsCombineMuonCouplingm150.Asymptotic.mH120.root", "higgsCombineMuonCouplingm200.Asymptotic.mH120.root", "higgsCombineMuonCouplingm400.Asymptotic.mH120.root", "higgsCombineMuonCouplingm600.Asymptotic.mH120.root", "higgsCombineMuonCouplingm800.Asymptotic.mH120.root", "higgsCombineMuonCouplingm1000.Asymptotic.mH120.root"};
	const float massPoints[nPoints] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 20, 30, 40, 50, 60, 80, 100, 130, 150, 200, 400, 600, 800, 1000};
	const float massCouplings[nPoints] = {0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001,  0.00001, 0.00001, 0.00001, 0.00001, 0.0001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
	/*
	const unsigned nPoints = 17; 
	const TString fileNames[nPoints] = {"higgsCombineMuonCouplingm5.Asymptotic.mH120.root", "higgsCombineMuonCouplingm6.Asymptotic.mH120.root",  "higgsCombineMuonCouplingm10.Asymptotic.mH120.root", "higgsCombineMuonCouplingm20.Asymptotic.mH120.root", "higgsCombineMuonCouplingm30.Asymptotic.mH120.root", "higgsCombineMuonCouplingm40.Asymptotic.mH120.root", "higgsCombineMuonCouplingm50.Asymptotic.mH120.root", "higgsCombineMuonCouplingm60.Asymptotic.mH120.root", "higgsCombineMuonCouplingm80.Asymptotic.mH120.root", "higgsCombineMuonCouplingm100.Asymptotic.mH120.root", "higgsCombineMuonCouplingm130.Asymptotic.mH120.root", "higgsCombineMuonCouplingm150.Asymptotic.mH120.root", "higgsCombineMuonCouplingm200.Asymptotic.mH120.root", "higgsCombineMuonCouplingm400.Asymptotic.mH120.root", "higgsCombineMuonCouplingm600.Asymptotic.mH120.root", "higgsCombineMuonCouplingm800.Asymptotic.mH120.root", "higgsCombineMuonCouplingm1000.Asymptotic.mH120.root"};
	const float massPoints[nPoints] = {5, 6, 10, 20, 30, 40, 50, 60, 80, 100, 130, 150, 200, 400, 600, 800, 1000};
	const float massCouplings[nPoints] = {0.00001, 0.00001, 0.00001, 0.00001,  0.00001, 0.00001, 0.00001, 0.00001, 0.0001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
	*/
	//Delphi limits
	double massPointsDelphi[38] = {0.260273972602739, 0.292237442922374, 0.342465753424657, 0.390410958904109, 0.447488584474885, 0.481735159817351, 0.534246575342466, 0.591324200913242, 0.625570776255707, 0.652968036529680, 0.694063926940639, 0.732876712328766, 0.776255707762557, 0.815068493150685, 0.860730593607306, 0.915525114155251, 0.972602739726027, 1.038812785388127, 1.098173515981735, 1.166666666666666, 1.230593607305935, 1.280821917808219, 1.331050228310502, 1.381278538812785, 1.456621004566209, 1.513698630136986, 1.586757990867580, 1.630136986301369, 1.666666666666666, 1.712328767123287, 1.751141552511415, 1.785388127853881, 1.815068493150684, 1.837899543378995, 1.856164383561643, 1.872146118721461, 1.888127853881278, 1.901826484018265};
	double limitsDelphi[38] = {-2.00000000000000, -2.19424460431654, -2.49640287769784, -2.78417266187050, -3.11510791366906, -3.33812949640287, -3.63309352517985, -3.95683453237410, -4.14388489208633, -4.30935251798561, -4.47482014388489, -4.59712230215827, -4.68345323741007, -4.72661870503597, -4.75539568345323, -4.74820143884892, -4.73381294964028, -4.72661870503597, -4.72661870503597,  -4.70503597122302, -4.68345323741007, -4.67625899280575, -4.66187050359712, -4.63309352517985, -4.60431654676259, -4.60431654676259, -4.60431654676259, -4.59712230215827, -4.56115107913669, -4.48920863309352, -4.38129496402877, -4.20143884892086, -3.98561151079136, -3.73381294964028, -3.50359712230215, -3.20863309352517, -2.92805755395683, -2.64748201438848};

	double massPointsAtlas[25] = {100.237932 ,110.246164 ,119.409396 ,140.781448 ,160.210482 ,179.085100 ,199.904443 ,221.280479 ,239.880451 ,257.926575 ,282.358730 ,298.739614 ,324.006808 ,345.941813 ,365.100470 ,386.202713 ,409.248828 ,427.297514 ,442.847060
 ,457.563558 ,471.447010 ,483.942430 ,491.994837 ,499.769752 ,499.769752};
	double limitsAtlas[25] = {0.003474 ,0.002831 ,0.003084 ,0.004687 ,0.006876 ,0.009806 ,0.013886 ,0.019115 ,0.025217 ,0.031881 ,0.043279 ,0.052810 ,0.068711 ,0.086890 ,0.106042 ,0.132205 ,0.167192 ,0.198336 ,0.230308 ,0.265536 ,0.303980 ,0.340654 ,0.368393
 ,0.395574 ,0.395574};

	double massPointsCMS[16] = {40.000000 ,50.000000 ,60.000000 ,80.000000 ,90.000000 ,100.000000 ,125.000000 ,150.000000 ,175.000000 ,200.000000 ,250.000000 ,300.000000 ,350.000000 ,400.000000 ,500.000000 ,500.000000};
	double limitsCMS[16] = {0.000024 ,0.000026 ,0.000068 ,0.001254 ,0.005507 ,0.002278 ,0.001838 ,0.003759 ,0.007691 ,0.011818 ,0.029263 ,0.046048 ,0.077838 ,0.155493 ,0.620506 ,0.620506};

	//const double lifeTimCorrection[9] = {1./0.00335705, 1./0.011472, 1./0.025203, 1./0.0561662, 1./0.0808458, 1./0.135414, 1./0.227908, 1/0.275631, 1/0.326697};//1./0.162945};
	const double lifeTimCorrection[12] = {297.241, 87.1449, 39.5894, 17.7988, 10.9343, 7.31146, 5.416, 3.54383, 2.84184, 2.40831, 1.90825, 1.61451};//1./0.162945};


	for(unsigned d = 0; d < 38; ++d){
		massPointsDelphi[d] = pow(10.0, massPointsDelphi[d]);
		limitsDelphi[d] = pow(10.0, limitsDelphi[d]);
	}
	//Make Delphi limit graph
	TGraphAsymmErrors* DelphiGraph;
	DelphiGraph = new TGraphAsymmErrors(38, massPointsDelphi, limitsDelphi);
	
	//Make ATLAS limit graph
	TGraphAsymmErrors* AtlasGraph;
	AtlasGraph = new TGraphAsymmErrors(25, massPointsAtlas, limitsAtlas);

	//Make CMS limit graph
	TGraphAsymmErrors* cmsGraph;
	cmsGraph = new TGraphAsymmErrors(16, massPointsCMS, limitsCMS);

	//Read the limits and uncertainties for every point
	TFile* limitFiles[nPoints];
	for(unsigned p = 0; p < nPoints; ++p){
		limitFiles[p] = TFile::Open("../limits/september2017/" + fileNames[p]);
	}
	float limits[7][nPoints]; //0: central, 1: 65% down 2: 65% up, 2: 95% down, 3: 95% up
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
	//Copy observed limit to display it with 
	for(unsigned p = 0; p < nPoints; ++p){
		limits[6][p] = limits[5][p];
	}
	for(unsigned p = 0; p < nPoints; ++p){
		limits[6][p] *= massCouplings[p];
	}
	//Correct for couplings
	for(unsigned p = 0; p < nPoints; ++p){
		for(unsigned i = 0; i < 6; ++i){
			limits[i][p] *= massCouplings[p];
			if(p < 12) limits[i][p]*=lifeTimCorrection[p];
		}
	}

	//make graph for errors and expected limit
	TGraphAsymmErrors* limitGraph[7];
	for(unsigned i = 0; i < 7; ++i){
		limitGraph[i] = new TGraphAsymmErrors(nPoints, massPoints, limits[i]);
	}
    //h->GetXaxis()->SetTitleOffset(1.4);
	for(unsigned i = 0; i < 6; ++i){
		limitGraph[i]->GetYaxis()->SetTitleOffset(1.4);
		limitGraph[i]->GetXaxis()->SetTitleOffset(1);
		limitGraph[i]->GetXaxis()->SetTitle("m_{N} (GeV)");
		limitGraph[i]->GetYaxis()->SetTitle("|V_{#muN}|^{2}");
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
		limitGraph[1]->GetPoint(p, xVal, yValDown);	
		limitGraph[2]->GetPoint(p, xVal, yValUp);
		std::cout << "for mass : " << massPoints[p]  << " exp = " <<  yValCenter << " + " << fabs(yValCenter - yValUp) << " - " << fabs(yValCenter - yValDown) << "   	obs = " << limits[5][p] << std::endl;
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

	limitGraph[6]->SetLineStyle(3);
	limitGraph[6]->SetLineWidth(2);
	limitGraph[6]->SetLineColor(kBlack);

	DelphiGraph->SetLineStyle(1);
	DelphiGraph->SetLineWidth(3);
	DelphiGraph->SetLineColor(kBlue);

	AtlasGraph->SetLineStyle(1);
	AtlasGraph->SetLineWidth(3);
	AtlasGraph->SetLineColor(kRed);

	cmsGraph->SetLineStyle(1);
	cmsGraph->SetLineWidth(3);
	cmsGraph->SetLineColor(kMagenta);
	
	TLegend* legend = new TLegend(0.6,0.3,0.9,0.6,NULL,"brNDC");
	//TLegend* legend = new TLegend(0.3,0.6,0.6,0.9,NULL,"brNDC");
	legend->SetFillStyle(0);
	//legend->AddEntry( (TObject*)nullptr, "95\% CL upper limits", "");
	legend->SetHeader("95\% CL upper limits");
	legend->AddEntry(limitGraph[0], "Expected", "l");
	legend->AddEntry(twoSigmaBand, "Expected \\pm 2 #sigma", "f");	
	legend->AddEntry(oneSigmaBand, "Expected \\pm 1 #sigma", "f");	
	legend->AddEntry(limitGraph[5], "Observed", "l");
	legend->AddEntry(limitGraph[6], "Observed, prompt N decays", "l");
	legend->AddEntry(DelphiGraph, "DELPHI", "l");
	legend->AddEntry(AtlasGraph, "ATLAS", "l");
	legend->AddEntry(cmsGraph, "CMS 8 TeV", "l");

	const TString extra = "_cutAt5";	
	for(unsigned i = 0; i < 2; ++i){
		std::shared_ptr<TCanvas> c = std::make_shared<TCanvas>("canv","canv",500,500);
		c->SetLogy();
		c->SetTopMargin(0.08);
		c->SetLeftMargin(0.18);
		c->SetBottomMargin(0.14);
		//if(i == 1) twoSigmaBand->SetMinimum(5.);
		twoSigmaBand->GetYaxis()->SetRangeUser(4*0.000001, 2);
		if(i == 0) twoSigmaBand->GetXaxis()->SetLimits(1, 1000);
		if(i == 1) twoSigmaBand->GetXaxis()->SetLimits(5, 1000);
		twoSigmaBand->Draw("CAE3");
		oneSigmaBand->Draw("CE3same");
		limitGraph[0]->Draw("Lsame");
		limitGraph[6]->Draw("Lsame");
		limitGraph[5]->Draw("Lsame");
		DelphiGraph->Draw("Lsame");
		AtlasGraph->Draw("Lsame");
		cmsGraph->Draw("Lsame");
		legend->Draw("same");
		drawLumi(c.get());
		gPad->RedrawAxis();
		c->SaveAs("plots/hnl_Limits_MuonCoupling" + ((i == 1) ? extra : "" ) + "_corrected.pdf");
		c->Close();
	}

	//Draw logarithmic plot

	TLegend* legend2 = new TLegend(0.3,0.6,0.7,0.9,NULL,"brNDC");
	legend2->SetFillStyle(0);
	legend2->SetHeader("95\% CL upper limits");
	legend2->AddEntry(limitGraph[0], "Expected", "l");
	legend2->AddEntry(twoSigmaBand, "Expected \\pm 2 #sigma", "f");	
	legend2->AddEntry(oneSigmaBand, "Expected \\pm 1 #sigma", "f");	
	legend2->AddEntry(limitGraph[5], "Observed", "l");
	legend2->AddEntry(limitGraph[6], "Observed, prompt N decays", "l");
	legend2->AddEntry(DelphiGraph, "DELPHI", "l");
	legend2->AddEntry(AtlasGraph, "ATLAS", "l");
	legend2->AddEntry(cmsGraph, "CMS 8 TeV", "l");
	legend2->Draw("same");


	for(unsigned i = 0; i < 2; ++i){
		std::shared_ptr<TCanvas> clog = std::make_shared<TCanvas>("canvlog","canvlog",500,500);
		clog->SetTopMargin(0.08);
		clog->SetLeftMargin(0.18);
		clog->SetBottomMargin(0.14);
		clog->SetLogy();
		clog->SetLogx();
		twoSigmaBand->GetYaxis()->SetRangeUser(4*0.000001, 2);
		if(i == 0) twoSigmaBand->GetXaxis()->SetLimits(1, 1000);
		if(i == 1) twoSigmaBand->GetXaxis()->SetLimits(5, 1000);
		twoSigmaBand->GetYaxis()->SetTitleOffset(1.4);
		twoSigmaBand->GetXaxis()->SetTitleOffset(1);
		twoSigmaBand->GetXaxis()->SetTitle("m_{N} (GeV)");
		twoSigmaBand->GetYaxis()->SetTitle("|V_{#muN}|^{2}");
		twoSigmaBand->Draw("CAE3");
		oneSigmaBand->Draw("CE3same");
		limitGraph[0]->Draw("Lsame");
		limitGraph[6]->Draw("Lsame");
		limitGraph[5]->Draw("Lsame");
		DelphiGraph->Draw("Lsame");
		AtlasGraph->Draw("Lsame");
		cmsGraph->Draw("Lsame");
		legend2->Draw("same");
		gPad->RedrawAxis();
		drawLumi(clog.get());
		clog->SaveAs("plots/hnl_Limits_MuonCoupling_logx" + ((i == 1) ? extra : "" ) + "_corrected.pdf");
	}
}


void plotMuonLimits(){
	plotLimitsMuon();
	plotLimitsMuonCorrected();
}














