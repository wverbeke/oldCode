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
#include <memory>
using std::cout;
using std::endl;
//Include other parts of the code
#include "drawLumi.h"
#include "tdrstyle.h"
void plotLimitsElectron(){
	setTDRStyle();
	gROOT->SetBatch(kTRUE);

	const unsigned nPoints = 26; 
	const TString fileNames[nPoints] = {"higgsCombineElectronCouplingm1.Asymptotic.mH120.root", "higgsCombineElectronCouplingm2.Asymptotic.mH120.root", "higgsCombineElectronCouplingm3.Asymptotic.mH120.root", "higgsCombineElectronCouplingm4.Asymptotic.mH120.root", "higgsCombineElectronCouplingm5.Asymptotic.mH120.root", "higgsCombineElectronCouplingm6.Asymptotic.mH120.root", "higgsCombineElectronCouplingm7.Asymptotic.mH120.root", "higgsCombineElectronCouplingm8.Asymptotic.mH120.root", "higgsCombineElectronCouplingm9.Asymptotic.mH120.root", "higgsCombineElectronCouplingm10.Asymptotic.mH120.root", "higgsCombineElectronCouplingm11.Asymptotic.mH120.root", "higgsCombineElectronCouplingm12.Asymptotic.mH120.root", "higgsCombineElectronCouplingm20.Asymptotic.mH120.root", "higgsCombineElectronCouplingm30.Asymptotic.mH120.root", "higgsCombineElectronCouplingm40.Asymptotic.mH120.root", "higgsCombineElectronCouplingm50.Asymptotic.mH120.root", "higgsCombineElectronCouplingm60.Asymptotic.mH120.root", "higgsCombineElectronCouplingm80.Asymptotic.mH120.root", "higgsCombineElectronCouplingm100.Asymptotic.mH120.root", "higgsCombineElectronCouplingm130.Asymptotic.mH120.root", "higgsCombineElectronCouplingm150.Asymptotic.mH120.root", "higgsCombineElectronCouplingm200.Asymptotic.mH120.root", "higgsCombineElectronCouplingm400.Asymptotic.mH120.root", "higgsCombineElectronCouplingm600.Asymptotic.mH120.root", "higgsCombineElectronCouplingm800.Asymptotic.mH120.root", "higgsCombineElectronCouplingm1000.Asymptotic.mH120.root"};
	const float massPoints[nPoints] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 20, 30, 40, 50, 60, 80, 100, 130, 150, 200, 400, 600, 800, 1000};
	const float massCouplings[nPoints] = {0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001,  0.00001, 0.00001, 0.00001, 0.00001, 0.0001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
	/*
	const unsigned nPoints = 17; 
	const TString fileNames[nPoints] = {"higgsCombineElectronCouplingm5.Asymptotic.mH120.root", "higgsCombineElectronCouplingm6.Asymptotic.mH120.root", "higgsCombineElectronCouplingm10.Asymptotic.mH120.root", "higgsCombineElectronCouplingm20.Asymptotic.mH120.root", "higgsCombineElectronCouplingm30.Asymptotic.mH120.root", "higgsCombineElectronCouplingm40.Asymptotic.mH120.root", "higgsCombineElectronCouplingm50.Asymptotic.mH120.root", "higgsCombineElectronCouplingm60.Asymptotic.mH120.root", "higgsCombineElectronCouplingm80.Asymptotic.mH120.root", "higgsCombineElectronCouplingm100.Asymptotic.mH120.root", "higgsCombineElectronCouplingm130.Asymptotic.mH120.root", "higgsCombineElectronCouplingm150.Asymptotic.mH120.root", "higgsCombineElectronCouplingm200.Asymptotic.mH120.root", "higgsCombineElectronCouplingm400.Asymptotic.mH120.root", "higgsCombineElectronCouplingm600.Asymptotic.mH120.root", "higgsCombineElectronCouplingm800.Asymptotic.mH120.root", "higgsCombineElectronCouplingm1000.Asymptotic.mH120.root"};
	const float massPoints[nPoints] = {5, 6,10, 20, 30, 40, 50, 60, 80, 100, 130, 150, 200, 400, 600, 800, 1000};
	const float massCouplings[nPoints] = {0.00001,0.00001, 0.00001, 0.00001,  0.00001, 0.00001, 0.00001, 0.00001, 0.0001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
	*/
	//Read the limits and uncertainties for every point

	//Delphi limits
	double massPointsDelphi[38] = {0.260273972602739, 0.292237442922374, 0.342465753424657, 0.390410958904109, 0.447488584474885, 0.481735159817351, 0.534246575342466, 0.591324200913242, 0.625570776255707, 0.652968036529680, 0.694063926940639, 0.732876712328766, 0.776255707762557, 0.815068493150685, 0.860730593607306, 0.915525114155251, 0.972602739726027, 1.038812785388127, 1.098173515981735, 1.166666666666666, 1.230593607305935, 1.280821917808219, 1.331050228310502, 1.381278538812785, 1.456621004566209, 1.513698630136986, 1.586757990867580, 1.630136986301369, 1.666666666666666, 1.712328767123287, 1.751141552511415, 1.785388127853881, 1.815068493150684, 1.837899543378995, 1.856164383561643, 1.872146118721461, 1.888127853881278, 1.901826484018265};
	double limitsDelphi[38] = {-2.00000000000000, -2.19424460431654, -2.49640287769784, -2.78417266187050, -3.11510791366906, -3.33812949640287, -3.63309352517985, -3.95683453237410, -4.14388489208633, -4.30935251798561, -4.47482014388489, -4.59712230215827, -4.68345323741007, -4.72661870503597, -4.75539568345323, -4.74820143884892, -4.73381294964028, -4.72661870503597, -4.72661870503597,  -4.70503597122302, -4.68345323741007, -4.67625899280575, -4.66187050359712, -4.63309352517985, -4.60431654676259, -4.60431654676259, -4.60431654676259, -4.59712230215827, -4.56115107913669, -4.48920863309352, -4.38129496402877, -4.20143884892086, -3.98561151079136, -3.73381294964028, -3.50359712230215, -3.20863309352517, -2.92805755395683, -2.64748201438848};


	//ATLAS limits
	double massPointsAtlas[20] = {1.00455861406709e+002, 1.09809766214851e+002, 1.38752503310975e+002, 1.78574156906282e+002, 1.98797339797022e+002, 2.20886448743258e+002, 2.38308779169597e+002, 2.57287535385483e+002, 2.79378072245287e+002, 3.02092250356086e+002, 3.21384433600234e+002, 3.43787683531516e+002, 3.63702080113091e+002, 3.82682978199330e+002, 3.99174308979077e+002, 4.19400704675346e+002, 4.35581821232361e+002, 4.49584298662402e+002, 4.64520988544563e+002, 4.80079534785759e+002};
	double limitsAtlas[20] = {4.04226688202133e-002, 2.87321017020558e-002, 3.33225809095133e-002, 5.51508368004247e-002, 6.92830498039127e-002, 9.05397993095326e-002, 1.11977313431971e-001, 1.38473609575320e-001, 1.75314576403352e-001, 2.16733719311297e-001, 2.55570350392377e-001, 3.13465143657095e-001, 3.72556547861380e-001, 4.39325988221587e-001, 5.10019175442642e-001, 5.96620239359251e-001, 6.76374403247621e-001, 7.60870825023857e-001, 8.49106552700366e-001, 9.62659703648251e-001};

	//L3 limits
	double massPointsL3[225] = {80.736041 ,84.756345 ,85.365482 ,85.852792 ,86.340102 ,86.827411 ,88.532995 ,89.142132 ,89.629442 ,90.238579 ,90.604061 ,91.091371 ,92.553299 ,94.380711 ,94.624365 ,94.989848 ,96.451777 ,96.939086 ,97.548223 ,98.157360 ,99.497462
 ,100.228426 ,100.472081 ,100.837563 ,101.446701 ,102.421320 ,103.274112 ,104.005076 ,104.614213 ,105.223350 ,105.710660 ,106.197970 ,106.563452 ,106.928934 ,107.172589 ,107.659898 ,108.390863 ,108.878173 ,109.365482 ,109.730964
 ,110.096447 ,110.218274 ,110.705584 ,111.314721 ,111.680203 ,112.167513 ,112.654822 ,113.873096 ,114.604061 ,114.969543 ,115.456853 ,115.822335 ,116.187817 ,116.553299 ,116.918782 ,117.162437 ,117.406091 ,118.015228 ,118.502538
 ,119.111675 ,119.355330 ,119.598985 ,119.964467 ,120.208122 ,120.451777 ,121.182741 ,122.035533 ,123.010152 ,123.375635 ,124.106599 ,124.715736 ,125.081218 ,125.568528 ,126.055838 ,126.421320 ,126.664975 ,127.030457 ,127.517766 
,128.005076 ,128.248731 ,128.614213 ,129.223350 ,129.954315 ,130.685279 ,131.294416 ,132.147208 ,132.634518 ,133.487310 ,134.218274 ,134.827411 ,135.680203 ,137.385787 ,138.238579 ,138.604061 ,138.969543 ,139.335025 ,139.822335
,140.553299 ,141.162437 ,141.893401 ,142.868020 ,143.598985 ,144.208122 ,144.573604 ,144.939086 ,145.304569 ,145.548223 ,146.035533 ,146.522843 ,147.010152 ,147.984772 ,148.593909 ,149.446701 ,150.055838 ,150.543147 ,151.395939
 ,151.883249 ,152.614213 ,153.345178 ,153.954315 ,154.563452 ,155.294416 ,155.659898 ,156.025381 ,156.512690 ,156.878173 ,157.487310 ,158.340102 ,159.192893 ,159.680203 ,160.289340 ,160.898477 ,161.507614 ,162.116751 ,162.604061
 ,163.335025 ,163.822335 ,164.187817 ,164.675127 ,165.406091 ,166.258883 ,166.989848 ,167.964467 ,168.451777 ,168.695431 ,169.182741 ,169.791878 ,170.401015 ,170.888325 ,171.131980 ,171.497462 ,171.862944 ,172.593909 ,173.324873
 ,173.690355 ,174.299492 ,174.543147 ,174.908629 ,175.274112 ,175.761421 ,176.126904 ,176.370558 ,176.736041 ,177.345178 ,177.588832 ,178.197970 ,178.441624 ,178.685279 ,179.050761 ,179.538071 ,180.390863 ,180.756345 ,181.243655
 ,181.609137 ,182.218274 ,182.705584 ,182.949239 ,183.436548 ,184.411168 ,185.873096 ,186.360406 ,186.604061 ,186.969543 ,187.335025 ,187.578680 ,187.944162 ,188.553299 ,188.796954 ,189.406091 ,189.771574 ,190.380711 ,190.624365
 ,191.355330 ,191.842640 ,192.208122 ,193.426396 ,193.791878 ,194.644670 ,195.375635 ,195.741117 ,196.228426 ,196.715736 ,197.203046 ,198.177665 ,198.664975 ,199.274112 ,199.639594 ,200.248731 ,200.614213 ,200.857868 ,201.223350
 ,201.588832 ,201.832487 ,201.954315 ,202.197970 ,202.441624 ,202.563452 ,202.685279 ,202.807107 ,203.050761 ,203.172589 ,203.294416 ,203.416244 ,203.416244 ,203.416244};
	
	double limitsL3[225] = {0.001642 ,0.002485 ,0.002736 ,0.002757 ,0.002736 ,0.002656 ,0.002290 ,0.002207 ,0.002240 ,0.002376 ,0.002412 ,0.002341 ,0.002049 ,0.001741 ,0.001741 ,0.001767 ,0.001902 ,0.002079 ,0.002340 ,0.002411 ,0.002393 ,0.002222 ,0.002189
 ,0.002238 ,0.002464 ,0.002654 ,0.002733 ,0.002900 ,0.003009 ,0.002878 ,0.002753 ,0.002693 ,0.002576 ,0.002288 ,0.001944 ,0.001727 ,0.001752 ,0.001752 ,0.001701 ,0.001752 ,0.001873 ,0.002002 ,0.002140 ,0.002047 ,0.002047 ,0.002155
 ,0.002304 ,0.002338 ,0.002426 ,0.002593 ,0.002793 ,0.002731 ,0.002593 ,0.002517 ,0.002632 ,0.002834 ,0.003120 ,0.003645 ,0.004258 ,0.005048 ,0.005476 ,0.005517 ,0.005316 ,0.004974 ,0.004619 ,0.004353 ,0.004226 ,0.004352 ,0.004550 
,0.004652 ,0.004550 ,0.004320 ,0.004163 ,0.004352 ,0.004583 ,0.004617 ,0.004449 ,0.004101 ,0.003589 ,0.003165 ,0.002874 ,0.002689 ,0.002669 ,0.002811 ,0.003072 ,0.003187 ,0.003187 ,0.003432 ,0.003807 ,0.004069 ,0.004285 ,0.004753
 ,0.005978 ,0.006438 ,0.006296 ,0.005890 ,0.005676 ,0.005846 ,0.006157 ,0.006295 ,0.006156 ,0.005717 ,0.004966 ,0.004646 ,0.004929 ,0.005387 ,0.005716 ,0.005632 ,0.005347 ,0.004680 ,0.004036 ,0.003947 ,0.004065 ,0.004749 ,0.005426
 ,0.005113 ,0.005038 ,0.005267 ,0.005345 ,0.005713 ,0.006106 ,0.005798 ,0.005755 ,0.006478 ,0.007568 ,0.008973 ,0.011976 ,0.013481 ,0.013381 ,0.012990 ,0.012894 ,0.013183 ,0.013579 ,0.013478 ,0.013280 ,0.013477 ,0.015060 ,0.016704
 ,0.017987 ,0.017333 ,0.015976 ,0.014191 ,0.011619 ,0.010631 ,0.010552 ,0.011705 ,0.012511 ,0.012235 ,0.012981 ,0.014722 ,0.017326 ,0.018383 ,0.020240 ,0.023469 ,0.023818 ,0.021472 ,0.021157 ,0.021954 ,0.023121 ,0.025646 ,0.027616
 ,0.030181 ,0.030858 ,0.029736 ,0.029735 ,0.032019 ,0.034994 ,0.035779 ,0.035778 ,0.035251 ,0.035249 ,0.033468 ,0.029950 ,0.028648 ,0.032011 ,0.033964 ,0.034215 ,0.033712 ,0.031771 ,0.029067 ,0.028218 ,0.028218 ,0.030839 ,0.033703
 ,0.035759 ,0.035759 ,0.034972 ,0.035757 ,0.041156 ,0.042392 ,0.042078 ,0.043022 ,0.049887 ,0.051385 ,0.049884 ,0.043981 ,0.043980 ,0.053709 ,0.077195 ,0.085624 ,0.087544 ,0.090173 ,0.097821 ,0.114269 ,0.132503 ,0.159439 ,0.186256
 ,0.230856 ,0.261814 ,0.273703 ,0.277780 ,0.286123 ,0.324495 ,0.346852 ,0.396292 ,0.459534 ,0.513506 ,0.548885 ,0.595453 ,0.700779 ,0.812618 ,0.881562 ,0.894706 ,0.970622 ,0.970622};
	     
    for(unsigned d = 0; d < 38; ++d){
		massPointsDelphi[d] = pow(10.0, massPointsDelphi[d]);
		limitsDelphi[d] = pow(10.0, limitsDelphi[d]);
	}
	//Make Delphi limit graph
	TGraphAsymmErrors* DelphiGraph;
	DelphiGraph = new TGraphAsymmErrors(38, massPointsDelphi, limitsDelphi);

	//Make ATLAS graph
	TGraphAsymmErrors* AtlasGraph = new TGraphAsymmErrors(20, massPointsAtlas, limitsAtlas);
	//Make L3 graph
	TGraphAsymmErrors* L3Graph = new TGraphAsymmErrors(225, massPointsL3, limitsL3);
	
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
			
			//if(p == 0) limits[i][p] *= (1./0.37);
			//if(p == 1) limits[i][p] *= (1./0.515808);
			
			//if(p == 2) limits[i][p] *= (1./0.580613);
			//cout << limits[i][p] << endl;
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
		limitGraph[i]->GetYaxis()->SetTitle("|V_{eN}|^{2}");
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
		std::cout << "for mass : " << massPoints[p]  << " exp = " <<  yValCenter << "   	obs = " << limits[5][p] << std::endl;
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

	L3Graph->SetLineStyle(1);
	L3Graph->SetLineWidth(3);
	L3Graph->SetLineColor(kMagenta);

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
	legend->AddEntry(L3Graph, "L3", "l");
	legend->AddEntry(AtlasGraph, "ATLAS", "l");
	
	const TString extra = "_cutAt5";
	for(unsigned i = 0; i < 2; ++i){
		std::shared_ptr<TCanvas> c = std::make_shared<TCanvas>("canv","canv",500,500);
		c->SetLogy();
		c->SetTopMargin(0.08);
		c->SetLeftMargin(0.18);
		c->SetBottomMargin(0.14);
		twoSigmaBand->GetYaxis()->SetRangeUser(4*0.000001, 2);
		if(i == 0) twoSigmaBand->GetXaxis()->SetLimits(1, 1000);
		if(i == 1) twoSigmaBand->GetXaxis()->SetLimits(5, 1000);
		twoSigmaBand->Draw("CAE3");
		oneSigmaBand->Draw("CE3same");
		limitGraph[0]->Draw("Lsame");
		limitGraph[5]->Draw("Lsame");
		DelphiGraph->Draw("Lsame");
		AtlasGraph->Draw("Lsame");
		L3Graph->Draw("Lsame");
		legend->Draw("same");
		drawLumi(c.get());
		gPad->RedrawAxis();
		c->SaveAs("plots/hnl_Limits_ElectronCoupling" + ((i == 1) ? extra : "" ) + ".pdf");
		c->Close();
	}

	//Draw logarithmic plot
	twoSigmaBand->GetYaxis()->SetTitleOffset(1.4);
	twoSigmaBand->GetXaxis()->SetTitleOffset(1);
	twoSigmaBand->GetXaxis()->SetTitle("m_{N} (GeV)");
	twoSigmaBand->GetYaxis()->SetTitle("|V_{eN}|^{2}");
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
	legend2->AddEntry(L3Graph, "L3", "l");
	legend2->AddEntry(AtlasGraph, "ATLAS", "l");

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
		AtlasGraph->Draw("Lsame");
		L3Graph->Draw("Lsame");
		legend2->Draw("same");
		drawLumi(clog.get());
		gPad->RedrawAxis();
		clog->SaveAs("plots/hnl_Limits_ElectronCoupling_logx" + ((i == 1) ? extra : "" ) + ".pdf");
    }
}


void plotLimitsElectronCorrected(){
	setTDRStyle();
	gROOT->SetBatch(kTRUE);

	const unsigned nPoints = 26; 
	const TString fileNames[nPoints] = {"higgsCombineElectronCouplingm1.Asymptotic.mH120.root", "higgsCombineElectronCouplingm2.Asymptotic.mH120.root", "higgsCombineElectronCouplingm3.Asymptotic.mH120.root", "higgsCombineElectronCouplingm4.Asymptotic.mH120.root", "higgsCombineElectronCouplingm5.Asymptotic.mH120.root", "higgsCombineElectronCouplingm6.Asymptotic.mH120.root", "higgsCombineElectronCouplingm7.Asymptotic.mH120.root", "higgsCombineElectronCouplingm8.Asymptotic.mH120.root", "higgsCombineElectronCouplingm9.Asymptotic.mH120.root", "higgsCombineElectronCouplingm10.Asymptotic.mH120.root", "higgsCombineElectronCouplingm11.Asymptotic.mH120.root", "higgsCombineElectronCouplingm12.Asymptotic.mH120.root", "higgsCombineElectronCouplingm20.Asymptotic.mH120.root", "higgsCombineElectronCouplingm30.Asymptotic.mH120.root", "higgsCombineElectronCouplingm40.Asymptotic.mH120.root", "higgsCombineElectronCouplingm50.Asymptotic.mH120.root", "higgsCombineElectronCouplingm60.Asymptotic.mH120.root", "higgsCombineElectronCouplingm80.Asymptotic.mH120.root", "higgsCombineElectronCouplingm100.Asymptotic.mH120.root", "higgsCombineElectronCouplingm130.Asymptotic.mH120.root", "higgsCombineElectronCouplingm150.Asymptotic.mH120.root", "higgsCombineElectronCouplingm200.Asymptotic.mH120.root", "higgsCombineElectronCouplingm400.Asymptotic.mH120.root", "higgsCombineElectronCouplingm600.Asymptotic.mH120.root", "higgsCombineElectronCouplingm800.Asymptotic.mH120.root", "higgsCombineElectronCouplingm1000.Asymptotic.mH120.root"};
	const float massPoints[nPoints] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 20, 30, 40, 50, 60, 80, 100, 130, 150, 200, 400, 600, 800, 1000};
	const float massCouplings[nPoints] = {0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001,  0.00001, 0.00001, 0.00001, 0.00001, 0.0001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
	/*
	const unsigned nPoints = 17; 
	const TString fileNames[nPoints] = {"higgsCombineElectronCouplingm5.Asymptotic.mH120.root", "higgsCombineElectronCouplingm6.Asymptotic.mH120.root", "higgsCombineElectronCouplingm10.Asymptotic.mH120.root", "higgsCombineElectronCouplingm20.Asymptotic.mH120.root", "higgsCombineElectronCouplingm30.Asymptotic.mH120.root", "higgsCombineElectronCouplingm40.Asymptotic.mH120.root", "higgsCombineElectronCouplingm50.Asymptotic.mH120.root", "higgsCombineElectronCouplingm60.Asymptotic.mH120.root", "higgsCombineElectronCouplingm80.Asymptotic.mH120.root", "higgsCombineElectronCouplingm100.Asymptotic.mH120.root", "higgsCombineElectronCouplingm130.Asymptotic.mH120.root", "higgsCombineElectronCouplingm150.Asymptotic.mH120.root", "higgsCombineElectronCouplingm200.Asymptotic.mH120.root", "higgsCombineElectronCouplingm400.Asymptotic.mH120.root", "higgsCombineElectronCouplingm600.Asymptotic.mH120.root", "higgsCombineElectronCouplingm800.Asymptotic.mH120.root", "higgsCombineElectronCouplingm1000.Asymptotic.mH120.root"};
	const float massPoints[nPoints] = {5, 6,10, 20, 30, 40, 50, 60, 80, 100, 130, 150, 200, 400, 600, 800, 1000};
	const float massCouplings[nPoints] = {0.00001,0.00001, 0.00001, 0.00001,  0.00001, 0.00001, 0.00001, 0.00001, 0.0001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
	*/
	//Read the limits and uncertainties for every point

	//Delphi limits
	double massPointsDelphi[38] = {0.260273972602739, 0.292237442922374, 0.342465753424657, 0.390410958904109, 0.447488584474885, 0.481735159817351, 0.534246575342466, 0.591324200913242, 0.625570776255707, 0.652968036529680, 0.694063926940639, 0.732876712328766, 0.776255707762557, 0.815068493150685, 0.860730593607306, 0.915525114155251, 0.972602739726027, 1.038812785388127, 1.098173515981735, 1.166666666666666, 1.230593607305935, 1.280821917808219, 1.331050228310502, 1.381278538812785, 1.456621004566209, 1.513698630136986, 1.586757990867580, 1.630136986301369, 1.666666666666666, 1.712328767123287, 1.751141552511415, 1.785388127853881, 1.815068493150684, 1.837899543378995, 1.856164383561643, 1.872146118721461, 1.888127853881278, 1.901826484018265};
	double limitsDelphi[38] = {-2.00000000000000, -2.19424460431654, -2.49640287769784, -2.78417266187050, -3.11510791366906, -3.33812949640287, -3.63309352517985, -3.95683453237410, -4.14388489208633, -4.30935251798561, -4.47482014388489, -4.59712230215827, -4.68345323741007, -4.72661870503597, -4.75539568345323, -4.74820143884892, -4.73381294964028, -4.72661870503597, -4.72661870503597,  -4.70503597122302, -4.68345323741007, -4.67625899280575, -4.66187050359712, -4.63309352517985, -4.60431654676259, -4.60431654676259, -4.60431654676259, -4.59712230215827, -4.56115107913669, -4.48920863309352, -4.38129496402877, -4.20143884892086, -3.98561151079136, -3.73381294964028, -3.50359712230215, -3.20863309352517, -2.92805755395683, -2.64748201438848};


	//ATLAS limits
	double massPointsAtlas[20] = {1.00455861406709e+002, 1.09809766214851e+002, 1.38752503310975e+002, 1.78574156906282e+002, 1.98797339797022e+002, 2.20886448743258e+002, 2.38308779169597e+002, 2.57287535385483e+002, 2.79378072245287e+002, 3.02092250356086e+002, 3.21384433600234e+002, 3.43787683531516e+002, 3.63702080113091e+002, 3.82682978199330e+002, 3.99174308979077e+002, 4.19400704675346e+002, 4.35581821232361e+002, 4.49584298662402e+002, 4.64520988544563e+002, 4.80079534785759e+002};
	double limitsAtlas[20] = {4.04226688202133e-002, 2.87321017020558e-002, 3.33225809095133e-002, 5.51508368004247e-002, 6.92830498039127e-002, 9.05397993095326e-002, 1.11977313431971e-001, 1.38473609575320e-001, 1.75314576403352e-001, 2.16733719311297e-001, 2.55570350392377e-001, 3.13465143657095e-001, 3.72556547861380e-001, 4.39325988221587e-001, 5.10019175442642e-001, 5.96620239359251e-001, 6.76374403247621e-001, 7.60870825023857e-001, 8.49106552700366e-001, 9.62659703648251e-001};

	//L3 limits
	double massPointsL3[225] = {80.736041 ,84.756345 ,85.365482 ,85.852792 ,86.340102 ,86.827411 ,88.532995 ,89.142132 ,89.629442 ,90.238579 ,90.604061 ,91.091371 ,92.553299 ,94.380711 ,94.624365 ,94.989848 ,96.451777 ,96.939086 ,97.548223 ,98.157360 ,99.497462
 ,100.228426 ,100.472081 ,100.837563 ,101.446701 ,102.421320 ,103.274112 ,104.005076 ,104.614213 ,105.223350 ,105.710660 ,106.197970 ,106.563452 ,106.928934 ,107.172589 ,107.659898 ,108.390863 ,108.878173 ,109.365482 ,109.730964
 ,110.096447 ,110.218274 ,110.705584 ,111.314721 ,111.680203 ,112.167513 ,112.654822 ,113.873096 ,114.604061 ,114.969543 ,115.456853 ,115.822335 ,116.187817 ,116.553299 ,116.918782 ,117.162437 ,117.406091 ,118.015228 ,118.502538
 ,119.111675 ,119.355330 ,119.598985 ,119.964467 ,120.208122 ,120.451777 ,121.182741 ,122.035533 ,123.010152 ,123.375635 ,124.106599 ,124.715736 ,125.081218 ,125.568528 ,126.055838 ,126.421320 ,126.664975 ,127.030457 ,127.517766 
,128.005076 ,128.248731 ,128.614213 ,129.223350 ,129.954315 ,130.685279 ,131.294416 ,132.147208 ,132.634518 ,133.487310 ,134.218274 ,134.827411 ,135.680203 ,137.385787 ,138.238579 ,138.604061 ,138.969543 ,139.335025 ,139.822335
,140.553299 ,141.162437 ,141.893401 ,142.868020 ,143.598985 ,144.208122 ,144.573604 ,144.939086 ,145.304569 ,145.548223 ,146.035533 ,146.522843 ,147.010152 ,147.984772 ,148.593909 ,149.446701 ,150.055838 ,150.543147 ,151.395939
 ,151.883249 ,152.614213 ,153.345178 ,153.954315 ,154.563452 ,155.294416 ,155.659898 ,156.025381 ,156.512690 ,156.878173 ,157.487310 ,158.340102 ,159.192893 ,159.680203 ,160.289340 ,160.898477 ,161.507614 ,162.116751 ,162.604061
 ,163.335025 ,163.822335 ,164.187817 ,164.675127 ,165.406091 ,166.258883 ,166.989848 ,167.964467 ,168.451777 ,168.695431 ,169.182741 ,169.791878 ,170.401015 ,170.888325 ,171.131980 ,171.497462 ,171.862944 ,172.593909 ,173.324873
 ,173.690355 ,174.299492 ,174.543147 ,174.908629 ,175.274112 ,175.761421 ,176.126904 ,176.370558 ,176.736041 ,177.345178 ,177.588832 ,178.197970 ,178.441624 ,178.685279 ,179.050761 ,179.538071 ,180.390863 ,180.756345 ,181.243655
 ,181.609137 ,182.218274 ,182.705584 ,182.949239 ,183.436548 ,184.411168 ,185.873096 ,186.360406 ,186.604061 ,186.969543 ,187.335025 ,187.578680 ,187.944162 ,188.553299 ,188.796954 ,189.406091 ,189.771574 ,190.380711 ,190.624365
 ,191.355330 ,191.842640 ,192.208122 ,193.426396 ,193.791878 ,194.644670 ,195.375635 ,195.741117 ,196.228426 ,196.715736 ,197.203046 ,198.177665 ,198.664975 ,199.274112 ,199.639594 ,200.248731 ,200.614213 ,200.857868 ,201.223350
 ,201.588832 ,201.832487 ,201.954315 ,202.197970 ,202.441624 ,202.563452 ,202.685279 ,202.807107 ,203.050761 ,203.172589 ,203.294416 ,203.416244 ,203.416244 ,203.416244};
	
	double limitsL3[225] = {0.001642 ,0.002485 ,0.002736 ,0.002757 ,0.002736 ,0.002656 ,0.002290 ,0.002207 ,0.002240 ,0.002376 ,0.002412 ,0.002341 ,0.002049 ,0.001741 ,0.001741 ,0.001767 ,0.001902 ,0.002079 ,0.002340 ,0.002411 ,0.002393 ,0.002222 ,0.002189
 ,0.002238 ,0.002464 ,0.002654 ,0.002733 ,0.002900 ,0.003009 ,0.002878 ,0.002753 ,0.002693 ,0.002576 ,0.002288 ,0.001944 ,0.001727 ,0.001752 ,0.001752 ,0.001701 ,0.001752 ,0.001873 ,0.002002 ,0.002140 ,0.002047 ,0.002047 ,0.002155
 ,0.002304 ,0.002338 ,0.002426 ,0.002593 ,0.002793 ,0.002731 ,0.002593 ,0.002517 ,0.002632 ,0.002834 ,0.003120 ,0.003645 ,0.004258 ,0.005048 ,0.005476 ,0.005517 ,0.005316 ,0.004974 ,0.004619 ,0.004353 ,0.004226 ,0.004352 ,0.004550 
,0.004652 ,0.004550 ,0.004320 ,0.004163 ,0.004352 ,0.004583 ,0.004617 ,0.004449 ,0.004101 ,0.003589 ,0.003165 ,0.002874 ,0.002689 ,0.002669 ,0.002811 ,0.003072 ,0.003187 ,0.003187 ,0.003432 ,0.003807 ,0.004069 ,0.004285 ,0.004753
 ,0.005978 ,0.006438 ,0.006296 ,0.005890 ,0.005676 ,0.005846 ,0.006157 ,0.006295 ,0.006156 ,0.005717 ,0.004966 ,0.004646 ,0.004929 ,0.005387 ,0.005716 ,0.005632 ,0.005347 ,0.004680 ,0.004036 ,0.003947 ,0.004065 ,0.004749 ,0.005426
 ,0.005113 ,0.005038 ,0.005267 ,0.005345 ,0.005713 ,0.006106 ,0.005798 ,0.005755 ,0.006478 ,0.007568 ,0.008973 ,0.011976 ,0.013481 ,0.013381 ,0.012990 ,0.012894 ,0.013183 ,0.013579 ,0.013478 ,0.013280 ,0.013477 ,0.015060 ,0.016704
 ,0.017987 ,0.017333 ,0.015976 ,0.014191 ,0.011619 ,0.010631 ,0.010552 ,0.011705 ,0.012511 ,0.012235 ,0.012981 ,0.014722 ,0.017326 ,0.018383 ,0.020240 ,0.023469 ,0.023818 ,0.021472 ,0.021157 ,0.021954 ,0.023121 ,0.025646 ,0.027616
 ,0.030181 ,0.030858 ,0.029736 ,0.029735 ,0.032019 ,0.034994 ,0.035779 ,0.035778 ,0.035251 ,0.035249 ,0.033468 ,0.029950 ,0.028648 ,0.032011 ,0.033964 ,0.034215 ,0.033712 ,0.031771 ,0.029067 ,0.028218 ,0.028218 ,0.030839 ,0.033703
 ,0.035759 ,0.035759 ,0.034972 ,0.035757 ,0.041156 ,0.042392 ,0.042078 ,0.043022 ,0.049887 ,0.051385 ,0.049884 ,0.043981 ,0.043980 ,0.053709 ,0.077195 ,0.085624 ,0.087544 ,0.090173 ,0.097821 ,0.114269 ,0.132503 ,0.159439 ,0.186256
 ,0.230856 ,0.261814 ,0.273703 ,0.277780 ,0.286123 ,0.324495 ,0.346852 ,0.396292 ,0.459534 ,0.513506 ,0.548885 ,0.595453 ,0.700779 ,0.812618 ,0.881562 ,0.894706 ,0.970622 ,0.970622};

	//const double lifeTimeCorrection[7] = {1./0.00310683, 1./0.00916059, 1./0.0267024, 1./0.0525266,  1./0.0808458 ,1./0.127422 , 1./0.227908};
	const double lifeTimeCorrection[12] = {306.376, 109.141, 37.4016, 19.0212,  9.12216 ,7.91352 , 4.35153, 3.42325, 2.89377, 2.44636, 1.60709, 1.55998};
	//const double lifeTimeCorrection[12] = {297.241, 87.1449, 39.5894, 17.7988, 10.9343, 7.31146, 5.416, 3.54383, 2.84184, 2.40831, 1.90825, 1.61451};//1./0.162945};
	     
    for(unsigned d = 0; d < 38; ++d){
		massPointsDelphi[d] = pow(10.0, massPointsDelphi[d]);
		limitsDelphi[d] = pow(10.0, limitsDelphi[d]);
	}
	//Make Delphi limit graph
	TGraphAsymmErrors* DelphiGraph;
	DelphiGraph = new TGraphAsymmErrors(38, massPointsDelphi, limitsDelphi);

	//Make ATLAS graph
	TGraphAsymmErrors* AtlasGraph = new TGraphAsymmErrors(20, massPointsAtlas, limitsAtlas);
	//Make L3 graph
	TGraphAsymmErrors* L3Graph = new TGraphAsymmErrors(225, massPointsL3, limitsL3);
	
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
			if(p < 12) limits[i][p]*=lifeTimeCorrection[p];
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
		limitGraph[i]->GetYaxis()->SetTitle("|V_{eN}|^{2}");
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
		//std::cout << "limit for mass : " << massPoints[p]  << " is " <<  yValCenter << std::endl;
		std::cout << "for mass : " << massPoints[p]  << " exp = " <<  yValCenter << "   	obs = " << limits[5][p] << std::endl;
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

	L3Graph->SetLineStyle(1);
	L3Graph->SetLineWidth(3);
	L3Graph->SetLineColor(kMagenta);

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
	legend->AddEntry(L3Graph, "L3", "l");
	legend->AddEntry(AtlasGraph, "ATLAS", "l");
	
	const TString extra = "_cutAt5";
	for(unsigned i = 0; i < 2; ++i){
		std::shared_ptr<TCanvas> c = std::make_shared<TCanvas>("canv","canv",500,500);
		c->SetLogy();
		c->SetTopMargin(0.08);
		c->SetLeftMargin(0.18);
		c->SetBottomMargin(0.14);
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
		L3Graph->Draw("Lsame");
		legend->Draw("same");
		drawLumi(c.get());
		gPad->RedrawAxis();
		c->SaveAs("plots/hnl_Limits_ElectronCoupling" + ((i == 1) ? extra : "" ) + "_corrected.pdf");
		c->Close();
	}

	//Draw logarithmic plot
	twoSigmaBand->GetYaxis()->SetTitleOffset(1.4);
	twoSigmaBand->GetXaxis()->SetTitleOffset(1);
	twoSigmaBand->GetXaxis()->SetTitle("m_{N} (GeV)");
	twoSigmaBand->GetYaxis()->SetTitle("|V_{eN}|^{2}");
	TLegend* legend2 = new TLegend(0.3,0.6,0.7,0.9,NULL,"brNDC");
	legend2->SetFillStyle(0);
	//legend->AddEntry( (TObject*)nullptr, "95\% CL upper limits", "");
	legend2->SetHeader("95\% CL upper limits");
	legend2->AddEntry(limitGraph[0], "Expected", "l");
	//legend2->AddEntry(twoSigmaBand, "Expected \\pm 2 std. deviation", "f");	
	//legend2->AddEntry(oneSigmaBand, "Expected \\pm 1 std. deviation", "f");	
	legend2->AddEntry(twoSigmaBand, "Expected \\pm 2 #sigma", "f");	
	legend2->AddEntry(oneSigmaBand, "Expected \\pm 1 #sigma", "f");
	legend2->AddEntry(limitGraph[5], "Observed", "l");
	legend2->AddEntry(limitGraph[6], "Observed, prompt N decays", "l");
	legend2->AddEntry(DelphiGraph, "DELPHI", "l");
	legend2->AddEntry(L3Graph, "L3", "l");
	legend2->AddEntry(AtlasGraph, "ATLAS", "l");

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
		limitGraph[6]->Draw("Lsame");
		limitGraph[5]->Draw("Lsame");
		DelphiGraph->Draw("Lsame");
		AtlasGraph->Draw("Lsame");
		L3Graph->Draw("Lsame");
		legend2->Draw("same");
		drawLumi(clog.get());
		gPad->RedrawAxis();
		clog->SaveAs("plots/hnl_Limits_ElectronCoupling_logx" + ((i == 1) ? extra : "" ) + "_corrected.pdf");
    }
}



void plotElectronLimits(){
	plotLimitsElectron();
	plotLimitsElectronCorrected();
}















