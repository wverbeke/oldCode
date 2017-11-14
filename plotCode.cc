//include c++ library classes
#include <fstream>
#include <set>
#include <math.h> 
//include Root classes
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
//Include other parts of the code
//#include "MultilepSUSYfunc.h"
#include "plotCode.h"
#include "drawLumi.h"
#include "hnlTools.h"



extern const double xPad;
extern const Color_t colors[];
extern const Color_t sigCols[];


void histcol(TH1D *h, const Color_t color){
	h->SetLineColor(color);
	h->SetMarkerColor(color);
	h->SetLineWidth(2);
}

TH1D *HistDiv(TH1D *h1, TH1D *h2, const bool abs){
	TH1D *h1c = (TH1D*) h1->Clone();
	TH1D *h2c = (TH1D*) h2->Clone();
	if(!abs){
		h1c->Scale(1/h1c->Integral(), "width");
		h2c->Scale(1/h2c->Integral(), "width");
	}
	h1c->Divide(h2c);
	return h1c;
}

void HistLabelSizes(TH1D *h, const double xlabel, const double xtitle, const double ylabel, const double ytitle){
	h->GetXaxis()->SetLabelSize(xlabel);
	h->GetXaxis()->SetTitleSize(xtitle);
	h->GetYaxis()->SetLabelSize(ylabel);
	h->GetYaxis()->SetTitleSize(ytitle);
}

void HistLabelSizes(TH2D *h, const double xlabel, const double xtitle, const double ylabel, const double ytitle){
	h->GetXaxis()->SetLabelSize(xlabel);
	h->GetXaxis()->SetTitleSize(xtitle);
	h->GetYaxis()->SetLabelSize(ylabel);
	h->GetYaxis()->SetTitleSize(ytitle);
}

void StackCol(TH1D *h, const Color_t color){
	histcol(h,color);
	h->SetFillColor(color);
	h->SetLineWidth(1);
	h->SetLineColor(kBlack);
}

void yieldOrder(TH1D** hists, unsigned* histInd, const unsigned nHist){
	unsigned ordered[nHist];
	std::set<unsigned> usedHist;
	for(unsigned h = 0; h < nHist; ++h){
		double maxYield = -99999.;
		for(unsigned k = 0; k <nHist; ++k){
			if(usedHist.find(k) == usedHist.end()){
				double yield = hists[k]->GetSumOfWeights();
				if(yield > maxYield){
					maxYield = yield;
					ordered[h] = k;
				}
			}
		}
		usedHist.insert(ordered[h]);
	}
	TH1D* histC[nHist];
	for(unsigned h = 0; h < nHist; ++h){
		histC[h] = (TH1D*) hists[ordered[h]]->Clone();
		histInd[h] = ordered[h];
	}
	for(unsigned h = 0; h < nHist; ++h){
		hists[h] = (TH1D*) histC[h]->Clone();
	}
}

Color_t bkgColor_EWK(const TString& bkgName){
	if(bkgName == "non-prompt") return kAzure + 1;
	else if(bkgName == "WZ") return kOrange;
	else if(bkgName == "ZZ/H") return  kGreen + 1;
	else if(bkgName == "TT/T + X") return  kViolet-3;
	else if(bkgName == "triboson") return kRed + 1;
	else if(bkgName == "X + #gamma") return kOrange + 7;
	else if(bkgName == "T + X") return kCyan + 1;
}

Color_t bkgColor_HNL(const TString& bkgName){
	if(bkgName == "non-prompt") return kAzure + 1;
	else if(bkgName == "nonprompt") return kAzure + 1;
	else if(bkgName == "WZ") return kRed - 7;
	else if(bkgName == "X + #gamma") return  kGreen + 1; //+ 1// -7 // -9
	else if(bkgName == "X#gamma^{(*)}") return  kGreen + 1; //+ 1// -7 // -9
	else if(bkgName == "TT/T + X" || bkgName == "t#bar{t}/t + X") return kMagenta -7;
	else if(bkgName == "ZZ/H") return kOrange + 6;
	//else if(bkgName == "T + X") return kCyan + 1; 
	else if(bkgName == "triboson") return kBlue - 7;

	//For MC fakes
	else if(bkgName == "Drell-Yan") return kAzure + 1;
	else if(bkgName == "TT") return kCyan + 1;
	else return kBlack;
}

/*
Color_t bkgColor_HNL(const TString& bkgName){
	if(bkgName == "non-prompt") return kOrange - 4;
	else if(bkgName == "nonprompt") return kOrange - 4;
	else if(bkgName == "WZ") return kAzure - 8;// kBlue - 6;
	else if(bkgName == "X + #gamma") return  kBlue - 10; // -9
	else if(bkgName == "X#gamma^{(*)}") return  kBlue - 10;  // -9
	else if(bkgName == "TT/T + X" || bkgName == "t#bar{t}/t + X") return kAzure + 10; //kCyan + 1//-3
	else if(bkgName == "ZZ/H") return kPink - 5;
	//else if(bkgName == "T + X") return kCyan + 1; 
	else if(bkgName == "triboson") return kMagenta - 10; //kPink + 1
	//For MC fakes
	else if(bkgName == "Drell-Yan") return kAzure + 1;
	else if(bkgName == "TT") return kCyan + 1;
	else return kBlack;
}
*/
/*
bool inBkgList(const TString& bkgname){
	if(bkgName == "non-prompt") return true;
	else if(bkgName == "WZ") return true;
	else if(bkgName == "X + #gamma") return true;
	else if(bkgName == "TT + X") return true;
	else if(bkgName == "ZZ") return true;
	else if(bkgName == "T + X") return true;
	else if(bkgName == "triboson") return true;
	return false;
}
*/
void plotDataVSMC(TH1D* data, TH1D** bkg, const TString* names, const unsigned nHist, const TString& file, const bool ylog,  const unsigned widthopt, const TString& analysis, TH1D** bkgSyst, const bool plotsig, TH1D** signal , const TString* signames, const unsigned nSig, const bool signorm){
	//Order background histograms in terms of yields 
	unsigned histI[nHist];
	TH1D* bkgC[nHist];
	for(unsigned b = 0; b < nHist; ++b){
		bkgC[b] = (TH1D*) bkg[b]->Clone();
	}
	yieldOrder(bkgC, histI, nHist);
	//Stack containing all background histograms
	THStack* bkgStack = new THStack("bkgStack", "bkgStack");
	for(int effsam = nHist -1; effsam > -1; --effsam){
		if(analysis == "EWKino"){
			StackCol(bkgC[effsam], bkgColor_EWK(names[histI[effsam] + 1]) );
		} else if(analysis == "HNL"){
			StackCol(bkgC[effsam], bkgColor_HNL(names[histI[effsam] + 1]) );
		} else{
			StackCol(bkgC[effsam], colors[histI[effsam]]);
		}
		bkgStack->Add(bkgC[effsam], "f");
	}
	//Total background
	TH1D* bkgTot = (TH1D*) bkgC[0]->Clone();
	for(unsigned i = 1; i <  nHist; ++i){
		bkgTot->Add(bkgC[i]);
	}
	//Background histograms with full uncertainty
	TH1D* bkgE[nHist];
	for(unsigned b = 0; b < nHist; ++b){
		bkgE[b] = (TH1D*) bkgC[b]->Clone();
		for(unsigned bin = 1; bin < bkgE[b]->GetNbinsX() + 1; ++bin){
			bkgE[b]->SetBinError(bin, sqrt(bkgE[b]->GetBinError(bin)*bkgE[b]->GetBinError(bin) + bkgSyst[histI[b]]->GetBinContent(bin)*bkgSyst[histI[b]]->GetBinContent(bin)) );
		}
	}
	//Total background with full unc.
	TH1D* bkgTotE = (TH1D*) bkgE[0]->Clone();
	for(unsigned bkg = 1; bkg < nHist; ++bkg){
		bkgTotE->Add(bkgE[bkg]);
	}

	//Set Poisonian errors to data
	data->SetBinErrorOption(TH1::kPoisson);
	//Replace data by TGRaphAsymmErrors for plotting
	TGraphAsymmErrors* obs = new TGraphAsymmErrors(data);
	for(unsigned b = 1; b < data->GetNbinsX() + 1; ++b){
		obs->SetPointError(b - 1, 0, 0, data->GetBinErrorLow(b), (data->GetBinContent(b) == 0 ) ? 0 : data->GetBinErrorUp(b) );
	}

	//Legend for data and all backgrounds
	//TLegend* legend = new TLegend(0.2,0.8,0.95,0.9,NULL,"brNDC");
	TLegend* legend;
	if(!plotsig) legend = new TLegend(0.2,0.8,0.95,0.9,NULL,"brNDC");
	else legend = new TLegend(0.2,0.7,0.95,0.9,NULL,"brNDC");
	if(!plotsig) legend-> SetNColumns(4); //2
	else legend->SetNColumns(3);
	//Avoid legend box
	legend->SetFillStyle(0);
	//Add data to legend
	//legend->AddEntry(data,names[0]);
	legend->AddEntry(obs,names[0], "pe1");
	//Add overlaid signals to legend
	if(plotsig){
		for(unsigned sig = 0; sig < nSig; ++sig){
			histcol(signal[sig], sigCols[sig]);
			signal[sig]->SetLineWidth(3);
			legend->AddEntry(signal[sig], signames[sig]);
		}
	}
	//Add backgrounds to the legend
	for(int effsam = nHist - 1; effsam > -1; --effsam){
    	legend->AddEntry(bkgE[effsam], names[histI[effsam] + 1], "f");
    }
	//Determine canvas size, depending on chosen option
	double width, height;
	if(widthopt == 0){
		width = 600*(1 - xPad);
		height = 600; //on request of Lesya
	} else if(widthopt == 1){
		width = 1200;
		height = 500;
	} else if(widthopt == 2){
		width = 300;
		height = 500;
	} else{
		std::cerr << "Incorrect width option given can't make plot" << std::endl;
		return;
	}
	//Make canvas to plot
	TCanvas *c =  new TCanvas(file,"",width,height);
    c->cd();
	//Make upper pad to draw main plot and lower pad for ratios
	TPad* p1, *p2;
	//Plot data and MC yields in first pad
    p1 = new TPad(file,"",0,xPad,1,1);
    p1->Draw();
    p1->cd();
    //p1->SetBottomMargin(0);
	p1->SetBottomMargin(0.03);
	//Make the total background uncertainty visible as a grey band
    bkgTotE->SetFillStyle(3244); //3005  3244
    bkgTotE->SetFillColor(kGray+2);
    bkgTotE->SetMarkerStyle(0); //1
	legend->AddEntry(bkgTotE, "total bkg. unc.", "f");
	//Set minimum slightly above 0 to avoid chopped off zero in the plots
  	if(!ylog) data->SetMinimum(0.0001);
    else if(ylog) p1->SetLogy();
	if(!ylog) bkgTotE->SetMinimum(0);
	bkgTotE->GetXaxis()->SetLabelSize(0);

	//Determine the maximum range of the histogram, depending on the maximum range of the bkg or data
	double dataMax = data->GetBinContent(data->GetMaximumBin()) + data->GetBinError(data->GetMaximumBin());
	double bkgMax = bkgTot->GetBinContent(bkgTot->GetMaximumBin()) + bkgTotE->GetBinError(bkgTot->GetMaximumBin());
	double totalMax = std::max(dataMax, bkgMax);
	if(plotsig && !signorm){
		double sigMax = 0;
		for(unsigned sig = 0; sig < nSig; ++sig){
			if(signal[sig]->GetSumOfWeights() == 0) continue;
			double localMax = signal[sig]->GetBinContent(signal[sig]->GetMaximumBin()) + signal[sig]->GetBinError(signal[sig]->GetMaximumBin());
			if( localMax > sigMax) sigMax = localMax;
		}
		totalMax = std::max(totalMax, sigMax);
	}
    //if(!ylog) data->SetMaximum(totalMax*1.3);
	//if(!plotsig && !ylog) data->SetMaximum(totalMax*1.3);
	//else if(!ylog) data->SetMaximum(totalMax*1.6);

	//Hack not to draw 0 points
	for(unsigned b = 1; b < data->GetNbinsX() + 1; ++b){
		if(obs->GetY()[b - 1] == 0)  obs->GetY()[b - 1] += totalMax*10;
	}
	if(!plotsig && !ylog) bkgTotE->SetMaximum(totalMax*1.3);
	else if(!ylog) bkgTotE->SetMaximum(totalMax*1.6);
    else if (!plotsig) {
		/*
		//NEW TEST CODE
		double minimum = totalMax;
		for(unsigned back = 0; back < nHist; ++back){
			for(unsigned b = 1; b < bkg[back]->GetNbinsX() + 1; ++b){
				if(bkg[back]->GetBinContent(b) != 0 &&  bkg[back]->GetBinContent(b) < minimum){
					minimum = bkg[back]->GetBinContent(b);
				}
			}
		}
		if(0.5*bkgTot->GetBinContent(bkgTot->GetMinimumBin()) > minimum*30) data->SetMinimum(minimum*30);
		else if(0.5*bkgTot->GetBinContent(bkgTot->GetMinimumBin()) < minimum*30) data->SetMinimum(0.5*bkgTot->GetBinContent(bkgTot->GetMinimumBin()) );
		*/
		double minimum = totalMax;
		for(unsigned b = 1; b < bkgTot->GetNbinsX() + 1; ++b){
			if(bkgTot->GetBinContent(b) != 0 && bkgTot->GetBinContent(b) < minimum){
				minimum = bkgTot->GetBinContent(b);
			}
		}
		//data->SetMinimum(minimum/5.); //10.
 		double SF = log10(std::max(10., totalMax/minimum));
		//data->SetMaximum(totalMax*3*SF); //3
		bkgTotE->SetMinimum(minimum/5.);
		bkgTotE->SetMaximum(totalMax*6*SF); //used to be 3
		//Hack not to draw 0 points
		for(unsigned b = 1; b < data->GetNbinsX() + 1; ++b){
			if(data->GetBinContent(b) == 0)  obs->GetY()[b - 1] += totalMax*30*SF;
		}
	}			
	else{
		/*
		double minimum = totalMax;
		for(unsigned back = 0; back < nHist; ++back){
			for(unsigned b = 1; b < bkg[back]->GetNbinsX() + 1; ++b){
				if(bkg[back]->GetBinContent(b) != 0 &&  bkg[back]->GetBinContent(b) < minimum){
					minimum = bkg[back]->GetBinContent(b);
				}
			}
		}
		if(0.5*bkgTot->GetBinContent(bkgTot->GetMinimumBin()) > minimum*30) data->SetMinimum(minimum*30);
		else if(0.5*bkgTot->GetBinContent(bkgTot->GetMinimumBin()) < minimum*30) data->SetMinimum(0.5*bkgTot->GetBinContent(bkgTot->GetMinimumBin()) );
		//double SF = log10(std::max(10., totalMax/bkgTot->GetBinContent(bkgTot->GetMinimumBin()) ));
		double SF = log10(std::max(10., totalMax/minimum));
		data->SetMaximum(totalMax*6*SF);
		*/
		double minimum = totalMax;
		for(unsigned b = 1; b < bkgTot->GetNbinsX() + 1; ++b){
			if(bkgTot->GetBinContent(b) != 0 && bkgTot->GetBinContent(b) < minimum){
				minimum = bkgTot->GetBinContent(b);
			}
		}
		//data->SetMinimum(minimum/5.); //10.
		double SF = log10(std::max(10., totalMax/minimum ));
		//data->SetMaximum(totalMax*6*SF);
		bkgTotE->SetMinimum(minimum/5.);
		bkgTotE->SetMaximum(totalMax*6*SF);
		//Hack not to draw 0 points
		for(unsigned b = 1; b < data->GetNbinsX() + 1; ++b){
			if(data->GetBinContent(b) == 0)  obs->GetY()[b - 1] += totalMax*30*SF;
		}
	}

	//Draw histograms and legends
	//First draw data to fix plot range
    //data->Draw("pe");
    bkgTotE->Draw("e2"); //e2same
    bkgStack->Draw("hist same");
	//Redraw data so it is overlaid on the background stack
    //data->Draw("pe same");
    legend->Draw("same");
    bkgTotE->Draw("e2same"); //e2same
	//data->Draw("pe1 same"); // NEW  //pesame
	obs->Draw("pe1 same");	
	//Draw signal plots
	if(plotsig){
		for(unsigned sig = 0; sig < nSig; ++sig){
			if(signorm && signal[sig]->GetSumOfWeights() != 0) signal[sig]->Scale(bkgTot->GetSumOfWeights()/ signal[sig]->GetSumOfWeights());
			signal[sig]->Draw("histsame");
		}
	}
	//redraw axis over histograms
    gPad->RedrawAxis();
	//Draw CMS header
	drawLumi(p1);
    c->cd(); 
	//Make ratio plot in second pad
	
	const unsigned nBins = data->GetNbinsX();
	//TH1D* dataErrors = new TH1D("dataerrors" + file, "dataerrors" + file, nBins, data->GetBinLowEdge(1), data->GetBinLowEdge(data->GetNbinsX()) + data->GetBinWidth(data->GetNbinsX()));
	TH1D* bkgStatErrors = new TH1D("bkgStaterrors" + file, "bkgStaterrors" + file, nBins, data->GetBinLowEdge(1), data->GetBinLowEdge(data->GetNbinsX()) + data->GetBinWidth(data->GetNbinsX()));
	//TH1D* bkgStatErrors = new TH1D("bkgStaterros" + file, "bkgStaterrors" + file, nBins, data->GetBinLowEdge(1), data->GetBinLowEdge(data->GetNbinsX()) + data->GetBinWidth(data->GetNbinsX()));
	for(unsigned b = 1; b < nBins + 1; ++b){
		/*
		if(data->GetBinContent(b) != 0){
			if(bkgTot->GetBinContent(b) != 0) dataErrors->SetBinContent(b, data->GetBinError(b)/bkgTot->GetBinContent(b));
			else dataErrors->SetBinContent(b, data->GetBinError(b)/data->GetBinContent(b));
		} else{
			dataErrors->SetBinContent(b, 0);
		}
		*/
		bkgStatErrors->SetBinContent(b, 1.);
		//bkgStatErrors->SetBinContent(b, 1.);
		if(bkgTot->GetBinContent(b) != 0){
			bkgStatErrors->SetBinError(b, bkgTot->GetBinError(b)/bkgTot->GetBinContent(b));
			//bkgStatErrors->SetBinError(b, sqrt( *bkgTot->GetBinContent(b) ));
		} else{
			bkgStatErrors->SetBinError(b, 0.);
			//bkgStatErrors->SetBinError(b, 0.);
		}			
	}
	TH1D* bkgErrors = (TH1D*) bkgTotE->Clone();//new TH1D("bkgerrors" + file, "bkgerrors" + file, nBins, data->GetBinLowEdge(1), data->GetBinLowEdge(data->GetNbinsX()) + data->GetBinWidth(data->GetNbinsX()));
	/*
	if(analysis == "HNL"){
		TString bkgNames[nHist];
		for(int b = nHist - 1; b > -1; --b){
    		bkgNames[b] = names[histI[b] + 1];
    	}
		hnl::setSystUnc(bkg, nHist - 1, bkgNames);
		bkgTot = (TH1D*) bkg[0]->Clone();
		for(int i = 1; i <  nHist; ++i){
			bkgTot->Add(bkg[i]);
		}
	*/
	for(unsigned b = 1; b < nBins + 1; ++b){
		bkgErrors->SetBinContent(b, 1.);
		//bkgStatErrors->SetBinContent(b, 1.);
		if(bkgTot->GetBinContent(b) != 0){
			bkgErrors->SetBinError(b, bkgTotE->GetBinError(b)/bkgTotE->GetBinContent(b));
			//bkgStatErrors->SetBinError(b, sqrt( *bkgTot->GetBinContent(b) ));
		} else{
			bkgErrors->SetBinError(b, 0.);
			//bkgStatErrors->SetBinError(b, 0.);
		}
	}			
	//}

	bkgStatErrors->SetFillStyle(1001);
	bkgErrors->SetFillStyle(1001);
	bkgStatErrors->SetFillColor(kCyan  - 4); //
	bkgErrors->SetFillColor(kOrange	- 4); //kOrange  //kOrange + 7 kYellow - 3   kOrange - 4
	//bkgStatErrors->SetFillColor(kOrange); 
    bkgStatErrors->SetMarkerStyle(1);
	bkgErrors->SetMarkerStyle(1);
	//bkgStatErrors->SetMarkerStyle(1);
		
    p2 = new TPad(file + "2","",0,0.0,1,xPad);
    p2->Draw();
    p2->cd();
    //p2->SetTopMargin(0);
	p2->SetTopMargin(0.01);
    p2->SetBottomMargin(0.4);
	
	//TH1D* dataC = (TH1D*) data->Clone();
	//TH1D* bkgTotC = (TH1D*) bkgTot->Clone();

  // dataC->Divide(bkgTotE);

	TGraphAsymmErrors* obsRatio = new TGraphAsymmErrors(data);
	for(unsigned b = 1; b < data->GetNbinsX() + 1; ++b){
		obsRatio->GetY()[b - 1] *= 1./bkgTotE->GetBinContent(b);
		obsRatio->SetPointError(b - 1, 0, 0, data->GetBinErrorLow(b)/bkgTotE->GetBinContent(b), data->GetBinErrorUp(b)/bkgTotE->GetBinContent(b));
		if(data->GetBinContent(b) == 0) obsRatio->GetY()[b - 1] += 5;
	}

	//Now reset bin errors 
	//for(unsigned b = 1; b < nBins + 1; ++b){
		//dataC->SetBinError(b, dataErrors->GetBinContent(b));
	//}
	//Legend for uncertainties
	TLegend* legend2 = new TLegend(0.18,0.85,0.8,0.98,NULL,"brNDC"); //0.18, 0.85, 0.65, 0.98
	legend2-> SetNColumns(3); //2
	//Avoid legend box
	legend2->SetFillStyle(0);
	//Add data to legend
	legend2->AddEntry(bkgStatErrors, "stat. pred. unc.", "f");
	legend2->AddEntry(bkgErrors, "total pred. unc.", "f");
	//legend2->AddEntry(dataC, "obs./pred. with total unc.", "pe12");
	legend2->AddEntry(obsRatio, "obs./pred.", "pe12");


	//legend2->AddEntry(bkgStatErrors, "stat. bkg. unc.", "f");
	/*
    dataC->SetMarkerColor(1);
    dataC->SetLineColor(1);
    dataC->GetYaxis()->SetRangeUser(0.,1.999);
    dataC->GetYaxis()->SetTitle("obs./pred.");
    dataC->GetYaxis()->SetTitleOffset(1.25/((1.-xPad)/xPad));
    dataC->GetYaxis()->SetTitleSize((1.-xPad)/xPad*0.06); //originally 0.06
    dataC->GetXaxis()->SetTitleSize((1.-xPad)/xPad*0.06); //originally 0.09
    dataC->GetYaxis()->SetLabelSize((1.-xPad)/xPad*0.05); //originally 0.05
    dataC->GetXaxis()->SetLabelSize((1.-xPad)/xPad*0.05); //originally 0.05
	*/
	

	bkgErrors->SetMarkerColor(1);
    bkgErrors->SetLineColor(1);
    bkgErrors->GetYaxis()->SetRangeUser(0.,1.999);
    bkgErrors->GetYaxis()->SetTitle("obs./pred.");
    bkgErrors->GetYaxis()->SetTitleOffset(1.25/((1.-xPad)/xPad)); //1.25
    bkgErrors->GetYaxis()->SetTitleSize((1.-xPad)/xPad*0.06); //originally 0.06
    bkgErrors->GetXaxis()->SetTitleSize((1.-xPad)/xPad*0.06); //originally 0.09
    bkgErrors->GetYaxis()->SetLabelSize((1.-xPad)/xPad*0.05); //originally 0.05
    bkgErrors->GetXaxis()->SetLabelSize((1.-xPad)/xPad*0.05); //originally 0.05
	//NeW 
	bkgErrors->GetXaxis()->SetLabelOffset((1-xPad)/xPad*0.009);
	bkgErrors->GetXaxis()->SetTitleOffset((1-xPad)/xPad*0.32);

    //dataC->Draw("pe");
	bkgErrors->Draw("e2"); //e2same
	//draw bkg errors
	bkgErrors->Draw("e2same"); //e2same
	bkgStatErrors->Draw("e2same");
	//bkgStatErrors->Draw("e2same");
	//dataC->Draw("pe1same"); //esame
	obsRatio->Draw("pe1same");
	legend2->Draw("same");
	gPad->RedrawAxis();
	//Draw line at 1 on ratio plot
    //double xmax = dataC->GetBinCenter(dataC->GetNbinsX()) + dataC->GetBinWidth(dataC->GetNbinsX())/2;
    //double xmin = dataC->GetBinCenter(0) + dataC->GetBinWidth(0)/2;
	double xmax = data->GetBinCenter(data->GetNbinsX()) + data->GetBinWidth(data->GetNbinsX())/2;
    double xmin = data->GetBinCenter(0) + data->GetBinWidth(0)/2;
    TLine *line = new TLine(xmin, 1, xmax, 1);
    line->SetLineStyle(2);
    line->Draw("same");
    c->SaveAs("plots/" + file + ".pdf");
	c->SaveAs("plots/" + file + ".png");
    //temp, remove later!!!
    c->SaveAs("plots/" + file + ".root");
    c->SaveAs("plots/" + file + ".C");
}

void plotHistRatio(TH1D *h1, TH1D *h2, const TString &entry1, const TString &entry2, const TString &file, const bool fit, const double fitmin, const double fitmax, const bool Ylog, const bool abs, const bool data, const bool Xlog, const bool UseTitle, const TString& Title, const bool WriteFit){
    const unsigned nHist = 2;
    TH1D *h_copies[nHist];
    h_copies[0] = (TH1D*) h1->Clone();
    h_copies[1] = (TH1D*) h2->Clone();
    for(unsigned i = 0; i < nHist; ++i){
        histcol(h_copies[i],colors[i]);
    }
    double maximum = 0;
    for(int i = 0; i < nHist; ++i){
        if(!abs) h_copies[i]->Scale(1/h_copies[i]->Integral(), "width");
        if(h_copies[i]->GetBinContent(h_copies[i]->GetMaximumBin()) > maximum) maximum = h_copies[i]->GetBinContent(h_copies[i]->GetMaximumBin());
    }
    h_copies[0]->SetMaximum(maximum*1.1);
    if(!Ylog) h_copies[0]->SetMinimum(0.0001);

    //const double xPad = 0.3;
    TCanvas *c =  new TCanvas(file,"",1920*(1-xPad),1080);				//used to be 2000 and 1000 instead of 1920 1080
	//TCanvas *c =  new TCanvas(file,"",600*(1-xPad),600);
    c->cd();

    TPad* p1, *p2;
    p1 = new TPad(file,"",0,xPad,1,1);						//what value to use for the last parameter?
    p1->Draw();
    p1->cd();				//Why this cd statement?
    //p1->SetTopMargin(0.1);
    p1->SetBottomMargin(0);
    if(Xlog) p1->SetLogx();	
    if(Ylog) p1->SetLogy();	
    //TLegend* leg = new TLegend(0.73,0.55,0.93,0.9,NULL,"brNDC");
	TLegend* leg = new TLegend(0.6,0.55,0.93,0.9,NULL,"brNDC");
    leg->SetFillStyle(0);			//NEW 23 april
    leg->SetTextSize(.08);
    //	leg->SetFillColor(0);
    leg->AddEntry(h_copies[0], entry1);
    leg->AddEntry(h_copies[1], entry2);
    if(!abs){
        h_copies[0]->DrawNormalized("histE");
        h_copies[1]->DrawNormalized("samehistE");
    }
    else{
        h_copies[0]->Draw("histE");
        h_copies[1]->Draw("samehistE");
    }
    leg->Draw("same");
	drawLumi(p1, "Simulation", false);

    if(UseTitle){
        //TLegend* leg2 = new TLegend(0.1, xPad + 0.1, 0.3, xPad + 0.3,NULL,"brNDC");
        //leg2->AddEntry(Title);
        //leg2->Draw("same");
        TLatex latex;
        latex.SetNDC();
        latex.SetTextAngle(0);
        latex.SetTextColor(kBlue); 
        latex.SetTextFont(42);
        //latex.SetTextAlign(31); 
        latex.SetTextSize(0.05);   
        latex.DrawLatex(0.17, 0.1, Title);
    }


    gPad->RedrawAxis();


    //code to try and make the cms logo
    //if(!data) SetCMSText("Simulation Preliminary", false);
    //else  SetCMSText("Preliminary", true);
   // CMS_lumi(p1,"Simulation", false);
    c->cd();

    p2 = new TPad(file + "2","",0,0.0,1,xPad);
    p2->Draw();
    p2->cd();
    p2->SetTopMargin(0);
    p2->SetBottomMargin(0.4);

    TH1D* histratio;
    histratio = (TH1D*) (HistDiv(h_copies[0],h_copies[1], abs))->Clone();
    histratio->SetMarkerColor(1);
    histratio->SetLineColor(1);
    histratio->GetYaxis()->SetRangeUser(0.,1.999);
    histratio->GetYaxis()->SetTitle("Ratio");//entry1 + "/" + entry2);
	//histratio->GetYaxis()->SetTitle("mll01WZ/old WZ");
    //histratio->GetXaxis()->SetTitle("");

    histratio->GetYaxis()->SetTitleOffset(1.25/((1.-xPad)/xPad));
    histratio->GetYaxis()->SetTitleSize((1.-xPad)/xPad*0.06); //originally 0.06
    histratio->GetXaxis()->SetTitleSize((1.-xPad)/xPad*0.06); //originally 0.09
    histratio->GetYaxis()->SetLabelSize((1.-xPad)/xPad*0.05); //originall 0.05 
    histratio->GetXaxis()->SetLabelSize((1.-xPad)/xPad*0.05); //originally 0.075

    if(histratio->GetBinContent(histratio->GetMaximumBin()) < 2) histratio->GetYaxis()->SetRangeUser(0.,1.999);
    else if(histratio->GetBinContent(histratio->GetMaximumBin()) < 5) histratio->SetMaximum((histratio->GetBinContent(histratio->GetMaximumBin()))*1.2);
    if(fit){
        TF1 *f = new TF1("f", "[0]",fitmin,fitmax);
        /*
           TF1 *f = new TF1("fa1","-0.5*[0]*(TMath::Erf((x-[1])/[2])-1)",fitmin,fitmax);
           f->SetParLimits(0,0.001,10);f->SetParLimits(1,40,150);f->SetParLimits(2,0.2,250); 
           f->SetLineColor(kPink);
         */
        histratio->Fit(f,"RS");
        histratio->Fit(f,"RS");
        if(WriteFit){
            std::ofstream fitfile;
            fitfile.open("fitratio.txt");
            fitfile <<  f->GetParameter(0);
            fitfile.close();
        }
    }
    histratio->Draw("pe");

    double xmax = histratio->GetBinCenter(histratio->GetNbinsX()) + histratio->GetBinWidth(histratio->GetNbinsX())/2;
    double xmin = histratio->GetBinCenter(0) + histratio->GetBinWidth(0)/2;
    TLine *line = new TLine(xmin, 1, xmax, 1);
    line->SetLineStyle(2);
    line->Draw("same");
    c->SaveAs("plots/" + file + ".pdf");
	c->SaveAs("plots/" + file + ".png");
}


void plotHist(TH1D* h, const TString &file, const bool data, const bool ylog, const bool abs, const bool xzoom, const double xmin, const double xmax){
    histcol(h, colors[0]);
    if(!ylog) h->SetMaximum(h->GetBinContent(h->GetMaximumBin())*1.3);
    else h->SetMaximum(h->GetBinContent(h->GetMaximumBin())*3);
    h->GetYaxis()->SetTitleOffset(1.4);
    h->GetXaxis()->SetTitleOffset(1.4);
    HistLabelSizes(h,0.05,0.05,0.05,0.05);
    if(xzoom){
        h->GetXaxis()->SetRangeUser(xmin, xmax);
    }
    TCanvas *c  = new TCanvas(file,file,1920*(1-xPad),1080);
    //c->SetTopMargin(0.12);
    c->SetBottomMargin(0.15);
    if(ylog) c->SetLogy();
    h->Draw("histE");
    gPad->RedrawAxis();
	drawLumi(c);
    //if(!data) SetCMSTextSmall("Simulation Preliminary", false);
    //else  SetCMSTextSmall("Preliminary", true);
  //  if(!data) CMS_lumi(c,"Simulation", false);
  //  else CMS_lumi(c,"Preliminary");
    c->cd();
    c->SaveAs("plots/" + file + ".pdf");
	c->SaveAs("plots/png/" + file + ".png");
}

void plotHist(TH2D* h, const TString &file, const bool data, const bool ylog, const bool xlog, const bool zlog){
    //histcol(h, colors[0]);
    h->GetYaxis()->SetTitleOffset(1.4);
    h->GetXaxis()->SetTitleOffset(1.4);
    h->SetMinimum(0);
    HistLabelSizes(h,0.05,0.05,0.05,0.05);
    //TCanvas *c  = new TCanvas(file,file,1920*(1-xPad),1080);
	TCanvas *c  = new TCanvas(file,file,600,600);
    c->SetRightMargin(0.12); //test
    //c->SetTopMargin(0.12);
	c->SetTopMargin(0.05);
    c->SetBottomMargin(0.15);
    if(ylog) c->SetLogy();
    if(xlog) c->SetLogx();
    if(zlog) c->SetLogz();
	//h->SetMarkerSize(1.5);	
    h->Draw("COLZtextE");
	//h->Draw("COLZ");
	if(!data) drawLumi(c, "Simulation", false);
	else drawLumi(c);
    gPad->RedrawAxis();
    //if(!data) SetCMSTextSmall("Simulation Preliminary", false);
    //else  SetCMSTextSmall("Preliminary", true);
  //  if(!data) CMS_lumi(c,"Simulation", false);
  //  else CMS_lumi(c,"Preliminary");
    c->cd();
    c->SaveAs("plots/" + file + ".pdf");
	c->SaveAs("plots/png/" + file + ".png");
}

void plotHist(std::vector<TH1D*> histvec, const std::vector<TString>& entries, const TString &file, const bool ylog, const bool abs, const bool xzoom, const double xmin, const double xmax){
    if(histvec.size() != entries.size()){
        std::cout << "ERROR: different amount of entries and histograms" << std::endl;
        return;
    }
    std::vector<TH1D*> histcopies;
    for(int i = 0; i < histvec.size(); ++i){
        histcopies.push_back( (TH1D*) histvec[i]->Clone() );
    }
    for(int i = 0; i != histcopies.size(); ++i){
        histcol(histcopies[i], colors[i]);
    }

    if(xzoom){
        (*histcopies.begin())->GetXaxis()->SetRangeUser(xmin, xmax);		//Perhaps this will have to be done for all histograms
    }
    double maximum = 0;
    for(int i = 0; i < histcopies.size(); ++i){
        if(!abs) histcopies[i]->Scale(1/histcopies[i]->Integral(), "width");
        if(histcopies[i]->GetBinContent(histcopies[i]->GetMaximumBin()) > maximum){
            maximum = histcopies[i]->GetBinContent(histcopies[i]->GetMaximumBin());
        }
    }

    if(!ylog) (*histcopies.begin())->SetMaximum(maximum*1.2);
	else (*histcopies.begin())->SetMaximum(maximum*3);
    (*histcopies.begin())->GetYaxis()->SetTitleOffset(1.4);
    (*histcopies.begin())->GetXaxis()->SetTitleOffset(1.4);
    HistLabelSizes(*histcopies.begin(),0.05,0.05,0.05,0.05);

    TCanvas *c  = new TCanvas(file,file,1920*(1-xPad),1080);
    c->SetTopMargin(0.12);
    c->SetBottomMargin(0.15);
    if(ylog) c->SetLogy();
    TLegend* leg = new TLegend(0.6,0.55,0.93,0.9,NULL,"brNDC"); 			//0.6 used to be 0.73
    leg->SetFillStyle(0);			//NEW 23 april
    leg->SetTextSize(.05);
    for(int i = 0; i < histcopies.size(); ++i){
        leg->AddEntry(histcopies[i], entries[i]);
    }
    if(!abs){
        (*histcopies.begin())->DrawNormalized("histE");
        for(std::vector<TH1D*>::iterator it = (histcopies.begin() + 1); it != histcopies.end(); ++it){
            (*it)->DrawNormalized("samehistE");
        }
    }
    else{
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		std::cout << (*histcopies.begin())->GetSumOfWeights() << std::endl;
        (*histcopies.begin())->Draw("histE");
        for(std::vector<TH1D*>::iterator it = (histcopies.begin() + 1); it != histcopies.end(); ++it){
            (*it)->Draw("samehistE");
			std::cout << (*it)->GetSumOfWeights() << std::endl;
        }
    }
    leg->Draw("same");

    gPad->RedrawAxis();
    //SetCMSTextSmall("Simulation Preliminary", false);
 	drawLumi(c,"Simulation", false);

    c->cd();
    //Set label sizes
    //HistLabelSizes(h_copies[0]);
    if(ylog) c->SetLogy();
    c->SaveAs("plots/" + file + ".pdf");
	c->SaveAs("plots/" + file + ".png");
}

void plotRoc(TH1D* sigeff, TH1D* bkgeff, const TString& file, const TString* pointNames, const bool legendUp, const TString& xName, const TString yName){
	if(sigeff->GetNbinsX() != bkgeff->GetNbinsX()){
		std::cerr << "ERROR: different number of background- and signal efficiency bins" << std::endl;
		return;
	}
	TCanvas* c = new TCanvas("roccanv","roccanv",500,500);
	c->SetTopMargin(0.08);
	const unsigned nPoints = sigeff->GetNbinsX();
	const Color_t dotColors[20]  = {kCyan +2, kYellow +1,  kBlue, kGreen +2, kMagenta, kOrange +1, kMagenta -1, kCyan,  kYellow, kAzure, kViolet, kTeal, kSpring, kGreen - 2, kRed - 7, kYellow + 3, kMagenta - 8, kBlack, kBlue - 7, kOrange + 2};
	/*
	for(unsigned p = 0; p < nPoints && p < 9; ++p){
		dotColors[p] = gROOT->GetColor(p);
	}
	if(nPoints > 9){
		for(unsigned p = 
	*/
    TGraphAsymmErrors* roc[nPoints];
	
	double sigeffRange[2] = {sigeff->GetBinContent(1), sigeff->GetBinContent(1)};
	double bkgeffRange[2] = {bkgeff	->GetBinContent(1), bkgeff->GetBinContent(1)};	
	for(unsigned p = 1; p < nPoints; ++p){
		if(sigeff->GetBinContent(p + 1) > sigeffRange[1]){
			sigeffRange[1] = sigeff->GetBinContent(p + 1);
		}
		if(sigeff->GetBinContent(p + 1) < sigeffRange[0]){
			sigeffRange[0] = sigeff->GetBinContent(p + 1);
		}
		if(bkgeff->GetBinContent(p + 1) > bkgeffRange[1]){
			bkgeffRange[1] = bkgeff->GetBinContent(p + 1);
		}
		if(bkgeff->GetBinContent(p + 1) < bkgeffRange[0]){
			bkgeffRange[0] = bkgeff->GetBinContent(p + 1);
		}
	}
	sigeffRange[0]*=0.8;
	sigeffRange[1]*=1.1;
	bkgeffRange[0]*=0.8;
	bkgeffRange[1]*=1.1;
	for(int n = nPoints - 1; n >= 0; --n){
		double* x = new double[n + 1];
		double* xerr = new double[n + 1];
		double* y = new double[n + 1];
		double* yerr = new double[n + 1];
		for(int p = 0; p <= n; ++p){
			x[p] = bkgeff->GetBinContent(p + 1);
			xerr[p] = bkgeff->GetBinError(p +1);
			y[p] = sigeff->GetBinContent(p + 1);
			yerr[p] = sigeff->GetBinError(p + 1);
		}
		roc[n] = new TGraphAsymmErrors(nPoints, x, y, xerr, xerr, yerr, yerr);
		roc[n]->SetMarkerColor(dotColors[n]);
		roc[n]->GetXaxis()->SetTitle(xName);
		roc[n]->GetYaxis()->SetTitle(yName);
		//roc[n]->GetXaxis()->SetNdivisions(10, false);
		//roc[n]->GetYaxis()->SetNdivisions(10, false);
		roc[n]->GetXaxis()->SetLimits(bkgeffRange[0], bkgeffRange[1]);
		roc[n]->GetYaxis()->SetRangeUser(sigeffRange[0], sigeffRange[1]);
		//roc->GetXaxis()->SetMoreLogLabels();
		//roc->GetYaxis()->SetMoreLogLabels();
	}
	double legXdown, legXup, legYdown, legYup;
	if(!legendUp){
		legXdown = 0.4;
		legXup =  0.7;
		legYdown = 0.2;
   		legYup = 0.5;
	} else{
		legXdown = 0.7;
		legXup = 1;
		legYdown = 0.6;
		legYup = 0.9;
	}
	TLegend* legend = new TLegend(legXdown,legYdown,legXup,legYup,NULL,"brNDC");
	//TLegend* legend = new TLegend(0.4,0.2,0.7,0.5,NULL,"brNDC");
    //legend->SetFillStyle(0);
    for(int p = nPoints - 1; p >= 0; --p){
    	legend->AddEntry(roc[p], pointNames[p], "p");
    }
	roc[nPoints - 1]->Draw("AP");
	for(int p = nPoints - 2; p >= 0; --p){
		roc[p]->Draw("Psame");
	}
	legend->Draw("same");
	//if(sigeffRange[0] < 0.1) c->SetLogy();
	//if(bkgeffRange[0] < 0.1) c->SetLogx();
	drawLumi(c, "Simulation", false);
    c->SaveAs("plots/" + file + ".pdf");
}





//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void makePaperPlotHNL(TH1D* dataYield, TH1D** bkgYields, TH1D** bkgErrors, const TString* bkgNames, const unsigned nBkg, const TString& title, const TString& fileName){
	const double xPadLocal = 0.35;
	//Order histograms
	const TString orderedNames[6] = {"nonprompt", "WZ", "X#gamma^{(*)}", "ZZ/H", "triboson", "t#bar{t}/t + X"};//"TT/T + X"};
	TH1D* orderedBkg[nBkg];
	TH1D* orderedBkgErrors[nBkg];
	std::cout << "nBkg = " << nBkg << std::endl;
	for(unsigned bkg = 0; bkg < nBkg; ++bkg){
		for(unsigned o = 0; o < nBkg; ++o){
			if(bkgNames[bkg] == orderedNames[o]){
				orderedBkg[o] = (TH1D*) bkgYields[bkg]->Clone();
				orderedBkgErrors[o] = (TH1D*) bkgErrors[bkg]->Clone();
			}
		}
	}
	//Histogram containing total bkg with stat unc.
	TH1D* totalBkg = (TH1D*) orderedBkg[0]->Clone();
	for(unsigned bkg = 1; bkg < nBkg; ++bkg){
		totalBkg->Add(orderedBkg[bkg]);
	}
	//Make histogram containing bkg yields with the total uncertainty
	TH1D* bkgYieldsE[nBkg];
	for(unsigned bkg = 0; bkg < nBkg; ++bkg){
		bkgYieldsE[bkg] = (TH1D*) orderedBkg[bkg]->Clone();
	}
	for(unsigned bkg = 0; bkg < nBkg; ++bkg){
		for(unsigned bin = 1; bin < bkgYieldsE[bkg]->GetNbinsX() + 1; ++bin){
			bkgYieldsE[bkg]->SetBinError(bin, sqrt(bkgYieldsE[bkg]->GetBinError(bin)*bkgYieldsE[bkg]->GetBinError(bin) + orderedBkgErrors[bkg]->GetBinContent(bin)*orderedBkgErrors[bkg]->GetBinContent(bin) ) );
		}
	}
	//Histogram containing total bkg with total unc.
	TH1D* totalBkgE = (TH1D*) bkgYieldsE[0]->Clone();
	for(unsigned bkg = 1; bkg < nBkg; ++bkg){
		totalBkgE->Add(bkgYieldsE[bkg]);
	}
	//Make bkg stack 
	THStack* bkgStack = new THStack("bkgStack", "bkgStack");
	for(int bkg = nBkg -1; bkg > -1; --bkg){
		StackCol(bkgYieldsE[bkg], bkgColor_HNL(orderedNames[bkg]) );
		bkgStack->Add(bkgYieldsE[bkg], "f");
	}
	//Calculate maximum and minimum to set plot ranges
	double dataMax = dataYield->GetBinContent(dataYield->GetMaximumBin()) + dataYield->GetBinError(dataYield->GetMaximumBin());
	double bkgMax = totalBkgE->GetBinContent(totalBkgE->GetMaximumBin()) + totalBkgE->GetBinError(totalBkgE->GetMaximumBin());
	double maximum = std::max(dataMax, bkgMax);
   	double minimum = maximum;
	for(unsigned b = 1; b < totalBkgE->GetNbinsX() + 1; ++b){
		if(totalBkg->GetBinContent(b) != 0 && totalBkg->GetBinContent(b) < minimum){
			minimum = totalBkgE->GetBinContent(b);
		}
	}
	minimum /= 5;
	double SF = log10(std::max(10., maximum/minimum ));
	dataYield->SetMinimum(minimum);
	//dataYield->SetMaximum(maximum*3*SF);
	dataYield->SetMaximum(maximum*9*SF);
	//Define plotting canvas
	TCanvas *c =  new TCanvas(fileName, "",1050,450);
    c->cd();

	//Make legend 
	//TLegend* legend = new TLegend(0.875,0,1,1,NULL,"brNDC");
	//Make upper pad to draw main plot and lower pad for ratios
	TPad* p1, *p2, *p3, *p4;
    p1 = new TPad(fileName,"",0,xPadLocal,0.85,1);
	p1->SetLogy();
    p1->Draw();
    p1->cd();
    //p1->SetBottomMargin(0);
	p1->SetBottomMargin(0.01);
	p1->SetRightMargin(0.01);
	p1->SetLeftMargin(1);
	p1->SetTopMargin(0.07);
	//Draw histograms on upper part of plot
	//Make the total background uncertainty visible as a grey band
    totalBkgE->SetFillStyle(3244); //3005  3244
    totalBkgE->SetFillColor(kGray+2);
    totalBkgE->SetMarkerStyle(0); //1

	//Set Correct offsets
	dataYield->GetYaxis()->SetTitleOffset(0.7);

	
	//Set bin labeles
	const TString labels[33] = {"0-10","10-20", "20-30", "> 30", "0-10","10-20", "20-30", "> 30", "0-100", "> 100", "0-100", "100-150", "150-250", ">250", "0-100", "> 100", "> 0", "0-100", "100-200", "200-300", "0-100", "100-200", "200-300", "300-400", "> 400", "0-100", "100-200", "200-300", "> 300", "0-100", "100-200", "200-300", ">300"};
	/*
	for(unsigned bin = 1; bin < dataYield->GetNbinsX() + 1; ++bin){
		dataYield->GetXaxis()->SetBinLabel(bin, labels[bin -1]);
	}
	dataYield->GetXaxis()->LabelsOption("v");
	dataYield->GetXaxis()->SetTitle("");
	dataYield->GetXaxis()->SetLabelSize(0);
	dataYield->Draw("pe1");
	bkgStack->Draw("hist same");
	totalBkgE->Draw("e2 same");
	dataYield->Draw("pe1 same");
	*/
	totalBkgE->SetMinimum(minimum);
	totalBkgE->SetMaximum(maximum*9*SF);
	totalBkgE->GetYaxis()->SetTitleOffset(0.7);
	TGraphAsymmErrors* obs = new TGraphAsymmErrors(dataYield);
	for(unsigned b = 1; b < dataYield->GetNbinsX() + 1; ++b){
		obs->SetPointError(b - 1, 0, 0, dataYield->GetBinErrorLow(b), (dataYield->GetBinContent(b) == 0 ) ? 0 : dataYield->GetBinErrorUp(b) );
		//obs->SetPointEXlow(b, 0);
		//obs->SetPointEXhigh(b, 0);
	}
	totalBkgE->Draw("e2");
	bkgStack->Draw("hist same");
	totalBkgE->Draw("e2 same");
	obs->Draw("pe1same");
	/*
	bkgStack->Draw("hist same");
	totalBkgE->Draw("e2 same");
	obs->Draw("pe1 same");
	*/
	//bin numbers where to draw solid line (inbetween categories)
	/*
	const unsigned solidBins[3] = {5, 9, 25};
	//bin numbers where to draw broken line (intbetween high mass minMos and M3l bins
	const unsigned brokenBins[6] = {12, 17, 21, 27, 31, 33};
	const TString brokenText[6] = {"M_{3l} > 100 GeV", "100 GeV < minM_{OS} < 200 GeV", "minM_{OS} > 200 GeV",  "M_{3l} > 100 GeV", "100 GeV < minM_{OS} < 200 GeV", "minM_{OS} > 200 GeV"};
	*/
	const unsigned solidBins[3] = {5, 9, 18};
	/*
	const unsigned brokenBins[6] = {11, 15, 17, 21, 26, 30};
	const TString brokenText[6] = {"M_{3l} > 100 GeV", "100 GeV < minM_{OS} < 200 GeV", "minM_{OS} > 200 GeV",  "M_{3l} > 100 GeV", "100 GeV < minM_{OS} < 200 GeV", "minM_{OS} > 200 GeV"};
	*/
	/*
	const unsigned nBroken = 8;
	const unsigned brokenBins[nBroken] = {11, 15, 17, 18, 21, 26, 30, 33};
	const TString brokenText[nBroken] = {"M_{3l} < 100 GeV", "minM_{OS} < 100 GeV", "100 GeV < minM_{OS} < 200 GeV", "minM_{OS} > 200 GeV",  "M_{3l} < 100 GeV", "minM_{OS} < 100 GeV", "100 GeV < minM_{OS} < 200 GeV", "minM_{OS} > 200 GeV"};
	*/
	const unsigned nBroken = 6;
	const unsigned brokenBins[nBroken] = {15, 17, 18, 26, 30, 33};
	const TString mVar = "M_{2lOS}^{min}";//minM_{OS};
	const TString brokenText[nBroken] = {mVar + " < 100 GeV", "100 GeV < " + mVar + " < 200 GeV", mVar + " > 200 GeV", mVar + " < 100 GeV", "100 GeV < " + mVar + " < 200 GeV", mVar + " > 200 GeV"};
	const unsigned nM3l = 2;
	const unsigned m3lBins[nBroken] = {11, 21};

	//Draw lines at correct bins
	for(unsigned sb = 0; sb < 3; ++sb){
		TLine* tempLine = new TLine(dataYield->GetBinLowEdge(solidBins[sb]), 0, dataYield->GetBinLowEdge(solidBins[sb]), maximum*SF*9 );//maximum*SF*3); //1000 //9?
		tempLine->SetLineStyle(1);
		tempLine->SetLineWidth(2);
		tempLine->SetLineColor(kBlue + 3);
		tempLine->Draw(); 
		//delete tempLine;
	}
	//find minimum of all bins with text
	double maxBin = 0;
	for(unsigned bb = 0; bb < nBroken; ++bb){
		if(std::max(dataYield->GetBinContent(brokenBins[bb] - 1) + dataYield->GetBinError(brokenBins[bb] - 1), totalBkg->GetBinContent(brokenBins[bb] - 1) + totalBkg->GetBinError(brokenBins[bb] - 1) ) > maxBin){
			maxBin = std::max(dataYield->GetBinContent(brokenBins[bb] - 1) + dataYield->GetBinError(brokenBins[bb] - 1), totalBkg->GetBinContent(brokenBins[bb] - 1) + totalBkg->GetBinError(brokenBins[bb] - 1) );
		}
	}
	for(unsigned bb = 0; bb < nBroken; ++bb){
		TLatex *xlabel = new TLatex(dataYield->GetBinLowEdge(brokenBins[bb]) - 0.4, maxBin*2, brokenText[bb]); //100  //maximum/5
		//xlabel->SetNDC();
		//xlabel->SetTextFont(1);
		xlabel->SetTextColor(kBlue + 2);
		xlabel->SetTextSize(0.035);
		//xlabel->SetTextAlign(22);
		xlabel->SetTextAlign(1);
		xlabel->SetTextAngle(90);
		//xlabel->PaintTextNDC(0.5, 0.5, "text");
		xlabel->Draw("same");
		//if(bb == 7) continue;
		if(bb == 2 || bb == 5) continue;
		TLine* tempLine = new TLine(dataYield->GetBinLowEdge(brokenBins[bb]), 0, dataYield->GetBinLowEdge(brokenBins[bb]), maximum*2); //500 //maximum*3*SF);
		tempLine->SetLineStyle(2);
		tempLine->SetLineWidth(2);
		tempLine->SetLineColor(kBlue + 3);
		tempLine->Draw(); 
		//delete tempLine;
	}

	for(unsigned bb = 0; bb < nM3l; ++bb){
		TLatex *xlabel = new TLatex(dataYield->GetBinLowEdge(m3lBins[bb]) - 0.4, maximum*1.3, "M_{3l} < 100 GeV"); //100  //maximum/5
		xlabel->SetTextColor(kBlue + 2);
		xlabel->SetTextSize(0.035);
		//xlabel->SetTextAlign(22);
		xlabel->SetTextAlign(1);
		xlabel->SetTextAngle(90);
		//xlabel->PaintTextNDC(0.5, 0.5, "text");
		xlabel->Draw("same");
		TLatex *xlabel2 = new TLatex(dataYield->GetBinLowEdge(m3lBins[bb] + 1) - 0.4, maximum*1.3, "M_{3l} > 100 GeV"); //100  //maximum/5
		xlabel2->SetTextColor(kBlue + 2);
		xlabel2->SetTextSize(0.035);
		//xlabel->SetTextAlign(22);
		xlabel2->SetTextAlign(1);
		xlabel2->SetTextAngle(90);
		//xlabel->PaintTextNDC(0.5, 0.5, "text");
		xlabel2->Draw("same");
		//if(bb == 7) continue;
		TLine* tempLine = new TLine(dataYield->GetBinLowEdge(m3lBins[bb]), 0, dataYield->GetBinLowEdge(m3lBins[bb]), maximum*7); //500 //maximum*3*SF);
		tempLine->SetLineStyle(7);
		tempLine->SetLineWidth(2);
		tempLine->SetLineColor(kBlue + 3);
		tempLine->Draw(); 
		//delete tempLine;
	}
	
	//Label all categories
	/*
	const TString categoryText[4] = {"#splitline{P_{T}^{lead} < 30 GeV}{#scale[1.3]{no OSSF}}", "#splitline{P_{T}^{lead} < 55 GeV}{#scale[1.3]{no OSSF}}", "#splitline{P_{T}^{lead} > 55 GeV}{#scale[1.3]{OSSF}}", "#splitline{P_{T}^{lead} > 55 GeV}{#scale[1.3]{no OSSF}}"};
	const double textPositions[4] = {1.3, 4.9, 12.1, 26.8};
	*/
	TString noOSSFstring = "no OSSF";
	TString OSSFstring = "OSSF";
	if(title == "3 leptons, #geq 2 muons"){
		noOSSFstring = "#bf{e^{#pm}#mu^{#mp}#mu^{#mp}}";
		OSSFstring = "#splitline{#bf{e^{#pm}#mu^{#pm}#mu^{#mp}}}{#bf{#mu^{#pm}#mu^{#pm}#mu^{#mp}}}";
	} else if(title == "3 leptons, #geq 2 electrons"){
		noOSSFstring = "#bf{e^{#pm}e^{#pm}#mu^{#mp}}";
		OSSFstring = "#splitline{#bf{e^{#pm}e^{#mp}#mu^{#mp}}}{#bf{e^{#pm}e^{#pm}e^{#mp}}}";
	}
	const TString categoryText[4] = {"#splitline{P_{T}^{lead} < 30 GeV}{#scale[1.5]{" + noOSSFstring + "}}", "#splitline{P_{T}^{lead} < 55 GeV}{#scale[1.5]{" + noOSSFstring + "}}", "#splitline{P_{T}^{lead} > 55 GeV}{#scale[1.5]{" + noOSSFstring + "}}", "#splitline{P_{T}^{lead} > 55 GeV}{#scale[1.5]{" + OSSFstring + "}}"};
	//const double textPositions[4] = {1.3, 4.9, 11.1, 21.5};
	const double textPositions[4] = {1.35, 4.9, 12, 21.8}; //1.3
	for(unsigned c = 0; c < 4; ++c){
		//TLatex *xlabel = new TLatex(textPositions[c], maximum*3, categoryText[c]); //1000
		TLatex *xlabel = new TLatex(textPositions[c], maximum*10, categoryText[c]); //1000
		//xlabel->SetNDC();
		//xlabel->SetTextFont(1);
		xlabel->SetTextAlign(1);
		xlabel->SetTextColor(kBlack);
		xlabel->SetTextSize(0.04); //0.035
		//xlabel->SetTextAlign(22);
		xlabel->SetTextAngle(0);
		//xlabel->PaintTextNDC(0.5, 0.5, "text");
		xlabel->Draw("same");
	}

	//redraw axis 
	gPad->RedrawAxis();
	
	/*
	//Make TLatex to draw extra info in plot
	TLatex latex(l,1+lumiTextOffset*t,"CMS");
	latex.SetNDC();
	latex.SetTextAngle(0);
	latex.SetTextColor(kBlack); 
	*/

	//Set individual bin labels
	c->cd();
	p2 = new TPad(fileName,"",0, 0, 0.85, xPadLocal);
	p2->Draw();
	p2->cd();
	//p2->SetTopMargin(0);
	p2->SetTopMargin(0.05);
	p2->SetLeftMargin(0.1);
	p2->SetRightMargin(0.01);
	p2->SetBottomMargin(0.5);
	//Make histogram with relative stat bkg. unc.
	const unsigned nBins = dataYield->GetNbinsX();
	TH1D* bkgStatErrors = new TH1D("bkgStatErrors" + fileName, "bkgStatErrors" + fileName, nBins, dataYield->GetBinLowEdge(1), dataYield->GetBinLowEdge(dataYield->GetNbinsX()) + dataYield->GetBinWidth(dataYield->GetNbinsX()));
	for(unsigned b = 1; b < nBins + 1; ++b){
		bkgStatErrors->SetBinContent(b, 1.);
		if(totalBkg->GetBinContent(b) != 0){
			bkgStatErrors->SetBinError(b, (totalBkg->GetBinError(b) )/ (totalBkg->GetBinContent(b) ) );
		} else{
			bkgStatErrors->SetBinError(b, 0.);
		}			
	}
	//Make histogtam with relative tot bkg unc.
	TH1D* bkgTotErrors = new TH1D("bkgTotErrors" + fileName, "bkgTotErrors" + fileName, nBins, dataYield->GetBinLowEdge(1), dataYield->GetBinLowEdge(dataYield->GetNbinsX()) + dataYield->GetBinWidth(dataYield->GetNbinsX()));
	for(unsigned b = 1; b < nBins + 1; ++b){
		bkgTotErrors->SetBinContent(b, 1.);
		if(totalBkgE->GetBinContent(b) != 0){
			bkgTotErrors->SetBinError(b, totalBkgE->GetBinError(b)/totalBkgE->GetBinContent(b));
		} else{
			bkgTotErrors->SetBinError(b, 0.);
		}
	}				
	bkgStatErrors->SetFillStyle(1001);
	bkgTotErrors->SetFillStyle(1001);
	bkgStatErrors->SetFillColor(kCyan  - 4);
	bkgTotErrors->SetFillColor(kOrange	- 4);
    bkgStatErrors->SetMarkerStyle(1);
	bkgTotErrors->SetMarkerStyle(1);
	
	TH1D* dataYieldC = (TH1D*) dataYield->Clone();
	/*
	//make sure the data points only have their stat unc plotted
	TH1D* totalBkgNoErrs = (TH1D*) totalBkgE->Clone();
	
	for(unsigned b = 1; b < totalBkgNoErrs->GetNbinsX(); ++b){
		totalBkgNoErrs->SetBinError(b, 0);
	}
	*/
	/*
	for(unsigned b = 1; b < dataYieldC->GetNbinsX() + 1; ++b){
		if(dataYieldDiv->GetBinContent(b) != 0 && totalBkgE->GetBinContent(b) != 0) dataYieldDiv->SetBinContent(b, dataYieldDiv->GetBinContent(b)/totalBkgE->GetBinContent(b));
	}
	*/

	TGraphAsymmErrors* obsRatio = new TGraphAsymmErrors(dataYieldC);
	for(unsigned b = 1; b < dataYieldC->GetNbinsX() + 1; ++b){
		obsRatio->GetY()[b - 1] *= 1./totalBkgE->GetBinContent(b);
		//obsRatio->SetPointEYlow(b-1 , dataYieldC->GetBinErrorLow(b)/totalBkgE->GetBinContent(b));
		//obsRatio->SetPointEYhigh(b -1, dataYieldC->GetBinErrorUp(b)/totalBkgE->GetBinContent(b));	
		obsRatio->SetPointError(b - 1, 0, 0, dataYieldC->GetBinErrorLow(b)/totalBkgE->GetBinContent(b), dataYieldC->GetBinErrorUp(b)/totalBkgE->GetBinContent(b));
	}
    //dataYieldC->Divide(totalBkgE);
	//dataYieldC->Divide(totalBkgNoErrs);
	//Legend for uncertainties
	//legend2->AddEntry(bkgStatErrors, "stat. bkg. unc.", "f");
	/*
    dataYieldC->SetMarkerColor(1);
    dataYieldC->SetLineColor(1);
    dataYieldC->GetYaxis()->SetRangeUser(0.,1.999);
    dataYieldC->GetYaxis()->SetTitle("obs./pred.");
    dataYieldC->GetYaxis()->SetTitleOffset(0.7/((1.-xPadLocal)/xPadLocal));
    dataYieldC->GetYaxis()->SetTitleSize((1.-xPadLocal)/xPadLocal*0.06);
    dataYieldC->GetXaxis()->SetTitleSize((1.-xPadLocal)/xPadLocal*0.06); 
    dataYieldC->GetYaxis()->SetLabelSize((1.-xPadLocal)/xPadLocal*0.05);
    dataYieldC->GetXaxis()->SetLabelSize((1.-xPadLocal)/xPadLocal*0.05);

	dataYieldC->GetXaxis()->SetLabelSize(0.1);
  	dataYieldC->Draw("pe");
	*/
	/*
	obsRatio->SetMarkerColor(1);
    obsRatio->SetLineColor(1);
    obsRatio->GetYaxis()->SetRangeUser(0.,1.999);
    obsRatio->GetYaxis()->SetTitle("obs./pred.");
    obsRatio->GetYaxis()->SetTitleOffset(0.7/((1.-xPadLocal)/xPadLocal));
    obsRatio->GetYaxis()->SetTitleSize((1.-xPadLocal)/xPadLocal*0.06);
    obsRatio->GetXaxis()->SetTitleSize((1.-xPadLocal)/xPadLocal*0.06); 
    obsRatio->GetYaxis()->SetLabelSize((1.-xPadLocal)/xPadLocal*0.05);
    obsRatio->GetXaxis()->SetLabelSize((1.-xPadLocal)/xPadLocal*0.05);
	obsRatio->GetXaxis()->SetLabelSize(0.1);
  	obsRatio->Draw("pe");
	*/
	//draw bkg errors
	for(unsigned bin = 1; bin < bkgTotErrors->GetNbinsX() + 1; ++bin){
		bkgTotErrors->GetXaxis()->SetBinLabel(bin, labels[bin -1]);
	}
	bkgTotErrors->SetLabelOffset(0.03);
	bkgTotErrors->GetXaxis()->LabelsOption("v");
	bkgTotErrors->GetXaxis()->SetTitle("");
	bkgTotErrors->GetXaxis()->SetLabelSize(0);
	bkgTotErrors->SetMarkerColor(1);
    bkgTotErrors->SetLineColor(1);
    bkgTotErrors->GetYaxis()->SetRangeUser(0.,3.999);//1.999
    bkgTotErrors->GetYaxis()->SetTitle("obs./pred.");
    bkgTotErrors->GetYaxis()->SetTitleOffset(0.7/((1.-xPadLocal)/xPadLocal));
    bkgTotErrors->GetYaxis()->SetTitleSize((1.-xPadLocal)/xPadLocal*0.06);
    bkgTotErrors->GetXaxis()->SetTitleSize((1.-xPadLocal)/xPadLocal*0.06); 
    bkgTotErrors->GetYaxis()->SetLabelSize((1.-xPadLocal)/xPadLocal*0.05);
    bkgTotErrors->GetXaxis()->SetLabelSize((1.-xPadLocal)/xPadLocal*0.05);

	dataYieldC->GetXaxis()->SetLabelSize(0.1);

	bkgTotErrors->Draw("e2same"); //e2same
	bkgStatErrors->Draw("e2same");
	obsRatio->Draw("pe01same");
	
	//Hack to avoid having dots drawn for 0 content bins
	for(unsigned b = 1; b < dataYieldC->GetNbinsX() + 1; ++b){
		if(dataYieldC->GetBinContent(b) == 0) dataYieldC->SetBinContent(b, 100.);
	}

	//dataYieldC->Draw("pe0same"); //esame   pe1same e01same
	//dataYieldC->Draw("p0e1same"); //esame   pe1same e01same
	//Draw lines at correct bins
	for(unsigned sb = 0; sb < 3; ++sb){
		TLine* tempLine = new TLine(dataYield->GetBinLowEdge(solidBins[sb]), 0, dataYield->GetBinLowEdge(solidBins[sb]), 4);
		tempLine->SetLineStyle(1);
		tempLine->SetLineWidth(2);
		tempLine->SetLineColor(kBlue + 3);
		tempLine->Draw(); 
		//delete tempLine;
	}
	/*
	for(unsigned bb = 0; bb < 6; ++bb){
		TLine* tempLine = new TLine(dataYield->GetBinLowEdge(brokenBins[bb]), 0, dataYield->GetBinLowEdge(brokenBins[bb]), 2);
		tempLine->SetLineStyle(2);
		tempLine->SetLineWidth(2);
		tempLine->SetLineColor(kBlue + 3);
		tempLine->Draw(); 
		//delete tempLine;
	}
	*/

	for(unsigned bb = 0; bb < nBroken; ++bb){
		if(bb == 2 || bb == 5) continue;
		TLine* tempLine = new TLine(dataYield->GetBinLowEdge(brokenBins[bb]), 0, dataYield->GetBinLowEdge(brokenBins[bb]), 4); //500 //maximum*3*SF);
		tempLine->SetLineStyle(2);
		tempLine->SetLineWidth(2);
		tempLine->SetLineColor(kBlue + 3);
		tempLine->Draw(); 
		//delete tempLine;
	}


	for(unsigned bb = 0; bb < nM3l; ++bb){
		TLine* tempLine = new TLine(dataYield->GetBinLowEdge(m3lBins[bb]), 0, dataYield->GetBinLowEdge(m3lBins[bb]), 4); //500 //maximum*3*SF);
		tempLine->SetLineStyle(7);
		tempLine->SetLineWidth(2);
		tempLine->SetLineColor(kBlue + 3);
		tempLine->Draw(); 
		//delete tempLine;
	}

	//Add extra text to X axis defining the variables 
	TLatex *xlabel = new TLatex(0.11, 0.15, mVar + " (GeV)");
	xlabel->SetNDC();
	xlabel->SetTextAlign(1);
	xlabel->SetTextColor(kBlack);
	xlabel->SetTextSize(0.08);
	xlabel->SetTextAngle(0);
	//xlabel->PaintTextNDC(0.5, 0.5, "text");
	xlabel->Draw("same");

	TLatex *xlabel2 = new TLatex(0.22, 0.15, mVar + " (GeV)");
	xlabel2->SetNDC();
	xlabel2->SetTextAlign(1);
	xlabel2->SetTextColor(kBlack);
	xlabel2->SetTextSize(0.08);
	xlabel2->SetTextAngle(0);
	//xlabel->PaintTextNDC(0.5, 0.5, "text");
	xlabel2->Draw("same");

	//TLatex *xlabel3 = new TLatex(0.48, 0.15, "M_{T} (GeV)");
	TLatex *xlabel3 = new TLatex(0.4, 0.15, "M_{T} (GeV)");
	xlabel3->SetNDC();
	xlabel3->SetTextAlign(1);
	xlabel3->SetTextColor(kBlack);
	xlabel3->SetTextSize(0.08);
	xlabel3->SetTextAngle(0);
	//xlabel->PaintTextNDC(0.5, 0.5, "text");
	xlabel3->Draw("same");

	//TLatex *xlabel4 = new TLatex(0.84, 0.15, "M_{T} (GeV)");
	TLatex *xlabel4 = new TLatex(0.75, 0.15, "M_{T} (GeV)");
	xlabel4->SetNDC();
	xlabel4->SetTextAlign(1);
	xlabel4->SetTextColor(kBlack);
	xlabel4->SetTextSize(0.08);
	xlabel4->SetTextAngle(0);
	//xlabel->PaintTextNDC(0.5, 0.5, "text");
	xlabel4->Draw("same");
	
	gPad->RedrawAxis();
	//Draw line at 1 on ratio plot
    double xmax = dataYieldC->GetBinCenter(dataYieldC->GetNbinsX()) + dataYieldC->GetBinWidth(dataYieldC->GetNbinsX())/2;
    double xmin = dataYieldC->GetBinCenter(0) + dataYieldC->GetBinWidth(0)/2;
    TLine *line = new TLine(xmin, 1, xmax, 1);
    line->SetLineStyle(2);
    line->Draw("same");


	c->cd();
	p3 = new TPad(fileName + "3","3",0.85,0,1,1);
	p3->Draw();	
	p3->cd();
	p3->SetRightMargin(0);
	TLegend* legend = new TLegend(0,0.4*xPadLocal,1,1,NULL,"brNDC"); //xPad
	legend->SetFillStyle(0);
	legend->SetTextFont(42); //62
	legend->SetTextSize(0.085);
	//legend->SetTextSize(0.06);
	legend->AddEntry(dataYield, "observed");
	legend->AddEntry(totalBkgE, "total pred. unc.", "f");
	for(unsigned bkg = 0; bkg < nBkg; ++bkg){
    	legend->AddEntry(bkgYieldsE[bkg], orderedNames[bkg], "f");
    }
	legend->AddEntry(bkgStatErrors, "stat. pred. unc.", "f");
	legend->AddEntry(bkgTotErrors, "total pred. unc.", "f");
	legend->AddEntry(dataYieldC, "#splitline{obs./pred. with}{stat. obs. unc.}", "pe1");
	legend->Draw("");
	
	c->cd();
	drawLumi(p1);
	p1->Draw();
	p1->cd();
	//Add TLatex to top of plot specifying the probed coupling
	TLatex *TitleText = new TLatex(0.5, 0.935, title);
	TitleText->SetNDC();
	TitleText->SetTextAlign(1);
	TitleText->SetTextColor(kBlack);
	TitleText->SetTextSize(0.06);
	TitleText->SetTextAngle(0);
	TitleText->Draw("same"); 
	gPad->RedrawAxis();
	/*
	gPad->RedrawAxis();
	c->cd();
	p4 = new TPad(fileName,"", 0.85, 0, 1, xPad);
	p4->SetRightMargin(0);
	p4->SetTopMargin(0);
	p4->Draw();
	p4->cd();
	
	TLegend* legend2 = new TLegend(0,0,1,1,NULL,"brNDC");
	legend2->SetTextSize(0.09);
	//Avoid legend box
	legend2->SetFillStyle(0);
	//Add data to legend
	legend2->AddEntry(bkgStatErrors, "stat. pred. unc.", "f");
	legend2->AddEntry(bkgTotErrors, "total pred. unc.", "f");
	legend2->AddEntry(dataYieldC, "obs./pred. with total unc.", "pe1");
	//Make canvas that is divided in 2 pads 
	legend2->Draw("");
	*/
	//Save plot
	c->SaveAs("plots/" + fileName + ".pdf");
	c->SaveAs("plots/" + fileName + ".eps");
	c->SaveAs("plots/" + fileName + ".png");
	c->SaveAs("plots/" + fileName + ".root");
	c->SaveAs("plots/" + fileName + ".C");
}


