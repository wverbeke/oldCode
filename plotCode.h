#ifndef PlotScript
#define PlotScript
 
//import ROOT classes
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TH2D.h"
 
//const double xPad = 0.35;
const double xPad = 0.25;
//Color order for plots and the background stack
//const Color_t colors[8] = {kCyan +2, kYellow +1,  kBlue, kGreen +2, kMagenta, kOrange +1, kMagenta -1, kCyan};
const Color_t colors[8] = {kAzure + 1, kGreen - 7, kMagenta -7, kRed - 7, kBlue -3, kOrange + 6, kCyan + 1, kMagenta +3};//kGreen - 7
//Color order for overlaid signals
const Color_t sigCols[7] = {kRed, kBlue, kGreen, kYellow, kCyan, kMagenta, kGreen + 3};
//const Color_t sigcolors[5] = 
//Set histogram colors and lines
void histcol(TH1D *, const Color_t);
//Return histogram divided by other histogram (both are normalized)
TH1D *HistDiv(TH1D *, TH1D *, const bool abs = false);
//Set Histogram labelsizes
void HistLabelSizes(TH1D *h, const double xlabel = 0.045, const double xtitle = 0.05, const double ylabel = 0.045, const double ytitle = 0.045);
void HistLabelSizes(TH2D *h, const double xlabel = 0.045, const double xtitle = 0.05, const double ylabel = 0.045, const double ytitle = 0.045);
//Set Stack colors
void StackCol(TH1D *h, const Color_t);
//Order histograms in terms of yields (biggest first)
void yieldOrder(TH1D**, unsigned*, const unsigned);
void plotDataVSMC(TH1D* data, TH1D** bkg, const TString* names, const unsigned nHist, const TString& file, const bool ylog = false, const unsigned widthopt = 0, const TString& analysis = "", TH1D** bkgSyst = nullptr, const bool plotsig = false, TH1D** signal = nullptr, const TString* signames =nullptr, const unsigned nSig = 0, const bool signorm = false);
void plotHistRatio(TH1D *h1, TH1D *h2, const TString &entry1, const TString &entry2, const TString &file, const bool fit = false, const double fitmin = 0, const double fitmax = 0, const bool Ylog = false, const bool abs = false, const bool data = false, const bool Xlog = false, const bool UseTitle = false, const TString& Title = "", const bool WriteFit = false);
 
void plotHist(TH1D* h, const TString &file, const bool data = false, const bool ylog = false, const bool abs = false, const bool xzoom = false, const double xmin = 0, const double xmax = 0);
void plotHist(std::vector<TH1D*> histvec,const std::vector<TString>& entries, const TString &file, const bool ylog = false, const bool abs = false, const bool xzoom = false, const double xmin = 0, const double xmax = 0);
void plotHist(TH2D* h, const TString &file, const bool data = false, const bool ylog = false, const bool xlog = false, const bool zlog = false);

//Draw a roc curve given two efficiency histograms
void plotRoc(TH1D*, TH1D*, const TString& file, const TString* pointNames, const bool legendUp = false, const TString& xName = "signal efficiency", const TString yName = "background efficiency");
//return color corresponding to bkg
Color_t bkgColor(const TString&); 
//code to make specific plots for the HNL paper
void makePaperPlotHNL(TH1D* dataYield, TH1D** bkgYields, TH1D** bkgErrors, const TString* bkgNames, const unsigned nBkg, const TString& title, const TString& fileName);
#endif
