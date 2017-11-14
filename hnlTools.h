#ifndef HNL_TOOLS
#define HNL_TOOLS

#include <string>
#include <vector>
#include <map>
#include <sstream>

#include "TString.h"
#include "TH1D.h"

namespace hnl{
	unsigned cat(const unsigned* ind, const int* flavors, const double* charges, const unsigned lCount, const double ptlead);//, const double m3l, const double met);
	unsigned sr(const double mt, const double minMos, const double m3l, const unsigned cat);
	unsigned sr_lowM_3lOSSF_lowPt(const double mt, const double minMos);
	unsigned sr_lowM_3lOSSF_highPt(const double mt, const double minMos);
	unsigned sr_lowM_3lnoOSSF_lowPt(const double minMos);
	unsigned sr_lowM_3lnoOSSF_highPt(const double minMos);
	unsigned sr_highM_3lOSSF(const double mt, const double minMos, const double m3l);
	unsigned sr_highM_3lnoOSSF(const double mt, const double minMos, const double m3l);
	unsigned sr_3lOSSF_mtminMos(const double mt, const double minMos, const double m3l);
	unsigned sr_3lnoOSSF_mtminMos(const double mt, const double minMos, const double m3l);
	unsigned flavorComposition(const unsigned* ind, const int* flavors, const double* charges, const unsigned lCount);
	//unsigned controlRegion(const unsigned* ind, const int* flavors, const double* charges, const unsigned lCount, const double m3l, const double mll, const double met, const unsigned nBJets);
	unsigned controlRegion(const unsigned* ind, const int* flavors, const double* charges, const unsigned lCount, const double m3l, const double mll);

	void printDataCard(const double obsYield, const double sigYield, const TString& sigName, const double* bkgYield, const unsigned nBkg, const TString* bkgNames, const std::vector<std::vector<double> >& systUnc, const unsigned nSyst, const TString* systNames, const TString* systDist, const TString& cardName, const bool shapeCard = false, const TString& shapeFileName = "");
	//void printDataCard
	//void makeAllDataCards
	void addSyst(TH1D*& hist, const double* syst, const unsigned nSyst);
	void setSystUnc(TH1D** bkgYields, const unsigned nBkg, const TString* bkgNames);
	void makeTable(const unsigned cat, TH1D** bkgYields, TH1D** bkgErrors, TH1D* observed, const TString* bkgNames, const unsigned nBkg, const TString& catName, const TString& extraName = "");
}
#endif
