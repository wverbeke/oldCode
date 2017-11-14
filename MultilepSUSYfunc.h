#ifndef MultilepSUSY_func
#define MultilepSUSY_func

//Import ROOT classes
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TString.h"

//Import C++ library classes

//Function to calculate a transverse mass
inline double transmass(const TLorentzVector& v1, const TLorentzVector& v2){
	return sqrt(2*v1.Pt()*v2.Pt()*( 1 - cos( v1.Phi()-v2.Phi() ) ) );	//this definition assumes the ultarelativistic limit
}
//Calculate MT2 given two leptons and the MET
double mt2ll(const TLorentzVector& l1, const TLorentzVector& l2, const TLorentzVector& metVec); 
//Calculate cone corrected Pt for ra7 selection
double PtCone_ra7(double Pt, double Pt_rel, double Pt_ratio, int flavor, double miniisolation);
//Calculate cone corrected Pt for EWKino selection
double PtCone(double, int, double, double);
//Determine search category
unsigned SR_EWK_cat(const unsigned*, const int*, const double*, unsigned);
//Determine Mll indices
void mllIndices(unsigned*, const unsigned*, const TLorentzVector*, const double*, const int*, const unsigned);
//find the right leptons for the MT2 calulcation in EWKino analysis, and calculate MT2
double find_mt2(const unsigned*, const int*, const double*, const TLorentzVector*, const TLorentzVector&, const unsigned); //NOTE: maybe rename function to mt2_EWKino(...) 
//Calculate MT2 of two highest Pt leptons
double mt2_maxPt(const unsigned*, const double*, const TLorentzVector*, const TLorentzVector&, const unsigned);
//void SR_EWK_mll(const unsigned*, const TLorentzVector* , const int*, const double*, const unsigned*, const unsigned);
//Determine MET binning index
unsigned metSR(const double);
//Determine the search region for the EWKino analysis
//old ICHEP search regions
unsigned SR_EWK_3lOSSF_ICHEP(const double mt, const double met, const double mll);
unsigned SR_EWK_3lnoOSSF_ICHEP(const double, const double, const double);
unsigned SR_EWK_2ltauOSSF_ICHEP(const double mt2, const double, const double);
unsigned SR_EWK_2ltauOSOF_ICHEP(const double, const double, const double);
unsigned SR_EWK_2ltauSS_ICHEP(const double, const double, const double);
unsigned SR_EWK_l2tau_ICHEP(const double, const double, const double);
unsigned SR_EWK_ICHEP(const double, const double, const double, const unsigned cat);
//new search regions for Moriond 2017
//trilepton search regions
unsigned SR_EWK_3lOSSF(const double, const double, const double);
unsigned SR_EWK_3lnoOSSF(const double, const double, const double);
unsigned SR_EWK_2ltauOSSF(const double, const double, const double);
unsigned SR_EWK_2ltauSS(const double, const double, const double);
unsigned SR_EWK_l2tau(const double, const double, const double);
unsigned SR_EWK_3lep(const double, const double, const double, const unsigned);
//four lepton search regions
unsigned SR_EWK_4l2OSSF(const double met);
unsigned SR_EWK_4l1OSSF(const double);
unsigned SR_EWK_3ltau(const double);
unsigned SR_EWK_2l2tau2OSSF(const double);
unsigned SR_EWK_2l2tau1OSSF(const double);
unsigned SR_EWK_4lep(const double, const unsigned);
//same-sign dilepton search regions
unsigned SR_EWK_2lSS(const unsigned nJets, const double mtmin, const double ptll, const double met, const bool plusplus);
//Function to determine the flavor composition of the event
unsigned tril_flavorComb(const unsigned*, const int*, const unsigned);
//reweighting factor for non-prompt control region
double fakeWeight(const unsigned*, const int*, const double*, const double*, const bool*, TH2D**, const unsigned);
//return btag SF for event
double bTagWeight(const int*, const double*, const double*, const double*, const double*, const unsigned, TH2D**, const bool);
//Add systematic uncertainties to a histogram
void addSyst(TH1D*&, const double*, const unsigned);
void addSyst(TH2D*&, const double*, const unsigned);
//Return error histogram
//TH1D* errorHist(const TH1D*, const TString& histName);
//New advanced function to calculate Mll indices
void mllIndNew(unsigned*, const unsigned*, const TLorentzVector*, const double*, const int*, const unsigned, const double, const double);

//Print the EWKino yields in latex tables
//void printEWKTables(TH1D**, const unsigned, const unsigned, const TString&);


#endif
