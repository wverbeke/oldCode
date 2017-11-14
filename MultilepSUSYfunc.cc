//include ROOT classes
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "THistPainter.h"
#include "THStack.h"
#include "TColor.h"
//#include "TGraphAsymmErrors.h"
#include "TF1.h"
//include C++ library classes
#include <sstream>
#include <iostream>
#include <cmath>
#include <fstream>
#include <set>
//include other parts of the code
#include "MultilepSUSYfunc.h"
//include code to calculate mt2
#include "mt2_bisect.h"


double mt2ll(const TLorentzVector& l1, const TLorentzVector& l2, const TLorentzVector& metVec){
	return  asymm_mt2_lester_bisect::get_mT2(l1.M(), l1.Px(), l1.Py(), l2.M(), l2.Px(), l2.Py(), metVec.Px(), metVec.Py(), 0, 0);
}


double PtCone_ra7(double Pt, double Pt_rel, double Pt_ratio, int flavor, double miniisolation){
	const double multiConst[5][3] = { {0.25, 0.67, 4.4}, {0.2, 0.69, 6.}, {0.16, 0.76, 7.2}, {0.12, 0.8, 7.2}, {0.09, 0.84, 7.2} };
	double pt_cone = Pt;
    if (Pt_rel > multiConst[flavor][2]){
    	pt_cone = Pt*(1+TMath::Max(0., miniisolation - multiConst[flavor][0]));
    } else if( multiConst[flavor][1]*Pt/Pt_ratio > pt_cone){
    	pt_cone = multiConst[flavor][1]*Pt/Pt_ratio;
    }
	return pt_cone;
}

double PtCone(double Pt, int flavor, double mva, double Pt_ratio){
	const double mvaCut[2] = {0.5, -0.2};
	if(flavor == 2){
		return Pt;
	}
	if(mva > mvaCut[flavor] && Pt > 10){
		return Pt;
	}
	double corr;
	if(flavor == 0) corr = 0.85;
	else corr = 0.75;
	return Pt*corr/Pt_ratio;
}

//Function to determine the search category
unsigned SR_EWK_cat(const unsigned* ind, const int* flavors, const double* charges, unsigned lCount){
	if(lCount > 4) lCount = 4;
	unsigned ntau = 0;
	unsigned nOSSF = 0;
	bool OSOF = false;
	std::set<unsigned> usedlep;
	for(unsigned l = 0; l < lCount; ++l){
		if(flavors[ind[l]] == 2) ++ntau;
		if(flavors[ind[l]] != 2 || lCount > 3){
			for(unsigned k = l + 1; k < lCount; ++k){
				if(flavors[ind[k]] != 2 || lCount > 3){
					if(charges[ind[l]] != charges[ind[k]]){
						if(flavors[ind[l]] == flavors[ind[k]]){
							std::set<unsigned>::iterator indl = usedlep.find(ind[l]);
							std::set<unsigned>::iterator indk = usedlep.find(ind[k]);
							if(indk == usedlep.end() && indl == usedlep.end()){
								++nOSSF;
								usedlep.insert(ind[l]);
								usedlep.insert(ind[k]);
							}
						} else OSOF = true;
					}
				}
			}
		}
	}
	if(lCount == 3){
		if(ntau == 0 && nOSSF){
			return 0;
		} else if(ntau == 0){
			return 1;
		} else if(ntau == 1 && nOSSF){
			return 2;
		} else if(ntau == 1 && OSOF){
			return 3;
		} else if(ntau == 1){
			return 4;
		} else if(ntau == 2){
			return 5;
		} else{
			return 999;		//event with 3 taus, skip this
		}
	} else if(lCount == 4){
		if(ntau == 0 && nOSSF > 1){
			return 6;
		} else if(ntau == 0){
			return 7;
		} else if(ntau == 1){
			return 8;
		} else if(nOSSF > 1 && ntau < 3){
			return 9;
		} else if(ntau < 3){
			return 10;
		} else{
			return 999;
		}
	} else if(lCount == 2 && !nOSSF && !OSOF){
		return 11;
	} else{
		return 9999;
	}
}

//Determine mll indices for events with category != 0
void mllIndices(unsigned* mllI, const unsigned* ind, const TLorentzVector* lepV, const double* charges, const int* flavors, const unsigned lCount){
	const double mZll = 91;
	const double mZlt = 60;
	const double mZem = 50;
	double mass;
	double minDiff = 99999.;
	double minDiff_OSSF = 99999.;
	bool usedOSSF = false;
	for(unsigned l = 0; l < lCount - 1; ++l){
		for(unsigned k = l + 1; k < lCount; ++k){
			bool OSSF = false;
			if(charges[ind[l]] != charges[ind[k]]){	
				if(flavors[ind[l]] == flavors[ind[k]] && flavors[ind[l]] != 2){
					OSSF = true;
					usedOSSF = true;
				} else if(flavors[ind[l]] != flavors[ind[k]] && flavors[ind[l]] != 2 && flavors[ind[k]] != 2){
					mass =mZem;
				} else if(flavors[ind[l]] != flavors[ind[k]]){ //one of the leptons is a tau, tau pair is never considered
					mass =mZlt;
				} else{		//skip tau pairs
					continue;
				}
				if( (!usedOSSF && fabs( (lepV[l] + lepV[k]).M()  - mass) < minDiff) || (OSSF && fabs( (lepV[l] + lepV[k]).M()  - mZll) < minDiff_OSSF) ){
					if(OSSF) minDiff_OSSF = fabs( (lepV[k] + lepV[l]).M()  - mZll);
					else minDiff = fabs( (lepV[k] + lepV[l]).M()  - mass);
					mllI[0] = ind[l];
					mllI[1] = ind[k];
				}
			}
		}
	}
}
//Find MT2 given event lepton info
double find_mt2(const unsigned* ind, const int* flavors, const double* charges, const TLorentzVector* lepV, const TLorentzVector& metV, const unsigned lCount){
	unsigned leadlep = 99;
	unsigned leadtau = 99;
	for(unsigned l = 0; l < lCount - 1; ++l){
		for(unsigned k = l +1; k < lCount; ++k){
			if(flavors[ind[l]] != 2 && flavors[ind[k]] != 2){
				if(charges[ind[l]] != charges[ind[k]]){
					return mt2ll(lepV[l], lepV[k], metV);
				}
			} else if(flavors[ind[l]] == 2 && flavors[ind[k]] != 2){
				if(k < leadlep){	//This function assumed ind is Pt ordered!
					leadlep = k;
				}
				if(l < leadtau){
					leadtau = l;
				}
			} else if(flavors[ind[k]] ==2 && flavors[ind[l]] != 2){
				if(l < leadlep){
					leadlep = l;
				}
				if(k < leadtau){
					leadtau = k;
				}
			}
		}
	}
	return mt2ll(lepV[leadlep], lepV[leadtau], metV);
}

//Caculate MT2 of highest Pt lepton pair
double mt2_maxPt(const unsigned* ind, const double* charges, const TLorentzVector* lepV, const TLorentzVector& metV, const unsigned lCount){
	unsigned subleadlep = 99;
	for(unsigned l = 1; l < lCount; ++l){
		if(charges[ind[l]] != charges[ind[0]]){
			subleadlep = l;
			break;
		}
	}
	if(subleadlep == 99) return 0;
	if(subleadlep == 0) std::cerr << "Bug in the function, the same lepton is tagged as leading and subleading" << std::endl;
	else return mt2ll(lepV[0], lepV[subleadlep], metV);
}

//Determine the MET SR index
unsigned metSR(const double met){
	if(met < 100){
		return 0;
	} else if(met < 150){
		return 1;
	} else if(met < 200){
		return 2;
	} else{
		return 3;
	}
}

//ICHEP SEARCH REGIONS! USE NEW DEFINITIONS IN CURRENT CODE
unsigned SR_EWK_3lOSSF_ICHEP(const double mt, const double met, const double mll){ //3 light leptons, OSSF
	unsigned metI, mtI, mllI;
	metI = metSR(met);
	if(mt < 120){
		mtI = 0;
	} else if(mt < 160){
		mtI = 1;
	} else{
		mtI = 2;
	}
	if(mll < 75){
		mllI = 0;
	} else if(mll < 105){
		mllI = 1;
	} else{
		mllI = 2;
	}
	return metI + 4*mtI + 12*mllI;
}

unsigned SR_EWK_3lnoOSSF_ICHEP(const double mt, const double met, const double mll){ 	//3 light leptons, no OSSF
	unsigned metI, mtI, mllI;
	if(mll < 100){
		mllI = 0;
	} else{
		mllI = 1;
	}
	if(mt < 120){
		mtI = 0;
		if(met < 100){
			metI = 0; 
		} else{
			metI = 1;
		}
	} else{
		mtI = 1;
		metI = 0;
	}
	return metI + 2*mtI + 3*mllI;
}

unsigned SR_EWK_2ltauOSSF_ICHEP(const double mt2, const double met, const double mll){ //OSSF + tau
	unsigned metI, mt2I, mllI;
	if(mll < 75 || mll > 105){
		if(mt2 < 100){
			mt2I = 0;
			if(mll < 75){
				mllI = 0;
			} else{
				mllI = 2;
			}
			metI = metSR(met);
		} else{
			mt2I = 1;
			mllI = 0;
			if(met < 200){
				metI = 0;
			} else{
				metI = 1;
			}
		}
	} else{
		mllI = 1;
		mt2I = 0;
		metI =metSR(met);
	}
	return metI +4*mllI + 12*mt2I;
}

unsigned SR_EWK_2ltauOSOF_ICHEP(const double mt2, const double met, const double mll){ //OSOF + tau
	unsigned metI, mt2I, mllI;
	if(mt2 < 100){
		mt2I = 0;
		metI = metSR(met);
		if(mll < 60){
			mllI = 0;
		} else if(mll < 100){
			mllI = 1;
		} else{
			mllI = 2;
		}
	} else{
		mt2I = 1;
		mllI = 0;
		if(met < 200){
			metI = 0;
		} else{
			metI = 1;
		}
	}
	return metI + 4*mllI + 12*mt2I;
}
	
unsigned SR_EWK_2ltauSS_ICHEP(const double mt2, const double met, const double mll){ //SS + tau
	unsigned metI, mt2I, mllI;
	if(mt2 < 100){
		mt2I = 0;
		if(mll < 100){
			metI = metSR(met);
			if(mll < 60){
				mllI = 0;
			} else{
				mllI = 1;
			}
		} else{
			mllI = 2;
			metI = 0;
		}
	} else{
		mt2I = 1;
		metI = 0;
		mllI = 0;
	}
	return metI + 4*mllI + 9*mt2I;
}

unsigned SR_EWK_l2tau_ICHEP(const double mt2, const double met, const double mll){ //l + 2 tau
	unsigned metI, mt2I, mllI;
	if(mt2 < 100){
		mt2I = 0;
		if(mll < 100){
			mllI = 0;
		} else{
			mllI = 1;
		}
		if(met < 100){
			metI = 0;
		} else if(met < 150){
			metI = 1;
		} else{
			metI = 2;
		}
	} else{
		mt2I = 1;
		mllI = 0;
		if(met < 200){
			metI = 0;
		} else{
			metI = 1;
		}
	}
	return metI + 3*mllI + 6*mt2I;
}

//Function to return the search region index for the EWKino analysis
unsigned SR_EWK(const double mt, const double met, const double mll, const unsigned cat){
	switch(cat){
		case 0: return SR_EWK_3lOSSF_ICHEP(mt, met, mll);
		case 1: return SR_EWK_3lnoOSSF_ICHEP(mt, met, mll);
		case 2: return SR_EWK_2ltauOSSF_ICHEP(mt, met, mll);
		case 3: return SR_EWK_2ltauOSOF_ICHEP(mt, met, mll);
		case 4: return SR_EWK_2ltauSS_ICHEP(mt, met, mll);
		case 5: return SR_EWK_l2tau_ICHEP(mt, met, mll);
		default:{
			std::cerr <<"ERROR: non-existing category index" << std::endl;
			return 999;
		}
	}
}

//New search region definitions for Moriond 2017
unsigned SR_EWK_3lOSSF(const double mt, const double met, const double mll){ //3 light leptons, OSSF
	unsigned metI = metSR(met), sr = 0;
	if(mt < 100 || mt > 160){
		if(met > 550){
			metI = 6;
		} else if(met > 400){
			metI = 5;
		} else if(met > 250){
			metI = 4;
		}
	}
	if(mll < 75){
		if(mt < 100){
			if(metI > 4) metI = 4;
		} else if(mt < 160){
			sr += 5;
		} else{
			sr += 9;
			if(metI > 4) metI = 4;
		}
	} else if(mll < 105){
		sr += 14;
		if(mt > 100 && mt <160){
			sr += 7;
		} else if(mt > 160){
			sr += 11;
			if(metI > 5) metI = 5;
		}
	} else{
		sr += 31;
		if(mt < 100){
			if(metI > 4) metI = 4;
		} else if(mt < 160){
			sr += 5;
		} else{
			sr += 9;
			if(metI > 3) metI = 3;
		}
	}
	return sr + metI;
}

unsigned SR_EWK_3lnoOSSF(const double mt, const double met, const double mll){ //3 light leptons, no OSSF
	unsigned metI, mtI, mllI;
	if(mll < 100){
		mllI = 0;
	} else{
		mllI = 1;
	}
	if(mt < 120){
		mtI = 0;
		if(met < 100){
			metI = 0; 
		} else{
			metI = 1;
		}
	} else{
		mtI = 1;
		metI = 0;
	}
	return metI + 2*mtI + 3*mllI;
}


unsigned SR_EWK_2ltauOSSF(const double mt2, const double met, const double mll){ //OSSF + tau
	unsigned metI, sr = 0;
	if(mll < 75 || mll > 105){
		if(mt2 < 100){
			if(mll > 105) sr += 11;
			metI = metSR(met);
			if(met > 250) metI = 4;
		}	
		if(mt2 > 100){
			sr += 16;
			if(met < 200){
				metI = 0;
			} else{
				metI = 1;
			}
		}
	} else{
		sr += 5;
		if(met < 100){
			metI = 0;
		} else if(met < 150){
			metI = 1;
		} else if(met < 200){
			metI = 2;
		} else if(met < 300){
			metI = 3;
		} else if(met < 400){
			metI = 4;
		} else{
			metI = 5;
		}
	}
	return sr + metI;	
}

unsigned SR_EWK_2ltauOSOF(const double mt2, const double met, const double mll){ //OSOF + tau
	unsigned metI = metSR(met), sr = 0;
	if(mt2 < 100){
		if(mll < 100 && met > 250) metI = 4;
		if(mll > 60 && mll < 100){
			sr += 5;
		} else if(mll >100){
			sr += 10;
		}
	} else{
		sr += 14;
		if(met < 200){
			metI = 0;
		} else{
			metI = 1;
		}
	}
	return sr + metI;
}

unsigned SR_EWK_2ltauSS(const double mt2, const double met, const double mll){ //SS + tau
	unsigned metI = 0, sr = 0;
	if(mt2 < 100){
		if(mll < 100){
			metI = metSR(met);
			if(met > 250) metI = 4;
			if(mll > 60) sr += 5;
		} else{
			sr += 10;
		}
	} else{
		sr += 11;
	}
	return sr + metI;
}

unsigned SR_EWK_l2tau(const double mt2, const double met, const double mll){ //l + 2 tau
	unsigned metI = 0, sr = 0;
	if(mt2 < 100){
		metI = metSR(met);
		if(mll < 100){
			if(met > 300){
				metI = 5;
			} else if(met > 250){
				metI = 4;
			}
		} else{
			sr += 6;
		}
	} else{
		sr += 10;
		if(met < 200){
			metI = 0;
		} else{
			metI = 1;
		}
	}
	return sr + metI;
}

unsigned SR_EWK_3lep(const double mt, const double met, const double mll, const unsigned cat){
	switch(cat){
		case 0: return SR_EWK_3lOSSF(mt, met, mll);
		case 1: return SR_EWK_3lnoOSSF(mt, met, mll);
		case 2: return SR_EWK_2ltauOSSF(mt, met, mll);
		case 3: return SR_EWK_2ltauOSOF(mt, met, mll);
		case 4: return SR_EWK_2ltauSS(mt, met, mll);
		case 5: return SR_EWK_l2tau(mt, met, mll);
		default:{
			std::cerr <<"ERROR: non-existing category index" << std::endl;
			return 999;
		}
	}
}



//four lepton search region definitions:
unsigned metSR4l(const double met){
	if(met < 50){
		return 0;
	} else if(met < 100){
		return 1;
	} else if(met < 150){
		return 2;
	} else{
		return 3;
	}
}
unsigned SR_EWK_4l2OSSF(const double met){
	if(met > 200){
		return 4;
	} else{
		return metSR4l(met);
	}
}

unsigned SR_EWK_4l1OSSF(const double met){
	return metSR4l(met);
}

unsigned SR_EWK_3ltau(const double met){
	return metSR4l(met);
}

unsigned SR_EWK_2l2tau2OSSF(const double met){
	return metSR4l(met);
}

unsigned SR_EWK_2l2tau1OSSF(const double met){
	if(met > 100){
		return 2;
	} else{
		return metSR4l(met);
	}
}

unsigned SR_EWK_4lep(const double met, const unsigned cat){
	switch(cat){
		case 6: return SR_EWK_4l2OSSF(met);
		case 7: return SR_EWK_4l1OSSF(met);
		case 8: return SR_EWK_3ltau(met);
		case 9: return SR_EWK_2l2tau2OSSF(met);
		case 10: return SR_EWK_2l2tau1OSSF(met);
		default:{
			std::cerr << "ERROR: non-existing category index" << std::endl;
		}
	}
}

//same-sign dilepton search regions
unsigned SR_EWK_2lSS(const unsigned nJets, const double mtmin, const double ptll, const double met, const bool plusplus){
	unsigned mtminI = (mtmin > 100);
	unsigned nJetsI = (nJets == 1);
	if(nJets > 1) return 999;
	unsigned metI = metSR(met);
	unsigned ptllI = (mtmin > 100 && ptll > 50);
	unsigned plusplusI = 0;
	if(metI  == 1 && !plusplus) plusplusI = 1;
	return plusplusI + (metI > 2 ? metI: metI + 1) + 5*ptllI + 10*mtminI + 15*nJetsI;
}


//Function to return info about the flavor composition of the event (useful for triggering)
unsigned tril_flavorComb(const unsigned* ind, const int* flavors, const unsigned lCount){
	unsigned flavCount[3] = {0,0,0};
	for(unsigned l = 0; l < lCount; ++l) ++flavCount[flavors[ind[l]]];
	if(flavCount[2] == 0){
		if(flavCount[1] == 0){
			return 0; //eee
		} else if(flavCount[1] == 1){
			return 1; //eem
		} else if(flavCount[1] == 2){
			return 2; //emm
		} else{
			return 3; //mmm
		}
	} else if(flavCount[2] == 1){
		if(flavCount[1] == 0){
			return 4; //eet
		} else if(flavCount[1] == 1){
			return 5; //emt
		} else{
			return 6; //mmt
		}
	} else if(flavCount[2] == 2){
		if(flavCount[1] == 0){
			return 7; //ett
		} else{
			return 8; //mtt
		}
	} else{
		return 9; //ttt
	}
}


//Calculate the weight to apply to events in the non-prompt control region:
double fakeWeight(const unsigned* ind, const int* flavors, const double* conePt, const double* eta, const bool* tight, TH2D** frMap, const unsigned lCount){
	/*
	static double maxPt[3];
	static unsigned count = 0;
	if(count < 1){
		TH2D* frMapC[3] = {(TH2D*) frMap[0]->Clone(), (TH2D*) frMap[1]->Clone(), (frMap[2] == nullptr ? nullptr : (TH2D*) frMap[2]->Clone())};
		maxPt[0] = (frMapC[0]->ProjectionX())->GetBinCenter((frMapC[0]->ProjectionX())->GetNbinsX());
		maxPt[1] = (frMapC[1]->ProjectionX())->GetBinCenter((frMapC[1]->ProjectionX())->GetNbinsX());
		//if(frMapC[2] != nullptr) maxPt[2] = (frMapC[2]->ProjectionX())->GetBinCenter((frMapC[2]->ProjectionX())->GetNbinsX());
		++count;
	}
	*/
	double weight = -1;
	for(unsigned l = 0; l < lCount; ++l){
		if(!tight[ind[l]]){
			double fr = frMap[flavors[ind[l]]]->GetBinContent(frMap[flavors[ind[l]]]->FindBin(std::min(conePt[l],49.), fabs(eta[ind[l]])));
			//double fr = frMap[flavors[ind[l]]]->GetBinContent(frMap[flavors[ind[l]]]->FindBin(TMath::Min(conePt[l], 65.), TMath::Min(fabs(eta[ind[l]]), 2.4) ) );
			weight *= -fr/(1-fr);
		}
	}
	return weight;
}

//Calculate the weight to apply 
double bTagWeight(const int* jetFlavor, const double* jetPt, const double* jetEta, const double* jetCsv, const double* bTagSF, const unsigned jetCount, TH2D** bTagEff, const bool looseBtag){
	//Warning: the given jet collection is assumed correct, i.e. precleaned if cleaned jets are used
	double Pdata = 1.;
	double Pmc = 1.;
	for(unsigned j = 0;  j < jetCount; ++j){
		if(jetFlavor[j] != 0 && jetFlavor[j] != 4 && jetFlavor[j] != 5) continue;
		unsigned flav = 0 + (jetFlavor[j] == 4) + 2*(jetFlavor[j] == 5);
		if(jetCsv[j] > (looseBtag ? 0.5426 : 0.8484) ){
			Pdata *= bTagEff[flav]->GetBinContent(bTagEff[flav]->FindBin(std::max(25., std::min(jetPt[j], 599.)), std::min(fabs(jetEta[j]), 2.399) ))*bTagSF[j];
			Pmc *= bTagEff[flav]->GetBinContent(bTagEff[flav]->FindBin(std::max(25., std::min(jetPt[j], 599.)), std::min(fabs(jetEta[j]), 2.399) ));
		}
		else{
			Pdata *= (1. - bTagEff[flav]->GetBinContent(bTagEff[flav]->FindBin(std::max(25., std::min(jetPt[j], 599.)), std::min(fabs(jetEta[j]), 2.399) ))*bTagSF[j] );
			Pmc *= (1. - bTagEff[flav]->GetBinContent(bTagEff[flav]->FindBin(std::max(25., std::min(jetPt[j], 599.)), std::min(fabs(jetEta[j]), 2.399) )) );
		}
	}
	return (Pdata/Pmc);
}

void addSyst(TH1D*& hist, const double* syst, const unsigned nSyst){
	for(unsigned b = 1; b < hist->GetNbinsX() + 1; ++b){
		for(unsigned s = 0; s < nSyst; ++s){
			hist->SetBinError(b, sqrt(hist->GetBinError(b)*hist->GetBinError(b) + hist->GetBinContent(b)*hist->GetBinContent(b)*syst[s]*syst[s]) );
		}
	}
}

void addSyst(TH2D*& hist, const double* syst, const unsigned nSyst){
	for(unsigned x = 1; x < hist->GetNbinsX() + 1; ++x){
		for(unsigned y = 1; y < hist->GetNbinsY() + 1; ++y){
			for(unsigned s = 0; s < nSyst; ++s){
				hist->SetBinError(x,y, sqrt(hist->GetBinError(x,y)*hist->GetBinError(x,y) + hist->GetBinContent(x,y)*hist->GetBinContent(x,y)*syst[s]*syst[s]) );
			}
		}
	}
}

/*
TH1D* errorHist(const TH1D* hist, const TString& histName){
	TH1D* errorHist = new TH1D(histName + "error", histName + "error", hist->GetNbinsX(), hist->GetBinLowEdge(1), hist->GetBinLowEdge(data->GetNbinsX()) + hist->GetBinWidth(data->GetNbinsX()));
	for(unsigned b = 1; b < hist->GetNbinsX() + 1; ++b){
		errorHist-SetBinContent(b, hist->GetBinError(b);
	}
	return errorHist;
}
*/



void mllIndNew(unsigned* mllI, const unsigned* ind, const TLorentzVector* lepV, const double* charges, const int* flavors, const unsigned lCount, const double met, const double metphi){
	const double mZ = 91;
	//const double mGS = 20;
	/*
	std::set<std::tuple<unsigned, unsigned, double>> pairs;
	for(unsigned l = 0; l < lCount; ++l){
		for(unsigned k = 0; k < lCount; ++k){
			pairs.insert(std::make_tuple(ind[l], ind[k], (lepV[l] + lepV[k]).M() ) );
		}
	}
	for(std::set<std::tuble<unsigned, unsigned, double>> iterator it = pairs.begin()
	*/
	//
	std::set<double> masses;
	for(unsigned l = 0; l < lCount -1 ; ++l){
		for(unsigned k = l + 1; k < lCount; ++k){
			masses.insert( (lepV[l] + lepV[k]).M() );
		}
	}
	//find min mas and mass closest to Z
	double min = 99999;
	double zcand = 99999.;
	for(std::set<double>::iterator it = masses.begin(); it != masses.end(); ++it){
		if(*it < min) min = *it;
		if(fabs(*it - mZ) < fabs(mZ - zcand) ) zcand = *it;
	}
	double best;
	//if(min > 25) best = mZ;
	//if(min < 30) best = min;
	//else best = mZ;
	//std::cout << best << std::endl;
	if(fabs(zcand - mZ) < 15){
		best = mZ;
		double minDiff = 99999;
		for(unsigned l = 0; l < lCount - 1; ++l){
			for(unsigned k = l + 1; k < lCount; ++k){
				if(charges[ind[l]] != charges[ind[k]]){
					if(flavors[ind[l]] != flavors[ind[k]]) continue;
					if(fabs(best - (lepV[k] + lepV[l]).M() ) < minDiff){
						mllI[0] = ind[k];
						mllI[1] = ind[l];
						minDiff = fabs(best - (lepV[k] + lepV[l]).M() );
					}
				}
			}
		}
	} else{
		double minMT = 0;
		TLorentzVector metV;
		unsigned lw = 999;
		metV.SetPtEtaPhiE(met, 0, metphi, met);
		for(unsigned l = 0; l < lCount; ++l){
			bool ossfcheck = false;
			for(unsigned k = 0; k < lCount -1; ++k){
				if(k == l) continue;
				for(unsigned j = k + 1; j < lCount; ++j){
					if(j == l) continue;
					if(charges[ind[k]] != charges[ind[j]] && flavors[ind[k]] == flavors[ind[j]]){
						ossfcheck = true;
						break;
					}
				}
				if(ossfcheck) break;
			}
			if(ossfcheck){
				if(transmass(lepV[l], metV) < minMT || minMT == 0){
					minMT = transmass(lepV[l], metV);
					lw = ind[l];
				}
			}
			else{
				best = mZ;
				double minDiff = 99999;
				for(unsigned l = 0; l < lCount - 1; ++l){
					for(unsigned k = l + 1; k < lCount; ++k){
						if(charges[ind[l]] != charges[ind[k]]){
							if(flavors[ind[l]] != flavors[ind[k]]) continue;
							if(fabs(best - (lepV[k] + lepV[l]).M() ) < minDiff){
								mllI[0] = ind[k];
								mllI[1] = ind[l];
								minDiff = fabs(best - (lepV[k] + lepV[l]).M() );
							}
						}
					}
				}
			}
		}
		unsigned zCount = 0;
		for(unsigned l = 0; l < lCount; ++l){
			if(ind[l] != lw){
				mllI[zCount] = ind[l];
				++zCount;
			}
		}
	}
	
	/*
	if(fabs(zcand - mZ) > 5){
		double minDiff = 99999;
		for(unsigned l = 0; l < lCount - 1; ++l){
			for(unsigned k = l + 1; k < lCount; ++k){
				if(charges[ind[l]] != charges[ind[k]]){
					if(flavors[ind[l]] != flavors[ind[k]]) continue;
					if( lepV[k].DeltaR(lepV[l]) < minDiff){
						mllI[0] = ind[k];
						mllI[1] = ind[l];
						minDiff = lepV[k].DeltaR(lepV[l]);
					}
				}
			}
		}
	} else{
		best = 91;	
		double minDiff = 99999;
		for(unsigned l = 0; l < lCount - 1; ++l){
			for(unsigned k = l + 1; k < lCount; ++k){
				if(charges[ind[l]] != charges[ind[k]]){
					if(flavors[ind[l]] != flavors[ind[k]]) continue;
					if(fabs(best - (lepV[k] + lepV[l]).M() ) < minDiff){
						mllI[0] = ind[k];
						mllI[1] = ind[l];
						minDiff = fabs(best - (lepV[k] + lepV[l]).M() );
					}
				}
			}
		}
	}
	*/
	/*
	double minDiff = 99999;
	for(unsigned l = 0; l < lCount - 1; ++l){
		for(unsigned k = l + 1; k < lCount; ++k){
			if(charges[ind[l]] != charges[ind[k]]){
				if(flavors[ind[l]] != flavors[ind[k]]) continue;
				if(fabs(best - (lepV[k] + lepV[l]).M() ) < minDiff){
					mllI[0] = ind[k];
					mllI[1] = ind[l];
					minDiff = fabs(best - (lepV[k] + lepV[l]).M() );
				}
			}
		}
	}
	*/
}




void printProgress(const float progress){
	
}




/*
void printEWKTables(TH1D** yields, const unsigned nBkg, const unsigned cat, const TString& name){
	const unsigned nSR[6] = {36, 6, 14, 14, 10, 8};
	const unsigned nMllb[6] = {3, 2, 3, 3, 3, 2};
	//Merge backgrounds for total yields
	TH1D* bkgTot = (TH1D*) yields[1]->Clone();
	for(unsigned bkg = 2; bkg < nBkg + 1; ++bkg){
		bkgTot->Add(yields[bkg]);
	}
	//Write latex table with all the yields to tex file
	std::ofstream table;
	table.open( "tables/" + name + "_yields.txt");
	table << "\\begin{tabular}{|c|c|";
	for(unsigned mll = 0; mll < nMllb[cat]; ++mll) table <<"c|";
	table <<"|}" << std::endl;
	table << "\\hline" << std::endl;
	if(cat < 2) table << "$M_{T}$ (GeV)";
	else table << "M_{T2}^{\\ell \\ell} (GeV)";
	table << " & $E_{T}^{miss}$ (GeV) & ";
	if(cat == 0){
		table << "$M_{\\ell \\ell} < 75 \\mathrm{GeV}$ & $75 < M_{\\ell \\ell} < 105 \\mathrm{GeV}$ & $M_{\\ell \\ell} > 105 \\mathrm{GeV}$ \\\\" << std::endl;
		table << " \\hline \\hline" << std::endl;
		TString metBins[4] = {"50-100", "100-150", "150-200", "> 200"};
		TString mtBins[3] = {"0-120", "120-160", "> 160"};
		unsigned j = 0;
		for(unsigned sr1 = 0; sr1 < 12; ++sr1){
			unsigned metI = sr1%4;
			unsigned mtI = sr1/4;
			if(sr1/4 >= j){
				table << "\\multirow{4}{*}{$" << mtBins[j] << "$}";
				++j;
			}
			table << "& $" << metBins[metI] << "$";
			unsigned i = 0;
			for(unsigned sr2 = 0; sr2 < nSR[cat]; ++sr2){
				unsigned mllI = sr2/12;
				unsigned bin = 1 +metI + 4*mtI + 12*mllI;
				if(sr2/12 >= i){
					table <<" & " << bkgTot->GetBinContent(bin) << " $\\pm$ " << bkgTot->GetBinError(bin) << "    " << yields[0]->GetBinContent(bin) << " $\\pm$ " <<  yields[0]->GetBinError(bin);
					++i;
				}
			}
			table << " \\\\ \\hline" << std::endl;
		}
	} else if(cat == 1){
		table << "$M_{\\ell \\ell} < 100 \\mathrm{GeV}$ & $M_{\\ell \\ell} > 100 \\mathrm{GeV}$ \\\\" << std::endl;
		table << "\\hline \\hline" << std::endl;
		table << "\\multirow{4}{*}{$ < 120$} & $50-100$ & " << bkgTot->GetBinContent(0) << " $\\pm$ " << bkgTot->GetBinError(0) << "    " << yields[0]->GetBinContent(0) << " $\\pm$ " <<  yields[0]->GetBinError(0) << " & ";
		table << bkgTot->GetBinContent(3) << " $\\pm$ " << bkgTot->GetBinError(3) << "    " << yields[0]->GetBinContent(3) << " $\\pm$ " <<  yields[0]->GetBinError(3) << " \\\\ \\hline" << std::endl;
		table << "& $	>100$ & " << bkgTot->GetBinContent(1) << " $\\pm$ " << bkgTot->GetBinError(1) << "    " << yields[0]->GetBinContent(1) << " $\\pm$ " <<  yields[0]->GetBinError(1) << " & ";
		table << bkgTot->GetBinContent(4) << " $\\pm$ " << bkgTot->GetBinError(4) << "    " << yields[0]->GetBinContent(4) << " $\\pm$ " <<  yields[0]->GetBinError(4) << " \\\\ \\hline" << std::endl;
		table << "\\multirow{1}{*}{$ > 120$} & $>50$ & " << bkgTot->GetBinContent(2) << " $\\pm$ " << bkgTot->GetBinError(2) << "    " << yields[0]->GetBinContent(2) << " $\\pm$ " <<  yields[0]->GetBinError(2) << " & ";
		table << "& " << bkgTot->GetBinContent(5) << " $\\pm$ " << bkgTot->GetBinError(5) << "    " << yields[0]->GetBinContent(5) << " $\\pm$ " <<  yields[0]->GetBinError(5) << " \\\\ \\hline" << std::endl;
	} else if(cat == 2 || cat == 3){
		if(cat == 2) table << "$M_{\\ell \\ell} < 75 \\mathrm{GeV}$ & $75 < M_{\\ell \\ell} < 105 \\mathrm{GeV}$ & $M_{\\ell \\ell} > 105 \\mathrm{GeV}$ \\\\" << std::endl;
		else table << "$M_{\\ell \\ell} < 60 \\mathrm{GeV}$ & $60 < M_{\\ell \\ell} < 100 \\mathrm{GeV}$ & $M_{\\ell \\ell} > 100 \\mathrm{GeV}$ \\\\" << std::endl;
		table << "\\hline \\hline" << std::endl;
		table << "\\multirow{4}{*}{$ < 100$} & $50-100$ & " << bkgTot->GetBinContent(0) << " $\\pm$ " << bkgTot->GetBinError(0) << "    " << yields[0]->GetBinContent(0) << " $\\pm$ " <<  yields[0]->GetBinError(0) << " & ";
		table << bkgTot->GetBinContent(4) << " $\\pm$ " << bkgTot->GetBinError(4) << "    " << yields[0]->GetBinContent(4) << " $\\pm$ " <<  yields[0]->GetBinError(4) << " & ";
		table << bkgTot->GetBinContent(8) << " $\\pm$ " << bkgTot->GetBinError(8) << "    " << yields[0]->GetBinContent(8) << " $\\pm$ " <<  yields[0]->GetBinError(8) << " \\\\ \\hline" << std::endl;
		table << "& $100-150$ & " << bkgTot->GetBinContent(1) << " $\\pm$ " << bkgTot->GetBinError(1) << "    " << yields[0]->GetBinContent(1) << " $\\pm$ " <<  yields[0]->GetBinError(1) << " & ";
		table << bkgTot->GetBinContent(5) << " $\\pm$ " << bkgTot->GetBinError(5) << "    " << yields[0]->GetBinContent(5) << " $\\pm$ " <<  yields[0]->GetBinError(5) << " & ";
		table << bkgTot->GetBinContent(9) << " $\\pm$ " << bkgTot->GetBinError(9) << "    " << yields[0]->GetBinContent(9) << " $\\pm$ " <<  yields[0]->GetBinError(9) << " \\\\ \\hline" << std::endl;
		table << "& $150-200$ & " << bkgTot->GetBinContent(2) << " $\\pm$ " << bkgTot->GetBinError(2) << "    " << yields[0]->GetBinContent(2) << " $\\pm$ " <<  yields[0]->GetBinError(2) << " & ";
		table << bkgTot->GetBinContent(6) << " $\\pm$ " << bkgTot->GetBinError(6) << "    " << yields[0]->GetBinContent(6) << " $\\pm$ " <<  yields[0]->GetBinError(6) << " & ";
		table << bkgTot->GetBinContent(10) << " $\\pm$ " << bkgTot->GetBinError(10) << "    " << yields[0]->GetBinContent(10) << " $\\pm$ " <<  yields[0]->GetBinError(10) << " \\\\ \\hline" << std::endl;
		table << "& $>200$ & " << bkgTot->GetBinContent(3) << " $\\pm$ " << bkgTot->GetBinError(3) << "    " << yields[0]->GetBinContent(3) << " $\\pm$ " <<  yields[0]->GetBinError(3) << " & ";
		table << bkgTot->GetBinContent(7) << " $\\pm$ " << bkgTot->GetBinError(7) << "    " << yields[0]->GetBinContent(7) << " $\\pm$ " <<  yields[0]->GetBinError(7) << " & ";
		table << bkgTot->GetBinContent(11) << " $\\pm$ " << bkgTot->GetBinError(11) << "    " << yields[0]->GetBinContent(11) << " $\\pm$ " <<  yields[0]->GetBinError(11) << " \\\\ \\hline" << std::endl;
		table << "\\multirow{2}{*}{$ > 100$} & $50-200$ & \\multicolumn{3}{c}{" << bkgTot->GetBinContent(12) << " $\\pm$ " << bkgTot->GetBinError(12) << "    " << yields[0]->GetBinContent(12) << " $\\pm$ " <<  yields[0]->GetBinError(12) << "} & ";
		table << "\\multicolumn{3}{c}{" << bkgTot->GetBinContent(13) << " $\\pm$ " << bkgTot->GetBinError(13) << "    " << yields[0]->GetBinContent(13) << " $\\pm$ " <<  yields[0]->GetBinError(13) << "} \\\\ \\hline" << std::endl;

	} else if(cat == 4){
		table << "$M_{\\ell \\ell} < 60 \\mathrm{GeV}$ & $60 < M_{\\ell \\ell} < 100 \\mathrm{GeV}$ & $M_{\\ell \\ell} > 100 \\mathrm{GeV}$ \\\\" << std::endl;
		table << "\\hline \\hline" << std::endl;
		table << "\\multirow{4}{*}{$ < 100$} & $50-100$ & " << bkgTot->GetBinContent(0) << " $\\pm$ " << bkgTot->GetBinError(0) << "    " << yields[0]->GetBinContent(0) << " $\\pm$ " <<  yields[0]->GetBinError(0) << " & ";
		table << bkgTot->GetBinContent(4) << " $\\pm$ " << bkgTot->GetBinError(4) << "    " << yields[0]->GetBinContent(4) << " $\\pm$ " <<  yields[0]->GetBinError(4) << " & ";
		table << "\\multirow{3}{*}{" << bkgTot->GetBinContent(8) << " $\\pm$ " << bkgTot->GetBinError(8) << "    " << yields[0]->GetBinContent(8) << " $\\pm$ " <<  yields[0]->GetBinError(8) << "} \\\\ \\hline" << std::endl;
		table << "& $100-150$ & " << bkgTot->GetBinContent(1) << " $\\pm$ " << bkgTot->GetBinError(1) << "    " << yields[0]->GetBinContent(1) << " $\\pm$ " <<  yields[0]->GetBinError(1) << " & ";
		table << bkgTot->GetBinContent(5) << " $\\pm$ " << bkgTot->GetBinError(5) << "    " << yields[0]->GetBinContent(5) << " $\\pm$ " <<  yields[0]->GetBinError(5) << " & ";
		table << "& $150-200$ & " << bkgTot->GetBinContent(2) << " $\\pm$ " << bkgTot->GetBinError(2) << "    " << yields[0]->GetBinContent(2) << " $\\pm$ " <<  yields[0]->GetBinError(2) << " & ";
		table << bkgTot->GetBinContent(6) << " $\\pm$ " << bkgTot->GetBinError(6) << "    " << yields[0]->GetBinContent(6) << " $\\pm$ " <<  yields[0]->GetBinError(6) << " & ";
		table << "& $>200$ & " << bkgTot->GetBinContent(3) << " $\\pm$ " << bkgTot->GetBinError(3) << "    " << yields[0]->GetBinContent(3) << " $\\pm$ " <<  yields[0]->GetBinError(3) << " & ";
		table << bkgTot->GetBinContent(7) << " $\\pm$ " << bkgTot->GetBinError(7) << "    " << yields[0]->GetBinContent(7) << " $\\pm$ " <<  yields[0]->GetBinError(7) << " & ";
		table << "\\multirow{1}{*}{$ > 100$} & $50-200$ & \\multicolumn{3}{c}{" << bkgTot->GetBinContent(9) << " $\\pm$ " << bkgTot->GetBinError(9) << "    " << yields[0]->GetBinContent(9) << " $\\pm$ " <<  yields[0]->GetBinError(9) << "} \\\\ \\hline " << std::endl;
	} else if(cat == 5){
		table << "$M_{\\ell \\ell} < 100 \\mathrm{GeV}$ & $M_{\\ell \\ell} > 100 \\mathrm{GeV}$ \\\\" << std::endl;
		table << "\\hline \\hline" << std::endl;
		table << "\\multirow{3}{*}{$ < 100$} & $50-100$ & " << bkgTot->GetBinContent(0) << " $\\pm$ " << bkgTot->GetBinError(0) << "    " << yields[0]->GetBinContent(0) << " $\\pm$ " <<  yields[0]->GetBinError(0) << " & ";
		table << bkgTot->GetBinContent(3) << " $\\pm$ " << bkgTot->GetBinError(3) << "    " << yields[0]->GetBinContent(3) << " $\\pm$ " <<  yields[0]->GetBinError(3) << " \\\\ \\hline" << std::endl;
		table << "& $>100$ & " << bkgTot->GetBinContent(1) << " $\\pm$ " << bkgTot->GetBinError(1) << "    " << yields[0]->GetBinContent(1) << " $\\pm$ " <<  yields[0]->GetBinError(1) << " & ";
		table << bkgTot->GetBinContent(4) << " $\\pm$ " << bkgTot->GetBinError(4) << "    " << yields[0]->GetBinContent(4) << " $\\pm$ " <<  yields[0]->GetBinError(4) << " \\\\ \\hline" << std::endl;
		table << "& $>150$ & " << bkgTot->GetBinContent(2) << " $\\pm$ " << bkgTot->GetBinError(2) << "    " << yields[0]->GetBinContent(2) << " $\\pm$ " <<  yields[0]->GetBinError(2) << " & ";
		table << bkgTot->GetBinContent(5) << " $\\pm$ " << bkgTot->GetBinError(5) << "    " << yields[0]->GetBinContent(5) << " $\\pm$ " <<  yields[0]->GetBinError(5) << " \\\\ \\hline" << std::endl;
		table << "\\multirow{2}{*}{$ > 100$} & $50-200$ & \\multicolumn{3}{c}{" << bkgTot->GetBinContent(6) << " $\\pm$ " << bkgTot->GetBinError(6) << "    " << yields[0]->GetBinContent(6) << " $\\pm$ " <<  yields[0]->GetBinError(6) << "} & ";
		table << "\\multicolumn{3}{c}{" << bkgTot->GetBinContent(7) << " $\\pm$ " << bkgTot->GetBinError(7) << "    " << yields[0]->GetBinContent(7) << " $\\pm$ " <<  yields[0]->GetBinError(7) << "} \\\\ \\hline" << std::endl;			
	} else{
		std::cerr <<"ERROR: non-existing category index" << std::endl;
		return;
	}
	table << "\\end{tabular}" << std::endl; 
	table.close();
}
*/			











