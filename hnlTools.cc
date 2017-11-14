#include "hnlTools.h"

//include c++ library tools
#include <fstream>
#include <iostream>

//Return search category
unsigned hnl::cat(const unsigned* ind, const int* flavors, const double* charges, const unsigned lCount, const double ptlead){
	if(lCount > 3) return 999;
	unsigned flavorComp = hnl::flavorComposition(ind, flavors, charges, lCount);
	if(flavorComp > 1) return 999;
	return flavorComp + 2*(ptlead >= 30 && ptlead < 55) + 4*(ptlead >= 55);
}

unsigned hnl::flavorComposition(const unsigned* ind, const int* flavors, const double* charges, const unsigned lCount){
	bool OSOF = false;
	for(unsigned l = 0; l < lCount - 1; ++l){
		if(flavors[ind[l]] == 2) continue;
		for(unsigned k = l + 1; k < lCount; ++k){
			if(flavors[ind[k]] == 2) continue;
			if(charges[ind[l]] != charges[ind[k]]){
				if(flavors[ind[l]] == flavors[ind[k]]){
					return 0;
				}
				else{
					OSOF = true;
				}
			}
		}
	}
	return (1 + !OSOF);	//0 = OSSF, 1 = OSOF, 2 = SSS
}
/*
unsigned hnl::controlRegion(const unsigned* ind, const int* flavors, const double* charges, const unsigned lCount, const double m3l, const double mll, const double met, const unsigned nBJets){
	unsigned flavorComp = hnl::flavorComposition(ind, flavors, charges, lCount);
	if(nBJets == 0){
		if(lCount == 3 && fabs(mll - 91.) < 15 && fabs(m3l - 91.) > 15 && flavorComp == 0){
			if(met > 50) return 0;
			else return 3;
		} else if(fabs(mll - 91.) < 15 && lCount == 4 && flavorComp == 0){
			return 1;
		} else if(fabs(mll - 91.) > 15 && fabs(m3l - 91.) < 15 && flavorComp == 0){
			return 2;
		}
	} else if(flavorComp != 0 || (fabs(mll - 91.) > 15 && fabs(m3l - 91.) > 15)){
		return 4;
	} return 999;
}
*/

unsigned hnl::controlRegion(const unsigned* ind, const int* flavors, const double* charges, const unsigned lCount, const double m3l, const double mll){
	unsigned flavorComp = hnl::flavorComposition(ind, flavors, charges, lCount);
	if(lCount == 3 && fabs(mll - 91.) < 15 && fabs(m3l - 91.) > 15 && flavorComp == 0){
		return 0;
	} else if(fabs(mll - 91.) < 15 && lCount == 4 && flavorComp == 0){
		return 1;
	} else if(fabs(mll - 91.) > 15 && fabs(m3l - 91.) < 15 && flavorComp == 0){
		return 2;
	}
	return 999;
}

//Function to return search region number in each category
unsigned hnl::sr_lowM_3lOSSF_lowPt(const double mt, const double minMos){
	unsigned sr = 0;
	if(mt < 40) sr = 0;
	else if(mt < 80) sr = 1;
	else sr = 2;
	if(minMos < 10) return sr;
	else if(minMos < 20) return 3 + sr;
	else if(minMos < 30) return 6 + sr;
	else return 9 + sr;
}
unsigned hnl::sr_lowM_3lOSSF_highPt(const double mt, const double minMos){
	unsigned sr = 0;
	if(mt < 40) sr = 0;
	else if(mt < 80) sr = 1;
	else sr = 2;
	if(minMos < 10) return sr;
	else if(minMos < 20) return 3 + sr;
	else if(minMos < 30) return 6 + sr;
	else return 9 + sr;
}

unsigned hnl::sr_lowM_3lnoOSSF_lowPt(const double minMos){
	if(minMos < 10) return 0;
	else if(minMos < 20) return 1;
	else if(minMos < 30) return 2;
	else return 3;
}

unsigned hnl::sr_lowM_3lnoOSSF_highPt(const double minMos){
	return sr_lowM_3lnoOSSF_lowPt(minMos);
}

unsigned hnl::sr_highM_3lOSSF(const double mt, const double minMos, const double m3l){
	return hnl::sr_3lOSSF_mtminMos(mt, minMos, m3l);
}

unsigned hnl::sr_highM_3lnoOSSF(const double mt, const double minMos, const double m3l){
	return hnl::sr_3lnoOSSF_mtminMos(mt, minMos, m3l);
}

unsigned hnl::sr_3lOSSF_mtminMos(const double mt, const double minMos, const double m3l){
	unsigned sr = 0;
	for(unsigned mti = 0; mti < 4; ++mti){
		if(mt < (100. + mti*100.)){
			sr += mti;
			break;
		}
	}
	if(mt > 400) sr += 4;
	if(m3l < 100){
		if(sr > 2) return 2;
		return sr;
	}
	if(minMos < 100) return 3 + sr;
	if(minMos > 100 && minMos < 200){
		if(sr > 3) sr = 3;
		return 8 + sr;
	}	
	else{
		/*
		if(mt >= 200 && mt < 300) sr = 2;
		else if(mt > 300) sr = 3;
		*/
		if(sr > 3) sr = 3;
		return sr + 12;
	}
}


unsigned hnl::sr_3lnoOSSF_mtminMos(const double mt, const double minMos, const double m3l){
	if(m3l < 100){
		if(mt < 100) return 0;
		else return 1;
	}
	if(minMos < 100){	
		if(mt < 100) return 2;
		if(mt < 150) return 3; //200
		if(mt < 250) return 4;
		return 5;
	}
	if(minMos < 200){
		if(mt < 100) return 6;
		else return 7;
	} else{
		return 8;
	}
}


unsigned hnl::sr(const double mt, const double minMos, const double m3l, const unsigned cat){
	switch(cat){
		case 0: return sr_lowM_3lOSSF_lowPt(mt, minMos);
		case 1: return sr_lowM_3lnoOSSF_lowPt(minMos);
		case 2: return sr_lowM_3lOSSF_highPt(mt, minMos);
		case 3: return sr_lowM_3lnoOSSF_highPt(minMos);
		case 4: return sr_highM_3lOSSF(mt, minMos, m3l);
		case 5: return sr_highM_3lnoOSSF(mt, minMos, m3l);
		default:{
			std::cerr << "ERROR: Undefined category index given" << std::endl;
			return 999;
		}
	}
}


//Function to print dataCard to be analysed by CMS combination tool
void hnl::printDataCard(const double obsYield, const double sigYield, const TString& sigName, const double* bkgYield, const unsigned nBkg, const TString* bkgNames, const std::vector<std::vector<double> >& systUnc, const unsigned nSyst, const TString* systNames, const TString* systDist, const TString& cardName, const bool shapeCard, const TString& shapeFileName){
	std::ofstream card;
	card.open(cardName + ".txt");
	//define number of channels, background sources and systematics
	card << "imax 1 number of channels \n";
	card << "jmax " << nBkg << " number of backgrounds \n";
	card << "kmax " << nSyst << " number of nuisance parameters (sources of systematical uncertainties) \n";
	card << "---------------------------------------------------------------------------------------- \n";
	//define the channels and the number of observed events
	card << "bin 1 \n";
	card << "observation " << obsYield << "\n";
	//define all backgrounds and their yields
	card << "---------------------------------------------------------------------------------------- \n";
	if(shapeCard){
		/*
		if(shapeFileName == "") shapeFileName = (std::string) cardName;
		shapeFileName = shapeFileName.substr(shapeFileName.find("/") + 1, shapeFileName.length());
		*/
		//card << "shapes * * " << shapeFileName + ".root  $PROCESS_$CHANNEL $PROCESS_$CHANNEL_$SYSTEMATIC\n";
		card << "shapes * * " << shapeFileName + ".root  $PROCESS $PROCESS_$SYSTEMATIC\n";
		card << "---------------------------------------------------------------------------------------- \n";
	}
	card << "bin	";
	for(unsigned proc = 0; proc < nBkg + 1; ++proc){
		card << "	" << 1;
	}
	card << "\n";
	card << "process";
	card << "	" << sigName;
	for(unsigned bkg = 0; bkg < nBkg; ++bkg){
		card << "	" << bkgNames[bkg];
	}
	card << "\n";
	card << "process";
	for(unsigned bkg = 0; bkg < nBkg + 1; ++bkg){
		card << "	" << bkg;
	}
	card << "\n";
	card <<	"rate";
	card << "	" << sigYield;
	for(unsigned bkg = 0; bkg < nBkg; ++bkg){
		if(bkgYield[bkg] <= 0) card << "	" << "0.00";
		else card << "	" << bkgYield[bkg];
	}
	card << "\n";
	card << "---------------------------------------------------------------------------------------- \n";
	//define sources of systematic uncertainty, what distibution they follow and how large their effect is
	for(unsigned syst = 0; syst < nSyst; ++syst){
		card << systNames[syst] << "	" << systDist[syst];
		/*
		if(systDist[syst] = "gmN"){
			for(unsigned bkg = 0; bkg < nBkg + 1; ++bkg){
				if(systUnc[syst][proc] != 0){
					card << systDist[syst] << " " << bkgYield[bkg];
				}
			}
		}
		*/
		for(unsigned proc = 0; proc < nBkg + 1; ++proc){
			card << "	";
			if(systUnc[syst][proc] == 0) card << "-";
			else card << systUnc[syst][proc];
		}
		card << "\n";
	}
	card.close();		
}


void hnl::addSyst(TH1D*& hist, const double* syst, const unsigned nSyst){
	for(unsigned b = 1; b < hist->GetNbinsX() + 1; ++b){
		for(unsigned s = 0; s < nSyst; ++s){
			hist->SetBinError(b, sqrt(hist->GetBinError(b)*hist->GetBinError(b) + hist->GetBinContent(b)*hist->GetBinContent(b)*syst[s]*syst[s]) );
		}
	}
}



void hnl::setSystUnc(TH1D** bkgYields, const unsigned nBkg, const TString* bkgNames){
	double syst[6] = {0.025,	//lumi
						0.04, //id eff
						0.02, //trig eff
						0.05, //JEC
						0.05, //PU
						0.    //Extra uncertainty varies per background
						};
	std::map<const TString, double> bkgExtraUnc;
	bkgExtraUnc["WZ"] = 0.1;
	bkgExtraUnc["X + gamma"] = 0.15;
	bkgExtraUnc["TT/T + X"] = 0.15;
	bkgExtraUnc["triboson"] = 0.5;
	bkgExtraUnc["ZZ/H"] = 0.25;
	bkgExtraUnc["non-prompt"] = 0.3;
	for(unsigned bkg = 0; bkg < nBkg; ++bkg){
		syst[5] = bkgExtraUnc[bkgNames[bkg]];
		hnl::addSyst(bkgYields[bkg], syst, 6);
	}
}



int precision(double x){
	int prec = 0;
	if(x == 0){;} 
	else if(x < 1.){
		--prec;
		while(x*10 < 1.){
			--prec;
			x*=10;
		}
	} else{
		while(x/10. > 1.){
			++prec;
			x/= 10.;
		}
	}
	return prec;
}



void hnl::makeTable(const unsigned cat, TH1D** bkgYields, TH1D** bkgErrors, TH1D* observed, const TString* bkgNames, const unsigned nBkg, const TString& catName, const TString& extraName){
	if(catName.Contains("lowM_3lOSSF")) return;
	TH1D* promptYield = (TH1D*) bkgYields[0]->Clone();
	TH1D* nonPromptYield = (TH1D*) bkgYields[nBkg -1]->Clone();
	for(unsigned bkg = 0; bkg < nBkg; ++bkg){ //nBkg -1 if split NP
		promptYield->Add(bkgYields[bkg]);
	}
	/*
	//Calculate Total stat unc
	TH1D* promptStat = (TH1D*) promptYield->Clone();
	TH1D* nonPromptStat = (TH1D*) nonPromptYield->Clone();
	for(unsigned b = 1; b < promptYield->GetNbinsX() + 1; ++b){
		promptStat->SetBinContent(b, promptYield->GetBinError(b));
		nonPromptStat->SetBinContent(b, nonPromptYield->GetBinError(b));
	}
	//Make copies of bkg histograms and remove their stat unc
	TH1D* bkgYieldsC[nBkg];
	for(unsigned bkg = 0; bkg < nBkg; ++bkg){
		bkgYieldsC[bkg] = (TH1D*) bkgYields[bkg]->Clone();
		for(unsigned b = 1; b < bkgYieldsC[bkg]->GetNbinsX() + 1; ++b){
			bkgYieldsC[bkg]->SetBinError(b, 0);
		}
	}
	//set systematic unertainties
	setSystUnc(bkgYieldsC, nBkg, bkgNames);
		promptYield = (TH1D*) bkgYieldsC[0]->Clone();
	for(unsigned bkg = 0; bkg < nBkg; ++bkg){ //nBkg -1 if split NP
		promptYield->Add(bkgYieldsC[bkg]);
	}
	*/
	TH1D* bkgYieldsC[nBkg];
	for(unsigned bkg = 0; bkg < nBkg; ++bkg) bkgYieldsC[bkg] = (TH1D*) bkgYields[bkg]->Clone();
	TH1D* totalBkg = (TH1D*) bkgYieldsC[0]->Clone(); //Total background with statistical uncertainty
	for(unsigned bkg = 1; bkg < nBkg; ++bkg){
		totalBkg->Add(bkgYieldsC[bkg]);
	}
	//Remove stat errors (set error to 0)
	//Add systematic uncertainty
	for(unsigned bkg = 0; bkg < nBkg; ++bkg){
		for(unsigned b = 1; b < bkgYields[bkg]->GetNbinsX() + 1; ++b){
			bkgYieldsC[bkg]->SetBinError(b, 0);
			bkgYieldsC[bkg]->SetBinError(b, bkgErrors[bkg]->GetBinContent(b));
		}
	}
	TH1D* totalBkgE = (TH1D*) bkgYieldsC[0]->Clone();
	for(unsigned bkg = 1; bkg < nBkg; ++bkg){
		totalBkgE->Add(bkgYieldsC[bkg]);
	}
	std::ofstream table;
	table.open("tables/table_" + catName + extraName + ".txt");
	const unsigned nSR = promptYield->GetNbinsX();
	TString srString[nSR];
	for(unsigned sr = 0; sr < nSR; ++sr){
		std::ostringstream strs; //strs << BinWidth; std::string Yaxis = strs.str();
		//strs << std::fixed << std::setprecision(2) << totalBkg->GetBinContent(sr + 1) << " $\\pm$ " << totalBkg->GetBinError(sr + 1) << " stat $\\pm$ " << totalBkgE->GetBinError(sr + 1) << " syst";
		int precYield = 0, precStat = 0, precSyst = 0;
		double x = totalBkg->GetBinContent(sr + 1);
		precYield = precision(x);
		x = totalBkg->GetBinError(sr + 1);
		precStat = precision(x);
		x = totalBkgE->GetBinError(sr + 1);
		precSyst =  precision(x);
		std::cout << "sr = " << sr << std::endl;
		std::cout << "precYield - std::max(precStat, precSyst) = " << precYield - std::max(precStat, precSyst) << std::endl;
		//Calculate observed precision
		int obsPrec = 1;
	 	double tempObs = observed->GetBinContent(sr + 1);
		while(tempObs > 10.){
			++obsPrec;
			tempObs /= 10.;
		}
		//strs << std::setprecision(2 + (precYield - std::max(precStat, precSyst) ) ) << totalBkg->GetBinContent(sr + 1) << " $\\pm$ " <<  std::setprecision( ((precStat < precSyst) ? 2 + (precStat - precSyst): 2 ) ) << totalBkg->GetBinError(sr + 1) << " stat. $\\pm$ " << std::setprecision( ((precSyst < precStat) ? 2 + (precSyst - precStat): 2 ) ) << totalBkgE->GetBinError(sr + 1) << " syst.";
		strs << std::setprecision(2 + (precYield - std::max(precStat, precSyst) ) ) << totalBkg->GetBinContent(sr + 1) << " $\\pm$ " <<  std::setprecision( ((precStat < precSyst) ? 2 + (precStat - precSyst): 2 ) ) << totalBkg->GetBinError(sr + 1) << " $\\pm$ " << std::setprecision( ((precSyst < precStat) ? 2 + (precSyst - precStat): 2 ) ) << totalBkgE->GetBinError(sr + 1) << " pred.   " << std::setprecision(obsPrec) << observed->GetBinContent(sr + 1) << " obs.";
		srString[sr] = strs.str();
	}
	if(catName.Contains("lowM")){
		TString ptBin;
		if(catName == "lowM_3lnoOSSF_highPt") ptBin = "30-55";
		else ptBin = "< 30";
		table << "\\begin{tabular}{|c|c|c|c|c|} \n";
		table << "\\hline \n";
		table << "\\multirow{2}{*}{$\\pt^\\text{leading}$ (GeV)} & \\multicolumn{4}{ c| } {$\\minMOS$  (GeV)} \\\\ \\cline{2-5} \n";
		table << " & $ < 10$ & $10 - 20$ & $20 - 30$ & $> 30$ \\\\ \n";
		table << "\\hline \\hline \n";
		table << "$" << ptBin << "$ & " << srString[0] << " & " << srString[1] << " & " << srString[2] << " & " << srString[3] << "\\\\ \\hline \n";
		table << "\\end{tabular} \n";
	} else if(catName == "highM_3lOSSF"){
		table << "\\begin{tabular}{|c|c|c|c|c|} \n";
		table << "\\hline \n";
		table << "\\multirow{2}{*}{$\\Mtril$ (GeV)} & \\multirow{2}{*}{$\\MT$ (GeV)} & \\multicolumn{3}{ c| } {$\\minMOS$  (GeV)} \\\\ \\cline{3-5} \n";
		table << "  & & $ < 100 $ & $100 - 200 $ & $ > 200$ \\\\ \n";
		table << "\\hline\\hline \n";
		table << "\\multirow{4}{*}{$ 0 - 100$}   & $ < 100$ & \\multicolumn{3}{c|}{" << srString[0] << "} \\\\ \\cline{2-5} \n";
		table << "& $100-200$ & \\multicolumn{3}{c|}{" << srString[1] << "} \\\\ \\cline{2-5} \n";
		table << "& $> 200$ & \\multicolumn{3}{c|}{" << srString[2] << "} \\\\ \\hline \n";
		table << "\\multirow{6}{*}{$> 100$} & $ < 100  $ &" << srString[3] << " & " <<  srString[8] << " & " <<  srString[12] << "\\\\ \\cline{2-5} \n";
		table << "& $100-200 $ & " << srString[4] << " & " << srString[9] << " & " << srString[13] << " \\\\ \\cline{2-5} \n";
		table << "& $200-300 $ & " << srString[5] << " & " << srString[10] << " & " << srString[14] << "\\\\ \\cline{2-5} \n";
		table << "& $300-400 $ & " << srString[6] << " & " << "\\multirow{3}{*}{" << srString[11] << "} & \\multirow{3}{*}{" << srString[15] << "} \\\\ \\cline{2-3} \n";
		table << "& $>400$ & " << srString[7] << " &  & \\\\ \\hline \n";
		table << "\\end{tabular} \n";
	} else if(catName == "highM_3lnoOSSF"){
		table << "\\begin{tabular}{|c|c|c|c|c|} \n";
		table << "\\hline \n";
		table << "\\multirow{2}{*}{$\\Mtril$ (GeV)} & \\multirow{2}{*}{$\\MT$ (GeV)} & \\multicolumn{3}{ c| } {$\\minMOS$  (GeV)} \\\\ \\cline{3-5} \n";
		table << "& & $ < 100 $ & $100 - 200 $ & $ > 200$ \\\\ \n";
		table << "\\hline \\hline \n";
		table << "\\multirow{2}{*}{$ 0 - 100$}  & $ < 100$ & \\multicolumn{3}{c|}{" << srString[0] << "} \\\\ \\cline{2-5} \n";
		table << "& $> 100$ & \\multicolumn{3}{c|}{" << srString[1] << "} \\\\ \\hline \n";
		table << "\\multirow{4}{*}{$> 100$} & $ < 100  $ & " << srString[2] << " & " << srString[6] <<" & \\multirow{4}{*}{" << srString[8] << "} \\\\ \\cline{2-4} \n";
		table << "& $100-150 $ & " << srString[3] << " & " << "\\multirow{3}{*}{" << srString[7] <<"} &  \\\\ \\cline{2-3} \n";
		table << "& $150-250 $ & " << srString[4] << " & & \\\\ \\cline{2-3} \n";
		table << "& $> 250$ & " << srString[5] << " & & \\\\ \\hline \n";
		table << "\\end{tabular} \n";
	}
	table.close();	
}





























