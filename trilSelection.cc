#include "trilTree.h"
#include "MultilepSUSYfunc.h"


bool trilTree::baseline(const bool multilepton, const bool looseBtag, const bool jetCleaning, const bool lowMllVeto){			
	//Baseline event selection:
	if(multilepton && _nL < 3) return false;
	if(!_METfilters_pass) return false;
	//Clean taus
	for(unsigned l = 0; l < _nL; ++l){
		if(_flavors[l] == 2){
			_isloose[l] = _isloose[l] && _lClean[l];
			_isFO[l] = _isloose[l] && _lClean[l];
			_istight[l] = _istight[l] && _lClean[l];
		}
	}
	//preliminary b veto on jets with underestimated JEC
	//if(nBJets(true, false, 2) != 0) return false;
	//Veto events with low mll lepton pairs
	if(lowMllVeto){
		if(!vetoLowMll(12.)) return false;
	}
	return true;
}


//Check if jet is clean
bool trilTree::jetIsClean(unsigned jetI){
	TLorentzVector jet;
	jet.SetPtEtaPhiE(_jetPt[jetI], _jetEta[jetI], _jetPhi[jetI], _jetE[jetI]);
	bool clean = true;
	for(unsigned l = 0; l < _nL; ++l){
		if(_isFO[l]){
			TLorentzVector lepton;
			lepton.SetPtEtaPhiE(_lPt[l], _lEta[l], _lPhi[l], _lE[l]);
			if(jet.DeltaR(lepton) < 0.4){
				return false;
			}
		}
	}
	return true;
}

//Make new cleaned jet collection
void trilTree::cleanJets(){
	unsigned nJ = 0;
	unsigned* jetI = new unsigned[_nJets];
	_HT = 0; 
	//Check which jets are cleaned
	for(unsigned j = 0; j < _nJets; ++j){
		if(jetIsClean(j)){
			jetI[nJ] = j;
			++nJ;	
			if(_jetPt[j] > 30)_HT+= _jetPt[j];
		}
	}
	//make new jet collection
	_nJets = nJ;
	for(unsigned j = 0; j < nJ; ++j){
		_jetPt[j] = _jetPt[jetI[j]];
		_jetEta[j] = _jetEta[jetI[j]];
		_jetPhi[j] = _jetPhi[jetI[j]];
		_jetE[j] = _jetE[jetI[j]];
		_csv[j] = _csv[jetI[j]];
		_jetFlavour[j] = _csv[jetI[j]];
	}
}

//Return number of cleaned jets
unsigned trilTree::nCleanedJets(){
	unsigned nJ = 0;
	//Check which jets are cleaned
	for(unsigned j = 0; j < _nJets; ++j){
		if(jetIsClean(j)){
			++nJ;
		}
	}
	return nJ;
}

void trilTree::selectJets(){
	unsigned nJ = 0;
	unsigned* jetI = new unsigned[_nJets];
	//Check which jets are cleaned
	for(unsigned j = 0; j < _nJets; ++j){
		if(_jetPt[j] > 25){
			jetI[nJ] = j;
			++nJ;	
		}
	}
	//make new jet collection
	_nJets = nJ;
	for(unsigned j = 0; j < nJ; ++j){
		_jetPt[j] = _jetPt[jetI[j]];
		_jetEta[j] = _jetEta[jetI[j]];
		_jetPhi[j] = _jetPhi[jetI[j]];
		_jetE[j] = _jetE[jetI[j]];
		_csv[j] = _csv[jetI[j]];
		_jetFlavour[j] = _csv[jetI[j]];
	}
}


unsigned trilTree::nBJets(const bool looseBtag, const bool jetCleaning, const unsigned jec){
	//Calculate number of bjets 
	unsigned nbJets = 0;
	for(unsigned j = 0; j < _nJets; ++j){
		//Clean jets
		if(jetCleaning && !jetIsClean(j)) continue;
		//Check jet Pt cuts
		if(jec == 0 && _jetPt[j] < 25) continue;
		else if(jec == 1 && _jetPtDown[j] < 25) continue;
		else if(jec == 2 && _jetPtUp[j] < 25) continue;
		if(_csv[j] > (looseBtag ? 0.5426 : 0.8484) ){
			++nbJets;
		}
	}	
	//_n_bJets = nbJets;
	return nbJets;
}


bool trilTree::vetoLowMll(const double mllCut){
	for(unsigned l1 = 0; l1 < _nL -1; ++l1){
		if(_isloose[l1]){
			for(unsigned l2 = l1 + 1; l2 < _nL; ++l2){
				if(_isloose[l2]){
					if(_flavors[l1] == _flavors[l2]){
						if(_charges[l1] != _charges[l2]){
							TLorentzVector lep1, lep2;
							lep1.SetPtEtaPhiE(_lPt[l1], _lEta[l1], _lPhi[l1], _lE[l1]);
							lep2.SetPtEtaPhiE(_lPt[l2], _lEta[l2], _lPhi[l2], _lE[l2]);
							if( (lep1 + lep2).M() < mllCut){
								return false;
							}
						}
					}
				}
			}
		}
	}
	return true;
}



unsigned trilTree::lepOrder(unsigned* ind, const unsigned lepCut, const bool skipTau, const bool cutBased){
	unsigned lCount = 0;	
	for(unsigned l = 0; l < _nL; ++l){
		if(_isFO[l] && (!skipTau || _flavors[l] != 2) ){
			ind[lCount] = l;
			++lCount;
		}
	}
	if(lCount < lepCut) return 0;
	unsigned* ordind = new unsigned[lCount];
	std::set<unsigned> usedLep;
	for(unsigned k =0; k < lCount; ++k){
		//unsigned maxI = 999;
		double maxConePt = 0;
		for(unsigned l = 0; l < lCount; ++l){
			if(usedLep.find(ind[l]) == usedLep.end()){
				double conePt;
				if(!cutBased) conePt = PtCone(_lPt[ind[l]], _flavors[ind[l]], _lepMVA[ind[l]], _ptratio[ind[l]]);
				else conePt = _lPt[ind[l]]*std::max(1., 1 + _isolation[ind[l]] - 0.1);
				if(conePt > maxConePt){
					maxConePt = conePt;
					ordind[k] = ind[l];
				}
			}
		}
		usedLep.insert(ordind[k]);
	}
	for(unsigned l = 0; l < lCount; ++l){
		ind[l] = ordind[l];
	}
	delete[] ordind;
	return lCount;
}

unsigned trilTree::tightCount(const unsigned* ind, const unsigned lCount){
	unsigned tightC = 0;
	for(unsigned l = 0; l < lCount; ++l){
		if(_istight[ind[l]]) ++tightC;
		else return tightC;
	}
	return tightC;
}

bool trilTree::ptCuts(unsigned* ind, unsigned& lCount){
	//MuonEG cut, code this in cleaner way
	if(_flavors[ind[0]] == 1){
		bool muoneg = false;
		if(lCount == 3){
			if(_flavors[ind[1]] == 0 && _flavors[ind[2]] != 2) muoneg = true;
		} else if(lCount > 3){
			if(_flavors[ind[2]] == 0 && _flavors[ind[3]] != 1) muoneg = true;
		}
		if(muoneg){
			if(PtCone(_lPt[ind[0]], _flavors[ind[0]], _lepMVA[ind[0]], _ptratio[ind[0]]) < 25){
				for(unsigned l = 1; l < lCount; ++l){
					ind[l-1] = ind[l];
				}
				lCount -= 1;
			}
		}
	}
	//Apply analysis Pt and eta cuts
	unsigned nTau = 0;
	for(unsigned l = 0; l < lCount; ++l){
		if(_flavors[ind[l]] == 2) ++nTau;
	}
	if(nTau < 2){
		if(_flavors[ind[1]] == 1 && _lPt[ind[1]] < 10) return false;
		else if(_flavors[ind[1]] == 0 && _lPt[ind[1]] < 10) return false;
		//if(_lPt[ind[2]] < 10) continue;
		if(_flavors[ind[0]] == 1 && _flavors[ind[1]] != 1 && _flavors[ind[2]] != 1){
			if(_lPt[ind[0]] < 20) return false; // USED TO BE 25
		} else if(_flavors[ind[0]] == 1){
			if(_lPt[ind[0]] < 20) return false;
		} else{
			if(_lPt[ind[0]] < 25) return false;
		}
		for(unsigned l = 2; l < lCount; ++l){
			if(_lPt[ind[l]] < 10) return false;
		}
	} else{
		for(unsigned l = 0; l < lCount; ++l){
			if(_flavors[ind[l]] == 2 && _lPt[ind[l]] < 20){
				return false;
			} else if (_flavors[ind[l]] == 1 && _lPt[ind[l]] < 25){
				return false;
			} else if (_flavors[ind[l]] == 0 && _lPt[ind[l]] < 30){
				return false;
			}
			if( fabs(_lEta[ind[l]]) > 2.1) return false;
		}
	}
	return true;
}

bool trilTree::promptMC(const unsigned* ind, const unsigned lCount){
	unsigned* match = new unsigned[_nL];
	matchGenLep(match);
	for(unsigned l = 0; l < lCount; ++l){
		if(_flavors[ind[l]] != 2){
			if(match[ind[l]] > 20){
				return false;
			}
			if(!_gen_isPromptl[match[ind[l]]] && fabs(_gen_lmompdg[match[ind[l]]]) != 15){
				return false;
			}
		}
	}
	return true;	
}

bool trilTree::trigger(const unsigned* ind, const unsigned lCount){
	if(lCount > 2){
		bool trigPass[9];
		trigPass[0] = _tril_trigger_eee || _tril_trigger_all;
		trigPass[1] = _tril_trigger_eem || _tril_trigger_all;
		trigPass[2] = _tril_trigger_emm || _tril_trigger_all;
		trigPass[3] = _tril_trigger_mmm || _tril_trigger_all;
		trigPass[4] = _tril_trigger_eee || _tril_trigger_all; //eet and eee events have same trigger
		trigPass[5] = _tril_trigger_emt || _tril_trigger_all;
		trigPass[6] = _tril_trigger_mmm || _tril_trigger_all; //mmt and mmm events have same trigger
		trigPass[7] = _tril_trigger_ett || _tril_trigger_all;
		trigPass[8] = _tril_trigger_mtt || _tril_trigger_all;
		trigPass[9] = false;   //we do not consider 3tau events at the time
		if(!trigPass[tril_flavorComb(ind, _flavors, lCount)]) return false;
	//same-sign dilepton trigger paths to be filled in!!!!
	} else{
	}
	return true;
}

//bool trilTree::ptCuts_2lSS(const unsigned* ind, const unsigned lCount){
//}

bool trilTree::ptCuts_hnl(const unsigned* ind, const unsigned lCount){	
	if(lCount < 3) return false;
	//baseline Pt cuts
	if(_lPt[ind[0]]*(1 + std::max(_isolation[ind[0]] - 0.1, 0.)) < 15) return false;
	if(_lPt[ind[1]]*(1 + std::max(_isolation[ind[1]] - 0.1, 0.)) < 10) return false;
	if(_lPt[ind[2]]*(1 + std::max(_isolation[ind[2]] - 0.1, 0.)) < (10 - 5*_flavors[ind[2]]) ) return false;
	//Additional Pt cuts required by trigger path
	if(tril_flavorComb(ind, _flavors, lCount) == 0){ //eee
		if(_lPt[ind[0]]*(1 + std::max(_isolation[ind[0]] - 0.1, 0.)) < 19) return false;
		if( (_lPt[ind[0]]*(1 + std::max(_isolation[ind[0]] - 0.1, 0.)) < 30) && (_lPt[ind[1]]*(1 + std::max(_isolation[ind[1]] - 0.1, 0.)) < 15) ) return false;
	} else if(tril_flavorComb(ind, _flavors, lCount) == 1){ //eem
		if(_flavors[ind[2]] == 0){
			if(_lPt[ind[2]]*(1 + std::max(_isolation[ind[2]] - 0.1, 0.)) < 15){
				if(_lPt[ind[0]]*(1 + std::max(_isolation[ind[0]] - 0.1, 0.)) < 23) return false; //23
			}
		} else if(_flavors[ind[2]] == 1){
			if(_lPt[ind[2]]*(1 + std::max(_isolation[ind[2]] - 0.1, 0.)) < 8){//8
				if(_lPt[ind[0]]*(1 + std::max(_isolation[ind[0]] - 0.1, 0.)) < 25) return false;
				if((_lPt[ind[0]]*(1 + std::max(_isolation[ind[0]] - 0.1, 0.)) < 30) && (_lPt[ind[1]]*(1 + std::max(_isolation[ind[1]] - 0.1, 0.)) < 15) ) return false;
			} else{
				if( (_lPt[ind[0]]*(1 + std::max(_isolation[ind[0]] - 0.1, 0.)) < 23) && (_lPt[ind[1]]*(1 + std::max(_isolation[ind[1]] - 0.1, 0.)) < 15)) return false;
			}
		}
	} else if(tril_flavorComb(ind, _flavors, lCount) == 2){//emm
		if(_flavors[ind[2]] == 1){
			if(_lPt[ind[2]]*(1 + std::max(_isolation[ind[2]] - 0.1, 0.)) < 9){
				if(_lPt[ind[0]]*(1 + std::max(_isolation[ind[0]] - 0.1, 0.)) < 23) return false; //23
			}
		}
	}
	//No cuts needed for mmm aside from the baseline
	return true;
}

double trilTree::getEventSF(const unsigned* ind, const unsigned lCount, const bool hnl){
	double eventSF = 1.;
	//Apply ID and reco SF to simulation
	for(unsigned l = 0; l < lCount; ++l){  //CHANGE BACK BACK BACK
		if(_istight[ind[l]]){
			if(_flavors[ind[l]] == 2){
				eventSF*=0.83;
			} else if(_flavors[ind[l]] == 0){
				eventSF*=idTightSFMap[0]->GetBinContent(idTightSFMap[0]->FindBin(std::min(_lPt[ind[l]], 199.), fabs(_lEta[ind[l]])));
				if(hnl) eventSF*=baseidSFMap[0]->GetBinContent(baseidSFMap[0]->FindBin(std::min(_lPt[ind[l]], 199.), fabs(_lEta[ind[l]])));
				eventSF*=recSFMap_ele->GetBinContent(recSFMap_ele->FindBin(_lEta[ind[l]]));
			} else if(_flavors[ind[l]] == 1){
				//double mapPt = _lPt[ind[l]] > 10 ? _lPt[ind[l]] : 10.1;
				double max = hnl ? 199. : 119.;
				if(hnl){ //Different binning order of maps
					eventSF*=idTightSFMap[1]->GetBinContent(idTightSFMap[1]->FindBin(fabs(_lEta[ind[l]]), std::min(_lPt[ind[l]], max) ) );
					eventSF*=baseidSFMap[1]->GetBinContent(baseidSFMap[1]->FindBin(fabs(_lEta[ind[l]]), std::min(_lPt[ind[l]], max) ) );
				} else {
					eventSF*=idTightSFMap[1]->GetBinContent(idTightSFMap[1]->FindBin(std::min(_lPt[ind[l]], max), fabs(_lEta[ind[l]])));
					eventSF*=baseidSFMap[1]->GetBinContent(baseidSFMap[1]->FindBin(std::min(_lPt[ind[l]], max), fabs(_lEta[ind[l]])));
				}
				if(_lPt[l] > 10) eventSF*=recSFMap_mu_ptAbove10->Eval(_lEta[ind[l]]);
				else eventSF*=recSFMap_mu_ptBelow10->Eval(_lEta[ind[l]]);
			}
		} else if(_isFO[ind[l]]){
			;
		} else if(_isloose[ind[l]]){
			;
		}
	}
	//Apply PU reweighing
	eventSF *= PUweights[0]->GetBinContent(PUweights[0]->FindBin(std::min(_n_trueInteractions, 49) ));
	//Apply btaf SF
	eventSF *= bTagSF(hnl, 0);
	return eventSF;
}

double trilTree::bTagSF(const bool hnl, const unsigned bTagSyst){//0 = central; 1 = down; 2 = up
	double* bTagSF = new double[_nJets];
	std::string syst;
 	if(bTagSyst == 0) syst = "central";
	else if(bTagSyst == 1) syst = "down";
	else if(bTagSyst == 2) syst = "up";
	BTagEntry::JetFlavor jetFlavor;
	for(unsigned j = 0; j < _nJets; ++j){
		if((bTagSyst == 0 && _jetPt[j] < 25) || (bTagSyst == 1 && _jetPtDown[j] < 25) || (bTagSyst == 2 && _jetPtUp[j] < 25) ){
			bTagSF[j] = 1;
			continue;
		}
		if(_jetFlavour[j] == 0) jetFlavor = BTagEntry::FLAV_UDSG;
		else if(_jetFlavour[j] == 4) jetFlavor = BTagEntry::FLAV_C;
		else if(_jetFlavour[j] == 5) jetFlavor = BTagEntry::FLAV_B;
		else{
			bTagSF[j] = 1;
			continue;
		}	
		bTagSF[j] = reader.eval_auto_bounds(syst, jetFlavor, _jetEta[j], _jetPt[j], _csv[j]);
	}
	return bTagWeight(_jetFlavour, _jetPt, _jetEta, _csv, bTagSF, _nJets, bTagEff, hnl);
}		


bool trilTree::wgSel(unsigned& ph, const unsigned* ind, const unsigned lCount, const double photonPtCut, const double deltaRCut, const double deltaPhiCut, const bool pixVeto, const bool eVeto, const bool onlyMu){
	//Cuts on the lepton
	//if(tightCount(ind, lCount) != 1) return false;
	if(_flavors[ind[0]] == 2) return false;
	if(onlyMu && _flavors[ind[0]] != 1) return false;
	if(_flavors[ind[0]] == 1 && _lPt[ind[0]] < 25) return false;
	if(_flavors[ind[0]] == 0 && _lPt[ind[0]] < 25) return false; //put else
	//Find a good photon
	TLorentzVector lep, phot;
	lep.SetPtEtaPhiE(_lPt[ind[0]], _lEta[ind[0]], _lPhi[ind[0]], _lE[ind[0]]);
	unsigned nPh = 0;
	//if(_nPh != 1) return false;
	double maxPhPt = 0;
	if(_nPh == 0) return false;
	for(unsigned p = 0; p < _nPh; ++p){
		if(_phmva_pass[p] && fabs(_phEta[p]) < 2.5){
			phot.SetPtEtaPhiE(_phPt[p], _phEta[p], _phPhi[p], _phE[p]);
			if(lep.DeltaR(phot) < deltaRCut) continue;
			if(fabs(_phEta[p]) > 1.4442 && fabs(_phEta[p]) < 1.566) continue;
			++nPh;
			if(_phPt[p] > maxPhPt){
				ph = p;
				maxPhPt = _phPt[p];
			}
		}
	}
	//if(nPh != 1) return false;
	for(unsigned ph = 0; ph < _nPh; ++ph){
		std::cout << "Photon " << ph << " has Pt = " << _phPt[ph] << std::endl;	
	}
	if(maxPhPt < photonPtCut) return false;
	std::cout << "Do we get here?" << std::endl;
	if(pixVeto && _phpixelseed_pass[ph]) return false;
	if(eVeto && !_pheveto_pass[ph]) return false;
	//Cuts on the distance between the lepton and photon
	//TLorentzVector lep, phot;
	//lep.SetPtEtaPhiE(_lPt[ind[0]], _lEta[ind[0]], _lPhi[ind[0]], _lE[ind[0]]);
	//phot.SetPtEtaPhiE(_phPt[ph], _phEta[ph], _phPhi[ph], _phE[ph]);
	if(lep.DeltaR(phot) < deltaRCut) return false;
	if(fabs(lep.DeltaPhi(phot)) < deltaPhiCut) return false;
	return true;

}

bool trilTree::wzSel(const unsigned* ind, unsigned* mllI, unsigned& lw, const unsigned lCount, const TLorentzVector* lepV, const bool flavorDiff, const bool onZ, const bool newalgo){
	if(tightCount(ind, lCount) != 3) return false;
	if(SR_EWK_cat(ind, _flavors, _charges, 3) != 0) return false;
	if(newalgo){
		mllIndNew(mllI, ind, lepV, _charges, _flavors, 3, _met, _met_phi);
	} else{
		mllIndices(mllI, ind, lepV, _charges, _flavors, 3);
	}
	lw = 99;
	for(unsigned l = 0; l < 3; ++l){
		if(ind[l] != mllI[0] && ind[l] != mllI[1]){
			lw  = ind[l];
		}
	}
	if(_lPt[lw] < 25) return false;

	if(flavorDiff && _flavors[lw] == _flavors[mllI[0]]) return false;
	//DELETE THIS LINE LATER
	//if(_flavors[lw] != _flavors[mllI[0]]) return false;
	double mll;
	if(mllI[0] == 99){
		mll = -1;
	} else{
		TLorentzVector lzV[2];
		for(unsigned l = 0; l < 2; ++l) lzV[l].SetPtEtaPhiE(_lPt[mllI[l]], _lEta[mllI[l]], _lPhi[mllI[l]], _lE[mllI[l]]);
		mll = (lzV[0] + lzV[1]).M();
	}
	if(onZ && fabs(mll - 91) > 15) return false;
	return true;
}
























