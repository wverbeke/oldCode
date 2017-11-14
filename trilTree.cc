#include "trilTree.h"

#include "../bTag/BTagCalibrationStandalone.h"

trilTree::trilTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
	if (tree != 0){
		Init(tree);
	}
}

Int_t trilTree::GetEntry(Long64_t entry)
{
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

/*
Long64_t fakeTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}
*/
void trilTree::Init(TTree *tree, const bool photonTree, const bool genTree)
{
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("_nL", &_nL, &b__nL);
    fChain->SetBranchAddress("_lPt", _lPt, &b__lPt);
    fChain->SetBranchAddress("_lEta", _lEta, &b__lEta);
    fChain->SetBranchAddress("_lPhi", _lPhi, &b__lPhi);
    fChain->SetBranchAddress("_lE", _lE, &b__lE);
    fChain->SetBranchAddress("_charges", _charges, &b__charges);
    fChain->SetBranchAddress("_flavors", _flavors, &b__flavors);
    fChain->SetBranchAddress("_isloose", _isloose, &b__isloose);
    fChain->SetBranchAddress("_isFO", _isFO, &b__isFO);
	fChain->SetBranchAddress("_istight", _istight, &b__istight);
	fChain->SetBranchAddress("_islooseCut", _islooseCut, &b__islooseCut);
	fChain->SetBranchAddress("_isFOCut", _isFOCut, &b__isFOCut);
	fChain->SetBranchAddress("_istightCut", _istightCut, &b__istightCut);
	fChain->SetBranchAddress("_tightCharge", _tightCharge, &b__tightCharge);
	fChain->SetBranchAddress("_mu_ismedium", _mu_ismedium, &b__mu_ismedium);
	fChain->SetBranchAddress("_hitsNumber", _hitsNumber, &b__hitsNumber);
	fChain->SetBranchAddress("_lClean", _lClean, &b__lClean);
	fChain->SetBranchAddress("_lCleanCut", _lCleanCut, &b__lCleanCut);
	fChain->SetBranchAddress("_met", &_met, &b__met);
	fChain->SetBranchAddress("_met_phi", &_met_phi, &b__met_phi);
	fChain->SetBranchAddress("_nJets", &_nJets, &b__nJets);
	fChain->SetBranchAddress("_jetEta", _jetEta, &b__jetEta);
	fChain->SetBranchAddress("_jetPhi", _jetPhi, &b__jetPhi);
	fChain->SetBranchAddress("_jetPt", _jetPt, &b__jetPt);
	fChain->SetBranchAddress("_jetE", _jetE, &b__jetE);
	fChain->SetBranchAddress("_jetFlavour", _jetFlavour, &b__jetFlavour);
	fChain->SetBranchAddress("_csv", _csv, &b__csv);
	fChain->SetBranchAddress("_closeJetCSVAll", _closeJetCSVAll, &b__closeJetCSVAll);
	fChain->SetBranchAddress("_closeJetPtAll", _closeJetPtAll, &b__closeJetPtAll);
	fChain->SetBranchAddress("_isolation", _isolation, &b__isolation);
	fChain->SetBranchAddress("_n_bJets", &_n_bJets, &b__n_bJets);
	fChain->SetBranchAddress("_eventNb", &_eventNb, &b__eventNb);
	fChain->SetBranchAddress("_runNb", &_runNb, &b__runNb);
	fChain->SetBranchAddress("_lumiBlock", &_lumiBlock, &b__lumiBlock);
	fChain->SetBranchAddress("_lepMVA", _lepMVA, &b__lepMVA);
	fChain->SetBranchAddress("_mvaValue", _mvaValue, &b__mvaValue);
	fChain->SetBranchAddress("_clusterpass", &_clusterpass, &b__clusterpass);
	fChain->SetBranchAddress("_vtxFitConversion", &_vtxFitConversion, &b__vtxFitConversion);
	fChain->SetBranchAddress("_n_PV", &_n_PV, &b__n_PV);
	fChain->SetBranchAddress("_tril_trigger_all", &_tril_trigger_all, &b__tril_trigger_all);
	fChain->SetBranchAddress("_tril_trigger_eee", &_tril_trigger_eee, &b__tril_trigger_eee);
	fChain->SetBranchAddress("_tril_trigger_mmm", &_tril_trigger_mmm, &b__tril_trigger_mmm);
	fChain->SetBranchAddress("_tril_trigger_eem", &_tril_trigger_eem, &b__tril_trigger_eem);
	fChain->SetBranchAddress("_tril_trigger_emm", &_tril_trigger_emm, &b__tril_trigger_emm);
	fChain->SetBranchAddress("_tril_trigger_emt", &_tril_trigger_emt, &b__tril_trigger_emt);
	fChain->SetBranchAddress("_tril_trigger_ett", &_tril_trigger_ett, &b__tril_trigger_ett);
	fChain->SetBranchAddress("_tril_trigger_mtt", &_tril_trigger_mtt, &b__tril_trigger_mtt);
	fChain->SetBranchAddress("_lowM_trigger_all", &_lowM_trigger_all, &b__lowM_trigger_all);
	fChain->SetBranchAddress("_lowM_trigger_mmm", &_lowM_trigger_mmm, &b__lowM_trigger_mmm);
	fChain->SetBranchAddress("_lowM_trigger_eee", &_lowM_trigger_eee, &b__lowM_trigger_eee);
	fChain->SetBranchAddress("_lowM_trigger_mme", &_lowM_trigger_mme, &b__lowM_trigger_mme);
	fChain->SetBranchAddress("_lowM_trigger_mee", &_lowM_trigger_mee, &b__lowM_trigger_mee);
	fChain->SetBranchAddress("_ptrel", _ptrel, &b__ptrel);
	fChain->SetBranchAddress("_miniisolation", _miniisolation, &b__miniisolation);
	fChain->SetBranchAddress("_ptratio", _ptratio, &b__ptratio);
	fChain->SetBranchAddress("_3dIP", _3dIP, &b__3dIP);
	fChain->SetBranchAddress("_ipPV", _ipPV, &b__ipPV);
	fChain->SetBranchAddress("_ipZPV", _ipZPV, &b__ipZPV);
	fChain->SetBranchAddress("_3dIPsig", _3dIPsig, &b__3dIPsig);
	fChain->SetBranchAddress("_METfilters_pass", &_METfilters_pass, &b__METfilters_pass);
	fChain->SetBranchAddress("_weight", &_weight, &b__weight);
	//Add an extra condition so that this isn't read for data files
	fChain->SetBranchAddress("_n_trueInteractions", &_n_trueInteractions, &b__n_trueInteractions);
	//Electron trigger emulation
	/*
	fChain->SetBranchAddress("_dEtaIn", _dEtaIn, &b__dEtaIn);
   	fChain->SetBranchAddress("_dPhiIn", _dPhiIn, &b__dPhiIn);
   	fChain->SetBranchAddress("_ele_hOverE", _ele_hOverE, &b__ele_hOverE);
  	fChain->SetBranchAddress("_full5x5SigmaIetaIeta", _full5x5SigmaIetaIeta, &b__full5x5SigmaIetaIeta);
   	fChain->SetBranchAddress("_ooEmooP", _ooEmooP, &b__ooEmooP);
	*/
	if(photonTree){
		fChain->SetBranchAddress("_nPh", &_nPh, &b__nPh);
		fChain->SetBranchAddress("_phPt", _phPt, &b__phPt);
		fChain->SetBranchAddress("_phEta", _phEta, &b__phEta);
		fChain->SetBranchAddress("_phPhi", _phPhi, &b__phPhi);
		fChain->SetBranchAddress("_phE", _phE, &b__phE);
		fChain->SetBranchAddress("_phmva", _phmva, &b__phmva);
		fChain->SetBranchAddress("_phmva_pass", _phmva_pass, &b__phmva_pass);
		fChain->SetBranchAddress("_pheveto_pass", _pheveto_pass, &b__pheveto_pass);
		fChain->SetBranchAddress("_phpixelseed_pass", _phpixelseed_pass, &b__phpixelseed_pass);
		//triggers used in Wgamma control region
		fChain->SetBranchAddress("_Mu17_Photon22_CaloIdL_L1ISO", &_Mu17_Photon22_CaloIdL_L1ISO, &b__Mu17_Photon22_CaloIdL_L1ISO);
		fChain->SetBranchAddress("_IsoMu24", &_IsoMu24, &b__IsoMu24);
		fChain->SetBranchAddress("_IsoTkMu24", &_IsoTkMu24, &b__IsoTkMu24);
	}
	if(genTree){
		fChain->SetBranchAddress("_gen_nL", &_gen_nL, &b__gen_nL);
		fChain->SetBranchAddress("_gen_lPt", _gen_lPt, &b__gen_lPt);
		fChain->SetBranchAddress("_gen_lEta", _gen_lEta, &b__gen_lEta);
		fChain->SetBranchAddress("_gen_lPhi", _gen_lPhi, &b__gen_lPhi);
		fChain->SetBranchAddress("_gen_lE", _gen_lE, &b__gen_lE);
		fChain->SetBranchAddress("_gen_lmompdg", _gen_lmompdg, &b__gen_lmompdg);
		fChain->SetBranchAddress("_gen_isPromptl", &_gen_isPromptl, &b__gen_isPromptl);
		fChain->SetBranchAddress("_gen_flavors", _gen_flavors, &b__gen_flavors);
		fChain->SetBranchAddress("_gen_charges", _gen_charges, &b__gen_charges);
		fChain->SetBranchAddress("_genmet", &_genmet, &b__genmet);
		fChain->SetBranchAddress("_genmet_phi", &_genmet_phi, &b__genmet_phi);
		//Branches for proper gen matching
		fChain->SetBranchAddress("_pdgmc", _pdgmc, &b__pdgmc);
		fChain->SetBranchAddress("_isPromptFinalState", _isPromptFinalState, &b__isPromptFinalState);
		fChain->SetBranchAddress("_fromHardProcessFinalState", _fromHardProcessFinalState, &b__fromHardProcessFinalState);
		fChain->SetBranchAddress("_origin", _origin, &b__origin);
		fChain->SetBranchAddress("_originDetailed", _originDetailed, &b__originDetailed);
		fChain->SetBranchAddress("_originPhot", _originPhot, &b__originPhot);
		fChain->SetBranchAddress("_mompdg", _mompdg, &b__mompdg);
		//Branches for MC uncertainties
		//MET uncertainties
		fChain->SetBranchAddress("_metJECDown", &_metJECDown, &b__metJECDown);
		fChain->SetBranchAddress("_metJECUp", &_metJECUp, &b__metJECUp);
		fChain->SetBranchAddress("_metOtherDown", &_metOtherDown, &b__metOtherDown);
		fChain->SetBranchAddress("_metOtherUp", &_metOtherUp, &b__metOtherUp);
		fChain->SetBranchAddress("_met_phiJECDown", &_met_phiJECDown, &b__met_phiJECDown);
		fChain->SetBranchAddress("_met_phiJECUp", &_met_phiJECUp, &b__met_phiJECUp);
		fChain->SetBranchAddress("_met_phiOtherDown", &_met_phiOtherDown, &b__met_phiOtherDown);
		fChain->SetBranchAddress("_met_phiOtherUp", &_met_phiOtherUp, &b__met_phiOtherUp);
		//JEC uncertainties
		fChain->SetBranchAddress("_jetPtDown", &_jetPtDown, &b__jetPtDown);
		fChain->SetBranchAddress("_jetPtUp", &_jetPtUp, &b__jetPtUp);
		//Scale and pdf uncertainties
		fChain->SetBranchAddress("_scaleWeight", &_scaleWeight, &b__scaleWeight);
		//lifetime
		fChain->SetBranchAddress("_ctau", &_ctau, &b__ctau);
	}
	if(photonTree && genTree){
		fChain->SetBranchAddress("_gen_nPh", &_gen_nPh, &b__gen_nPh);
		fChain->SetBranchAddress("_gen_phPt", _gen_phPt, &b__gen_phPt);
		fChain->SetBranchAddress("_gen_phEta", _gen_phEta, &b__gen_phEta);
		fChain->SetBranchAddress("_gen_phPhi", _gen_phPhi, &b__gen_phPhi);
		fChain->SetBranchAddress("_gen_phE", _gen_phE, &b__gen_phE);
		fChain->SetBranchAddress("_gen_phmompdg", _gen_phmompdg, &b__gen_phmompdg);
		fChain->SetBranchAddress("_gen_isPromptPh", _gen_isPromptPh, &b__gen_isPromptPh);
	}
	//hcounter[sam] = _hCounter->GetBinContent(1);	
}

void trilTree::readSF(const bool hnl){
	//Read FR maps
	//TFile* light_FRfile = TFile::Open("../weights/fr_comb.root");
	if(hnl){
		TFile* hnl_FRfile[3];
		hnl_FRfile[0] = TFile::Open("../weights/MartinaFRMap_nominal.root");
		hnl_FRfile[1] = TFile::Open("../weights/MartinaFRMap_EWKDown.root");
		hnl_FRfile[2] = TFile::Open("../weights/MartinaFRMap_EWKUp.root");
		for(unsigned ewk = 0; ewk < 3; ++ewk){
			frMap[ewk][0] = (TH2D*) hnl_FRfile[ewk]->Get("frMapEle");
			frMap[ewk][1] = (TH2D*) hnl_FRfile[ewk]->Get("frMapMu");
		}
	} else{	
		TFile* light_FRfile=TFile::Open("../weights/fr_susy_data_v2.2_310117.root");
		TFile* tau_FRfile = TFile::Open("../weights/DY_tauFR.root");
		frMap[0][0] = (TH2D*) light_FRfile->Get("FR_susy_wpM_el_data_comb");
		frMap[0][1] = (TH2D*) light_FRfile->Get("FR_susy_wpM_mu_data_comb");
		frMap[0][2] = (TH2D*) tau_FRfile->Get("taufakesNew_2");
		frMap[1][0] = (TH2D*) light_FRfile->Get("FR_susy_wpM_el_data_comb_down");
		frMap[1][1] = (TH2D*) light_FRfile->Get("FR_susy_wpM_mu_data_comb_down");
		frMap[2][0] = (TH2D*) light_FRfile->Get("FR_susy_wpM_el_data_comb_up");
		frMap[2][1] = (TH2D*) light_FRfile->Get("FR_susy_wpM_mu_data_comb_up");
		//frMap[0] = (TH2D*) light_FRfile->Get("FR_susy_wpV_el_data_comb");
		//frMap[1] = (TH2D*) light_FRfile->Get("FR_susy_wpV_mu_data_comb");
	}

	//Read SF maps
	//Tight ID efficiencies
	/*
	TFile* idSF_file_mu = TFile::Open("../weights/TnP_MuonID_NUM_mvaSUSYM_DENOM_mvaPreSel_VAR_map_pt_eta.root");
	TFile* idSF_file_ele = TFile::Open("../weights/eleScaleFactors.root");
	idTightSFMap[0] = (TH2D*) idSF_file_ele->Get("GsfElectronToLeptonMvaMIDEmuTightIP2DSIP3D8miniIso04");
	idTightSFMap[1] = (TH2D*) idSF_file_mu->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_tag_combRelIsoPF04dBeta_bin0_&_tag_pt_bin0_&_mvaPreSel_pass");
	*/
	TFile* baseidSF_file_mu = TFile::Open("../weights/TnP_NUM_MediumID_DENOM_generalTracks_VAR_map_pt_eta.root");
	TFile* idSF_file_ele = TFile::Open("../weights/scaleFactors_electronMoriond17.root");
	if(hnl){
		//TFile* idSF_file_mu = TFile::Open("../weights/ratio_NUM_RelIsoVTight_DENOM_MediumID_VAR_map_pt_eta_v2.root");
		//baseidSFMap[1] = (TH2D*) baseidSF_file_mu->Get("SF");
		//idTightSFMap[1] = (TH2D*) idSF_file_mu->Get("pt_abseta_ratio"); //Old SUSY SF
		baseidSFMap[0] = (TH2D*) idSF_file_ele->Get("GsfElectronToMVATightTightIP2DSIP3D4");
		idTightSFMap[0] = (TH2D*) idSF_file_ele->Get("MVATightElectronToRelIso010");
		TFile* idSF_file_mu = TFile::Open("../weights/ScaleFactors_ID_ISO.root");
		baseidSFMap[1] = (TH2D*) idSF_file_mu->Get("sf_id");
		idTightSFMap[1] = (TH2D*) idSF_file_mu->Get("sf_iso");
	} else{
		TFile* idSF_file_mu = TFile::Open("../weights/TnP_NUM_mvaSUSYM_DENOM_mvaPreSel_VAR_map_pt_eta.root");
		//baseidSFMap[0] = idSF_file_ele->Get("GsfElectronToMVATightTightIP2DSIP3D4");
		baseidSFMap[1] = (TH2D*) baseidSF_file_mu->Get("SF");
		idTightSFMap[0] = (TH2D*) idSF_file_ele->Get("GsfElectronToLeptonMvaMIDEmuTightIP2DSIP3D8mini04");
		idTightSFMap[1] = (TH2D*) idSF_file_mu->Get("SF");
		/*
		TFile* idSF_file_mu = TFile::Open("../weights/TnP_NUM_mvaSUSYVT_DENOM_mvaPreSel_VAR_map_pt_eta.root");
		idTightSFMap[0] = (TH2D*) idSF_file_ele->Get("GsfElectronToLeptonMvaVTIDEmuTightIP2DSIP3D8mini04");
		idTightSFMap[1] = (TH2D*) idSF_file_mu->Get("SF");
		*/
	}

	/*
	//Loose ID efficiencies
	TFile* idSF_file_muLoose = TFile::Open("../weights/TnP_MuonID_NUM_MediumID_DENOM_generalTracks_VAR_map_pt_eta.root");
	idLooseSFMap[0] = (TH2D*) idSF_file_ele->Get("GsfElectronToLoose2D");
	idLooseSFMap[1] = (TH2D*) idSF_file_muLoose->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0");
	*/
	//Tracking efficiencies
	TFile* recSF_file_mu_ptAbove10 = TFile::Open("../weights/muonTrackingEff_above10.root");
	TFile* recSF_file_mu_ptBelow10 = TFile::Open("../weights/muonTrackingEff_below10.root");
	//TFile* recSF_file_ele = TFile::Open("../weights/egammaEffi.txt_SF2D.root");
	TFile* recSF_file_ele = TFile::Open("../weights/egammaEffi.txt_EGM2D.root");
	recSFMap_ele = (TH1D*) ((TH2D*) recSF_file_ele->Get("EGamma_SF2D"))->ProjectionX();
	recSFMap_mu_ptAbove10 = (TGraph*) recSF_file_mu_ptAbove10->Get("ratio_eta");
	recSFMap_mu_ptBelow10 = (TGraph*) recSF_file_mu_ptBelow10->Get("ratio_eta");
	//Trigger SF
	/*
	TFile* triggerEffFile = TFile::Open("../weights/triggerSF.root");
	triggerEffMap3l[0] = (TH2D*) triggerEffFile->Get("eff_3l_ele");
	triggerEffMap3l[1] = (TH2D*) triggerEffFile->Get("eff_3l_mu");
	triggerEffMap2l[0] = (TH2D*) triggerEffFile->Get("eff_2l_ele");
	triggerEffMap2l[1] = (TH2D*) triggerEffFile->Get("eff_2l_mu");
	//Electron leg  for ett events
	TFile* triggerEleLegFile = TFile::Open("../weights/triggerSF_Ele27_EWKino_fullsim_ICHEP2016_12p9fb.root");
	triggerEffEleLeg = (TH2D*) triggerEleLegFile->Get("hist2dnum_Ele27_WPLoose_Gsf__HLT_Ele27_WPLoose_Gsf");
	*/
	//Set up btag SF reader
	calib = *(new BTagCalibration("csvv2", "../bTag/CSVv2_Moriond17_B_H.csv"));
	if(hnl) reader = *(new BTagCalibrationReader(BTagEntry::OP_LOOSE, "central", {"up", "down"}));
	else reader =  *(new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up", "down"}));
	reader.load(calib, BTagEntry::FLAV_B, "comb");
	reader.load(calib, BTagEntry::FLAV_C, "comb");
	reader.load(calib, BTagEntry::FLAV_UDSG, "incl");
	//Read btag efficiencies
	TFile* bTagEffFile;
	if(hnl) bTagEffFile =  TFile::Open("../weights/bTagEff_hnl.root");
	else bTagEffFile =  TFile::Open("../weights/bTagEff_ewkino.root");
	bTagEff[0] = (TH2D*) bTagEffFile->Get("btagEff_udsg");
	bTagEff[1] = (TH2D*) bTagEffFile->Get("btagEff_charm");
	bTagEff[2] = (TH2D*) bTagEffFile->Get("btagEff_beauty");
	//PU reweighing
	TFile* PUfile = TFile::Open("../weights/puw_nTrueInt_Moriond2017_36p5fb_Summer16_central.root");
	PUweights[0] = (TH1D*) PUfile->Get("puw");
	PUfile = TFile::Open("../weights/puw_nTrueInt_Moriond2017_36p5fb_Summer16_down.root");
	PUweights[1] = (TH1D*) PUfile->Get("puw");
	PUfile = TFile::Open("../weights/puw_nTrueInt_Moriond2017_36p5fb_Summer16_up.root");
	PUweights[2] = (TH1D*) PUfile->Get("puw");
}

void trilTree::matchGenLep(unsigned* match){
	if(_gen_nL == 0 || _nL == 0) return;
	TLorentzVector* Lep = new TLorentzVector[_nL];
    TLorentzVector* GenLep = new TLorentzVector[_gen_nL];
    for(unsigned l = 0; l < _nL; ++l){
    	Lep[l].SetPtEtaPhiE(_lPt[l], _lEta[l], _lPhi[l], _lE[l]);
    }
    for(unsigned l = 0; l < _gen_nL; ++l){
    	GenLep[l].SetPtEtaPhiE(_gen_lPt[l], _gen_lEta[l], _gen_lPhi[l], _gen_lE[l]);
    }
    for(unsigned l = 0; l < _nL; ++l){
    	match[l] = 9999; 
    }
	std::set<unsigned> usedGen;
	std::set<unsigned> usedReco;
    for(unsigned n = 0; n < _nL; ++n){
		double minDeltaR = 9999.;
		unsigned matched = 9999;
    	for(unsigned l = 0; l < _nL; ++l){
        	if(_lPt[l] == 0) continue;
			if(usedReco.find(l) == usedReco.end()){
				for(unsigned g = 0; g < _gen_nL; ++g){
            		if(_gen_lPt[g] == 0) continue;	
					if(std::isnan(_gen_lPt[g]) || std::isnan(_gen_lEta[g]) || std::isnan(_gen_lPhi[g]) || std::isnan(_gen_lE[g])) continue;
					if(std::isinf(_gen_lPt[g]) || std::isinf(_gen_lEta[g]) || std::isinf(_gen_lPhi[g]) || std::isinf(_gen_lE[g])) continue;
            		//if(_charges[l] != _gen_charges[g] || _flavors[l] != _gen_flavors[g]) continue;
					if(_flavors[l] != _gen_flavors[g]) continue;
					if(usedGen.find(g) == usedGen.end()){
						if(Lep[l].DeltaR(GenLep[g]) < minDeltaR){
							matched = l;
							match[l] = g;
							minDeltaR = Lep[l].DeltaR(GenLep[g]);
						}
					}
				}
			}
		}	
		if(matched != 9999){
			usedReco.insert(matched);
			usedGen.insert(match[matched]);
		}
	}
}

void trilTree::cutBased(){
	//const double MVA_cuts_pt15[3] = {-0.86, -0.85, -0.81};
	//const double MVA_cuts_pt25[3] = {-0.96, -0.96, -0.95};
	const double MVA_cuts_tuned[3] = {-0.02, -0.52, -0.52};
	//const double MVA_cuts_tuned[3] = {0.48, 0.10, -0.36};
	//const double MVA_cuts_tuned[3] = {0.5, 0.1, -0.02};
	for(unsigned l = 0; l < _nL; ++l){
		if(_flavors[l] != 2){
			_isloose[l] = _islooseCut[l];
			if(_flavors[l] == 1)_isFO[l] = _islooseCut[l] && _mu_ismedium[l] && _3dIPsig[l] < 4;
			else if(_flavors[l] == 0){
				 _isFO[l] = _islooseCut[l] && _hitsNumber[l] == 0 && _3dIPsig[l] < 4;
				 int eta = -1;
                 if(fabs(_lEta[l]) < 0.8 ) {
                 	eta = 0;
                 } else  if(fabs(_lEta[l]) < 1.479 ) {
                 	eta = 1;
                 } else{
                 	eta = 2;
                 }
				bool passedMVA = _mvaValue[l] > MVA_cuts_tuned[eta];
                // bool passedMVA = _mvaValue[l] >  std::min( MVA_cuts_pt15[eta], std::max(MVA_cuts_pt25[eta] , MVA_cuts_pt15[eta] + (MVA_cuts_pt25[eta] - MVA_cuts_pt15[eta])*0.1 *(_lPt[l]-15) ) );
			 	_isFO[l] = _isFO[l] && passedMVA && _clusterpass[l];
				_istightCut[l] = _istightCut[l] && _clusterpass[l];
				//_istight[l] = _isFO[l] && (_isolation[l] < 0.1);
				//_istightCut[l] = _istightCut[l] && passedMVA;
			}
			//_istight[l] = _isFO[l] && _isolation[l] < 0.1;
			_istight[l] = _istightCut[l] && _3dIPsig[l] < 4;
 			//_istight[l] = _isFO[l] && (_isolation[l] < 0.1);
		}
		if(_flavors[l] == 2){
			_isFO[l] == false;
			_istight[l] == false;
		}
	}
}
