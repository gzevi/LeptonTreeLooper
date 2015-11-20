#include "LeptonTreeLooper.h"
#include "Tools/goodrun.h"
#include "Tools/goodrun.cc"

int LeptonTreeLooper( TChain* chain, TString output_name , int nEvents );
void makePlots(std::map<std::string, TH1*> & h_1, TString sel, float weight);
void makeDilepPlots(std::map<std::string, TH1*> & h_1, TString sel, float weight);

int LeptonTreeLooper( TChain* chain, TString output_name , int nEvents ) {

  //bool doPUweight = false;
  bool fast = true;
  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");


  // Loop over events to Analyze
  unsigned int nEventsTotal = 0;
  unsigned int nEventsChain = chain->GetEntries();
  if( nEvents >= 0 ) nEventsChain = nEvents;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  std::map<std::string, TH1*> h_1d;

  //global var for checking if event already has dilep mass stored
  int lastEventSaved_ = -1;
  
  //set lumi in fb
  const float lumi = .010;

  //set goodrun file
  //set_goodrun_file("goodRunList/private_json_and_DCS_150716_snt.txt");
  set_goodrun_file("goodRunList/shilpi.txt");


  //vector for trigger names
  //need same number of entries and same order as trigDecisions!!!
  std::vector<TString> trigNames;
  trigNames.push_back("global");
  trigNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  trigNames.push_back("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300");
  trigNames.push_back("HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p5PF");
  trigNames.push_back("HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30");
  trigNames.push_back("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30");
  trigNames.push_back("HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30");
  trigNames.push_back("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30");
  trigNames.push_back("HLT_Ele33_CaloIdM_TrackIdM_PFJet30");
  trigNames.push_back("HLT_Ele23_CaloIdM_TrackIdM_PFJet30");
  trigNames.push_back("HLT_Ele18_CaloIdM_TrackIdM_PFJet30");
  trigNames.push_back("HLT_Ele12_CaloIdM_TrackIdM_PFJet30");
  trigNames.push_back("HLT_Ele8_CaloIdM_TrackIdM_PFJet30");
  trigNames.push_back("tag_HLT_Ele25WP60_Ele8_Mass55_LeadingLeg");
  trigNames.push_back("tag_HLT_Ele25WP60_SC4_Mass55_LeadingLeg");
  trigNames.push_back("HLT_Ele32_eta2p1_WPTight_Gsf");
  trigNames.push_back("HLT_Ele32_eta2p1_WPLoose_Gsf");
  trigNames.push_back("HLT_Ele27_eta2p1_WPTight_Gsf");
  trigNames.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf");
  // trigNames.push_back("HLT_Ele22_eta2p1_WPLoose_Gsf");
  // trigNames.push_back("HLT_Ele22_eta2p1_WPTight_Gsf");
  trigNames.push_back("HLT_Ele23_WPLoose_Gsf");

  //load PU-ratio histogram for reweighting
  TFile * f_pu = new TFile("puWeight.root","READ");
  TH1D * h_ratio = (TH1D*) f_pu->Get("h_dataOverMC_nvtxEB");

  // File Loop
  while ( (currentFile = (TFile*)fileIter.Next()) ) {

    // Get File Content
    TFile *file = new TFile( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("t");
    if(fast) TTreeCache::SetLearnEntries(10);
    if(fast) tree->SetCacheSize(128*1024*1024);
    lepton_tree_obj.Init(tree);
    
    // Loop over Events in current file
    if( nEventsTotal >= nEventsChain ) continue;
    unsigned int nEventsTree = tree->GetEntriesFast();
    for( unsigned int event = 0; event < nEventsTree; ++event) {
    
      // Get Event Content
      if( nEventsTotal >= nEventsChain ) continue;
      if(fast) tree->LoadTree(event);
      lepton_tree_obj.GetEntry(event);
      ++nEventsTotal;
    
      // Progress
      LeptonTree::progress( nEventsTotal, nEventsChain );

      // Analysis Code

      //good run list
      //if (evt_isRealData() && !goodrun(evt_run(), evt_lumiBlock())) continue; 
      //poor mans goodrunlist
      if (evt_isRealData() && evt_run()!=254231 && evt_run()!=254232 && evt_run()!=254790 && evt_run()!=254852 && evt_run()!=254879) continue; //25ns
      //if (evt_isRealData() && evt_run()!=254833) continue; //50ns
      
      if ( fabs(id()) != 11) continue;
      LorentzVector el_p4 = p4();
      float pt = el_p4.pt();
      float eta =  el_p4.eta();
      float phi =  el_p4.phi();

      //easy variables
      plot1D("h_pt", pt,  1, h_1d, "pT [GeV]", 150, 0, 150);      
      plot1D("h_eta", eta,  1, h_1d, "#eta", 150, -3, 3);
      plot1D("h_phi", phi,  1, h_1d, "#phi", 150, -3.5, 3.5);

      //vector for trigger decisions
      //need same number of entries and same order as trigNames!!!
      std::vector<int> trigDecision;      std::vector<int> tagLL; std::vector<float> trigPtPlat;
      trigDecision.push_back(1);       tagLL.push_back(0);
      trigDecision.push_back(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ());                             tagLL.push_back(/*tag_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_ElectronLeg()*/ 0);                trigPtPlat.push_back(0);
      trigDecision.push_back(HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300());                         tagLL.push_back(/*tag_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_ElectronLeg()*/ 0);            trigPtPlat.push_back(0);
      trigDecision.push_back(HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p5PF());                tagLL.push_back(/*tag_HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p5PF_ElectronLeg()*/ 0);   trigPtPlat.push_back(0);
      trigDecision.push_back(HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30());                              tagLL.push_back(tag_HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_ElectronLeg());                       trigPtPlat.push_back(38);
      trigDecision.push_back(HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30());                              tagLL.push_back(tag_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_ElectronLeg());                       trigPtPlat.push_back(30);
      trigDecision.push_back(HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30());                              tagLL.push_back(tag_HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30_ElectronLeg());                       trigPtPlat.push_back(0);
      trigDecision.push_back(HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30());                              tagLL.push_back(tag_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_ElectronLeg());                       trigPtPlat.push_back(0);
      trigDecision.push_back(HLT_Ele33_CaloIdM_TrackIdM_PFJet30());                                    tagLL.push_back(tag_HLT_Ele33_CaloIdM_TrackIdM_PFJet30_ElectronLeg());                             trigPtPlat.push_back(38);
      trigDecision.push_back(HLT_Ele23_CaloIdM_TrackIdM_PFJet30());                                    tagLL.push_back(tag_HLT_Ele23_CaloIdM_TrackIdM_PFJet30_ElectronLeg());                             trigPtPlat.push_back(28);
      trigDecision.push_back(HLT_Ele18_CaloIdM_TrackIdM_PFJet30());                                    tagLL.push_back(tag_HLT_Ele18_CaloIdM_TrackIdM_PFJet30_ElectronLeg());                             trigPtPlat.push_back(0);
      trigDecision.push_back(HLT_Ele12_CaloIdM_TrackIdM_PFJet30());                                    tagLL.push_back(tag_HLT_Ele12_CaloIdM_TrackIdM_PFJet30_ElectronLeg());                             trigPtPlat.push_back(0);
      trigDecision.push_back(HLT_Ele8_CaloIdM_TrackIdM_PFJet30());                                     tagLL.push_back(tag_HLT_Ele8_CaloIdM_TrackIdM_PFJet30_ElectronLeg());                              trigPtPlat.push_back(0);
      trigDecision.push_back(tag_HLT_Ele25WP60_Ele8_Mass55_LeadingLeg());                              tagLL.push_back(tag_HLT_Ele25WP60_Ele8_Mass55_LeadingLeg());                                       trigPtPlat.push_back(0);
      trigDecision.push_back(tag_HLT_Ele25WP60_SC4_Mass55_LeadingLeg());                               tagLL.push_back(tag_HLT_Ele25WP60_SC4_Mass55_LeadingLeg());                                        trigPtPlat.push_back(0);
      trigDecision.push_back(HLT_Ele32_eta2p1_WPTight_Gsf());                                          tagLL.push_back(tag_HLT_Ele32_eta2p1_WPTight_Gsf());                                               trigPtPlat.push_back(46);
      trigDecision.push_back(HLT_Ele32_eta2p1_WPLoose_Gsf());                                          tagLL.push_back(tag_HLT_Ele32_eta2p1_WPLoose_Gsf());                                               trigPtPlat.push_back(46);
      trigDecision.push_back(HLT_Ele27_eta2p1_WPTight_Gsf());                                          tagLL.push_back(tag_HLT_Ele27_eta2p1_WPTight_Gsf());                                               trigPtPlat.push_back(40);
      trigDecision.push_back(HLT_Ele27_eta2p1_WPLoose_Gsf());                                          tagLL.push_back(tag_HLT_Ele27_eta2p1_WPLoose_Gsf());                                               trigPtPlat.push_back(40);
      // trigDecision.push_back(HLT_Ele22_eta2p1_WPLoose_Gsf());                                       tagLL.push_back(/* HLT_Ele22_eta2p1_WPLoose_Gsf()*/0);                                             trigPtPlat.push_back(0);
      // trigDecision.push_back(HLT_Ele22_eta2p1_WPTight_Gsf());                                       tagLL.push_back(/*HLT_Ele22_eta2p1_WPTight_Gsf()*/0);                                              trigPtPlat.push_back(0);
      trigDecision.push_back(HLT_Ele23_WPLoose_Gsf());                                                 tagLL.push_back(/*HLT_Ele23_WPLoose_Gsf()*/0);                                                     trigPtPlat.push_back(25);


      //make sure all trigger vectors are same size
      if ( trigNames.size() != trigDecision.size() && trigNames.size() != tagLL.size() ){ cout << "Trigger vectors are not equal size! Continuing." << endl; continue; }

      //if real data, set lumiScale to 1
      //float lumiScale = scale1fb() * lumi;
      float lumiScale = scale1fb() * lumi;
      if (evt_isRealData()) lumiScale = 1;      

      //if not data, calculate weight based on PUreweighting
      if (!evt_isRealData()){
	const int vtx_bin = h_ratio->GetXaxis()->FindBin(nvtx());
	float puWeight = h_ratio->GetBinContent(vtx_bin);
	//lumiScale *= puWeight;
      }
      
      if (evt_isRealData()) {
	for ( unsigned int trigIdx=0; trigIdx < trigNames.size(); trigIdx++) {
	  if ( trigDecision[trigIdx] != 0) { makePlots( h_1d, trigNames[trigIdx]+"_noTagReq", lumiScale); }
	}//trig loop
      }

      //skip leptons where a (different) tag doesn't exist in event
      const float mll = dilep_mass();
      if ( mll == -1 ) continue;

      //tag requirements
      if (tag_p4().pt() < 30 ) continue;
      if (fabs(tag_p4().eta()) > 2.1) continue;

      LorentzVector el_seedCl = el_p4;
      el_seedCl *= eSeed()/el_p4.E();
      LorentzVector el_SC = el_p4;
      LorentzVector el_SCraw = el_p4;
      el_SCraw *= eSCRaw() / el_p4.E();

      LorentzVector t_seedCl = tag_p4();
      t_seedCl *= tag_eSeed()/tag_p4().E();
      LorentzVector t_SC = tag_p4();
      LorentzVector t_SCraw = tag_p4();
      t_SCraw *= tag_eSCraw() / tag_p4().E();
      
      LorentzVector dilep_seedCL_p4 = t_seedCl + el_seedCl;
      LorentzVector dilep_SCraw_p4 = t_SCraw + el_SCraw;
      
      float mll_seedCL = dilep_seedCL_p4.M();
      float mll_SCraw = dilep_SCraw_p4.M();
      

      bool bothEB = false;
      bothEB = (fabs(tag_p4().eta()) < 1.479) && (fabs(eta) < 1.479);
      
      if((evt_isRealData() && tag_HLT_Ele27_eta2p1_WPLoose_Gsf() > 0) ||  (!evt_isRealData() ) ){
	if (mll < 100 && mll > 80){ makePlots( h_1d, "Zprobe", lumiScale); }
	if (mll < 100 && mll > 80 && mva()>0){ makePlots( h_1d, "ZprobePassMVA", lumiScale); }
	if (mll < 100 && mll > 80 && (pfChargedHadronIso() +pfPhotonIso() +pfNeutralHadronIso())/pt < 0.1){ makePlots( h_1d, "ZprobePassIso", lumiScale); }
	if (mll < 100 && mll > 80 && passes_HAD_veto_noiso_v3()){ makePlots( h_1d, "ZprobePassID", lumiScale); }
	if ( fabs(eta) > 1.77 ) {
	  if (mll < 100 && mll > 80){ makePlots( h_1d, "ZprobeEta177", lumiScale); }
	  if (mll < 100 && mll > 80 && mva()>0){ makePlots( h_1d, "ZprobePassMVAEta177", lumiScale); }
	  if (mll < 100 && mll > 80 && (pfChargedHadronIso() +pfPhotonIso() +pfNeutralHadronIso())/pt < 0.1){ makePlots( h_1d, "ZprobePassIsoEta177", lumiScale); }
	  if (mll < 100 && mll > 80 && passes_HAD_veto_noiso_v3()){ makePlots( h_1d, "ZprobePassIDEta177", lumiScale); }
	}
	//if (mll < 100 && mll > 80 && conv_vtx_flag() == 0 ){ makePlots( h_1d, "ZprobePixelVeto", lumiScale); }
	//pt 10t20
	if (mll < 100 && mll > 80 && pt > 10 && pt < 20 ){ makePlots( h_1d, "ZprobePt10t20", lumiScale); }
	if (mll < 100 && mll > 80 && mva()>0 && pt > 10 && pt < 20 ){ makePlots( h_1d, "ZprobePassMVAPt10t20", lumiScale); }
	if (mll < 100 && mll > 80 && (pfChargedHadronIso() +pfPhotonIso() +pfNeutralHadronIso())/pt < 0.1 && pt > 10 && pt < 20 ){ makePlots( h_1d, "ZprobePassIsoPt10t20", lumiScale); }
	if (mll < 100 && mll > 80 && passes_HAD_veto_noiso_v3() && pt > 10 && pt < 20 ){ makePlots( h_1d, "ZprobePassIDPt10t20", lumiScale); }
	//pt 20t30
	if (mll < 100 && mll > 80 && pt > 20 && pt < 30 ){ makePlots( h_1d, "ZprobePt20t30", lumiScale); }
	if (mll < 100 && mll > 80 && mva()>0 && pt > 20 && pt < 30 ){ makePlots( h_1d, "ZprobePassMVAPt20t30", lumiScale); }
	if (mll < 100 && mll > 80 && (pfChargedHadronIso() +pfPhotonIso() +pfNeutralHadronIso())/pt < 0.1 && pt > 20 && pt < 30 ){ makePlots( h_1d, "ZprobePassIsoPt20t30", lumiScale); }
	if (mll < 100 && mll > 80 && passes_HAD_veto_noiso_v3() && pt > 20 && pt < 30 ){ makePlots( h_1d, "ZprobePassIDPt20t30", lumiScale); }
	//pt 30t40
	if (mll < 100 && mll > 80 && pt > 30 && pt < 40 ){ makePlots( h_1d, "ZprobePt30t40", lumiScale); }
	if (mll < 100 && mll > 80 && mva()>0 && pt > 30 && pt < 40 ){ makePlots( h_1d, "ZprobePassMVAPt30t40", lumiScale); }
	if (mll < 100 && mll > 80 && (pfChargedHadronIso() +pfPhotonIso() +pfNeutralHadronIso())/pt < 0.1 && pt > 30 && pt < 40 ){ makePlots( h_1d, "ZprobePassIsoPt30t40", lumiScale); }
	if (mll < 100 && mll > 80 && passes_HAD_veto_noiso_v3() && pt > 30 && pt < 40 ){ makePlots( h_1d, "ZprobePassIDPt30t40", lumiScale); }
	//pt 40tInf
	if (mll < 100 && mll > 80 && pt > 40 ){ makePlots( h_1d, "ZprobePt40tInf", lumiScale); }
	if (mll < 100 && mll > 80 && mva()>0 && pt > 40 ){ makePlots( h_1d, "ZprobePassMVAPt40tInf", lumiScale); }
	if (mll < 100 && mll > 80 && (pfChargedHadronIso() +pfPhotonIso() +pfNeutralHadronIso())/pt < 0.1 && pt > 40 ){ makePlots( h_1d, "ZprobePassIsoPt40tInf", lumiScale); }
	if (mll < 100 && mll > 80 && passes_HAD_veto_noiso_v3() && pt > 40 ){ makePlots( h_1d, "ZprobePassIDPt40tInf", lumiScale); }
	//fake
	if (mll < 70 && mll > 30){ makePlots( h_1d, "Fake", lumiScale);}	
      }// isRealData?
      
      //things to only fill once per event
      const int evt = evt_event();
      if( lastEventSaved_ != evt ){
	
	if((evt_isRealData() && tag_HLT_Ele27_eta2p1_WPLoose_Gsf() > 0) ||  (!evt_isRealData() ) ){
	  
	  //inclusive mll
	  plot1D("hZprobe_mll_all",mll, lumiScale, h_1d, "mll",150,0,150);
	  if (bothEB){	plot1D("hZprobe_mllEB_all",mll, lumiScale, h_1d, "mll",150,0,150); }
	  else { plot1D("hZprobe_mllEE_all",mll, lumiScale, h_1d, "mll",150,0,150); }
	  
	  //for failing mva
	  if ( mva() < -0.4) { plot1D("hZprobe_mll_mvaFail_all",mll, lumiScale, h_1d, "mll",150,0,150); }
	  
	  //for passing pogMed id
	  if ( passes_POG_mediumID() ) {
	    plot1D("hZprobe_mll_pogMedium_all",mll, lumiScale, h_1d, "mll",150,0,150);
	    if (bothEB){
	      plot1D("hZprobe_mllEB_pogMedium_all",mll, lumiScale, h_1d, "mll",150,0,150);
	      plot1D("hZprobe_mllEB_SCraw_pogMedium_all",mll_SCraw, lumiScale, h_1d, "mll",150,0,150);
	      plot1D("hZprobe_mllEB_seedCL_pogMedium_all",mll_seedCL, lumiScale, h_1d, "mll",150,0,150); 
	    }
	    else {
	      plot1D("hZprobe_mllEE_pogMedium_all",mll, lumiScale, h_1d, "mll",150,0,150);
	      plot1D("hZprobe_mllEE_SCraw_pogMedium_all",mll_SCraw, lumiScale, h_1d, "mll",150,0,150);
	      plot1D("hZprobe_mllEE_seedCL_pogMedium_all",mll_seedCL, lumiScale, h_1d, "mll",150,0,150); 
	    }	  
	    
	  }//pass med ID
	  
	}// isRealData?

	//------trigger turn ons-------
	
	//denominator for trigger eff
	//require tag matched to HLT_Ele27_eta2p1_WPLoose_Gsf AND passes_POG_mediumID on probe
	if ( (tag_HLT_Ele27_eta2p1_WPLoose_Gsf() > 0 || !evt_isRealData()) && passes_POG_mediumID() ) {
	  makeDilepPlots( h_1d, "tagIsLead_probeMedPog", lumiScale);
	  if ( fabs(eta)<2 ) makeDilepPlots( h_1d, "tagIsLead_probeMedPog_etaLess2", lumiScale) ;
	}

	//make plots for individual triggers
	for ( unsigned int trigIdx=0; trigIdx < trigNames.size(); trigIdx++) {
	  if ( trigDecision[trigIdx] != 0) {
	    
	    //require ptPlat in denominator, for eta efficiency
	    if ( (tag_HLT_Ele27_eta2p1_WPLoose_Gsf() > 0 || !evt_isRealData()) && passes_POG_mediumID() && pt>trigPtPlat[trigIdx] ) makeDilepPlots( h_1d, "tagIsLead_probeMedPog_"+ trigNames[trigIdx] +"_ptPlat", lumiScale) ;

	    //EE + EB plots
	    makeDilepPlots( h_1d, trigNames[trigIdx], lumiScale);
	    
	    //EB or EE plots
	    if ( abs(tag_p4().eta())<1.44 && abs(eta) < 1.44 ) { makeDilepPlots( h_1d, trigNames[trigIdx]+"_EB", lumiScale); }
	    else { makeDilepPlots( h_1d, trigNames[trigIdx]+"_EE", lumiScale); }

	    //require tag matched to HLT_Ele27_eta2p1_WPLoose_Gsf AND passes_POG_mediumID on probe AND probe triggered trigger
	    if ( tag_HLT_Ele27_eta2p1_WPLoose_Gsf() > 0 && passes_POG_mediumID() && trigDecision[trigIdx] > 0 ) {
	      makeDilepPlots( h_1d, trigNames[trigIdx]+"_tagIsLead_probeMedPog_probeTrig", lumiScale);
	      if ( fabs(eta)<2 ) makeDilepPlots( h_1d, trigNames[trigIdx]+"tagIsLead_probeMedPog_probeTrig_etaLess2", lumiScale) ;
	      if ( pt>trigPtPlat[trigIdx] ) makeDilepPlots( h_1d, trigNames[trigIdx]+"tagIsLead_probeMedPog_ptPlat", lumiScale) ;
	    }

	    if ( tag_HLT_Ele27_eta2p1_WPLoose_Gsf() > 0 && passes_POG_mediumID() && probe_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_LeadingLeg() > 0 ) {
	      makeDilepPlots( h_1d, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_LeadingLeg_tagIsLead_probeMedPog_probeTrig", lumiScale);
	      if ( fabs(eta)<2 ) makeDilepPlots( h_1d, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_LeadingLeg_tagIsLead_probeMedPog_probeTrig_etaLess2", lumiScale) ;
	      if ( pt>trigPtPlat[trigIdx] ) makeDilepPlots( h_1d, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_LeadingLeg_tagIsLead_probeMedPog_ptPlat", lumiScale) ;
	    }
	    if ( tag_HLT_Ele27_eta2p1_WPLoose_Gsf() > 0 && passes_POG_mediumID() && probe_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_TrailingLeg() > 0 ) {
	      makeDilepPlots( h_1d, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_TrailingLeg_tagIsLead_probeMedPog_probeTrig", lumiScale);
	      if ( fabs(eta)<2 ) makeDilepPlots( h_1d, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_TrailingLeg_tagIsLead_probeMedPog_probeTrig_etaLess2", lumiScale) ;
	      if ( pt>trigPtPlat[trigIdx] ) makeDilepPlots( h_1d, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_TrailingLeg_tagIsLead_probeMedPog_ptPlat", lumiScale) ;
	    }
	    plot1D("h_trigs",trigIdx-0.5, lumiScale, h_1d, "trigger",trigNames.size(),-0.5,trigNames.size()-0.5); //counts how often each trigger fired

	  }//trigDecision
	}//trig loop
	lastEventSaved_ = evt;
      }//lastEvtSaved
      
    } // end of event loop
    
    // Clean Up
    delete tree;
    file->Close();
    delete file;
  } // end of file loop
  
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }
  
  cout<<"\nWriting file"<<endl;
  savePlots(h_1d, output_name+".root");

  
  // return
  bmark->Stop("benchmark");
  cout << endl;
  cout << nEventsTotal << " Events Processed" << endl;
  cout << "------------------------------" << endl;
  cout << "CPU  Time:	" << Form( "%.01f", bmark->GetCpuTime("benchmark")  ) << endl;
  cout << "Real Time:	" << Form( "%.01f", bmark->GetRealTime("benchmark") ) << endl;
  cout << endl;
  delete bmark;
  return 0;
}


void makePlots(std::map<std::string, TH1*> & h_1, TString sel, float weight = 1) {

  const float pt = p4().pt();
  const float phi = p4().phi();
  const float eta = etaSC();
  const float seedE = eSeed();
  const float theta = 2.0*TMath::ATan(TMath::Exp(-1.*eta));
  const float seedEt = seedE * TMath::Sin(theta);
  const float SCrawEt = eSCRaw() * TMath::Sin(theta);

  const float mll = dilep_mass();
  const float numberVtx = nvtx();

  const float t_pt = tag_p4().pt();
  const float t_phi = tag_p4().phi();
  const float t_eta = tag_p4().eta();

  const float z_pt = dilep_p4().pt();
  const float z_phi = dilep_p4().phi();
  const float z_eta = dilep_p4().eta();
  const float z_y = dilep_p4().Rapidity();
  
  const float probe_dxyPV = dxyPV();
  const float probe_dxyPV_err = dxyPV_err();
  const float probe_dZ = dZ();
  const float probe_RelIso03 = RelIso03();
  const float probe_RelIso03EA = RelIso03EA();
  const float probe_RelIso03DB = RelIso03DB();
  const float probe_pfChargedHadronIso = pfChargedHadronIso();
  const float probe_pfPhotonIso = pfPhotonIso();
  const float probe_pfNeutralHadronIso = pfNeutralHadronIso();
  const float probe_tkIso = tkIso();
  const float probe_sumPUPt = sumPUPt();
  const float probe_sigmaIEtaIEta_full5x5 = sigmaIEtaIEta_full5x5();
  const float probe_dEtaIn = dEtaIn();
  const float probe_dPhiIn = dPhiIn();
  const float probe_hOverE = hOverE();
  const float probe_ip3d = ip3d();
  const float probe_ip3derr = ip3derr();
  const float probe_ecalEnergy = ecalEnergy();
  const float probe_eOverPIn = eOverPIn();
  const float probe_conv_vtx_flag = conv_vtx_flag();
  const float probe_exp_innerlayers = exp_innerlayers();
  const float probe_exp_outerlayers = exp_outerlayers();
  const float probe_charge = charge();
  const float probe_sccharge = sccharge();
  const float probe_ckf_charge = ckf_charge();
  //const float probe_trk_charge = trk_charge();
  const float probe_threeChargeAgree = threeChargeAgree();
  const float probe_mva = mva();
  const float probe_type = type();
  const float probe_ecalIso = ecalIso();
  const float probe_hcalIso = hcalIso();
  const float probe_sigmaIEtaIEta = sigmaIEtaIEta();
  const float probe_ecalPFClusterIso = ecalPFClusterIso();
  const float probe_hcalPFClusterIso = hcalPFClusterIso();
  //additional vars used in mva
  const float probe_ckf_laywithmeas = ckf_laywithmeas();
  const float probe_sigmaIPhiIPhi_full5x5 = sigmaIPhiIPhi_full5x5();
  const float probe_e1x5_full5x5 = e1x5_full5x5();
  const float probe_e5x5_full5x5 = e5x5_full5x5();
  const float probe_r9_full5x5 = r9_full5x5();
  const float probe_etaSCwidth = etaSCwidth();
  const float probe_phiSCwidth = phiSCwidth();
  const float probe_eSeed = eSeed();
  const float probe_eSCRaw = eSCRaw();
  const float probe_eSCPresh = eSCPresh();
  const float probe_ckf_chi2 = ckf_chi2();
  const float probe_ckf_ndof = ckf_ndof();
  const float probe_chi2 = chi2();
  const float probe_ndof = ndof();
  const float probe_fbrem = fbrem();
  const float probe_eOverPOut = eOverPOut();
  const float probe_dEtaOut = dEtaOut();
  const float probe_dPhiOut = dPhiOut();

  //find leading/trailing leptons
  float pt1 = -999;
  float eta1 = -999;
  float phi1 = -999;
  float pt2 = -999;
  float eta2 = -999;
  float phi2 = -999;
  if (pt > t_pt){
    pt1 = pt;
    eta1 = eta;
    phi1 = phi;
    pt2 = t_pt;
    eta2 = t_eta;
    phi2 = t_phi;
  } else {
    pt1 = t_pt;
    eta1 = t_eta;
    phi1 = t_phi;
    pt2 = pt;
    eta2 = eta;
    phi2 = phi;
  }
  
  TString EBEE = "";
  TString EBEEsign = "";
  
  //endcap
  if (fabs(eta) > 1.479) {
    EBEE = "EE";
    // if (fabs(eta) > 2.2) EBEE="EE3";
    // if (fabs(eta) > 1.8 && fabs(eta) < 2.2) EBEE="EE2";
    // if (fabs(eta) < 1.8) EBEE="EE1";
    if (eta > 0) EBEEsign = "EEpos";
    if (eta < 0) EBEEsign = "EEneg";
  }
  //barrel
  else if (fabs(eta)<1.479){
    EBEE = "EB";
    // if (fabs(eta) > 1.) EBEE="EB3";
    // if (fabs(eta) > 0.5 && fabs(eta) <1.0 ) EBEE="EB2";
    // if (fabs(eta) <0.5) EBEE="EB1";
    EBEEsign = "EB";
  }
  else return;
  
  plot1D(("h"+sel+"_pt"+EBEE).Data(), pt,  weight, h_1, "pt", 100, 0, 200);
  plot1D(("h"+sel+"_leading_pt"+EBEE).Data(), pt1,  weight, h_1, "pt", 50, 0, 100);
  plot1D(("h"+sel+"_trailing_pt"+EBEE).Data(), pt2,  weight, h_1, "pt", 50, 0, 100);
  plot1D(("h"+sel+"_seedEt"+EBEE).Data(), seedEt,  weight, h_1, "seed ET", 50, 0, 100);
  plot1D(("h"+sel+"_SCrawEt"+EBEE).Data(), SCrawEt,  weight, h_1, "raw SC ET", 50, 0, 100);
  plot1D(("h"+sel+"_eta").Data(), eta,  weight, h_1, "eta", 50, -2.5, 2.5);
  plot1D(("h"+sel+"_phi").Data(), phi,  weight, h_1, "phi", 50, -3.5, 3.5);
  
  plot1D(("h"+sel+"_relchiso"+EBEE).Data(), pfChargedHadronIso()/pt,  weight, h_1, "PFCh", 100, 0, 1);
  plot1D(("h"+sel+"_relemiso"+EBEE).Data(), pfPhotonIso()/pt,  weight, h_1, "PFEM", 100, 0, 1);
  plot1D(("h"+sel+"_relnhiso"+EBEE).Data(), pfNeutralHadronIso()/pt,  weight, h_1, "PFNh", 100, 0, 1);
  
  plot1D(("h"+sel+"_relECALiso"+EBEE).Data(), ecalIso()/pt,  weight, h_1, "ECAL RelIso", 100, 0, 1);
  plot1D(("h"+sel+"_relHCALiso"+EBEE).Data(), hcalIso()/pt,  weight, h_1, "HCAL RelIso", 100, 0, 1);
  plot1D(("h"+sel+"_relECALHCALiso"+EBEE).Data(), (ecalIso()+hcalIso())/pt,  weight, h_1, "HCAL RelIso", 100, 0, 1);
  
  plot1D(("h"+sel+"_sieie"+EBEE).Data(), sigmaIEtaIEta_full5x5(),  weight, h_1, "sieie_5x5", 100, 0, EBEE=="EB" ? 0.025 : 0.06);
  plot1D(("h"+sel+"_sipip"+EBEE).Data(), sigmaIPhiIPhi_full5x5(),  weight, h_1, "sipip_5x5", 100, 0, EBEE=="EB" ? 0.05 : 0.1);
  plot1D(("h"+sel+"_deta"+EBEEsign).Data(), dEtaIn(),  weight, h_1, "DeltaEta", 50, -0.04, 0.04);
  if (fabs(eta) > 1.57) plot1D(("h"+sel+"_deta"+EBEE).Data(), dEtaIn(),  weight, h_1, "DeltaEta", 50, -0.04, 0.04);
  plot1D(("h"+sel+"_dphi"+EBEEsign).Data(), dPhiIn(),  weight, h_1, "DeltaPhi", 50, -0.2, 0.2);
  if (fabs(eta) > 1.57) plot1D(("h"+sel+"_dphi"+EBEE).Data(), dPhiIn(),  weight, h_1, "DeltaPhi", 50, -0.2, 0.2);
  plot1D(("h"+sel+"_HoverE"+EBEE).Data(), hOverE(),  weight, h_1, "H/E", 50, 0, 0.5);
  if (fabs(eta) > 1.57) plot1D(("h"+sel+"_detaSeed"+EBEE).Data(), dEtaOut(),  weight, h_1, "DeltaEta", 50, -0.1, 0.1);
  if (fabs(eta) > 1.57) plot1D(("h"+sel+"_dphiSeed"+EBEE).Data(), dPhiOut(),  weight, h_1, "DeltaPhi", 50, -0.3, 0.3);
  plot1D(("h"+sel+"_detaSeed"+EBEEsign).Data(), dEtaOut(),  weight, h_1, "DeltaEta", 50, -0.1, 0.1);
  plot1D(("h"+sel+"_dphiSeed"+EBEEsign).Data(), dPhiOut(),  weight, h_1, "DeltaPhi", 50, -0.3, 0.3);
  
  plot1D(("h"+sel+"_PSEoverRawE"+EBEE).Data(), eSCPresh() / eSCRaw() ,  weight, h_1, "EPSoverERawSC", 50, 0, 1);

  // ---------

  //  plot1D(("h"+sel+"_VAR"+EBEE).Data(), VAR,  weight, h_1, "VAR", 50, 0, 100);

  plot1D(("h"+sel+"_dxyPV"+EBEE).Data(), probe_dxyPV,  weight, h_1, "dxyPV", 40, -0.1, 0.1);
  plot1D(("h"+sel+"_dxyPV_err"+EBEE).Data(), probe_dxyPV_err,  weight, h_1, "dxyPV_err", 50, 0, 0.04);
  plot1D(("h"+sel+"_dZ"+EBEE).Data(), probe_dZ,  weight, h_1, "dZ", 50,  EBEE=="EB" ? -0.05 : -0.1 , EBEE=="EB" ? 0.05 : 0.1);
  plot1D(("h"+sel+"_RelIso03"+EBEE).Data(), probe_RelIso03,  weight, h_1, "RelIso03", 50, 0, 1);
  plot1D(("h"+sel+"_RelIso03EA"+EBEE).Data(), probe_RelIso03EA,  weight, h_1, "RelIso03EA", 50, 0, 1);
  plot1D(("h"+sel+"_RelIso03DB"+EBEE).Data(), probe_RelIso03DB,  weight, h_1, "RelIso03DB", 50, 0, 1);
  plot1D(("h"+sel+"_pfChargedHadronIso"+EBEE).Data(), probe_pfChargedHadronIso,  weight, h_1, "pfChargedHadronIso", 40, 0, 20);
  plot1D(("h"+sel+"_pfPhotonIso"+EBEE).Data(), probe_pfPhotonIso,  weight, h_1, "pfPhotonIso", 40, 0, 20);
  plot1D(("h"+sel+"_pfNeutralHadronIso"+EBEE).Data(), probe_pfNeutralHadronIso,  weight, h_1, "pfNeutralHadronIso", 40, 0, 20);
  plot1D(("h"+sel+"_tkIso"+EBEE).Data(), probe_tkIso,  weight, h_1, "tkIso", 40, 0, 20);
  plot1D(("h"+sel+"_sumPUPt"+EBEE).Data(), probe_sumPUPt,  weight, h_1, "sumPUPt", 40, 0, 20);
  // plot1D(("h"+sel+"_sigmaIEtaIEta_full5x5"+EBEE).Data(), probe_sigmaIEtaIEta_full5x5,  weight, h_1, "sigmaIEtaIEta_full5x5", 50, 0, 100);
  // plot1D(("h"+sel+"_etaSC"+EBEE).Data(), probe_etaSC,  weight, h_1, "etaSC", 50, 0, 100);
  // plot1D(("h"+sel+"_dEtaIn"+EBEE).Data(), probe_dEtaIn,  weight, h_1, "dEtaIn", 50, 0, 100);
  // plot1D(("h"+sel+"_dPhiIn"+EBEE).Data(), probe_dPhiIn,  weight, h_1, "dPhiIn", 50, 0, 100);
  // plot1D(("h"+sel+"_hOverE"+EBEE).Data(), probe_hOverE,  weight, h_1, "hOverE", 50, 0, 100);
  if ( EBEE == "EE" ) plot1D(("h"+sel+"_ip3d"+EBEE).Data(), probe_ip3d,  weight, h_1, "ip3d", 50, -0.1, 0.1);
  plot1D(("h"+sel+"_ip3d"+EBEEsign).Data(), probe_ip3d,  weight, h_1, "ip3d", 50, -0.1, 0.1);
  plot1D(("h"+sel+"_ip3derr"+EBEE).Data(), probe_ip3derr,  weight, h_1, "ip3derr", 50, 0, 0.05);
  plot1D(("h"+sel+"_ecalEnergy"+EBEE).Data(), probe_ecalEnergy,  weight, h_1, "ecalEnergy", 50, 0, EBEE=="EB" ? 200 : 400);
  plot1D(("h"+sel+"_eOverPIn"+EBEE).Data(), probe_eOverPIn,  weight, h_1, "eOverPIn", 50, 0, 10);
  plot1D(("h"+sel+"_conv_vtx_flag"+EBEE).Data(), probe_conv_vtx_flag,  weight, h_1, "conv_vtx_flag", 2, 0, 2);
  plot1D(("h"+sel+"_exp_innerlayers"+EBEE).Data(), probe_exp_innerlayers,  weight, h_1, "exp_innerlayers", 4, 0, 4);
  plot1D(("h"+sel+"_exp_outerlayers"+EBEE).Data(), probe_exp_outerlayers,  weight, h_1, "exp_outerlayers", 5, 0, 5);
  plot1D(("h"+sel+"_charge"+EBEE).Data(), probe_charge,  weight, h_1, "charge", 3, -1, 2);
  plot1D(("h"+sel+"_sccharge"+EBEE).Data(), probe_sccharge,  weight, h_1, "sccharge", 3, -1, 2);
  plot1D(("h"+sel+"_ckf_charge"+EBEE).Data(), probe_ckf_charge,  weight, h_1, "ckf_charge", 3, -1, 2);
  //plot1D(("h"+sel+"_trk_charge"+EBEE).Data(), probe_trk_charge,  weight, h_1, "trk_charge", 3, -1, 2);
  plot1D(("h"+sel+"_threeChargeAgree"+EBEE).Data(), probe_threeChargeAgree,  weight, h_1, "threeChargeAgree", 2, 0, 1);
  plot1D(("h"+sel+"_mva"+EBEE).Data(), probe_mva,  weight, h_1, "mva", 50, -1, 1);
  plot1D(("h"+sel+"_type"+EBEE).Data(), probe_type,  weight, h_1, "type", 80, -20, 60);
  plot1D(("h"+sel+"_ecalIso"+EBEE).Data(), probe_ecalIso,  weight, h_1, "ecalIso", 20, 0, 10);
  plot1D(("h"+sel+"_hcalIso"+EBEE).Data(), probe_hcalIso,  weight, h_1, "hcalIso", 20, 0, 10);
  plot1D(("h"+sel+"_sigmaIEtaIEta"+EBEE).Data(), probe_sigmaIEtaIEta,  weight, h_1, "sigmaIEtaIEta", 100, 0, EBEE=="EB" ? 0.05 : 0.1);
  plot1D(("h"+sel+"_ecalPFClusterIso"+EBEE).Data(), probe_ecalPFClusterIso,  weight, h_1, "ecalPFClusterIso", 20, 0, 10);
  plot1D(("h"+sel+"_hcalPFClusterIso"+EBEE).Data(), probe_hcalPFClusterIso,  weight, h_1, "hcalPFClusterIso", 20, 0, 10);
  plot1D(("h"+sel+"_ckf_laywithmeas"+EBEE).Data(), probe_ckf_laywithmeas,  weight, h_1, "ckf_laywithmeas", 20, 0, 20);
  // plot1D(("h"+sel+"_sigmaIPhiIPhi_full5x5"+EBEE).Data(), probe_sigmaIPhiIPhi_full5x5,  weight, h_1, "sigmaIPhiIPhi_full5x5", 50, 0, 100);
  plot1D(("h"+sel+"_e1x5_full5x5"+EBEE).Data(), probe_e1x5_full5x5,  weight, h_1, "e1x5_full5x5", 50, 0, EBEE=="EB" ? 200 : 400);
  plot1D(("h"+sel+"_e5x5_full5x5"+EBEE).Data(), probe_e5x5_full5x5,  weight, h_1, "e5x5_full5x5", 50, 0, EBEE=="EB" ? 200 : 400);
  plot1D(("h"+sel+"_r9_full5x5"+EBEE).Data(), probe_r9_full5x5,  weight, h_1, "r9_full5x5", 96, 0, 1.2); 
  plot1D(("h"+sel+"_r9_full5x5_fineBin"+EBEE).Data(), probe_r9_full5x5,  weight, h_1, "r9_full5x5", 60, 0.4, 1); 
  plot1D(("h"+sel+"_etaSCwidth"+EBEE).Data(), probe_etaSCwidth,  weight, h_1, "etaSCwidth", 50, 0, 0.05);
  plot1D(("h"+sel+"_phiSCwidth"+EBEE).Data(), probe_phiSCwidth,  weight, h_1, "phiSCwidth", 50, 0, 0.5);
  // plot1D(("h"+sel+"_eSeed"+EBEE).Data(), probe_eSeed,  weight, h_1, "eSeed", 50, 0, 500);
  plot1D(("h"+sel+"_eSCRaw"+EBEE).Data(), probe_eSCRaw,  weight, h_1, "eSCRaw", 50, 0,  EBEE=="EB" ? 200 : 400);
  plot1D(("h"+sel+"_eSCPresh"+EBEE).Data(), probe_eSCPresh,  weight, h_1, "eSCPresh", 40, 0, 40);
  plot1D(("h"+sel+"_ckf_chi2"+EBEE).Data(), probe_ckf_chi2,  weight, h_1, "ckf_chi2", 50, 0, 100);
  plot1D(("h"+sel+"_ckf_ndof"+EBEE).Data(), probe_ckf_ndof,  weight, h_1, "ckf_ndof", 50, 0, 50);
  plot1D(("h"+sel+"_chi2"+EBEE).Data(), probe_chi2,  weight, h_1, "chi2", 50, 0, 200);
  plot1D(("h"+sel+"_ndof"+EBEE).Data(), probe_ndof,  weight, h_1, "ndof", 50, 0, 50);
  plot1D(("h"+sel+"_fbrem"+EBEE).Data(), probe_fbrem,  weight, h_1, "fbrem", 50, -1, 1); 
  plot1D(("h"+sel+"_eOverPOut"+EBEE).Data(), probe_eOverPOut,  weight, h_1, "eOverPOut", 25, 0, 10);
  // plot1D(("h"+sel+"_dEtaOut"+EBEE).Data(), probe_dEtaOut,  weight, h_1, "dEtaOut", 50, 0, 100);
  // plot1D(("h"+sel+"_dPhiOut"+EBEE).Data(), probe_dPhiOut,  weight, h_1, "dPhiOut", 50, 0, 100);

  //additional stuff
  if (probe_r9_full5x5 > 0.9) plot1D(("h"+sel+"_sigmaIPhiIPhi_full5x5_r9High"+EBEE).Data(), probe_sigmaIPhiIPhi_full5x5,  weight, h_1, "sigmaIPhiIPhi_full5x5", 100, 0, EBEE=="EB" ? 0.05 : 0.1);
  if (probe_r9_full5x5 < 0.9) plot1D(("h"+sel+"_sigmaIPhiIPhi_full5x5_r9Low"+EBEE).Data(), probe_sigmaIPhiIPhi_full5x5,  weight, h_1, "sigmaIPhiIPhi_full5x5", 100, 0, EBEE=="EB" ? 0.05 : 0.1);
  plot1D(("h"+sel+"_SIP3D"+EBEE).Data(), probe_ip3d/probe_ip3derr,  weight, h_1, "SIP3D", 60, -10, 10);
  plot1D(("h"+sel+"_mll"+EBEE).Data(), mll,  weight, h_1, "mll", 40, 80, 100); 
  if (probe_mva < -0.4) plot1D(("h"+sel+"_mll_mvaFail"+EBEE).Data(), mll,  weight, h_1, "mll", 40, 80, 100); 
  plot1D(("h"+sel+"_nvtx"+EBEE).Data(), numberVtx,  weight, h_1, "nvtx", 50, 0, 50); 
  plot1D(("h"+sel+"_eSCRawVSpt"+EBEE).Data(), probe_eSCRaw/pt,  weight, h_1, "eSCRaw/pt", 50, 0.5, 1.5);
  plot1D(("h"+sel+"_ckf_chi2_over_ckf_ndof"+EBEE).Data(), probe_chi2/probe_ndof,  weight, h_1, "chi2/ndof", 50, 0, 50);
  plot1D(("h"+sel+"_chi2_over_ndof"+EBEE).Data(), probe_chi2/probe_ndof,  weight, h_1, "ckf_chi2/ckf_ndof", 50, 0, 50);
  plot1D(("h"+sel+"_tagpt"+EBEE).Data(), t_pt,  weight, h_1, "tag pt", 50, 0, 100);
  plot1D(("h"+sel+"_tageta").Data(), t_eta,  weight, h_1, "tag eta", 50, -2.5, 2.5);
  plot1D(("h"+sel+"_tagphi").Data(), t_phi,  weight, h_1, "tag phi", 50, -3.5, 3.5);
  plot1D(("h"+sel+"_Zpt"+EBEE).Data(), z_pt,  weight, h_1, "Z pt", 50, 0, 100);
  plot1D(("h"+sel+"_Zeta").Data(), z_eta,  weight, h_1, "Z eta", 100, -5, 5);
  plot1D(("h"+sel+"_Zy").Data(), z_y,  weight, h_1, "Z rapidity", 100, -5, 5);
  plot1D(("h"+sel+"_Zphi").Data(), z_phi,  weight, h_1, "Z phi", 50, -3.5, 3.5);
  plot1D(("h"+sel+"_ooemoop"+EBEE).Data(), fabs( (1.0/ecalEnergy()) - (eOverPIn()/ecalEnergy())) ,  weight, h_1, "|1/E - 1/p|", 50, 0, 1);
  plot1D(("h"+sel+"_e1x5_over_e5x5"+EBEE).Data(), probe_e1x5_full5x5 / probe_e5x5_full5x5,  weight, h_1, "e1x5/e5x5",50, 0, 1);
  // plot1D(("h"+sel+"_e3x3"+EBEE).Data(), probe_r9_full5x5 * probe_eSCRaw,  weight, h_1, "e3x3",100, 0, EBEE=="EB" ? 200 : 400);
  // plot1D(("h"+sel+"_e3x3_over_e5x5"+EBEE).Data(), probe_r9_full5x5 * probe_eSCRaw / probe_e5x5_full5x5,  weight, h_1, "e3x3/e5x5",50, 0, 1);
  // plot1D(("h"+sel+"_e1x5_over_eSCRaw"+EBEE).Data(), probe_e1x5_full5x5 / probe_eSCRaw,  weight, h_1, "e1x5/eSCRaw",50, 0, 1);

  
  return;
  
}

void makeDilepPlots(std::map<std::string, TH1*> & h_1, TString sel, float weight = 1) {

  const float mll = dilep_mass();

  const LorentzVector el_p4 = p4();
  const float pt = el_p4.pt();
  const float eta =  el_p4.eta();
  const float phi =  el_p4.phi();

  const LorentzVector t_p4 = tag_p4();
  const float t_pt = t_p4.pt();
  const float t_eta =  t_p4.eta();
  const float t_phi =  t_p4.phi();

  const LorentzVector z_p4 = dilep_p4();
  const float z_pt = z_p4.pt();
  const float z_eta =  z_p4.eta();
  const float z_phi =  z_p4.phi();

  const bool inZWindow = (mll > 80) && (mll < 100);

  //find leading/trailing leptons
  float pt1 = -999;
  float eta1 = -999;
  float phi1 = -999;
  float pt2 = -999;
  float eta2 = -999;
  float phi2 = -999;
  if (pt > t_pt){
    pt1 = pt;
    eta1 = eta;
    phi1 = phi;
    pt2 = t_pt;
    eta2 = t_eta;
    phi2 = t_phi;
  } else {
    pt1 = t_pt;
    eta1 = t_eta;
    phi1 = t_phi;
    pt2 = pt;
    eta2 = eta;
    phi2 = phi;
  }
    
  TString EBEE = "";
  TString EBEEsign = "";
  
   if (fabs(eta) > 1.479) {
    EBEE = "EE";
    if (eta > 0) EBEEsign = "EEpos";
    if (eta < 0) EBEEsign = "EEneg";
  }
  else if (fabs(eta)<1.479){
    EBEE = "EB";
    EBEEsign = "EB";
  }
  else return;
  
  plot1D(("h"+sel+"_mll").Data(), mll,  weight, h_1, "m_{ll} [GeV]", 150, 0, 150); //dilepton mass
  if (tag_p4().pt()>25 && pt>15) {
    plot1D(("h"+sel+"_doubleE_mll").Data(), mll,  weight, h_1, "m_{ll} [GeV]", 150, 0, 150);//dilepton mass, for els past doubleE threshold
  }

  if (inZWindow) {
    plot1D(("h"+sel+"_Zprobe_pt").Data(), pt,  weight, h_1, "p_{T} [GeV]", 50, 0, 100); //probe pt
    plot1D(("h"+sel+"_Zprobe_pt"+EBEE).Data(), pt,  weight, h_1, "p_{T} [GeV]", 50, 0, 100); //probe pt
    plot1D(("h"+sel+"_Zprobe_eta").Data(), eta,  weight, h_1, "eta", 50, -2.5, 2.5); //probe eta
    plot1D(("h"+sel+"_Zprobe_phi").Data(), phi,  weight, h_1, "phi", 50, -3.5, 3.5); //probe phi
    
    plot1D(("h"+sel+"_Z_pt").Data(), z_pt,  weight, h_1, "Z p_{T} [GeV]", 50, 0, 100); //probe pt
    plot1D(("h"+sel+"_Z_eta").Data(), z_eta,  weight, h_1, "Z eta", 50, -2.5, 2.5); //probe eta
    plot1D(("h"+sel+"_Z_phi").Data(), z_phi,  weight, h_1, "Z phi", 50, -3.5, 3.5); //probe phi
    
    plot1D(("h"+sel+"_lead_pt").Data(), pt1,  weight, h_1, "p_{T} [GeV]", 100, 0, 200); //lead pt
    plot1D(("h"+sel+"_lead_pt"+EBEE).Data(), pt1,  weight, h_1, "p_{T} [GeV]", 100, 0, 200); //lead pt
    plot1D(("h"+sel+"_lead_eta").Data(), eta1,  weight, h_1, "eta", 50, -2.5, 2.5); //lead eta
    plot1D(("h"+sel+"_lead_phi").Data(), phi1,  weight, h_1, "phi", 50, -3.5, 3.5); //lead phi

    plot1D(("h"+sel+"_trail_pt").Data(), pt2,  weight, h_1, "p_{T} [GeV]", 100, 0, 200); //trail pt
    plot1D(("h"+sel+"_trail_pt"+EBEE).Data(), pt2,  weight, h_1, "p_{T} [GeV]", 100, 0, 200); //trail pt	
    plot1D(("h"+sel+"_trail_eta").Data(), eta2,  weight, h_1, "eta", 50, -2.5, 2.5); //trail eta
    plot1D(("h"+sel+"_trail_phi").Data(), phi2,  weight, h_1, "phi", 50, -3.5, 3.5); //trail phi
    
  }
  
}
