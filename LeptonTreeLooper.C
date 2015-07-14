#include "LeptonTreeLooper.h"

int LeptonTreeLooper( TChain* chain, TString output_name , int nEvents );
void makePlots(std::map<std::string, TH1*> & h_1, TString sel, float weight);


int LeptonTreeLooper( TChain* chain, TString output_name , int nEvents ) {

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
  const float lumi = .006;
  
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
      std::vector<int> trigDecision;
      trigDecision.push_back(1);
      trigDecision.push_back(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ());
      trigDecision.push_back(HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300());
      trigDecision.push_back(HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p5PF());
      trigDecision.push_back(HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30());
      trigDecision.push_back(HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30());
      trigDecision.push_back(HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30());
      trigDecision.push_back(HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30());
      trigDecision.push_back(HLT_Ele33_CaloIdM_TrackIdM_PFJet30());
      trigDecision.push_back(HLT_Ele23_CaloIdM_TrackIdM_PFJet30());
      trigDecision.push_back(HLT_Ele18_CaloIdM_TrackIdM_PFJet30());
      trigDecision.push_back(HLT_Ele12_CaloIdM_TrackIdM_PFJet30());
      trigDecision.push_back(HLT_Ele8_CaloIdM_TrackIdM_PFJet30());
      trigDecision.push_back(tag_HLT_Ele25WP60_Ele8_Mass55_LeadingLeg());
      trigDecision.push_back(tag_HLT_Ele25WP60_SC4_Mass55_LeadingLeg());

      //vector for tag LL
      //need same number of entries and same order as trigNames!!!
      std::vector<int> tagLL;
      tagLL.push_back(0);
      tagLL.push_back(/*tag_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_ElectronLeg()*/ 0);
      tagLL.push_back(/*tag_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_ElectronLeg()*/ 0);
      tagLL.push_back(/*tag_HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p5PF_ElectronLeg()*/ 0);
      tagLL.push_back(tag_HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_ElectronLeg());
      tagLL.push_back(tag_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_ElectronLeg());
      tagLL.push_back(tag_HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30_ElectronLeg());
      tagLL.push_back(tag_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_ElectronLeg());
      tagLL.push_back(tag_HLT_Ele33_CaloIdM_TrackIdM_PFJet30_ElectronLeg());
      tagLL.push_back(tag_HLT_Ele23_CaloIdM_TrackIdM_PFJet30_ElectronLeg());
      tagLL.push_back(tag_HLT_Ele18_CaloIdM_TrackIdM_PFJet30_ElectronLeg());
      tagLL.push_back(tag_HLT_Ele12_CaloIdM_TrackIdM_PFJet30_ElectronLeg());
      tagLL.push_back(tag_HLT_Ele8_CaloIdM_TrackIdM_PFJet30_ElectronLeg());
      tagLL.push_back(tag_HLT_Ele25WP60_Ele8_Mass55_LeadingLeg());
      tagLL.push_back(tag_HLT_Ele25WP60_SC4_Mass55_LeadingLeg());

      //skip leptons where a (different) tag doesn't exist in event
      float mll = dilep_mass();
      if ( mll == -1 ) continue;

      const float lumiScale = scale1fb() * lumi;
      //loop over triggers
      for ( int trigIdx=0; trigIdx < trigNames.size(); trigIdx++) {
	if ( trigDecision[trigIdx] != 0) {
	  makePlots( h_1d, trigNames[trigIdx], 1);
	  makePlots( h_1d, trigNames[trigIdx]+"_scaled", lumiScale);
	  if ( tagLL[trigIdx] != 0 ) makePlots( h_1d, trigNames[trigIdx]+"_tagLL", 1);
	}//trigDecision	
      }//trig loop

      
      //things to only fill once per event
      int evt = evt_event();
      if( lastEventSaved_ != evt ){
	for ( int trigIdx=0; trigIdx < trigNames.size(); trigIdx++) {
	  if ( trigDecision[trigIdx] != 0) {
	    plot1D(("h"+trigNames[trigIdx]+"_mll").Data(), mll,  1, h_1d, "m_{ll} [GeV]", 150, 0, 150); //dilepton mass
	    plot1D(("h"+trigNames[trigIdx]+"_mll_pt").Data(), pt,  1, h_1d, "p_{T} [GeV]", 50, 0, 100); //probe pt
	    plot1D(("h"+trigNames[trigIdx]+"_mll_eta").Data(), eta,  1, h_1d, "eta", 50, -2.5, 2.5); //probe eta
	    plot1D(("h"+trigNames[trigIdx]+"_scaled_mll").Data(), mll,  lumiScale, h_1d, "m_{ll} [GeV]", 150, 0, 150); //dilepton mass
	    plot1D("h_trigs",trigIdx-0.5, 1, h_1d, "trigger",trigNames.size(),-0.5,trigNames.size()-0.5); //counts how often each trigger fired
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
  
  float pt = p4().pt();
  float phi = p4().phi();
  float eta = etaSC();
  //float seedE = eSeed();
  float theta = 2.0*TMath::ATan(TMath::Exp(-1.*eta));
  //float seedEt = seedE * TMath::Sin(theta);
  float SCrawEt = eSCRaw() * TMath::Sin(theta);

  //float mll = dilep_mass();

  TString EBEE = "";
  TString EBEEsign = "";
  
  if (fabs(eta) > 1.57) {
    EBEE = "EE";
    if (eta > 0) EBEEsign = "EEpos";
    if (eta < 0) EBEEsign = "EEneg";
  }
  else if (fabs(eta)<1.44){
    EBEE = "EB";
    EBEEsign = "EB";
  }
  else return;
  
  plot1D(("h"+sel+"_pt"+EBEE).Data(), pt,  weight, h_1, "pt", 50, 0, 100);
  //plot1D(("h"+sel+"_seedEt"+EBEE).Data(), seedEt,  weight, h_1, "seed ET", 50, 0, 100);
  plot1D(("h"+sel+"_SCrawEt"+EBEE).Data(), SCrawEt,  weight, h_1, "raw SC ET", 50, 0, 100);
  plot1D(("h"+sel+"_eta").Data(), eta,  weight, h_1, "eta", 50, -2.5, 2.5);
  plot1D(("h"+sel+"_phi").Data(), phi,  weight, h_1, "phi", 50, -3.5, 3.5);
  
  // plot1D(("h"+sel+"_relchiso"+EBEE).Data(), pfChargedHadronIso()/seedEt,  weight, h_1, "PFCh", 100, 0, 1);
  // plot1D(("h"+sel+"_relemiso"+EBEE).Data(), pfPhotonIso()/seedEt,  weight, h_1, "PFEM", 100, 0, 1);
  // plot1D(("h"+sel+"_relnhiso"+EBEE).Data(), pfNeutralHadronIso()/seedEt,  weight, h_1, "PFNh", 100, 0, 1);
  
  // plot1D(("h"+sel+"_relECALiso"+EBEE).Data(), ecalIso()/seedEt,  weight, h_1, "ECAL RelIso", 100, 0, 1);
  // plot1D(("h"+sel+"_relHCALiso"+EBEE).Data(), hcalIso()/seedEt,  weight, h_1, "HCAL RelIso", 100, 0, 1);
  // plot1D(("h"+sel+"_relECALHCALiso"+EBEE).Data(), (ecalIso()+hcalIso())/seedEt,  weight, h_1, "HCAL RelIso", 100, 0, 1);
  
  plot1D(("h"+sel+"_sieie"+EBEE).Data(), sigmaIEtaIEta_full5x5(),  weight, h_1, "sieie_5x5", 100, 0, EBEE=="EB" ? 0.05 : 0.1);
  plot1D(("h"+sel+"_sipip"+EBEE).Data(), sigmaIPhiIPhi_full5x5(),  weight, h_1, "sipip_5x5", 100, 0, EBEE=="EB" ? 0.05 : 0.1);
  plot1D(("h"+sel+"_deta"+EBEEsign).Data(), dEtaIn(),  weight, h_1, "DeltaEta", 50, -0.1, 0.1);
  if (fabs(eta) > 1.57) plot1D(("h"+sel+"_deta"+EBEE).Data(), dEtaIn(),  weight, h_1, "DeltaEta", 50, -0.1, 0.1);
  plot1D(("h"+sel+"_dphi"+EBEEsign).Data(), dPhiIn(),  weight, h_1, "DeltaPhi", 50, -0.1, 0.1);
  if (fabs(eta) > 1.57) plot1D(("h"+sel+"_dphi"+EBEE).Data(), dPhiIn(),  weight, h_1, "DeltaPhi", 50, -0.1, 0.1);
  plot1D(("h"+sel+"_HoverE"+EBEE).Data(), hOverE(),  weight, h_1, "H/E", 50, 0, 0.5);
  if (fabs(eta) > 1.57) plot1D(("h"+sel+"_detaSeed"+EBEE).Data(), dEtaOut(),  weight, h_1, "DeltaEta", 50, -0.1, 0.1);
  //if (fabs(eta) > 1.57) plot1D(("h"+sel+"_dphiSeed"+EBEE).Data(), dPhiOut(),  weight, h_1, "DeltaPhi", 50, -0.1, 0.1);
  plot1D(("h"+sel+"_detaSeed"+EBEEsign).Data(), dEtaOut(),  weight, h_1, "DeltaEta", 50, -0.1, 0.1);
  //plot1D(("h"+sel+"_dphiSeed"+EBEEsign).Data(), dPhiOut(),  weight, h_1, "DeltaPhi", 50, -0.1, 0.1);
  
  // plot1D(("h"+sel+"_PSEoverRawE"+EBEE).Data(), eSCPresh() / eSCRaw() ,  weight, h_1, "EPSoverERawSC", 50, 0, 1);

  //plot1D("h"+sel+"_mll", mll,  weight, h_1, "m_{ll} [GeV]", 150, 0, 150);

  return;
  
}
