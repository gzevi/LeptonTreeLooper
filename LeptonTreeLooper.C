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
  const float lumi = .042;

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
  // trigNames.push_back("HLT_Ele23_WPLoose_Gsf");

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
      if (evt_isRealData() && !goodrun(evt_run(), evt_lumiBlock())) continue; 

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
      // trigDecision.push_back(HLT_Ele23_WPLoose_Gsf());                                              tagLL.push_back(/*HLT_Ele23_WPLoose_Gsf()*/0);                                                     trigPtPlat.push_back(0);


      //make sure all trigger vectors are same size
      if ( trigNames.size() != trigDecision.size() && trigNames.size() != tagLL.size() ){ cout << "Trigger vectors are not equal size! Continuing." << endl; continue; }

      //if real data, set lumiScale to 1
      float lumiScale = scale1fb() * lumi;
      //float lumiScale = scale1fb() * lumi * 20; //multiply by 20 for DYtest2
      if (evt_isRealData()) lumiScale = 1;      

      //if not data, calculate weight based on PUreweighting
      if (!evt_isRealData()){
	const int vtx_bin = h_ratio->GetXaxis()->FindBin(nvtx());
	float puWeight = h_ratio->GetBinContent(vtx_bin);
	lumiScale *= puWeight;
      }
      
      //skip leptons where a (different) tag doesn't exist in event
      const float mll = dilep_mass();
      if ( mll == -1 ) continue;

      //tag requirements
      if (tag_p4().pt() < 30 ) continue;
      if (tag_p4().eta() > 2.1) continue;      
      
      //things to only fill once per event
      const int evt = evt_event();
      if( lastEventSaved_ != evt ){
	
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
    
    plot1D(("h"+sel+"_lead_pt").Data(), pt1,  weight, h_1, "p_{T} [GeV]", 50, 0, 100); //lead pt
    plot1D(("h"+sel+"_lead_pt"+EBEE).Data(), pt1,  weight, h_1, "p_{T} [GeV]", 50, 0, 100); //lead pt
    plot1D(("h"+sel+"_lead_eta").Data(), eta1,  weight, h_1, "eta", 50, -2.5, 2.5); //lead eta
    plot1D(("h"+sel+"_lead_phi").Data(), phi1,  weight, h_1, "phi", 50, -3.5, 3.5); //lead phi

    plot1D(("h"+sel+"_trail_pt").Data(), pt2,  weight, h_1, "p_{T} [GeV]", 50, 0, 100); //trail pt
    plot1D(("h"+sel+"_trail_pt"+EBEE).Data(), pt2,  weight, h_1, "p_{T} [GeV]", 50, 0, 100); //trail pt	
    plot1D(("h"+sel+"_trail_eta").Data(), eta2,  weight, h_1, "eta", 50, -2.5, 2.5); //trail eta
    plot1D(("h"+sel+"_trail_phi").Data(), phi2,  weight, h_1, "phi", 50, -3.5, 3.5); //trail phi
    
  }
  
}
