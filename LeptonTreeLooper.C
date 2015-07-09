
#include "LeptonTreeLooper.h"

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

      float mll = dilep_mass();
      if ( mll == -1 ) continue;

      int evt = evt_event();

      if( lastEventSaved_ != evt ){

      	//global dilep mass
      	//save only the first dilep mass hyp for event.
      	//Less than 2% efficiency loss for Z window in DY MC by not taking best Z mll hyp
      	plot1D("h_global_mll", mll,  1, h_1d, "m_{ll} [GeV]", 150, 0, 150);
	
      	//trigger specific dilep mass
      	if ( HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() != 0 ){ plot1D("h_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_mll", mll,  1, h_1d, "m_{ll} [GeV]", 150, 0, 150); }
      	if ( HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300() != 0 ){ plot1D("h_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_mll", mll,  1, h_1d, "m_{ll} [GeV]", 150, 0, 150); }
      	if ( HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p5PF() != 0 ){ plot1D("h_HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p5PF_mll", mll,  1, h_1d, "m_{ll} [GeV]", 150, 0, 150); }
      	if ( HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30() != 0 ){ plot1D("h_HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_mll", mll,  1, h_1d, "m_{ll} [GeV]", 150, 0, 150); }
      	if ( HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30() != 0 ){ plot1D("h_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_mll", mll,  1, h_1d, "m_{ll} [GeV]", 150, 0, 150); }
      	if ( HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30() != 0 ){ plot1D("h_HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30_mll", mll,  1, h_1d, "m_{ll} [GeV]", 150, 0, 150); }
      	if ( HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30() != 0 ){ plot1D("h_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_mll", mll,  1, h_1d, "m_{ll} [GeV]", 150, 0, 150); }
      	if ( HLT_Ele33_CaloIdM_TrackIdM_PFJet30() != 0 ){ plot1D("h_HLT_Ele33_CaloIdM_TrackIdM_PFJet30_mll", mll,  1, h_1d, "m_{ll} [GeV]", 150, 0, 150); }
      	if ( HLT_Ele23_CaloIdM_TrackIdM_PFJet30() != 0 ){ plot1D("h_HLT_Ele23_CaloIdM_TrackIdM_PFJet30_mll", mll,  1, h_1d, "m_{ll} [GeV]", 150, 0, 150); }
      	if ( HLT_Ele18_CaloIdM_TrackIdM_PFJet30() != 0 ){ plot1D("h_HLT_Ele18_CaloIdM_TrackIdM_PFJet30_mll", mll,  1, h_1d, "m_{ll} [GeV]", 150, 0, 150); }
      	if ( HLT_Ele12_CaloIdM_TrackIdM_PFJet30() != 0 ){ plot1D("h_HLT_Ele12_CaloIdM_TrackIdM_PFJet30_mll", mll,  1, h_1d, "m_{ll} [GeV]", 150, 0, 150); }
      	if ( HLT_Ele8_CaloIdM_TrackIdM_PFJet30() != 0 ){ plot1D("h_HLT_Ele8_CaloIdM_TrackIdM_PFJet30_mll", mll,  1, h_1d, "m_{ll} [GeV]", 150, 0, 150); }
	
      	lastEventSaved_ = evt;	
      } //evtSaved
      
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
