{

  gROOT->ProcessLine(".L LeptonTree.cc+");
  gROOT->ProcessLine(".L PlotUtilities.cc+");
  gROOT->ProcessLine(".L LeptonTreeLooper.C+");

  TChain *ch = new TChain("t");

  // ch->Add("/nfs-7/userdata/leptonTree/v2.00TP/DY_madgraph80X.root");
  // LeptonTreeLooper(ch, "DY_jun8_puReweight", -1);
  // ch->Add("/nfs-7/userdata/leptonTree/v2.01TP/2016SingleEl.root");
  // LeptonTreeLooper(ch, "SingleEl_jun8_puReweight", -1);

  // ch->Add("/nfs-7/userdata/leptonTree/v2.04TP/DY_madgraph80X*");
  // LeptonTreeLooper(ch, "DY_Jul19_puReweight_Iso05", -1);
  
  // ch->Add("/nfs-7/userdata/leptonTree/v2.05TPPF/DY_madgraph80X*");
  // LeptonTreeLooper(ch, "DY_Jul19_puReweight_newBabies_XXXXXX", -1);

  // ch->Add("/nfs-7/userdata/leptonTree/v2.05TPPF/DY_madgraph80X.root");
  // LeptonTreeLooper(ch, "DY_Jul19_puReweight_newBabies_subset_test", -1);
  
  // ch->Add("/nfs-7/userdata/leptonTree/v2.05TPPF/2016SingleEl*");
  // LeptonTreeLooper(ch, "SingleEl_Jul19_RunC", -1);

  // ch->Add("/nfs-7/userdata/leptonTree/v2.07EGM20GeV/DY_madgraph80X*.root");
  // LeptonTreeLooper(ch, "DY_Jul22_puReweight_207_Iso05", -1);  
  // ch->Add("/nfs-7/userdata/leptonTree/v2.07EGM20GeV/2016SingleEl*");
  // LeptonTreeLooper(ch, "SingleEl_Jul22_207_Iso05", -1);

  // ch->Add("/nfs-7/userdata/leptonTree/v2.08EGM20GeV/DY_madgraph80X*.root");
  // LeptonTreeLooper(ch, "DY_Jul27_puReweight_208", -1);  
  // ch->Add("/nfs-7/userdata/leptonTree/v2.09EGM20GeV/2016SingleEl*");
  // LeptonTreeLooper(ch, "SingleEl_Jul27_209", -1);

  // ch->Add("/nfs-7/userdata/leptonTree/v2.09EGM20GeV/QCD*");
  // LeptonTreeLooper(ch, "QCD_Jul27_209", -1);
  // ch->Add("/nfs-7/userdata/leptonTree/v2.09EGM20GeV/WJets*");
  // LeptonTreeLooper(ch, "WJets_Jul27_209", -1);

  
  // ch->Add("/nfs-7/userdata/leptonTree/v2.08EGM20GeV/DY_madgraph80X.root");
  // LeptonTreeLooper(ch, "DY_Jul26_test", -1);
  // ch->Add("/nfs-7/userdata/leptonTree/v2.08EGM20GeV/2016SingleEl.root");
  // LeptonTreeLooper(ch, "SingleEl_Jul26_test", -1);


  //-----------------------
  //november rereco testing
  //-----------------------
  
  // ch->Add("/nfs-7/userdata/leptonTree/v2.08EGM20GeV/DY_madgraph80X*.root");
  // LeptonTreeLooper(ch, "DY_Nov13_puReweight2016H_208", -1);  

  // ch->Add("/nfs-7/userdata/leptonTree/v2.10EGM20GeV/2016*SingleEl_P*root");
  // LeptonTreeLooper(ch, "Nov13_2016-SingleEl_P", -1);
  // ch->Add("/nfs-7/userdata/leptonTree/v2.10EGM20GeV/2016*SingleEl_R*root");
  // LeptonTreeLooper(ch, "Nov13_2016-SingleEl_R", -1);

  // ch->Add("/nfs-7/userdata/leptonTree/v2.11EGM20GeV/2016C*SingleEl_R*root");
  // LeptonTreeLooper(ch, "Nov13_v211_2016C-SingleEl_R", -1);
  // ch->Reset();
  // ch->Add("/nfs-7/userdata/leptonTree/v2.11EGM20GeV/2016D*SingleEl_R*root");
  // LeptonTreeLooper(ch, "Nov13_v211_2016D-SingleEl_R", -1);
  // ch->Reset();
  // ch->Add("/nfs-7/userdata/leptonTree/v2.11EGM20GeV/2016E*SingleEl_R*root");
  // LeptonTreeLooper(ch, "Nov13_v211_2016E-SingleEl_R", -1);
  // ch->Reset();
  // ch->Add("/nfs-7/userdata/leptonTree/v2.11EGM20GeV/2016F*SingleEl_R*root");
  // LeptonTreeLooper(ch, "Nov13_v211_2016F-SingleEl_R", -1);
  // ch->Reset();
  // ch->Add("/nfs-7/userdata/leptonTree/v2.11EGM20GeV/2016G*SingleEl_R*root");
  // LeptonTreeLooper(ch, "Nov13_v211_2016G-SingleEl_R", -1);

  // ch->Add("/nfs-7/userdata/leptonTree/v2.11EGM20GeV/2016B*SingleEl_P*root");
  // LeptonTreeLooper(ch, "Nov13_v211_2016B-SingleEl_P", -1);
  // ch->Reset();
  // ch->Add("/nfs-7/userdata/leptonTree/v2.11EGM20GeV/2016C*SingleEl_P*root");
  // LeptonTreeLooper(ch, "Nov13_v211_2016C-SingleEl_P", -1);
  // ch->Reset();
  // ch->Add("/nfs-7/userdata/leptonTree/v2.11EGM20GeV/2016D*SingleEl_P*root");
  // LeptonTreeLooper(ch, "Nov13_v211_2016D-SingleEl_P", -1);
  // ch->Reset();
  // ch->Add("/nfs-7/userdata/leptonTree/v2.11EGM20GeV/2016E*SingleEl_P*root");
  // LeptonTreeLooper(ch, "Nov13_v211_2016E-SingleEl_P", -1);
  // ch->Reset();
  // ch->Add("/nfs-7/userdata/leptonTree/v2.11EGM20GeV/2016F*SingleEl_P*root");
  // LeptonTreeLooper(ch, "Nov13_v211_2016F-SingleEl_P", -1);
  // ch->Reset();
  // ch->Add("/nfs-7/userdata/leptonTree/v2.11EGM20GeV/2016G*SingleEl_P*root");
  // LeptonTreeLooper(ch, "Nov13_v211_2016G-SingleEl_P", -1);
  // ch->Reset();
  // ch->Add("/nfs-7/userdata/leptonTree/v2.11EGM20GeV/2016H*SingleEl_P*root");
  // LeptonTreeLooper(ch, "Nov13_v211_2016H-SingleEl_P", -1);

  //-----------------------
  //February moriond testing
  //-----------------------
  
  // ch->Add("/nfs-7/userdata/leptonTree/ReRecoDataMoriondMC_5PF_17Feb17/DY_madgraph*.root");
  // LeptonTreeLooper(ch, "DY_Feb23_puReweight_fixScaleSmear_moriond", -1);  

  ch->Add("/nfs-7/userdata/leptonTree/ReRecoDataMoriondMC_5PF_17Feb17/2016SingleEl*.root");
  LeptonTreeLooper(ch, "2016SingleEl_Feb23_moriond", -1);  
 
}
