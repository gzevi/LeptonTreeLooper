{

  gROOT->ProcessLine(".L LeptonTree.cc+");
  gROOT->ProcessLine(".L PlotUtilities.cc+");
  gROOT->ProcessLine(".L LeptonTreeLooper.C+");

  TChain *ch = new TChain("t");

 
  // ch->Reset();
  // ch->Add("/nfs-7/userdata/leptonTree/v0.09/SingleElectronAll.root");
  // LeptonTreeLooper(ch, "fullSingleElectron_splitEta", -1);


  // ch->Reset();
  // ch->Add("../LeptonBabyMaker/WW_test.root");
  // LeptonTreeLooper(ch, "WW_test", -1);
  ch->Reset();
  ch->Add("../LeptonBabyMaker/WW_old.root");
  LeptonTreeLooper(ch, "WW", -1);

}
