{

  gROOT->ProcessLine(".L LeptonTree.cc++");
  gROOT->ProcessLine(".L PlotUtilities.cc++");
  gROOT->ProcessLine(".L LeptonTreeLooper.C++");

  TChain *ch = new TChain("t");
  ch->Add("LowMassDY_50ns.root");
  LeptonTreeLooper(ch, "plots", -1);
  
}