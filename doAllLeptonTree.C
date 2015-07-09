{

  gROOT->ProcessLine(".L LeptonTree.cc+");
  gROOT->ProcessLine(".L PlotUtilities.cc+");
  gROOT->ProcessLine(".L LeptonTreeLooper.C+");

  TChain *ch = new TChain("t");
  ch->Add("/hadoop/cms/store/user/mderdzinski/leptonBabies/commissioning/merged_files/DrellYan_test.root");
  LeptonTreeLooper(ch, "plots", -1);
  
}
