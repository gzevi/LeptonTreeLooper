#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TString.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TEfficiency.h"
#include "tdrstyle_SUSY.C"


void makeEff( TChain* t_ele, TString var, TCut den, TCut num, int nbins, float min, float max, TString legend, TString xaxis, TString yaxis, TString saveas, bool log = true) {
  
  TH1D* h_pt_den = new TH1D("h_pt_den",";"+xaxis,nbins,min,max);
  TH1D* h_pt_num = (TH1D*) h_pt_den->Clone("h_pt_num");
  TCanvas* c = new TCanvas("c","c");
  c->SetGrid(1,1);
  c->cd();
  
  t_ele->Draw(var+">>h_pt_den",den);
  t_ele->Draw(var+">>h_pt_num",num);
  
  TH2F* h_axis = new TH2F("h_axis",";"+xaxis+";"+yaxis,nbins,min,max,20,0,1);
  h_axis->Draw();
  
  TEfficiency* h_pt_eff = new TEfficiency(*h_pt_num, *h_pt_den);
  h_pt_eff->SetLineColor(kBlue);
  h_pt_eff->SetMarkerColor(kBlue);
  
  h_pt_eff->Draw("pe same");
  
  TLegend* leg = new TLegend(0.6,0.2,0.9,0.5);
  leg->AddEntry(h_pt_eff,legend,"pe");
  leg->Draw("same");
  
  c->SaveAs(saveas+".pdf");
  c->Close();
  h_axis->Delete();
  leg->Delete();
  
  c = new TCanvas("c","c");
  c->SetGrid(1,1);
  c->cd();
  if (log) h_axis = new TH2F("h_axis",";"+xaxis+";Events",nbins,min,max,20,0.01,h_pt_den->GetMaximum()*1.3);
  else h_axis = new TH2F("h_axis",";"+xaxis+";"+yaxis,nbins,min,max,20,0,h_pt_den->GetMaximum()*1.3);
  h_axis->Draw();
  if (log) c->SetLogy();
  h_pt_num->SetLineColor(kBlue);
  h_pt_num->SetMarkerColor(kBlue);
  h_pt_den->Draw("pe same");
  h_pt_num->Draw("pe same");
  leg = new TLegend(0.65,0.3,0.96,0.5);
  leg->AddEntry(h_pt_den,"Pass Denominator","pe");
  leg->AddEntry(h_pt_num,legend,"pe");
  leg->Draw("same");
  c->SaveAs("Overlay"+saveas+".pdf");
  c->Close();
  
  h_pt_den->Delete();
  h_pt_num->Delete();
  h_pt_eff->Delete();

}


void makeOverlay( TChain* t_ele, TString var, TCut den, TCut num, int nbins, float min, float max, TString legend, TString xaxis, TString yaxis, TString saveas) {
  
  TH1D* h_pt_den = new TH1D("h_pt_den",";"+xaxis,nbins,min,max);
  TH1D* h_pt_num = (TH1D*) h_pt_den->Clone("h_pt_num");
  TCanvas* c = new TCanvas("c","c");
  c->SetGrid(1,1);
  c->cd();
  
  t_ele->Draw(var+">>h_pt_den",den);
  t_ele->Draw(var+">>h_pt_num",num);
  
  TH2F* h_axis = new TH2F("h_axis",";"+xaxis+";"+yaxis,nbins,min,max,20,0,h_pt_den->GetMaximum()*1.3);
  h_axis->Draw();
  
  h_pt_num->SetLineColor(kBlue);
  h_pt_num->SetMarkerColor(kBlue);
  
  h_pt_den->Draw("pe same");
  h_pt_num->Draw("pe same");
  
  TLegend* leg = new TLegend(0.65,0.3,0.96,0.5);
  leg->AddEntry(h_pt_den,"Pass Denominator","pe");
  leg->AddEntry(h_pt_num,legend,"pe");
  leg->Draw("same");
  
  c->SaveAs(saveas+".pdf");
  c->Close();
  h_pt_den->Delete();
  h_pt_num->Delete();
  
}

void plot_trigeff_Ele (const TString& indir = "babies/") {

//  gROOT->ProcessLine(".L ./tdrstyle_SUSY.C");
  
  setTDRStyle();
  
  TH1::SetDefaultSumw2();
  
  TChain* t_ele = new TChain("t");

  t_ele->Add(Form("%s/2016SingleEl273017full.root", indir.Data()));

  TFile* f_out = new TFile("trigeff_Ele.root","RECREATE");
  
  TH1D* h_pt_denom_ele25WPTight = new TH1D("h_pt_denom_ele25WPTight",";electron p_{T} [GeV]",50,10,60);
  TH1D* h_pt_denom_eb_ele25WPTight = (TH1D*) h_pt_denom_ele25WPTight->Clone("h_pt_denom_eb_ele25WPTight");
  TH1D* h_pt_denom_ee_ele25WPTight = (TH1D*) h_pt_denom_ele25WPTight->Clone("h_pt_denom_ee_ele25WPTight");
  TH1D* h_pt_num_ele25WPTight = (TH1D*) h_pt_denom_ele25WPTight->Clone("h_pt_num_ele25WPTight");
  TH1D* h_pt_num_eb_ele25WPTight = (TH1D*) h_pt_denom_ele25WPTight->Clone("h_pt_num_eb_ele25WPTight");
  TH1D* h_pt_num_ee_ele25WPTight = (TH1D*) h_pt_denom_ele25WPTight->Clone("h_pt_num_ee_ele25WPTight");

  TCanvas* c = new TCanvas("c","c");
  c->SetGrid(1,1);
  c->cd();

  TCut tag = " p4.pt()>20 && abs(p4.eta())<2.5 && tag_p4.pt()>25 && tag_HLT_Ele23_WPLoose_Gsf>0";
  TCut base = "p4.pt()>20 && abs(p4.eta())<2.5 && dilep_mass>70 && dilep_mass<110 && tag_p4.pt()>25 && tag_HLT_Ele23_WPLoose_Gsf>0";
  TCut ID = "passes_POG_tightID";
  TCut eb = base+"abs(p4.eta()) < 1.412";
  TCut ee = base+"abs(p4.eta()) > 1.479 && abs(p4.eta()) < 2.1";

  t_ele->Draw("p4.pt()>>h_pt_denom_ele25WPTight",base);
  t_ele->Draw("p4.pt()>>h_pt_num_ele25WPTight",base+"HLT_Ele25_eta2p1_WPTight_Gsf>0");

  t_ele->Draw("p4.pt()>>h_pt_denom_eb_ele25WPTight",base+eb);
  t_ele->Draw("p4.pt()>>h_pt_num_eb_ele25WPTight",base+eb+"HLT_Ele25_eta2p1_WPTight_Gsf>0");

  t_ele->Draw("p4.pt()>>h_pt_denom_ee_ele25WPTight",base+ee);
  t_ele->Draw("p4.pt()>>h_pt_num_ee_ele25WPTight",base+ee+"HLT_Ele25_eta2p1_WPTight_Gsf>0");

  TH2F* h_axis = new TH2F("h_axis",";electron p_{T} [GeV];Efficiency of Ele25_eta2p1_WPTight",50,10,60,20,0,1);
  h_axis->Draw();
  
  TEfficiency* h_pt_eff_ele25WPTight = new TEfficiency(*h_pt_num_ele25WPTight, *h_pt_denom_ele25WPTight);
  h_pt_eff_ele25WPTight->SetLineColor(kBlue);
  h_pt_eff_ele25WPTight->SetMarkerColor(kBlue);

  h_pt_eff_ele25WPTight->Draw("pe same");

  TLegend* leg = new TLegend(0.6,0.2,0.9,0.5);
  leg->AddEntry(h_pt_eff_ele25WPTight,"Ele25_eta2p1_WPTight","pe");
  leg->Draw("same");

  c->SaveAs("trigeff_Ele25_eta2p1_WPTight.pdf");
  //c->SaveAs("trigeff_Ele25_eta2p1_WPTight.eps");

  // split into EB and EE
  
  TCanvas* c_ebee = new TCanvas("c_ebee","c_ebee");
  c_ebee->SetGrid(1,1);
  c_ebee->cd();

  h_axis->Draw();

  TEfficiency* h_pt_eff_eb_ele25WPTight = new TEfficiency(*h_pt_num_eb_ele25WPTight, *h_pt_denom_eb_ele25WPTight);
  h_pt_eff_eb_ele25WPTight->SetLineColor(kGreen+2);
  h_pt_eff_eb_ele25WPTight->SetMarkerColor(kGreen+2);

  TEfficiency* h_pt_eff_ee_ele25WPTight = new TEfficiency(*h_pt_num_ee_ele25WPTight, *h_pt_denom_ee_ele25WPTight);
  h_pt_eff_ee_ele25WPTight->SetLineColor(kMagenta);
  h_pt_eff_ee_ele25WPTight->SetMarkerColor(kMagenta);

  h_pt_eff_eb_ele25WPTight->Draw("pe same");
  h_pt_eff_ee_ele25WPTight->Draw("pe same");

  TLegend* leg_ebee = new TLegend(0.6,0.2,0.9,0.5);
  leg_ebee->AddEntry(h_pt_eff_eb_ele25WPTight,"Ele25_eta2p1_WPTight, EB only","pe");
  leg_ebee->AddEntry(h_pt_eff_ee_ele25WPTight,"Ele25_eta2p1_WPTight, EE only","pe");
  leg_ebee->Draw("same");

  c_ebee->SaveAs("trigeff_Ele25_eta2p1_WPTight_EBEE.pdf");
//  c_ebee->SaveAs("trigeff_Photon165_HE10_EBEE.eps");

  
  //Add more plots
  TCut passEle25WPTight = "HLT_Ele25_eta2p1_WPTight_Gsf>0";
  TCut passEle27WPLoose = "HLT_Ele27_eta2p1_WPLoose_Gsf>0";
  TCut sieieEB15 = "sigmaIEtaIEta_full5x5<0.01";
  TCut hoeEB15   = "hOverE<0.01";
//  TCut detaEB15  = "dEtaIn<0.002";
  TCut detaEB15  = "dEtaOut<0.004";
  TCut dphiEB15  = "dPhiIn<0.05";
  TCut eOpEB15   = "abs((1./ecalEnergy) - (eOverPIn/ecalEnergy))<0.004";
  TCut eIsoEB15  = "ecalPFClusterIso/p4.pt()<0.02";
  TCut hIsoEB15  = "hcalPFClusterIso/p4.pt()<0.02";
  TCut tIsoEB15  = "tkIso/p4.pt()<0.04";

  TCut ptCut ="p4.pt()>30";
  TCut pt30 ="p4.pt()>30";
  TCut base30 = base+ptCut;
  TCut EB15 = base30+eb+sieieEB15+hoeEB15+detaEB15+dphiEB15+eOpEB15+eIsoEB15+hIsoEB15+tIsoEB15;
  TCut EB15nopt = base+eb+sieieEB15+hoeEB15+detaEB15+dphiEB15+eOpEB15+eIsoEB15+hIsoEB15+tIsoEB15;
  TCut EB15noS = base30+eb+hoeEB15+detaEB15+dphiEB15+eOpEB15+eIsoEB15+hIsoEB15+tIsoEB15;
  TCut EB15noH = base30+eb+sieieEB15+detaEB15+dphiEB15+eOpEB15+eIsoEB15+hIsoEB15+tIsoEB15;
  TCut EB15noE = base30+eb+sieieEB15+hoeEB15+dphiEB15+eOpEB15+eIsoEB15+hIsoEB15+tIsoEB15;
  TCut EB15noP = base30+eb+sieieEB15+hoeEB15+detaEB15+eOpEB15+eIsoEB15+hIsoEB15+tIsoEB15;
  TCut EB15noO = base30+eb+sieieEB15+hoeEB15+detaEB15+dphiEB15+eIsoEB15+hIsoEB15+tIsoEB15;
  TCut EB15noEI = base30+eb+sieieEB15+hoeEB15+detaEB15+dphiEB15+eOpEB15+hIsoEB15+tIsoEB15;
  TCut EB15noHI = base30+eb+sieieEB15+hoeEB15+detaEB15+dphiEB15+eOpEB15+eIsoEB15+tIsoEB15;
  TCut EB15noTI = base30+eb+sieieEB15+hoeEB15+detaEB15+dphiEB15+eOpEB15+eIsoEB15+hIsoEB15;
  
  //void makeEff( TCut den, TCut num, int nbins, float min, float max, TString legend, TString xaxis, TString yaxis, TString saveas) {
  makeEff(t_ele, "p4.pt()", base+eb,base+eb+passEle25WPTight,60,20,80, "Ele25_eta2p1_WPTight, EB only", "el p_{T} [GeV]", "Efficiency of Ele25_eta2p1_WPTight", "ptEB_Ele25WPTight", false);
  makeEff(t_ele, "p4.pt()", base+ID+eb,base+ID+eb+passEle25WPTight,60,20,80, "Ele25_eta2p1_WPTight, EB only", "el p_{T} [GeV]", "Efficiency of Ele25_eta2p1_WPTight", "ptEB_Ele25WPTight_ID", false);
  makeEff(t_ele, "p4.pt()", base+ID+ee,base+ID+ee+passEle25WPTight,30,20,80, "Ele25_eta2p1_WPTight, EB only", "el p_{T} [GeV]", "Efficiency of Ele25_eta2p1_WPTight", "ptEE_Ele25WPTight_ID", false);
  
  makeEff(t_ele, "p4.pt()", base+ID+eb,base+ID+eb+passEle27WPLoose,60,20,80, "Ele27_eta2p1_WPLoose, EB only", "el p_{T} [GeV]", "Efficiency of Ele27_eta2p1_WPLoose", "ptEB_Ele27WPLoose_ID", false);
  makeEff(t_ele, "p4.pt()", base+ID+ee,base+ID+ee+passEle27WPLoose,30,20,80, "Ele27_eta2p1_WPLoose, EB only", "el p_{T} [GeV]", "Efficiency of Ele27_eta2p1_WPLoose", "ptEE_Ele27WPLoose_ID", false);

  
  makeEff(t_ele, "sigmaIEtaIEta_full5x5", EB15noS,EB15noS+passEle25WPTight,30,0,0.030, "Ele25_eta2p1_WPTight, EB only", "Sieie", "Efficiency of Ele25_eta2p1_WPTight", "IDsieieEB_Ele25WPTight_ID15");
  makeEff(t_ele, "hOverE", EB15noH,EB15noH+passEle25WPTight,21,0,0.21, "Ele25_eta2p1_WPTight, EB only", "H/E", "Efficiency of Ele25_eta2p1_WPTight", "IDhOverEEB_Ele25WPTight_ID15");
  makeEff(t_ele, "dEtaIn", EB15noE,EB15noE+passEle25WPTight,25,0,0.025, "Ele25_eta2p1_WPTight, EB only", "dEtaIn", "Efficiency of Ele25_eta2p1_WPTight", "IDdEtaInEB_Ele25WPTight_ID15");
  makeEff(t_ele, "dEtaOut", EB15noE,EB15noE+passEle25WPTight,30,0,0.015, "Ele25_eta2p1_WPTight, EB only", "dEtaOut", "Efficiency of Ele25_eta2p1_WPTight", "IDdEtaOutEB_Ele25WPTight_ID15");
  makeEff(t_ele, "dPhiIn", EB15noP,EB15noP+passEle25WPTight,20,0,0.1, "Ele25_eta2p1_WPTight, EB only", "dPhiIn", "Efficiency of Ele25_eta2p1_WPTight", "IDdPhiInEB_Ele25WPTight_ID15");
  makeEff(t_ele, "abs((1./ecalEnergy) - (eOverPIn/ecalEnergy))", EB15noO,EB15noO+passEle25WPTight,25,0,0.025, "Ele25_eta2p1_WPTight, EB only", "1/E-1/P", "Efficiency of Ele25_eta2p1_WPTight", "IDeOpEB_Ele25WPTight_ID15");
  makeEff(t_ele, "ecalPFClusterIso/p4.pt()", EB15noEI,EB15noEI+passEle25WPTight,20,0,0.4, "Ele25_eta2p1_WPTight, EB only", "ecalIso", "Efficiency of Ele25_eta2p1_WPTight", "IDeIsoEB_Ele25WPTight_ID15");
  makeEff(t_ele, "hcalPFClusterIso/p4.pt()", EB15noHI,EB15noHI+passEle25WPTight,20,0,0.4, "Ele25_eta2p1_WPTight, EB only", "hcalIso", "Efficiency of Ele25_eta2p1_WPTight", "IDhIsoEB_Ele25WPTight_ID15");
  makeEff(t_ele, "tkIso/p4.pt()", EB15noTI,EB15noTI+passEle25WPTight,20,0,0.4, "Ele25_eta2p1_WPTight, EB only", "tkIso", "Efficiency of Ele25_eta2p1_WPTight", "IDtIsoEB_Ele25WPTight_ID15");

  
  makeEff(t_ele, "sigmaIEtaIEta_full5x5", EB15noS,EB15noS+passEle27WPLoose,30,0,0.030, "Ele27_eta2p1_WPLoose, EB only", "Sieie", "Efficiency of Ele27_eta2p1_WPLoose", "IDsieieEB_Ele27WPLoose_ID15");
  makeEff(t_ele, "hOverE", EB15noH,EB15noH+passEle27WPLoose,21,0,0.21, "Ele27_eta2p1_WPLoose, EB only", "H/E", "Efficiency of Ele27_eta2p1_WPLoose", "IDhOverEEB_Ele27WPLoose_ID15");
  makeEff(t_ele, "dEtaIn", EB15noE,EB15noE+passEle27WPLoose,25,0,0.025, "Ele27_eta2p1_WPLoose, EB only", "dEtaIn", "Efficiency of Ele27_eta2p1_WPLoose", "IDdEtaInEB_Ele27WPLoose_ID15");
  makeEff(t_ele, "dEtaOut", EB15noE,EB15noE+passEle27WPLoose,30,0,0.015, "Ele27_eta2p1_WPLoose, EB only", "dEtaOut", "Efficiency of Ele27_eta2p1_WPLoose", "IDdEtaOutEB_Ele27WPLoose_ID15");
  makeEff(t_ele, "dPhiIn", EB15noP,EB15noP+passEle27WPLoose,20,0,0.1, "Ele27_eta2p1_WPLoose, EB only", "dPhiIn", "Efficiency of Ele27_eta2p1_WPLoose", "IDdPhiInEB_Ele27WPLoose_ID15");
  makeEff(t_ele, "abs((1./ecalEnergy) - (eOverPIn/ecalEnergy))", EB15noO,EB15noO+passEle27WPLoose,25,0,0.025, "Ele27_eta2p1_WPLoose, EB only", "1/E-1/P", "Efficiency of Ele27_eta2p1_WPLoose", "IDeOpEB_Ele27WPLoose_ID15");
  makeEff(t_ele, "ecalPFClusterIso/p4.pt()", EB15noEI,EB15noEI+passEle27WPLoose,20,0,0.4, "Ele27_eta2p1_WPLoose, EB only", "ecalIso", "Efficiency of Ele27_eta2p1_WPLoose", "IDeIsoEB_Ele27WPLoose_ID15");
  makeEff(t_ele, "hcalPFClusterIso/p4.pt()", EB15noHI,EB15noHI+passEle27WPLoose,20,0,0.4, "Ele27_eta2p1_WPLoose, EB only", "hcalIso", "Efficiency of Ele27_eta2p1_WPLoose", "IDhIsoEB_Ele27WPLoose_ID15");
  makeEff(t_ele, "tkIso/p4.pt()", EB15noTI,EB15noTI+passEle27WPLoose,20,0,0.4, "Ele27_eta2p1_WPLoose, EB only", "tkIso", "Efficiency of Ele27_eta2p1_WPLoose", "IDtIsoEB_Ele27WPLoose_ID15");

  //makeEff(t_ele, "p4.pt()", EB15nopt,EB15nopt+passEle25WPTight,60,20,80, "Ele25_eta2p1_WPTight, EB only", "el p_{T} [GeV]", "Efficiency of Ele25_eta2p1_WPTight", "ptEB_Ele25WPTight_ID15");
  makeOverlay(t_ele, "dilep_mass", tag+pt30+ID,tag+pt30+ID+passEle25WPTight,60,60,120, "Ele25_eta2p1_WPTight", "m_{ee} [GeV]", "Events", "Overlay_mEE_Ele25WPTight_ID");

  
  
  f_out->Write();
  f_out->Close();
  
  return;
}











