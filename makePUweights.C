void makePUweights(TString outName = "puWeight_moriond_Prompt.root", TString inData = "2016SingleEl_Feb22_moriond.root", TString inMC = "DY_Feb22_noPuReweight_moriond.root")
{

  TFile * fOut = new TFile(outName, "RECREATE");
  TFile * fData = new TFile(inData);
  TFile * fMC = new TFile(inMC);

  fOut->cd();

  TH1D* hOut   = (TH1D*) fMC->Get("hZprobe_nvtxEB")->Clone();
  TH1D* hMC    = (TH1D*) fMC->Get("hZprobe_nvtxEB")->Clone();  
  TH1D* hData  = (TH1D*) fData->Get("hZprobe_nvtxEB")->Clone();

  hOut->Reset();
  hOut->SetName("h_dataOverMC_nvtxEB");

  hOut->Divide(hData,hMC);

  hOut->Write();

  fOut->Write();
}
