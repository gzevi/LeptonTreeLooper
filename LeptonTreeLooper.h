#ifndef LEPTONTREELOOPER_H
#define LEPTONTREELOOPER_H

// C++
#include <iostream>
#include <vector>

// ROOT
#include "TBenchmark.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTreeCache.h"

#include "TLorentzVector.h"
#include "Math/Vector4D.h"
#include "Math/LorentzVector.h"

// LeptonTree
#include "LeptonTree.h"

#include "PlotUtilities.h"


using namespace std;
using namespace lepton_tree;

//int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test");

int LeptonTreeLooper(TChain* chain, TString output_name = "test", int nEvents = -1);
  
float pt_corrected_;
float E_corrected_;

#endif
