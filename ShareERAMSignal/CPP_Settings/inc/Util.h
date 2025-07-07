#ifndef Util_H
#define Util_H

#include "Includes.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

#include <TStyle.h>
#include <TCanvas.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TLine.h>

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TMultiGraph.h>

TGraphErrors *Convert_TH2_TGE(const TH2 *pTH2);

int MakeMyDir(const std::string &DirName);

struct RankedValue {
   double Value;
   int Rank;

   bool operator<(const RankedValue &rhs) const { return Value < rhs.Value; }
};

#endif
