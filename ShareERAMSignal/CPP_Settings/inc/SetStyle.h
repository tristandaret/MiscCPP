#ifndef SETSTYLE_H
#define SETSTYLE_H

#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TLegend.h"

// --- T2K style ---

TStyle *SetMyStyle();

// Set graphic settings of a TGraphErrors
void Graphic_setup(TGraphErrors *ptge, Size_t markersize, Style_t markerstyle,
                   Color_t markercolor, Width_t linewidth, Color_t linecolor);

// Set graphic settings of an histogram
void Graphic_setup(TH1 *th1, Size_t markersize, Style_t markerstyle, Color_t markercolor,
                   Width_t linewidth, Color_t linecolor, Color_t fillcolor = 0,
                   Float_t alpha = 1.0);

// Set graphic settings of an TH2
void Graphic_setup(TH2 *th2, Size_t markersize, Style_t markerstyle, Color_t markercolor);

// Set graphic settings for a TF1
void Graphic_setup(TF1 *tf1, Width_t linewidth = 2, Color_t linecolor = kBlack,
                   Style_t linestyle = 1);

// Set graphic settings for a TLegend
void Graphic_setup(TLegend *leg, Float_t textsize = 0.04, Color_t textcolor = kBlack,
                   Style_t fillStyle = 1001, Color_t fillcolor = kWhite,
                   Int_t bordersize = 0, Int_t ncolumns = 1);

// Set the stats box position
void SetStatBoxPosition(TH1 *pTH, const double &xmin, const double &xmax,
                        const double &ymin, const double &ymax);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
TStyle *SetT2KStyle(Int_t WhichStyle = 3, TString styleName = "T2K");

void CenterHistoTitles(TH1 *thisHisto);

void AddGridLinesToPad(TPad *thisPad);

#endif