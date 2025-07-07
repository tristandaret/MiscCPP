#include <TLegend.h>
#include "TAxis.h"
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TGraph.h>
#include <TMath.h>
#include <TF1Convolution.h>
#include <TVirtualFFT.h>
#include "TROOT.h"
#include <TStyle.h>
#include <TPaveText.h>

using namespace std;

Double_t derivated_ER_ltspice(Double_t *x, Double_t *par)
{
    // double t0 = 0;   par[0]
     double Q = 0.6368;
   //  double Amp = 3000;   par[1]
     double ws = 0.1951;

     double Q_factor = (120.e-15)/(1.6e-19);
     double ADC_factor = 4096;
     double fmax1 = 0.1357;

     Double_t ADC = 0;
     Double_t prefactor = x[0]>par[0];
     Double_t t = x[0]-par[0];
     Double_t exp1 = TMath::Exp(-ws*t);
     Double_t exp2 = TMath::Exp(-ws*t/(2*Q));
     Double_t Qsqrt = TMath::Sqrt((2*Q-1)/(2*Q+1));
     Double_t arg = ws*0.5*TMath::Sqrt(4-1./pow(Q, 2));

     Double_t term1 = -ws*exp1;
     Double_t term2 = -ws/(2*Q) * exp2 * (Qsqrt * TMath::Sin(arg*t) - TMath::Cos(arg*t));
     Double_t term3 = exp2 * arg * (Qsqrt * TMath::Cos(arg*t*prefactor) + TMath::Sin(arg*t*prefactor));

     if (Q != 0 && Q != 0.5 && Q != -0.5)
       {
         ADC = prefactor * (ADC_factor/(Q_factor*fmax1)) * par[1] * (term1 + term2 + term3);  //(ADC_factor/(Q_factor*fmax1))
       }

     return ADC;
}

//-----------------------------------------------------------

Double_t Q_Lead(Double_t *x, Double_t *par)
{
    // double t0 ---   par[0]
    // double Q ---    par[1]
    // double x0 ---   par[2]
    // double y0 ---   par[3]
    // double RC ---   par[4]
    double LenX = 1.128;    // cm
    double LenY = 1.019;    // cm

    double Qtotal = 0;
    Double_t prefactor = x[0]>par[0];
    Double_t t = x[0]-par[0];
    double sigma = TMath::Sqrt(2*t*prefactor/par[4]);

    Qtotal = prefactor * (par[1]/4) * ( TMath::Erf((LenX/2 - par[2])/(sigma*TMath::Sqrt(2))) - TMath::Erf((-LenX/2 - par[2])/(sigma*TMath::Sqrt(2))) )
            * ( TMath::Erf((LenY/2 - par[3])/(sigma*TMath::Sqrt(2))) - TMath::Erf((-LenY/2 - par[3])/(sigma*TMath::Sqrt(2))) );
    return Qtotal;
}
// Double_t Q_LeftNeigh(Double_t *x, Double_t *par)
// {
//     double LenX = 1.128;    // cm
//     double LenY = 1.019;    // cm

//     double Qtotal = 0;
//     Double_t prefactor = x[0]>par[0];
//     Double_t t = x[0]-par[0];
//     double sigma = TMath::Sqrt(2*t*prefactor/par[4]);

//     Qtotal = prefactor * (par[1]/4) * ( TMath::Erf((-LenX/2 - par[2])/(sigma*TMath::Sqrt(2))) - TMath::Erf((-3*LenX/2 - par[2])/(sigma*TMath::Sqrt(2))) )
//             * ( TMath::Erf((LenY/2 - par[3])/(sigma*TMath::Sqrt(2))) - TMath::Erf((-LenY/2 - par[3])/(sigma*TMath::Sqrt(2))) );
//     return Qtotal;
// }
// Double_t Q_RightNeigh(Double_t *x, Double_t *par)
// {
//     double LenX = 1.128;    // cm
//     double LenY = 1.019;    // cm

//     double Qtotal = 0;
//     Double_t prefactor = x[0]>par[0];
//     Double_t t = x[0]-par[0];
//     double sigma = TMath::Sqrt(2*t*prefactor/par[4]);

//     Qtotal = prefactor * (par[1]/4) * ( TMath::Erf((3*LenX/2 - par[2])/(sigma*TMath::Sqrt(2))) - TMath::Erf((LenX/2 - par[2])/(sigma*TMath::Sqrt(2))) )
//             * ( TMath::Erf((LenY/2 - par[3])/(sigma*TMath::Sqrt(2))) - TMath::Erf((-LenY/2 - par[3])/(sigma*TMath::Sqrt(2))) );
//     return Qtotal;
// }
// Double_t Q_TopNeigh(Double_t *x, Double_t *par)
// {
//     double LenX = 1.128;    // cm
//     double LenY = 1.019;    // cm

//     double Qtotal = 0;
//     Double_t prefactor = x[0]>par[0];
//     Double_t t = x[0]-par[0];
//     double sigma = TMath::Sqrt(2*t*prefactor/par[4]);

//     Qtotal = prefactor * (par[1]/4) * ( TMath::Erf((LenX/2 - par[2])/(sigma*TMath::Sqrt(2))) - TMath::Erf((-LenX/2 - par[2])/(sigma*TMath::Sqrt(2))) )
//             * ( TMath::Erf((3*LenY/2 - par[3])/(sigma*TMath::Sqrt(2))) - TMath::Erf((LenY/2 - par[3])/(sigma*TMath::Sqrt(2))) );
//     return Qtotal;
// }
// Double_t Q_BottomNeigh(Double_t *x, Double_t *par)
// {
//     double LenX = 1.128;    // cm
//     double LenY = 1.019;    // cm

//     double Qtotal = 0;
//     Double_t prefactor = x[0]>par[0];
//     Double_t t = x[0]-par[0];
//     double sigma = TMath::Sqrt(2*t*prefactor/par[4]);

//     Qtotal = prefactor * (par[1]/4) * ( TMath::Erf((LenX/2 - par[2])/(sigma*TMath::Sqrt(2))) - TMath::Erf((-LenX/2 - par[2])/(sigma*TMath::Sqrt(2))) )
//             * ( TMath::Erf((-LenY/2 - par[3])/(sigma*TMath::Sqrt(2))) - TMath::Erf((-3*LenY/2 - par[3])/(sigma*TMath::Sqrt(2))) );
//     return Qtotal;
// }
// Double_t Q_TLNeigh(Double_t *x, Double_t *par)
// {
//     double LenX = 1.128;    // cm
//     double LenY = 1.019;    // cm

//     double Qtotal = 0;
//     Double_t prefactor = x[0]>par[0];
//     Double_t t = x[0]-par[0];
//     double sigma = TMath::Sqrt(2*t*prefactor/par[4]);

//     Qtotal = prefactor * (par[1]/4) * ( TMath::Erf((-LenX/2 - par[2])/(sigma*TMath::Sqrt(2))) - TMath::Erf((-3*LenX/2 - par[2])/(sigma*TMath::Sqrt(2))) )
//             * ( TMath::Erf((3*LenY/2 - par[3])/(sigma*TMath::Sqrt(2))) - TMath::Erf((LenY/2 - par[3])/(sigma*TMath::Sqrt(2))) );
//     return Qtotal;
// }
// Double_t Q_TRNeigh(Double_t *x, Double_t *par)
// {
//     double LenX = 1.128;    // cm
//     double LenY = 1.019;    // cm

//     double Qtotal = 0;
//     Double_t prefactor = x[0]>par[0];
//     Double_t t = x[0]-par[0];
//     double sigma = TMath::Sqrt(2*t*prefactor/par[4]);

//     Qtotal = prefactor * (par[1]/4) * ( TMath::Erf((3*LenX/2 - par[2])/(sigma*TMath::Sqrt(2))) - TMath::Erf((LenX/2 - par[2])/(sigma*TMath::Sqrt(2))) )
//             * ( TMath::Erf((3*LenY/2 - par[3])/(sigma*TMath::Sqrt(2))) - TMath::Erf((LenY/2 - par[3])/(sigma*TMath::Sqrt(2))) );
//     return Qtotal;
// }
// Double_t Q_BLNeigh(Double_t *x, Double_t *par)
// {
//     double LenX = 1.128;    // cm
//     double LenY = 1.019;    // cm

//     double Qtotal = 0;
//     Double_t prefactor = x[0]>par[0];
//     Double_t t = x[0]-par[0];
//     double sigma = TMath::Sqrt(2*t*prefactor/par[4]);

//     Qtotal = prefactor * (par[1]/4) * ( TMath::Erf((-LenX/2 - par[2])/(sigma*TMath::Sqrt(2))) - TMath::Erf((-3*LenX/2 - par[2])/(sigma*TMath::Sqrt(2))) )
//             * ( TMath::Erf((-LenY/2 - par[3])/(sigma*TMath::Sqrt(2))) - TMath::Erf((-3*LenY/2 - par[3])/(sigma*TMath::Sqrt(2))) );
//     return Qtotal;
// }
// Double_t Q_BRNeigh(Double_t *x, Double_t *par)
// {
//     double LenX = 1.128;    // cm
//     double LenY = 1.019;    // cm

//     double Qtotal = 0;
//     Double_t prefactor = x[0]>par[0];
//     Double_t t = x[0]-par[0];
//     double sigma = TMath::Sqrt(2*t*prefactor/par[4]);

//     Qtotal = prefactor * (par[1]/4) * ( TMath::Erf((3*LenX/2 - par[2])/(sigma*TMath::Sqrt(2))) - TMath::Erf((LenX/2 - par[2])/(sigma*TMath::Sqrt(2))) )
//             * ( TMath::Erf((-LenY/2 - par[3])/(sigma*TMath::Sqrt(2))) - TMath::Erf((-3*LenY/2 - par[3])/(sigma*TMath::Sqrt(2))) );
//     return Qtotal;
// }

// // // // //  // // // // //  // // // // //  // // // // //  // // // // //  // // // // //

void toyMC()
{
    TCanvas* c1 = new TCanvas("c1", "c1", 1000, 1000);
    c1->Divide(3,3);

    TCanvas* c4 = new TCanvas("c4", "c4", 1000, 1000);
    c4->Divide(3,3);

    TCanvas* c2 = new TCanvas();
    TCanvas* c3 = new TCanvas();
    TCanvas* c5 = new TCanvas();

    TH1D *hlead = new TH1D("hlead", "hlead", 510, 0, 510);
    TH1D *hLN = new TH1D("hLN", "hLN", 510, 0, 510);
    TH1D *hRN = new TH1D("hRN", "hRN", 510, 0, 510);
    TH1D *hTN = new TH1D("hTN", "hTN", 510, 0, 510);
    TH1D *hBN = new TH1D("hBN", "hBN", 510, 0, 510);
    TH1D *hTLN = new TH1D("hTLN", "hTLN", 510, 0, 510);
    TH1D *hTRN = new TH1D("hTRN", "hTRN", 510, 0, 510);
    TH1D *hBLN = new TH1D("hBLN", "hBLN", 510, 0, 510);
    TH1D *hBRN = new TH1D("hBRN", "hBRN", 510, 0, 510);

    TH1D *hsum = new TH1D("hsum", "hsum", 510, 0, 510);

    double LenX = 1.128;    // cm
    double LenY = 1.019;    // cm
    double Qinc = 336000;
    Double_t RC = 10000/40;    //time bin/cm2
    Double_t x0 = LenX/3;
    Double_t y0 = -LenY/3;

    TF1 *d_ER_f1 = new TF1("d_ER_f1", derivated_ER_ltspice, 0, 50, 2);
    d_ER_f1->SetParameter(0, 0);
    d_ER_f1->SetParameter(1, 1);

  //-------------------------------------------------------------------------------
  //Charge in the leading Pad

    TF1 *Q_Leading = new TF1("Q_Leading", Q_Lead, 0., 300, 5);
    Q_Leading->SetParameter(0, 0);
    Q_Leading->SetParameter(1, Qinc);
    Q_Leading->SetParameter(2, x0);
    Q_Leading->SetParameter(3, y0);
    Q_Leading->SetParameter(4, RC);

    c4->cd(5)->SetGrid();
    Q_Leading->SetLineColor(2);
    Q_Leading->Draw();

    TF1Convolution *conv_Q_dER_Lead = new TF1Convolution("d_ER_f1", "Q_Leading", -510, 510, true);
    conv_Q_dER_Lead->SetRange(-510, 510);
    conv_Q_dER_Lead->SetNofPointsFFT(10000);

    TF1 *fLead = new TF1("Leading pad", *conv_Q_dER_Lead, 0., 100., conv_Q_dER_Lead->GetNpar());
    fLead->SetParameter(0, 0);    // t0 --- ER
    fLead->SetParameter(1, 1);    // Amp --- ER
    fLead->SetParameter(2, 0);   // t0 --- Q_pad
    fLead->SetParameter(3, Qinc);   // Charge --- Q_pad
    fLead->SetParameter(4, x0);   // x0 --- Q_pad
    fLead->SetParameter(5, y0);   // y0 --- Q_pad
    fLead->SetParameter(6, RC);   // RC --- Q_pad

    c1->cd(5)->SetGrid();
    fLead->SetLineColor(2);
    fLead->GetXaxis()->SetTitle("time bins (40ns)");
    fLead->GetYaxis()->SetTitle("Signal amplitude (ADC)");
  //  fLead->GetYaxis()->SetRangeUser(-25, 100);
    fLead->Draw();

    cout << fLead->GetMaximum() << endl;

    c2->cd();
    fLead->SetLineColor(2);
    fLead->Draw();
    //--------------------------------------------------------------------------------------------------------------------------

//     TF1 *Q_LeftN = new TF1("Q_LeftN", Q_LeftNeigh, 0, 300, 5);
//     Q_LeftN->SetParameter(0, 0);
//     Q_LeftN->SetParameter(1, Qinc);
//     Q_LeftN->SetParameter(2, x0);
//     Q_LeftN->SetParameter(3, y0);
//     Q_LeftN->SetParameter(4, RC);

//     c4->cd(4)->SetGrid();
//     Q_LeftN->SetLineColor(4);
//     Q_LeftN->Draw();

//     TF1Convolution *conv_Q_dER_LN = new TF1Convolution("d_ER_f1", "Q_LeftN", -510, 510, true);
//     conv_Q_dER_LN->SetRange(-510, 510);
//     conv_Q_dER_LN->SetNofPointsFFT(10000);

//     TF1 *fLN = new TF1("Left neigh.", *conv_Q_dER_LN, 0., 100., conv_Q_dER_LN->GetNpar());
//     fLN->SetParameter(0, 0);
//     fLN->SetParameter(1, 1);
//     fLN->SetParameter(2, 0);
//     fLN->SetParameter(3, Qinc);
//     fLN->SetParameter(4, x0);
//     fLN->SetParameter(5, y0);
//     fLN->SetParameter(6, RC);

//     c2->cd();
//     fLN->SetLineColor(4);
//     fLN->Draw("SAME");

//     c1->cd(4)->SetGrid();
//     fLN->SetLineColor(4);
//     fLN->GetXaxis()->SetTitle("time bins (40ns)");
//     fLN->GetYaxis()->SetTitle("Signal amplitude (ADC)");
//   //  fLN->GetYaxis()->SetRangeUser(-25, 100);
//     fLN->Draw();



//     //--------------------------------------------------------------------------------------------------------------------------

//     TF1 *Q_RightN = new TF1("Q_RightN", Q_RightNeigh, 0, 300, 5);
//     Q_RightN->SetParameter(0, 0);
//     Q_RightN->SetParameter(1, Qinc);
//     Q_RightN->SetParameter(2, x0);
//     Q_RightN->SetParameter(3, y0);
//     Q_RightN->SetParameter(4, RC);

//     c4->cd(6)->SetGrid();
//     Q_RightN->SetLineColor(4);
//     Q_RightN->Draw();

//     TF1Convolution *conv_Q_dER_RN = new TF1Convolution("d_ER_f1", "Q_RightN", -510, 510, true);
//     conv_Q_dER_RN->SetRange(-510, 510);
//     conv_Q_dER_RN->SetNofPointsFFT(10000);

//     TF1 *fRN = new TF1("Right neigh.", *conv_Q_dER_RN, 0., 100., conv_Q_dER_RN->GetNpar());
//     fRN->SetParameter(0, 0);
//     fRN->SetParameter(1, 1);
//     fRN->SetParameter(2, 0);
//     fRN->SetParameter(3, Qinc);
//     fRN->SetParameter(4, x0);
//     fRN->SetParameter(5, y0);
//     fRN->SetParameter(6, RC);

//     c1->cd(6)->SetGrid();
//     fRN->SetLineColor(4);
//     fRN->GetXaxis()->SetTitle("time bins (40ns)");
//     fRN->GetYaxis()->SetTitle("Signal amplitude (ADC)");
//  //   fRN->GetYaxis()->SetRangeUser(-25, 100);
//     fRN->Draw();

//     cout << fRN->GetMaximum() << endl;

//     c2->cd();
//     fRN->Draw("SAME");

//     //--------------------------------------------------------------------------------------------------------------------------

//     TF1 *Q_TopN = new TF1("Q_TopN", Q_TopNeigh, 0, 300, 5);
//     Q_TopN->SetParameter(0, 0);
//     Q_TopN->SetParameter(1, Qinc);
//     Q_TopN->SetParameter(2, x0);
//     Q_TopN->SetParameter(3, y0);
//     Q_TopN->SetParameter(4, RC);

//     c4->cd(2)->SetGrid();
//     Q_TopN->SetLineColor(4);
//     Q_TopN->Draw();

//     TF1Convolution *conv_Q_dER_TN = new TF1Convolution("d_ER_f1", "Q_TopN", -510, 510, true);
//     conv_Q_dER_TN->SetRange(-510, 510);
//     conv_Q_dER_TN->SetNofPointsFFT(10000);

//     TF1 *fTN = new TF1("Top neigh.", *conv_Q_dER_TN, 0., 100., conv_Q_dER_TN->GetNpar());
//     fTN->SetParameter(0, 0);
//     fTN->SetParameter(1, 1);
//     fTN->SetParameter(2, 0);
//     fTN->SetParameter(3, Qinc);
//     fTN->SetParameter(4, x0);
//     fTN->SetParameter(5, y0);
//     fTN->SetParameter(6, RC);

//     c1->cd(2)->SetGrid();
//     fTN->SetLineColor(4);
//     fTN->GetXaxis()->SetTitle("time bins (40ns)");
//     fTN->GetYaxis()->SetTitle("Signal amplitude (ADC)");
//   //  fTN->GetYaxis()->SetRangeUser(-25, 100);
//     fTN->Draw();

//     c2->cd();
//     fTN->Draw("SAME");

//     //--------------------------------------------------------------------------------------------------------------------------

//     TF1 *Q_BottomN = new TF1("Q_BottomN", Q_BottomNeigh, 0, 300, 5);
//     Q_BottomN->SetParameter(0, 0);
//     Q_BottomN->SetParameter(1, Qinc);
//     Q_BottomN->SetParameter(2, x0);
//     Q_BottomN->SetParameter(3, y0);
//     Q_BottomN->SetParameter(4, RC);

//     c4->cd(8)->SetGrid();
//     Q_BottomN->SetLineColor(4);
//     Q_BottomN->Draw();

//     TF1Convolution *conv_Q_dER_BN = new TF1Convolution("d_ER_f1", "Q_BottomN", -510, 510, true);
//     conv_Q_dER_BN->SetRange(-510, 510);
//     conv_Q_dER_BN->SetNofPointsFFT(10000);

//     TF1 *fBN = new TF1("Bottom neigh.", *conv_Q_dER_BN, 0., 100., conv_Q_dER_BN->GetNpar());
//     fBN->SetParameter(0, 0);
//     fBN->SetParameter(1, 1);
//     fBN->SetParameter(2, 0);
//     fBN->SetParameter(3, Qinc);
//     fBN->SetParameter(4, x0);
//     fBN->SetParameter(5, y0);
//     fBN->SetParameter(6, RC);

//     c1->cd(8)->SetGrid();
//     fBN->SetLineColor(4);
//     fBN->GetXaxis()->SetTitle("time bins (40ns)");
//     fBN->GetYaxis()->SetTitle("Signal amplitude (ADC)");
//  //   fBN->GetYaxis()->SetRangeUser(-25, 100);
//     fBN->Draw();

//     cout << fBN->GetMaximum() << endl;

//     c2->cd();
//     fBN->Draw("SAME");

//     //--------------------------------------------------------------------------------------------------------------------------

//     TF1 *Q_TopLeftN = new TF1("Q_TopLeftN", Q_TLNeigh, 0, 300, 5);
//     Q_TopLeftN->SetParameter(0, 0);
//     Q_TopLeftN->SetParameter(1, Qinc);
//     Q_TopLeftN->SetParameter(2, x0);
//     Q_TopLeftN->SetParameter(3, y0);
//     Q_TopLeftN->SetParameter(4, RC);

//     c4->cd(1)->SetGrid();
//     Q_TopLeftN->SetLineColor(3);
//     Q_TopLeftN->Draw();

//     TF1Convolution *conv_Q_dER_TLN = new TF1Convolution("d_ER_f1", "Q_TopLeftN", -510, 510, true);
//     conv_Q_dER_TLN->SetRange(-510, 510);
//     conv_Q_dER_TLN->SetNofPointsFFT(10000);

//     TF1 *fTLN = new TF1("Top-Left neigh.", *conv_Q_dER_TLN, 0., 100., conv_Q_dER_TLN->GetNpar());
//     fTLN->SetParameter(0, 0);
//     fTLN->SetParameter(1, 1);
//     fTLN->SetParameter(2, 0);
//     fTLN->SetParameter(3, Qinc);
//     fTLN->SetParameter(4, x0);
//     fTLN->SetParameter(5, y0);
//     fTLN->SetParameter(6, RC);

//     c1->cd(1)->SetGrid();
//     fTLN->SetLineColor(3);
//     fTLN->GetXaxis()->SetTitle("time bins (40ns)");
//     fTLN->GetYaxis()->SetTitle("Signal amplitude (ADC)");
//  //   fTLN->GetYaxis()->SetRangeUser(-25, 100);
//     fTLN->Draw();

//     c2->cd();
//     fTLN->SetLineColor(3);
//     fTLN->Draw("SAME");

//     //--------------------------------------------------------------------------------------------------------------------------

//     TF1 *Q_TopRightN = new TF1("Q_TopRightN", Q_TRNeigh, 0, 300, 5);
//     Q_TopRightN->SetParameter(0, 0);
//     Q_TopRightN->SetParameter(1, Qinc);
//     Q_TopRightN->SetParameter(2, x0);
//     Q_TopRightN->SetParameter(3, y0);
//     Q_TopRightN->SetParameter(4, RC);

//     c4->cd(3)->SetGrid();
//     Q_TopRightN->SetLineColor(3);
//     Q_TopRightN->Draw();

//     TF1Convolution *conv_Q_dER_TRN = new TF1Convolution("d_ER_f1", "Q_TopRightN", -510, 510, true);
//     conv_Q_dER_TRN->SetRange(-510, 510);
//     conv_Q_dER_TRN->SetNofPointsFFT(10000);

//     TF1 *fTRN = new TF1("Top-Right neigh.", *conv_Q_dER_TRN, 0., 100., conv_Q_dER_TRN->GetNpar());
//     fTRN->SetParameter(0, 0);
//     fTRN->SetParameter(1, 1);
//     fTRN->SetParameter(2, 0);
//     fTRN->SetParameter(3, Qinc);
//     fTRN->SetParameter(4, x0);
//     fTRN->SetParameter(5, y0);
//     fTRN->SetParameter(6, RC);

//     c1->cd(3)->SetGrid();
//     fTRN->SetLineColor(3);
//     fTRN->GetXaxis()->SetTitle("time bins (40ns)");
//     fTRN->GetYaxis()->SetTitle("Signal amplitude (ADC)");
//  //   fTRN->GetYaxis()->SetRangeUser(-25, 100);
//     fTRN->Draw();

//     c2->cd();
//     fTRN->Draw("SAME");

//     //--------------------------------------------------------------------------------------------------------------------------

//     TF1 *Q_BottomLeftN = new TF1("Q_BottomLeftN", Q_BLNeigh, 0, 300, 5);
//     Q_BottomLeftN->SetParameter(0, 0);
//     Q_BottomLeftN->SetParameter(1, Qinc);
//     Q_BottomLeftN->SetParameter(2, x0);
//     Q_BottomLeftN->SetParameter(3, y0);
//     Q_BottomLeftN->SetParameter(4, RC);

//     c4->cd(7)->SetGrid();
//     Q_BottomLeftN->SetLineColor(3);
//     Q_BottomLeftN->Draw();

//     TF1Convolution *conv_Q_dER_BLN = new TF1Convolution("d_ER_f1", "Q_BottomLeftN", -510, 510, true);
//     conv_Q_dER_BLN->SetRange(-510, 510);
//     conv_Q_dER_BLN->SetNofPointsFFT(10000);

//     TF1 *fBLN = new TF1("Bottom-Left neigh.", *conv_Q_dER_BLN, 0., 100., conv_Q_dER_BLN->GetNpar());
//     fBLN->SetParameter(0, 0);
//     fBLN->SetParameter(1, 1);
//     fBLN->SetParameter(2, 0);
//     fBLN->SetParameter(3, Qinc);
//     fBLN->SetParameter(4, x0);
//     fBLN->SetParameter(5, y0);
//     fBLN->SetParameter(6, RC);

//     c1->cd(7)->SetGrid();
//     fBLN->SetLineColor(3);
//     fBLN->GetXaxis()->SetTitle("time bins (40ns)");
//     fBLN->GetYaxis()->SetTitle("Signal amplitude (ADC)");
//  //   fBLN->GetYaxis()->SetRangeUser(-25, 100);
//     fBLN->Draw();

//     c2->cd();
//     fBLN->Draw("SAME");

//     //--------------------------------------------------------------------------------------------------------------------------

//     TF1 *Q_BottomRightN = new TF1("Q_BottomRightN", Q_BRNeigh, 0, 300, 5);
//     Q_BottomRightN->SetParameter(0, 0);
//     Q_BottomRightN->SetParameter(1, Qinc);
//     Q_BottomRightN->SetParameter(2, x0);
//     Q_BottomRightN->SetParameter(3, y0);
//     Q_BottomRightN->SetParameter(4, RC);

//     c4->cd(9)->SetGrid();
//     Q_BottomRightN->SetLineColor(3);
//     Q_BottomRightN->Draw();

//     TF1Convolution *conv_Q_dER_BRN = new TF1Convolution("d_ER_f1", "Q_BottomRightN", -510, 510, true);
//     conv_Q_dER_BRN->SetRange(-510, 510);
//     conv_Q_dER_BRN->SetNofPointsFFT(10000);

//     TF1 *fBRN = new TF1("Bottom-Right neigh.", *conv_Q_dER_BRN, 0., 100., conv_Q_dER_BRN->GetNpar());
//     fBRN->SetParameter(0, 0);
//     fBRN->SetParameter(1, 1);
//     fBRN->SetParameter(2, 0);
//     fBRN->SetParameter(3, Qinc);
//     fBRN->SetParameter(4, x0);
//     fBRN->SetParameter(5, y0);
//     fBRN->SetParameter(6, RC);

//     c1->cd(9)->SetGrid();
//     fBRN->SetLineColor(3);
//     fBRN->GetXaxis()->SetTitle("time bins (40ns)");
//     fBRN->GetYaxis()->SetTitle("Signal amplitude (ADC)");
//  //   fBRN->GetYaxis()->SetRangeUser(-25, 100);
//     fBRN->Draw();

//     cout << fBRN->GetMaximum() << endl;

//     //--------------------------------------------------------------------------------------------------------------------------
//     //--------------------------------------------------------------------------------------------------------------------------
//     //--------------------------------------------------------------------------------------------------------------------------

//     for(int i=0; i<510; i++)
//     {
//         hsum->SetBinContent(i+1, 0);

//         hlead->SetBinContent(i+1, fLead->Eval(i+0.5));
//         hLN->SetBinContent(i+1, fLN->Eval(i+0.5));
//         hRN->SetBinContent(i+1, fRN->Eval(i+0.5));
//         hTN->SetBinContent(i+1, fTN->Eval(i+0.5));
//         hBN->SetBinContent(i+1, fBN->Eval(i+0.5));
//         hTLN->SetBinContent(i+1, fTLN->Eval(i+0.5));
//         hTRN->SetBinContent(i+1, fTRN->Eval(i+0.5));
//         hBLN->SetBinContent(i+1, fBLN->Eval(i+0.5));
//         hBRN->SetBinContent(i+1, fBRN->Eval(i+0.5));
//     }

//     hsum->Add(hlead, 1);
//     hsum->Add(hLN, 1);
//     hsum->Add(hRN, 1);
//     hsum->Add(hTN, 1);
//     hsum->Add(hBN, 1);
//     hsum->Add(hTLN, 1);
//     hsum->Add(hTRN, 1);
//     hsum->Add(hBLN, 1);
//     hsum->Add(hBRN, 1);

//    // c2->cd();
//   //  fBRN->Draw("SAME");

//     c3->cd();
//     hsum->Draw();

//     c5->cd()->SetGrid();
//     Q_Leading->SetLineColor(2);
//  //   Q_Leading->SetLineWidth(3);
//     Q_Leading->Draw();
//     Q_TopN->SetLineColor(4);
//  //   Q_TopN->SetLineWidth(3);
//     Q_TopN->Draw("same");
//     Q_TopRightN->SetLineColor(3);
//  //   Q_TopRightN->SetLineWidth(3);
//     Q_TopRightN->Draw("same");

//     cout << fLead->GetMaximum() + fBN->GetMaximum() + fRN->GetMaximum()  + fTN->GetMaximum() + fLN->GetMaximum()
//             + fBRN->GetMaximum() + fBLN->GetMaximum() + fTRN->GetMaximum() + fTLN->GetMaximum() << endl;

//     cout << "Gain = " << hsum->GetMaximum()*183/224 << endl;


}


