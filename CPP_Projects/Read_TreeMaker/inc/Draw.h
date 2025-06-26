#ifndef DRAW_H
#define DRAW_H

#include "Util.h"
#include "Process.h"

class Draw {
public:
   Draw();
   ~Draw();

   // Setters
   void SetOutputComparisonFolder(const std::string &outputComparisonFolder)
   {
      foutputComparisonFolder = outputComparisonFolder;
   }
   void SetOutputPDFFolder(const std::string &outputPDFfolder)
   {
      foutputPDFFolder = outputPDFfolder;
   }
   void SetOutputPDFrealpath(const std::string &outputPDFrealpath)
   {
      fRealpathPDF = outputPDFrealpath;
   }

   // Member functions
   void Run(const std::string &filepath);
   void Compare(const std::vector<std::string> &v_filepaths, const std::string &type);
   void CompareComments(const std::vector<std::string> &v_filepaths);

private:
   std::string frun;
   std::string ftag;
   std::string foutputPDFFolder;
   std::string foutputComparisonFolder;
   std::string fRealpathPDF;

   TFile fpInputTFile;
   TTree fpInputTree;

   // Settings
   std::vector<int> v_fvalues;
   TStyle *fpStyle;
   TCanvas *fpCanvas;
   TLegend *fpLegend;
   std::vector<int> colorsXP = {kOrange + 1, kBlue + 1, kOrange + 2, kBlue + 1};
   std::vector<int> markers = {20, 34, 23, 34};
   std::vector<int> colorsWF = {kRed + 1, kGreen + 1, kRed + 3, kGreen + 1};
   std::vector<int> colors = {kOrange + 1,  kBlue + 1, kRed + 1,    kGreen + 1,
                              kMagenta + 1, kCyan + 1, kYellow + 1, kViolet + 1};
   // std::vector<int> markers = {47, 34, 20, 21, 22, 23, 33, 43, 29};
   // std::vector<int> markers = {20, 34, 20, 21, 22, 23, 33, 43, 29};

   float resomin = 4;
   float resomax = 11;
   int meanmin = 300;
   int meanmax = 600;
   int stdmin = 0;
   int stdmax = 100;

   // Margins
   float marginTop = 0.06;
   float marginBottom = 0.15;
   float marginLeft = 0.14;
   float marginRighth1 = 0.045;
   float marginRighth2 = 0.15;
   void SetMarginH1(TPad *pad)
   {
      pad->SetMargin(marginLeft, marginRighth1, marginBottom, marginTop);
   }
   void SetMarginH2(TPad *pad)
   {
      pad->SetMargin(marginLeft, marginRighth2, marginBottom, marginTop);
   }
};

#endif