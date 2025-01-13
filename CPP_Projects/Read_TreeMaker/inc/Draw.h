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
   void SetOutputPDFFolder(const std::string &outputPDFfolder) { foutputPDFFolder = outputPDFfolder; }
   void SetOutputPDFrealpath(const std::string &outputPDFrealpath) { fRealpathPDF = outputPDFrealpath; }

   // Member functions
   void Run(const std::string &filepath);
   void CompareRuns(const std::vector<std::string> &v_filepaths, const std::string &type);

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
   std::vector<int> colors = {kRed + 2,     kBlue + 2, kGreen + 3,  kOrange + 7,
                              kMagenta + 2, kCyan + 2, kYellow + 2, kViolet + 2};
   std::vector<int> markers = {47, 34, 20, 21, 22, 23, 33, 43, 29};

   int resomin = 0;
   int resomax = 30;
   int meanmin = 200;
   int meanmax = 800;
   int stdmin = 0;
   int stdmax = 100;
};

#endif