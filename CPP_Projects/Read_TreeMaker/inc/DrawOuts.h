#include "Util.h"

class Draw {
public:
   // Constructor
   Draw();
   // Destructor
   ~Draw();
   // Member functions
   void Run(const std::string &filepath);
   void Scan(const std::vector<std::string> &v_filepath);
   void CompareScans(const std::vector<std::string> &v_tags, const std::string &method);

   std::string comment;
   std::vector<std::string> v_comment;
   std::string type;
   std::string run;
   std::string fileName;

private:
   // Data file
   TFile *fpInputFile;
   TTree *fpInputTree;

   // Histograms
   // Scans
   TH1F *fph1f_WF;
   TH1F *fph1f_XP;
   TF1 *fptf1_WF;
   TF1 *fptf1_XP;
   TGraphErrors *fpTGE_mean_WF;
   TGraphErrors *fpTGE_mean_XP;
   TGraphErrors *fpTGE_std_WF;
   TGraphErrors *fpTGE_std_XP;
   TGraphErrors *fpTGE_reso_WF;
   TGraphErrors *fpTGE_reso_XP;

   // Settings
   std::vector<int> v_fvalues;
   TStyle *fpStyle;
   TCanvas *fpCanvas;
   TLegend *fpLegend;
   std::vector<int> colors = {kRed + 2, kMagenta + 2, kBlue + 2};
   std::vector<int> markers = {22, 47, 23};
   std::vector<int> colors2 = {kBlue + 2, kRed + 2};
};