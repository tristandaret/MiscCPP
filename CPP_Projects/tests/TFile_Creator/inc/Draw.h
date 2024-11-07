#ifndef DRAW_H
#define DRAW_H

#include "Util.h"
#include "Process.h"

class Draw
{
    public:
		Draw();
		~Draw();

      // Setters
      void SetOutputPDFfolder(const std::string &outputPDFfolder){foutputPDFfolder = outputPDFfolder;}
      void SetOutputPDFrealpath(const std::string &outputPDFrealpath){fRealpathPDF = outputPDFrealpath;}

		// Member functions
		void CompareRuns(const std::vector<std::string> &v_filepaths);

	private:
		std::string frun;
		std::string ftag;
		std::string foutputPDFfolder;
		std::string fRealpathPDF;

		TFile fpInputTFile;
		TTree fpInputTree;

		// Settings
		std::vector<int> v_fvalues;
		TStyle *fpStyle;
		TCanvas *fpCanvas;
		TLegend *fpLegend;
		std::vector<int> colors = {kRed+2, kBlue+2, kGreen+2, kOrange+2, kMagenta+2, kCyan+2, kYellow+2, kViolet+2};
		std::vector<int> markers = {47, 34, 20, 21, 22, 23, 33, 43, 29};

		int resomin = 0;
		int resomax = 25;
		int meanmin = 200;
		int meanmax = 800;
		int stdmin = 0;
		int stdmax = 100;
};

#endif