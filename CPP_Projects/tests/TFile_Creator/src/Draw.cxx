#include "Draw.h"
#include "Misc_Functions.h"
#include "TFrame.h"
#include "SetStyle.h"

#include <fstream>

// Constructor
Draw::Draw()
{
   fpStyle = SetMyStyle();
   gROOT->SetStyle(fpStyle->GetName());
   fpCanvas = new TCanvas("fpCanvas", "fpCanvas", 1800, 1350);
   gPad->UseCurrentStyle();
}

// Destructor
Draw::~Draw() {}

void Draw::CompareRuns(const std::vector<std::string> &v_filepaths)
{

   // Load processes
   std::vector<Process *> v_processes;
   for (int i = 0; i < v_filepaths.size(); i++) {
      TFile inputFile(v_filepaths[i].c_str(), "READ");
      TObject *obj = inputFile.Get("Process"); // delete dealt with by ROOT
      Process *process = dynamic_cast<Process *>(obj);
      v_processes.push_back(process);
   }

   for (int i = 0; i < v_filepaths.size(); i++) {
      std::cout << v_filepaths[i] << std::endl;
      std::cout << v_processes[i]->frun << " " << v_processes[i]->ftag << " " << v_processes[i]->fcomment << " "
                << v_processes[i]->fph1f_WF->GetMean() << std::endl;
   }

   for (Process *pr : v_processes)
      delete pr;
}
